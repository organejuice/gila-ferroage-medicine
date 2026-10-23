"""Reproduce Chhatbar source-level G6 audit from immutable public files.

Run the companion R scripts to export Seurat metadata, BA17 state detection,
per-donor QC, independent-marker diagnostics and geometry linkage. This script
never computes AD–control effects.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
from collections import defaultdict
from pathlib import Path

from openpyxl import load_workbook


ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data/raw/G6_chhatbar_metadata_2026-09-29"
OUT = ROOT / "results/2026-09-29-g6-chhatbar-metadata-audit"
REG = ROOT / "results/2026-09-24-disease-side-v3.5-g6-input-audit/registries"
CONFIG = ROOT / "configs/g6_chhatbar_source_audit_014.json"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def hash_file(path: Path, algorithm: str) -> str:
    digest = hashlib.new(algorithm)
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def normalized_donor(sample: str, suffixes: list[str]) -> str:
    for suffix in suffixes:
        if sample.endswith(suffix):
            return sample[: -len(suffix)]
    return sample


def norm_diagnosis(value: str) -> str:
    if value in {"AD", "AD_early", "AD_late", "AD late"}:
        return "AD"
    if value in {"Control", "Controls"}:
        return "Control"
    return value


def main() -> None:
    cfg = json.loads(CONFIG.read_text(encoding="utf-8"))
    OUT.mkdir(parents=True, exist_ok=True)
    record = json.loads((RAW / "zenodo_record_16938034.json").read_text(encoding="utf-8"))
    source_rows = []
    for item in record["files"]:
        path = RAW / item["key"]
        if not path.exists():
            continue
        md5 = hash_file(path, "md5")
        assert path.stat().st_size == item["size"], item["key"]
        assert md5 == item["checksum"].split(":", 1)[1], item["key"]
        source_rows.append({"file": item["key"], "bytes": path.stat().st_size,
                            "md5": md5, "sha256": hash_file(path, "sha256"),
                            "source": item["links"]["self"]})
    for name, url in [
        ("zenodo_record_16938034.json", "https://zenodo.org/api/records/16938034"),
        ("41590_2026_2472_MOESM4_ESM.xlsx",
         "https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41590-026-02472-z/MediaObjects/41590_2026_2472_MOESM4_ESM.xlsx"),
    ]:
        path = RAW / name
        source_rows.append({"file": name, "bytes": path.stat().st_size,
                            "md5": hash_file(path, "md5"), "sha256": hash_file(path, "sha256"),
                            "source": url})
    write_tsv(OUT / "input_checksums.tsv", source_rows,
              ["file", "bytes", "md5", "sha256", "source"])
    assert {"Nanostring_CosMX_v1.RDS", "Human_VisiumHD_compressed_v1.RDS"} <= {
        row["file"] for row in source_rows
    }

    workbook = load_workbook(RAW / "41590_2026_2472_MOESM4_ESM.xlsx", read_only=True,
                             data_only=True)
    sample_rows = list(workbook["Sup_table_3_Human_samples"].values)
    headers = sample_rows[0]
    supplement = [dict(zip(headers, row)) for row in sample_rows[1:] if row[0] is not None]
    by_donor_region = defaultdict(list)
    by_donor = defaultdict(list)
    for row in supplement:
        donor = str(row["Subject_ID"])
        by_donor_region[(donor, row["Region"])].append(row)
        by_donor[donor].append(row)
    panel = {str(row[0]).upper() for row in workbook["Sup_table_4_CosMX_gene_panel"].values
             if row[0] is not None}
    assert len(panel) == 1008

    with gzip.open(REG / "G1_microglia_signed_sender_244.tsv.gz", "rt", encoding="utf-8") as handle:
        g1 = {row["gene"].upper() for row in csv.DictReader(handle, delimiter="\t")}
    with gzip.open(REG / "G5_receiver_one_to_one_human_ortholog_registry.tsv.gz", "rt",
                   encoding="utf-8") as handle:
        g5_rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(g1) == 244
    modules: dict[str, set[str]] = defaultdict(set)
    status = {}
    for row in g5_rows:
        mod = row["module"]
        status[mod] = f"{row['G5_primary_status']}/{row['G5_module_gate']}"
        if row["one_to_one_mapping_status"].startswith("accepted") and row["human_one_to_one_symbol"]:
            modules[mod].add(row["human_one_to_one_symbol"].upper())
    signatures = {"G1_Microglia": g1, **modules}
    feature_sets = {"supplement_cosmx_panel": panel}
    for platform in ("cosmx", "visium_hd"):
        feature_sets[platform + "_object"] = {
            row["gene"].upper() for row in read_tsv(OUT / f"{platform}_features.tsv")
        }
    overlap_rows = []
    for source, features in feature_sets.items():
        for mod, genes in signatures.items():
            overlap = sorted(features & genes)
            overlap_rows.append({"source": source, "source_features": len(features),
                                 "module": mod, "upstream_status": status.get(mod, "frozen_G1"),
                                 "frozen_genes": len(genes), "overlap_n": len(overlap),
                                 "overlap_genes": ";".join(overlap)})
    write_tsv(OUT / "frozen_feature_overlap.tsv", overlap_rows,
              ["source", "source_features", "module", "upstream_status", "frozen_genes",
               "overlap_n", "overlap_genes"])
    discrepancy = [
        {"comparison": "panel_only", "gene": gene} for gene in sorted(panel - feature_sets["cosmx_object"])
    ] + [
        {"comparison": "object_only", "gene": gene} for gene in sorted(feature_sets["cosmx_object"] - panel)
    ]
    write_tsv(OUT / "cosmx_panel_object_discrepancy.tsv", discrepancy, ["comparison", "gene"])

    cosmx = read_tsv(OUT / "cosmx_sample_summary.tsv")
    visium = read_tsv(OUT / "visium_hd_sample_summary.tsv")
    assert sum(int(row["cell_count"]) for row in cosmx) == 68327
    assert sum(int(row["cell_count"]) for row in visium) == 96409
    assert len({row["Subject_ID"] for row in visium}) == len(visium)
    design = []
    for row in visium:
        donor = normalized_donor(row["Subject_ID"], cfg["donor_suffixes"])
        region = row["Tissue_Region"]
        exact = by_donor_region.get((donor, region), [])
        donor_rows = by_donor.get(donor, [])
        match = "exact_donor_region" if exact else ("donor_other_region_only" if donor_rows else "absent")
        supplement_row = exact[0] if exact else (donor_rows[0] if donor_rows else {})
        dx = norm_diagnosis(row["Clinical_Brain_Diagnosis_Short"])
        supplement_dx = norm_diagnosis(str(supplement_row.get("Clinical_Brain_Diagnosis_Short", "")))
        design.append({"platform": "visium_hd", "sample_id": row["Subject_ID"],
                       "derived_donor_id": donor, "diagnosis": dx,
                       "detailed_diagnosis": row["Clinical_Brain_Diagnosis"],
                       "region": row["Region"], "tissue_region": region,
                       "tma": row["orig.ident"], "cells": row["cell_count"],
                       "object_age": row["Age"], "object_sex": row["Reported_Sex"],
                       "supplement_match": match, "supplement_diagnosis": supplement_dx,
                       "supplement_age": str(supplement_row.get("Age", "")),
                       "supplement_sex": str(supplement_row.get("Reported_Sex", "")),
                       "supplement_matching_region_rows": len(exact)})
        if exact:
            assert supplement_dx == dx, row["Subject_ID"]
            assert str(supplement_row["Age"]) == row["Age"], row["Subject_ID"]
            assert supplement_row["Reported_Sex"] == row["Reported_Sex"], row["Subject_ID"]
    write_tsv(OUT / "visium_donor_region_design.tsv", design,
              ["platform", "sample_id", "derived_donor_id", "diagnosis", "detailed_diagnosis",
               "region", "tissue_region", "tma", "cells", "object_age", "object_sex",
               "supplement_match", "supplement_diagnosis", "supplement_age", "supplement_sex",
               "supplement_matching_region_rows"])
    region_groups = defaultdict(lambda: {"AD": set(), "Control": set(), "tmas": defaultdict(set)})
    for row in design:
        if row["diagnosis"] not in ("AD", "Control"):
            continue
        group = region_groups[row["region"]]
        group[row["diagnosis"]].add(row["derived_donor_id"])
        group["tmas"][row["tma"]].add(row["diagnosis"])
    region_rows = []
    for region, group in sorted(region_groups.items()):
        region_rows.append({"region": region, "ad_donors": len(group["AD"]),
                            "control_donors": len(group["Control"]),
                            "tmas_with_both_groups": sum(len(v) == 2 for v in group["tmas"].values()),
                            "screen": "SIZE_PASS_PENDING_FULL_AUDIT" if
                            len(group["AD"]) >= cfg["minimum_independent_ad_donors_per_region"] and
                            len(group["Control"]) >= cfg["minimum_independent_control_donors_per_region"]
                            else "AUXILIARY_ONLY_DONOR_N"})
    write_tsv(OUT / "visium_region_counts.tsv", region_rows,
              ["region", "ad_donors", "control_donors", "tmas_with_both_groups", "screen"])
    ba17 = [row for row in design if row["region"] == cfg["candidate_region"] and
            row["tissue_region"] == cfg["candidate_tissue_region"]]
    assert len(ba17) == 10
    assert len({row["derived_donor_id"] for row in ba17 if row["diagnosis"] == "AD"}) == 6
    assert len({row["derived_donor_id"] for row in ba17 if row["diagnosis"] == "Control"}) == 4
    detection = read_tsv(OUT / "visium_hd_ba17_detection.tsv")
    mg_zero = [row["sample_id"] for row in detection
               if row["module"] == "G1_Microglia" and int(row["target_cells"]) == 0]
    qc = [row for row in read_tsv(OUT / "visium_hd_ba17_donor_qc.tsv")
          if row["final_celltype"] == "All_nuclei"]
    detection_rate = read_tsv(OUT / "visium_hd_ba17_detection_rate.tsv")
    identity = read_tsv(OUT / "visium_hd_ba17_identity_marker_summary.tsv")
    fov = read_tsv(OUT / "visium_hd_fov_geometry_linkage.tsv")
    fov_cells = read_tsv(OUT / "visium_hd_fov_ba17_cell_linkage.tsv")
    duplicate_fov = read_tsv(OUT / "visium_hd_fov_duplicate_summary.tsv")
    assert len(qc) == 10 and len(identity) == 20 and len(fov) == 14
    assert len({row["cell_id"] for row in fov_cells}) == 23480
    qc_feature_fail = sum(int(row["n_below_author_5_feature_cutoff"]) for row in qc)
    qc_mito_fail = sum(int(row["n_above_author_75pct_mt_cutoff"]) for row in qc)
    g1_rates = [row for row in detection_rate if row["module"] == "G1_Microglia"]
    g1_at_1pct = [int(row["genes_at_or_above_1pct"]) for row in g1_rates]

    def donors_with_any_1pct(module: str) -> int:
        return len({row["sample_id"] for row in detection_rate
                    if row["module"] == module and int(row["genes_at_or_above_1pct"]) > 0})

    mg_identity = [row for row in identity if row["target_celltype"] == "Microglia"
                   and int(row["target_cells"]) > 0]
    asc_identity = [row for row in identity if row["target_celltype"] == "Astrocytes"]
    exact_ad = sum(row["supplement_match"] == "exact_donor_region" and row["diagnosis"] == "AD"
                   for row in ba17)
    exact_control = sum(row["supplement_match"] == "exact_donor_region" and row["diagnosis"] == "Control"
                        for row in ba17)
    duplicate_n = sum(int(row["duplicate_cell_ids"]) for row in duplicate_fov)
    ba17_fov = [row for row in fov if int(row["centroid_ids_matching_ba17_metadata"]) > 0]
    stage_tma_counts: dict[str, dict[str, int]] = defaultdict(lambda: defaultdict(int))
    for row in ba17:
        stage_tma_counts[row["detailed_diagnosis"]][row["tma"]] += 1
    decision = {
        "audit_id": cfg["audit_id"],
        "g6_status": "G6_NOT_EXECUTABLE",
        "visium_ba17_screen": "SIZE_PASS_PENDING_FULL_AUDIT",
        "visium_ba17_followup_status": "HOLD_DONOR_CROSSCHECK_IDENTITY_AND_SPATIAL_INPUTS",
        "visium_ba17_object_ad_donor_labels": 6,
        "visium_ba17_object_control_donor_labels": 4,
        "visium_ba17_supplement_exact_matches": sum(row["supplement_match"] == "exact_donor_region" for row in ba17),
        "visium_ba17_supplement_exact_ad_matches": exact_ad,
        "visium_ba17_supplement_exact_control_matches": exact_control,
        "visium_ba17_supplement_donor_other_region_only": sum(row["supplement_match"] == "donor_other_region_only" for row in ba17),
        "visium_ba17_supplement_absent": sum(row["supplement_match"] == "absent" for row in ba17),
        "visium_ba17_zero_microglia_samples": mg_zero,
        "visium_ba17_author_qc_cutoff_violations": {
            "nFeature_RNA_lt_5": qc_feature_fail,
            "percent_mt_gt_75": qc_mito_fail,
        },
        "visium_ba17_per_donor_median_nFeature_range": [
            min(int(float(row["nFeature_median"])) for row in qc),
            max(int(float(row["nFeature_median"])) for row in qc),
        ],
        "visium_ba17_G1_genes_at_1pct_range": [min(g1_at_1pct), max(g1_at_1pct)],
        "visium_ba17_state_detection_reference_donors": {
            "G1_any_gene_at_1pct": sum(n > 0 for n in g1_at_1pct),
            "G1_at_least_61_genes_at_1pct_historical_reference": sum(n >= 61 for n in g1_at_1pct),
            "M5_Iron_Uptake": donors_with_any_1pct("M5_Iron_Uptake"),
            "M5_Iron_Export": donors_with_any_1pct("M5_Iron_Export"),
            "M6_LPO_ACSL4": donors_with_any_1pct("M6_LPO_ACSL4"),
        },
        "visium_ba17_independent_marker_median_score_higher": {
            "Microglia": sum(row["target_marker_score_higher"] == "TRUE" for row in mg_identity),
            "Astrocytes": sum(row["target_marker_score_higher"] == "TRUE" for row in asc_identity),
        },
        "visium_ba17_unique_nuclei_linked_to_fov": len({row["cell_id"] for row in fov_cells}),
        "visium_ba17_fov_coordinate_rows_linked_including_representations": sum(int(row["centroid_ids_matching_ba17_metadata"]) for row in ba17_fov),
        "visium_ba17_duplicate_fov_cell_ids_same_coordinates": duplicate_n,
        "visium_ba17_detailed_diagnosis_by_tma": {diag: dict(counts) for diag, counts in stage_tma_counts.items()},
        "visium_ba17_coordinate_unit_confirmed": all(row["coordinate_unit_confirmed"] == "TRUE" for row in fov),
        "visium_ba17_tissue_mask_in_object": any(row["has_tissue_mask"] == "TRUE" for row in fov),
        "cosmx_parenchymal_ad_sample_labels": sum(row["Pathology"] == "AD" and "Meninges" not in row["Sample"] for row in cosmx),
        "cosmx_source_screen": "AUXILIARY_ONLY_DONOR_N_AND_PANEL_STRUCTURAL",
        "open_input_hold": "original_H_and_E_or_tissue_mask_not_in_zenodo_file_list;KaroSpace_route_not_validated_for_image_mask_or_scale;coordinate_unit_not_verified",
        "effect_analysis_performed": False,
    }
    (OUT / "decision.json").write_text(json.dumps(decision, indent=2), encoding="utf-8")
    print(json.dumps(decision, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
