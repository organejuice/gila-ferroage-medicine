"""Source-level SEA-AD G6 audit; never reads spatial expression or effects.

Run ``python src/103_audit_g6_seaad_metadata.py --fetch`` once to retrieve only
public small metadata/panel files, then omit ``--fetch`` for local reproduction.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import re
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from collections import Counter, defaultdict
from pathlib import Path

from openpyxl import load_workbook


ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data/raw/G6_sea_ad_metadata_2026-09-29"
OUT = ROOT / "results/2026-09-29-g6-sea-ad-metadata-audit"
REG = ROOT / "results/2026-09-24-disease-side-v3.5-g6-input-audit/registries"
CONFIG = ROOT / "configs/g6_seaad_source_audit_013.json"
BUCKET = "https://sea-ad-spatial-transcriptomics.s3.amazonaws.com"
NS = {"s3": "http://s3.amazonaws.com/doc/2006-03-01/"}
DONOR_URL = ("https://cdn.prod.website-files.com/689cfbd308fa7373b604d290/"
             "68debdfdd1b8e9f8fd64dab0_sea-ad_cohort_donor_metadata_072524.xlsx")
SUPP_URL = ("https://static-content.springer.com/esm/art%3A10.1038%2Fs41593-024-01774-5/"
            "MediaObjects/41593_2024_1774_MOESM3_ESM.xlsx")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def download(url: str, destination: Path) -> None:
    if destination.exists():
        return
    destination.parent.mkdir(parents=True, exist_ok=True)
    with urllib.request.urlopen(url, timeout=45) as response:
        if int(response.headers.get("Content-Length", "0")) > 2_000_000:
            raise ValueError(f"Refusing large metadata download: {url}")
        data = response.read(2_000_001)
    if len(data) > 2_000_000:
        raise ValueError(f"Refusing large metadata download: {url}")
    destination.write_bytes(data)


def list_prefix(prefix: str) -> tuple[list[str], list[dict]]:
    query = urllib.parse.urlencode({"list-type": "2", "delimiter": "/", "prefix": prefix})
    with urllib.request.urlopen(f"{BUCKET}/?{query}", timeout=45) as response:
        root = ET.fromstring(response.read())
    if root.findtext("s3:IsTruncated", namespaces=NS) != "false":
        raise ValueError(f"S3 listing is truncated: {prefix}")
    folders = [e.findtext("s3:Prefix", namespaces=NS)
               for e in root.findall("s3:CommonPrefixes", NS)]
    files = [{"key": e.findtext("s3:Key", namespaces=NS),
              "size": int(e.findtext("s3:Size", namespaces=NS))}
             for e in root.findall("s3:Contents", NS)]
    return folders, files


def fetch_small_sources() -> None:
    if (RAW / "s3_listings.json").exists() and (RAW / "cah_run_index.json").exists():
        return  # Immutable raw snapshot already captured; use a new directory for a new source version.
    download(DONOR_URL, RAW / "sea-ad_cohort_donor_metadata_072524.xlsx")
    download(SUPP_URL, RAW / "Gabitto_2024_supplementary_tables_1-8.xlsx")
    listings: dict[str, object] = {}
    region_dirs = ["middle-temporal-gyrus/", "caudate_nucleus/",
                   "hippocampus/", "medial-entorhinal-cortex/"]
    for region in region_dirs:
        folders, files = list_prefix(region)
        listings[region] = {"folders": folders, "files": files}
    runs: list[dict] = []
    for donor_prefix in listings["caudate_nucleus/"]["folders"]:
        donor = donor_prefix.rstrip("/").split("/")[-1]
        if not re.fullmatch(r"H\d{2}\.33\.\d{3}", donor):
            continue
        specimen_folders, _ = list_prefix(donor_prefix)
        for specimen_prefix in specimen_folders:
            specimen = specimen_prefix.rstrip("/").split("/")[-1]
            if not specimen.isdigit():
                continue
            output_parent = specimen_prefix + "xenium_output_and_alignment/"
            output_folders, _ = list_prefix(output_parent)
            if len(output_folders) != 1:
                raise ValueError(f"Expected one Xenium output for {specimen_prefix}: {output_folders}")
            output_prefix = output_folders[0]
            _, contents = list_prefix(output_prefix)
            by_name = {Path(item["key"]).name: item for item in contents}
            for name in ("gene_panel.json", "experiment.xenium"):
                if name not in by_name:
                    raise ValueError(f"Missing {name}: {output_prefix}")
                if by_name[name]["size"] > 2_000_000:
                    raise ValueError(f"Unexpectedly large {name}: {output_prefix}")
                download(f"{BUCKET}/{by_name[name]['key']}", RAW / "cah_runs" / donor / specimen / name)
            runs.append({"donor_id": donor, "specimen_id": specimen, "output_prefix": output_prefix,
                         "gene_panel_key": by_name["gene_panel.json"]["key"],
                         "experiment_key": by_name["experiment.xenium"]["key"],
                         "available_objects": {name: by_name[name]
                                               for name in ("cell_feature_matrix.h5", "cells.parquet",
                                                            "cell_boundaries.parquet", "nucleus_boundaries.parquet",
                                                            "transcripts.parquet", "morphology.ome.tif")
                                               if name in by_name}})
    (RAW / "s3_listings.json").write_text(json.dumps(listings, indent=2), encoding="utf-8")
    (RAW / "cah_run_index.json").write_text(json.dumps(runs, indent=2), encoding="utf-8")


def verify_access() -> None:
    runs = json.loads((RAW / "cah_run_index.json").read_text(encoding="utf-8"))
    representative = runs[0]
    keys = [
        "middle-temporal-gyrus/all_donors-h5ad/SEAAD_MTG_MERFISH.2024-12-11.h5ad",
        "caudate_nucleus/combined_anndata_object/CaH_Xenium.2026-01-07.h5ad",
    ] + [item["key"] for item in representative["available_objects"].values()]
    records = []
    for key in keys:
        url = f"{BUCKET}/{key}"
        request = urllib.request.Request(url, method="HEAD")
        try:
            with urllib.request.urlopen(request, timeout=45) as response:
                records.append({"key": key, "http_status": response.status,
                                "content_length": response.headers.get("Content-Length", ""), "url": url})
        except Exception as exc:
            records.append({"key": key, "http_status": f"ERROR: {type(exc).__name__}",
                            "content_length": "", "url": url})
    write_tsv(OUT / "access_head_checks.tsv", records,
              ["key", "http_status", "content_length", "url"])


def workbook_records(path: Path, sheet: str) -> list[dict]:
    rows = load_workbook(path, read_only=True, data_only=True)[sheet].values
    fields = next(rows)
    return [dict(zip(fields, values)) for values in rows if values[0]]


def signature_sets() -> tuple[set[str], dict[str, set[str]], dict[str, str]]:
    with gzip.open(REG / "G1_microglia_signed_sender_244.tsv.gz", "rt", encoding="utf-8") as handle:
        g1 = {row["gene"] for row in csv.DictReader(handle, delimiter="\t")}
    with gzip.open(REG / "G5_receiver_one_to_one_human_ortholog_registry.tsv.gz", "rt", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    g5: dict[str, set[str]] = defaultdict(set)
    statuses = {}
    for row in rows:
        statuses[row["module"]] = f"{row['G5_primary_status']}/{row['G5_module_gate']}"
        if row["one_to_one_mapping_status"].startswith("accepted") and row["human_one_to_one_symbol"]:
            g5[row["module"]].add(row["human_one_to_one_symbol"])
    assert len(g1) == 244
    return g1, dict(g5), statuses


def classify(row: dict, config: dict) -> dict:
    ad = row[config["clinical_ad_column"]] == config["clinical_positive_value"]
    control = row[config["clinical_control_column"]] == config["clinical_positive_value"]
    assert not (ad and control), row["Donor ID"]
    return {"clinical_ad": ad, "clinical_control": control,
            "adnc_high": row[config["adnc_column"]] == config["adnc_high_value"],
            "adnc_not_ad": row[config["adnc_column"]] == config["adnc_not_ad_value"]}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fetch", action="store_true", help="Download public small sources and S3 listings")
    parser.add_argument("--verify-access", action="store_true", help="HEAD representative large spatial inputs")
    args = parser.parse_args()
    if args.fetch:
        fetch_small_sources()
    OUT.mkdir(parents=True, exist_ok=True)
    if args.verify_access:
        verify_access()
    config = json.loads(CONFIG.read_text(encoding="utf-8"))
    supp = RAW / "Gabitto_2024_supplementary_tables_1-8.xlsx"
    donor_book = RAW / "sea-ad_cohort_donor_metadata_072524.xlsx"
    listings = json.loads((RAW / "s3_listings.json").read_text(encoding="utf-8"))
    runs = json.loads((RAW / "cah_run_index.json").read_text(encoding="utf-8"))
    st1 = {r["Donor ID"]: r for r in workbook_records(supp, "ST1")}
    donors = {r["Donor ID"]: r for r in workbook_records(donor_book, "SEA-AD_Cohort_Metadata")}
    assert set(st1) == set(donors) and len(donors) == 84
    mtg_flags = {donor for donor, row in st1.items()
                 if row[config["mtg_merfish_assayed_flag_column"]] == config["mtg_merfish_assayed_flag_value"]}
    mtg_folders = {folder.rstrip("/").split("/")[-1] for folder in listings["middle-temporal-gyrus/"]["folders"]
                   if re.fullmatch(r"H\d{2}\.33\.\d{3}", folder.rstrip("/").split("/")[-1])}
    # The bucket can hold raw sections outside the paper's assayed MERFISH set.
    assert mtg_flags <= mtg_folders, {"flag_only": sorted(mtg_flags - mtg_folders)}
    mtg_bucket_extra = sorted(mtg_folders - mtg_flags)
    cah_ids = {run["donor_id"] for run in runs}
    assert len(cah_ids) == 10 and cah_ids <= set(donors)
    assert len(runs) == len({(r["donor_id"], r["specimen_id"]) for r in runs})

    cohort_rows = []
    for region, ids in (("MTG_MERFISH", mtg_flags), ("CaH_XENIUM", cah_ids)):
        for donor in sorted(ids):
            row = donors[donor]
            flags = classify(row, config)
            cohort_rows.append({"region_platform": region, "donor_id": donor,
                                "cognitive_status": row["Cognitive Status"],
                                "adnc": row["Overall AD neuropathological Change"],
                                "consensus_clinical_ad": int(flags["clinical_ad"]),
                                "consensus_clinical_control": int(flags["clinical_control"]),
                                "adnc_high": int(flags["adnc_high"]),
                                "adnc_not_ad": int(flags["adnc_not_ad"]),
                                "n_spatial_rois": sum(r["donor_id"] == donor for r in runs) if region == "CaH_XENIUM" else ""})
    write_tsv(OUT / "donor_region_design.tsv", cohort_rows,
              ["region_platform", "donor_id", "cognitive_status", "adnc", "consensus_clinical_ad",
               "consensus_clinical_control", "adnc_high", "adnc_not_ad", "n_spatial_rois"])

    g1, g5, g5_status = signature_sets()
    merfish = {r["Gene"] for r in workbook_records(supp, "ST3")}
    assert len(merfish) == 140
    run_rows = []
    xenium_panels = []
    required_files = ("cell_feature_matrix.h5", "cells.parquet", "cell_boundaries.parquet",
                      "nucleus_boundaries.parquet", "transcripts.parquet", "morphology.ome.tif")
    for run in runs:
        local = RAW / "cah_runs" / run["donor_id"] / run["specimen_id"]
        panel = json.loads((local / "gene_panel.json").read_text(encoding="utf-8"))
        experiment = json.loads((local / "experiment.xenium").read_text(encoding="utf-8"))
        genes = {target["type"]["data"]["name"] for target in panel["payload"]["targets"]
                 if target["type"]["descriptor"] == "gene"}
        assert len(genes) == panel["payload"]["panel"]["num_gene_targets"]
        xenium_panels.append(genes)
        # Multiple regions can share a Xenium run/cassette; region_name is the ROI barcode.
        assert experiment["region_name"] == run["specimen_id"]
        run_rows.append({"donor_id": run["donor_id"], "specimen_id": run["specimen_id"],
                         "run_name": experiment.get("run_name", ""),
                         "cassette_name": experiment.get("cassette_name", ""),
                         "slide_id": experiment.get("slide_id", ""),
                         "instrument_sn": experiment.get("instrument_sn", ""),
                         "run_start_time": experiment.get("run_start_time", ""),
                         "analysis_sw_version": experiment.get("analysis_sw_version", ""),
                         "panel_name": experiment.get("panel_name", ""),
                         "panel_design_id": experiment.get("panel_design_id", ""),
                         "n_genes": len(genes), "n_cells_reported": experiment.get("num_cells", ""),
                         "all_six_spatial_objects_listed": int(all(k in run["available_objects"] for k in required_files)),
                         "panel_sha256": sha256(local / "gene_panel.json"),
                         "experiment_sha256": sha256(local / "experiment.xenium"),
                         "output_prefix": run["output_prefix"]})
    assert xenium_panels and all(genes == xenium_panels[0] for genes in xenium_panels)
    write_tsv(OUT / "cah_run_design.tsv", run_rows,
              ["donor_id", "specimen_id", "run_name", "cassette_name", "slide_id", "instrument_sn",
               "run_start_time", "analysis_sw_version", "panel_name", "panel_design_id", "n_genes",
               "n_cells_reported", "all_six_spatial_objects_listed", "panel_sha256",
               "experiment_sha256", "output_prefix"])
    panels = {"MTG_MERFISH_ST3": merfish, "CaH_XENIUM_all_rois": xenium_panels[0]}
    overlap = []
    for name, panel in panels.items():
        for signature, genes, status in [("G1_signed_microglia", g1, "FROZEN_G1")] + [
            (module, genes, g5_status[module]) for module, genes in sorted(g5.items())
        ]:
            overlap.append({"panel": name, "panel_size": len(panel), "signature": signature,
                            "upstream_status": status, "frozen_genes": len(genes),
                            "panel_genes": len(panel & genes),
                            "genes_on_panel": ";".join(sorted(panel & genes))})
    write_tsv(OUT / "frozen_panel_overlap.tsv", overlap,
              ["panel", "panel_size", "signature", "upstream_status", "frozen_genes",
               "panel_genes", "genes_on_panel"])

    input_files = [CONFIG, donor_book, supp, RAW / "s3_listings.json", RAW / "cah_run_index.json"] + [
        p for run in runs for p in ((RAW / "cah_runs" / run["donor_id"] / run["specimen_id"] / "gene_panel.json"),
                                  (RAW / "cah_runs" / run["donor_id"] / run["specimen_id"] / "experiment.xenium"))]
    write_tsv(OUT / "input_sha256.tsv",
              [{"file": str(p.relative_to(ROOT)).replace("\\", "/"), "bytes": p.stat().st_size,
                "sha256": sha256(p)} for p in input_files], ["file", "bytes", "sha256"])
    flags_by_donor = {donor: classify(row, config) for donor, row in donors.items()}
    summary = {"scope": "source_metadata_and_panel_only", "g6_status": "G6_NOT_EXECUTABLE",
               "source_level_min_independent_donors_per_group": config["source_level_min_independent_donors_per_group"],
               "mtg_merfish": {"donors": len(mtg_flags),
                               "consensus_clinical_ad": sum(flags_by_donor[d]["clinical_ad"] for d in mtg_flags),
                               "consensus_clinical_control": sum(flags_by_donor[d]["clinical_control"] for d in mtg_flags),
                               "adnc_high": sum(flags_by_donor[d]["adnc_high"] for d in mtg_flags),
                               "adnc_not_ad": sum(flags_by_donor[d]["adnc_not_ad"] for d in mtg_flags),
                               "panel_genes": len(merfish),
                               "bucket_donors_not_in_ST1_MERFISH2": mtg_bucket_extra},
               "cah_xenium": {"donors": len(cah_ids), "spatial_roi_outputs": len(runs),
                              "distinct_run_names": len({row["run_name"] for row in run_rows}),
                              "consensus_clinical_ad": sum(flags_by_donor[d]["clinical_ad"] for d in cah_ids),
                              "consensus_clinical_control": sum(flags_by_donor[d]["clinical_control"] for d in cah_ids),
                              "adnc_high": sum(flags_by_donor[d]["adnc_high"] for d in cah_ids),
                              "adnc_not_ad": sum(flags_by_donor[d]["adnc_not_ad"] for d in cah_ids),
                              "panel_genes": len(xenium_panels[0]),
                              "panels_identical_across_runs": True,
                              "rois_with_six_spatial_objects_listed": sum(row["all_six_spatial_objects_listed"] for row in run_rows),
                              "slide_by_clinical_group": {
                                  slide: dict(Counter("AD" if flags_by_donor[row["donor_id"]]["clinical_ad"]
                                                      else "Control" if flags_by_donor[row["donor_id"]]["clinical_control"]
                                                      else "Other" for row in run_rows if row["slide_id"] == slide))
                                  for slide in sorted({row["slide_id"] for row in run_rows})}},
               "other_spatial_regions": {region: listings[region]
                                         for region in ("hippocampus/", "medial-entorhinal-cortex/")},
               "large_expression_downloaded": False, "diagnosis_group_effects_computed": False}
    (OUT / "audit.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
