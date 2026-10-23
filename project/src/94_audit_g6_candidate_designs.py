#!/usr/bin/env python
"""Metadata-only design screen for G6-FP-2026-09-27-001 WP3.

This runner reads donor/sample metadata and platform/panel descriptions only.
It never opens an expression matrix or computes expression/diagnosis effects.
"""
from __future__ import annotations

import csv
import gzip
import hashlib
import json
import sys
from collections import Counter
from pathlib import Path

import numpy as np


def read_tsv(path: Path, selected: set[str] | None = None) -> tuple[list[str], list[dict[str, str]]]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames or []
        rows = []
        for row in reader:
            rows.append({k: row.get(k, "") for k in (selected or set(fields)) if k in row})
    return fields, rows


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def crosstab(rows: list[dict[str, str]], group: str, factor: str) -> list[dict[str, object]]:
    counts = Counter((r.get(group, ""), r.get(factor, "")) for r in rows)
    return [
        {"group": g, "factor": factor, "level": level, "n_donors": n}
        for (g, level), n in sorted(counts.items())
    ]


def design_rank(rows: list[dict[str, str]], group: str, factors: list[str]) -> tuple[int, int]:
    n = len(rows)
    cols = [np.ones(n, dtype=float)]
    for variable in [group, *factors]:
        levels = sorted({r.get(variable, "") for r in rows if r.get(variable, "") not in ("", "NA", "NaN")})
        for level in levels[1:]:
            cols.append(np.asarray([float(r.get(variable, "") == level) for r in rows]))
    x = np.column_stack(cols)
    return int(np.linalg.matrix_rank(x)), int(x.shape[1])


def factor_perfectly_aliases(rows: list[dict[str, str]], group: str, factor: str) -> bool:
    levels = {}
    for row in rows:
        level = row.get(factor, "")
        if not level:
            continue
        levels.setdefault(level, set()).add(row.get(group, ""))
    groups_by_level = list(levels.values())
    group_levels = {r.get(group, "") for r in rows}
    return len(group_levels) > 1 and len(levels) > 1 and all(len(x) == 1 for x in groups_by_level)


def main() -> int:
    if len(sys.argv) != 2:
        raise SystemExit("Usage: 94_audit_g6_candidate_designs.py <project_root>")
    root = Path(sys.argv[1])
    audit = root / "results/2026-09-24-disease-side-v3.5-g6-input-audit"
    out = root / "results/2026-09-27-g6-next-feasibility/candidates"
    out.mkdir(parents=True, exist_ok=False)
    source_rows: list[dict[str, object]] = []
    candidate_rows: list[dict[str, object]] = []
    design_rows: list[dict[str, object]] = []

    # GSE233208: only the author-corrected late sporadic AD/control donor set.
    gse233_path = audit / "registries/GSE233208_donor_section_manifest.tsv"
    _, all_233 = read_tsv(gse233_path)
    fields_233 = {"donor_id", "Diagnosis", "role_status", "CaptureArea", "Slide", "LibraryDay1", "LibraryDay2", "Sectioning", "RIN", "region", "platform", "section_sample_id"}
    rows_233 = [{k: row.get(k, "") for k in fields_233} for row in all_233
                if row.get("role_status") in {"validation_eligible_control", "validation_eligible_late_sporadic_AD"}]
    groups_233 = Counter(r["Diagnosis"] for r in rows_233)
    if len({r["donor_id"] for r in rows_233}) != len(rows_233):
        raise RuntimeError("GSE233208 selected donor IDs are not unique")
    tech_factors = ["Slide", "CaptureArea", "LibraryDay1", "LibraryDay2", "Sectioning"]
    aliases = [f for f in tech_factors if factor_perfectly_aliases(rows_233, "Diagnosis", f)]
    rank, ncol = design_rank(rows_233, "Diagnosis", tech_factors)
    for factor in tech_factors:
        design_rows.extend({"candidate": "GSE233208", **r} for r in crosstab(rows_233, "Diagnosis", factor))
    overlap_path = audit / "GSE174367_GSE233208_overlap_audit.json"
    overlap = json.loads(overlap_path.read_text(encoding="utf-8"))
    overlap_exact_crosswalk = False
    primary_gate = (
        groups_233.get("AD", 0) >= 10 and groups_233.get("Control", 0) >= 10
        and rank == ncol and not aliases and overlap_exact_crosswalk
    )
    candidate_rows.append({
        "candidate": "GSE233208", "status": "PRIMARY_CANDIDATE" if primary_gate else "NOT_ELIGIBLE_CONFIRMATORY",
        "n_AD": groups_233.get("AD", 0), "n_Control": groups_233.get("Control", 0),
        "region": "author metadata labels capture areas; anatomical-region common support not established",
        "platform": "10x Visium whole-transcriptome; 36,601 features in processed-object metadata audit",
        "design_rank": rank, "design_columns": ncol, "perfect_group_technical_alias": ";".join(aliases),
        "reference_independence": f"NOT_ESTABLISHED; {overlap['full_metadata_candidate_pairs']} demographic/pathology-matched candidates ({overlap['candidate_GSE233208_Diagnosis_counts']}); no exact donor crosswalk",
        "G1_G5_panel_coverage": "whole-transcriptome platform; frozen panel coverage not independently tabulated in this metadata-only run",
        "coordinates_mask_counts": "author-corrected metadata, images, coordinates and count-layer existence previously audited",
        "decision_reason": "A1 contains all selected controls and C1 all late sporadic AD; diagnosis and CaptureArea are perfectly aliased; exact reference donor independence unresolved",
    })
    for p in (gse233_path, overlap_path, audit / "tables/GSE233208_metadata_correction_audit.txt"):
        source_rows.append({"candidate": "GSE233208", "source_file": str(p.relative_to(root)), "source_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE233208", "sha256": sha256(p), "retrieved_or_audited": "2026-09-27"})

    # GSE220442: group counts are reported in the author table, but GEO-to-patient crosswalk is unresolved.
    patient_path = audit / "registries/GSE220442_patientID_donor_manifest.tsv"
    gsm_path = audit / "registries/GSE220442_GSM_source_sample_manifest.tsv"
    donor_path = audit / "registries/GSE220442_donor_metadata_manifest.tsv"
    _, patients = read_tsv(patient_path)
    _, gsms = read_tsv(gsm_path)
    _, donors_220 = read_tsv(donor_path)
    counts_220 = Counter(r.get("source_diagnosis", "") for r in patients)
    direct_matches = sum(r.get("mapping_status", "") == "exact_one_to_one_patientID_match" for r in gsms)
    mapped_220 = direct_matches == 6 and len({r.get("donor_id", "") for r in donors_220}) == 6
    candidate_rows.append({
        "candidate": "GSE220442", "status": "AUXILIARY_ONLY" if mapped_220 else "HOLD_METADATA_INCOMPLETE",
        "n_AD": counts_220.get("AD", 0), "n_Control": counts_220.get("control", 0),
        "region": "middle temporal gyrus", "platform": "10x Visium whole-transcriptome",
        "design_rank": "NOT_ESTIMABLE_FROM_GEO_TO_PATIENT_MAPPING", "design_columns": "",
        "perfect_group_technical_alias": "unknown: GEO sample-to-patient crosswalk unresolved",
        "reference_independence": "UNRESOLVED; no auditable crosswalk in this audit",
        "G1_G5_panel_coverage": "whole-transcriptome platform; frozen panel coverage not independently tabulated in this metadata-only run",
        "coordinates_mask_counts": "GEO archive lists spot matrices/images; exact GSM/patient identity link missing",
        "decision_reason": "3 AD + 3 control is below 10+10 primary threshold; currently cannot map any GSM one-to-one to patientID, so even auxiliary execution is held",
    })
    for p in (patient_path, gsm_path, donor_path, audit / "GSE220442_GSM_sample_mapping_audit.json"):
        source_rows.append({"candidate": "GSE220442", "source_file": str(p.relative_to(root)), "source_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE220442", "sha256": sha256(p), "retrieved_or_audited": "2026-09-27"})

    # Zenodo records: metadata/author descriptions only. No RDS or expression layer is opened.
    zenodo = [
        {
            "candidate": "Zenodo_16938034_CosMx", "status": "HOLD_METADATA_INCOMPLETE", "n_AD": "UNKNOWN", "n_Control": "UNKNOWN",
            "region": "multiple CNS disease contexts; exact per-donor CosMx map unavailable in local metadata audit",
            "platform": "Nanostring CosMx Human Universal Cell Characterization panel, 1,000-plex plus 30 GBM-only genes",
            "design_rank": "NOT_ESTIMABLE", "design_columns": "", "perfect_group_technical_alias": "UNKNOWN",
            "reference_independence": "UNRESOLVED", "G1_G5_panel_coverage": "NOT_COMPUTED: Supplementary Table 4 XLSX was not machine-readable in source inspection; 1,000-plex panel is targeted",
            "coordinates_mask_counts": "processed Seurat RDS listed by Zenodo; no matrix downloaded; per-donor coordinate/mask manifest unavailable",
            "decision_reason": "per-donor AD/control counts and exact panel-to-frozen-signature overlap not established; cannot rank as primary or auxiliary",
        },
        {
            "candidate": "Zenodo_16938034_VisiumHD", "status": "HOLD_METADATA_INCOMPLETE", "n_AD": "UNKNOWN", "n_Control": "UNKNOWN",
            "region": "multiple CNS regions/disease contexts; per-donor sample mapping not available locally",
            "platform": "Visium HD processed Seurat object; broad transcriptome",
            "design_rank": "NOT_ESTIMABLE", "design_columns": "", "perfect_group_technical_alias": "UNKNOWN",
            "reference_independence": "UNRESOLVED", "G1_G5_panel_coverage": "WHOLE_TRANSCRIPTOME_EXPECTED; not quantified from a source feature catalog",
            "coordinates_mask_counts": "processed Seurat RDS listed by Zenodo; no matrix downloaded; per-donor coordinate/mask manifest unavailable",
            "decision_reason": "per-donor AD/control counts, brain-region common support and technical design rank are not auditable from record metadata",
        },
    ]
    candidate_rows.extend(zenodo)
    zenodo_snapshot = [
        {"candidate": "Zenodo_16938034_CosMx", "record": "https://zenodo.org/records/16938034", "object": "Nanostring_CosMX_v1.RDS", "file_size": "162.3 MB", "platform_description": "CosMx Human Universal Cell Characterization panel; 1,000 genes; source article Supplementary Table 4 has human metadata and panel", "doi": "10.5281/zenodo.16938034", "source_article_doi": "10.1038/s41590-026-02472-z", "snapshot_date": "2026-09-27"},
        {"candidate": "Zenodo_16938034_VisiumHD", "record": "https://zenodo.org/records/16938034", "object": "Human_VisiumHD_compressed_v1.RDS", "file_size": "270.6 MB", "platform_description": "processed Visium HD object; record does not give per-donor AD/control crosswalk", "doi": "10.5281/zenodo.16938034", "source_article_doi": "10.1038/s41590-026-02472-z", "snapshot_date": "2026-09-27"},
    ]
    for row in zenodo_snapshot:
        source_rows.append({"candidate": row["candidate"], "source_file": "web metadata snapshot (see zenodo_source_snapshot.tsv)", "source_url": row["record"], "sha256": "", "retrieved_or_audited": "2026-09-27"})

    def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
        if not rows:
            path.write_text("", encoding="utf-8")
            return
        fields = list(dict.fromkeys(key for row in rows for key in row))
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
            writer.writeheader()
            writer.writerows(rows)

    write_tsv(out / "candidate_design_summary.tsv", candidate_rows)
    write_tsv(out / "technical_factor_crosstabs.tsv", design_rows)
    write_tsv(out / "source_hash_manifest.tsv", source_rows)
    write_tsv(out / "zenodo_source_snapshot.tsv", zenodo_snapshot)
    no_primary = all(r["status"] != "PRIMARY_CANDIDATE" for r in candidate_rows)
    report = [
        "# G6 alternative candidate screen (metadata only)",
        "",
        "Study: G6-FP-2026-09-27-001. No expression matrix was opened or downloaded, and no diagnosis-linked expression, state, or neighborhood effect was calculated.",
        "",
        f"GSE233208 selected late sporadic AD/control donors: AD={groups_233.get('AD', 0)}, Control={groups_233.get('Control', 0)}; design matrix rank={rank}/{ncol}; perfectly aliased factors={','.join(aliases) or 'none'}.",
        f"Reference independence remains unresolved: {overlap['full_metadata_candidate_pairs']} demographic/pathology-matched possible overlaps ({overlap['candidate_GSE233208_Diagnosis_counts']}); no exact shared donor ID crosswalk.",
        f"GSE220442 author donor metadata: AD={counts_220.get('AD', 0)}, Control={counts_220.get('control', 0)}; exact GEO GSM to patientID matches={direct_matches}/{len(gsms)}.",
        "CosMx and Visium HD processed objects are listed in Zenodo 16938034. The CosMx source article specifies a 1,000-plex targeted panel and links human sample metadata/panel in Supplementary Tables 3/4; the spreadsheet could not be machine-read in this metadata-only pass, so exact frozen G1/G5 coverage and donor design remain on hold.",
        "",
        f"No candidate qualifies as a primary candidate: {no_primary}.",
        "G6 remains G6_NOT_EXECUTABLE unless a future metadata correction resolves the listed gates. This screen does not authorize any group-effect model.",
    ]
    (out / "report.md").write_text("\n".join(report) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
