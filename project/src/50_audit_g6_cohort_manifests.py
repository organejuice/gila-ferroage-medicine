#!/usr/bin/env python3
"""Build effect-blind G6 donor, sample, and input-acquisition manifests."""

from __future__ import annotations

import gzip
import hashlib
import json
import re
import zipfile
import xml.etree.ElementTree as ET
from collections import Counter
from datetime import date
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RAW_META = ROOT / "data/raw/G6_source_metadata_2026-09-23"
OUT = ROOT / "results/2026-09-24-disease-side-v3.5-g6-input-audit"
REG = OUT / "registries"
TABLES = OUT / "tables"
for directory in (REG, TABLES):
    directory.mkdir(parents=True, exist_ok=True)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def write_tsv(rows: list[dict], path: Path) -> None:
    frame = pd.DataFrame(rows)
    frame.to_csv(path, sep="\t", index=False, compression="gzip" if path.suffix == ".gz" else None, na_rep="")


def parse_xlsx_sheet(xlsx_path: Path) -> list[dict[str, str]]:
    """Read the deposited workbook's shared-string worksheet without loss."""
    ns = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
    with zipfile.ZipFile(xlsx_path) as archive:
        strings_root = ET.fromstring(archive.read("xl/sharedStrings.xml"))
        strings = ["".join(item.itertext()) for item in strings_root.findall("m:si", ns)]
        sheet = ET.fromstring(archive.read("xl/worksheets/sheet1.xml"))

    def column_index(cell_ref: str) -> int:
        letters = re.match(r"[A-Z]+", cell_ref).group(0)
        value = 0
        for letter in letters:
            value = value * 26 + ord(letter) - ord("A") + 1
        return value - 1

    table: dict[int, list[str]] = {}
    for row in sheet.findall(".//m:sheetData/m:row", ns):
        values: list[str] = []
        for cell in row.findall("m:c", ns):
            value = cell.find("m:v", ns)
            text = "" if value is None else value.text or ""
            if cell.get("t") == "s" and text:
                text = strings[int(text)]
            col = column_index(cell.get("r", "A1"))
            while len(values) <= col:
                values.append("")
            values[col] = text
        table[int(row.get("r", "0"))] = values

    headers = table[2]
    records = []
    for row_number in sorted(table):
        if row_number < 3:
            continue
        values = table[row_number]
        row = {header: values[i] if i < len(values) else "" for i, header in enumerate(headers) if header}
        if row.get("Sample ID") and row.get("Diagnosis") in {"AD", "Control"}:
            records.append(row)
    return records


def parse_geo_filelist(path: Path) -> dict[str, int]:
    result: dict[str, int] = {}
    with path.open("rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) >= 4 and fields[0] == "File":
                try:
                    result[fields[1]] = int(fields[3])
                except ValueError:
                    continue
    return result


def parse_soft_samples(path: Path) -> list[dict[str, str]]:
    records: list[dict[str, str]] = []
    current: dict[str, str] | None = None
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.rstrip("\r\n")
            if line.startswith("^SAMPLE = "):
                if current is not None:
                    records.append(current)
                current = {"GSM": line.split("=", 1)[1].strip()}
            elif current is not None and line.startswith("!Sample_") and " = " in line:
                field, value = line.split(" = ", 1)
                key = field.removeprefix("!Sample_")
                value = value.strip().strip('"')
                current[key] = current.get(key, "") + (";" if key in current else "") + value
    if current is not None:
        records.append(current)
    return records


def classify_nia_aa_adnc(thal: str, braak: str, cerad: str) -> str:
    """Apply the published NIA-AA 2012 Thal/Braak/CERAD ADNC table."""
    try:
        thal_value, braak_value, cerad_value = int(float(thal)), int(float(braak)), int(float(cerad))
    except (TypeError, ValueError):
        return "not_classifiable_missing_or_noninteger_score"
    if thal_value not in range(6) or braak_value not in range(7) or cerad_value not in range(4):
        return "not_classifiable_out_of_range_score"
    a_score = 0 if thal_value == 0 else 1 if thal_value in (1, 2) else 2 if thal_value == 3 else 3
    b_score = 0 if braak_value == 0 else 1 if braak_value in (1, 2) else 2 if braak_value in (3, 4) else 3
    c_score = cerad_value
    if a_score == 0 and c_score == 0:
        return "Not_AD_neuropathologic_change"
    if a_score == 1 and c_score <= 1:
        return "Low_AD_neuropathologic_change"
    if a_score == 1 and c_score >= 2:
        return "Low_AD_neuropathologic_change" if b_score <= 1 else "Intermediate_AD_neuropathologic_change"
    if a_score == 2:
        return "Low_AD_neuropathologic_change" if b_score <= 1 else "Intermediate_AD_neuropathologic_change"
    if a_score == 3 and c_score <= 1:
        return "Low_AD_neuropathologic_change" if b_score <= 1 else "Intermediate_AD_neuropathologic_change"
    if a_score == 3 and c_score >= 2:
        return "Low_AD_neuropathologic_change" if b_score <= 1 else "Intermediate_AD_neuropathologic_change" if b_score == 2 else "High_AD_neuropathologic_change"
    return "not_classifiable"


# GSE269906: link each GEO section to the author's donor record and exact file sizes.
gse269_auth = parse_xlsx_sheet(RAW_META / "GSE269906_Supplementary_Data1.xlsx")
gse269_sizes = parse_geo_filelist(RAW_META / "GSE269906_filelist.txt")
gse269_soft = parse_soft_samples(RAW_META / "GSE269906_family.soft.gz")
soft_by_gsm = {record.get("GSM", ""): record for record in gse269_soft}
cellbin_dir = ROOT / "data/raw/GSE269906_cellbin"
tissue_dir_candidates = [ROOT / "data/raw/GSE269906_tissue", ROOT / "data/raw/GSE269906"]
tissue_dir = next((path for path in tissue_dir_candidates if path.exists()), None)
tissue_cn_manifest_path = REG / "GSE269906_tissue_CNGB_file_manifest.tsv"
tissue_cn_manifest = (
    pd.read_csv(tissue_cn_manifest_path, sep="\t", dtype={"GSM": str})
    if tissue_cn_manifest_path.exists() else pd.DataFrame()
)
tissue_cn_by_gsm = {
    row["GSM"]: row for row in tissue_cn_manifest.to_dict(orient="records")
}

gse269_rows: list[dict] = []
primary_files: list[dict] = []
for donor in gse269_auth:
    donor_id = donor["Sample ID"].strip()
    cell_name = next((name for name in gse269_sizes if name.endswith(".cellbin.gef") and donor_id in name), "")
    tissue_name = next((name for name in gse269_sizes if name.endswith(".tissue.gef") and donor_id in name), "")
    gsm_match = re.match(r"GSM\d+", cell_name)
    gsm = gsm_match.group(0) if gsm_match else ""
    soft = soft_by_gsm.get(gsm, {})
    expected_cell = gse269_sizes.get(cell_name, 0)
    cell_path = cellbin_dir / cell_name if cell_name else None
    actual_cell = cell_path.stat().st_size if cell_path and cell_path.exists() else 0
    cell_complete = bool(expected_cell and actual_cell == expected_cell)
    tissue_path = tissue_dir / tissue_name if tissue_dir and tissue_name else None
    actual_tissue = tissue_path.stat().st_size if tissue_path and tissue_path.exists() else 0
    expected_tissue = gse269_sizes.get(tissue_name, 0)
    tissue_complete = bool(expected_tissue and actual_tissue == expected_tissue)
    diagnosis = donor.get("Diagnosis", "")
    adnc = classify_nia_aa_adnc(donor.get("Thal", ""), donor.get("Braak", ""), donor.get("CERAD", ""))
    gse269_rows.append({
        "accession": "GSE269906",
        "role": "primary",
        "GSM": gsm,
        "donor_id": donor_id,
        "source_diagnosis": diagnosis,
        "source_subgroup": donor.get("AD levels", ""),
        "region": "BA10 prefrontal cortex",
        "platform": "Stereo-seq cell-bin",
        "sex": donor.get("Sex", ""),
        "age_band": donor.get("Age", ""),
        "Thal": donor.get("Thal", ""),
        "CERAD": donor.get("CERAD", ""),
        "Braak": donor.get("Braak", ""),
        "NIA_AA_2012_ADNC_level": adnc,
        "NIA_AA_2012_classification_rule": "Thal A0-A3, Braak B0-B3, and CERAD C0-C3 mapping per Montine et al. 2012 Table 3; descriptive pathology audit, not a post hoc diagnosis relabel",
        "race": donor.get("Race", ""),
        "sample_title": soft.get("title", ""),
        "cellbin_file": cell_name,
        "cellbin_expected_bytes": expected_cell,
        "cellbin_actual_bytes": actual_cell,
        "cellbin_status": "complete" if cell_complete else ("partial" if actual_cell else "missing"),
        "tissue_gef_file": tissue_name,
        "tissue_gef_expected_bytes": expected_tissue,
        "tissue_gef_actual_bytes": actual_tissue,
        "tissue_gef_status": "complete" if tissue_complete else ("partial" if actual_tissue else "not_downloaded"),
        "control_pathology_note": (
            "author-labelled Normal Aging; Braak II-III/CERAD 0-1 and Thal 0-4 are recorded; not pathology-free"
            if diagnosis.lower() == "control" else ""
        ),
    })
    for filename, expected, actual, status, local_path in (
        (cell_name, expected_cell, actual_cell, "complete" if cell_complete else ("partial" if actual_cell else "missing"), cell_path),
        (tissue_name, expected_tissue, actual_tissue, "complete" if tissue_complete else ("partial" if actual_tissue else "not_downloaded"), tissue_path),
    ):
        if filename:
            cn_tissue = tissue_cn_by_gsm.get(gsm, {}) if filename.endswith(".tissue.gef") else {}
            primary_files.append({
                "accession": "GSE269906",
                "GSM": gsm,
                "donor_id": donor_id,
                "file": filename,
                "file_type": "cellbin" if filename.endswith(".cellbin.gef") else "tissue",
                "source_url": f"https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM8330nnn/{gsm}/suppl/{filename}",
                "mirror_url": cn_tissue.get("mirror_url", ""),
                "mirror_source_database": cn_tissue.get("source_database", ""),
                "expected_bytes": expected,
                "actual_bytes": actual,
                "status": status,
                "sha256": (
                    cn_tissue.get("sha256", "")
                    if cn_tissue and status == "complete" and actual == int(cn_tissue.get("CNGB_actual_bytes", 0))
                    else sha256(local_path) if status == "complete" and local_path else ""
                ),
            })

if len(gse269_rows) != 12 or Counter(row["source_diagnosis"] for row in gse269_rows) != {"AD": 6, "Control": 6}:
    raise RuntimeError("GSE269906 author donor manifest failed the frozen 6/6 source-label check.")
if len({row["donor_id"] for row in gse269_rows}) != 12 or len({row["GSM"] for row in gse269_rows}) != 12:
    raise RuntimeError("GSE269906 donor/GSM mapping is not one-to-one.")
write_tsv(gse269_rows, REG / "GSE269906_donor_section_manifest.tsv")
write_tsv(primary_files, REG / "GSE269906_file_acquisition_manifest.tsv")

# GSE233208: consume the separately saved, author-corrected sample table.
gse233_path = TABLES / "GSE233208_author_corrected_sample_donor_manifest.tsv.gz"
if not gse233_path.exists():
    raise RuntimeError("Run src/49_prepare_g6_metadata_correction.R first.")
gse233 = pd.read_csv(gse233_path, sep="\t", compression="gzip", dtype={"Seq.ID": str})
gse233["role_status"] = gse233["Diagnosis"].map({
    "AD": "validation_eligible_late_sporadic_AD",
    "Control": "validation_eligible_control",
}).fillna("excluded_early_or_AD_DS")
gse233["accession"] = "GSE233208"
gse233["platform"] = "10x Visium"
gse233["region"] = gse233.get("CaptureArea", "")
gse233["donor_id"] = gse233["Case.Year"].astype(str) + "_" + gse233["Case.Num"].astype(str)
gse233["section_sample_id"] = gse233["Seq.ID"]
write_tsv(gse233.to_dict(orient="records"), REG / "GSE233208_donor_section_manifest.tsv")
late = gse233[gse233["Diagnosis"].isin(["AD", "Control"])]
if late["Diagnosis"].value_counts().to_dict() != {"AD": 10, "Control": 10}:
    raise RuntimeError("GSE233208 corrected late sporadic AD/control donor counts are not 10/10.")
if late["donor_id"].nunique() != 20 or late["Seq.ID"].nunique() != 20:
    raise RuntimeError("GSE233208 eligible late-AD/control donors are not unique.")
overlap_audit_path = OUT / "GSE174367_GSE233208_overlap_audit.json"
overlap_audit = json.loads(overlap_audit_path.read_text(encoding="utf-8")) if overlap_audit_path.exists() else {}
overlap_counts = overlap_audit.get("candidate_GSE233208_Diagnosis_counts", {})
candidate_overlap_ad = int(overlap_counts.get("AD", 0))
candidate_overlap_control = int(overlap_counts.get("Control", 0))

# GSE220442: only spot metadata is staged; keep the donor key and label audit separate.
gse220_metadata = pd.read_csv(RAW_META / "GSE220442_metadata.csv.gz", compression="gzip")
gse220_rows = []
for donor_id, group in gse220_metadata.groupby("patientID", dropna=False, sort=True):
    categories = sorted(str(value) for value in group["category"].dropna().unique())
    stages = sorted(str(value) for value in group["stage"].dropna().unique())
    gse220_rows.append({
        "accession": "GSE220442",
        "role": "cross_region_supportive_only",
        "GSM": "",
        "donor_id": donor_id,
        "source_diagnosis": ";".join(categories),
        "source_stage": ";".join(stages),
        "region": "middle temporal gyrus",
        "platform": "10x Visium",
        "n_spots": len(group),
        "layers_present": ";".join(sorted(str(value) for value in group["Layer"].dropna().unique())),
        "GSM_to_patientID_mapping": "not established from staged GEO metadata",
        "expression_matrix_and_image_available_locally": False,
    })
write_tsv(gse220_rows, REG / "GSE220442_donor_metadata_manifest.tsv")
if len(gse220_rows) != 6 or Counter(row["source_diagnosis"] for row in gse220_rows) != {"control": 3, "AD": 3}:
    raise RuntimeError("GSE220442 staged metadata did not resolve to the reported 3/3 donor groups.")

# Hash source-level manifests and small metadata already staged.
source_specs = [
    (RAW_META / "GSE269906_family.soft.gz", "GSE269906", "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE269nnn/GSE269906/soft/GSE269906_family.soft.gz", "GEO source family metadata"),
    (RAW_META / "GSE269906_filelist.txt", "GSE269906", "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE269nnn/GSE269906/", "GEO series file-size listing"),
    (RAW_META / "GSE269906_Supplementary_Data1.xlsx", "GSE269906", "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE269nnn/GSE269906/suppl/GSE269906_Supplementary_Data1.xlsx", "author donor neuropathology table"),
    (RAW_META / "GSE233208_family.soft.gz", "GSE233208", "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE233nnn/GSE233208/soft/GSE233208_family.soft.gz", "GEO source family metadata"),
    (RAW_META / "GSE233208_SupplementaryTables.xlsx", "GSE233208", "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE233nnn/GSE233208/suppl/GSE233208_SupplementaryTables.xlsx", "complete publisher supplementary workbook"),
    (RAW_META / "GSE233208_author_Supplementary_Table1.tsv.gz", "GSE233208", "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE233nnn/GSE233208/suppl/", "author Supplementary Table 1 extracted without modification from publisher workbook"),
    (ROOT / "data/raw/GSE233208/GSE233208_Human_visium_ADDS_seurat_processed.rds.gz", "GSE233208", "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE233nnn/GSE233208/suppl/GSE233208_Human_visium_ADDS_seurat_processed.rds.gz", "immutable publisher-deposited compressed Seurat object"),
    (ROOT / "data/interim/GSE233208/GSE233208_Human_visium_ADDS_seurat_processed.rds", "GSE233208", "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE233nnn/GSE233208/suppl/GSE233208_Human_visium_ADDS_seurat_processed.rds.gz", "one outer gzip layer removed for readRDS; derived interim object"),
    (RAW_META / "GSE220442_family.soft.gz", "GSE220442", "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE220nnn/GSE220442/soft/GSE220442_family.soft.gz", "GEO source family metadata"),
    (RAW_META / "GSE220442_metadata.csv.gz", "GSE220442", "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE220nnn/GSE220442/suppl/GSE220442_metadata.csv.gz", "GEO spot metadata"),
]
source_rows = []
for path, accession, source_url, processing_level in source_specs:
    source_rows.append({
        "accession": accession,
        "database": "NCBI GEO",
        "source_url": source_url,
        "download_date": date.today().isoformat(),
        "file": path.relative_to(ROOT).as_posix(),
        "bytes": path.stat().st_size,
        "sha256": sha256(path),
        "processing_level": processing_level,
        "license_access_conditions": "Public NCBI GEO; no separate file-level license was recorded in the staged series metadata",
    })
write_tsv(source_rows, REG / "G6_source_metadata_manifest.tsv")

primary_cell_complete = sum(row["cellbin_status"] == "complete" for row in gse269_rows)
primary_tissue_complete = sum(row["tissue_gef_status"] == "complete" for row in gse269_rows)
validation_images = TABLES / "GSE233208_spatial_image_coordinate_manifest.tsv.gz"
validation_images_df = pd.read_csv(validation_images, sep="\t", compression="gzip")
validation_counts = TABLES / "GSE233208_expression_layer_audit.tsv.gz"
validation_counts_df = pd.read_csv(validation_counts, sep="\t", compression="gzip")
audit = {
    "run_id": "2026-09-24-disease-side-v3.5-g6-input-audit",
    "protocol": "docs/protocol/disease_side_v3.5_latest_frozen.md",
    "protocol_status": "G6_ANALYSIS_LOCK_PENDING_INPUT_AUDIT",
    "date": date.today().isoformat(),
    "primary": {
        "accession": "GSE269906",
        "donors": len(gse269_rows),
        "AD_donors": sum(row["source_diagnosis"] == "AD" for row in gse269_rows),
        "Control_donors": sum(row["source_diagnosis"] == "Control" for row in gse269_rows),
        "complete_cellbin_files": primary_cell_complete,
        "total_cellbin_files": len(gse269_rows),
        "complete_tissue_gef_files": primary_tissue_complete,
        "total_tissue_gef_files": len(gse269_rows),
        "NIA_AA_2012_ADNC_levels_in_author_labelled_controls": pd.DataFrame(gse269_rows).query("source_diagnosis == 'Control'")["NIA_AA_2012_ADNC_level"].value_counts().to_dict(),
        "author_control_label_retained": "Normal Aging; neuropathology is described per donor and is not silently called pathology-free",
        "input_admissible": False,
        "reason": "Raw CellBin and tissue GEF files are audited for donor mapping and integrity; the effect-blind bin50, segmentation-mask, spatial-support, and frozen-state coverage audit must finish before the primary input can be released.",
    },
    "independent_validation": {
        "accession": "GSE233208",
        "author_correction_reproduced": True,
        "unique_samples": int(gse233["Seq.ID"].nunique()),
        "unique_author_donors": int(gse233["donor_id"].nunique()),
        "late_sporadic_AD_donors": int((gse233["Diagnosis"] == "AD").sum()),
        "control_donors": int((gse233["Diagnosis"] == "Control").sum()),
        "earlyAD_excluded_donors": int((gse233["Diagnosis"] == "earlyAD").sum()),
        "AD_DS_excluded_donors": int((gse233["Diagnosis"] == "AD_DS").sum()),
        "all_spot_metadata_and_image_coordinates_mapped": bool(
            validation_images_df["n_coordinates_unmatched_to_metadata"].sum() == 0
            and validation_images_df["n_coordinates_matched_to_metadata"].sum() == 115451
        ),
        "integer_count_layer": int(validation_counts_df.loc[0, "noninteger_nonzero_entries"]) == 0,
        "cell_type_inference_provenance": "processed object includes MG_deconv and ASC_deconv fields but the author's public workflow does not document their method/reference; an independent local human GSE174367 reference matrix and donor manifest have been audited",
        "GSE174367_candidate_overlap_donors": overlap_audit.get("candidate_GSE233208_donors", []),
        "candidate_overlap_AD_donors": candidate_overlap_ad,
        "candidate_overlap_Control_donors": candidate_overlap_control,
        "eligible_independent_AD_donors_after_fail_closed_exclusion": max(0, 10 - candidate_overlap_ad),
        "eligible_independent_Control_donors_after_fail_closed_exclusion": max(0, 10 - candidate_overlap_control),
        "independent_validation_eligible": bool(candidate_overlap_ad == 0 and candidate_overlap_control == 0),
        "independence_audit_rule": "Public age, sex, PMI, diagnosis, tangle-stage, and plaque-stage concordance is treated as candidate donor overlap; candidate donors are not counted as independent absent an author crosswalk.",
    },
    "supportive": {
        "accession": "GSE220442",
        "donors": len(gse220_rows),
        "AD_donors": sum(row["source_diagnosis"] == "AD" for row in gse220_rows),
        "control_donors": sum(row["source_diagnosis"] == "control" for row in gse220_rows),
        "available_input": "spot metadata only",
        "GSM_to_patientID_mapping_resolved": False,
    },
    "execution_lock_issued": False,
    "g6_group_effects_calculated": False,
    "decision": "G6_NOT_EXECUTABLE",
    "blockers": [],
}
if primary_cell_complete != 12:
    audit["blockers"].append("GSE269906 primary cellbin files are not complete for all 12 donors.")
audit["blockers"].append("Primary per-section QC, segmentation-mask overlays, physical coordinate scale, and state-gene coverage have not yet been released by the effect-blind audit.")
if not (REG / "G5_receiver_one_to_one_human_ortholog_registry.tsv.gz").exists():
    audit["blockers"].append("One-to-one human ortholog mapping for the G5 receiver family is not available.")
if overlap_audit and not audit["independent_validation"]["independent_validation_eligible"]:
    audit["blockers"].append("GSE233208 does not provide an independent AD validation set after fail-closed exclusion of candidate donors overlapping the G1 GSE174367 reference cohort.")
if not overlap_audit:
    audit["blockers"].append("Cross-cohort donor-identity audit against the G1 reference cohort has not yet run.")
audit["primary"]["tissue_gef_files_required"] = True
audit["primary"]["mask_source_pending_direct_audit"] = "CellBin/cellBorder segmentation polygons are available; polygon-derived bin50 tissue masks and their agreement with tissue-GEF occupied bins remain under effect-blind audit."
(OUT / "G6_input_audit.json").write_text(json.dumps(audit, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
print(json.dumps({
    "primary_cellbin_complete": primary_cell_complete,
    "primary_tissue_complete": primary_tissue_complete,
    "GSE233208_samples": len(gse233),
    "GSE220442_donors": len(gse220_rows),
    "decision": audit["decision"],
}, ensure_ascii=False))
