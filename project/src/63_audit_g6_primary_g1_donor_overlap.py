#!/usr/bin/env python3
"""Effect-blind donor candidate audit: GSE269906 primary versus GSE174367 G1."""

from __future__ import annotations

import json
import re
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/2026-09-24-disease-side-v3.5-g6-input-audit"
TABLES = OUT / "tables"
REG = OUT / "registries"

primary = pd.read_csv(REG / "GSE269906_donor_section_manifest.tsv", sep="\t", dtype=str).fillna("")
d8 = pd.read_csv(ROOT / "data/raw/D8_GSE174367/GSE174367_snRNA-seq_cell_meta.csv.gz", dtype=str).fillna("")
d8 = (
    d8.groupby("SampleID", as_index=False, sort=True)
    .agg({
        "Diagnosis": "first", "Age": "first", "Sex": "first", "PMI": "first",
        "Tangle.Stage": "first", "Plaque.Stage": "first", "Barcode": "size",
    })
    .rename(columns={
        "SampleID": "GSE174367_donor_id", "Diagnosis": "D8_diagnosis",
        "Age": "D8_age", "Sex": "D8_sex", "PMI": "D8_PMI",
        "Tangle.Stage": "D8_Tangle.Stage", "Plaque.Stage": "D8_Plaque.Stage",
        "Barcode": "D8_labelled_nuclei",
    })
)

def stage_number(value: str) -> str:
    found = re.search(r"\d+", value)
    return found.group(0) if found else ""

def age_interval(value: str) -> tuple[int, int] | None:
    found = re.search(r"(\d+)\s*[-–]\s*(\d+)", value)
    return (int(found.group(1)), int(found.group(2))) if found else None

def cerad_label(value: str) -> str:
    try:
        code = int(float(value))
    except ValueError:
        return ""
    return {0: "None", 1: "Stage A", 2: "Stage B", 3: "Stage C"}.get(code, "")

rows: list[dict[str, str | bool]] = []
for _, p in primary.iterrows():
    interval = age_interval(p["age_band"])
    for _, d in d8.iterrows():
        try:
            age = int(float(d["D8_age"]))
        except ValueError:
            continue
        diagnosis_sex_age_candidate = (
            p["source_diagnosis"].casefold() == d["D8_diagnosis"].casefold()
            and p["sex"].casefold() == d["D8_sex"].casefold()
            and interval is not None and interval[0] <= age <= interval[1]
        )
        if not diagnosis_sex_age_candidate:
            continue
        braak_match = stage_number(p["Braak"]) == stage_number(d["D8_Tangle.Stage"])
        cerad_match = cerad_label(p["CERAD"]).casefold() == d["D8_Plaque.Stage"].casefold()
        rows.append({
            "GSE269906_donor_id": p["donor_id"],
            "GSE269906_GSM": p["GSM"],
            "GSE269906_diagnosis": p["source_diagnosis"],
            "GSE269906_age_band": p["age_band"],
            "GSE269906_sex": p["sex"],
            "GSE269906_Braak": p["Braak"],
            "GSE269906_CERAD": p["CERAD"],
            "GSE174367_donor_id": d["GSE174367_donor_id"],
            "GSE174367_diagnosis": d["D8_diagnosis"],
            "GSE174367_age": d["D8_age"],
            "GSE174367_sex": d["D8_sex"],
            "GSE174367_PMI": d["D8_PMI"],
            "GSE174367_Tangle_Stage": d["D8_Tangle.Stage"],
            "GSE174367_Plaque_Stage": d["D8_Plaque.Stage"],
            "diagnosis_sex_age_band_candidate": True,
            "Braak_stage_number_match": braak_match,
            "CERAD_plaque_stage_match_under_explicit_crosswalk": cerad_match,
            "both_available_pathology_fields_match": braak_match and cerad_match,
        })

matches = pd.DataFrame(rows)
if matches.empty:
    matches = pd.DataFrame(columns=[
        "GSE269906_donor_id", "GSE269906_GSM", "GSE269906_diagnosis", "GSE269906_age_band",
        "GSE269906_sex", "GSE269906_Braak", "GSE269906_CERAD", "GSE174367_donor_id",
        "GSE174367_diagnosis", "GSE174367_age", "GSE174367_sex", "GSE174367_PMI",
        "GSE174367_Tangle_Stage", "GSE174367_Plaque_Stage", "diagnosis_sex_age_band_candidate",
        "Braak_stage_number_match", "CERAD_plaque_stage_match_under_explicit_crosswalk",
        "both_available_pathology_fields_match",
    ])
matches.to_csv(TABLES / "GSE269906_GSE174367_candidate_donor_overlap.tsv", sep="\t", index=False)

summary = {
    "primary_accession": "GSE269906",
    "G1_reference_accession": "GSE174367",
    "primary_donors": int(primary["donor_id"].nunique()),
    "G1_reference_donors": int(d8["GSE174367_donor_id"].nunique()),
    "demographic_diagnosis_candidate_pairs": int(len(matches)),
    "GSE269906_candidate_donors": sorted(matches["GSE269906_donor_id"].unique().tolist()) if len(matches) else [],
    "GSE174367_candidate_donors": sorted(matches["GSE174367_donor_id"].unique().tolist()) if len(matches) else [],
    "pairs_also_matching_braak_and_plaque_stage": int(matches["both_available_pathology_fields_match"].sum()) if len(matches) else 0,
    "candidate_identity_basis": "GSE269906 reports five-year age bands; candidates match diagnosis, sex, and a GSE174367 exact age within that band. Available Braak and CERAD versus Tangle/Plaque stage fields are reported as corroboration only, not as proof of identity. No shared donor crosswalk or GSE269906 PMI was found in staged source metadata.",
    "donor_identity_crosswalk_found": False,
    "effects_used": False,
}
(OUT / "GSE269906_GSE174367_overlap_audit.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
print(json.dumps(summary, ensure_ascii=False))
