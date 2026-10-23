#!/usr/bin/env python3
"""Audit GSE233208 versus the GSE174367 G1 reference using public donor fields."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/2026-09-24-disease-side-v3.5-g6-input-audit"
TABLES = OUT / "tables"
D8_META = ROOT / "data/raw/D8_GSE174367/GSE174367_snRNA-seq_cell_meta.csv.gz"
GSE233 = TABLES / "GSE233208_author_corrected_sample_donor_manifest.tsv.gz"

d8 = pd.read_csv(D8_META, dtype=str).fillna("")
d8_donors = (
    d8.groupby("SampleID", as_index=False, sort=True)
    .agg({
        "Diagnosis": "first",
        "Age": "first",
        "Sex": "first",
        "PMI": "first",
        "Tangle.Stage": "first",
        "Plaque.Stage": "first",
        "Barcode": "size",
    })
    .rename(columns={
        "SampleID": "D8_donor_id",
        "Barcode": "D8_labelled_nuclei",
        "Tangle.Stage": "D8_Tangle_stage",
        "Plaque.Stage": "D8_Plaque_stage",
    })
)
v = pd.read_csv(GSE233, sep="\t", dtype=str).fillna("")
v["GSE233_author_donor_id"] = v["Case.Year"].astype(str) + "_" + v["Case.Num"].astype(str)
v["role_status"] = v["Diagnosis"].map({
    "AD": "validation_eligible_late_sporadic_AD",
    "Control": "validation_eligible_control",
}).fillna("excluded_early_or_AD_DS")

shared = []
for _, d in d8_donors.iterrows():
    for _, g in v.iterrows():
        diagnosis_match = d["Diagnosis"] == g["Diagnosis"]
        age_match = d["Age"].strip() == g["Age"].strip()
        sex_match = d["Sex"].strip().upper() == g["Sex"].strip().upper()
        pmi_match = d["PMI"].strip() == g["PMI"].strip()
        tangle_match = d["D8_Tangle_stage"].strip() == g["Tangle.Stage"].strip()
        plaque_match = d["D8_Plaque_stage"].strip() == g["Plaque.Stage"].strip()
        # Require the full independent-field tuple to avoid identifying people
        # from age/sex alone; exact repeated PMI and both pathology stages add
        # strong source-level corroboration.
        exact = diagnosis_match and age_match and sex_match and pmi_match and tangle_match and plaque_match
        broad = diagnosis_match and age_match and sex_match and pmi_match
        if exact or broad:
            shared.append({
                "D8_donor_id": d["D8_donor_id"],
                "GSE233208_donor_id": g["GSE233_author_donor_id"],
                "GSE233208_Seq.ID": g["Seq.ID"],
                "diagnosis": d["Diagnosis"],
                "Age": d["Age"],
                "Sex": d["Sex"],
                "PMI_D8": d["PMI"],
                "PMI_GSE233208": g["PMI"],
                "Tangle_stage_D8": d["D8_Tangle_stage"],
                "Tangle_stage_GSE233208": g["Tangle.Stage"],
                "Plaque_stage_D8": d["D8_Plaque_stage"],
                "Plaque_stage_GSE233208": g["Plaque.Stage"],
                "age_sex_PMI_diagnosis_match": broad,
                "pathology_stages_match": tangle_match and plaque_match,
                "full_metadata_match": exact,
                "GSE233208_validation_role_status": g["role_status"],
            })

matches = pd.DataFrame(shared)
if matches.empty:
    matches = pd.DataFrame(columns=[
        "D8_donor_id", "GSE233208_donor_id", "GSE233208_Seq.ID", "diagnosis",
        "Age", "Sex", "PMI_D8", "PMI_GSE233208", "Tangle_stage_D8",
        "Tangle_stage_GSE233208", "Plaque_stage_D8", "Plaque_stage_GSE233208",
        "age_sex_PMI_diagnosis_match", "pathology_stages_match", "full_metadata_match",
        "GSE233208_validation_role_status",
    ])
matches.to_csv(TABLES / "GSE174367_GSE233208_candidate_donor_overlap.tsv", sep="\t", index=False)

summary = {
    "GSE174367_G1_reference_donors": int(d8_donors["D8_donor_id"].nunique()),
    "GSE233208_author_donors": int(v["GSE233_author_donor_id"].nunique()),
    "candidate_pairs_matching_age_sex_PMI_diagnosis": int(matches["age_sex_PMI_diagnosis_match"].sum()) if not matches.empty else 0,
    "candidate_pairs_matching_both_pathology_stages": int(matches["pathology_stages_match"].sum()) if not matches.empty else 0,
    "full_metadata_candidate_pairs": int(matches["full_metadata_match"].sum()) if not matches.empty else 0,
    "candidate_GSE233208_donors": sorted(matches["GSE233208_donor_id"].unique().tolist()) if not matches.empty else [],
    "candidate_GSE233208_Diagnosis_counts": matches.drop_duplicates("GSE233208_donor_id")["diagnosis"].value_counts().to_dict() if not matches.empty else {},
    "identity_key_source": "No shared donor ID crosswalk was found; candidate identity is based on exact author-provided age, sex, PMI, diagnosis, and pathology-stage concordance.",
    "G6_effects_used": False,
}
(OUT / "GSE174367_GSE233208_overlap_audit.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
print(json.dumps(summary, ensure_ascii=False))
