#!/usr/bin/env python3
"""Invariant audit for the disease-side completion package."""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    tables = root / "results" / "tables"
    interim = root / "data" / "interim"
    checks = {}

    d1 = pd.read_csv(tables / "D1_annotation_audit.tsv", sep="\t")
    checks["D1_peak_to_gene_rate"] = float(d1["mapping_rate"].min()) >= 0.90
    d2 = pd.read_csv(tables / "D2_transcript_gene_audit.tsv", sep="\t")
    checks["D2_transcript_mapping_rate"] = float(d2["transcript_mapping_rate"].iloc[0]) >= 0.80
    d3 = pd.read_csv(tables / "D3_module_effects_all.tsv", sep="\t")
    checks["D3_all_module_rows"] = len(d3) == 18
    checks["D3_no_module_fdr_hit"] = int((d3["FDR"] < 0.05).sum()) == 0
    d4 = pd.read_csv(tables / "D4_spatial_support_summary.tsv", sep="\t")
    checks["D4_two_samples"] = len(d4) == 2 and int(d4["n_cells"].sum()) > 40000
    i3 = pd.read_csv(tables / "I3_disease_module_vector.tsv", sep="\t")
    d6_h5 = root / "data" / "raw" / "D6_GSE214979" / "GSE214979_filtered_feature_bc_matrix.h5"
    d6_module = pd.read_csv(tables / "D6_module_effects.tsv", sep="\t")
    d6_concordance = pd.read_csv(tables / "D6_signature_concordance.tsv", sep="\t")
    checks["I3_modules_present"] = set(i3["module"]) >= {"M1", "M2", "M3", "M4", "M5_Iron_Uptake", "M5_Iron_Storage", "M5_Iron_Export", "M5_Ferritinophagy", "M6", "M7"}
    checks["I3_includes_D6_layer"] = "D6_human_independent_validation" in set(i3["layer"])
    checks["D6_h5_complete"] = d6_h5.stat().st_size == 1369492123 and d6_h5.open("rb").read(8) == bytes([137, 72, 68, 70, 13, 10, 26, 10])
    checks["D6_validation_outputs"] = len(d6_module) == 12 and len(d6_concordance) == 2 and d6_module["FDR"].notna().all()
    c1 = pd.read_csv(tables / "C1_communication_gate.tsv", sep="\t")
    checks["C1_no_go"] = c1.loc[0, "decision"] == "NO_GO" and int(c1.loc[0, "observed_D3_gene_FDR_hits"]) == 0 and int(c1.loc[0, "observed_D3_module_FDR_hits"]) == 0
    checks["C1_d6_fields"] = {"observed_D6_module_FDR_hits", "D6_microglia_signature_spearman", "D6_astrocyte_signature_spearman"}.issubset(set(c1.columns))
    assessment = pd.read_csv(tables / "MechanismEvidenceAssessment.tsv", sep="\t")
    checks["mechanism_level_1_of_5"] = assessment.loc[0, "level"] == "1/5"
    completion = pd.read_csv(tables / "DiseaseSideCompletion.tsv", sep="\t")
    checks["completion_has_D1_to_D6"] = set(completion["stage"]) >= {"D1", "D2", "D3", "D4", "D5", "D6", "I1", "I2", "I3", "C1"}
    checks = {key: bool(value) for key, value in checks.items()}
    passed = all(checks.values())
    audit = {"status": "PASS" if passed else "FAIL", "checks": checks}
    (interim / "disease_side_audit.json").write_text(json.dumps(audit, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(audit, ensure_ascii=False, indent=2))
    if not passed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
