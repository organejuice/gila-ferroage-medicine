#!/usr/bin/env python3
"""Audit drug-side gate invariants and generate a machine-readable report."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def run(root: Path) -> None:
    tables = root / "results" / "tables"
    interim = root / "data" / "interim"
    provenance = json.loads((interim / "drug_side_prior_provenance.json").read_text(encoding="utf-8"))
    gate = json.loads((interim / "drug_side_gate_summary.json").read_text(encoding="utf-8"))
    ilincs = json.loads((interim / "ilincs_candidate_summary.json").read_text(encoding="utf-8"))
    interface = json.loads((interim / "disease_interface_package_v1.3.json").read_text(encoding="utf-8"))
    execution_dir = root / "results" / "2026-09-20-drug-side-execution"
    map_meta = json.loads((execution_dir / "map_reference_candidates" / "meta.json").read_text(encoding="utf-8")) if (execution_dir / "map_reference_candidates" / "meta.json").exists() else {}
    map_audit = json.loads((execution_dir / "map_drug_sensitivity_audit.json").read_text(encoding="utf-8")) if (execution_dir / "map_drug_sensitivity_audit.json").exists() else {}
    map_fp32 = json.loads((execution_dir / "map_transformer_fp32_sensitivity.json").read_text(encoding="utf-8")) if (execution_dir / "map_transformer_fp32_sensitivity.json").exists() else {}
    matrix = pd.read_csv(tables / "drug_side_evidence_matrix.tsv", sep="\t")
    next_experiment = pd.read_csv(tables / "drug_side_next_experiment.tsv", sep="\t")
    modules = pd.read_csv(tables / "DiseaseModuleEffects_fast.tsv", sep="\t")
    qc = pd.read_csv(tables / "D3_fast_qc_audit.tsv", sep="\t")
    pseudo_meta = pd.concat([
        pd.read_csv(interim / "D3_pseudobulk_metadata_microglia.tsv", sep="\t"),
        pd.read_csv(interim / "D3_pseudobulk_metadata_astrocyte.tsv", sep="\t"),
    ], ignore_index=True)
    checks: list[dict] = []

    def check(name: str, passed: bool, detail: str) -> None:
        checks.append({"check": name, "passed": bool(passed), "detail": detail})

    check("formula_claim_blocked", gate["formula_level_eligible"] == 0 and interface["normalization_claim_allowed"] is False, "formula-level and strict-normalization claims remain blocked")
    check("real_perturbation_not_fabricated", bool(matrix["drug_effect_status"].astype(str).str.startswith("not_available").all()), "all drug-effect rows are explicit unavailable statuses")
    check("iLINCS_query_complete", ilincs["query_records"] == 17 and ilincs["exact_signature_rows"] == 0 and ilincs["unique_matched_cids"] == 0, "all 17 local compound records were queried and no exact CID match was promoted")
    check("exposure_not_fabricated", bool((matrix["free_exposure_status"] == "not_available").all()), "all free-exposure rows are explicit unavailable statuses")
    vc_allowed_status = {"not_available", "completed_reference_run_but_drug_invariant"}
    check("VC_not_fabricated", bool(matrix["VC_MAP_status"].isin(vc_allowed_status).all()), "virtual-cell outputs are either unavailable or explicitly marked invariant and not used for ranking")
    check("MAP_reference_runs_complete", map_meta.get("n_unique_pubchem_cids", 0) == 16 and map_audit.get("drug_vector_max_abs_diff", 0) > 0 and map_audit.get("pred_hvg_max_abs_diff", 0) == 0, "16 unique compounds were run with the official checkpoint; drug encodings vary but the released perturbation head is invariant in the reference probe")
    check("MAP_transformer_fp32_probe", map_fp32.get("dtype") == "float32" and map_fp32.get("transformer_cls_max_abs_diff") == 0, "independent FP32 probe confirms the released transformer CLS output is invariant to a synthetic drug-token perturbation")
    check("no_candidate_ranking", bool((next_experiment["primary_rank_status"] == "not_ranked").all() and (matrix["ranking_status"] == "not_ranked").all()), "no compound was ranked or locked")
    check("claim_ceiling_C0", bool((matrix["claim_level"] == "C0").all() and (next_experiment["claim_ceiling"] == "C0").all()), "current evidence is limited to database prior/hypothesis level")
    check("M7_included", "M7_FerroAging_Author" in set(modules["module"]), "D5 Table S2 ferro-aging gene set is represented")
    check("D3_sample_count", int(qc.shape[0]) == 24 and int(pseudo_meta.shape[0]) == 48, "24 samples and 48 glial pseudobulk groups are present")
    qc_max = int(qc["n_cells_after"].max())
    group_ok = True
    for _, row in pseudo_meta.iterrows():
        if int(row["n_cells"]) <= 0 or int(row["n_cells"]) > qc_max:
            group_ok = False
    check("pseudobulk_cell_counts_plausible", group_ok, "pseudobulk n_cells are positive and bounded by sample QC counts")
    source_hash_ok = True
    for source, record in provenance["sources"].items():
        current = sha256(Path(source))
        source_hash_ok = source_hash_ok and current == record["sha256"]
    check("local_source_hash_recheck", source_hash_ok, "local prior source hashes match the hashes recorded at ingestion")
    check("raw_sources_declared_read_only", provenance["raw_sources_modified"] is False and gate["raw_sources_modified"] is False, "source caches were read without modification")

    report = {
        "audit_status": "PASS_WITH_PREDECLARED_BLOCKS" if all(x["passed"] for x in checks) else "FAIL",
        "checks": checks,
        "failed_checks": [x["check"] for x in checks if not x["passed"]],
        "blocked_claims": ["strict_normalization", "compound_directional_reversal", "exposure_potency_match", "VC_prediction", "formula_level_state_modulation"],
        "evidence_counts": {"compound_prior_rows": int(next_experiment.shape[0]), "matrix_rows": int(matrix.shape[0]), "target_prior_overlap_rows": int(gate["candidate_compounds_with_database_target_prior_overlap"])},
        "reproducibility": {"seed": 42, "raw_sources_modified": False, "source_hash_rechecked": source_hash_ok},
    }
    (interim / "drug_side_audit.json").write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    pd.DataFrame(checks).to_csv(tables / "drug_side_audit_checks.tsv", sep="\t", index=False)
    print(json.dumps(report, ensure_ascii=False, indent=2))
    if report["audit_status"] == "FAIL":
        raise SystemExit(1)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    args = parser.parse_args()
    run(args.root.resolve())
