#!/usr/bin/env python3
"""Execute the drug-side state/gate matrix without fabricating perturbations.

Disease-side signed effects are usable as a discovery reference.  The current
workspace contains no matched real compound perturbation or free-exposure
measurements.  Official MAP reference runs are retained, but the released
perturbation head is invariant across different drug encodings, so every
drug-effect and compatibility field remains explicitly non-estimable and no
compound is ranked or locked.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd
import yaml


MODULE_TARGETS = {
    "M1_Glycolysis": {"HK1", "HK2", "PFKFB3", "PFKFB4", "ALDOA", "ALDOC", "GAPDH", "PGK1", "PGAM1", "ENO1", "ENO2", "PKM", "LDHA", "PDK1", "SLC2A1", "HIF1A", "PFKP"},
    "M3_SASP_Senescence": {"CDKN1A", "CDKN2A", "IL6", "IL1B", "TNF", "CCL2", "CCL3", "CCL4", "CSF1", "CXCL10", "MMP9", "MMP14", "SERPINE1", "IGFBP3", "IGFBP7", "NFKB1", "RELA", "TGFB1"},
    "M5_Iron_Uptake": {"TFRC", "SLC11A2", "STEAP3", "DMT1"},
    "M5_Iron_Storage": {"FTH1", "FTL"},
    "M5_Iron_Export": {"SLC40A1", "CP"},
    "M5_Ferritinophagy": {"NCOA4", "ATG3", "ATG5", "ATG7", "MAP1LC3B", "SQSTM1"},
    "M6_LPO_ACSL4": {"ACSL4", "ALOX15", "ALOX5", "ALOX12", "NOX4", "GPX4", "LPCAT3", "SAT1"},
}


def normalise_cid(value: object) -> str:
    if pd.isna(value):
        return ""
    try:
        return str(int(float(value)))
    except (TypeError, ValueError):
        return str(value)


def run(root: Path) -> None:
    config = yaml.safe_load((root / "configs" / "drug_side_params.yaml").read_text(encoding="utf-8"))["drug_side"]
    tables = root / "results" / "tables"
    interim = root / "data" / "interim"
    tables.mkdir(parents=True, exist_ok=True)
    interim.mkdir(parents=True, exist_ok=True)

    master = pd.read_csv(tables / "drug_compound_master.tsv", sep="\t")
    targets = pd.read_csv(tables / "drug_target_evidence.tsv", sep="\t")
    modules = pd.read_csv(tables / "DiseaseModuleEffects_fast.tsv", sep="\t")
    package = json.loads((interim / "disease_interface_package_v1.3.json").read_text(encoding="utf-8"))
    ilincs_summary = json.loads((interim / "ilincs_candidate_summary.json").read_text(encoding="utf-8")) if (interim / "ilincs_candidate_summary.json").exists() else {"query_records": 0, "exact_signature_rows": 0, "unique_matched_cids": 0}
    ilincs_exact = pd.read_csv(interim / "ilincs_exact_candidate_signatures.tsv", sep="\t") if (interim / "ilincs_exact_candidate_signatures.tsv").exists() else pd.DataFrame(columns=["PubChem_CID"])
    map_audit_path = root / "results" / "2026-09-20-drug-side-execution" / "map_drug_sensitivity_audit.json"
    map_fp32_path = root / "results" / "2026-09-20-drug-side-execution" / "map_transformer_fp32_sensitivity.json"
    map_meta_path = root / "results" / "2026-09-20-drug-side-execution" / "map_reference_candidates" / "meta.json"
    map_outputs_path = root / "results" / "2026-09-20-drug-side-execution" / "map_reference_candidates" / "drug_side_map_reference_outputs.tsv"
    map_audit = json.loads(map_audit_path.read_text(encoding="utf-8")) if map_audit_path.exists() else {}
    map_fp32 = json.loads(map_fp32_path.read_text(encoding="utf-8")) if map_fp32_path.exists() else {}
    map_meta = json.loads(map_meta_path.read_text(encoding="utf-8")) if map_meta_path.exists() else {}
    map_outputs = pd.read_csv(map_outputs_path, sep="\t") if map_outputs_path.exists() else pd.DataFrame()
    map_reference_complete = int(map_meta.get("n_unique_pubchem_cids", 0)) == int(master["PubChem_CID"].nunique()) if not master.empty else False
    map_drug_sensitive = bool(map_audit.get("pred_hvg_max_abs_diff", 0.0) > 0.0)
    map_status = "completed_reference_run_but_drug_invariant" if map_reference_complete and not map_drug_sensitive else "not_available"
    hubs = {str(x).upper() for x in config["target_prior"]["literature_hubs"]}

    target_sets = {name: {gene.upper() for gene in genes} for name, genes in MODULE_TARGETS.items()}
    targets["Target_Symbol"] = targets["Target_Symbol"].astype(str).str.upper()
    master["PubChem_CID_norm"] = master["PubChem_CID"].map(normalise_cid)
    targets["PubChem_CID_norm"] = targets["PubChem_CID"].map(normalise_cid)
    ilincs_exact["PubChem_CID_norm"] = ilincs_exact["PubChem_CID"].map(normalise_cid) if "PubChem_CID" in ilincs_exact.columns else ""
    ilincs_counts = ilincs_exact.groupby("PubChem_CID_norm").size().to_dict() if not ilincs_exact.empty else {}

    target_by_compound = targets.groupby("PubChem_CID_norm")["Target_Symbol"].apply(lambda x: sorted(set(x))).to_dict()
    primary_endpoint_names = list(MODULE_TARGETS)
    effect_rows = []
    next_rows = []
    for _, compound in master.sort_values(["PubChem_CID_norm", "Ingredient"], na_position="last").iterrows():
        cid = compound["PubChem_CID_norm"]
        compound_id = f"CID:{cid}|InChIKey:{compound.get('InChIKey', '')}|Ingredient:{compound.get('Ingredient', '')}"
        compound_targets = set(target_by_compound.get(cid, []))
        hub_matches = sorted(compound_targets & hubs)
        module_matches = sorted(module for module, genes in target_sets.items() if compound_targets & genes)
        target_prior_status = "hub_or_disease_module_overlap_database_prior_only" if (hub_matches or module_matches) else "no_frozen_target_prior_overlap"
        ilincs_count = int(ilincs_counts.get(cid, 0))
        ilincs_status = "exact_pubchem_cid_signatures_available" if ilincs_count > 0 else "no_exact_pubchem_cid_match_in_iLINCS"

        for cell_type in ("Microglia", "Astrocyte"):
            for endpoint in primary_endpoint_names:
                hit = modules[(modules["cell_type"] == cell_type) & (modules["module"] == endpoint)]
                if hit.empty:
                    beta = None
                    fdr = None
                    disease_direction = "not_available"
                    coverage = None
                    disease_status = "endpoint_not_estimated"
                else:
                    row = hit.iloc[0]
                    beta = None if pd.isna(row.get("beta")) else float(row["beta"])
                    fdr = None if pd.isna(row.get("FDR")) else float(row["FDR"])
                    disease_direction = "AD_up" if beta is not None and beta > 0 else "AD_down" if beta is not None and beta < 0 else "no_direction"
                    coverage = None if pd.isna(row.get("coverage_ratio")) else float(row["coverage_ratio"])
                    disease_status = str(row.get("status", "estimated"))
                effect_rows.append({
                    "compound_id": compound_id,
                    "herb_prior": compound.get("Herb", ""),
                    "ingredient": compound.get("Ingredient", ""),
                    "PubChem_CID": cid,
                    "cell_type": cell_type,
                    "endpoint": endpoint,
                    "disease_beta_discovery": beta,
                    "disease_FDR_discovery": fdr,
                    "disease_direction": disease_direction,
                    "disease_endpoint_status": disease_status,
                    "disease_gene_coverage_ratio": coverage,
                    "drug_effect_status": "not_available_no_exact_pubchem_cid_match_in_iLINCS",
                    "drug_effect_direction": "not_available",
                    "free_exposure_status": "not_available",
                    "exposure_vs_potency_status": "not_estimable",
                    "VC_MAP_status": map_status,
                    "iLINCS_exact_signature_count": ilincs_count,
                    "iLINCS_status": ilincs_status,
                    "compatibility_status": "not_estimable",
                    "overshoot_status": "not_estimable",
                    "claim_level": "C0",
                    "formula_level_status": "blocked_formula_identity_exposure_and_whole_formula_perturbation",
                    "target_prior_role": target_prior_status,
                    "target_prior_hub_matches": ";".join(hub_matches),
                    "target_prior_module_matches": ";".join(module_matches),
                    "ranking_status": "not_ranked",
                })

        next_rows.append({
            "compound_id": compound_id,
            "herb_prior": compound.get("Herb", ""),
            "ingredient": compound.get("Ingredient", ""),
            "PubChem_CID": cid,
            "InChIKey": compound.get("InChIKey", ""),
            "database_target_count": int(len(compound_targets)),
            "target_prior_hub_matches": ";".join(hub_matches),
            "target_prior_module_matches": ";".join(module_matches),
            "chemical_space_status": "local_identity_database_prior_only",
            "observed_formula_status": "not_observed",
            "free_exposure_status": "not_available",
            "real_compound_perturbation_status": "not_available",
            "iLINCS_exact_signature_count": ilincs_count,
            "iLINCS_status": ilincs_status,
            "VC_MAP_status": map_status,
            "primary_rank_status": "not_ranked",
            "candidate_lock_status": "not_eligible",
            "claim_ceiling": "C0",
            "next_required_experiment": "matched_microglia_or_astrocyte_real_compound_perturbation_plus_free_exposure_and_QC",
        })

    evidence_matrix = pd.DataFrame(effect_rows)
    next_experiments = pd.DataFrame(next_rows).drop_duplicates(subset=["compound_id"])
    evidence_matrix.to_csv(tables / "drug_side_evidence_matrix.tsv", sep="\t", index=False)
    next_experiments.to_csv(tables / "drug_side_next_experiment.tsv", sep="\t", index=False)

    summary = {
        "formula_name": config["formula"]["name"],
        "formula_evidence_mode": config["formula"]["evidence_mode"],
        "formula_equivalence_status": config["formula"]["equivalence_status"],
        "disease_interface_package": package["package"],
        "disease_interface_discovery_status": package["discovery_status"],
        "strict_normalization_allowed": bool(package["normalization_claim_allowed"]),
        "candidate_compounds_in_local_prior": int(next_experiments.shape[0]),
        "candidate_compounds_with_database_target_prior_overlap": int((next_experiments[["target_prior_hub_matches", "target_prior_module_matches"]].fillna("").astype(str).agg("".join, axis=1).str.len() > 0).sum()),
        "iLINCS_query_records": int(ilincs_summary.get("query_records", 0)),
        "iLINCS_exact_signature_rows": int(ilincs_summary.get("exact_signature_rows", 0)),
        "iLINCS_exact_matched_cids": int(ilincs_summary.get("unique_matched_cids", 0)),
        "real_compound_perturbation_search_status": "iLINCS_queried_exact_pubchem_cid_no_match",
        "real_compound_perturbation_available": False,
        "free_exposure_available": False,
        "VC_MAP_available": bool(map_reference_complete and map_drug_sensitive),
        "VC_MAP_reference_runs": int(map_meta.get("n_unique_pubchem_cids", 0)),
        "VC_MAP_reference_output_rows": int(map_outputs.shape[0]),
        "VC_MAP_drug_sensitive": map_drug_sensitive,
        "VC_MAP_status": map_status,
        "VC_MAP_transformer_fp32_probe_max_abs_diff": float(map_fp32.get("transformer_cls_max_abs_diff", -1.0)),
        "directional_compatibility_eligible": 0,
        "primary_ranked_candidates": 0,
        "formula_level_eligible": 0,
        "claim_ceiling": "C0",
        "raw_sources_modified": False,
    }
    (interim / "drug_side_gate_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    args = parser.parse_args()
    run(args.root.resolve())
