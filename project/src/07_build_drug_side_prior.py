#!/usr/bin/env python3
"""Build the traceable drug-side chemical and target prior.

The two external CSVs are read-only local evidence caches.  Their records are
never promoted to observed formula composition, exposure, or compound
perturbation evidence without the corresponding primary data.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd
import yaml


FORMULA_HERBS = {
    "熟地黄": "\u719f\u5730\u9ec4",
    "丹参": "\u4e39\u53c2",
    "石菖蒲": "\u77f3\u83d6\u84b2",
    "紫河车": "\u7d2b\u6cb3\u8f66",
    "龟甲胶": "\u9f9f\u7532\u80f6",
    "茯神": "\u830f\u795e",
    "益智仁": "\u76ca\u667a\u4ec1",
}


def read_csv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, encoding="utf-8-sig", low_memory=False)


def norm_cid(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce").astype("Int64")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def run(root: Path) -> None:
    config = yaml.safe_load((root / "configs" / "drug_side_params.yaml").read_text(encoding="utf-8"))["drug_side"]
    pubchem_path = Path(config["local_prior"]["pubchem_identity_csv"])
    mapping_path = Path(config["local_prior"]["target_mapping_csv"])
    pubchem = read_csv(pubchem_path)
    mapping = read_csv(mapping_path)

    herb_values = set(FORMULA_HERBS.values())
    pubchem = pubchem[pubchem["Herb"].astype(str).isin(herb_values)].copy()
    mapping = mapping[mapping["Herb"].astype(str).isin(herb_values)].copy()
    pubchem["PubChem_CID_norm"] = norm_cid(pubchem["PubChem_CID"])
    mapping["PubChem_CID_norm"] = norm_cid(mapping["PubChem_CID"])

    out_tables = root / "results" / "tables"
    out_interim = root / "data" / "interim"
    out_tables.mkdir(parents=True, exist_ok=True)
    out_interim.mkdir(parents=True, exist_ok=True)

    formula_rows = []
    for label, herb in FORMULA_HERBS.items():
        pub_rows = pubchem[pubchem["Herb"].astype(str) == herb]
        map_rows = mapping[mapping["Herb"].astype(str) == herb]
        formula_rows.append({
            "formula_name": config["formula"]["name"],
            "evidence_mode": config["formula"]["evidence_mode"],
            "formula_equivalence_status": config["formula"]["equivalence_status"],
            "candidate_herb": label,
            "local_identity_records": int(pub_rows.shape[0]),
            "local_target_records": int(map_rows.shape[0]),
            "formula_observed_status": "not_observed_in_current_project",
            "exposure_status": "not_available",
            "inference_role": "candidate_formula_context_only; not proof of current fixed prescription",
        })
    pd.DataFrame(formula_rows).to_csv(out_tables / "drug_formula_candidate_herb_coverage.tsv", sep="\t", index=False)

    formula_identity = pd.DataFrame([{
        "formula_name": config["formula"]["name"],
        "evidence_mode": config["formula"]["evidence_mode"],
        "formula_equivalence_status": config["formula"]["equivalence_status"],
        "exact_prescription_available": False,
        "batch_and_preparation_qc_available": False,
        "whole_formula_perturbation_available": False,
        "formula_level_claim_allowed": False,
        "claim_ceiling_before_new_data": "C0_database_prior_or_hypothesis; no formula-level state claim",
        "decision": "continue_constituent_prior_audit_but_block_formula_level_inference",
    }])
    formula_identity.to_csv(out_tables / "drug_formula_identity.tsv", sep="\t", index=False)

    identity_cols = [
        "Herb", "Ingredient", "PubChem_CID_norm", "Compound_Title", "SMILES", "Canonical_SMILES",
        "Isomeric_SMILES", "InChI", "InChIKey", "Molecular_Formula", "Molecular_Weight", "XLogP",
        "TPSA", "Query_Status", "Query_Method",
    ]
    identity = pubchem[[c for c in identity_cols if c in pubchem.columns]].copy()
    identity = identity.rename(columns={"PubChem_CID_norm": "PubChem_CID"})
    for column in ["SMILES", "Canonical_SMILES", "Isomeric_SMILES", "InChIKey"]:
        if column not in identity.columns:
            identity[column] = pd.NA
    identity["Formula_Evidence_Mode"] = config["formula"]["evidence_mode"]
    identity["Formula_Equivalence"] = config["formula"]["equivalence_status"]
    identity["Observed_in_Formula"] = False
    identity["Parent_or_Metabolite"] = "unknown"
    identity["Exposure_Evidence"] = "not_available"
    identity["Mechanism_Track"] = "unassigned"
    identity["Chemical_Evidence_Class"] = "local_PubChem_identity_cache_database_prior"

    map_identity = mapping[[c for c in ["Herb", "Ingredient", "PubChem_CID_norm", "SMILES", "Canonical_SMILES", "Isomeric_SMILES", "InChIKey"] if c in mapping.columns]].copy()
    map_identity = map_identity.rename(columns={"PubChem_CID_norm": "PubChem_CID"})
    identity = pd.concat([identity, map_identity], ignore_index=True, sort=False)
    identity = identity.drop_duplicates(subset=["PubChem_CID", "InChIKey", "Ingredient"], keep="first")
    identity = identity.sort_values(["Herb", "PubChem_CID", "Ingredient"], na_position="last")
    identity.to_csv(out_tables / "drug_compound_master.tsv", sep="\t", index=False)

    evidence_cols = [
        "Herb", "Ingredient", "PubChem_CID_norm", "Target_Symbol", "Target_EntrezID", "Target_Source",
        "Target_Source_Tier", "Source", "Evidence_DOI", "Target_Species", "Target_Species_Status",
        "Final_Target_Source_DB", "Target_DB_Provider", "Target_DB_Notes",
    ]
    evidence = mapping[[c for c in evidence_cols if c in mapping.columns]].copy()
    evidence = evidence.rename(columns={"PubChem_CID_norm": "PubChem_CID"})
    evidence["Target_Direction"] = "not_available"
    evidence["Potency_Value"] = pd.NA
    evidence["Potency_Unit"] = pd.NA
    evidence["Exposure_Free_Concentration"] = pd.NA
    evidence["Exposure_vs_Potency"] = "not_estimable"
    evidence["Evidence_Independence_Count"] = 1
    evidence["Evidence_Role"] = "database_target_prior_not_compound_reversal_evidence"
    evidence = evidence.drop_duplicates(subset=["PubChem_CID", "Target_Symbol", "Source"], keep="first")
    evidence = evidence.sort_values(["Herb", "PubChem_CID", "Target_Symbol"], na_position="last")
    evidence.to_csv(out_tables / "drug_target_evidence.tsv", sep="\t", index=False)

    provenance = {
        "formula_identity": "Mode B; exact current prescription and preparation not resolved",
        "sources": {
            str(pubchem_path): {"sha256": sha256(pubchem_path), "records_used": int(pubchem.shape[0]), "role": "local chemical identity cache"},
            str(mapping_path): {"sha256": sha256(mapping_path), "records_used": int(mapping.shape[0]), "role": "local target mapping cache"},
        },
        "formula_candidate_herbs": list(FORMULA_HERBS),
        "local_evidence_herbs": sorted({str(x) for x in mapping["Herb"].dropna().unique()}),
        "raw_sources_modified": False,
        "promotion_blockers": ["exact_prescription_and_FE-A/B verification", "in-vivo exposure evidence", "real compound perturbation", "whole-formula perturbation"],
    }
    (out_interim / "drug_side_prior_provenance.json").write_text(json.dumps(provenance, ensure_ascii=False, indent=2), encoding="utf-8")
    summary = {
        "formula_candidate_herbs": len(FORMULA_HERBS),
        "herbs_with_local_identity": int(sum(int(x > 0) for x in formula_identity["exact_prescription_available"]) if False else len(set(pubchem["Herb"].dropna()))),
        "chemical_identity_rows": int(identity.shape[0]),
        "unique_compound_cids": int(identity["PubChem_CID"].dropna().nunique()),
        "target_evidence_rows": int(evidence.shape[0]),
        "unique_target_symbols": int(evidence["Target_Symbol"].dropna().nunique()),
        "formula_level_claim_allowed": False,
        "raw_sources_modified": False,
    }
    (out_interim / "drug_side_prior_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    args = parser.parse_args()
    run(args.root.resolve())
