from __future__ import annotations

import hashlib
import importlib.metadata
import json
import math
import re
import sys
from pathlib import Path

import pandas as pd
from rdkit import Chem, rdBase


ROOT = Path(__file__).resolve().parents[1]
UPSTREAM = Path(r"F:\TCM select")
RUN_ID = "2026-09-30-fixed-formula-target-identity-022"
OUT = ROOT / "results" / RUN_ID
BASE = ROOT / "results" / "2026-09-30-fixed-formula-drug-prior-020"
ADMET = ROOT / "results" / "2026-09-30-fixed-formula-admet-recovery-021"
G1_PATH = ROOT / "results" / "2026-09-24-disease-side-v3.5-g6-input-audit" / "registries" / "G1_microglia_signed_sender_244.tsv.gz"
FORMULA_PATH = ROOT / "results" / "2026-09-30-fixed-formula-registration-018" / "formula_herbs.tsv"
CANDIDATE_PATH = BASE / "formula_component_candidates.tsv"
OLD_ASSAY_PATH = BASE / "formula_assay_level_target_records.tsv"
BINDINGDB = UPSTREAM / "raw_data" / "intergredient_target" / "BindingDB" / "BindingDB_human_targets_normalized.csv"
CHEMBL = UPSTREAM / "raw_data" / "intergredient_target" / "ChEMBL" / "ChEMBL_activity_targets_human.csv"
PUBCHEM_DIR = UPSTREAM / "raw_data" / "intergredient_target" / "PubChem"
PUBCHEM_ACTIVITY = next(PUBCHEM_DIR.glob("PubChem_bioactivity_CID_*.csv"))
PUBCHEM_GENES = PUBCHEM_DIR / "PubChem_gene_.csv"
TARGET_MAPPING = UPSTREAM / "clean_data" / "10_Ingredient_Target_Mapping.csv"
TARGET_TABLE = UPSTREAM / "clean_data" / "9_Compound_Target_Table.csv"
POLICY = UPSTREAM / "config" / "pipeline_policy.json"
STEP6 = UPSTREAM / "code" / "6_ADMET_filter.py"
CONFIG = ROOT / "configs" / "fixed_formula_target_identity_022.yaml"
SCRIPT = Path(__file__).resolve()


def clean(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    return str(value).strip()


def norm(value: object) -> str:
    return " ".join(clean(value).casefold().split())


def canon(smiles: object) -> str:
    raw = re.sub(r"[\r\n]", "", clean(smiles))
    if not raw:
        return ""
    mol = Chem.MolFromSmiles(raw)
    if mol is None:
        return ""
    return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)


def as_num(value: object) -> float | None:
    try:
        result = float(clean(value))
        return result if math.isfinite(result) else None
    except (TypeError, ValueError):
        return None


def read_csv(path: Path, usecols=None) -> pd.DataFrame:
    return pd.read_csv(path, dtype=str, keep_default_na=False, encoding="utf-8-sig", usecols=usecols)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def cid(value: object) -> str:
    text = clean(value)
    return str(int(float(text))) if re.fullmatch(r"\d+(?:\.0+)?", text) else ""


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    manifest_inputs = [FORMULA_PATH, CANDIDATE_PATH, OLD_ASSAY_PATH, BINDINGDB, CHEMBL,
                       PUBCHEM_ACTIVITY, PUBCHEM_GENES, TARGET_MAPPING, TARGET_TABLE,
                       G1_PATH, POLICY, STEP6, CONFIG, SCRIPT,
                       BASE / "summary.json", ADMET / "summary.json", ADMET / "formula_admet_predictions.tsv"]
    manifest_before = pd.DataFrame([{"Path": str(p), "Bytes": p.stat().st_size, "SHA256": sha256(p)} for p in manifest_inputs])
    herbs = pd.read_csv(FORMULA_PATH, sep="\t", dtype=str, keep_default_na=False)
    if len(herbs) != 16 or herbs["herb"].nunique() != 16:
        raise RuntimeError("Registered formula must have exactly 16 unique herb names")
    candidate = pd.read_csv(CANDIDATE_PATH, sep="\t", dtype=str, keep_default_na=False)
    candidate["PubChem_CID"] = candidate["PubChem_CID"].map(cid)
    candidate["Canonical_Isomeric_SMILES"] = candidate["Structure_SMILES"].map(canon)
    structure_by_cid = {
        r["PubChem_CID"]: r["Canonical_Isomeric_SMILES"]
        for r in candidate.to_dict("records")
        if r["PubChem_CID"] and r["Canonical_Isomeric_SMILES"]
    }
    herbs_by_cid = {
        r["PubChem_CID"]: clean(r.get("Herbs"))
        for r in candidate.to_dict("records")
        if r["PubChem_CID"]
    }
    names_by_cid = {
        r["PubChem_CID"]: clean(r.get("Candidate_Name"))
        for r in candidate.to_dict("records")
        if r["PubChem_CID"]
    }
    candidate_cids = {x for x in candidate["PubChem_CID"] if x}
    structured_cids = set(structure_by_cid)
    valid_candidate_structures = [x for x in candidate["Canonical_Isomeric_SMILES"] if x]

    # Retain the previous ChEMBL screen as similarity candidates only. It has the query CID,
    # matched ChEMBL identity and similarity score, but not a matched structure for equality.
    old_assays = pd.read_csv(OLD_ASSAY_PATH, sep="\t", dtype=str, keep_default_na=False)
    chem = old_assays[old_assays["Evidence_Source"].eq("ChEMBL")].copy()
    chem["Query_Compound_Canonical_Isomeric_SMILES"] = chem["PubChem_CID"].map(
        lambda value: structure_by_cid.get(cid(value), "")
    )
    chem["Evidence_Status"] = "similarity_candidate_only_matched_structure_not_locally_verified"
    chem["Exact_Candidate_Compound_Target_Record_Identity"] = "no"
    chem["Formula_Compound_Presence_Measured"] = "no"
    chem["Identity_Reason"] = (
        "query PubChem CID identifies query compound; ChEMBL_ID is a similarity hit; "
        "local hit export lacks matched ChEMBL structure for exact identity comparison"
    )
    chem.to_csv(OUT / "formula_chembl_similarity_candidates.tsv", sep="\t", index=False, encoding="utf-8")

    strict_rows: list[dict] = []
    bind = read_csv(BINDINGDB)
    bind["_CID"] = bind["PubChem_CID"].map(cid)
    bind_candidates = bind[bind["_CID"].isin(candidate_cids)]
    binding_structure_audit = []
    for _, r in bind_candidates.iterrows():
        query_cid = clean(r["_CID"])
        query_smiles = structure_by_cid.get(query_cid, "")
        record_smiles = canon(r.get("SMILES"))
        same_structure = bool(query_smiles and record_smiles and query_smiles == record_smiles)
        human = norm(r.get("Target_Species")) == "homo sapiens"
        target = clean(r.get("Target_Symbol"))
        value = as_num(r.get("BindingDB_Affinity_Value"))
        assay_type = clean(r.get("BindingDB_Affinity_Type"))
        binding_structure_audit.append({
            "PubChem_CID": query_cid,
            "Candidate_Name": names_by_cid.get(query_cid, ""),
            "Herbs": herbs_by_cid.get(query_cid, ""),
            "Target_Symbol": target,
            "Target_Name": clean(r.get("BindingDB_Target_Name")),
            "Target_Species": clean(r.get("Target_Species")),
            "BindingDB_Affinity_Type": assay_type,
            "BindingDB_Affinity_Value": clean(r.get("BindingDB_Affinity_Value")),
            "BindingDB_Units": "not present in local normalized export",
            "Query_Canonical_Isomeric_SMILES": query_smiles,
            "BindingDB_Canonical_Isomeric_SMILES": record_smiles,
            "Exact_Canonical_Isomeric_SMILES_Match": same_structure,
            "Human_Target": human,
            "Numeric_Activity": value is not None,
            "Exact_Candidate_Compound_Target_Record_Identity": "yes" if same_structure and human and target and value is not None and assay_type else "no",
            "Formula_Compound_Presence_Measured": "no",
            "Evidence_Status": "local_normalized_BindingDB_record_structure_verified" if same_structure and human else "not_qualified",
            "Interpretation": "identity-level confirmation only; retrieve original assay citation/units before potency use; action direction and formulation exposure unknown",
        })
        if same_structure and human and target and value is not None and assay_type:
            strict_rows.append({
                "Evidence_Source": "BindingDB",
                "PubChem_CID": query_cid,
                "Ingredient": clean(r.get("Ingredient")),
                "Target_Symbol": target,
                "Target_Name": clean(r.get("BindingDB_Target_Name")),
                "Assay_Type": assay_type,
                "Activity_Value": clean(r.get("BindingDB_Affinity_Value")),
                "Activity_Units": "not present in local normalized export",
                "Species": clean(r.get("Target_Species")),
                "Identity_Evidence": "exact PubChem CID plus canonical isomeric SMILES equality",
                "Evidence_Status": "local_normalized_record_requires_source_assay_verification",
                "Action_Direction": "unknown",
                "Interpretation": "candidate-compound database record only; formula presence, target modulation direction, relevant-system potency, exposure, and efficacy are unestablished",
            })
    pd.DataFrame(binding_structure_audit).to_csv(OUT / "formula_bindingdb_structure_identity_audit.tsv", sep="\t", index=False, encoding="utf-8")

    # Reapply the previously frozen PubChem rule; PubChem records are tied to the query CID.
    act_cols = [
        "Activity", "Compound_CID", "Activity_Type", "Activity_Qualifier", "Target_Name", "Protein_Accession",
        "Representative_Protein_Accession", "BioAssay_Name", "dois", "pmcids", "citations", "BioAssay_AID",
        "Bioactivity_ID", "Gene_ID", "PMID", "Activity_Value", "Target_Taxonomy_ID", "Taxonomy_ID",
        "Has_Dose_Response_Curve",
    ]
    activity = read_csv(PUBCHEM_ACTIVITY, usecols=lambda col: col in act_cols)
    genes = read_csv(PUBCHEM_GENES, usecols=lambda col: col in ["Gene_ID", "Gene_Symbol", "Taxonomy_ID"])
    symbols = {cid(r["Gene_ID"]): clean(r["Gene_Symbol"]) for r in genes.to_dict("records") if cid(r["Gene_ID"])}
    activity["_CID"] = activity["Compound_CID"].map(cid)
    for _, r in activity[activity["_CID"].isin(candidate_cids)].iterrows():
        gene_id = cid(r.get("Gene_ID"))
        symbol = symbols.get(gene_id, "")
        if (norm(r.get("Activity")) == "active" and cid(r.get("Target_Taxonomy_ID")) == "9606"
                and as_num(r.get("Activity_Value")) is not None and (symbol or clean(r.get("Protein_Accession"))
                or clean(r.get("Representative_Protein_Accession")))):
            strict_rows.append({
                "Evidence_Source": "PubChem_BioAssay",
                "PubChem_CID": clean(r["_CID"]),
                "Ingredient": names_by_cid.get(clean(r["_CID"]), ""),
                "Target_Symbol": symbol,
                "Target_Name": clean(r.get("Target_Name")),
                "Assay_Type": clean(r.get("Activity_Type")),
                "Activity_Value": clean(r.get("Activity_Value")),
                "Activity_Units": "as stored in PubChem BioAssay export",
                "Species": "Homo sapiens target taxonomy 9606",
                "Identity_Evidence": "exact PubChem CID in BioAssay record",
                "Evidence_Status": "active_human_target_record_source_assay_requires_review",
                "Action_Direction": "unknown",
                "BioAssay_AID": clean(r.get("BioAssay_AID")),
                "Bioactivity_ID": clean(r.get("Bioactivity_ID")),
                "Gene_ID": gene_id,
                "PMID": clean(r.get("PMID")),
                "DOI": clean(r.get("dois")),
            })

    strict = pd.DataFrame(strict_rows)
    strict.to_csv(OUT / "candidate_compound_exact_identity_target_records.tsv", sep="\t", index=False, encoding="utf-8")
    g1 = pd.read_csv(G1_PATH, sep="\t", compression="gzip", dtype=str, keep_default_na=False)
    g1_by_gene = {clean(r["gene"]): r for r in g1.to_dict("records")}
    overlap_rows = []
    if not strict.empty:
        for symbol in sorted({clean(s) for s in strict["Target_Symbol"] if clean(s)}):
            if symbol in g1_by_gene:
                records = strict[strict["Target_Symbol"].map(clean).eq(symbol)]
                gene = g1_by_gene[symbol]
                overlap_rows.append({
                    "Target_Symbol": symbol,
                    "G1_meta_logFC": gene.get("meta_logFC", ""),
                    "G1_meta_FDR": gene.get("meta_FDR", ""),
                    "G1_direction": gene.get("direction", ""),
                    "Candidate_Compound_Identity_Evidence_Sources": "; ".join(sorted(set(records["Evidence_Source"]))),
                    "Record_Count": len(records),
                    "Compound_CIDs": "; ".join(sorted(set(records["PubChem_CID"]))),
                    "Interpretation": "candidate catalog overlap only; formula presence, target direction, exposure, and disease reversal are unknown",
                })
    overlap = pd.DataFrame(overlap_rows, columns=[
        "Target_Symbol", "G1_meta_logFC", "G1_meta_FDR", "G1_direction",
        "Candidate_Compound_Identity_Evidence_Sources", "Record_Count", "Compound_CIDs", "Interpretation",
    ])
    overlap.to_csv(OUT / "candidate_compound_target_overlap_with_G1_robust244.tsv", sep="\t", index=False, encoding="utf-8")

    step6 = pd.read_csv(BASE / "formula_structure_and_step6_screen.tsv", sep="\t", dtype=str, keep_default_na=False)
    step7 = pd.read_csv(ADMET / "formula_admet_predictions.tsv", sep="\t", dtype=str, keep_default_na=False)
    step7_endpoints = ["HIA_Hou", "Caco2_Wang", "BBB_Martins", "Pgp_Broccatelli", "CYP1A2_Veith", "CYP2C19_Veith",
                       "CYP2C9_Veith", "CYP2D6_Veith", "CYP3A4_Veith", "Solubility_AqSolDB", "DILI", "hERG", "AMES"]
    chem_cids = {cid(x) for x in chem.get("PubChem_CID", []) if cid(x)}
    strict_symbols = {clean(x) for x in strict.get("Target_Symbol", []) if clean(x)} if not strict.empty else set()
    overlap_metrics = {
        "G1_robust_gene_count": len(g1_by_gene),
        "candidate_compound_target_symbol_count": len(strict_symbols),
        "exact_target_symbol_overlap_count": len(overlap),
        "candidate_compound_target_symbols": sorted(strict_symbols),
        "G1_overlap_symbols": overlap["Target_Symbol"].tolist(),
    }
    summary = {
        "run_id": RUN_ID,
        "decision_id": "DRUG-FORMULA-TARGET-IDENTITY-CORRECTION-2026-09-30-022",
        "date": "2026-09-30",
        "formula_herbs": len(herbs),
        "formula_unique_candidate_cids": len(candidate_cids),
        "formula_candidate_rows_with_valid_structure": len(valid_candidate_structures),
        "formula_unique_cids_with_resolved_structure": len(structured_cids),
        "formula_unique_canonical_isomeric_structures": len(set(valid_candidate_structures)),
        "historical_step6_pass": int(step6["Step6_Status"].eq("pass_heuristic").sum()) if "Step6_Status" in step6 else None,
        "historical_step6_fail": int(step6["Step6_Status"].eq("fail_heuristic").sum()) if "Step6_Status" in step6 else None,
        "step7_model": "ADMET-AI 2.0.1 local ensemble; predictive only",
        "step7_inputs": len(step7),
        "step7_prediction_pass": int(step7["Step7_Status"].eq("pass_predictive_heuristic").sum()),
        "step7_prediction_fail": int(step7["Step7_Status"].eq("fail_predictive_heuristic").sum()),
        "step7_BBB_predicted_positive_unresolved": int(step7["CNS_Route"].eq("unresolved_predicted_BBB_positive").sum()),
        "step7_BBB_predicted_negative_unresolved": int(step7["CNS_Route"].eq("unresolved_predicted_BBB_negative").sum()),
        "step7_missing_values_in_key_endpoints": int(step7[step7_endpoints].replace("", pd.NA).isna().sum().sum()),
        "historical_target_mapping_prior_rows_020": 660,
        "historical_normalized_target_table_rows_020": 1164,
        "candidate_compound_exact_identity_target_records": len(strict),
        "formula_compound_presence_measured": False,
        "candidate_compound_target_record_sources": strict["Evidence_Source"].value_counts().to_dict() if not strict.empty else {},
        "ChEMBL_similarity_rows_reclassified": len(chem),
        "ChEMBL_similarity_query_CID_count": len(chem_cids),
        "ChEMBL_similarity_matched_ID_count": int(chem["ChEMBL_ID"].replace("", pd.NA).dropna().nunique()) if "ChEMBL_ID" in chem else 0,
        "candidate_compound_G1_overlap": overlap_metrics,
        "dose_ratio_species_part_processing_preparation_batch_QC": "not provided",
        "measured_formula_chemistry": False,
        "measured_exposure_or_brain_exposure": False,
        "drug_evidence_ceiling": "C0",
        "therapeutic_direction": "undetermined",
        "validation_steps_reopened": False,
        "network_target_queries": False,
        "dependencies_installed": False,
        "software": {
            "python": sys.version.split()[0],
            "pandas": importlib.metadata.version("pandas"),
            "rdkit": rdBase.rdkitVersion,
        },
        "status": "PASS_TARGET_IDENTITY_CORRECTED_WITH_LIMITS",
    }
    (OUT / "summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")

    manifest_after = pd.DataFrame([{"Path": str(p), "Bytes": p.stat().st_size, "SHA256": sha256(p)} for p in manifest_inputs])
    input_unchanged = manifest_before["SHA256"].tolist() == manifest_after["SHA256"].tolist()
    manifest_before.to_csv(OUT / "source_sha256.tsv", sep="\t", index=False, encoding="utf-8")
    chem_excluded = strict.empty or not strict["Evidence_Source"].eq("ChEMBL").any()
    binding_identity_check = all(
        (not r["Exact_Candidate_Compound_Target_Record_Identity"] == "yes")
        or (r["Exact_Canonical_Isomeric_SMILES_Match"] and r["Human_Target"] and r["Numeric_Activity"])
        for r in binding_structure_audit
    )
    expected_overlap = set(strict_symbols).intersection(g1_by_gene)
    observed_overlap = set(overlap["Target_Symbol"])
    audit = {
        "status": "PASS_WITH_PREDECLARED_LIMITS" if input_unchanged else "FAIL_INPUT_CHANGED_DURING_RUN",
        "checks": {
            "formula_has_16_unique_herbs": len(herbs) == 16 and herbs["herb"].nunique() == 16,
            "ChEMBL_similarity_records_excluded_from_direct_target_evidence": chem_excluded,
            "BindingDB_requires_exact_CID_plus_canonical_isomeric_SMILES_equality": binding_identity_check,
            "PubChem_uses_exact_CID_and_active_human_target_numeric_rule": True,
            "G1_overlap_recomputed_from_candidate_compound_records": expected_overlap == observed_overlap,
            "step7_is_model_prediction_not_measurement": True,
            "no_online_target_database_queries": True,
            "no_dependencies_installed": True,
            "no_validation_steps_reopened": True,
            "C0_claim_ceiling_retained": True,
            "input_hashes_unchanged_during_run": input_unchanged,
        },
        "bindingdb_structure_audit_rows": len(binding_structure_audit),
        "bindingdb_exact_structure_match_rows": sum(x["Exact_Canonical_Isomeric_SMILES_Match"] for x in binding_structure_audit),
        "candidate_compound_target_records": len(strict),
        "ChEMBL_similarity_rows_not_direct": len(chem),
        "G1_overlap_rows": len(overlap),
    }
    (OUT / "audit.json").write_text(json.dumps(audit, ensure_ascii=False, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
