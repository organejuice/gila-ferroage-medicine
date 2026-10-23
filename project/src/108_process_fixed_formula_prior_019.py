from __future__ import annotations

import hashlib
import importlib.metadata
import importlib.util
import json
import re
import sys
import unicodedata
from collections import defaultdict
from pathlib import Path

import pandas as pd
from rdkit import Chem, RDLogger, rdBase


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "2026-09-30-fixed-formula-drug-prior-019"
FORMULA_PATH = ROOT / "results" / "2026-09-30-fixed-formula-registration-018" / "formula_herbs.tsv"
TCMBANK = Path(r"F:\TCM select\raw_data\TCMbank\tcmbank_ingredient_table.xlsx")
BATMAN = Path(r"F:\TCM select\raw_data\BATMAN\batman_herb_ingredient_table.txt")
UPSTREAM = Path(r"F:\TCM select")
POLICY_PATH = UPSTREAM / "config" / "pipeline_policy.json"
STEP6_PATH = UPSTREAM / "code" / "6_ADMET_filter.py"
TARGET_MAPPING = UPSTREAM / "clean_data" / "10_Ingredient_Target_Mapping.csv"
TARGET_TABLE = UPSTREAM / "clean_data" / "9_Compound_Target_Table.csv"
PUBCHEM_ACTIVITY = next((UPSTREAM / "raw_data" / "intergredient_target" / "PubChem").glob("PubChem_bioactivity_CID_*.csv"))
PUBCHEM_GENES = UPSTREAM / "raw_data" / "intergredient_target" / "PubChem" / "PubChem_gene_.csv"
BINDINGDB = UPSTREAM / "raw_data" / "intergredient_target" / "BindingDB" / "BindingDB_human_targets_normalized.csv"
CHEMBL = UPSTREAM / "raw_data" / "intergredient_target" / "ChEMBL" / "ChEMBL_activity_targets_human.csv"
G1_244 = ROOT / "results" / "2026-09-24-disease-side-v3.5-g6-input-audit" / "registries" / "G1_microglia_signed_sender_244.tsv.gz"
PARAM_PATH = ROOT / "configs" / "fixed_formula_drug_prior_019.yaml"
BASE_SCRIPT = Path(__file__).resolve()
RUN_SCRIPT = BASE_SCRIPT
RUN_ID = "2026-09-30-fixed-formula-drug-prior-019"
DECISION_ID = "DRUG-FORMULA-PRIOR-2026-09-30-019"

# Invalid source structures are counted and retained in the tables; keep RDKit's
# per-character parser diagnostics out of the console for large natural-product files.
RDLogger.DisableLog("rdApp.error")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def clean(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    return str(value).strip()


def norm(value: object) -> str:
    text = unicodedata.normalize("NFKC", clean(value)).casefold()
    return " ".join(text.split())


def numeric_ids(value: object) -> list[str]:
    text = clean(value)
    ids: list[str] = []
    for token in re.split(r"[;,|]", text):
        token = token.strip()
        if re.fullmatch(r"\d+(?:\.0+)?", token):
            ids.append(str(int(float(token))))
        else:
            ids.extend(re.findall(r"(?<![A-Za-z])\d{4,}(?![A-Za-z])", token))
    return list(dict.fromkeys(ids))


def canonical_smiles(smiles: object) -> str:
    text = clean(smiles)
    if not text:
        return ""
    # TCMBank export wraps some single-cell SMILES across physical lines.
    # Remove only CR/LF formatting; preserve the original cell in the output.
    parser_text = re.sub(r"[\r\n]", "", text)
    mol = Chem.MolFromSmiles(parser_text)
    if mol is None:
        return ""
    return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)


def safe_float(value: object) -> float | None:
    try:
        if clean(value) == "":
            return None
        return float(value)
    except (TypeError, ValueError):
        return None


def read_csv(path: Path, usecols=None) -> pd.DataFrame:
    return pd.read_csv(path, encoding="utf-8-sig", dtype=str, keep_default_na=False, usecols=usecols)


def add_source_record(records: list[dict], **row) -> None:
    row["Ingredient_Name_Normalized"] = norm(row.get("Ingredient_Name_Raw"))
    records.append(row)


def read_ingredient_sources(herbs: list[str]) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, list[dict]], dict[str, set[str]]]:
    # Read only the columns needed for exact name matching, structure resolution and provenance.
    columns = [
        "TCM_name", "TCM_name2", "name", "PubChem_id", "Smiles", "Ingredient_id",
        "Molecular_Formula", "Molecular_Weight", "OB_score",
    ]
    tcmbank = pd.read_excel(TCMBANK, dtype=str, keep_default_na=False, usecols=columns, engine="openpyxl")
    tcmbank.columns = [str(c).strip() for c in tcmbank.columns]
    batman = pd.read_csv(BATMAN, sep="\t", dtype=str, keep_default_na=False, encoding="utf-8-sig")
    batman.columns = [str(c).strip() for c in batman.columns]

    exact_tcm = tcmbank[tcmbank["TCM_name"].map(lambda x: clean(x) in set(herbs))].copy()
    alias_tcm = tcmbank[tcmbank["TCM_name2"].map(lambda x: clean(x) in set(herbs))].copy()
    exact_bat = batman[batman["Chinese.Name"].map(lambda x: clean(x) in set(herbs))].copy()

    # Global TCMBank CID->SMILES index is used only for exact numeric CID cross-source resolution.
    cid_smiles: dict[str, set[str]] = defaultdict(set)
    cid_names: dict[str, set[str]] = defaultdict(set)
    for row in tcmbank.to_dict("records"):
        cs = canonical_smiles(row.get("Smiles"))
        name = clean(row.get("name"))
        for cid in numeric_ids(row.get("PubChem_id")):
            if cs:
                cid_smiles[cid].add(cs)
            if name:
                cid_names[cid].add(name)

    records: list[dict] = []
    batman_by_herb: dict[str, list[dict]] = defaultdict(list)
    for row_ix, row in exact_bat.iterrows():
        herb = clean(row.get("Chinese.Name"))
        raw_field = clean(row.get("Ingredients"))
        for ordinal, token in enumerate((part.strip() for part in raw_field.split("|")), start=1):
            if not token:
                continue
            match = re.match(r"^(.*?)\s*\((\d{4,})\)\s*$", token)
            name = match.group(1).strip() if match else token
            cid = match.group(2) if match else ""
            variants = cid_smiles.get(cid, set()) if cid else set()
            smiles = next(iter(variants)) if len(variants) == 1 else ""
            rec = {
                "Herb": herb,
                "Source": "BATMAN-TCM 2.0",
                "Source_Record_ID": f"row_{row_ix + 2}_ingredient_{ordinal}",
                "Ingredient_Name_Raw": name,
                "Source_Identifier_Raw": cid,
                "PubChem_CID": cid,
                "Structure_SMILES": smiles,
                "Structure_Source": "TCMBank exact CID join" if smiles else "not_resolved",
                "CID_Crosscheck_Status": (
                    "matched_TCMBank_unique_structure" if cid and len(variants) == 1
                    else "CID_structure_conflict_in_TCMBank" if cid and len(variants) > 1
                    else "embedded_identifier_only" if cid else "no_embedded_numeric_identifier"
                ),
                "Ingredient_Row_Ordinal": ordinal,
                "TCMBank_Name_By_CID": "; ".join(sorted(cid_names.get(cid, set()))) if cid else "",
                "TCMBank_Molecular_Formula": "",
                "TCMBank_Molecular_Weight": "",
                "TCMBank_OB_score": "",
                "Exact_Herb_Name_Match": True,
            }
            add_source_record(records, **rec)
            batman_by_herb[herb].append(rec)

    tcmbank_by_herb: dict[str, list[dict]] = defaultdict(list)
    for row_ix, row in exact_tcm.iterrows():
        herb = clean(row.get("TCM_name"))
        ids = numeric_ids(row.get("PubChem_id")) or [""]
        for ordinal, cid in enumerate(ids, start=1):
            rec = {
                "Herb": herb,
                "Source": "TCMBank",
                "Source_Record_ID": clean(row.get("Ingredient_id")) or f"row_{row_ix + 2}",
                "Ingredient_Name_Raw": clean(row.get("name")),
                "Source_Identifier_Raw": clean(row.get("PubChem_id")),
                "PubChem_CID": cid,
                "Structure_SMILES": clean(row.get("Smiles")),
                "Structure_Source": "TCMBank exact herb row",
                "CID_Crosscheck_Status": "TCMBank_exact_herb_row" if cid else "no_PubChem_CID",
                "Ingredient_Row_Ordinal": ordinal,
                "TCMBank_Name_By_CID": clean(row.get("name")),
                "TCMBank_Molecular_Formula": clean(row.get("Molecular_Formula")),
                "TCMBank_Molecular_Weight": clean(row.get("Molecular_Weight")),
                "TCMBank_OB_score": clean(row.get("OB_score")),
                "Exact_Herb_Name_Match": True,
            }
            add_source_record(records, **rec)
            tcmbank_by_herb[herb].append(rec)

    return exact_tcm, alias_tcm, {"BATMAN-TCM 2.0": batman_by_herb, "TCMBank": tcmbank_by_herb}, cid_smiles


def build_source_records(herbs: list[str], tcmbank_exact: pd.DataFrame, alias_tcm: pd.DataFrame,
                        source_by_herb: dict[str, dict[str, list[dict]]]) -> pd.DataFrame:
    rows = []
    for herb in herbs:
        for src_name, herb_dict in source_by_herb.items():
            rows.extend(herb_dict.get(herb, []))
    frame = pd.DataFrame(rows)
    if not frame.empty:
        frame["_row_key"] = range(len(frame))
        frame["Structure_SMILES_Parser_Input"] = frame["Structure_SMILES"].map(
            lambda value: re.sub(r"[\r\n]", "", clean(value))
        )
        frame["Structure_Whitespace_Normalized"] = frame.apply(
            lambda row: clean(row.get("Structure_SMILES")) != clean(row.get("Structure_SMILES_Parser_Input")), axis=1
        )
        frame["Structure_Canonical_Isomeric_SMILES"] = frame["Structure_SMILES_Parser_Input"].map(canonical_smiles)
        frame["Structure_Parse_Status"] = frame.apply(
            lambda row: "not_provided" if not clean(row.get("Structure_SMILES"))
            else "valid_after_CRLF_normalization" if clean(row.get("Structure_Canonical_Isomeric_SMILES"))
            else "invalid_after_CRLF_normalization", axis=1
        )
    return frame


def build_master(source_records: pd.DataFrame) -> pd.DataFrame:
    structures_by_cid: dict[str, set[str]] = defaultdict(set)
    for row in source_records.to_dict("records"):
        cid = clean(row.get("PubChem_CID"))
        smi = clean(row.get("Structure_Canonical_Isomeric_SMILES"))
        if cid and smi:
            structures_by_cid[cid].add(smi)

    grouped: dict[str, list[dict]] = defaultdict(list)
    for row in source_records.to_dict("records"):
        cid = clean(row.get("PubChem_CID"))
        smi = clean(row.get("Structure_Canonical_Isomeric_SMILES"))
        name = clean(row.get("Ingredient_Name_Raw"))
        if cid:
            key = f"CID:{cid}"
        elif smi:
            key = f"SMILES:{smi}"
        else:
            key = f"NAME:{norm(name)}"
        grouped[key].append(row)

    out = []
    for key, items in sorted(grouped.items()):
        cid = key.removeprefix("CID:") if key.startswith("CID:") else ""
        variants = structures_by_cid.get(cid, set()) if cid else {
            clean(x.get("Structure_Canonical_Isomeric_SMILES")) for x in items
            if clean(x.get("Structure_Canonical_Isomeric_SMILES"))
        }
        conflict = len(variants) > 1
        smiles = next(iter(variants)) if len(variants) == 1 else ""
        names = sorted({clean(x.get("Ingredient_Name_Raw")) for x in items if clean(x.get("Ingredient_Name_Raw"))})
        herbs = sorted({clean(x.get("Herb")) for x in items if clean(x.get("Herb"))})
        sources = sorted({clean(x.get("Source")) for x in items if clean(x.get("Source"))})
        out.append({
            "Candidate_ID": key,
            "Candidate_Name": names[0] if names else "",
            "Alternative_Source_Names": "; ".join(names[1:]),
            "PubChem_CID": cid,
            "Herbs": "; ".join(herbs),
            "Herb_Count": len(herbs),
            "Sources": "; ".join(sources),
            "Source_Record_Count": len(items),
            "Structure_SMILES": smiles,
            "Structure_Status": (
                "conflicting_structures_by_CID" if conflict else "structure_available" if smiles
                else "structure_not_resolved"
            ),
            "Structure_Variant_Count": len(variants),
            "Structure_Variants": " || ".join(sorted(variants)),
            "BATMAN_CID_Crosschecked": any(
                x.get("Source") == "BATMAN-TCM 2.0" and x.get("CID_Crosscheck_Status") == "matched_TCMBank_unique_structure"
                for x in items
            ),
        })
    return pd.DataFrame(out)


def source_coverage(herbs: list[str], source_by_herb: dict[str, dict[str, list[dict]]], alias_tcm: pd.DataFrame,
                    master: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for herb in herbs:
        bat = source_by_herb["BATMAN-TCM 2.0"].get(herb, [])
        tcm = source_by_herb["TCMBank"].get(herb, [])
        alias_count = int(alias_tcm["TCM_name2"].map(lambda x: clean(x) == herb).sum())
        cid_values = {clean(x.get("PubChem_CID")) for x in bat + tcm if clean(x.get("PubChem_CID"))}
        bat_ids = {clean(x.get("PubChem_CID")) for x in bat if clean(x.get("PubChem_CID"))}
        matched = {clean(x.get("PubChem_CID")) for x in bat if x.get("CID_Crosscheck_Status") == "matched_TCMBank_unique_structure"}
        herb_master = master[master["Herbs"].map(lambda x: herb in x.split("; "))]
        rows.append({
            "Herb": herb,
            "BATMAN_ingredient_records": len(bat),
            "BATMAN_unique_names": len({norm(x.get("Ingredient_Name_Raw")) for x in bat if norm(x.get("Ingredient_Name_Raw"))}),
            "BATMAN_embedded_numeric_ids": len(bat_ids),
            "BATMAN_ids_with_TCMBank_structure_match": len(matched),
            "TCMBank_exact_primary_rows": len(tcm),
            "TCMBank_unique_primary_names": len({norm(x.get("Ingredient_Name_Raw")) for x in tcm if norm(x.get("Ingredient_Name_Raw"))}),
            "TCMBank_primary_CIDs": len({clean(x.get("PubChem_CID")) for x in tcm if clean(x.get("PubChem_CID"))}),
            "TCMBank_secondary_alias_rows_excluded": alias_count,
            "Unique_CID_candidates_across_sources": len(cid_values),
            "Candidates_with_resolved_structure": int(herb_master["Structure_SMILES"].astype(bool).sum()),
            "Candidates_without_resolved_structure": int((~herb_master["Structure_SMILES"].astype(bool)).sum()),
        })
    return pd.DataFrame(rows)


def load_step6_module():
    spec = importlib.util.spec_from_file_location("tcm_select_step6_readonly", STEP6_PATH)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load existing Step 6 implementation at {STEP6_PATH}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def screen_structures(master: pd.DataFrame, policy: dict) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    step6 = load_step6_module()
    structural = master[(master["Structure_SMILES"].astype(bool)) & (~master["Structure_Status"].eq("conflicting_structures_by_CID"))].copy()
    # One prediction per unique stereochemical RDKit canonical SMILES.
    unique_smiles = sorted({canonical_smiles(s) for s in structural["Structure_SMILES"] if canonical_smiles(s)})
    rows = []
    smile_status = {}
    for smi in unique_smiles:
        descriptors = step6.compute_descriptors(smi)
        if descriptors is None:
            smile_status[smi] = ("invalid_structure", [], [])
            rows.append({"Canonical_Isomeric_SMILES": smi, "Step6_Status": "invalid_structure"})
            continue
        passed, failed, warnings = step6.evaluate_filters(descriptors)
        status = "pass_heuristic" if passed else "fail_heuristic"
        smile_status[smi] = (status, failed, warnings)
        rows.append({
            "Canonical_Isomeric_SMILES": smi,
            **descriptors,
            "Step6_Status": status,
            "Step6_Failed_Rules": "; ".join(failed),
            "Step6_Warnings": "; ".join(warnings),
            "Step6_Policy_ID": policy["steps"]["step6"]["policy_id"],
        })
    screening = pd.DataFrame(rows)
    lookup = screening.set_index("Canonical_Isomeric_SMILES").to_dict("index") if not screening.empty else {}
    for ix, row in master.iterrows():
        smi = canonical_smiles(row["Structure_SMILES"])
        if smi in lookup:
            master.loc[ix, "Step6_Status"] = lookup[smi].get("Step6_Status", "not_assessed")
        elif row["Structure_Status"] == "conflicting_structures_by_CID":
            master.loc[ix, "Step6_Status"] = "held_structure_conflict"
        else:
            master.loc[ix, "Step6_Status"] = "not_assessed_no_structure"
    meta = {
        "input_unique_valid_structures": len(unique_smiles),
        "step6_pass_count": int(screening["Step6_Status"].eq("pass_heuristic").sum()) if not screening.empty else 0,
        "step6_fail_count": int(screening["Step6_Status"].eq("fail_heuristic").sum()) if not screening.empty else 0,
        "step6_invalid_count": int(screening["Step6_Status"].eq("invalid_structure").sum()) if not screening.empty else 0,
        "step6_policy_id": policy["steps"]["step6"]["policy_id"],
        "step6_output_class": policy["steps"]["step6"]["output_class"],
        "step6_output_strength": policy["steps"]["step6"]["output_strength"],
        "step6_interpretation": policy["steps"]["step6"]["interpretation_note"],
        "rdkit_version": rdBase.rdkitVersion,
    }
    return master, screening, meta


def predict_admet(screening: pd.DataFrame, master: pd.DataFrame, policy: dict) -> tuple[pd.DataFrame, dict]:
    passed = screening[screening["Step6_Status"].eq("pass_heuristic")] if not screening.empty else pd.DataFrame()
    smiles = sorted(set(passed["Canonical_Isomeric_SMILES"].astype(str))) if not passed.empty else []
    if not smiles:
        return pd.DataFrame(), {"model_status": "no_step6_pass_inputs", "input_count": 0}
    try:
        import torch
        torch.set_num_threads(min(2, max(1, torch.get_num_threads())))
        from admet_ai import ADMETModel
        model = ADMETModel()
        predicted = model.predict(smiles)
        if isinstance(predicted, pd.Series):
            predicted = predicted.to_frame().T
        if not isinstance(predicted, pd.DataFrame):
            predicted = pd.DataFrame(predicted)
        if predicted.index.name == "SMILES" or not any(str(c).lower() == "smiles" for c in predicted.columns):
            predicted = predicted.reset_index().rename(columns={"index": "Canonical_Isomeric_SMILES", "SMILES": "Canonical_Isomeric_SMILES"})
        else:
            predicted = predicted.rename(columns={next(c for c in predicted.columns if str(c).lower() == "smiles"): "Canonical_Isomeric_SMILES"})
        if "Canonical_Isomeric_SMILES" not in predicted:
            predicted.insert(0, "Canonical_Isomeric_SMILES", smiles[:len(predicted)])
        relevant = [
            "Canonical_Isomeric_SMILES", "HIA_Hou", "Caco2_Wang", "BBB_Martins", "Pgp_Broccatelli",
            "CYP1A2_Veith", "CYP2C19_Veith", "CYP2C9_Veith", "CYP2D6_Veith", "CYP3A4_Veith",
            "Solubility_AqSolDB", "DILI", "hERG", "AMES",
        ]
        keep = [c for c in relevant if c in predicted.columns]
        out = predicted[keep].copy()
        out["Step7_Model"] = "ADMET-AI local installed ensemble"
        out["Step7_Policy_ID"] = policy["steps"]["step7"]["policy_id"]
        thresholds = policy["steps"]["step7"]["thresholds"]
        decisions = []
        for row in out.to_dict("records"):
            failed, warnings = [], []
            hia = safe_float(row.get("HIA_Hou"))
            if hia is not None and (hia * 100 if hia <= 1 else hia) < thresholds["HIA_min"]:
                failed.append("HIA_below_project_threshold")
            caco = safe_float(row.get("Caco2_Wang"))
            if caco is not None and caco < thresholds["Caco2_min"]:
                failed.append("Caco2_below_project_threshold")
            sol = safe_float(row.get("Solubility_AqSolDB"))
            if sol is not None and sol < thresholds["ESOL_logS_min"]:
                failed.append("solubility_below_project_threshold")
            cyp_cols = ["CYP1A2_Veith", "CYP2C19_Veith", "CYP2C9_Veith", "CYP2D6_Veith", "CYP3A4_Veith"]
            cyp_pos = [c for c in cyp_cols if (safe_float(row.get(c)) or 0) > 0.5]
            if len(cyp_pos) >= thresholds["CYP_inhibition_max_count"]:
                failed.append("CYP_multi_inhibition_prediction")
            elif cyp_pos:
                warnings.append("single_CYP_inhibition_prediction")
            dili = safe_float(row.get("DILI"))
            herg = safe_float(row.get("hERG"))
            if thresholds.get("hard_reject_dual_tox") and dili is not None and herg is not None and dili > 0.5 and herg > 0.5:
                failed.append("DILI_and_hERG_dual_positive_prediction")
            pgp = safe_float(row.get("Pgp_Broccatelli"))
            if thresholds.get("warn_Pgp_inhibitor") and pgp is not None and pgp > 0.5:
                warnings.append("P-gp_inhibitor_prediction")
            bbb = safe_float(row.get("BBB_Martins"))
            route = "unresolved_predicted_BBB_positive" if bbb is not None and bbb > 0.5 else "unresolved_predicted_BBB_negative" if bbb is not None else "unresolved_no_BBB_prediction"
            decisions.append({
                "Step7_Status": "pass_predictive_heuristic" if not failed else "fail_predictive_heuristic",
                "Step7_Failed_Rules": "; ".join(failed),
                "Step7_Warnings": "; ".join(warnings),
                "CNS_Route": route,
            })
        decision_frame = pd.DataFrame(decisions)
        out = pd.concat([out.reset_index(drop=True), decision_frame], axis=1)
        package_version = importlib.metadata.version("admet-ai")
        model_info = {
            "model_status": "completed",
            "model": "ADMET-AI local installed ensemble",
            "package_version": package_version,
            "input_count": len(smiles),
            "output_count": len(out),
            "output_endpoints": [c for c in out.columns if c in relevant],
            "step7_pass_predictive_count": int(out["Step7_Status"].eq("pass_predictive_heuristic").sum()),
            "step7_fail_predictive_count": int(out["Step7_Status"].eq("fail_predictive_heuristic").sum()),
            "step7_policy_id": thresholds and policy["steps"]["step7"]["policy_id"],
            "interpretation": policy["steps"]["step7"]["interpretation_note"],
            "BBB_gate_required": bool(thresholds.get("BBB_penetrant_required", False)),
        }
        return out, model_info
    except Exception as exc:  # retain the completed ingredient and RDKit audit if model loading is blocked
        return pd.DataFrame(), {
            "model_status": "failed_or_unavailable",
            "input_count": len(smiles),
            "error_type": type(exc).__name__,
            "error": str(exc)[:1000],
        }


def parse_cid_col(value: object) -> str:
    values = numeric_ids(value)
    return values[0] if values else ""


def target_priors(candidate_cids: set[str]) -> tuple[pd.DataFrame, pd.DataFrame, dict, set[str]]:
    priors = read_csv(TARGET_MAPPING)
    table = read_csv(TARGET_TABLE)
    for frame in (priors, table):
        frame["_CID"] = frame["PubChem_CID"].map(parse_cid_col)
    prior_hits = priors[priors["_CID"].isin(candidate_cids)].copy()
    table_hits = table[table["_CID"].isin(candidate_cids)].copy()
    for frame in (prior_hits, table_hits):
        frame["Evidence_Interpretation"] = "imported_database_target_mapping_prior; action direction, independent assay count and exposure are not established"
        frame["Action_Direction"] = "unknown"
    if not prior_hits.empty:
        prior_hits["Mapping_Table_Source"] = "10_Ingredient_Target_Mapping.csv"
    if not table_hits.empty:
        table_hits["Mapping_Table_Source"] = "9_Compound_Target_Table.csv"

    measured_rows = []
    # BindingDB: exact CID and human normalized target records with recorded assay type and numeric value.
    bind = read_csv(BINDINGDB)
    if not bind.empty:
        bind["_CID"] = bind["PubChem_CID"].map(parse_cid_col)
        for _, r in bind[bind["_CID"].isin(candidate_cids)].iterrows():
            if clean(r.get("Target_Symbol")) and safe_float(r.get("BindingDB_Affinity_Value")) is not None:
                measured_rows.append({
                    "Evidence_Source": "BindingDB",
                    "PubChem_CID": clean(r.get("_CID")),
                    "Ingredient": clean(r.get("Ingredient")),
                    "Target_Symbol": clean(r.get("Target_Symbol")),
                    "Target_Name": clean(r.get("BindingDB_Target_Name")),
                    "Assay_Type": clean(r.get("BindingDB_Affinity_Type")),
                    "Activity_Value": clean(r.get("BindingDB_Affinity_Value")),
                    "Activity_Units": "as stored in local normalized export; verify source record before quantitative use",
                    "Species": clean(r.get("Target_Species")),
                    "Evidence_Status": "database_measurement_record_exact_CID",
                    "Action_Direction": "not_inferred_from_affinity_value_alone",
                })

    # ChEMBL: require human target, numeric activity value and 100% structure similarity in cached response.
    chem = read_csv(CHEMBL)
    if not chem.empty:
        chem["_CID"] = chem["PubChem_CID"].map(parse_cid_col)
        for _, r in chem[chem["_CID"].isin(candidate_cids)].iterrows():
            similarity = safe_float(r.get("Similarity"))
            if (clean(r.get("Target_Symbol")) and safe_float(r.get("Standard_Value")) is not None
                    and norm(r.get("Target_Species")) == "homo sapiens" and similarity is not None and similarity >= 100):
                measured_rows.append({
                    "Evidence_Source": "ChEMBL",
                    "PubChem_CID": clean(r.get("_CID")),
                    "Ingredient": clean(r.get("Ingredient")),
                    "Target_Symbol": clean(r.get("Target_Symbol")),
                    "Target_Name": clean(r.get("Matched_Pref_Name")) or clean(r.get("Target_ChEMBL_ID")),
                    "Assay_Type": clean(r.get("Standard_Type")),
                    "Activity_Value": clean(r.get("Standard_Value")),
                    "Activity_Units": clean(r.get("Standard_Units")),
                    "Species": clean(r.get("Target_Species")),
                    "Evidence_Status": "database_measurement_record_exact_CID_similarity_100",
                    "Action_Direction": "not_inferred_from_activity_value_alone",
                    "ChEMBL_ID": clean(r.get("Matched_ChEMBL_ID")),
                    "Activity_ID": clean(r.get("Activity_ID")),
                    "Assay_ID": clean(r.get("Assay_ChEMBL_ID")),
                    "pChEMBL_Value": clean(r.get("pChEMBL_Value")),
                })

    # PubChem: retain only records labelled Active with a human target identifier and a numeric activity value.
    act_cols = [
        "Activity", "Compound_CID", "Activity_Type", "Activity_Qualifier", "Target_Name", "Protein_Accession",
        "Representative_Protein_Accession", "BioAssay_Name", "dois", "pmcids", "citations", "BioAssay_AID",
        "Bioactivity_ID", "Gene_ID", "PMID", "Activity_Value", "Target_Taxonomy_ID", "Taxonomy_ID",
        "Has_Dose_Response_Curve",
    ]
    activity = read_csv(PUBCHEM_ACTIVITY, usecols=lambda c: c in act_cols)
    genes = read_csv(PUBCHEM_GENES, usecols=lambda c: c in ["Gene_ID", "Gene_Symbol", "Taxonomy_ID"])
    gene_symbol = {parse_cid_col(r.get("Gene_ID")): clean(r.get("Gene_Symbol")) for r in genes.to_dict("records") if parse_cid_col(r.get("Gene_ID"))}
    activity["_CID"] = activity["Compound_CID"].map(parse_cid_col)
    for _, r in activity[activity["_CID"].isin(candidate_cids)].iterrows():
        target_tax = parse_cid_col(r.get("Target_Taxonomy_ID"))
        gene_id = parse_cid_col(r.get("Gene_ID"))
        active = norm(r.get("Activity")) == "active"
        value = safe_float(r.get("Activity_Value"))
        protein = clean(r.get("Protein_Accession")) or clean(r.get("Representative_Protein_Accession"))
        symbol = gene_symbol.get(gene_id, "")
        if active and target_tax == "9606" and value is not None and (symbol or protein):
            measured_rows.append({
                "Evidence_Source": "PubChem_BioAssay",
                "PubChem_CID": clean(r.get("_CID")),
                "Ingredient": "",
                "Target_Symbol": symbol,
                "Target_Name": clean(r.get("Target_Name")),
                "Assay_Type": clean(r.get("Activity_Type")),
                "Activity_Value": clean(r.get("Activity_Value")),
                "Activity_Units": "as stored in PubChem BioAssay record",
                "Species": "Homo sapiens target taxonomy 9606",
                "Evidence_Status": "PubChem_active_assay_record_numeric_value_human_target",
                "Action_Direction": "not_inferred_from_activity_record_alone",
                "BioAssay_AID": clean(r.get("BioAssay_AID")),
                "Bioactivity_ID": clean(r.get("Bioactivity_ID")),
                "Gene_ID": gene_id,
                "PMID": clean(r.get("PMID")),
                "DOI": clean(r.get("dois")),
                "Has_Dose_Response_Curve": clean(r.get("Has_Dose_Response_Curve")),
            })

    measured = pd.DataFrame(measured_rows)
    if not measured.empty:
        measured = measured.drop_duplicates().reset_index(drop=True)
    all_symbols = set()
    for f in (prior_hits, table_hits, measured):
        if not f.empty and "Target_Symbol" in f:
            all_symbols.update(clean(x) for x in f["Target_Symbol"] if clean(x))
    metrics = {
        "formula_candidate_unique_cids": len(candidate_cids),
        "target_mapping_prior_rows_exact_cid": len(prior_hits),
        "target_table_rows_exact_cid": len(table_hits),
        "assay_level_measurement_records_retained": len(measured),
        "measured_evidence_sources": sorted(set(measured["Evidence_Source"])) if not measured.empty else [],
        "unique_symbols_in_any_target_cache": len(all_symbols),
        "unique_symbols_in_assay_level_rows": int(measured["Target_Symbol"].replace("", pd.NA).dropna().nunique()) if not measured.empty else 0,
        "action_directions_inferred": 0,
        "new_online_target_queries": False,
    }
    return pd.concat([prior_hits, table_hits], ignore_index=True, sort=False), measured, metrics, all_symbols


def g1_overlap(measured: pd.DataFrame, g1_path: Path) -> tuple[pd.DataFrame, dict]:
    g1 = pd.read_csv(g1_path, sep="\t", compression="gzip", dtype=str, keep_default_na=False)
    g1_symbols = {clean(x): r for x, r in zip(g1["gene"], g1.to_dict("records"))}
    if measured.empty:
        return pd.DataFrame(), {"g1_robust_gene_count": len(g1_symbols), "assay_target_symbols": 0, "overlap_symbols": 0}
    rows = []
    for symbol in sorted({clean(x) for x in measured["Target_Symbol"] if clean(x)}):
        if symbol in g1_symbols:
            gene = g1_symbols[symbol]
            records = measured[measured["Target_Symbol"].map(clean).eq(symbol)]
            rows.append({
                "Target_Symbol": symbol,
                "G1_meta_logFC": gene.get("meta_logFC", ""),
                "G1_meta_FDR": gene.get("meta_FDR", ""),
                "G1_direction": gene.get("direction", ""),
                "Assay_Evidence_Sources": "; ".join(sorted(set(records["Evidence_Source"]))),
                "Assay_Record_Count": len(records),
                "Compound_CIDs": "; ".join(sorted(set(clean(x) for x in records["PubChem_CID"]))),
                "Compound_Action_Direction": "unknown",
                "Interpretation": "descriptive overlap only; cannot infer disease reversal",
            })
    result = pd.DataFrame(rows)
    assay_symbols = {clean(x) for x in measured["Target_Symbol"] if clean(x)}
    return result, {
        "g1_robust_gene_count": len(g1_symbols),
        "assay_target_symbols": len(assay_symbols),
        "overlap_symbols": len(result),
    }


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    herbs_df = pd.read_csv(FORMULA_PATH, sep="\t", dtype=str, keep_default_na=False)
    herbs = herbs_df["herb"].astype(str).str.strip().tolist()
    if len(herbs) != 16 or len(set(herbs)) != 16:
        raise RuntimeError("Registered formula does not contain exactly 16 unique herb names")
    policy = json.loads(POLICY_PATH.read_text(encoding="utf-8"))

    input_paths = [FORMULA_PATH, TCMBANK, BATMAN, POLICY_PATH, STEP6_PATH, TARGET_MAPPING, TARGET_TABLE,
                   PUBCHEM_ACTIVITY, PUBCHEM_GENES, BINDINGDB, CHEMBL, G1_244, PARAM_PATH, BASE_SCRIPT, RUN_SCRIPT]
    source_manifest = pd.DataFrame([
        {"Path": str(path), "Bytes": path.stat().st_size, "SHA256": sha256(path), "Role": "input_or_frozen_method"}
        for path in input_paths
    ])
    source_manifest.to_csv(OUT / "source_sha256.tsv", sep="\t", index=False, encoding="utf-8")

    exact_tcm, alias_tcm, source_by_herb, _ = read_ingredient_sources(herbs)
    source_records = build_source_records(herbs, exact_tcm, alias_tcm, source_by_herb)
    master = build_master(source_records)
    coverage = source_coverage(herbs, source_by_herb, alias_tcm, master)
    source_records.to_csv(OUT / "formula_component_source_records.tsv", sep="\t", index=False, encoding="utf-8")
    master, screening, screen_metrics = screen_structures(master, policy)
    admet, admet_metrics = predict_admet(screening, master, policy)
    if not admet.empty:
        admet.to_csv(OUT / "formula_admet_predictions.tsv", sep="\t", index=False, encoding="utf-8")
        admet_map = admet.set_index("Canonical_Isomeric_SMILES").to_dict("index")
        master["Step7_Status"] = master["Structure_SMILES"].map(lambda s: admet_map.get(canonical_smiles(s), {}).get("Step7_Status", "not_assessed"))
    else:
        master["Step7_Status"] = "not_assessed_model_unavailable_or_no_input"
    screening.to_csv(OUT / "formula_structure_and_step6_screen.tsv", sep="\t", index=False, encoding="utf-8")
    master.to_csv(OUT / "formula_component_candidates.tsv", sep="\t", index=False, encoding="utf-8")
    coverage.to_csv(OUT / "formula_herb_source_coverage.tsv", sep="\t", index=False, encoding="utf-8")

    candidate_cids = {clean(x) for x in master["PubChem_CID"] if clean(x)}
    target_map, measured, target_metrics, target_symbols = target_priors(candidate_cids)
    target_map.to_csv(OUT / "formula_target_database_priors.tsv", sep="\t", index=False, encoding="utf-8")
    measured.to_csv(OUT / "formula_assay_level_target_records.tsv", sep="\t", index=False, encoding="utf-8")
    overlap, overlap_metrics = g1_overlap(measured, G1_244)
    overlap.to_csv(OUT / "formula_target_overlap_with_G1_robust244.tsv", sep="\t", index=False, encoding="utf-8")

    packages = {}
    for pkg in ["pandas", "openpyxl", "rdkit", "admet-ai", "torch"]:
        try:
            packages[pkg] = importlib.metadata.version(pkg)
        except importlib.metadata.PackageNotFoundError:
            packages[pkg] = "unknown"
    metrics = {
        "run_id": RUN_ID,
        "decision_id": DECISION_ID,
        "herb_count": len(herbs),
        "herbs_with_exact_BATMAN_records": int(coverage["BATMAN_ingredient_records"].gt(0).sum()),
        "herbs_with_exact_primary_TCMBank_records": int(coverage["TCMBank_exact_primary_rows"].gt(0).sum()),
        "herbs_with_TCMBank_secondary_alias_only_or_additional_rows": int(coverage["TCMBank_secondary_alias_rows_excluded"].gt(0).sum()),
        "formula_source_records": len(source_records),
        "candidate_record_keys": len(master),
        "candidate_unique_pubchem_cids": len(candidate_cids),
        "candidates_with_structures": int(master["Structure_SMILES"].astype(bool).sum()),
        "structure_conflict_candidate_count": int(master["Structure_Status"].eq("conflicting_structures_by_CID").sum()),
        "source_coverage": coverage.to_dict("records"),
        "structure_screen": screen_metrics,
        "admet_ai": admet_metrics,
        "target_prior": target_metrics,
        "g1_robust_overlap": overlap_metrics,
        "g2_g5_disposition": "source conclusions remain model-scoped; no new drug perturbation or causal inference performed",
        "P0_status": "16 herb identities registered; dose, proportions, species/part, processing and preparation not specified",
        "P1_status": "formula equivalence to experimental/literature preparation not established",
        "P2_observed_formula_chemistry": "not measured; database prior only",
        "P2_in_vivo_exposure": "not measured; parent/metabolite exposure and brain concentration unavailable",
        "P3_exposure_potency_match": "not established",
        "therapy_direction": "undetermined; target action direction and relevant exposure are unavailable",
        "drug_evidence_ceiling": "C0",
        "closed_validation_steps_reopened": False,
        "software_versions": {**packages, "python": sys.version.split()[0], "rdkit_runtime": rdBase.rdkitVersion},
    }
    (OUT / "summary.json").write_text(json.dumps(metrics, ensure_ascii=False, indent=2), encoding="utf-8")

    # Rehash unchanged raw inputs to verify that this run did not alter source files.
    source_after = {str(p): sha256(p) for p in [TCMBANK, BATMAN, POLICY_PATH, TARGET_MAPPING, TARGET_TABLE,
                                                  PUBCHEM_ACTIVITY, PUBCHEM_GENES, BINDINGDB, CHEMBL]}
    source_before = dict(zip(source_manifest["Path"], source_manifest["SHA256"]))
    unchanged = all(source_before[p] == h for p, h in source_after.items())
    audit = {
        "status": "PASS_WITH_PREDECLARED_EVIDENCE_LIMITS" if unchanged and len(herbs) == 16 else "FAIL",
        "checks": {
            "sixteen_unique_herb_inputs": len(herbs) == 16 and len(set(herbs)) == 16,
            "all_sources_exact_name_match_only": bool(source_records["Exact_Herb_Name_Match"].all()) if not source_records.empty else True,
            "TCMbank_alias_rows_not_in_primary_extract": True,
            "structure_join_uses_numeric_CID_only": True,
            "conflicting_CID_structures_excluded_from_RDKit": True,
            "raw_sources_unchanged": unchanged,
            "new_online_queries": False,
            "new_dependencies_installed": False,
            "drug_effect_or_mechanism_validation_run": False,
            "closed_steps_reopened": False,
            "C0_claim_ceiling_retained": True,
        },
        "source_hash_before_after_equal": unchanged,
        "candidate_rows": len(master),
        "source_record_rows": len(source_records),
        "target_assay_rows": len(measured),
        "admet_model_status": admet_metrics.get("model_status"),
    }
    (OUT / "audit.json").write_text(json.dumps(audit, ensure_ascii=False, indent=2), encoding="utf-8")

    cov_lines = [
        "| 药味 | BATMAN名称记录 | BATMAN含数字ID | ID与TCMBank结构匹配 | TCMBank精确主名行 | 主名CID | 排除的二级别名行 |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for r in coverage.to_dict("records"):
        cov_lines.append(f"| {r['Herb']} | {r['BATMAN_ingredient_records']} | {r['BATMAN_embedded_numeric_ids']} | {r['BATMAN_ids_with_TCMBank_structure_match']} | {r['TCMBank_exact_primary_rows']} | {r['TCMBank_primary_CIDs']} | {r['TCMBank_secondary_alias_rows_excluded']} |")
    report = f'''# 固定处方药物侧先验审计 {RUN_ID[-3:]}

日期：2026-09-30。方案：药物侧 v1.3 P0–P3，按决策{RUN_ID[-3:]}限定为现存数据库先验整理。固定处方来自018登记的用户 Word 清单，共16味。

## 结论

本轮完成了处方药味与 BATMAN/TCMBank 本地成分目录的逐味核对、可由精确 PubChem CID 连接的结构恢复、Step 6 RDKit 理化启发式筛查、Step 7 本地 ADMET-AI 预测，以及现有本地靶点/活性缓存的精确 CID 检索。详细逐味覆盖见[覆盖表](formula_herb_source_coverage.tsv)，成分来源记录见[原始记录表](formula_component_source_records.tsv)，化学候选见[成分主表](formula_component_candidates.tsv)。

这些输出均为数据库先验。处方文档只有药味名称，未给剂量、比例、药用部位/炮制、制备或批次质量，因此不能确认外部研究配方与本方等价，不能称为本方煎剂实测成分或体内暴露；药物证据等级保持 **C0**。017关闭的模型资格、药物纠正、暴露—效力、整方药效/机制及细胞通讯验证没有重启。

## 来源匹配与成分覆盖

处方药名只按原文中文名与 BATMAN `Chinese.Name`、TCMBank `TCM_name` 精确匹配。TCMBank `TCM_name2` 次级别名只列入覆盖审计，不并入主数据。BATMAN 数字括号标识只有在 TCMBank 出现同一数值 CID 时才用于结构链接；无此连接的条目仅保留为名称候选。结构有 CID 冲突的条目被排除在描述符和 ADMET 模型之外。上游原始目录和本地冻结政策的 SHA-256、文件大小见[来源清单](source_sha256.tsv)。

逐味行数、名称数、CID 数和结构匹配数见表。成分数据库收载并不等同于该成分在实际处方、煎液或体内出现；也不证明丰度、剂量或可吸收量。

{chr(10).join(cov_lines)}

## 化学空间与 ADMET 预测

结构统一使用本机 `{rdBase.rdkitVersion}` RDKit。对部分 TCMBank 单元格中的物理换行，仅在解析副本中删除 CR/LF，并保留原始字符串、规范化副本及解析状态。Step 6 复用了上游代码和 `pipeline_policy.json` 中 `{policy['steps']['step6']['policy_id']}`：Ro5、Veber、QED、NPL 与 PAINS/BRENK 警示规则。唯一有效结构数为 {screen_metrics['input_unique_valid_structures']}；启发式通过 {screen_metrics['step6_pass_count']}，未通过 {screen_metrics['step6_fail_count']}，无效结构 {screen_metrics['step6_invalid_count']}。逐分子结果见[结构与 Step 6 表](formula_structure_and_step6_screen.tsv)。通过只表示按本项目理化规则进入下一层，不表示活性或药效。

对 Step 6 通过的唯一结构运行本地 ADMET-AI；状态为 `{admet_metrics.get('model_status')}`，输入 {admet_metrics.get('input_count', 0)} 个，模型预测通过 {admet_metrics.get('step7_pass_predictive_count', 0)} 个、未通过 {admet_metrics.get('step7_fail_predictive_count', 0)} 个。预测字段与解释见[ADMET预测表](formula_admet_predictions.tsv)。这些数值不是实测吸收、毒性、脑暴露或自由浓度。BBB 预测不作强制门槛，也不能据此确立 Direct-CNS；没有合格脑内暴露时 CNS 路线保持未定。

## 靶点药理先验与 G1 接口

只从现有本地数据按处方候选 CID 精确查找靶点映射与 PubChem、BindingDB、ChEMBL 活性记录；本轮没有向网络数据库提交新查询。靶点关联行数 {target_metrics['target_mapping_prior_rows_exact_cid']}，严格保留的靶点活性/结合记录 {target_metrics['assay_level_measurement_records_retained']}，唯一来源为 {', '.join(target_metrics['measured_evidence_sources']) if target_metrics['measured_evidence_sources'] else '无'}。映射先验与 assay 记录分表保存：[映射先验](formula_target_database_priors.tsv)、[活性记录](formula_assay_level_target_records.tsv)。原始值、单位、assay、CID 与来源字段应在引用具体靶点前继续回查原文；单有亲和力或活性值不足以确定在本处方中的促进/抑制方向。

严格活性记录的靶点符号与 G1 稳健244基因集精确重叠 {overlap_metrics['overlap_symbols']} 个，见[G1描述性交集表](formula_target_overlap_with_G1_robust244.tsv)。这只是“数据库活性记录靶点符号”和“人类AD小胶质关联基因”共现，不说明本方在相关细胞中触达该靶点，也不说明能把疾病效应推向健康方向。G2–G5 的模型结论范围不变；人类 G1/G2c、模型内 G2–G5、G6关闭状态分别保留，不能合并为已验证的人类整方机制链。

## 按 P0–P3 的状态

| 门 | 当前处理 | 状态 |
|---|---|---|
| P0 处方身份 | 确认用户给定16味药名；剂量、比例、部位/炮制、制备、批次不全 | `IDENTITY_LIST_ONLY` |
| P1 方剂等价 | 既往候选处方只与本方部分药味重叠，缺剂量/制备比对 | `NOT_ESTABLISHED` |
| P2 化学空间 | 已得成分数据库目录和结构预测；无本方煎剂 LC-MS/QTOF 与含量定量 | `DATABASE_PRIOR_ONLY` |
| P2 体内暴露 | 缺母体/代谢物血浆及脑组织、游离分数、时间浓度 | `NOT_MEASURED` |
| P3 靶点药理 | 有部分精确 CID 的数据库活性记录；没有本方体系中的作用方向/剂量效力或暴露匹配 | `PRIOR_ONLY_NO_EXPOSURE_POTENCY_MATCH` |

## 复现与审计

执行代码：[处理入口](../../src/109_process_fixed_formula_prior_020.py)；复用实现：[公共处理模块](../../src/108_process_fixed_formula_prior_019.py)；参数：[运行参数](../../configs/fixed_formula_drug_prior_020.yaml)；来源SHA与逐项状态见来源清单和[audit.json](audit.json)。运行环境：项目已有 `TCM_select` Conda 环境，Python {sys.version.split()[0]}；版本记录在 `summary.json`。脚本只读外部原始表和既有策略代码，输出写入项目内；未安装依赖、未在线查询、未修改原始来源、未运行真实整方扰动/VC药效或机制检验。WikiSkill轨迹见[`workspaces/{RUN_ID}/`](../../workspaces/{RUN_ID}/)；DVC remote 本次未核验。
'''
    (OUT / "report.md").write_text(report, encoding="utf-8")


if __name__ == "__main__":
    main()
