#!/usr/bin/env python3
"""Query iLINCS/LINCS for exact PubChem-CID matched compound signatures.

Names are used only as a search handle.  A record is accepted as an exact
candidate perturbation only when the returned iLINCS pubChemID equals the
local PubChem CID.  Fuzzy name matches are retained in the audit but never
promoted to perturbation evidence.
"""
from __future__ import annotations

import argparse
import json
import re
import time
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import requests


API = "https://www.ilincs.org/api/SignatureMeta/findTermWithSynonyms"


def norm_cid(value: object) -> str:
    if pd.isna(value):
        return ""
    try:
        return str(int(float(value)))
    except (TypeError, ValueError):
        return str(value).strip()


def norm_name(value: object) -> str:
    return re.sub(r"[^a-z0-9]", "", str(value).lower())


def run(root: Path) -> None:
    master = pd.read_csv(root / "results" / "tables" / "drug_compound_master.tsv", sep="\t")
    unique = master[["Herb", "Ingredient", "PubChem_CID", "InChIKey"]].copy()
    unique["cid"] = unique["PubChem_CID"].map(norm_cid)
    unique = unique[unique["cid"].astype(bool)].drop_duplicates(subset=["cid", "Ingredient"])
    raw_dir = root / "data" / "raw" / "DrugSide_iLINCS"
    interim = root / "data" / "interim"
    raw_dir.mkdir(parents=True, exist_ok=True)
    interim.mkdir(parents=True, exist_ok=True)
    raw_path = raw_dir / "metadata_queries.jsonl"
    audit_rows: list[dict] = []
    exact_rows: list[dict] = []
    with raw_path.open("w", encoding="utf-8") as raw_handle:
        for record in unique.to_dict(orient="records"):
            term = str(record["Ingredient"])
            payload = {"query_time_utc": datetime.now(timezone.utc).isoformat(), "term": term, "requested_cid": record["cid"]}
            try:
                response = requests.get(API, params={"term": term, "library": "LIB_5"}, timeout=90)
                payload["http_status"] = response.status_code
                data = response.json().get("data", []) if response.ok else []
                payload["response_count"] = len(data)
                payload["response"] = data
                exact = [row for row in data if norm_cid(row.get("pubChemID")) == record["cid"]]
                name_term = norm_name(term)
                name_matches = [row for row in data if name_term and (name_term == norm_name(row.get("compound", "")) or name_term == norm_name(row.get("treatment", "")))]
                audit_rows.append({
                    "ingredient": term,
                    "PubChem_CID": record["cid"],
                    "http_status": response.status_code,
                    "response_count": len(data),
                    "exact_cid_match_count": len(exact),
                    "exact_signature_count": len({str(row.get("signatureid")) for row in exact if row.get("signatureid")}),
                    "name_match_count": len(name_matches),
                    "error": "",
                })
                for row in exact:
                    out = {
                        "Herb": record["Herb"],
                        "Ingredient": term,
                        "PubChem_CID": record["cid"],
                        "InChIKey": record["InChIKey"],
                        "signatureid": row.get("signatureid"),
                        "lincsSigID": row.get("lincsSigID"),
                        "cellline": row.get("cellline"),
                        "tissue": row.get("tissue"),
                        "concentration": row.get("concentration"),
                        "time": row.get("time"),
                        "compound": row.get("compound"),
                        "treatment": row.get("treatment"),
                        "pubChemID_returned": norm_cid(row.get("pubChemID")),
                        "perturbagenID": row.get("perturbagenID"),
                        "GeneTargets": row.get("GeneTargets"),
                        "integratedMoas": row.get("integratedMoas"),
                        "is_exemplar": row.get("is_exemplar"),
                        "platform": row.get("platform"),
                        "evidence_status": "exact_pubchem_cid_match_iLINCS_metadata",
                    }
                    exact_rows.append(out)
            except Exception as exc:
                payload["http_status"] = None
                payload["response_count"] = 0
                payload["error"] = f"{type(exc).__name__}: {exc}"
                audit_rows.append({"ingredient": term, "PubChem_CID": record["cid"], "http_status": None, "response_count": 0, "exact_cid_match_count": 0, "exact_signature_count": 0, "name_match_count": 0, "error": payload["error"]})
            raw_handle.write(json.dumps(payload, ensure_ascii=False) + "\n")
            raw_handle.flush()
            time.sleep(0.2)
    audit = pd.DataFrame(audit_rows).sort_values(["PubChem_CID", "ingredient"])
    exact = pd.DataFrame(exact_rows)
    audit.to_csv(interim / "ilincs_candidate_query_audit.tsv", sep="\t", index=False)
    if exact.empty:
        exact = pd.DataFrame(columns=["Herb", "Ingredient", "PubChem_CID", "InChIKey", "signatureid", "lincsSigID", "cellline", "tissue", "concentration", "time", "compound", "treatment", "pubChemID_returned", "perturbagenID", "GeneTargets", "integratedMoas", "is_exemplar", "platform", "evidence_status"])
    exact.to_csv(interim / "ilincs_exact_candidate_signatures.tsv", sep="\t", index=False)
    summary = {
        "query_records": int(unique.shape[0]),
        "http_success_records": int((audit["http_status"] == 200).sum()),
        "exact_cid_matched_compound_records": int(audit["exact_cid_match_count"].gt(0).sum()),
        "exact_signature_rows": int(exact.shape[0]),
        "unique_matched_cids": int(exact["PubChem_CID"].nunique()) if not exact.empty else 0,
        "fuzzy_name_matches_not_promoted": int((audit["name_match_count"] > 0).sum()),
        "source": API,
        "matching_rule": "returned pubChemID equals local PubChem CID; name-only matches are excluded",
    }
    (interim / "ilincs_candidate_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    args = parser.parse_args()
    run(args.root.resolve())
