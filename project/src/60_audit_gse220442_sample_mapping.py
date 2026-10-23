#!/usr/bin/env python3
"""Resolve GEO sample accessions to the GSE220442 spot-level donor IDs."""

from __future__ import annotations

import gzip
import json
import re
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RAW_META = ROOT / "data/raw/G6_source_metadata_2026-09-23"
OUT = ROOT / "results/2026-09-24-disease-side-v3.5-g6-input-audit"
REG = OUT / "registries"
TABLES = OUT / "tables"


def parse_soft(path: Path) -> list[dict[str, str]]:
    records: list[dict[str, str]] = []
    current: dict[str, str] | None = None
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.rstrip("\r\n")
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


soft = pd.DataFrame(parse_soft(RAW_META / "GSE220442_family.soft.gz")).fillna("")
metadata = pd.read_csv(RAW_META / "GSE220442_metadata.csv.gz", compression="gzip", dtype=str).fillna("")
filelist = (RAW_META / "GSE220442_filelist.txt")
if filelist.exists():
    lines = filelist.read_text(encoding="utf-8", errors="replace").splitlines()
    files = []
    for line in lines[1:]:
        values = line.split("\t")
        if len(values) >= 4 and values[0] == "File":
            files.append({"GEO_name": values[1], "bytes": values[3], "type": values[4] if len(values) > 4 else ""})
    file_table = pd.DataFrame(files)
else:
    file_table = pd.DataFrame(columns=["GEO_name", "bytes", "type"])

donors = (
    metadata.groupby("patientID", dropna=False, sort=True)
    .agg(source_diagnosis=("category", lambda x: ";".join(sorted(set(x)))), n_spots=("patientID", "size"))
    .reset_index()
)
gsm_rows = []
for row in soft.to_dict(orient="records"):
    attributes = " | ".join(
        f"{key}={value}" for key, value in row.items()
        if key != "GSM" and value and ("title" in key.lower() or "characteristic" in key.lower() or "source" in key.lower() or "supplement" in key.lower() or "relation" in key.lower())
    )
    title = row.get("title", "")
    matches = []
    for patient in donors["patientID"].astype(str):
        if patient and patient.lower() in (title + " " + attributes).lower():
            matches.append(patient)
    gsm_rows.append({
        "GSM": row.get("GSM", ""),
        "title": title,
        "source_name_ch1": row.get("source_name_ch1", ""),
        "characteristics_ch1": row.get("characteristics_ch1", ""),
        "supplementary_file": row.get("supplementary_file", ""),
        "relations": row.get("relation", ""),
        "metadata_fields_for_audit": attributes,
        "direct_patientID_match_n": len(matches),
        "direct_patientID_match": ";".join(matches),
        "mapping_status": "unique_exact_patientID_in_GEO_sample_metadata" if len(matches) == 1 else "unresolved_no_or_ambiguous_patientID_in_GEO_sample_metadata",
    })

gsm = pd.DataFrame(gsm_rows)
gsm.to_csv(REG / "GSE220442_GSM_source_sample_manifest.tsv", sep="\t", index=False)
donors.to_csv(REG / "GSE220442_patientID_donor_manifest.tsv", sep="\t", index=False)
file_table.to_csv(REG / "GSE220442_GEO_file_manifest.tsv", sep="\t", index=False)

summary = {
    "GEO_sample_records": int(len(gsm)),
    "spot_metadata_patientIDs": int(donors["patientID"].nunique()),
    "unique_GSM_to_patientID_matches": int((gsm["mapping_status"] == "unique_exact_patientID_in_GEO_sample_metadata").sum()),
    "unresolved_GSM_records": int((gsm["mapping_status"] != "unique_exact_patientID_in_GEO_sample_metadata").sum()),
    "GEO_supplementary_files": int(len(file_table)),
    "source_rows": gsm.to_dict(orient="records"),
    "G6_effects_used": False,
}
(OUT / "GSE220442_GSM_sample_mapping_audit.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
print(json.dumps({key: value for key, value in summary.items() if key != "source_rows"}, ensure_ascii=False))
