#!/usr/bin/env python3
"""Create a compact, auditable table from an NCBI GEO esummary response."""
from __future__ import annotations

import csv
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "data" / "interim" / "G5_dataset_search" / "ncbi_gds_esummary_50.json"
OUTPUT = ROOT / "data" / "interim" / "G5_dataset_search" / "ncbi_gds_candidates_50.tsv"


def as_text(value: object) -> str:
    if value is None:
        return ""
    if isinstance(value, list):
        values: list[str] = []
        for item in value:
            if isinstance(item, dict):
                values.append(str(item.get("name") or item.get("acc") or item))
            else:
                values.append(str(item))
        return " | ".join(values)
    if isinstance(value, dict):
        return json.dumps(value, ensure_ascii=False, sort_keys=True)
    return str(value)


def main() -> None:
    payload = json.loads(INPUT.read_text(encoding="utf-8"))
    result = payload["result"]
    rows: list[dict[str, str]] = []
    for uid in result.get("uids", []):
        record = result[str(uid)]
        rows.append(
            {
                "uid": str(uid),
                "accession": as_text(record.get("accession")),
                "title": as_text(record.get("title")),
                "summary": as_text(record.get("summary")),
                "taxon": as_text(record.get("taxon")),
                "entrytype": as_text(record.get("entrytype")),
                "gds_type": as_text(record.get("gdsType") or record.get("gdstype")),
                "n_samples": as_text(record.get("n_samples") or record.get("nsamples")),
                "pubmedids": as_text(record.get("pubmedids")),
                "suppfile": as_text(record.get("suppfile")),
                "ftp_link": as_text(record.get("ftplink")),
                "sample_titles": as_text(record.get("samples")),
            }
        )
    fields = list(rows[0]) if rows else []
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(json.dumps({"records": len(rows), "output": str(OUTPUT)}, ensure_ascii=False))


if __name__ == "__main__":
    main()
