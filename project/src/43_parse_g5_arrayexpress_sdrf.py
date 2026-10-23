#!/usr/bin/env python3
"""Parse ArrayExpress SDRF files while preserving duplicate columns."""
from __future__ import annotations

import csv
import json
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
INPUT_DIR = ROOT / "data" / "interim" / "G5_dataset_search"
OUTPUT = INPUT_DIR / "g5_arrayexpress_sample_metadata.tsv"
ACCESSIONS = ("E-MTAB-10985", "E-MTAB-13937")


def unique_headers(headers: list[str]) -> list[str]:
    seen: Counter[str] = Counter()
    output: list[str] = []
    for header in headers:
        seen[header] += 1
        output.append(header if seen[header] == 1 else f"{header}__{seen[header]}")
    return output


def main() -> None:
    records: list[dict[str, str]] = []
    summaries: dict[str, dict[str, object]] = {}
    for accession in ACCESSIONS:
        path = INPUT_DIR / f"{accession}.sdrf.txt"
        with path.open("r", encoding="utf-8-sig", newline="") as handle:
            reader = csv.reader(handle, delimiter="\t")
            headers = unique_headers(next(reader))
            accession_records = [dict(zip(headers, row, strict=False)) for row in reader]
        for record in accession_records:
            records.append({"accession": accession, **record})
        factor_fields = [field for field in headers if field.startswith("Factor Value[")]
        summaries[accession] = {
            "rows": len(accession_records),
            "factors": {
                field: dict(sorted(Counter(row.get(field, "") for row in accession_records).items()))
                for field in factor_fields
            },
            "individuals": dict(
                sorted(Counter(row.get("Characteristics[individual]", "") for row in accession_records).items())
            ),
        }
    fields = ["accession"] + sorted({key for record in records for key in record if key != "accession"})
    with OUTPUT.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(records)
    print(json.dumps({"output": str(OUTPUT), "summaries": summaries}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
