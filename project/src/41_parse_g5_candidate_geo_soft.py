#!/usr/bin/env python3
"""Parse downloaded GEO family SOFT metadata for G5 candidate screening."""
from __future__ import annotations

import csv
import gzip
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
IN_DIR = ROOT / "data" / "interim" / "G5_dataset_search" / "geo_soft"
OUT_DIR = ROOT / "data" / "interim" / "G5_dataset_search"


def add_value(record: dict[str, object], key: str, value: str) -> None:
    if key not in record:
        record[key] = value
    elif isinstance(record[key], list):
        record[key].append(value)
    else:
        record[key] = [record[key], value]


def joined(value: object) -> str:
    if isinstance(value, list):
        return " | ".join(str(item) for item in value)
    return "" if value is None else str(value)


def parse_soft(path: Path) -> tuple[dict[str, object], list[dict[str, str]]]:
    series: dict[str, object] = {}
    samples: list[dict[str, str]] = []
    current: dict[str, object] | None = None
    section = ""
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.rstrip("\r\n")
            if line.startswith("^SERIES = "):
                section = "series"
                series["accession"] = line.split("=", 1)[1].strip()
                continue
            if line.startswith("^SAMPLE = "):
                if current is not None:
                    samples.append({key: joined(value) for key, value in current.items()})
                section = "sample"
                current = {"sample_accession": line.split("=", 1)[1].strip()}
                continue
            if not line.startswith("!") or " = " not in line:
                continue
            key, value = line[1:].split(" = ", 1)
            if section == "series" and key.startswith("Series_"):
                add_value(series, key.removeprefix("Series_"), value)
            elif section == "sample" and current is not None and key.startswith("Sample_"):
                short_key = key.removeprefix("Sample_")
                if short_key == "characteristics_ch1" and ":" in value:
                    char_key, char_value = value.split(":", 1)
                    current[f"characteristic_{char_key.strip().lower().replace(' ', '_')}"] = char_value.strip()
                else:
                    add_value(current, short_key, value)
    if current is not None:
        samples.append({key: joined(value) for key, value in current.items()})
    return series, samples


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    series_rows: list[dict[str, str]] = []
    sample_rows: list[dict[str, str]] = []
    for path in sorted(IN_DIR.glob("GSE*_family.soft.gz")):
        series, samples = parse_soft(path)
        accession = str(series.get("accession", path.name.split("_", 1)[0]))
        series_rows.append({
            "accession": accession,
            "title": joined(series.get("title")),
            "summary": joined(series.get("summary")),
            "overall_design": joined(series.get("overall_design")),
            "supplementary_files": joined(series.get("supplementary_file")),
            "citation": joined(series.get("citation")),
            "sample_count": str(len(samples)),
            "soft_file": path.relative_to(ROOT).as_posix(),
        })
        for sample in samples:
            sample_rows.append({"series_accession": accession, **sample})

    series_fields = ["accession", "title", "summary", "overall_design", "supplementary_files", "citation", "sample_count", "soft_file"]
    with (OUT_DIR / "g5_candidate_series_summary.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=series_fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(series_rows)
    sample_fields = sorted({key for row in sample_rows for key in row}, key=lambda key: (key != "series_accession", key != "sample_accession", key))
    with (OUT_DIR / "g5_candidate_sample_metadata.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=sample_fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(sample_rows)
    print(json.dumps({"series": len(series_rows), "samples": len(sample_rows), "output": str(OUT_DIR)}, ensure_ascii=False))


if __name__ == "__main__":
    main()
