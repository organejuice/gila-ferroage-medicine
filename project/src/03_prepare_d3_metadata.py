#!/usr/bin/env python3
"""Parse GEO family SOFT metadata for the GSE160936 D3 sample manifest."""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import pandas as pd


def parse_soft(path: Path) -> list[dict[str, str]]:
    records: list[dict[str, str]] = []
    current: dict[str, str] | None = None
    characteristics: list[str] = []
    for raw in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if raw.startswith("^SAMPLE = "):
            if current is not None:
                current["characteristics"] = " | ".join(characteristics)
                records.append(current)
            current = {"geo_sample": raw.split(" = ", 1)[1].strip()}
            characteristics = []
            continue
        if current is None or not raw.startswith("!Sample_") or " = " not in raw:
            continue
        key, value = raw.split(" = ", 1)
        key = key.removeprefix("!Sample_")
        value = value.strip()
        if key == "characteristics_ch1":
            characteristics.append(value)
        elif key == "supplementary_file_1":
            current["supplementary_file"] = value
        elif key in {"title", "description", "source_name_ch1"}:
            current[key] = value
    if current is not None:
        current["characteristics"] = " | ".join(characteristics)
        records.append(current)
    return records


def characteristics_map(text: str) -> dict[str, str]:
    values: dict[str, str] = {}
    for item in text.split(" | "):
        if ":" not in item:
            continue
        key, value = item.split(":", 1)
        values[key.strip().lower()] = value.strip()
    return values


def parse_numeric(value: str | None) -> float | None:
    if value is None:
        return None
    match = re.search(r"[-+]?\d+(?:\.\d+)?", value)
    return float(match.group(0)) if match else None


def prepare(root: Path) -> pd.DataFrame:
    soft = root / "data" / "interim" / "GSE160936_family_soft.txt"
    raw_dir = root / "data" / "raw" / "D3_GSE160936"
    records = parse_soft(soft)
    rows = []
    for record in records:
        chars = characteristics_map(record.get("characteristics", ""))
        state_raw = chars.get("disease state", "")
        if state_raw.lower().startswith("non-disease"):
            disease = "Control"
        elif state_raw.lower().startswith("alzheimer") or state_raw.upper() == "AD":
            disease = "AD"
        else:
            disease = "Unresolved"
        region_raw = chars.get("brain region (somatosensory or entorhinal cortex)", "")
        region = region_raw.upper().strip()
        title = record.get("title", "")
        donor = re.sub(r"\s+(SSC|EC)$", "", title).strip()
        archive_name = Path(record.get("supplementary_file", "")).name
        archive = raw_dir / archive_name
        rows.append({
            "geo_sample": record["geo_sample"],
            "sample_title": title,
            "donor_id": donor,
            "disease": disease,
            "disease_state_raw": state_raw,
            "region": region,
            "age": parse_numeric(chars.get("age")),
            "sex": chars.get("sex", chars.get("Sex", "")),
            "rin": parse_numeric(chars.get("rin")),
            "braak_tangle_stage": parse_numeric(chars.get("braak tangle stage")),
            "ptau_percent_positive_cells": parse_numeric(chars.get("ptau immunostaining (%ptau positive cells)")),
            "amyloid_percent_area": parse_numeric(chars.get("amyloid immunostaining (%area stained)")),
            "archive": str(archive.relative_to(root)),
            "archive_exists": archive.exists(),
            "archive_size": archive.stat().st_size if archive.exists() else 0,
            "source": "GSE160936 family SOFT",
        })
    table = pd.DataFrame(rows).sort_values(["donor_id", "region", "geo_sample"]).reset_index(drop=True)
    if len(table) != 24:
        raise ValueError(f"Expected 24 D3 samples, found {len(table)}")
    if set(table["disease"]) != {"AD", "Control"}:
        raise ValueError(f"Unexpected D3 disease labels: {sorted(table['disease'].unique())}")
    if table["archive_exists"].sum() != 24:
        raise ValueError("Not all D3 archives are present")
    out = root / "data" / "interim" / "D3_sample_metadata.tsv"
    out.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(out, sep="\t", index=False)
    summary = {
        "sample_count": int(len(table)),
        "donor_count": int(table["donor_id"].nunique()),
        "samples_by_disease": table.groupby("disease").size().to_dict(),
        "samples_by_region": table.groupby("region").size().to_dict(),
        "donors_by_disease": table[["donor_id", "disease"]].drop_duplicates().groupby("disease").size().to_dict(),
        "metadata_table": str(out.relative_to(root)),
    }
    (root / "data" / "interim" / "D3_metadata_audit.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    return table


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    args = parser.parse_args()
    table = prepare(args.root.resolve())
    print(table[["geo_sample", "sample_title", "donor_id", "disease", "region", "age", "sex", "rin"]].to_string(index=False))


if __name__ == "__main__":
    main()