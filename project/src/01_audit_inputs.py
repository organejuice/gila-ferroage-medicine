#!/usr/bin/env python3
"""Audit the locally staged inputs for the AD disease-side v2.0 plan.

This stage is deliberately conservative. It validates file structure, sample
labels, donor/cell metadata, and the frozen D5 reference gene set. It does not
perform peak-to-gene annotation or infer a disease direction from a single
cell-level observation.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import re
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable

import pandas as pd


EXPECTED_DATASETS = {
    "D1": "D1_GSE188765",
    "D2": "D2_GSE227221",
    "D5": "D5_CellMetab2026",
    "D6": "D6_GSE214979",
}


def sha256_file(path: Path, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def parse_soft_series(path: Path) -> dict[str, list[str]]:
    """Read the series-level portion of a GEO SOFT text export."""

    values: dict[str, list[str]] = defaultdict(list)
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.rstrip("\r\n")
            if not line.startswith("!Series_") or " = " not in line:
                continue
            key, value = line.split(" = ", 1)
            values[key.removeprefix("!Series_")].append(value)
    return dict(values)


def parse_d1_name(path: Path) -> dict[str, str]:
    match = re.match(r"(?P<gsm>GSM\d+)_(?P<context>[A-Za-z]+)_(?P<condition>WT|AD)-(?P<replicate>\d+)_peaks\.txt\.gz$", path.name)
    if not match:
        raise ValueError(f"Unexpected D1 peak filename: {path.name}")
    return match.groupdict()


def read_d1_peaks(path: Path) -> tuple[int, set[tuple[str, int, int]]]:
    total_rows = 0
    peaks: set[tuple[str, int, int]] = set()
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) < 3:
                raise ValueError(f"D1 row has fewer than 3 columns: {path.name}:{line_number}")
            chrom, start, end = fields[:3]
            try:
                start_i, end_i = int(start), int(end)
            except ValueError as exc:
                raise ValueError(f"D1 coordinates are not integers: {path.name}:{line_number}") from exc
            if start_i < 0 or end_i <= start_i:
                raise ValueError(f"D1 interval is invalid: {path.name}:{line_number}")
            total_rows += 1
            peaks.add((chrom, start_i, end_i))
    return total_rows, peaks


def audit_d1(root: Path, tables_dir: Path) -> dict:
    directory = root / "data" / "raw" / EXPECTED_DATASETS["D1"]
    rows: list[dict] = []
    peak_sets: dict[str, set[tuple[str, int, int]]] = {}
    for path in sorted(directory.glob("*.txt.gz")):
        labels = parse_d1_name(path)
        total_rows, peaks = read_d1_peaks(path)
        peak_sets[labels["gsm"]] = peaks
        rows.append(
            {
                **labels,
                "file": str(path.relative_to(root)),
                "total_rows": total_rows,
                "unique_peaks": len(peaks),
                "duplicate_rows": total_rows - len(peaks),
                "duplicate_fraction": (total_rows - len(peaks)) / total_rows if total_rows else None,
                "sha256": sha256_file(path),
            }
        )
    if len(rows) != 8:
        raise ValueError(f"Expected 8 D1 sample peak files, found {len(rows)}")

    pd.DataFrame(rows).to_csv(tables_dir / "D1_peak_summary.tsv", sep="\t", index=False)

    overlap_rows: list[dict] = []
    for context in sorted({row["context"] for row in rows}):
        by_condition = {
            condition: set().union(*(peak_sets[row["gsm"]] for row in rows if row["context"] == context and row["condition"] == condition))
            for condition in ("WT", "AD")
        }
        intersection = by_condition["WT"] & by_condition["AD"]
        union = by_condition["WT"] | by_condition["AD"]
        overlap_rows.append(
            {
                "context": context,
                "wt_union_peaks": len(by_condition["WT"]),
                "ad_union_peaks": len(by_condition["AD"]),
                "shared_peaks": len(intersection),
                "union_peaks": len(union),
                "jaccard": len(intersection) / len(union) if union else None,
                "note": "descriptive interval overlap only; no peak-to-gene annotation",
            }
        )
    pd.DataFrame(overlap_rows).to_csv(tables_dir / "D1_peak_overlap.tsv", sep="\t", index=False)
    return {"sample_count": len(rows), "summary": rows, "overlap": overlap_rows}


def audit_d2(root: Path) -> dict:
    paths = sorted((root / "data" / "raw" / EXPECTED_DATASETS["D2"]).glob("*.gz"))
    if len(paths) != 1:
        raise ValueError(f"Expected one D2 TPM file, found {len(paths)}")
    path = paths[0]
    with gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="") as handle:
        header = handle.readline().rstrip("\r\n").split("\t")
        if len(header) < 2:
            raise ValueError("D2 TPM header has fewer than two columns")
        sample_ids = header
        rows = 0
        bad_rows = 0
        transcript_ids: list[str] = []
        for line in handle:
            fields = line.rstrip("\r\n").split("\t")
            if not fields or not fields[0]:
                continue
            rows += 1
            transcript_ids.append(fields[0])
            # The first data column is the transcript ID; the header contains only sample IDs.
            if len(fields) != len(sample_ids) + 1:
                bad_rows += 1
    versioned = sum(1 for item in transcript_ids[:1000] if re.match(r"^ENS[A-Z]*\d+\.\d+$", item))
    return {
        "file": str(path.relative_to(root)),
        "sha256": sha256_file(path),
        "sample_count": len(sample_ids),
        "sample_ids": sample_ids,
        "feature_rows": rows,
        "rows_with_wrong_column_count": bad_rows,
        "identifier_type": "versioned Ensembl transcript IDs" if versioned else "unresolved",
        "first_transcript_ids": transcript_ids[:5],
    }


def write_d2_sample_metadata(root: Path, tables_dir: Path) -> dict:
    """Freeze the 12-sample iMGL time-course labels from the GEO series record."""
    d2_path = next((root / "data" / "raw" / EXPECTED_DATASETS["D2"]).glob("*.gz"))
    with gzip.open(d2_path, "rt", encoding="utf-8", errors="replace") as handle:
        sample_ids = handle.readline().rstrip("\r\n").split("\t")
    titles = ["BaselineA", "BaselineB", "BaselineC", "Treat4HrsA", "Treat4HrsB", "Treat4HrsC", "Treat24HrsA", "Treat24HrsB", "Treat24HrsC", "Treat96HrsA", "Treat96HrsB", "Treat96HrsC"]
    times = [0, 0, 0, 4, 4, 4, 24, 24, 24, 96, 96, 96]
    if len(sample_ids) != len(titles):
        raise ValueError("D2 sample count does not match the GEO time-course manifest")
    rows = []
    for index, (sample_id, title, time_hours) in enumerate(zip(sample_ids, titles, times), start=1):
        rows.append({"sample_id_in_expression_table": sample_id, "geo_sample": f"GSM{7093346 + index}", "sample_title": title, "treatment": "Control" if time_hours == 0 else "Amyloid_beta_1_42_fibrils", "time_hours": time_hours, "replicate": title[-1], "source": "GEO GSE227223 sample labels; processed matrix from GSE227221"})
    path = tables_dir / "D2_sample_metadata.tsv"
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)
    return {"sample_count": len(rows), "metadata_table": str(path.relative_to(root)), "time_hours": sorted(set(times))}


def audit_d3_downloads(root: Path, tables_dir: Path) -> dict:
    """Track sample-level D3 processed matrices without accepting partial archives."""
    filelist = root / "data" / "interim" / "GSE160936_filelist.txt"
    raw_dir = root / "data" / "raw" / "D3_GSE160936"
    rows = []
    for line in filelist.read_text(encoding="utf-8", errors="replace").splitlines():
        parts = line.split("\t")
        if len(parts) < 4 or parts[0] != "File":
            continue
        name, expected_size = parts[1], int(parts[3])
        path = raw_dir / name
        actual_size = path.stat().st_size if path.exists() else 0
        status = "complete" if actual_size == expected_size else ("partial" if actual_size else "missing")
        rows.append({"dataset": "D3_GSE160936", "file": str(path.relative_to(root)), "expected_size": expected_size, "actual_size": actual_size, "status": status})
    if len(rows) != 24:
        raise ValueError(f"Expected 24 D3 sample files in the GEO file list, found {len(rows)}")
    pd.DataFrame(rows).to_csv(tables_dir / "D3_download_manifest.tsv", sep="\t", index=False)
    counts = Counter(row["status"] for row in rows)
    return {"expected_sample_files": len(rows), "complete_sample_files": counts.get("complete", 0), "partial_sample_files": counts.get("partial", 0), "missing_sample_files": counts.get("missing", 0), "manifest": str((tables_dir / "D3_download_manifest.tsv").relative_to(root))}

def audit_d5(root: Path, tables_dir: Path) -> dict:
    paths = sorted((root / "data" / "raw" / EXPECTED_DATASETS["D5"]).glob("Table_S2*.xlsx"))
    if len(paths) != 1:
        raise ValueError(f"Expected one D5 Table S2 workbook, found {len(paths)}")
    path = paths[0]
    raw = pd.read_excel(path, sheet_name="Ferroaging gene set", header=None)
    if raw.shape[1] < 2:
        raise ValueError("D5 ferro-aging sheet does not contain gene symbols")
    symbols = (
        raw.iloc[1:, 1]
        .dropna()
        .astype(str)
        .str.strip()
        .loc[lambda series: series.str.upper().ne("SYMBOL")]
        .loc[lambda series: series.ne("")]
        .drop_duplicates()
        .tolist()
    )
    if len(symbols) < 50:
        raise ValueError(f"D5 ferro-aging reference unexpectedly small: {len(symbols)} genes")
    d5_path = tables_dir / "D5_ferroaging_reference.tsv"
    pd.DataFrame(
        {
            "module": "M7",
            "gene_symbol": symbols,
            "source_dataset": "D5_CellMetab2026",
            "reference_type": "author_ferroaging_gene_set",
        }
    ).to_csv(d5_path, sep="\t", index=False)
    return {
        "file": str(path.relative_to(root)),
        "sha256": sha256_file(path),
        "sheet": "Ferroaging gene set",
        "gene_count": len(symbols),
        "first_symbols": symbols[:10],
        "reference_table": str(d5_path.relative_to(root)),
    }


def audit_d6(root: Path, tables_dir: Path) -> dict:
    paths = sorted((root / "data" / "raw" / EXPECTED_DATASETS["D6"]).glob("*.csv.gz"))
    if len(paths) != 1:
        raise ValueError(f"Expected one D6 metadata file, found {len(paths)}")
    path = paths[0]
    metadata = pd.read_csv(path, low_memory=False)
    required = {"id", "Status", "predicted.id", "Age", "Sex", "PMI..h.", "RIN", "structure"}
    missing = sorted(required.difference(metadata.columns))
    if missing:
        raise ValueError(f"D6 metadata missing required columns: {missing}")
    donor = metadata[["id", "Status", "Age", "Sex", "PMI..h.", "RIN"]].drop_duplicates()
    donor_counts = donor.groupby("Status", dropna=False).size().to_dict()
    cell_counts = metadata.groupby(["Status", "predicted.id"], dropna=False).size().reset_index(name="cell_count")
    cell_counts.to_csv(tables_dir / "D6_celltype_status_counts.tsv", sep="\t", index=False)
    donor.to_csv(tables_dir / "D6_donor_metadata.tsv", sep="\t", index=False)
    return {
        "file": str(path.relative_to(root)),
        "sha256": sha256_file(path),
        "cell_count": int(len(metadata)),
        "donor_count": int(len(donor)),
        "donors_by_status": {str(key): int(value) for key, value in donor_counts.items()},
        "cell_types": sorted(metadata["predicted.id"].dropna().astype(str).unique()),
        "statuses": sorted(metadata["Status"].dropna().astype(str).unique()),
        "missing_required_values": {column: int(metadata[column].isna().sum()) for column in sorted(required)},
        "metadata_columns": list(metadata.columns),
    }


def audit_geo_series(root: Path, tables_dir: Path) -> list[dict]:
    rows: list[dict] = []
    for dataset_id in ("GSE160936", "GSE214979", "GSE227221"):
        path = root / "data" / "interim" / f"{dataset_id}_geo_full.txt"
        series = parse_soft_series(path)
        rows.append(
            {
                "dataset": dataset_id,
                "title": series.get("title", [""])[0],
                "status": series.get("status", [""])[0],
                "series_sample_count": len(series.get("sample_id", [])),
                "platform": series.get("platform_id", [""])[0],
                "organism": series.get("platform_organism", [""])[0],
                "supplementary_files": " | ".join(series.get("supplementary_file", [])),
                "source_file": str(path.relative_to(root)),
            }
        )
    pd.DataFrame(rows).to_csv(tables_dir / "GEO_series_audit.tsv", sep="\t", index=False)
    return rows


def write_signature_registry(root: Path, tables_dir: Path, m7_gene_count: int) -> None:
    """Write the auditable module registry without inventing unfinished gene sets."""
    rows = [
        {"module": "M1", "status": "planned", "feature_count": "pending", "source": "HALLMARK_GLYCOLYSIS plus Reactome support", "note": "Freeze after human gene-symbol harmonization and sensitivity review."},
        {"module": "M2", "status": "pending_D1_annotation", "feature_count": "pending", "source": "D1_GSE188765 H4K12la peaks", "note": "Requires pinned mm10 peak-to-gene annotation; do not call this a gene set yet."},
        {"module": "M3", "status": "planned", "feature_count": "pending", "source": "SenMayo, CellAge/SASP, cell-cycle arrest references", "note": "Use donor-level score and preserve component scores."},
        {"module": "M4", "status": "planned", "feature_count": "pending", "source": "C1 state association, then LIANA/CellChat/NicheNet", "note": "Communication analysis is conditional on reproducible state association."},
        {"module": "M5", "status": "planned", "feature_count": "pending", "source": "Iron handling subcomponents", "note": "Retain subcomponents if a composite is unstable."},
        {"module": "M6", "status": "planned", "feature_count": "pending", "source": "ACSL4-centered lipid peroxidation program", "note": "Freeze after source and species harmonization."},
        {"module": "M7", "status": "frozen_reference", "feature_count": str(m7_gene_count), "source": "D5_CellMetab2026 Table S2 Ferroaging gene set", "note": "Author reference signature; interpret as ferro-aging-related if component concordance fails."},
    ]
    pd.DataFrame(rows).to_csv(tables_dir / "SignatureRegistry.tsv", sep="\t", index=False)

def audit_all(root: Path) -> dict:
    tables_dir = root / "results" / "tables"
    interim_dir = root / "data" / "interim"
    tables_dir.mkdir(parents=True, exist_ok=True)
    interim_dir.mkdir(parents=True, exist_ok=True)
    result = {
        "project_root": str(root),
        "datasets": {
            "D1": audit_d1(root, tables_dir),
            "D2": audit_d2(root),
            "D2_sample_metadata": write_d2_sample_metadata(root, tables_dir),
            "D3": audit_d3_downloads(root, tables_dir),
            "D5": audit_d5(root, tables_dir),
            "D6": audit_d6(root, tables_dir),
        },
        "geo_series": audit_geo_series(root, tables_dir),
        "limitations": [
            "D1 is audited at the interval level; peak-to-gene annotation is intentionally not performed in this stage.",
            "The local GEO text exports contain series-level records only; sample-level metadata must be retrieved separately before D3 donor pseudobulk modeling.",
            "D2 is transcript-level TPM and needs an explicit transcript-to-gene aggregation policy before differential analysis.",
        ],
    }
    write_signature_registry(root, tables_dir, result["datasets"]["D5"]["gene_count"])
    out = interim_dir / "input_audit.json"
    out.write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
    return result


def run_checks(root: Path) -> None:
    audit_path = root / "data" / "interim" / "input_audit.json"
    registry_path = root / "results" / "tables" / "D5_ferroaging_reference.tsv"
    if not audit_path.exists() or not registry_path.exists():
        raise AssertionError("Input audit outputs are missing; run the audit first")
    audit = json.loads(audit_path.read_text(encoding="utf-8"))
    assert audit["datasets"]["D1"]["sample_count"] == 8
    assert audit["datasets"]["D2"]["sample_count"] == 12
    assert audit["datasets"]["D2"]["rows_with_wrong_column_count"] == 0
    assert audit["datasets"]["D5"]["gene_count"] >= 90
    assert audit["datasets"]["D6"]["cell_count"] == 105332
    assert set(audit["datasets"]["D6"]["statuses"]) == {"AD", "Ctrl"}
    assert audit["datasets"]["D6"]["donor_count"] == 15
    print("input audit checks: PASS")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--check", action="store_true", help="check previously generated outputs")
    args = parser.parse_args()
    root = args.root.resolve()
    if args.check:
        run_checks(root)
        return
    result = audit_all(root)
    print(json.dumps({
        "D1_samples": result["datasets"]["D1"]["sample_count"],
        "D2_samples": result["datasets"]["D2"]["sample_count"],
        "D2_features": result["datasets"]["D2"]["feature_rows"],
        "D5_genes": result["datasets"]["D5"]["gene_count"],
        "D6_cells": result["datasets"]["D6"]["cell_count"],
        "D6_donors": result["datasets"]["D6"]["donor_count"],
    }, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
