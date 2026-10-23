#!/usr/bin/env python3
"""Create animal-level astrocyte pseudobulk counts from GSE143758.

The GEO supplementary matrices are already restricted to astrocyte nuclei.
This script applies a transparent cell-level QC rule, aggregates retained
cells by the source sample/animal, and preserves the source-defined region,
age, and 5xFAD/WT labels. Raw files are never modified.
"""
from __future__ import annotations

import csv
import gzip
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw" / "GSE143758"
OUT = ROOT / "data" / "interim" / "GSE143758"

FILES = [
    RAW / "GSE143758_Admouse_Hippocampus_3TimeCourses_Astrocytes_WT_AD_Aging.csv.gz",
    RAW / "GSE143758_Admouse_Crtx_7-10m_Astrocytes_UMIcounts.csv.gz",
]


def parse_gsm_metadata(path: Path) -> dict[str, dict[str, str]]:
    records: dict[str, dict[str, str]] = {}
    current: dict[str, str] | None = None
    current_id: str | None = None
    for raw in path.read_text(encoding="utf-8", errors="replace").splitlines():
        line = raw.strip("\r")
        if line.startswith("^SAMPLE = "):
            if current_id and current is not None:
                records[current_id] = current
            current_id = line.split("=", 1)[1].strip()
            current = {"geo_sample": current_id}
        elif current is not None and line.startswith("!Sample_title = "):
            current["title"] = line.split("=", 1)[1].strip()
        elif current is not None and line.startswith("!Sample_characteristics_ch1 = "):
            value = line.split("=", 1)[1].strip()
            if ":" in value:
                key, val = value.split(":", 1)
                current[key.strip().lower()] = val.strip()
    if current_id and current is not None:
        records[current_id] = current
    return records


def sample_prefix(cell: str) -> str:
    match = re.match(r"^(.*)_([ACGTN]{16})(?:-1)?$", str(cell))
    if not match:
        raise ValueError(f"Cannot parse cell barcode: {cell}")
    return match.group(1)


def numeric_age(value: str | None, title: str) -> float:
    if value:
        match = re.search(r"\d+(?:\.\d+)?", value)
        if match:
            return float(match.group(0))
    match = re.search(r"(\d+(?:\.\d+)?)m", title, flags=re.IGNORECASE)
    if match:
        return float(match.group(1))
    return float("nan")


def animal_id(prefix: str) -> str:
    match = re.search(r"(?:WT|AD)(\d+)$", prefix, flags=re.IGNORECASE)
    if match:
        return match.group(1)
    return prefix


def source_meta(prefix: str, gsm: dict[str, dict[str, str]]) -> dict[str, str | float]:
    aliases = {
        "HipR-AD-G1-8w": "HipR-AD-G1-2m",
        "HipR-AD-G2-8w": "HipR-AD-G2-2m",
        "HipR-WT-G1-8w": "HipR-WT-G1-2m",
        "HipR-WT-G2-8w": "HipR-WT-G2-2m",
        "Untreated-Hip-10m-S135": "HIP_10mAD135",
        "Untreated-Hip-4m-S186": "HIP_4mAD186",
        "Untreated-Hip-6w-S331": "HIP_6wAD331",
        "Untreated-Hip-7m-S343": "HIP_7mAD343",
        "WT-Hip-4m-S185": "HIP_4mWT185",
        "WT-Hip-6w-S331": "HIP_6wWT330",
        "Wt-Hip-10m-S124": "HIP_10mWT124",
        "Wt-Hip-7m-S341": "HIP_7mWT341",
        "M15": "AD_M15_4.7m",
        "M206": "WT_M206_19.6m",
        "M208": "WT_M208_19.6m",
        "M268": "WT_M768_13.3m",
        "M680": "WT_M680_8.1m",
        "M704": "AD_M704_7.3m",
        "M721": "WT_M721_4.7m",
        "M766": "WT_M766_13.3m",
    }
    lookup = prefix.removeprefix("EZ.Batch1_")
    lookup = aliases.get(lookup, lookup)
    matches = [row for row in gsm.values() if row.get("title") == lookup]
    row = matches[0] if matches else {}
    strain = row.get("strain", "")
    if strain.lower() == "5xfad":
        disease = "AD_model"
    elif strain:
        disease = "WT_control"
    else:
        disease = "unknown"
    region = row.get("region", "")
    if not region:
        region = "Cortex" if prefix.upper().startswith("CRTX") else "Hippocampus"
    return {
        "sample_prefix": prefix,
        "geo_sample": row.get("geo_sample", ""),
        "title": row.get("title", prefix),
        "strain": strain,
        "disease": disease,
        "region": region,
        "age_months": numeric_age(row.get("age"), prefix),
        "animal_id": animal_id(prefix),
    }


def get_columns(path: Path) -> list[str]:
    header = pd.read_csv(path, sep=";", compression="gzip", nrows=0)
    return [str(x) for x in header.columns[1:]]


def read_qc(path: Path, columns: list[str], groups: list[str]) -> tuple[pd.DataFrame, dict[str, np.ndarray]]:
    n = len(columns)
    n_counts = np.zeros(n, dtype=np.float64)
    n_features = np.zeros(n, dtype=np.int64)
    mt_counts = np.zeros(n, dtype=np.float64)
    group_index = {group: np.flatnonzero(np.asarray(groups) == group) for group in sorted(set(groups))}
    first = True
    for chunk in pd.read_csv(path, sep=";", compression="gzip", chunksize=500, low_memory=False):
        genes = chunk.iloc[:, 0].astype(str).str.upper().to_numpy()
        values = chunk.iloc[:, 1:].to_numpy(dtype=np.float64, copy=False)
        if values.shape[1] != n:
            raise ValueError(f"Column count changed while reading {path}")
        n_counts += values.sum(axis=0)
        n_features += (values > 0).sum(axis=0)
        mt_counts += values[np.char.startswith(genes.astype(str), "MT-")].sum(axis=0)
        first = False
    if first:
        raise ValueError(f"Empty matrix: {path}")
    mt_pct = 100.0 * mt_counts / np.maximum(n_counts, 1.0)
    keep = (n_counts >= 500) & (n_features >= 200) & (mt_pct <= 20)
    cell_qc = pd.DataFrame(
        {
            "cell_id": columns,
            "sample_prefix": groups,
            "n_count": n_counts,
            "n_feature": n_features,
            "mt_pct": mt_pct,
            "keep": keep,
        }
    )
    return cell_qc, group_index


def aggregate(path: Path, columns: list[str], groups: list[str], keep: np.ndarray) -> tuple[pd.DataFrame, list[str]]:
    unique_groups = sorted(set(groups))
    group_indices = {group: np.flatnonzero((np.asarray(groups) == group) & keep) for group in unique_groups}
    chunks: list[pd.DataFrame] = []
    for chunk in pd.read_csv(path, sep=";", compression="gzip", chunksize=500, low_memory=False):
        genes = chunk.iloc[:, 0].astype(str).to_numpy()
        values = chunk.iloc[:, 1:].to_numpy(dtype=np.float64, copy=False)
        out = {"gene": genes}
        for group in unique_groups:
            idx = group_indices[group]
            out[group] = values[:, idx].sum(axis=1) if idx.size else np.zeros(len(genes))
        chunks.append(pd.DataFrame(out))
    result = pd.concat(chunks, ignore_index=True)
    result = result.groupby("gene", sort=False, as_index=False)[unique_groups].sum()
    return result, unique_groups


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    gsm_path = RAW / "GSE143758_gsm_metadata.txt"
    gsm = parse_gsm_metadata(gsm_path)
    all_qc: list[pd.DataFrame] = []
    all_meta: list[dict] = []
    matrices: list[pd.DataFrame] = []
    matrix_groups: list[str] = []
    for path in FILES:
        if not path.exists() or path.stat().st_size == 0:
            raise FileNotFoundError(f"Missing or incomplete raw matrix: {path}")
        columns = get_columns(path)
        groups = [sample_prefix(cell) for cell in columns]
        cell_qc, _ = read_qc(path, columns, groups)
        cell_qc["source_file"] = path.name
        all_qc.append(cell_qc)
        keep = cell_qc["keep"].to_numpy(bool)
        matrix, unique_groups = aggregate(path, columns, groups, keep)
        region = "Cortex" if "Crtx" in path.name else "Hippocampus"
        renamed = []
        for group in unique_groups:
            meta = source_meta(group, gsm)
            meta.update(
                {
                    "group_id": f"GSE143758__{region}__{group}__Astrocyte",
                    "dataset_id": "GSE143758",
                    "cell_type": "Astrocyte",
                    "source_file": path.name,
                    "n_cells_total": int((np.asarray(groups) == group).sum()),
                    "n_cells_used": int(((np.asarray(groups) == group) & keep).sum()),
                }
            )
            if meta["disease"] == "unknown":
                raise ValueError(f"Could not resolve strain for sample {group}")
            all_meta.append(meta)
            renamed.append(meta["group_id"])
        matrix = matrix.rename(columns=dict(zip(unique_groups, renamed)))
        matrices.append(matrix)
        matrix_groups.extend(renamed)
    combined = matrices[0]
    for matrix in matrices[1:]:
        if matrix["gene"].duplicated().any() or combined["gene"].duplicated().any():
            raise ValueError("GSE143758 matrix contains duplicated gene identifiers")
        combined = combined.merge(matrix, on="gene", how="outer", validate="one_to_one")
    combined = combined[["gene"] + matrix_groups].fillna(0)
    combined.to_csv(OUT / "GSE143758_astrocyte_pseudobulk_counts.tsv.gz", sep="\t", index=False, compression="gzip")
    metadata = pd.DataFrame(all_meta).sort_values(["region", "age_months", "disease", "sample_prefix"])
    metadata.to_csv(OUT / "GSE143758_astrocyte_pseudobulk_metadata.tsv", sep="\t", index=False)
    qc = pd.concat(all_qc, ignore_index=True)
    qc.to_csv(OUT / "GSE143758_astrocyte_cell_qc.tsv.gz", sep="\t", index=False, compression="gzip")
    summary = {
        "dataset": "GSE143758",
        "input_files": [p.name for p in FILES],
        "sample_count": int(len(metadata)),
        "animal_count": int(metadata["animal_id"].nunique()),
        "regions": sorted(metadata["region"].unique().tolist()),
        "cells_total": int(len(qc)),
        "cells_used": int(qc["keep"].sum()),
        "qc_rule": "n_count>=500; n_feature>=200; mt_pct<=20",
        "raw_data_modified": False,
    }
    (OUT / "GSE143758_preparation_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
