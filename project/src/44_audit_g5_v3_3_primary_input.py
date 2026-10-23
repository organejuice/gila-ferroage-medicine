#!/usr/bin/env python3
"""Audit frozen G5 v3.3 primary inputs without altering GEO source files."""
from __future__ import annotations

import gzip
import hashlib
import json
import re
import tarfile
from pathlib import Path

import numpy as np
import pandas as pd
import yaml


ROOT = Path(__file__).resolve().parents[1]
CFG = yaml.safe_load((ROOT / "configs/g5_v3.3_execution.yaml").read_text(encoding="utf-8"))
RUN = ROOT / CFG["outputs"]["result_dir"]
PRIMARY = CFG["inputs"]["primary"]
TOLERANCE = float(CFG["input_gate"]["fractional_tolerance"])
TIME_TO_HOUR = {
    "12am": 0, "2am": 2, "4am": 4, "6am": 6, "8am": 8, "10am": 10,
    "12pm": 12, "2pm": 14, "4pm": 16, "6pm": 18, "8pm": 20, "10pm": 22,
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def source_path(key: str) -> Path:
    return ROOT / PRIMARY[key]


def load_counts(path: Path) -> tuple[pd.DataFrame, list[str], pd.Series]:
    frame = pd.read_csv(path, compression="gzip", low_memory=False)
    if "gene_id" not in frame.columns:
        raise ValueError(f"gene_id column missing from {path}")
    samples = [column for column in frame.columns if column != "gene_id"]
    values = frame.loc[:, samples].apply(pd.to_numeric, errors="coerce")
    return values, samples, frame["gene_id"].astype(str)


def audit_matrix(label: str, values: pd.DataFrame, samples: list[str], gene_ids: pd.Series) -> dict:
    matrix = values.to_numpy(dtype=float, copy=False)
    finite = np.isfinite(matrix)
    nonnegative = bool(finite.all() and (matrix >= 0).all())
    fractional = finite & (np.abs(matrix - np.rint(matrix)) > TOLERANCE)
    return {
        "source_group": label,
        "n_genes": int(values.shape[0]),
        "n_samples": int(values.shape[1]),
        "gene_ids_unique": bool(gene_ids.is_unique),
        "sample_ids_unique": bool(len(samples) == len(set(samples))),
        "finite_values": bool(finite.all()),
        "nonnegative_values": nonnegative,
        "fractional_value_count": int(fractional.sum()),
        "total_value_count": int(matrix.size),
        "fractional_value_fraction": float(fractional.mean()),
        "minimum_fractional_value": float(matrix[fractional].min()) if fractional.any() else None,
        "maximum_fractional_value": float(matrix[fractional].max()) if fractional.any() else None,
        "minimum_library_size": float(matrix.sum(axis=0).min()),
        "maximum_library_size": float(matrix.sum(axis=0).max()),
    }


def parse_metadata(samples: list[str]) -> pd.DataFrame:
    pattern = re.compile(PRIMARY["sample_id_pattern"])
    rows: list[dict] = []
    for sample_id in samples:
        matched = pattern.fullmatch(sample_id)
        if not matched:
            rows.append({"sample_id": sample_id, "metadata_parse_pass": False})
            continue
        source_group, numeric_id, time_label, sex, batch = matched.groups()
        hour = TIME_TO_HOUR[time_label]
        rows.append({
            "sample_id": sample_id,
            "metadata_parse_pass": True,
            "source_group": source_group,
            "disease": "APP_PS1" if source_group == "APP" else "WT",
            "animal_id": f"{source_group}{numeric_id}",
            "clock_label": time_label,
            "clock_hour": hour,
            "sin_time": float(np.sin(2 * np.pi * hour / 24)),
            "cos_time": float(np.cos(2 * np.pi * hour / 24)),
            "sex": sex,
            "batch": f"batch{batch}",
        })
    return pd.DataFrame(rows)


def main() -> None:
    if RUN.exists():
        raise FileExistsError(f"Refusing to overwrite existing run directory: {RUN}")
    for directory in (RUN / "tables", RUN / "gates", RUN / "registries", RUN / "logs"):
        directory.mkdir(parents=True, exist_ok=True)

    app_path, wt_path = source_path("app_counts"), source_path("wt_counts")
    app_values, app_samples, app_genes = load_counts(app_path)
    wt_values, wt_samples, wt_genes = load_counts(wt_path)
    app_audit = audit_matrix("APP", app_values, app_samples, app_genes)
    wt_audit = audit_matrix("WT", wt_values, wt_samples, wt_genes)
    metadata = parse_metadata(app_samples + wt_samples)

    raw_tar_path = source_path("geo_raw_tar")
    with tarfile.open(raw_tar_path, "r") as archive:
        raw_tar_members = archive.getnames()
    raw_tar_inventory = pd.DataFrame({"member": raw_tar_members})
    raw_tar_inventory.to_csv(RUN / "registries/GSE261698_RAW_tar_inventory.tsv", sep="\t", index=False)

    expected = PRIMARY
    total_fractional = app_audit["fractional_value_count"] + wt_audit["fractional_value_count"]
    total_values = app_audit["total_value_count"] + wt_audit["total_value_count"]
    checks = {
        "input_files_exist": all(source_path(key).is_file() for key in ("app_counts", "wt_counts", "family_soft", "geo_raw_tar")),
        "source_hashes_match": all(
            sha256(source_path(key)) == expected["expected_sha256"][key]
            for key in ("app_counts", "wt_counts", "family_soft", "geo_raw_tar")
        ),
        "expected_gene_rows": app_audit["n_genes"] == expected["expected_n_genes"] == wt_audit["n_genes"],
        "identical_gene_order": bool(app_genes.equals(wt_genes)),
        "unique_gene_ids": bool(app_audit["gene_ids_unique"] and wt_audit["gene_ids_unique"]),
        "expected_samples": len(app_samples) == expected["expected_n_model"] and len(wt_samples) == expected["expected_n_control"],
        "unique_sample_ids": bool(app_audit["sample_ids_unique"] and wt_audit["sample_ids_unique"] and len(set(app_samples + wt_samples)) == len(app_samples + wt_samples)),
        "nonnegative_finite_counts": bool(app_audit["finite_values"] and wt_audit["finite_values"] and app_audit["nonnegative_values"] and wt_audit["nonnegative_values"]),
        "strict_integer_counts": total_fractional == 0,
        "complete_sample_metadata": bool(metadata["metadata_parse_pass"].all()),
        "expected_disease_counts": bool((metadata.get("disease", pd.Series(dtype=str)) == "APP_PS1").sum() == expected["expected_n_model"] and (metadata.get("disease", pd.Series(dtype=str)) == "WT").sum() == expected["expected_n_control"]),
        "complete_24_hour_coverage": set(metadata.get("clock_hour", pd.Series(dtype=float)).dropna().tolist()) == set(TIME_TO_HOUR.values()),
        "complete_raw_matrix_alternative": len(raw_tar_members) >= 48 * 3,
    }
    core_checks = [
        "input_files_exist", "source_hashes_match", "expected_gene_rows", "identical_gene_order",
        "unique_gene_ids", "expected_samples", "unique_sample_ids", "nonnegative_finite_counts",
        "strict_integer_counts", "complete_sample_metadata", "expected_disease_counts", "complete_24_hour_coverage",
    ]
    primary_input_pass = all(checks[key] for key in core_checks)

    metadata.to_csv(RUN / "tables/GSE261698_sample_metadata_audit.tsv", sep="\t", index=False)
    integrity_rows = []
    for key, value in checks.items():
        integrity_rows.append({"criterion": key, "observed": bool(value), "required_for_primary_input": key in core_checks})
    integrity_rows.extend([
        {"criterion": "total_fractional_values", "observed": total_fractional, "required_for_primary_input": True},
        {"criterion": "total_values", "observed": total_values, "required_for_primary_input": True},
        {"criterion": "fractional_value_fraction", "observed": total_fractional / total_values, "required_for_primary_input": True},
    ])
    integrity = pd.DataFrame(integrity_rows)
    integrity.to_csv(RUN / "tables/GSE261698_input_integrity.tsv", sep="\t", index=False)

    audit = {
        "run_id": CFG["run_id"],
        "dataset": PRIMARY["accession"],
        "raw_data_modified": False,
        "matrix_audits": {"APP": app_audit, "WT": wt_audit},
        "source_hashes": {key: sha256(source_path(key)) for key in ("app_counts", "wt_counts", "family_soft", "geo_raw_tar")},
        "metadata_summary": {
            "n_samples": int(len(metadata)),
            "n_APP_PS1": int((metadata.get("disease", pd.Series(dtype=str)) == "APP_PS1").sum()),
            "n_WT": int((metadata.get("disease", pd.Series(dtype=str)) == "WT").sum()),
            "time_points": sorted(int(value) for value in metadata.get("clock_hour", pd.Series(dtype=float)).dropna().unique()),
            "batches": sorted(metadata.get("batch", pd.Series(dtype=str)).dropna().unique().tolist()),
            "sexes": sorted(metadata.get("sex", pd.Series(dtype=str)).dropna().unique().tolist()),
        },
        "raw_tar_inventory": {"n_members": len(raw_tar_members), "members": raw_tar_members},
        "checks": checks,
        "fractional_value_count_total": total_fractional,
        "fractional_value_fraction_total": total_fractional / total_values,
        "primary_input_pass": primary_input_pass,
        "failure_reason": None if primary_input_pass else "G0_FAIL_STRICT_INTEGER_COUNT_REQUIREMENT",
    }
    (RUN / CFG["outputs"]["audit_json"]).write_text(json.dumps(audit, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps({"run": str(RUN.relative_to(ROOT)), "primary_input_pass": primary_input_pass, "fractional_value_count": total_fractional, "fractional_value_fraction": total_fractional / total_values}, ensure_ascii=False))


if __name__ == "__main__":
    main()
