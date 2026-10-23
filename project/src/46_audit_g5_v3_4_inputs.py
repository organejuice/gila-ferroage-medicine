#!/usr/bin/env python3
"""Audit all G5 v3.4 inputs before any model is fitted."""
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
CFG_PATH = ROOT / "configs/g5_v3.4_execution.yaml"
CFG = yaml.safe_load(CFG_PATH.read_text(encoding="utf-8"))
RUN = ROOT / CFG["outputs"]["result_dir"]
PRIMARY = CFG["inputs"]["primary"]
TOLERANCE = float(CFG["analysis"]["fractional_tolerance"])
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


def load_counts(path: Path) -> tuple[pd.DataFrame, list[str], pd.Series]:
    frame = pd.read_csv(path, compression="gzip", low_memory=False)
    if "gene_id" not in frame.columns:
        raise ValueError(f"gene_id column missing from {path}")
    samples = [column for column in frame.columns if column != "gene_id"]
    values = frame.loc[:, samples].apply(pd.to_numeric, errors="coerce")
    return values, samples, frame["gene_id"].astype(str)


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
            "time_bin": "AM" if hour < 12 else "PM",
            "sex": sex,
            "batch": f"batch{batch}",
        })
    return pd.DataFrame(rows)


def matrix_audit(label: str, values: pd.DataFrame, samples: list[str], gene_ids: pd.Series) -> dict:
    matrix = values.to_numpy(dtype=float, copy=False)
    finite = np.isfinite(matrix)
    fractional = finite & (np.abs(matrix - np.rint(matrix)) > TOLERANCE)
    return {
        "source_group": label,
        "n_genes": int(values.shape[0]),
        "n_samples": int(values.shape[1]),
        "gene_ids_unique": bool(gene_ids.is_unique),
        "sample_ids_unique": bool(len(samples) == len(set(samples))),
        "finite_values": bool(finite.all()),
        "nonnegative_values": bool(finite.all() and (matrix >= 0).all()),
        "fractional_value_count": int(fractional.sum()),
        "total_value_count": int(matrix.size),
        "fractional_value_fraction": float(fractional.mean()),
        "minimum_value": float(np.nanmin(matrix)),
        "maximum_value": float(np.nanmax(matrix)),
        "minimum_library_size": float(np.nansum(matrix, axis=0).min()),
        "maximum_library_size": float(np.nansum(matrix, axis=0).max()),
    }


def registered_files() -> list[tuple[str, Path, str]]:
    rows: list[tuple[str, Path, str]] = []
    for key, expected in PRIMARY["expected_sha256"].items():
        rows.append((f"GSE261698:{key}", ROOT / PRIMARY[key], expected))
    for accession, entry in CFG["inputs"]["validations"].items():
        for key, expected in entry["expected_sha256"].items():
            rows.append((f"{accession}:{key}", ROOT / entry[key], expected))
    return rows


def main() -> None:
    if RUN.exists():
        raise FileExistsError(f"Refusing to overwrite existing run directory: {RUN}")
    for directory in (RUN / "tables", RUN / "gates", RUN / "figures", RUN / "registries", RUN / "logs"):
        directory.mkdir(parents=True, exist_ok=True)

    file_rows = []
    for label, path, expected in registered_files():
        observed = sha256(path) if path.is_file() else None
        file_rows.append({
            "source": label, "file": path.relative_to(ROOT).as_posix(),
            "size_bytes": path.stat().st_size if path.is_file() else None,
            "expected_sha256": expected, "observed_sha256": observed,
            "exists": path.is_file(), "hash_match": observed == expected,
        })
    file_manifest = pd.DataFrame(file_rows)
    file_manifest.to_csv(RUN / "registries/source_file_manifest.tsv", sep="\t", index=False)

    app_values, app_samples, app_genes = load_counts(ROOT / PRIMARY["app_counts"])
    wt_values, wt_samples, wt_genes = load_counts(ROOT / PRIMARY["wt_counts"])
    app_audit = matrix_audit("APP", app_values, app_samples, app_genes)
    wt_audit = matrix_audit("WT", wt_values, wt_samples, wt_genes)
    metadata = parse_metadata(app_samples + wt_samples)

    with tarfile.open(ROOT / PRIMARY["geo_raw_tar"], "r") as archive:
        raw_members = archive.getnames()
    pd.DataFrame({"member": raw_members}).to_csv(
        RUN / "registries/GSE261698_RAW_tar_inventory.tsv", sep="\t", index=False
    )

    total_fractional = app_audit["fractional_value_count"] + wt_audit["fractional_value_count"]
    total_values = app_audit["total_value_count"] + wt_audit["total_value_count"]
    g4 = json.loads((ROOT / CFG["inputs"]["g4_decision"]).read_text(encoding="utf-8"))
    checks = {
        "upstream_G4_pass": g4.get("decision") in {"G4_PASS", "G4_PASS_MODEL_SCOPED"},
        "all_registered_files_exist": bool(file_manifest["exists"].all()),
        "all_registered_hashes_match": bool(file_manifest["hash_match"].all()),
        "expected_gene_rows": app_audit["n_genes"] == PRIMARY["expected_n_genes"] == wt_audit["n_genes"],
        "identical_gene_order": bool(app_genes.equals(wt_genes)),
        "unique_gene_ids": bool(app_audit["gene_ids_unique"] and wt_audit["gene_ids_unique"]),
        "expected_samples": len(app_samples) == PRIMARY["expected_n_model"] and len(wt_samples) == PRIMARY["expected_n_control"],
        "unique_sample_ids": len(set(app_samples + wt_samples)) == len(app_samples + wt_samples),
        "finite_nonnegative_count_like_values": bool(
            app_audit["finite_values"] and wt_audit["finite_values"]
            and app_audit["nonnegative_values"] and wt_audit["nonnegative_values"]
        ),
        "complete_sample_metadata": bool(metadata["metadata_parse_pass"].all()),
        "expected_disease_counts": bool(
            (metadata["disease"] == "APP_PS1").sum() == PRIMARY["expected_n_model"]
            and (metadata["disease"] == "WT").sum() == PRIMARY["expected_n_control"]
        ),
        "complete_24_hour_coverage": set(metadata["clock_hour"].tolist()) == set(TIME_TO_HOUR.values()),
        "fractional_values_preserved_not_rounded": total_fractional > 0,
    }
    required = [
        "upstream_G4_pass", "all_registered_files_exist", "all_registered_hashes_match",
        "expected_gene_rows", "identical_gene_order", "unique_gene_ids", "expected_samples",
        "unique_sample_ids", "finite_nonnegative_count_like_values", "complete_sample_metadata",
        "expected_disease_counts", "complete_24_hour_coverage", "fractional_values_preserved_not_rounded",
    ]
    primary_input_pass = all(checks[key] for key in required)

    metadata.to_csv(RUN / "tables/GSE261698_sample_metadata.tsv", sep="\t", index=False)
    gate_rows = [
        {"criterion": key, "required": key in required, "observed": value,
         "detail": "v3.4 accepts source-provided fractional count-like values" if key == "fractional_values_preserved_not_rounded" else ""}
        for key, value in checks.items()
    ]
    gate_rows.extend([
        {"criterion": "fractional_value_count", "required": False, "observed": total_fractional, "detail": "descriptive"},
        {"criterion": "fractional_value_fraction", "required": False, "observed": total_fractional / total_values, "detail": "descriptive"},
        {"criterion": "rounding_allowed", "required": True, "observed": False, "detail": "must remain false"},
    ])
    pd.DataFrame(gate_rows).to_csv(RUN / CFG["outputs"]["input_gate_table"], sep="\t", index=False)

    audit = {
        "run_id": CFG["run_id"], "dataset": "GSE261698", "raw_data_modified": False,
        "input_mode": "source_confirmed_nonnegative_expected_or_estimated_gene_counts",
        "rounding_performed": False, "DESeq2_used": False,
        "matrix_audits": {"APP": app_audit, "WT": wt_audit},
        "fractional_value_count_total": total_fractional,
        "fractional_value_fraction_total": total_fractional / total_values,
        "metadata_summary": {
            "n_samples": int(len(metadata)),
            "n_APP_PS1": int((metadata["disease"] == "APP_PS1").sum()),
            "n_WT": int((metadata["disease"] == "WT").sum()),
            "time_points": sorted(int(x) for x in metadata["clock_hour"].unique()),
            "batches": sorted(metadata["batch"].unique().tolist()),
            "sexes": sorted(metadata["sex"].unique().tolist()),
        },
        "checks": checks, "primary_input_pass": primary_input_pass,
        "failure_reason": None if primary_input_pass else "G0_V3_4_INPUT_ADMISSIBILITY_FAILURE",
    }
    (RUN / CFG["outputs"]["audit_json"]).write_text(
        json.dumps(audit, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps({
        "run": RUN.relative_to(ROOT).as_posix(), "primary_input_pass": primary_input_pass,
        "fractional_value_count": total_fractional,
        "fractional_value_fraction": total_fractional / total_values,
    }, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
