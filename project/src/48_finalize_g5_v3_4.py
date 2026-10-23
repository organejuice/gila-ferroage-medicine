#!/usr/bin/env python3
"""Finalize and audit the frozen AD disease-side G5 v3.4 run."""
from __future__ import annotations

import hashlib
import json
import platform
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import yaml


ROOT = Path(__file__).resolve().parents[1]
CFG = yaml.safe_load((ROOT / "configs/g5_v3.4_execution.yaml").read_text(encoding="utf-8"))
RUN = ROOT / CFG["outputs"]["result_dir"]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def rel(path: Path) -> str:
    return path.relative_to(ROOT).as_posix()


def read_table(name: str) -> pd.DataFrame:
    path = RUN / "tables" / name
    if not path.is_file():
        raise FileNotFoundError(path)
    return pd.read_csv(path, sep="\t")


def numeric(frame: pd.DataFrame, columns: list[str]) -> pd.DataFrame:
    for column in columns:
        if column in frame:
            frame[column] = pd.to_numeric(frame[column], errors="coerce")
    return frame


def file_record(path: Path) -> dict:
    return {"file": rel(path), "sha256": sha256(path), "size_bytes": path.stat().st_size}


def runtime_version() -> str:
    try:
        result = subprocess.run([r"D:\R\R-4.3.3\bin\Rscript.exe", "--version"], capture_output=True, text=True, check=False)
        return (result.stdout + result.stderr).strip().splitlines()[-1]
    except Exception as exc:  # pragma: no cover - reporting fallback
        return f"unavailable: {exc}"


def bh_fixed(values: pd.Series, n_tests: int = 7) -> pd.Series:
    """BH correction with a fixed family size, matching R p.adjust(..., n=7)."""
    result = pd.Series(np.nan, index=values.index, dtype=float)
    valid = values.dropna()
    if valid.empty:
        return result
    order = np.argsort(valid.to_numpy())
    sorted_p = valid.to_numpy()[order]
    ranks = np.arange(1, len(sorted_p) + 1, dtype=float)
    adjusted = np.minimum.accumulate((sorted_p * n_tests / ranks)[::-1])[::-1]
    adjusted = np.minimum(adjusted, 1.0)
    ordered_index = valid.index.to_numpy()[order]
    result.loc[ordered_index] = adjusted
    return result


def main() -> None:
    decision_path = RUN / "g5_v3.4_decision.json"
    audit = json.loads((RUN / CFG["outputs"]["audit_json"]).read_text(encoding="utf-8"))
    g4 = json.loads((ROOT / CFG["inputs"]["g4_decision"]).read_text(encoding="utf-8"))
    primary = numeric(read_table("GSE261698_primary_module_effects.tsv"), ["logFC", "SE", "P.Value", "FDR_BH_7", "CI_low", "CI_high", "n_present", "n_total", "coverage"])
    expected_modules = [
        "M3_SASP_Senescence", "M5_Iron_Uptake", "M5_Iron_Storage", "M5_Iron_Export",
        "M5_Ferritinophagy", "M6_LPO_ACSL4", "M7_FerroAging_Author",
    ]
    if sorted(primary["module"].tolist()) != sorted(expected_modules):
        raise ValueError("Primary result does not contain exactly the seven frozen G5 modules")
    if primary["module"].duplicated().any():
        raise ValueError("Duplicate primary module rows")

    primary["FDR_recomputed"] = bh_fixed(primary["P.Value"], n_tests=7)
    if not primary["FDR_BH_7"].fillna(-1).round(12).equals(primary["FDR_recomputed"].fillna(-1).round(12)):
        raise ValueError("Primary FDR_BH_7 does not match recomputed BH correction")

    upstream_pass = g4.get("decision") in {"G4_PASS", "G4_PASS_MODEL_SCOPED"}
    primary_pass = bool(((primary["status"] == "estimable") & (primary["FDR_BH_7"] < 0.05)).any())
    passing_modules = primary.loc[(primary["status"] == "estimable") & (primary["FDR_BH_7"] < 0.05), "module"].tolist()
    decision = "G5_PASS_MODEL_SCOPED" if upstream_pass and primary_pass else "G5_NOT_YET_PASS"

    gate_rows = [
        {"criterion": "G0_v3_4_primary_input_admissibility", "required": True, "observed": bool(audit.get("primary_input_pass")), "detail": "source-confirmed finite nonnegative count-like matrices; fractional values preserved"},
        {"criterion": "upstream_G4_pass", "required": True, "observed": upstream_pass, "detail": g4.get("decision", "missing")},
        {"criterion": "primary_has_exactly_seven_modules", "required": True, "observed": len(primary) == 7, "detail": "BH family is exactly the seven frozen modules"},
        {"criterion": "primary_module_FDR_lt_0_05", "required": True, "observed": primary_pass, "detail": ";".join(passing_modules) if passing_modules else "none"},
        {"criterion": "validation_can_rescue_primary_failure", "required": False, "observed": False, "detail": "frozen rule: false"},
        {"criterion": "G5_decision", "required": True, "observed": decision == "G5_PASS_MODEL_SCOPED", "detail": decision},
    ]
    pd.DataFrame(gate_rows).to_csv(RUN / "gates/G5_v3.4_gate.tsv", sep="\t", index=False)

    validation_frames = []
    for dataset, filename, role in [
        ("GSE291239", "GSE291239_module_effects.tsv", "independent_animal_validation"),
        ("GSE74614", "GSE74614_module_effects.tsv", "orthogonal_microarray_validation"),
        ("GSE290101", "GSE290101_module_effects.tsv", "animal_aware_repeated_measure_sensitivity"),
    ]:
        frame = numeric(read_table(filename), ["logFC", "SE", "P.Value", "FDR_BH_7", "CI_low", "CI_high", "n_present", "n_total", "coverage"])
        for column in ("n_present", "n_total", "coverage"):
            mapped = f"{column}_map"
            if mapped in frame.columns:
                frame[column] = frame[column].fillna(frame[mapped])
        frame["role"] = role
        frame["primary_direction_concordant"] = frame["module"].map(primary.set_index("module")["logFC"]).mul(frame["logFC"]) > 0
        validation_frames.append(frame[["dataset", "role", "module", "n_present", "n_total", "coverage", "n_units", "n_case", "n_control", "logFC", "SE", "P.Value", "FDR_BH_7", "primary_direction_concordant", "status", "note"]])
    validation_summary = pd.concat(validation_frames, ignore_index=True)
    validation_summary.to_csv(RUN / "tables/G5_validation_summary.tsv", sep="\t", index=False)

    sensitivity = numeric(read_table("GSE261698_primary_sensitivities.tsv"), ["logFC", "SE", "P.Value", "FDR_BH_7"])
    sensitivity_summary = (
        sensitivity.groupby("analysis", dropna=False)
        .agg(n_modules=("module", "size"), n_estimable=("status", lambda x: int((x == "estimable").sum())),
             min_FDR=("FDR_BH_7", "min"), n_FDR_lt_0_05=("FDR_BH_7", lambda x: int((x < 0.05).sum())))
        .reset_index()
    )
    sensitivity_summary.to_csv(RUN / "tables/G5_sensitivity_summary.tsv", sep="\t", index=False)

    primary_view = primary[["module", "n_present", "n_total", "coverage", "n_units", "n_case", "n_control", "logFC", "SE", "CI_low", "CI_high", "P.Value", "FDR_BH_7", "status", "note"]]
    decision_record = {
        "run_id": CFG["run_id"],
        "protocol": CFG["protocol"],
        "protocol_checksum_file": CFG["protocol_checksum"],
        "decision": decision,
        "upstream_G4": g4.get("decision"),
        "primary_dataset": "GSE261698",
        "primary_biological_unit": "animal",
        "primary_formula": CFG["analysis"]["primary_formula"],
        "primary_coefficient": CFG["analysis"]["primary_coefficient"],
        "primary_module_family_size": 7,
        "passing_modules": passing_modules,
        "primary_effects": primary_view.to_dict(orient="records"),
        "validation_roles": {"GSE291239": "independent_animal_validation", "GSE74614": "orthogonal_microarray_validation", "GSE290101": "animal_aware_repeated_measure_sensitivity"},
        "validation_can_rescue_primary_failure": False,
        "claim_boundary": {"allowed": "the tested APP/PS1-21 animal astrocyte response and model-scoped replication context", "not_allowed": "direct human confirmation, spatial causality, or causal microglia-to-astrocyte communication"},
        "downstream": "G6 may proceed only within the model-scoped boundary" if decision == "G5_PASS_MODEL_SCOPED" else "G6 remains locked",
        "runtime": {"python": platform.python_version(), "R": runtime_version(), "remote_computation": False},
        "created_at": datetime.now(timezone.utc).isoformat(),
    }
    decision_path.write_text(json.dumps(decision_record, ensure_ascii=False, indent=2), encoding="utf-8")

    report_lines = [
        "# Disease-side v3.4 G5 execution report",
        "",
        f"## Decision: **{decision}**",
        "",
        "The frozen v3.4 primary route was applied without rounding or integer coercion. The primary BH family contains exactly seven frozen modules; validation datasets are supportive and non-rescuing.",
        "",
        "## Primary result",
        "",
        "| Module | n present/total | logFC | 95% CI | P | BH-FDR | Status |",
        "|---|---:|---:|---|---:|---:|---|",
    ]
    for _, row in primary.iterrows():
        report_lines.append(f"| {row['module']} | {row['n_present']:.0f}/{row['n_total']:.0f} | {row['logFC']:.6g} | [{row['CI_low']:.6g}, {row['CI_high']:.6g}] | {row['P.Value']:.6g} | {row['FDR_BH_7']:.6g} | {row['status']} |")
    report_lines += [
        "",
        f"Passing primary module(s): {', '.join(passing_modules) if passing_modules else 'none'}.",
        "",
        "## Validation and sensitivity context",
        "",
        "- GSE291239: independent AppNL-G-F versus WT astrocyte animal validation, with sex adjustment when estimable.",
        "- GSE74614: GFAP+/+ Vim+/+ cortical astrocyte microarray validation; probe-to-gene mapping was parsed from the registered GPL10333 platform table embedded in the family SOFT file.",
        "- GSE290101: animal-blocked repeated-measures sensitivity across STAT/NFKB/DOUBLE reporter-defined astrocyte libraries.",
        f"- Primary sensitivity package: {len(sensitivity)} module-analysis rows, including circular-time interaction, AM/PM, leave-one-time-point-out, leave-one-batch-out, sex-stratified, leave-one-animal-out, and rounded-input non-primary sensitivity.",
        "",
        "## Claim boundary",
        "",
        "This result is model-scoped. It does not establish direct human tissue confirmation, spatial communication, or causal microglia-to-astrocyte direction. G6 remains governed by its own spatial AD/control donor-replicated gate.",
        "",
        "## Reproducibility",
        "",
        f"- Run ID: `{CFG['run_id']}`",
        f"- R runtime: `{runtime_version()}`",
        "- Raw files were not modified; registered source hashes were checked by the G0 audit.",
        "- Historical v3.2/v3.3 result directories were not overwritten.",
    ]
    (RUN / "execution_report.md").write_text("\n".join(report_lines) + "\n", encoding="utf-8")

    manifest_paths = [
        ROOT / "docs/protocol/disease_side_v3.4_latest_frozen.md",
        ROOT / "docs/protocol/disease_side_v3.4_freeze.sha256",
        ROOT / "configs/disease_side_v3.4_frozen.yaml",
        ROOT / "configs/g5_v3.4_execution.yaml",
        ROOT / "docs/decisions/2026-09-21-disease-side-v3.4-noninteger-count-amendment.md",
        ROOT / "docs/decisions/2026-09-22-disease-side-v3.4-g5-execution.md",
        ROOT / "src/46_audit_g5_v3_4_inputs.py",
        ROOT / "src/47_run_g5_v3_4.R",
        ROOT / "src/48_finalize_g5_v3_4.py",
        ROOT / CFG["inputs"]["g4_decision"],
        ROOT / CFG["inputs"]["m7_reference"],
    ]
    for entry in CFG["inputs"]["primary"]["expected_sha256"]:
        manifest_paths.append(ROOT / CFG["inputs"]["primary"][entry])
    for dataset in CFG["inputs"]["validations"].values():
        for key in dataset:
            if key != "expected_sha256":
                manifest_paths.append(ROOT / dataset[key])
    manifest = {
        "run_id": CFG["run_id"], "decision": decision, "raw_data_modified": False,
        "remote_computation": False, "inputs": [file_record(p) for p in manifest_paths if p.is_file()],
    }
    (RUN / "run_manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")

    checksum_lines = []
    for path in sorted(RUN.rglob("*")):
        if path.is_file() and path.name != "checksums.sha256":
            checksum_lines.append(f"{sha256(path)}  {rel(path)}")
    (RUN / "checksums.sha256").write_text("\n".join(checksum_lines) + "\n", encoding="utf-8")
    print(json.dumps({"run": rel(RUN), "decision": decision, "passing_modules": passing_modules, "primary_fdr_min": float(primary["FDR_BH_7"].min(skipna=True))}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
