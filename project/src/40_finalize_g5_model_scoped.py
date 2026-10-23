#!/usr/bin/env python3
"""Finalize G5 under the frozen model-scoped rule: animal OR cell model."""
from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import yaml


ROOT = Path(__file__).resolve().parents[1]
CFG_PATH = ROOT / "configs" / "g5_model_scoped_v3.2_execution.yaml"
CFG = yaml.safe_load(CFG_PATH.read_text(encoding="utf-8"))
RUN = ROOT / "results" / CFG["run_id"]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def rel(path: Path) -> str:
    return path.relative_to(ROOT).as_posix()


def numeric(frame: pd.DataFrame, columns: list[str]) -> pd.DataFrame:
    for column in columns:
        frame[column] = pd.to_numeric(frame[column], errors="coerce")
    return frame


def main() -> None:
    RUN.joinpath("gates").mkdir(parents=True, exist_ok=True)
    g4 = json.loads((ROOT / CFG["inputs"]["g4_decision"]).read_text(encoding="utf-8"))
    cell_cfg_path = ROOT / "configs" / "g5_model_scoped_v3.2_cell_rescue_execution.yaml"
    cell_cfg = yaml.safe_load(cell_cfg_path.read_text(encoding="utf-8"))
    threshold = float(CFG["analysis"]["gate"]["primary_region_module_fdr_max"])
    minimum = int(CFG["analysis"]["gate"]["minimum_passing_modules"])
    primary_region = CFG["analysis"]["primary_region"]

    animal_path = RUN / "tables" / "GSE143758_animal_region_module_effects.tsv"
    animal = numeric(pd.read_csv(animal_path, sep="\t"), ["logFC", "SE", "P.Value", "FDR_within_stratum"])
    animal_primary = animal[animal["region"] == primary_region].copy()
    animal_primary["route"] = "animal_model"
    animal_primary["model_dataset"] = "GSE143758"
    animal_primary["gate_FDR"] = animal_primary["FDR_within_stratum"]
    animal_primary["pass_model_gate"] = animal_primary["gate_FDR"] < threshold

    cell_path = RUN / "tables" / "GSE319380_cell_module_effects.tsv"
    cell = numeric(pd.read_csv(cell_path, sep="\t"), ["logFC", "SE", "P.Value", "FDR_within_contrast"])
    cell["route"] = "cell_model"
    cell["model_dataset"] = "GSE319380"
    cell["gate_FDR"] = cell["FDR_within_contrast"]
    cell["pass_model_gate"] = cell["gate_FDR"] < threshold
    cell["region"] = "cell_model"

    module_gate = pd.concat([
        animal_primary[["route", "model_dataset", "region", "module", "n_present", "n_total", "coverage", "n_animals", "n_AD_model", "n_WT_control", "logFC", "SE", "P.Value", "gate_FDR", "pass_model_gate"]].rename(columns={"n_animals": "n_biological_units", "n_AD_model": "n_case_units", "n_WT_control": "n_control_units"}),
        cell[["route", "model_dataset", "region", "module", "n_present", "n_total", "coverage", "n_units", "n_sAD_lines", "n_control_lines", "logFC", "SE", "P.Value", "gate_FDR", "pass_model_gate"]].rename(columns={"n_units": "n_biological_units", "n_sAD_lines": "n_case_units", "n_control_lines": "n_control_units"}),
    ], ignore_index=True)
    module_gate.to_csv(RUN / "gates" / "G5_model_module_gate.tsv", sep="\t", index=False)

    animal_passing = animal_primary.loc[animal_primary["pass_model_gate"], "module"].tolist()
    cell_passing = cell.loc[cell["pass_model_gate"], "module"].tolist()
    upstream_pass = g4.get("decision") in {"G4_PASS_MODEL_SCOPED", "G4_PASS"}
    model_pass = bool(upstream_pass and (len(animal_passing) + len(cell_passing) >= minimum))
    decision = "G5_PASS_MODEL_SCOPED" if model_pass else "G5_NOT_YET_PASS"
    gate_rows = [
        {"criterion": "upstream_G4_pass", "required": True, "observed": upstream_pass, "detail": g4.get("decision", "missing")},
        {"criterion": "animal_model_primary_hippocampus", "required": True, "observed": len(animal_passing) >= minimum, "detail": f"GSE143758; passing={';'.join(animal_passing) if animal_passing else 'none'}"},
        {"criterion": "cell_model_rescue", "required": True, "observed": len(cell_passing) >= minimum, "detail": f"GSE319380; passing={';'.join(cell_passing) if cell_passing else 'none'}"},
        {"criterion": "animal_or_cell_model_gate", "required": True, "observed": model_pass, "detail": "at least one predeclared module with FDR<0.05 in either qualified model route"},
        {"criterion": "human_translation", "required": False, "observed": False, "detail": "HUMAN_TRANSLATION_UNRESOLVED; not a hard gate"},
    ]
    pd.DataFrame(gate_rows).to_csv(RUN / "gates" / "G5_model_gate_v3.2.tsv", sep="\t", index=False)

    decision_record = {
        "run_id": CFG["run_id"],
        "protocol": "AD disease-side v3.2 model-scoped correction",
        "decision": decision,
        "upstream_G4": g4.get("decision"),
        "human_translation": "HUMAN_TRANSLATION_UNRESOLVED",
        "model_route": "animal_or_cell_model_sufficient",
        "animal_route": {"dataset": "GSE143758", "primary_region": primary_region, "passing_modules": animal_passing},
        "cell_route": {"dataset": "GSE319380", "passing_modules": cell_passing},
        "passing_modules": sorted(set(animal_passing + cell_passing)),
        "claim_boundary": {
            "allowed": "only the tested animal or iPSC-astrocyte model response",
            "not_allowed": "human tissue confirmation, human spatial causality, or microglia-to-astrocyte causality",
        },
        "downstream": "G6 may proceed only as model-scoped" if model_pass else "G6 remains locked",
        "created_at": datetime.now(timezone.utc).isoformat(),
    }
    (RUN / "g5_model_decision.json").write_text(json.dumps(decision_record, ensure_ascii=False, indent=2), encoding="utf-8")

    report = f"""# Disease-side v3.2 G5 model-scoped execution

## Decision

**{decision}**

The frozen model-scoped rule was applied: a qualified animal **or** cell model
can satisfy G5. Human D3/D7/D8 results remain supportive translation context
and are labeled `HUMAN_TRANSLATION_UNRESOLVED`.

### Animal route

- Dataset: GSE143758 5xFAD versus WT source-defined astrocytes.
- Biological unit: animal/sample, not nucleus.
- Primary region: hippocampus; passing modules: {', '.join(animal_passing) if animal_passing else 'none'}.

### Cell route

- Dataset: GSE319380 sAD versus control iPSC-derived astrocytes.
- Biological unit: cell line × stressor after summing three technical replicates.
- Model: disease adjusted for stressor, blocked by cell line; passing modules: {', '.join(cell_passing) if cell_passing else 'none'}.

### Boundary

The model-scoped decision does not establish human tissue confirmation, human
spatial causality, or a causal microglia-to-astrocyte direction. G6 is
{'available only for model-scoped execution' if model_pass else 'locked because G5 did not pass'}.
"""
    (RUN / "execution_report.md").write_text(report, encoding="utf-8")

    manifest_values = list(CFG["inputs"].values()) + list(cell_cfg["inputs"].values()) + [
        "configs/g5_model_scoped_v3.2_execution.yaml",
        "configs/g5_model_scoped_v3.2_cell_rescue_execution.yaml",
        "src/36_prepare_gse143758_astrocyte_pseudobulk.py",
        "src/37_run_g5_animal_model.R",
        "src/39_run_g5_cell_model.R",
        "src/40_finalize_g5_model_scoped.py",
    ]
    manifest_inputs = []
    for value in manifest_values:
        path = ROOT / value
        if path.exists():
            manifest_inputs.append({"file": rel(path), "sha256": sha256(path), "size": path.stat().st_size})
    manifest = {"run_id": CFG["run_id"], "raw_data_modified": False, "remote_computation": False, "result": decision, "inputs": manifest_inputs}
    (RUN / "run_manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")
    checksum_lines = []
    for path in sorted(RUN.rglob("*")):
        if path.is_file() and path.name != "checksums.sha256":
            checksum_lines.append(f"{sha256(path)}  {rel(path)}")
    (RUN / "checksums.sha256").write_text("\n".join(checksum_lines) + "\n", encoding="utf-8")
    print(json.dumps({"run": rel(RUN), "decision": decision, "animal_passing": animal_passing, "cell_passing": cell_passing}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
