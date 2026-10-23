#!/usr/bin/env python3
"""Finalize the model-scoped G5 decision without treating human G5 as hard."""
from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
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


def main() -> None:
    g4 = json.loads((ROOT / CFG["inputs"]["g4_decision"]).read_text(encoding="utf-8"))
    effects_path = RUN / "tables" / "GSE143758_animal_region_module_effects.tsv"
    effects = pd.read_csv(effects_path, sep="\t")
    effects["logFC"] = pd.to_numeric(effects["logFC"], errors="coerce")
    effects["SE"] = pd.to_numeric(effects["SE"], errors="coerce")
    effects["P.Value"] = pd.to_numeric(effects["P.Value"], errors="coerce")
    effects["FDR_within_stratum"] = pd.to_numeric(effects["FDR_within_stratum"], errors="coerce")

    primary_region = CFG["analysis"]["primary_region"]
    primary = effects[effects["region"] == primary_region].copy()
    sensitivity = effects[effects["region"] != primary_region].copy()
    gate_rows: list[dict] = []
    modules = sorted(primary["module"].unique().tolist())
    module_rows: list[dict] = []
    fdr_limit = float(CFG["analysis"]["gate"]["primary_region_module_fdr_max"])
    for module in modules:
        p = primary[primary["module"] == module].iloc[0]
        s = sensitivity[sensitivity["module"] == module]
        support_direction = bool(len(s) == 0 or (np.sign(s["logFC"]) == np.sign(p["logFC"])).all())
        module_rows.append({
            "module": module,
            "primary_region": primary_region,
            "primary_beta": float(p["logFC"]),
            "primary_SE": float(p["SE"]),
            "primary_P": float(p["P.Value"]),
            "primary_FDR": float(p["FDR_within_stratum"]),
            "sensitivity_regions": ";".join(sorted(s["region"].astype(str).unique().tolist())),
            "sensitivity_direction_concordant": support_direction,
            "pass_primary_model_gate": bool(float(p["FDR_within_stratum"]) < fdr_limit),
        })
    module_gate = pd.DataFrame(module_rows)
    module_gate.to_csv(RUN / "gates" / "G5_model_module_gate.tsv", sep="\t", index=False)

    upstream_pass = g4.get("decision") in {"G4_PASS_MODEL_SCOPED", "G4_PASS"}
    passing = module_gate.loc[module_gate["pass_primary_model_gate"], "module"].tolist()
    model_pass = bool(upstream_pass and len(passing) >= int(CFG["analysis"]["gate"]["minimum_passing_modules"]))
    decision = "G5_PASS_MODEL_SCOPED" if model_pass else "G5_NOT_YET_PASS"
    gate_rows.extend([
        {"criterion": "upstream_G4_pass", "required": True, "observed": upstream_pass, "detail": g4.get("decision", "missing")},
        {"criterion": "animal_model_astrocyte_endpoint", "required": True, "observed": True, "detail": "GSE143758 source-defined astrocyte matrices aggregated by animal/sample"},
        {"criterion": "primary_hippocampus_module_FDR", "required": True, "observed": model_pass, "detail": ";".join(f"{m}={m in passing}" for m in modules)},
        {"criterion": "cortex_sensitivity", "required": False, "observed": bool(len(sensitivity) > 0), "detail": "supportive only; not a hard gate"},
        {"criterion": "human_translation", "required": False, "observed": False, "detail": "HUMAN_TRANSLATION_UNRESOLVED; not a hard gate"},
    ])
    pd.DataFrame(gate_rows).to_csv(RUN / "gates" / "G5_model_gate_v3.2.tsv", sep="\t", index=False)

    decision_record = {
        "run_id": CFG["run_id"],
        "protocol": "AD disease-side v3.2 model-scoped correction",
        "decision": decision,
        "upstream_G4": g4.get("decision"),
        "human_translation": "HUMAN_TRANSLATION_UNRESOLVED",
        "model_route": "animal_or_cell_model_sufficient",
        "primary_model": "GSE143758 5xFAD versus WT astrocytes",
        "passing_modules": passing,
        "claim_boundary": {
            "allowed": "GSE143758 animal-model astrocyte response in the analyzed source-defined regions",
            "not_allowed": "human astrocyte confirmation, human spatial causality, or microglia-to-astrocyte causality",
        },
        "downstream": "G6 may proceed only as model-scoped" if model_pass else "G6 remains locked",
        "created_at": datetime.now(timezone.utc).isoformat(),
    }
    (RUN / "g5_model_decision.json").write_text(json.dumps(decision_record, ensure_ascii=False, indent=2), encoding="utf-8")
    report = f"""# Disease-side v3.2 G5 model-scoped execution

## Decision

**{decision}**

This execution applies the established model-scoped rule: a qualified animal
or cell model can satisfy G5; the human D3/D7/D8 result is supportive
translation context and is labeled `HUMAN_TRANSLATION_UNRESOLVED`.

The primary model was GSE143758 (5xFAD versus WT source-defined astrocytes),
with hippocampus as the prespecified primary region and cortex as sensitivity
support. The biological unit was the animal/sample, not the nucleus.

- Passing modules: {', '.join(passing) if passing else 'none'}
- Human translation: `HUMAN_TRANSLATION_UNRESOLVED`
- G6: {'model-scoped route available' if model_pass else 'locked'}

## Boundary

`{decision}` supports only the model-level astrocyte response tested here. It
does not establish a human astrocyte response, spatial causality, or a causal
microglia-to-astrocyte direction.
"""
    (RUN / "execution_report.md").write_text(report, encoding="utf-8")

    manifest_inputs = []
    for value in list(CFG["inputs"].values()) + ["src/36_prepare_gse143758_astrocyte_pseudobulk.py", "src/37_run_g5_animal_model.R", "src/38_finalize_g5_model_scoped.py", "configs/g5_model_scoped_v3.2_execution.yaml"]:
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
    print(json.dumps({"run": rel(RUN), "decision": decision, "passing_modules": passing}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
