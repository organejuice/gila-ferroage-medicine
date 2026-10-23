#!/usr/bin/env python3
"""Finalize G5 v3.3 when a frozen primary-input hard gate fails."""
from __future__ import annotations

import hashlib
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import yaml


ROOT = Path(__file__).resolve().parents[1]
CFG = yaml.safe_load((ROOT / "configs/g5_v3.3_execution.yaml").read_text(encoding="utf-8"))
RUN = ROOT / CFG["outputs"]["result_dir"]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def rel(path: Path) -> str:
    return path.relative_to(ROOT).as_posix()


def git_revision() -> str | None:
    try:
        return subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    except (OSError, subprocess.CalledProcessError):
        return None


def main() -> None:
    audit_path = RUN / CFG["outputs"]["audit_json"]
    if not audit_path.is_file():
        raise FileNotFoundError(f"Required input audit is missing: {audit_path}")
    audit = json.loads(audit_path.read_text(encoding="utf-8"))
    g4 = json.loads((ROOT / CFG["inputs"]["g4_decision"]).read_text(encoding="utf-8"))
    upstream_pass = g4.get("decision") in {"G4_PASS", "G4_PASS_MODEL_SCOPED"}
    if audit["primary_input_pass"]:
        raise RuntimeError("Primary input passed; run the preregistered G5 model rather than this input-stop finalizer.")

    fractional = int(audit["fractional_value_count_total"])
    fraction = float(audit["fractional_value_fraction_total"])
    raw_tar_members = int(audit["raw_tar_inventory"]["n_members"])
    gates = pd.DataFrame([
        {"criterion": "upstream_G4_pass", "required": True, "observed": upstream_pass, "detail": g4.get("decision", "missing")},
        {"criterion": "source_defined_astrocyte_animal_endpoint", "required": True, "observed": True, "detail": "GSE261698 APP/PS1-21 and WT astrocyte TRAP-RNA-seq"},
        {"criterion": "animal_metadata_and_24_hour_design", "required": True, "observed": audit["checks"]["complete_sample_metadata"] and audit["checks"]["expected_disease_counts"], "detail": "48 distinct animals; metadata parsed from GEO count-matrix column names"},
        {"criterion": "strict_integer_processed_counts", "required": True, "observed": audit["checks"]["strict_integer_counts"], "detail": f"{fractional} fractional values of {sum(row['total_value_count'] for row in audit['matrix_audits'].values())} ({fraction:.6%}); frozen maximum is 0"},
        {"criterion": "official_RAW_tar_complete_alternative", "required": False, "observed": audit["checks"]["complete_raw_matrix_alternative"], "detail": f"RAW tar contains {raw_tar_members} members, not a complete 48-animal matrix"},
        {"criterion": "primary_module_model", "required": True, "observed": False, "detail": "not run because the strict-integer primary input gate failed"},
        {"criterion": "validation_datasets", "required": False, "observed": False, "detail": "not downloaded or analyzed; preregistered validation cannot rescue a failed primary input gate"},
    ])
    gate_path = RUN / CFG["outputs"]["input_gate_table"]
    gates.to_csv(gate_path, sep="\t", index=False)

    decision = "G5_NOT_YET_PASS"
    record = {
        "run_id": CFG["run_id"],
        "protocol": "AD disease-side v3.3",
        "decision": decision,
        "stop_stage": "G0_primary_input_admissibility",
        "upstream_G4": g4.get("decision"),
        "primary_dataset": "GSE261698",
        "reason": "The downloaded supplementary matrices are nonnegative and sample-mapped, but contain fractional values while the frozen v3.3 primary input gate requires strictly integer processed counts.",
        "fractional_value_count": fractional,
        "fractional_value_fraction": fraction,
        "raw_tar_alternative": {"member_count": raw_tar_members, "complete_48_animal_matrix": audit["checks"]["complete_raw_matrix_alternative"]},
        "model_statistics_executed": False,
        "validation_statistics_executed": False,
        "human_translation": "SUPPORTIVE_NOT_HARD_GATE",
        "downstream": "G6 remains locked",
        "claim_boundary": {
            "allowed": "GSE261698 supplementary-matrix admissibility audit only",
            "not_allowed": "a conclusion about astrocyte response, module direction, or a failed biological mechanism",
        },
        "created_at": datetime.now(timezone.utc).isoformat(),
    }
    decision_path = RUN / "g5_v3.3_decision.json"
    decision_path.write_text(json.dumps(record, ensure_ascii=False, indent=2), encoding="utf-8")

    report = f"""# Disease-side v3.3 G5 execution

## Decision

**{decision}** at `G0_primary_input_admissibility`.

GSE261698's APP and WT astrocyte supplementary matrices have matching gene order,
48 parsable animal-level samples, full two-hour circadian coverage, and nonnegative
finite values. However, they contain {fractional:,} fractional values
({fraction:.6%} of all matrix values). The frozen v3.3 protocol requires strict
integer processed counts for the primary G5 gate, with no noninteger values allowed.

The GEO `RAW.tar` archive contains only {raw_tar_members} MTX/TSV members and does
not supply a complete 48-animal alternative matrix. No rounding, threshold change,
or post hoc acceptance rule was applied.

## Analyses not run

The edgeR TMM/limma-voom primary module model and all GSE291239, GSE74614, and
GSE290101 validation analyses were not run. Under v3.3, validation and time-bin
effects cannot rescue a failed primary input gate.

## Downstream status

G4 remains `{g4.get('decision')}`. G6 remains locked. This is an input-admissibility
result, not evidence for or against the astrocyte-response mechanism.
"""
    (RUN / "execution_report.md").write_text(report, encoding="utf-8")

    manifest_paths = [
        "configs/g5_v3.3_execution.yaml", "configs/disease_side_v3.3_frozen.yaml",
        "docs/protocol/disease_side_v3.3_latest_frozen.md", "docs/protocol/disease_side_v3.3_freeze.sha256",
        "src/44_audit_g5_v3_3_primary_input.py", "src/45_finalize_g5_v3_3_input_stop.py",
        "data/raw/GSE261698.dvc",
        CFG["inputs"]["g4_decision"],
        CFG["inputs"]["primary"]["app_counts"], CFG["inputs"]["primary"]["wt_counts"],
        CFG["inputs"]["primary"]["family_soft"], CFG["inputs"]["primary"]["geo_raw_tar"],
    ]
    manifest = {
        "run_id": CFG["run_id"], "result": decision, "stop_stage": "G0_primary_input_admissibility",
        "raw_data_modified": False, "remote_computation": False, "git_revision": git_revision(),
        "inputs": [{"file": value, "sha256": sha256(ROOT / value), "size": (ROOT / value).stat().st_size} for value in manifest_paths],
    }
    (RUN / "run_manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")
    checksum_lines = [f"{sha256(path)}  {rel(path)}" for path in sorted(RUN.rglob("*")) if path.is_file() and path.name != "checksums.sha256"]
    (RUN / "checksums.sha256").write_text("\n".join(checksum_lines) + "\n", encoding="utf-8")
    print(json.dumps({"run": rel(RUN), "decision": decision, "stop_stage": record["stop_stage"]}, ensure_ascii=False))


if __name__ == "__main__":
    main()
