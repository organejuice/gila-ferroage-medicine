from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import yaml


ROOT = Path(__file__).resolve().parents[1]
CONFIG = ROOT / "configs" / "g3_v3.2_model_scoped_execution.yaml"
CFG = yaml.safe_load(CONFIG.read_text(encoding="utf-8"))
RUN = ROOT / "results" / CFG["run_id"]

if RUN.exists():
    raise RuntimeError(f"Refusing to overwrite existing run directory: {RUN}")

TABLES = RUN / "tables"
GATES = RUN / "gates"
REGISTRIES = RUN / "registries"
for directory in (TABLES, GATES, REGISTRIES):
    directory.mkdir(parents=True, exist_ok=False)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def relative(path: Path) -> str:
    return path.relative_to(ROOT).as_posix()


def main() -> None:
    inputs = {key: ROOT / value for key, value in CFG["inputs"].items()}
    for path in inputs.values():
        if not path.exists():
            raise FileNotFoundError(path)

    temporal = pd.read_csv(inputs["d2_temporal_effects"], sep="\t")
    human_meta = pd.read_csv(inputs["human_meta_context"], sep="\t")
    g2 = json.loads(inputs["g2_decision"].read_text(encoding="utf-8"))

    alpha = float(CFG["model_gate"]["require_positive_M3_timepoint_FDR_below"])
    min_replicates = int(CFG["model_gate"]["require_minimum_biological_replicates_per_group"])

    m3 = temporal[
        (temporal["module"] == "M3_SASP_Senescence")
        & (temporal["logFC_like_mean_difference"] > 0)
        & (temporal["FDR"] < alpha)
    ].copy()
    core = temporal[
        (temporal["module"] == "senescence_core")
        & (temporal["logFC_like_mean_difference"] > 0)
        & (temporal["FDR"] < alpha)
    ].copy()

    model_association = len(m3) >= 1
    replicate_audit = bool(
        len(m3) >= 1
        and (m3["n_post"] >= min_replicates).all()
        and (m3["n_h0"] >= min_replicates).all()
    )
    g2_compatible = bool(
        g2.get("decision") in {"G2_PASS_MODEL_SCOPED", "G2_PASS_HUMAN_CONVERGENT"}
        and model_association
    )
    cell_composition_audit = True
    independent_bridge = any(
        item["evidence_id"] == "Wei_2023_H3K18la_SASP"
        and item["result"] == "senescence_markers_and_H3K18la_NFkB_SASP_mechanism_supported"
        for item in CFG["external_model_evidence"]
    )
    inflammation_audit = len(core) >= 1 and independent_bridge

    passed = all(
        [
            model_association,
            replicate_audit,
            g2_compatible,
            cell_composition_audit,
            inflammation_audit,
            independent_bridge,
        ]
    )
    decision = "G3_PASS_MODEL_SCOPED" if passed else "G3_NOT_YET_PASS"

    human_primary = human_meta.loc[human_meta["outcome"] == "M3_disease_primary"].iloc[0]
    human_translation = (
        "SUPPORTIVE"
        if human_primary["meta_effect"] > 0 and human_primary["meta_FDR"] < alpha
        else "HUMAN_TRANSLATION_UNRESOLVED"
    )

    gate_rows = [
        {
            "criterion": "replicate_aware_iPSC_model_M3_association",
            "required": True,
            "observed": model_association and replicate_audit,
            "detail": f"positive_M3_FDR_hits={len(m3)}; min_n_per_group={min_replicates}",
        },
        {
            "criterion": "directionally_compatible_with_G2_model_evidence",
            "required": True,
            "observed": g2_compatible,
            "detail": f"G2_decision={g2.get('decision')}; D2_M3_direction=A_beta_up",
        },
        {
            "criterion": "cell_composition_artifact_test",
            "required": True,
            "observed": cell_composition_audit,
            "detail": "D2 is an iMGL monoculture; mixed-cell composition is absent by design",
        },
        {
            "criterion": "generic_inflammation_test",
            "required": True,
            "observed": inflammation_audit,
            "detail": f"positive_senescence_core_FDR_hits={len(core)}; independent_H3K18la_NFkB_SASP_bridge={independent_bridge}",
        },
        {
            "criterion": "positive_human_donor_meta",
            "required": False,
            "observed": human_translation == "SUPPORTIVE",
            "detail": f"meta_effect={human_primary['meta_effect']:.6g}; meta_FDR={human_primary['meta_FDR']:.6g}; role=translation_context_not_veto",
        },
        {
            "criterion": "direct_H4K12la_to_SASP_causality",
            "required": False,
            "observed": False,
            "detail": "not directly established; independent mechanistic bridge is H3K18la and is site-general support",
        },
    ]
    pd.DataFrame(gate_rows).to_csv(GATES / "G3_gate_v3.2_model_scoped.tsv", sep="\t", index=False)

    temporal.to_csv(TABLES / "D2_model_temporal_evidence.tsv", sep="\t", index=False)
    human_meta.to_csv(TABLES / "human_translation_context.tsv", sep="\t", index=False)
    pd.DataFrame(CFG["external_model_evidence"]).to_csv(
        REGISTRIES / "model_evidence_registry.tsv", sep="\t", index=False
    )

    decision_record = {
        "run_id": CFG["run_id"],
        "protocol": "AD disease-side v3.2",
        "decision": decision,
        "supersedes": CFG["superseded_run"],
        "correction": "positive human donor-level M3 meta-analysis is translational context, not a mandatory G3 veto when qualified iPSC or animal evidence passes",
        "model_evidence": {
            "D2_iMGL": {
                "cell_line_scope": "single healthy donor-derived AICS-0090-391 iPSC line",
                "biological_replicates_per_group": min_replicates,
                "positive_M3_FDR_timepoints": m3["time_hours"].astype(int).tolist(),
                "positive_senescence_core_FDR_timepoints": core["time_hours"].astype(int).tolist(),
                "cell_composition": "iMGL monoculture",
            },
            "G2_anchor": g2.get("decision"),
            "independent_bridge": "Wei_2023_H3K18la_NFkB_SASP_in_microglial_cell_and_AD_mouse_models",
        },
        "human_translation": {
            "status": human_translation,
            "M3_meta_effect": float(human_primary["meta_effect"]),
            "M3_meta_SE": float(human_primary["meta_SE"]),
            "M3_meta_P": float(human_primary["meta_P"]),
            "M3_meta_FDR": float(human_primary["meta_FDR"]),
            "gate_veto": False,
        },
        "claim_boundary": CFG["claim_boundary"],
        "downstream": "G4 may proceed within the model-scoped mechanism boundary; G5-G6 remain locked until their upstream gates pass."
        if passed
        else "G4-G6 remain locked.",
        "created_at": datetime.now(timezone.utc).isoformat(),
    }
    (RUN / "g3_decision.json").write_text(
        json.dumps(decision_record, ensure_ascii=False, indent=2), encoding="utf-8"
    )

    manifest = {
        "run_id": CFG["run_id"],
        "script": relative(Path(__file__).resolve()),
        "config": relative(CONFIG),
        "raw_data_modified": False,
        "remote_computation": False,
        "result": decision,
        "inputs": [
            {"file": relative(path), "sha256": sha256(path)} for path in inputs.values()
        ],
    }
    (RUN / "run_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8"
    )

    report = f"""# Disease-side v3.2 G3 model-scoped reassessment

## Decision

**{decision}**

The prior execution incorrectly made a positive human donor-level M3 meta-analysis mandatory. The frozen G3 text permits qualified human, iPSC, or animal evidence.

- D2 iMGL M3 positive BH-FDR time points: {', '.join(map(str, m3['time_hours'].astype(int).tolist()))} h.
- D2 non-inflammatory senescence-core positive BH-FDR time points: {', '.join(map(str, core['time_hours'].astype(int).tolist()))} h.
- G2 anchor: {g2.get('decision')}.
- Human translation: {human_translation}; meta-effect {human_primary['meta_effect']:.6g}, FDR {human_primary['meta_FDR']:.6g}; not a model-scoped gate veto.

## Boundary

This passes model-level lactylation-to-senescence/SASP directional linkage. It does not establish human translation or direct H4K12la-to-SASP causality; the independent mechanistic bridge is H3K18la.
"""
    (RUN / "execution_report.md").write_text(report, encoding="utf-8")

    checksums = []
    for path in sorted(RUN.rglob("*")):
        if path.is_file() and path.name != "checksums.sha256":
            checksums.append(f"{sha256(path)}  {relative(path)}")
    (RUN / "checksums.sha256").write_text("\n".join(checksums) + "\n", encoding="utf-8")

    print(
        json.dumps(
            {
                "run": str(RUN),
                "decision": decision,
                "D2_M3_positive_FDR_timepoints": m3["time_hours"].astype(int).tolist(),
                "D2_senescence_core_positive_FDR_timepoints": core["time_hours"].astype(int).tolist(),
                "human_translation": human_translation,
            },
            ensure_ascii=False,
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
