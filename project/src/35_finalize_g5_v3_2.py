from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from scipy.stats import chi2, norm


ROOT = Path(__file__).resolve().parents[1]
CFG_PATH = ROOT / "configs" / "g5_v3.2_execution.yaml"
CFG = yaml.safe_load(CFG_PATH.read_text(encoding="utf-8"))
RUN = ROOT / "results" / CFG["run_id"]
if not RUN.exists():
    raise FileNotFoundError("Run directory missing; execute src/34_run_g5_region_effects.R first")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def rel(path: Path) -> str:
    return path.relative_to(ROOT).as_posix()


def random_one(group: pd.DataFrame) -> dict:
    x = group[np.isfinite(group["logFC"]) & np.isfinite(group["SE"]) & (group["SE"] > 0)].copy()
    k = len(x)
    if k == 0:
        return {"k": 0, "beta": np.nan, "SE": np.nan, "z": np.nan, "P.Value": np.nan, "tau2": np.nan, "I2": np.nan, "Q_p": np.nan}
    y = x["logFC"].to_numpy(float)
    se = x["SE"].to_numpy(float)
    w = 1.0 / (se**2)
    fixed = np.sum(w * y) / np.sum(w)
    if k == 1:
        z = y[0] / se[0]
        return {"k": 1, "beta": y[0], "SE": se[0], "z": z, "P.Value": 2 * norm.sf(abs(z)), "tau2": 0.0, "I2": 0.0, "Q_p": np.nan}
    q = float(np.sum(w * (y - fixed) ** 2))
    c = float(np.sum(w) - np.sum(w**2) / np.sum(w))
    tau2 = max(0.0, (q - (k - 1)) / c) if c > 0 else 0.0
    wr = 1.0 / (se**2 + tau2)
    beta = float(np.sum(wr * y) / np.sum(wr))
    se_beta = float(np.sqrt(1.0 / np.sum(wr)))
    z = beta / se_beta
    p = float(2 * norm.sf(abs(z)))
    i2 = max(0.0, (q - (k - 1)) / q) if q > 0 else 0.0
    return {"k": k, "beta": beta, "SE": se_beta, "z": z, "P.Value": p, "tau2": tau2, "I2": i2, "Q_p": float(chi2.sf(q, k - 1)) if q > 0 else 1.0}


def meta_table(effects: pd.DataFrame, by: list[str]) -> pd.DataFrame:
    rows: list[dict] = []
    for keys, group in effects.groupby(by, sort=True):
        if not isinstance(keys, tuple):
            keys = (keys,)
        result = {key: value for key, value in zip(by, keys)}
        result.update(random_one(group))
        rows.append(result)
    out = pd.DataFrame(rows)
    out["FDR"] = np.nan
    valid = np.isfinite(out["P.Value"])
    if valid.any():
        p = out.loc[valid, "P.Value"].to_numpy(float)
        order = np.argsort(p)
        adj = p[order] * len(p) / np.arange(1, len(p) + 1)
        adj = np.minimum.accumulate(adj[::-1])[::-1]
        restored = np.empty(len(p))
        restored[order] = np.minimum(adj, 1.0)
        out.loc[valid, "FDR"] = restored
    return out


def main() -> None:
    g4_path = ROOT / CFG["inputs"]["g4_decision"]
    g4 = json.loads(g4_path.read_text(encoding="utf-8"))
    effects_path = RUN / "tables" / "astrocyte_region_module_effects.tsv"
    effects = pd.read_csv(effects_path, sep="\t")
    effects["P.Value"] = pd.to_numeric(effects["P.Value"], errors="coerce")
    effects["logFC"] = pd.to_numeric(effects["logFC"], errors="coerce")
    effects["SE"] = pd.to_numeric(effects["SE"], errors="coerce")

    cohort_meta = meta_table(effects, ["cohort", "module"])
    overall = meta_table(effects, ["module"])
    loo_rows = []
    cohorts = sorted(effects["cohort"].dropna().unique().tolist())
    for excluded in cohorts:
        reduced = effects.loc[effects["cohort"] != excluded].copy()
        x = meta_table(reduced, ["module"])
        x["excluded_cohort"] = excluded
        loo_rows.append(x)
    loo = pd.concat(loo_rows, ignore_index=True) if loo_rows else pd.DataFrame()

    cohort_meta.to_csv(RUN / "tables" / "astrocyte_cohort_random_effects.tsv", sep="\t", index=False)
    overall.to_csv(RUN / "tables" / "astrocyte_primary_random_effects.tsv", sep="\t", index=False)
    loo.to_csv(RUN / "tables" / "astrocyte_leave_one_cohort_out.tsv", sep="\t", index=False)

    primary = ["D3_GSE160936", "D7_GSE148822", "D8_GSE174367"]
    estimability_rows = []
    for cohort in primary:
        x = effects.loc[effects["cohort"] == cohort]
        estimability_rows.append({
            "cohort": cohort,
            "estimable": bool(len(x) > 0 and x["region"].nunique() >= 1 and (x["n_AD_donors"] >= 2).all() and (x["n_Control_donors"] >= 2).all()),
            "region_strata": int(x["region"].nunique()),
            "regions": ";".join(sorted(x["region"].dropna().unique().tolist())),
        })
    estimability = pd.DataFrame(estimability_rows)
    estimability.to_csv(RUN / "tables" / "astrocyte_estimability.tsv", sep="\t", index=False)

    module_rows = []
    for module in overall["module"].tolist():
        om = overall.loc[overall["module"] == module].iloc[0]
        cm = cohort_meta.loc[cohort_meta["module"] == module]
        concordant = int(np.sum(np.sign(cm["beta"].to_numpy(float)) == np.sign(float(om["beta"]))))
        lm = loo.loc[loo["module"] == module]
        stable = bool(len(lm) == len(cohorts) and np.isfinite(lm["beta"]).all() and (np.sign(lm["beta"]) == np.sign(float(om["beta"]))).all())
        module_rows.append({
            "module": module,
            "overall_beta": float(om["beta"]),
            "overall_SE": float(om["SE"]),
            "overall_P": float(om["P.Value"]),
            "overall_FDR": float(om["FDR"]),
            "concordant_cohorts": concordant,
            "minimum_concordant": int(CFG["analysis"]["gate"]["minimum_concordant_cohorts"]),
            "leave_one_cohort_out_direction_stable": stable,
        })
    module_gate = pd.DataFrame(module_rows)
    module_gate["pass"] = (
        (module_gate["overall_FDR"] < float(CFG["analysis"]["gate"]["random_effects_meta_fdr_max"]))
        & (module_gate["concordant_cohorts"] >= module_gate["minimum_concordant"])
        & module_gate["leave_one_cohort_out_direction_stable"]
    )
    module_gate.to_csv(RUN / "gates" / "G5_module_gate.tsv", sep="\t", index=False)

    upstream_pass = g4.get("decision") in {"G4_PASS_MODEL_SCOPED", "G4_PASS"}
    all_estimable = bool(estimability["estimable"].all())
    region_stratified = bool((estimability["region_strata"] >= 1).all())
    module_pass = bool(module_gate["pass"].any())
    decision = "G5_PASS_MODEL_SCOPED" if upstream_pass and all_estimable and region_stratified and module_pass else "G5_NOT_YET_PASS"
    stable_pass = bool((module_gate["pass"] & module_gate["leave_one_cohort_out_direction_stable"]).any())

    gate = pd.DataFrame([
        {"criterion": "upstream_G4_pass", "required": True, "observed": upstream_pass, "detail": g4.get("decision", "missing")},
        {"criterion": "all_primary_cohorts_estimable", "required": True, "observed": all_estimable, "detail": ";".join(f"{r.cohort}={r.estimable}" for r in estimability.itertuples())},
        {"criterion": "donor_aware_region_stratified_analysis", "required": True, "observed": region_stratified, "detail": ";".join(f"{r.cohort}={r.regions}" for r in estimability.itertuples())},
        {"criterion": "random_effects_meta_FDR_and_cohort_concordance", "required": True, "observed": module_pass, "detail": ";".join(f"{module}={passed}" for module, passed in zip(module_gate["module"], module_gate["pass"]))},
        {"criterion": "leave_one_cohort_out_direction_stability", "required": True, "observed": stable_pass, "detail": "checked for every response module"},
        {"criterion": "spatial_colocalization_primary_gate", "required": False, "observed": False, "detail": "not used; G5 requires expression direction before G6"},
    ])
    gate.to_csv(RUN / "gates" / "G5_gate_v3.2.tsv", sep="\t", index=False)

    decision_record = {
        "run_id": CFG["run_id"],
        "protocol": "AD disease-side v3.2",
        "decision": decision,
        "upstream_G4": g4.get("decision"),
        "estimability": estimability.to_dict(orient="records"),
        "passing_modules": module_gate.loc[module_gate["pass"], "module"].tolist(),
        "claim_boundary": {
            "allowed": "donor-aware region-stratified astrocyte transcriptional response",
            "not_allowed": "spatial causality, ligand-receptor causality, or microglia-to-astrocyte direction without G6",
            "spatial_support": "not used as a primary G5 gate",
        },
        "downstream": "G6 may proceed to its separate spatial-support gate" if decision == "G5_PASS_MODEL_SCOPED" else "G6 remains locked",
        "created_at": datetime.now(timezone.utc).isoformat(),
    }
    (RUN / "g5_decision.json").write_text(json.dumps(decision_record, ensure_ascii=False, indent=2), encoding="utf-8")

    report = f"""# Disease-side v3.2 G5 execution

## Decision

**{decision}**

G5 used donor-level astrocyte pseudobulk with region-stratified AD-Control models in D3/GSE160936, D7/GSE148822, and D8/GSE174367. D7 astrocytes were selected by the existing marker/QC rule; D8 used the source ASC annotation with the frozen QC thresholds. Random-effects meta-analysis and leave-one-cohort-out direction checks were performed for all predeclared response modules. Spatial colocalization was not used as a primary directional gate.

- Estimability: {'; '.join(f'{r.cohort}={r.estimable}' for r in estimability.itertuples())}
- Regions: {'; '.join(f'{r.cohort}={r.regions}' for r in estimability.itertuples())}
- Passing response modules: {', '.join(module_gate.loc[module_gate['pass'], 'module'].tolist()) if module_pass else 'none'}

## Boundary

This gate tests an astrocyte transcriptional response. It does not establish spatial directionality, ligand-receptor causality, or a microglia-to-astrocyte causal effect. G6 is {'unlocked for its separate spatial-support gate' if decision == 'G5_PASS_MODEL_SCOPED' else 'locked because G5 did not pass'}.
"""
    (RUN / "execution_report.md").write_text(report, encoding="utf-8")

    input_values = list(CFG["inputs"].values()) + ["src/34_run_g5_region_effects.R", "src/35_finalize_g5_v3_2.py", "configs/g5_v3.2_execution.yaml"]
    manifest_inputs = []
    for value in input_values:
        path = ROOT / value
        if path.exists():
            manifest_inputs.append({"file": rel(path), "sha256": sha256(path), "size": path.stat().st_size})
    manifest = {"run_id": CFG["run_id"], "script": "src/35_finalize_g5_v3_2.py", "region_effect_script": "src/34_run_g5_region_effects.R", "raw_data_modified": False, "remote_computation": False, "result": decision, "inputs": manifest_inputs}
    (RUN / "run_manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")

    checksums = []
    for path in sorted(RUN.rglob("*")):
        if path.is_file() and path.name != "checksums.sha256":
            checksums.append(f"{sha256(path)}  {rel(path)}")
    (RUN / "checksums.sha256").write_text("\n".join(checksums) + "\n", encoding="utf-8")
    print(json.dumps({"run": rel(RUN), "decision": decision, "passing_modules": module_gate.loc[module_gate["pass"], "module"].tolist(), "estimability": estimability.to_dict(orient="records")}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
