from __future__ import annotations

import hashlib
import json
import math
import shutil
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from scipy.stats import t as student_t
from scipy.stats import ttest_ind


ROOT = Path(__file__).resolve().parents[1]
CONFIG = ROOT / "configs" / "g3_v3.2_execution.yaml"
CFG = yaml.safe_load(CONFIG.read_text(encoding="utf-8"))
RUN = ROOT / "results" / CFG["run_id"]
if RUN.exists():
    raise RuntimeError(f"Refusing to overwrite existing run directory: {RUN}")
TABLES = RUN / "tables"
META = RUN / "meta"
GATES = RUN / "gates"
REGISTRIES = RUN / "registries"
for d in [TABLES, META, GATES, REGISTRIES]:
    d.mkdir(parents=True, exist_ok=False)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def rel(path: Path) -> str:
    return path.relative_to(ROOT).as_posix()


def write_tsv(path: Path, rows: list[dict]) -> None:
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)


def bh(values: list[float]) -> np.ndarray:
    p = np.asarray(values, dtype=float)
    out = np.full(len(p), np.nan)
    ok = np.isfinite(p)
    if not ok.any():
        return out
    order = np.argsort(p[ok])
    ranked = p[ok][order]
    adj = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    tmp = np.empty(len(ranked))
    tmp[order] = np.minimum(adj, 1.0)
    out[ok] = tmp
    return out


def read_count_matrix(path: Path, metadata_path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    counts = pd.read_csv(path, sep="\t", compression="infer")
    gene_col = counts.columns[0]
    counts[gene_col] = counts[gene_col].astype(str).str.upper()
    counts = counts.rename(columns={gene_col: "gene"}).set_index("gene")
    counts = counts.groupby(level=0, sort=False).sum(numeric_only=True)
    meta = pd.read_csv(metadata_path, sep="\t")
    if "group_id" not in meta.columns:
        raise ValueError(f"group_id missing from {metadata_path}")
    missing = [x for x in meta["group_id"] if x not in counts.columns]
    if missing:
        raise ValueError(f"metadata/count mismatch in {metadata_path}: {missing[:3]}")
    counts = counts.loc[:, meta["group_id"].tolist()]
    return counts.astype(float), meta


def normalize_human(counts: pd.DataFrame) -> pd.DataFrame:
    lib = counts.sum(axis=0).replace(0, np.nan)
    return np.log2(counts.div(lib, axis=1) * 1_000_000 + 2)


def zscore_rows(expr: pd.DataFrame) -> pd.DataFrame:
    means = expr.mean(axis=1)
    sds = expr.std(axis=1, ddof=1).replace(0, np.nan)
    return expr.sub(means, axis=0).div(sds, axis=0).fillna(0.0)


def module_score(expr: pd.DataFrame, genes: list[str]) -> tuple[pd.Series, int]:
    present = [g.upper() for g in genes if g.upper() in expr.index]
    if len(present) < 2:
        return pd.Series(np.nan, index=expr.columns), len(present)
    return expr.loc[present].mean(axis=0), len(present)


def fit_ols(meta: pd.DataFrame, y: pd.Series, predictor: str = "disease", extra_numeric: list[str] | None = None) -> dict:
    extra_numeric = extra_numeric or []
    m = meta.copy()
    m["y"] = pd.to_numeric(y, errors="coerce").values
    keep = m["y"].notna()
    for col in ["age", "RIN"] + extra_numeric:
        if col in m.columns:
            m[col] = pd.to_numeric(m[col], errors="coerce")
            keep &= m[col].notna()
    if predictor == "disease":
        keep &= m["disease"].isin(["Control", "AD"])
    elif predictor not in m.columns:
        return {"estimable": False, "reason": f"missing predictor {predictor}"}
    m = m.loc[keep].copy()
    if len(m) < 4:
        return {"estimable": False, "reason": "fewer than four observations"}
    cols = ["intercept"]
    xparts = [np.ones(len(m))]
    if predictor == "disease":
        xparts.append((m["disease"].astype(str) == "AD").astype(float).to_numpy())
        cols.append("diseaseAD")
    else:
        xparts.append(pd.to_numeric(m[predictor], errors="coerce").to_numpy())
        cols.append(predictor)
        # The G2-to-G3 association is adjusted for disease status so that
        # cross-sectional AD-control separation is not mistaken for a
        # within-cohort module association.
        if "disease" in m.columns and m["disease"].isin(["Control", "AD"]).all() and m["disease"].nunique() == 2:
            xparts.append((m["disease"].astype(str) == "AD").astype(float).to_numpy())
            cols.append("diseaseAD")
    for col in ["age", "RIN"] + extra_numeric:
        if col in m.columns and m[col].nunique(dropna=True) > 1:
            xparts.append(m[col].to_numpy(dtype=float))
            cols.append(col)
    for col in ["sex", "batch"]:
        if col in m.columns and m[col].notna().all() and m[col].nunique() > 1:
            dummies = pd.get_dummies(m[col].astype(str), prefix=col, drop_first=True, dtype=float)
            for dcol in dummies.columns:
                xparts.append(dummies[dcol].to_numpy())
                cols.append(dcol)
    X = np.column_stack(xparts)
    yv = m["y"].to_numpy(dtype=float)
    beta, _, rank, _ = np.linalg.lstsq(X, yv, rcond=None)
    fitted = X @ beta
    df_resid = len(yv) - int(rank)
    if df_resid <= 0:
        return {"estimable": False, "reason": "non-positive residual degrees of freedom"}
    rss = float(np.sum((yv - fitted) ** 2))
    sigma2 = rss / df_resid
    cov = sigma2 * np.linalg.pinv(X.T @ X)
    idx = 1
    se = math.sqrt(max(float(cov[idx, idx]), 0.0))
    effect = float(beta[idx])
    if se <= 0 or not np.isfinite(se):
        p = 0.0 if effect != 0 else 1.0
        tstat = np.inf if effect != 0 else 0.0
    else:
        tstat = effect / se
        p = float(2 * student_t.sf(abs(tstat), df_resid))
    return {
        "estimable": True,
        "n": int(len(yv)),
        "rank": int(rank),
        "df_resid": int(df_resid),
        "effect": effect,
        "SE": se,
        "t": float(tstat),
        "P.Value": p,
        "CI95_low": effect - float(student_t.ppf(0.975, df_resid)) * se,
        "CI95_high": effect + float(student_t.ppf(0.975, df_resid)) * se,
        "model_columns": ";".join(cols),
    }


def random_effects(rows: pd.DataFrame, effect_col: str = "effect", se_col: str = "SE") -> dict:
    x = rows[[effect_col, se_col]].copy()
    x = x[np.isfinite(x[effect_col]) & np.isfinite(x[se_col]) & (x[se_col] > 0)]
    k = len(x)
    if k == 0:
        return {"k": 0, "meta_effect": np.nan, "meta_SE": np.nan, "meta_P": np.nan, "tau2": np.nan, "I2": np.nan}
    y = x[effect_col].to_numpy(float)
    se = x[se_col].to_numpy(float)
    w = 1 / se**2
    fixed = np.sum(w * y) / np.sum(w)
    if k == 1:
        meta, meta_se, tau2, q_p, i2 = y[0], se[0], 0.0, np.nan, 0.0
    else:
        Q = np.sum(w * (y - fixed) ** 2)
        c = np.sum(w) - np.sum(w**2) / np.sum(w)
        tau2 = max(0.0, (Q - (k - 1)) / max(c, np.finfo(float).eps))
        wr = 1 / (se**2 + tau2)
        meta = np.sum(wr * y) / np.sum(wr)
        meta_se = math.sqrt(1 / np.sum(wr))
        q_p = float(__import__("scipy").stats.chi2.sf(Q, k - 1))
        i2 = max(0.0, (Q - (k - 1)) / Q) if Q > 0 else 0.0
    z = meta / meta_se if meta_se > 0 else 0.0
    p = float(2 * student_t.sf(abs(z), max(k - 1, 1))) if np.isfinite(z) else 0.0
    return {"k": int(k), "meta_effect": float(meta), "meta_SE": float(meta_se), "meta_P": p, "tau2": float(tau2), "I2": float(i2), "Q_p": q_p}


M3 = ["CDKN1A", "CDKN2A", "IL6", "IL1B", "TNF", "CCL2", "CCL3", "CCL4", "CSF1", "CXCL10", "MMP9", "MMP14", "SERPINE1", "IGFBP3", "IGFBP7", "NFKB1", "RELA", "TGFB1"]
INFLAMMATION = ["IL6", "IL1B", "TNF", "CCL2", "CCL3", "CCL4", "CSF1", "CXCL10", "NFKB1", "RELA"]
SEN_CORE = ["CDKN1A", "CDKN2A", "MMP9", "MMP14", "SERPINE1", "IGFBP3", "IGFBP7", "TGFB1"]


def load_m2_weights(path: Path) -> dict[str, float]:
    tab = pd.read_csv(path, sep="\t")
    tab = tab[(tab["context"] == "M") & tab["delta_presence"].notna() & (tab["delta_presence"] != 0)]
    weights: dict[str, float] = {}
    for _, row in tab.iterrows():
        symbols = str(row["human_one_to_one_symbols"]).split(";")
        for symbol in symbols:
            symbol = symbol.strip().upper()
            if symbol and symbol != "NAN":
                weights[symbol] = weights.get(symbol, 0.0) + float(row["delta_presence"])
    return weights


def score_human_dataset(cohort: str, counts_path: Path, meta_path: Path, m2_weights: dict[str, float]) -> tuple[pd.DataFrame, dict]:
    counts, meta = read_count_matrix(counts_path, meta_path)
    expr = normalize_human(counts)
    scores = pd.DataFrame(index=meta["group_id"].astype(str))
    scores["cohort"] = cohort
    for col in ["donor_id", "disease", "region", "age", "sex", "RIN", "rin", "batch", "n_cells_used", "n_cells"]:
        if col in meta.columns:
            scores[col] = meta.set_index("group_id")[col]
    if "RIN" not in scores.columns and "rin" in scores.columns:
        scores["RIN"] = scores["rin"]
    if "n_cells_used" not in scores.columns and "n_cells" in scores.columns:
        scores["n_cells_used"] = scores["n_cells"]
    scores["M3"] , n_m3 = module_score(expr, M3)
    scores["inflammation"] , n_inf = module_score(expr, INFLAMMATION)
    scores["senescence_core"] , n_core = module_score(expr, SEN_CORE)
    m2_genes = [g for g in m2_weights if g in expr.index]
    if len(m2_genes) >= 2:
        z = zscore_rows(expr.loc[m2_genes])
        w = pd.Series({g: m2_weights[g] for g in m2_genes})
        scores["M2_H4K12la_proxy"] = (z.mul(w, axis=0).sum(axis=0) / np.abs(w).sum()).reindex(scores.index).to_numpy()
    else:
        scores["M2_H4K12la_proxy"] = np.nan
    scores.reset_index(names="group_id").to_csv(TABLES / f"{cohort}_G3_donor_region_scores.tsv", sep="\t", index=False)
    return scores.reset_index(names="group_id"), {"n_genes": int(expr.shape[0]), "n_M3_present": n_m3, "n_inflammation_present": n_inf, "n_senescence_core_present": n_core, "n_M2_present": len(m2_genes), "n_groups": int(len(meta))}


def temporal_analysis() -> tuple[pd.DataFrame, dict]:
    tpm_path = ROOT / CFG["inputs"]["d2_gene_tpm"]
    meta_path = ROOT / CFG["inputs"]["d2_metadata"]
    tpm = pd.read_csv(tpm_path, sep="\t", index_col=0)
    tpm.index = tpm.index.astype(str).str.upper()
    tpm = tpm.groupby(level=0, sort=False).sum()
    meta = pd.read_csv(meta_path, sep="\t")
    cols = meta["sample_id_in_expression_table"].tolist()
    expr = np.log2(tpm.loc[:, cols].astype(float) + 1)
    z = zscore_rows(expr)
    rows = []
    for module, genes in [("M3_SASP_Senescence", M3), ("inflammation_control", INFLAMMATION), ("senescence_core", SEN_CORE)]:
        present = [g for g in genes if g in z.index]
        score = z.loc[present].mean(axis=0)
        for time in [4, 24, 96]:
            x = score.loc[meta.loc[meta["time_hours"] == time, "sample_id_in_expression_table"]].to_numpy(float)
            y = score.loc[meta.loc[meta["time_hours"] == 0, "sample_id_in_expression_table"]].to_numpy(float)
            tt = ttest_ind(x, y, equal_var=True)
            rows.append({"module": module, "time_hours": time, "n_post": len(x), "n_h0": len(y), "logFC_like_mean_difference": float(np.mean(x) - np.mean(y)), "t": float(tt.statistic), "P.Value": float(tt.pvalue), "n_present": len(present), "direction": "A_beta_up" if np.mean(x) > np.mean(y) else "A_beta_down"})
    out = pd.DataFrame(rows)
    out["FDR"] = bh(out["P.Value"].tolist())
    out.to_csv(TABLES / "D2_G3_temporal_module_effects.tsv", sep="\t", index=False)
    return out, {"input": rel(tpm_path), "n_samples": int(len(meta)), "modules": sorted(out["module"].unique())}


def main() -> None:
    m2_weights = load_m2_weights(ROOT / CFG["inputs"]["m2_weights"])
    temporal, temporal_info = temporal_analysis()
    sources = {
        "D3_GSE160936": (ROOT / CFG["inputs"]["d3_counts"], ROOT / CFG["inputs"]["d3_metadata"]),
        "D7_GSE148822": (ROOT / CFG["inputs"]["d7_counts"], ROOT / CFG["inputs"]["d7_metadata"]),
        "D8_GSE174367": (ROOT / CFG["inputs"]["d8_counts"], ROOT / CFG["inputs"]["d8_metadata"]),
    }
    all_scores = []
    input_audit = []
    for cohort, (counts_path, meta_path) in sources.items():
        scores, info = score_human_dataset(cohort, counts_path, meta_path, m2_weights)
        all_scores.append(scores)
        for path in [counts_path, meta_path]:
            input_audit.append({"file": rel(path), "sha256": sha256(path), "role": cohort})
        input_audit.append({"file": rel(counts_path), "sha256": sha256(counts_path), "role": f"{cohort}_counts_reused_processed"})
    scores = pd.concat(all_scores, ignore_index=True)

    effect_rows = []
    association_rows = []
    for cohort, group in scores.groupby("cohort", sort=True):
        for region, r in group.groupby("region", sort=True):
            r = r.copy()
            for name, extra, predictor in [
                ("M3_disease_primary", [], "disease"),
                ("M3_disease_cell_count_adjusted", ["log_cells"], "disease"),
                ("M3_disease_inflammation_adjusted", ["inflammation"], "disease"),
            ]:
                if "log_cells" in extra:
                    r["log_cells"] = np.log1p(pd.to_numeric(r["n_cells_used"], errors="coerce"))
                fit = fit_ols(r, r["M3"], predictor="disease", extra_numeric=extra)
                if fit.get("estimable"):
                    effect_rows.append({"cohort": cohort, "region": region, "outcome": name, **fit, "n_AD": int((r["disease"] == "AD").sum()), "n_Control": int((r["disease"] == "Control").sum()), "n_donor_region": int(len(r))})
            fit_core = fit_ols(r, r["senescence_core"], predictor="disease")
            if fit_core.get("estimable"):
                effect_rows.append({"cohort": cohort, "region": region, "outcome": "senescence_core_disease", **fit_core, "n_AD": int((r["disease"] == "AD").sum()), "n_Control": int((r["disease"] == "Control").sum()), "n_donor_region": int(len(r))})
            r["M2_H4K12la_proxy"] = pd.to_numeric(r["M2_H4K12la_proxy"], errors="coerce")
            fit_assoc = fit_ols(r, r["M3"], predictor="M2_H4K12la_proxy", extra_numeric=[])
            if fit_assoc.get("estimable"):
                association_rows.append({"cohort": cohort, "region": region, "outcome": "M2_to_M3_association", **fit_assoc, "n_AD": int((r["disease"] == "AD").sum()), "n_Control": int((r["disease"] == "Control").sum()), "n_donor_region": int(len(r))})
    effects = pd.DataFrame(effect_rows).drop_duplicates(subset=["cohort", "region", "outcome"])
    associations = pd.DataFrame(association_rows)
    effects.to_csv(TABLES / "human_G3_region_effects.tsv", sep="\t", index=False)
    associations.to_csv(TABLES / "human_G2_G3_region_associations.tsv", sep="\t", index=False)

    meta_rows = []
    for outcome, subset in pd.concat([effects, associations], ignore_index=True).groupby("outcome"):
        re = random_effects(subset, "effect", "SE")
        meta_rows.append({"outcome": outcome, **re, "direction": "AD_up_or_positive" if re["meta_effect"] > 0 else "AD_down_or_negative"})
    meta_summary = pd.DataFrame(meta_rows)
    meta_summary["meta_FDR"] = bh(meta_summary["meta_P"].tolist())
    meta_summary.to_csv(META / "G3_random_effects_meta.tsv", sep="\t", index=False)

    loo_rows = []
    for outcome in meta_summary["outcome"]:
        for excluded in sources:
            subset = pd.concat([effects, associations], ignore_index=True)
            subset = subset[(subset["outcome"] == outcome) & (subset["cohort"] != excluded)]
            re = random_effects(subset, "effect", "SE")
            full = meta_summary.loc[meta_summary["outcome"] == outcome].iloc[0]
            loo_rows.append({"outcome": outcome, "excluded_cohort": excluded, **re, "direction_stable": bool(np.sign(re["meta_effect"]) == np.sign(full["meta_effect"]) and re["meta_effect"] != 0 and full["meta_effect"] != 0)})
    pd.DataFrame(loo_rows).to_csv(META / "G3_leave_one_cohort_out.tsv", sep="\t", index=False)

    m3 = meta_summary.loc[meta_summary["outcome"] == "M3_disease_primary"].iloc[0]
    m3_cell = meta_summary.loc[meta_summary["outcome"] == "M3_disease_cell_count_adjusted"].iloc[0]
    m3_inf = meta_summary.loc[meta_summary["outcome"] == "M3_disease_inflammation_adjusted"].iloc[0]
    assoc = meta_summary.loc[meta_summary["outcome"] == "M2_to_M3_association"].iloc[0]
    temporal_m3 = temporal[temporal["module"] == "M3_SASP_Senescence"]
    temporal_m3_positive = int(((temporal_m3["logFC_like_mean_difference"] > 0) & (temporal_m3["FDR"] < CFG["analysis"]["fdr"])).sum())
    temporal_compat = temporal_m3_positive >= 1
    human_association = bool(m3["meta_effect"] > 0 and m3["meta_FDR"] < CFG["analysis"]["fdr"])
    g2_association = bool(assoc["meta_effect"] > 0 and assoc["meta_FDR"] < CFG["analysis"]["fdr"])
    cell_audit = bool(np.sign(m3_cell["meta_effect"]) == np.sign(m3["meta_effect"]))
    inflammation_audit = bool(np.sign(m3_inf["meta_effect"]) == np.sign(m3["meta_effect"]))
    g3_pass = human_association and (temporal_compat or g2_association) and cell_audit and inflammation_audit
    adequately_contradictory = bool(m3["meta_effect"] < 0 and m3["meta_FDR"] < CFG["analysis"]["fdr"])
    decision = "G3_PASS" if g3_pass else ("G3_FAIL" if adequately_contradictory else "G3_NOT_YET_PASS")
    gate_rows = [
        {"criterion": "positive_donor_aware_M3_AD_association", "observed": human_association, "detail": f"meta_effect={m3['meta_effect']:.6g}; meta_FDR={m3['meta_FDR']:.6g}"},
        {"criterion": "directional_G2_compatibility", "observed": temporal_compat or g2_association, "detail": f"D2_M3_positive_FDR_hits={temporal_m3_positive}; human_M2_to_M3_meta_effect={assoc['meta_effect']:.6g}; meta_FDR={assoc['meta_FDR']:.6g}"},
        {"criterion": "cell_composition_audit", "observed": cell_audit, "detail": f"primary_effect={m3['meta_effect']:.6g}; cell_count_adjusted_effect={m3_cell['meta_effect']:.6g}"},
        {"criterion": "generic_inflammation_audit", "observed": inflammation_audit, "detail": f"primary_effect={m3['meta_effect']:.6g}; inflammation_adjusted_effect={m3_inf['meta_effect']:.6g}"},
    ]
    write_tsv(GATES / "G3_gate_v3.2.tsv", gate_rows)
    decision_json = {
        "run_id": CFG["run_id"],
        "protocol": "AD disease-side v3.2",
        "decision": decision,
        "human_primary_cohorts": ["D3_GSE160936", "D7_GSE148822", "D8_GSE174367"],
        "d2_temporal": {"positive_M3_post_A_beta_FDR_hits": temporal_m3_positive, "compatible": temporal_compat},
        "meta_summary": meta_summary.to_dict(orient="records"),
        "gate_rows": gate_rows,
        "downstream": "G4-G6 remain locked because G3 did not pass." if decision != "G3_PASS" else "G4 may proceed; G5-G6 remain fail-closed until their upstream gates pass.",
        "claim_boundary": "transcriptional SASP/senescence proxy; not direct senescence, secreted SASP, lactate, H4K12la, iron, or lipid-peroxidation measurement",
        "limitations": ["M3 is an operational 18-gene proxy", "human effects are cohort-region units and D3 regions share donors", "D2 is iPSC A_beta time-course support, not human AD-Control evidence", "cell-level data were aggregated to donor-region pseudobulk", "D6 excluded from primary gate"],
        "created_at": datetime.now(timezone.utc).isoformat(),
    }
    (RUN / "g3_decision.json").write_text(json.dumps(decision_json, ensure_ascii=False, indent=2), encoding="utf-8")
    manifest = {"run_id": CFG["run_id"], "protocol": CFG["protocol"], "config": rel(CONFIG), "script": "src/30_run_g3_v3_2.py", "raw_data_read_only": True, "raw_reads_downloaded": False, "processed_inputs_reused": True, "result": decision, "input_files": input_audit}
    (RUN / "run_manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")
    report = f"""# Disease-side v3.2 G3 execution report

Run ID: `{CFG['run_id']}`  
Protocol: `{CFG['protocol']}`

## Decision

**{decision}**

G3 tests donor-aware Senescence/SASP disease association, directional compatibility with G2, cell-composition sensitivity, and generic-inflammation sensitivity. The independent unit is donor within brain region; D2 replicates are used only for the Aβ time-course.

- D2 M3 positive post-Aβ timepoints passing BH-FDR < 0.05: `{temporal_m3_positive}`.
- Human primary M3 random-effects meta: effect `{m3['meta_effect']:.6g}`, P `{m3['meta_P']:.6g}`, FDR `{m3['meta_FDR']:.6g}`, I2 `{m3['I2']:.6g}`.
- Human M2-to-M3 association meta: effect `{assoc['meta_effect']:.6g}`, P `{assoc['meta_P']:.6g}`, FDR `{assoc['meta_FDR']:.6g}`.
- Cell-count-adjusted M3 meta: effect `{m3_cell['meta_effect']:.6g}`.
- Inflammation-adjusted M3 meta: effect `{m3_inf['meta_effect']:.6g}`.

## Claim boundary

M3 is an operational transcriptional Senescence/SASP proxy. It cannot be reported as direct cellular senescence, secreted SASP, lactate, H4K12la, iron, or lipid-peroxidation measurement. D6 was excluded from the primary gate.

## Downstream policy

{decision_json['downstream']}
"""
    (RUN / "execution_report.md").write_text(report, encoding="utf-8")
    checksum_rows = []
    for path in sorted(RUN.rglob("*")):
        if path.is_file() and path.name != "checksums.sha256":
            checksum_rows.append(f"{sha256(path)}  {rel(path)}")
    (RUN / "checksums.sha256").write_text("\n".join(checksum_rows) + "\n", encoding="utf-8")
    print(json.dumps({"run": str(RUN), "decision": decision, "m3_meta_effect": float(m3["meta_effect"]), "m3_meta_FDR": float(m3["meta_FDR"]), "m2_m3_meta_effect": float(assoc["meta_effect"]), "m2_m3_meta_FDR": float(assoc["meta_FDR"]), "d2_m3_positive_FDR_hits": temporal_m3_positive}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
