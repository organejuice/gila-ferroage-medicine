#!/usr/bin/env python3
"""Conservative D4 CosMx spatial support analysis.

This is a descriptive support layer only: GSE263791 contains two AD CosMx
sections, so no AD-control disease contrast or causal communication model is
estimated here. Cell identities are provisional marker-score labels because
the downloaded processed files do not contain an authoritative cell-type label.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.stats import spearmanr

SAMPLES = [
    {
        "sample_id": "GSM8199188",
        "prefix": "GSM8199188_ID61-ID62_S10",
        "subject_label": "AD2_AD3_CosMx",
        "disease": "AD",
    },
    {
        "sample_id": "GSM8199189",
        "prefix": "GSM8199189_ID67-ID68_S18",
        "subject_label": "AD4_AD5_CosMx",
        "disease": "AD",
    },
]

MICROGLIA_MARKERS = ["Hexb", "C1qa", "C1qb", "C1qc", "P2ry12", "Tmem119", "Cx3cr1", "Aif1", "Tyrobp", "Lpl"]
ASTROCYTE_MARKERS = ["Aldoc", "Aqp4", "Gfap", "Slc1a2", "Slc1a3", "Sox9", "Gja1", "C3"]


def score_genes(log_expr: pd.DataFrame, genes: list[str]) -> tuple[pd.Series, int]:
    present = [g for g in genes if g in log_expr.columns]
    if not present:
        return pd.Series(0.0, index=log_expr.index), 0
    return log_expr[present].mean(axis=1), len(present)


def load_sample(root: Path, sample: dict) -> pd.DataFrame:
    raw = root / "data" / "raw" / "D4_GSE263791"
    expr_path = raw / f"{sample['prefix']}_exprMat_file.csv.gz"
    meta_path = raw / f"{sample['prefix']}_metadata_file.csv.gz"
    expr = pd.read_csv(expr_path)
    meta = pd.read_csv(meta_path)
    key = ["fov", "cell_ID"]
    if expr[key].duplicated().any() or meta[key].duplicated().any():
        raise ValueError(f"duplicate fov/cell_ID in {sample['sample_id']}")
    feature_cols = [c for c in expr.columns if c not in key]
    neg_cols = [c for c in feature_cols if c.startswith("NegPrb")]
    gene_cols = [c for c in feature_cols if c not in neg_cols]
    counts = expr[gene_cols].apply(pd.to_numeric, errors="coerce").fillna(0.0)
    totals = counts.sum(axis=1)
    n_features = (counts > 0).sum(axis=1)
    neg_counts = expr[neg_cols].apply(pd.to_numeric, errors="coerce").fillna(0.0).sum(axis=1) if neg_cols else pd.Series(0.0, index=expr.index)
    log_expr = np.log1p(counts.div(totals.replace(0, np.nan), axis=0) * 1e4).fillna(0.0)
    micro, n_micro = score_genes(log_expr, MICROGLIA_MARKERS)
    astro, n_astro = score_genes(log_expr, ASTROCYTE_MARKERS)
    margin = micro - astro
    label = np.where((micro >= astro + 0.25) & (micro >= 0.2), "Microglia_provisional", np.where((astro >= micro + 0.25) & (astro >= 0.2), "Astrocyte_provisional", "Unresolved"))
    out = expr[key].copy()
    out["sample_id"] = sample["sample_id"]
    out["subject_label"] = sample["subject_label"]
    out["disease"] = sample["disease"]
    out["n_genes"] = len(gene_cols)
    out["nCount"] = totals
    out["nFeature"] = n_features
    out["negative_probe_count"] = neg_counts
    out["negative_probe_fraction"] = neg_counts / totals.replace(0, np.nan)
    out["microglia_marker_score"] = micro
    out["astrocyte_marker_score"] = astro
    out["marker_score_margin_micro_minus_astro"] = margin
    out["provisional_state"] = label
    out["marker_genes_present_microglia"] = n_micro
    out["marker_genes_present_astrocyte"] = n_astro
    # The expression export contains a small number of cell_ID=0 rows that are
    # not represented in the official metadata table. Keep only coordinate-
    # resolved cells and record the input audit in main().
    out = out.merge(meta, on=key, how="inner", validate="one_to_one", suffixes=("", "_metadata"))
    if out[["CenterX_global_px", "CenterY_global_px"]].isna().any().any():
        raise ValueError(f"missing spatial coordinates after metadata-resolved merge in {sample['sample_id']}")
    coords = out[["CenterX_global_px", "CenterY_global_px"]].to_numpy(float)
    k = min(7, len(out))
    tree = cKDTree(coords)
    _, neighbors = tree.query(coords, k=k)
    if k == 1:
        neighbors = neighbors[:, None]
    states = out["provisional_state"].to_numpy()
    micro_frac = []
    astro_frac = []
    for i, inds in enumerate(neighbors):
        inds = np.atleast_1d(inds)
        inds = inds[inds != i]
        vals = states[inds]
        micro_frac.append(float(np.mean(vals == "Microglia_provisional")) if len(vals) else np.nan)
        astro_frac.append(float(np.mean(vals == "Astrocyte_provisional")) if len(vals) else np.nan)
    out["k6_neighbor_microglia_fraction"] = micro_frac
    out["k6_neighbor_astrocyte_fraction"] = astro_frac
    out["beta_amyloid_percentile_within_sample"] = out["Mean.BetaAmyloid"].rank(pct=True)
    out["descriptive_qc_flag"] = np.where(out["nCount"] > 0, "retained_for_descriptive_support", "zero_expression")
    return out


def summarize(cell: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for sample_id, g in cell.groupby("sample_id", sort=True):
        micro = g[g.provisional_state == "Microglia_provisional"]
        astro = g[g.provisional_state == "Astrocyte_provisional"]
        rho_micro = spearmanr(g["beta_amyloid_percentile_within_sample"], g["microglia_marker_score"], nan_policy="omit")
        rho_astro = spearmanr(g["beta_amyloid_percentile_within_sample"], g["astrocyte_marker_score"], nan_policy="omit")
        rows.append(
            {
                "sample_id": sample_id,
                "subject_label": g["subject_label"].iloc[0],
                "disease": g["disease"].iloc[0],
                "n_cells": len(g),
                "n_microglia_provisional": len(micro),
                "n_astrocyte_provisional": len(astro),
                "n_unresolved": int((g.provisional_state == "Unresolved").sum()),
                "median_nCount": float(g.nCount.median()),
                "median_nFeature": float(g.nFeature.median()),
                "median_negative_probe_fraction": float(g.negative_probe_fraction.median()),
                "mean_microglia_score": float(g.microglia_marker_score.mean()),
                "mean_astrocyte_score": float(g.astrocyte_marker_score.mean()),
                "mean_beta_amyloid_microglia": float(micro["Mean.BetaAmyloid"].mean()) if len(micro) else np.nan,
                "mean_beta_amyloid_astrocyte": float(astro["Mean.BetaAmyloid"].mean()) if len(astro) else np.nan,
                "mean_neighbor_astrocyte_fraction_for_microglia": float(micro["k6_neighbor_astrocyte_fraction"].mean()) if len(micro) else np.nan,
                "mean_neighbor_microglia_fraction_for_astrocyte": float(astro["k6_neighbor_microglia_fraction"].mean()) if len(astro) else np.nan,
                "spearman_beta_vs_microglia_score": float(rho_micro.statistic),
                "spearman_beta_vs_microglia_p": float(rho_micro.pvalue),
                "spearman_beta_vs_astrocyte_score": float(rho_astro.statistic),
                "spearman_beta_vs_astrocyte_p": float(rho_astro.pvalue),
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    args = parser.parse_args()
    root = args.root.resolve()
    tables = root / "results" / "tables"
    interim = root / "data" / "interim"
    tables.mkdir(parents=True, exist_ok=True)
    interim.mkdir(parents=True, exist_ok=True)
    input_audit = []
    loaded = []
    for sample in SAMPLES:
        raw_dir = root / "data" / "raw" / "D4_GSE263791"
        expr_path = raw_dir / f"{sample['prefix']}_exprMat_file.csv.gz"
        meta_path = raw_dir / f"{sample['prefix']}_metadata_file.csv.gz"
        expr_keys = pd.read_csv(expr_path, usecols=["fov", "cell_ID"])
        meta_keys = pd.read_csv(meta_path, usecols=["fov", "cell_ID"])
        matched = expr_keys.merge(meta_keys, on=["fov", "cell_ID"], how="inner").shape[0]
        input_audit.append({"sample_id": sample["sample_id"], "expr_rows": len(expr_keys), "metadata_rows": len(meta_keys), "matched_coordinate_rows": matched, "unmatched_expr_rows": len(expr_keys) - matched, "disease": sample["disease"]})
        loaded.append(load_sample(root, sample))
    cells = pd.concat(loaded, ignore_index=True)
    summary = summarize(cells)
    pd.DataFrame(input_audit).to_csv(tables / "D4_input_audit.tsv", sep="\t", index=False)
    cells.to_csv(tables / "D4_cell_spatial_qc.tsv", sep="\t", index=False)
    summary.to_csv(tables / "D4_spatial_support_summary.tsv", sep="\t", index=False)
    provenance = {
        "dataset": "GSE263791",
        "input": "official processed CosMx expression, metadata and FOV position files for GSM8199188 and GSM8199189",
        "samples": SAMPLES,
        "input_audit": input_audit,
        "cells": int(len(cells)),
        "features_excluding_negative_probes": int(cells["n_genes"].iloc[0]),
        "marker_sets": {"microglia": MICROGLIA_MARKERS, "astrocyte": ASTROCYTE_MARKERS},
        "state_rule": "provisional marker-score labels; margin >= 0.25 and winning score >= 0.2; unresolved otherwise",
        "spatial_rule": "within-sample k=6 nearest-neighbor descriptive fractions using global pixel coordinates",
        "boundaries": [
            "both available sections are AD; no AD-control contrast is possible",
            "processed metadata has no authoritative cell-type labels; provisional labels are not cell annotations",
            "no ligand-receptor or causal communication inference is performed",
            "no hard QC deletion threshold is applied; QC fields are descriptive",
        ],
    }
    (interim / "D4_spatial_provenance.json").write_text(json.dumps(provenance, indent=2), encoding="utf-8")
    print(json.dumps({"cells": len(cells), "samples": len(summary), "summary": summary.to_dict(orient="records")}, indent=2))


if __name__ == "__main__":
    main()