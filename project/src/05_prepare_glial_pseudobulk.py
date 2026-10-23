#!/usr/bin/env python3
"""Prepare sample-aware glial pseudobulk matrices without a full Seurat merge.

This is an execution-oriented fallback for the D3 discovery input.  It keeps
the locked QC thresholds and marker definitions, but reads one Cell Ranger
matrix at a time and aggregates raw counts at donor x region x cell-type
level.  No raw files are modified.
"""
from __future__ import annotations

import argparse
import gzip
import json
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.io import mmread
import yaml


MARKERS = {
    "Microglia": ["AIF1", "CX3CR1", "P2RY12", "TMEM119", "C1QA", "C1QB", "C1QC", "TYROBP", "LST1", "SPI1"],
    "Astrocyte": ["AQP4", "GFAP", "SLC1A2", "SLC1A3", "ALDH1L1", "SOX9", "GJA1", "CLU"],
    "Neuron": ["SNAP25", "RBFOX3", "SYT1", "STMN2", "GRIN1"],
    "Oligodendrocyte": ["MBP", "PLP1", "MOG", "MOBP", "OLIG1", "OLIG2"],
    "OPC": ["PDGFRA", "CSPG4", "VCAN"],
    "Endothelial": ["CLDN5", "KDR", "EMCN", "FLT1"],
    "Pericyte": ["RGS5", "PDGFRB", "CSPG4", "MCAM"],
}


def read_lines(path: Path) -> list[str]:
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        return [line.rstrip("\n") for line in handle]


def read_10x(matrix_dir: Path) -> tuple[sp.csc_matrix, list[str]]:
    features = read_lines(matrix_dir / "features.tsv.gz")
    genes = [row.split("\t")[1] if "\t" in row else row for row in features]
    genes = pd.Index(genes).astype(str)
    with gzip.open(matrix_dir / "matrix.mtx.gz", "rb") as handle:
        matrix = mmread(handle).tocsr()
    if matrix.shape[0] != len(genes):
        raise ValueError(f"feature/matrix mismatch in {matrix_dir}: {len(genes)} vs {matrix.shape[0]}")
    return matrix.tocsc(), genes.tolist()


def marker_scores(counts: sp.csc_matrix, genes: list[str], keep: np.ndarray) -> tuple[pd.DataFrame, pd.DataFrame]:
    gene_index = {gene.upper(): i for i, gene in enumerate(genes)}
    sub = counts[:, keep]
    totals = np.asarray(sub.sum(axis=0)).ravel()
    totals[totals <= 0] = 1.0
    lognorm = sub.multiply(10000.0 / totals).log1p().tocsr()
    scores: dict[str, np.ndarray] = {}
    available: dict[str, int] = {}
    for name, markers in MARKERS.items():
        idx = [gene_index[g] for g in markers if g.upper() in gene_index]
        available[name] = len(idx)
        scores[name] = np.asarray(lognorm[idx, :].mean(axis=0)).ravel() if idx else np.zeros(sub.shape[1])
    score_df = pd.DataFrame(scores)
    order = np.argsort(-score_df.to_numpy(), axis=1)
    best = score_df.columns.to_numpy()[order[:, 0]]
    second = score_df.to_numpy()[np.arange(len(score_df)), order[:, 1]]
    margin = score_df.to_numpy()[np.arange(len(score_df)), order[:, 0]] - second
    labels = np.where(margin >= 0.05, best, "Unresolved")
    label_df = pd.DataFrame({"cell_type": labels, "best_score": score_df.max(axis=1), "margin": margin})
    return score_df, label_df


def sample_qc(counts: sp.csc_matrix, genes: list[str], config: dict) -> tuple[np.ndarray, dict]:
    n_features = np.asarray((counts > 0).sum(axis=0)).ravel()
    n_counts = np.asarray(counts.sum(axis=0)).ravel()
    mt_idx = np.array([g.upper().startswith("MT-") for g in genes], dtype=bool)
    mt_counts = np.asarray(counts[mt_idx, :].sum(axis=0)).ravel() if mt_idx.any() else np.zeros(counts.shape[1])
    mt_pct = 100.0 * mt_counts / np.maximum(n_counts, 1.0)
    max_features = float(np.quantile(n_features, float(config["max_feature_quantile"])))
    max_mt = float(np.quantile(mt_pct, float(config["max_mt_quantile"])))
    keep = (n_features >= int(config["min_features"])) & (n_features <= max_features) & (mt_pct <= max_mt)
    return keep, {
        "n_cells_before": int(counts.shape[1]),
        "n_cells_after": int(keep.sum()),
        "retention_fraction": float(keep.mean()),
        "median_features": float(np.median(n_features)),
        "median_umis": float(np.median(n_counts)),
        "median_mt_percent": float(np.median(mt_pct)),
        "max_features": max_features,
        "max_mt_percent": max_mt,
    }


def write_matrix(path: Path, matrix: sp.csc_matrix, genes: list[str], columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
        handle.write("gene\t" + "\t".join(columns) + "\n")
        for i, gene in enumerate(genes):
            vals = np.asarray(matrix.getrow(i).todense()).ravel()
            handle.write(gene + "\t" + "\t".join(str(int(v)) for v in vals) + "\n")


def run(root: Path) -> None:
    config = yaml.safe_load((root / "configs" / "d3_params.yaml").read_text(encoding="utf-8"))["d3"]
    metadata = pd.read_csv(root / "data" / "interim" / "D3_sample_metadata.tsv", sep="\t")
    out_tables = root / "results" / "tables"
    out_interim = root / "data" / "interim"
    out_tables.mkdir(parents=True, exist_ok=True)
    out_interim.mkdir(parents=True, exist_ok=True)
    qc_rows: list[dict] = []
    cell_rows: list[dict] = []
    aggregate: dict[str, dict[str, sp.csc_matrix]] = {"Microglia": {}, "Astrocyte": {}}
    aggregate_genes: list[str] | None = None
    score_availability: dict[str, int] = {}

    for row in metadata.to_dict(orient="records"):
        matrix_dir = root / "data" / "interim" / "D3_matrices" / row["geo_sample"] / "filtered_feature_bc_matrix"
        counts, genes = read_10x(matrix_dir)
        keep, qc = sample_qc(counts, genes, config)
        qc_rows.append({"sample_id": row["geo_sample"], **qc})
        score_df, labels = marker_scores(counts, genes, keep)
        score_availability = {key: int(value) for key, value in score_availability.items()}
        for name, markers in MARKERS.items():
            score_availability[name] = len([g for g in markers if g.upper() in {x.upper() for x in genes}])
        label_counts = labels["cell_type"].value_counts().to_dict()
        for cell_type, n_cells in label_counts.items():
            cell_rows.append({"sample_id": row["geo_sample"], "donor_id": row["donor_id"], "disease": row["disease"], "region": row["region"], "cell_type": cell_type, "n_cells": int(n_cells)})
        if aggregate_genes is None:
            aggregate_genes = genes
        elif aggregate_genes != genes:
            raise ValueError("D3 feature ordering differs across samples; harmonization is required before aggregation")
        kept_indices = np.flatnonzero(keep)
        for cell_type in ("Microglia", "Astrocyte"):
            selected = kept_indices[labels["cell_type"].to_numpy() == cell_type]
            if selected.size == 0:
                continue
            group = f"{row['donor_id']}__{row['region']}__{cell_type}"
            summed = sp.csc_matrix(counts[:, selected].sum(axis=1))
            aggregate[cell_type][group] = summed

    if aggregate_genes is None:
        raise ValueError("No D3 matrices were processed")
    pd.DataFrame(qc_rows).to_csv(out_tables / "D3_fast_qc_audit.tsv", sep="\t", index=False)
    pd.DataFrame(cell_rows).to_csv(out_tables / "D3_fast_celltype_summary.tsv", sep="\t", index=False)
    for cell_type, groups in aggregate.items():
        group_names = sorted(groups)
        matrix = sp.hstack([groups[g] for g in group_names], format="csc")
        write_matrix(out_interim / f"D3_pseudobulk_counts_{cell_type.lower()}.tsv.gz", matrix, aggregate_genes, group_names)
        rows = []
        for group in group_names:
            donor, region, _ = group.split("__")
            meta = metadata[(metadata["donor_id"] == donor) & (metadata["region"] == region)].iloc[0]
            n_cells = int(next(item["n_cells"] for item in cell_rows if item["sample_id"] == meta["geo_sample"] and item["cell_type"] == cell_type))
            rows.append({"group_id": group, "donor_id": donor, "disease": meta["disease"], "region": region, "cell_type": cell_type, "age": meta["age"], "sex": meta["sex"], "rin": meta["rin"], "ptau_percent_positive_cells": meta["ptau_percent_positive_cells"], "amyloid_percent_area": meta["amyloid_percent_area"], "n_cells": n_cells})
        pd.DataFrame(rows).to_csv(out_interim / f"D3_pseudobulk_metadata_{cell_type.lower()}.tsv", sep="\t", index=False)
    summary = {"samples": int(len(metadata)), "qc_rows": int(len(qc_rows)), "pseudobulk_groups": {key: len(value) for key, value in aggregate.items()}, "marker_availability": score_availability, "annotation_method": "sample-wise QC plus marker-score assignment on log1p library-size normalized counts; unresolved cells excluded from primary glial pseudobulk", "raw_data_modified": False}
    (out_interim / "D3_fast_pseudobulk_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    args = parser.parse_args()
    run(args.root.resolve())
