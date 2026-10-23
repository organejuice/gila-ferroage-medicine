#!/usr/bin/env python3
"""Build the auditable D0-D2 disease-facing package for drug-side analysis."""
from __future__ import annotations

import argparse
import gzip
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import t as t_dist


MODULES = {
    "M1_Glycolysis": ["HK1", "HK2", "PFKFB3", "PFKFB4", "ALDOA", "ALDOC", "GAPDH", "PGK1", "PGAM1", "ENO1", "ENO2", "PKM", "LDHA", "PDK1", "SLC2A1", "HIF1A", "PFKP"],
    "M3_SASP_Senescence": ["CDKN1A", "CDKN2A", "IL6", "IL1B", "TNF", "CCL2", "CCL3", "CCL4", "CSF1", "CXCL10", "MMP9", "MMP14", "SERPINE1", "IGFBP3", "IGFBP7", "NFKB1", "RELA", "TGFB1"],
    "M5_Iron_Uptake": ["TFRC", "SLC11A2", "STEAP3", "DMT1"],
    "M5_Iron_Storage": ["FTH1", "FTL"],
    "M5_Iron_Export": ["SLC40A1", "CP"],
    "M5_Ferritinophagy": ["NCOA4", "ATG3", "ATG5", "ATG7", "MAP1LC3B", "SQSTM1"],
    "M6_LPO_ACSL4": ["ACSL4", "ALOX15", "ALOX5", "ALOX12", "NOX4", "GPX4", "LPCAT3", "SAT1"],
}


def read_matrix(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", compression="gzip", index_col=0)


def bh(p_values: np.ndarray) -> np.ndarray:
    p = np.asarray(p_values, dtype=float)
    out = np.full_like(p, np.nan)
    ok = np.isfinite(p)
    if not ok.any():
        return out
    order = np.argsort(p[ok])
    ranked = p[ok][order]
    adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    values = np.empty_like(adjusted)
    values[order] = np.minimum(adjusted, 1.0)
    out[ok] = values
    return out


def design_matrix(meta: pd.DataFrame) -> tuple[np.ndarray, list[str]]:
    disease = (meta["disease"].astype(str) == "AD").astype(float).to_numpy()
    region = (meta["region"].astype(str) == "SSC").astype(float).to_numpy()
    sex = (meta["sex"].astype(str).str.lower() == "male").astype(float).to_numpy()
    age = meta["age"].astype(float).to_numpy()
    rin = meta["rin"].astype(float).to_numpy()
    age = (age - np.nanmean(age)) / np.nanstd(age)
    rin = (rin - np.nanmean(rin)) / np.nanstd(rin)
    X = np.column_stack([np.ones(len(meta)), disease, region, age, sex, rin])
    return X, ["Intercept", "diseaseAD", "regionSSC", "age_z", "sexMale", "rin_z"]


def fit_effect(values: np.ndarray, X: np.ndarray, coef_name: str = "diseaseAD") -> dict:
    n, p = X.shape
    beta = np.linalg.lstsq(X, values, rcond=None)[0]
    residual = values - X @ beta
    dof = max(n - p, 1)
    sigma2 = np.sum(residual**2, axis=0) / dof
    cov = np.linalg.pinv(X.T @ X)
    se = np.sqrt(np.maximum(sigma2 * cov[1, 1], 1e-15))
    t_value = beta[1] / se
    p_value = 2.0 * t_dist.sf(np.abs(t_value), dof)
    return {"beta": beta[1], "se": se, "t": t_value, "p": p_value, "dof": dof}


def module_effects(log_expr: pd.DataFrame, meta: pd.DataFrame, X: np.ndarray, cell_type: str, modules: dict[str, list[str]]) -> pd.DataFrame:
    rows = []
    gene_upper = {str(g).upper(): g for g in log_expr.index}
    for module, genes in modules.items():
        available = [gene_upper[g] for g in genes if g in gene_upper]
        if len(available) < 2:
            rows.append({"cell_type": cell_type, "module": module, "n_total": len(genes), "n_usable": len(available), "coverage_ratio": len(available) / len(genes), "beta": np.nan, "se": np.nan, "t": np.nan, "p": np.nan, "FDR": np.nan, "status": "insufficient_gene_coverage"})
            continue
        score = log_expr.loc[available].mean(axis=0).to_numpy()
        score = (score - np.mean(score)) / max(np.std(score), 1e-12)
        effect = fit_effect(score[:, None], X)
        rows.append({"cell_type": cell_type, "module": module, "n_total": len(genes), "n_usable": len(available), "coverage_ratio": len(available) / len(genes), **{k: float(np.asarray(v).ravel()[0]) for k, v in effect.items() if k != "dof"}, "dof": effect["dof"], "status": "estimated"})
    out = pd.DataFrame(rows)
    out["FDR"] = bh(out["p"].to_numpy())
    return out


def run(root: Path) -> None:
    out_tables = root / "results" / "tables"
    interim = root / "data" / "interim"
    out_tables.mkdir(parents=True, exist_ok=True)
    interface: dict = {"package": "Disease Interface Package v1.3", "source_dataset": "D3_GSE160936", "discovery_status": "D3_discovery_only; D6_expression_validation_unavailable", "normalization_claim_allowed": False, "modules": {}}
    module_sets = dict(MODULES)
    d5_path = root / "data" / "raw" / "D5_CellMetab2026" / "Table_S2_ferro_aging_reference.xlsx"
    if d5_path.exists():
        try:
            d5 = pd.read_excel(d5_path, sheet_name="Ferroaging gene set", header=1)
            symbol_col = next((c for c in d5.columns if str(c).lower() in {"gene", "symbol", "gene_symbol", "genesymbol"}), None)
            if symbol_col is not None:
                symbols = sorted({str(x).strip().upper() for x in d5[symbol_col].dropna() if str(x).strip()})
                if symbols:
                    module_sets["M7_FerroAging_Author"] = symbols
        except Exception as exc:
            interface["M7_reference_warning"] = f"D5_Table_S2_read_failed: {type(exc).__name__}: {exc}"
    module_rows = []
    signature_paths = {}
    for cell_type in ("microglia", "astrocyte"):
        matrix_path = interim / f"D3_pseudobulk_counts_{cell_type}.tsv.gz"
        meta_path = interim / f"D3_pseudobulk_metadata_{cell_type}.tsv"
        counts = read_matrix(matrix_path)
        meta = pd.read_csv(meta_path, sep="\t").set_index("group_id").loc[counts.columns]
        lib = counts.sum(axis=0).replace(0, np.nan)
        log_expr = np.log2(counts.div(lib, axis=1) * 1_000_000.0 + 1.0)
        X, names = design_matrix(meta)
        keep_genes = (counts >= 10).sum(axis=1).to_numpy() >= 3
        expr = log_expr.loc[keep_genes]
        effects = fit_effect(expr.to_numpy().T, X)
        sig = pd.DataFrame({"gene": expr.index.astype(str), "beta": effects["beta"], "se": effects["se"], "t": effects["t"], "P.Value": effects["p"]})
        sig["FDR"] = bh(sig["P.Value"].to_numpy())
        sig["direction"] = np.where(sig["beta"] > 0, "AD_up", "AD_down")
        sig["cell_type"] = cell_type.capitalize()
        sig["n_groups"] = len(meta)
        sig["n_AD"] = int((meta["disease"] == "AD").sum())
        sig["n_Control"] = int((meta["disease"] == "Control").sum())
        sig = sig.sort_values(["FDR", "t"], ascending=[True, False])
        path = out_tables / f"DiseaseSignature_{cell_type}_fast.tsv"
        sig.to_csv(path, sep="\t", index=False)
        signature_paths[cell_type] = str(path.relative_to(root))
        effects_modules = module_effects(expr, meta, X, cell_type.capitalize(), module_sets)
        effects_modules.to_csv(out_tables / f"DiseaseModuleEffects_{cell_type}_fast.tsv", sep="\t", index=False)
        module_rows.append(effects_modules)

    all_modules = pd.concat(module_rows, ignore_index=True)
    all_modules.to_csv(out_tables / "DiseaseModuleEffects_fast.tsv", sep="\t", index=False)
    for module in sorted(all_modules["module"].unique()):
        sub = all_modules[all_modules["module"] == module]
        interface["modules"][module] = {"effects": sub.to_dict(orient="records"), "therapeutic_qualification": "TQ-B_or_TQ-D_pending_external_validation", "desired_state": "directional_compatibility_only; D3_Control_reference; strict_normalization_not_allowed"}
    interface["interfaces"] = {"I1_signed_disease_signature": {"status": "available_D3_discovery", "files": signature_paths, "statistical_unit": "donor x region x glial cell type pseudobulk"}, "I2_signed_regulatory_activity": {"status": "not_executed", "reason": "regulon/pathway activity resources and independent validation not yet available"}, "I3_7D_module_vector": {"status": "partial", "missing": ["M2_H4K12la_target_from_D1_peak_to_gene", "M4_same-system_communication_validation"], "files": ["results/tables/DiseaseModuleEffects_fast.tsv"], "M7_source": "D5_CellMetab2026_Table_S2_if_available"}}
    interface["D1_D2_freeze"] = {"freeze_before_drug_ranking": True, "TQ_rule": "M1/M3/M5/M6/M7 remain context-dependent or descriptive until functional/external validation; no automatic reversal of TQ-C/D programs", "reference": "D3 Control groups in the same discovery dataset; not a strict normalization reference", "primary_metric": "signed disease beta and module standardized beta; coverage and FDR reported"}
    (interim / "disease_interface_package_v1.3.json").write_text(json.dumps(interface, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps({"signatures": signature_paths, "module_rows": len(all_modules), "strict_normalization_allowed": False}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    args = parser.parse_args()
    run(args.root.resolve())
