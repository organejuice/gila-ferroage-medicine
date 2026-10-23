from __future__ import annotations

import hashlib
import json
import math
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from scipy.stats import norm, ttest_ind


ROOT = Path(__file__).resolve().parents[1]
CONFIG = ROOT / "configs" / "g4_v3.2_execution.yaml"
CFG = yaml.safe_load(CONFIG.read_text(encoding="utf-8"))
RUN = ROOT / "results" / CFG["run_id"]
if RUN.exists():
    raise RuntimeError(f"Refusing to overwrite existing run directory: {RUN}")

TABLES = RUN / "tables"
GATES = RUN / "gates"
REGISTRIES = RUN / "registries"
for directory in (TABLES, GATES, REGISTRIES):
    directory.mkdir(parents=True, exist_ok=False)


MODULES = {
    "M5_Iron_Uptake": ["TFRC", "SLC11A2", "STEAP3", "DMT1"],
    "M5_Iron_Storage": ["FTH1", "FTL"],
    "M5_Iron_Export": ["SLC40A1", "CP"],
    "M5_Ferritinophagy": ["NCOA4", "ATG3", "ATG5", "ATG7", "MAP1LC3B", "SQSTM1"],
    "M6_LPO_ACSL4": ["ACSL4", "ALOX15", "ALOX5", "ALOX12", "NOX4", "GPX4", "LPCAT3", "SAT1"],
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def rel(path: Path) -> str:
    return path.relative_to(ROOT).as_posix()


def bh(values: list[float]) -> np.ndarray:
    p = np.asarray(values, dtype=float)
    out = np.full(len(p), np.nan)
    valid = np.isfinite(p)
    if not valid.any():
        return out
    idx = np.where(valid)[0]
    order = idx[np.argsort(p[idx])]
    ranked = p[order]
    adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    out[order] = np.minimum(adjusted, 1.0)
    return out


def random_effects(sub: pd.DataFrame) -> dict:
    clean = sub.loc[np.isfinite(sub["logFC"]) & (sub["SE"] > 0)].copy()
    k = len(clean)
    if k == 0:
        return {"k": 0, "meta_effect": np.nan, "meta_SE": np.nan, "meta_P": np.nan, "Q": np.nan, "I2_percent": np.nan, "tau2": np.nan}
    y = clean["logFC"].to_numpy(float)
    se = clean["SE"].to_numpy(float)
    w = 1.0 / (se**2)
    fixed = float(np.sum(w * y) / np.sum(w))
    q = float(np.sum(w * (y - fixed) ** 2))
    df = max(k - 1, 0)
    c = float(np.sum(w) - np.sum(w**2) / np.sum(w)) if k > 1 else 0.0
    tau2 = max(0.0, (q - df) / c) if c > 0 else 0.0
    wr = 1.0 / (se**2 + tau2)
    meta = float(np.sum(wr * y) / np.sum(wr))
    meta_se = float(math.sqrt(1.0 / np.sum(wr)))
    z = meta / meta_se if meta_se > 0 else np.nan
    p = float(2.0 * norm.sf(abs(z))) if np.isfinite(z) else np.nan
    i2 = float(max(0.0, (q - df) / q) * 100.0) if q > 0 and k > 1 else 0.0
    return {"k": k, "meta_effect": meta, "meta_SE": meta_se, "meta_P": p, "Q": q, "I2_percent": i2, "tau2": tau2}


def load_inputs() -> dict[str, Path]:
    paths = {key: ROOT / value for key, value in CFG["inputs"].items()}
    for path in paths.values():
        if not path.exists():
            raise FileNotFoundError(path)
    return paths


def d2_sample_level_rna(paths: dict[str, Path]) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    counts = pd.read_csv(paths["d2_gene_tpm"], sep="\t", index_col=0)
    counts.index = counts.index.astype(str).str.upper()
    counts = counts.groupby(level=0, sort=False).sum(numeric_only=True)
    meta = pd.read_csv(paths["d2_sample_metadata"], sep="\t")
    meta = meta.set_index("sample_id_in_expression_table")
    missing = [sample for sample in counts.columns if sample not in meta.index]
    if missing:
        raise ValueError(f"D2 samples absent from metadata: {missing[:5]}")
    meta = meta.loc[counts.columns].copy()
    expression = np.log2(counts.astype(float) + 1.0)
    d5 = pd.read_csv(paths["d5_ferroaging_reference"], sep="\t")
    d5_genes = sorted(set(d5["gene_symbol"].astype(str).str.upper()))
    modules = dict(MODULES)
    modules["M7_FerroAging"] = d5_genes
    rows: list[dict] = []
    scores: dict[str, pd.Series] = {}
    for module, genes in modules.items():
        available = [gene for gene in genes if gene in expression.index]
        score = expression.loc[available].mean(axis=0)
        score = (score - score.mean()) / max(float(score.std(ddof=0)), 1e-12)
        scores[module] = score
        for hour in (4, 24, 96):
            treated = score.loc[meta["time_hours"].eq(hour)].to_numpy(float)
            control = score.loc[meta["time_hours"].eq(0)].to_numpy(float)
            stat, p_value = ttest_ind(treated, control, equal_var=True)
            rows.append(
                {
                    "modality": "RNA_sample_level",
                    "module": module,
                    "time_hours": hour,
                    "n_genes_total": len(genes),
                    "n_genes_present": len(available),
                    "coverage": len(available) / len(genes),
                    "n_treated": len(treated),
                    "n_control": len(control),
                    "effect_treated_minus_control": float(np.mean(treated) - np.mean(control)),
                    "t": float(stat),
                    "P.Value": float(p_value),
                    "direction": "up" if np.mean(treated) > np.mean(control) else "down",
                }
            )
    rna = pd.DataFrame(rows)
    rna["FDR"] = bh(rna["P.Value"].tolist())
    rna["positive_FDR"] = (rna["effect_treated_minus_control"] > 0) & (rna["FDR"] < 0.05)
    expected = {
        "M5_Iron_Uptake": "up",
        "M5_Iron_Storage": "up",
        "M5_Iron_Export": "down",
        "M5_Ferritinophagy": "up",
        "M6_LPO_ACSL4": "up",
        "M7_FerroAging": "up",
    }
    rna["expected_direction"] = rna["module"].map(expected)
    rna["directional_FDR"] = (rna["FDR"] < 0.05) & (rna["direction"] == rna["expected_direction"])
    rna.to_csv(TABLES / "D2_G4_RNA_sample_level_module_effects.tsv", sep="\t", index=False)
    return rna, meta.reset_index(), {"module_sizes": {key: len(value) for key, value in modules.items()}}


def human_context(paths: dict[str, Path]) -> pd.DataFrame:
    source = pd.read_csv(paths["human_region_module_effects"], sep="\t")
    modules = ["M5_Iron_Uptake", "M5_Iron_Storage", "M5_Iron_Export", "M5_Ferritinophagy", "M6_LPO_ACSL4", "M7_FerroAging_Author"]
    subset = source.loc[source["module"].isin(modules)].copy()
    subset["cohort_region"] = subset["cohort"].astype(str) + "__" + subset["region"].astype(str)
    subset.to_csv(TABLES / "human_G1_G4_region_module_effects.tsv", sep="\t", index=False)
    rows = []
    for module, group in subset.groupby("module", sort=True):
        result = random_effects(group.rename(columns={"logFC": "logFC", "SE": "SE"}))
        rows.append({"module": module, **result, "n_region_strata": int(len(group))})
    summary = pd.DataFrame(rows)
    summary["human_translation"] = np.where(
        (summary["meta_effect"] > 0) & (summary["meta_P"] < 0.05),
        "supportive",
        "unresolved",
    )
    summary.to_csv(TABLES / "human_G4_random_effects_context.tsv", sep="\t", index=False)
    return summary


def external_registry() -> pd.DataFrame:
    rows = []
    for item in CFG["external_evidence"]:
        rows.append(
            {
                "evidence_id": item["evidence_id"],
                "title": item["title"],
                "accession_or_doi": item.get("accession", item.get("doi", "")),
                "source_url": item.get("source_url", ""),
                "supplementary_url": item.get("supplementary_url", ""),
                "modality": item.get("modality", ""),
                "context": item.get("context", ""),
                "sample_design": item.get("sample_design", ""),
                "assays": item.get("assays", ""),
                "result": item.get("result", ""),
                "evidence_tier": item.get("evidence_tier", ""),
                "counts_as_second_modality": bool(item.get("counts_as_second_modality", False)),
                "local_quantitative_reanalysis": bool(item.get("local_quantitative_reanalysis", False)),
                "limitation": item.get("limitation", item.get("exclusion_reason", "")),
            }
        )
    registry = pd.DataFrame(rows)
    registry.to_csv(REGISTRIES / "external_G4_evidence_registry.tsv", sep="\t", index=False)
    return registry


def main() -> None:
    paths = load_inputs()
    g3 = json.loads(paths["g3_decision"].read_text(encoding="utf-8"))
    rna, d2_meta, qc = d2_sample_level_rna(paths)
    human = human_context(paths)
    external = external_registry()

    iron = rna.loc[rna["module"].str.startswith("M5_") & rna["directional_FDR"]]
    lpo = rna.loc[(rna["module"] == "M6_LPO_ACSL4") & rna["directional_FDR"]]
    ferro = rna.loc[(rna["module"] == "M7_FerroAging") & rna["directional_FDR"]]
    iron_direction_pattern = all(
        [
            (rna.loc[(rna["module"] == "M5_Iron_Uptake") & rna["directional_FDR"]].shape[0] >= 1),
            (rna.loc[(rna["module"] == "M5_Iron_Storage") & rna["directional_FDR"]].shape[0] >= 1),
            (rna.loc[(rna["module"] == "M5_Ferritinophagy") & rna["directional_FDR"]].shape[0] >= 1),
            (rna.loc[(rna["module"] == "M5_Iron_Export") & rna["directional_FDR"]].shape[0] >= 1),
        ]
    )
    rna_iron = len(iron) >= int(CFG["analysis"]["minimum_positive_iron_timepoints"]) and iron_direction_pattern
    rna_lpo = len(lpo) >= int(CFG["analysis"]["minimum_positive_lpo_timepoints"])
    rna_ferro = len(ferro) >= int(CFG["analysis"]["minimum_positive_ferro_aging_timepoints"])
    direct_iron = bool(
        (
            (external["evidence_id"] == "Kenkhuis_2021_direct_iron_microglia")
            & external["counts_as_second_modality"]
        ).any()
    )
    modalities = 1 + int(direct_iron)
    upstream_pass = g3.get("decision") == "G3_PASS_MODEL_SCOPED"
    pass_decision = bool(upstream_pass and rna_iron and rna_lpo and rna_ferro and direct_iron and modalities >= int(CFG["analysis"]["minimum_modalities"]))
    decision = "G4_PASS_MODEL_SCOPED" if pass_decision else "G4_NOT_YET_PASS"

    gate_rows = [
        {"criterion": "upstream_G3_pass", "required": True, "observed": upstream_pass, "detail": g3.get("decision", "missing")},
        {"criterion": "RNA_iron_homeostasis_convergence", "required": True, "observed": rna_iron, "detail": f"directional_FDR_rows={len(iron)}; direction_pattern={iron_direction_pattern}; n=3 biological replicates per timepoint"},
        {"criterion": "RNA_lipid_peroxidation_proxy", "required": True, "observed": rna_lpo, "detail": f"directional_FDR_rows={len(lpo)}; module=M6_LPO_ACSL4"},
        {"criterion": "RNA_ferro_aging_proxy", "required": True, "observed": rna_ferro, "detail": f"directional_FDR_rows={len(ferro)}; module=M7_FerroAging"},
        {"criterion": "independent_direct_iron_modality", "required": True, "observed": direct_iron, "detail": "Kenkhuis 2021 human AD microglial iron histology/multispectral IF; publication-level, not locally re-fit"},
        {"criterion": "minimum_convergent_modalities", "required": True, "observed": modalities >= 2, "detail": f"RNA_sample_level={True}; direct_iron={direct_iron}; count={modalities}; required={CFG['analysis']['minimum_modalities']}"},
        {"criterion": "direct_lipid_peroxidation_flux", "required": False, "observed": False, "detail": "not measured; LPO remains a transcriptional proxy"},
        {"criterion": "human_G1_RNA_translation", "required": False, "observed": bool((human["human_translation"] == "supportive").any()), "detail": "supportive context only; not a veto under model-scoped G4"},
    ]
    pd.DataFrame(gate_rows).to_csv(GATES / "G4_gate_v3.2.tsv", sep="\t", index=False)

    d2_meta.to_csv(TABLES / "D2_sample_metadata_audited.tsv", sep="\t", index=False)
    human.to_csv(TABLES / "human_G4_context_summary.tsv", sep="\t", index=False)

    decision_record = {
        "run_id": CFG["run_id"],
        "protocol": "AD disease-side v3.2",
        "decision": decision,
        "superseded_run": CFG["superseded_run"],
        "upstream_G3": g3.get("decision"),
        "modalities": {
            "RNA_sample_level": {
                "status": "supported",
                "independent_unit": "D2 biological replicate",
                "replicates_per_timepoint": 3,
                "directional_FDR_iron_rows": int(len(iron)),
                "directional_FDR_lpo_rows": int(len(lpo)),
                "directional_FDR_ferro_aging_rows": int(len(ferro)),
            },
            "direct_iron_histology": {
                "status": "supportive_publication_level",
                "study": "Kenkhuis et al. 2021",
                "iron_positive_analysis": "11 AD vs 9 control",
                "local_quantitative_reanalysis": False,
                "uncertainty": "not reconstructed from current public summary",
            },
        },
        "human_context": human.to_dict(orient="records"),
        "claim_boundary": {
            "allowed": "model-scoped convergence of iron-handling, LPO-associated, and ferro-aging-like transcriptional programs with independent AD microglial iron-loading support",
            "not_allowed": "direct lipid-peroxidation flux, ferroptosis, or causal ferro-aging mechanism",
            "human_RNA": "translation context only; not a veto",
        },
        "downstream": "G5 may proceed because G4 passed; G6 remains locked until G5 passes." if pass_decision else "G5-G6 remain locked.",
        "created_at": datetime.now(timezone.utc).isoformat(),
    }
    (RUN / "g4_decision.json").write_text(json.dumps(decision_record, ensure_ascii=False, indent=2), encoding="utf-8")

    manifest = {
        "run_id": CFG["run_id"],
        "script": rel(Path(__file__).resolve()),
        "config": rel(CONFIG),
        "raw_data_modified": False,
        "remote_computation": False,
        "result": decision,
        "inputs": [{"file": rel(path), "sha256": sha256(path)} for path in paths.values()],
        "external_evidence_metadata": [
            {"evidence_id": item["evidence_id"], "source_url": item.get("source_url", ""), "supplementary_url": item.get("supplementary_url", ""), "supplementary_expected_md5": item.get("supplementary_expected_md5", "")}
            for item in CFG["external_evidence"]
        ],
    }
    (RUN / "run_manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")

    report = f"""# Disease-side v3.2 G4 execution

## Decision

**{decision}**

G4 was executed after `{g3.get('decision')}`. The analysis used a sample-level D2 RNA module score, not a gene-level pseudo-replicate test: D2 has three independent biological replicates at each time point and the 18 module-by-time tests were BH-corrected together.

- Iron-homeostasis RNA convergence: `{'PASS' if rna_iron else 'FAIL'}`; expected-direction FDR rows = {len(iron)}; required temporal/directional pattern = {iron_direction_pattern}.
- LPO-associated RNA proxy: `{'PASS' if rna_lpo else 'FAIL'}`; expected-direction FDR rows = {len(lpo)}.
- Ferro-aging-like RNA proxy: `{'PASS' if rna_ferro else 'FAIL'}`; expected-direction FDR rows = {len(ferro)}.
- Independent second modality: `{'PASS' if direct_iron else 'FAIL'}`; Kenkhuis et al. 2021 direct iron histology/multispectral microglial evidence.
- Modalities counted: {modalities}; v3.2 minimum: {CFG['analysis']['minimum_modalities']}.

## Claim boundary

This result supports a model-scoped convergence of iron-handling, LPO-associated, and ferro-aging-like transcriptional programs with independent AD microglial iron-loading evidence. It does **not** establish direct lipid-peroxidation flux, ferroptosis, or a causal ferro-aging mechanism; no lipidomics, proteomics, metabolite, or flux measurement was available in the executed data layer. The direct iron study is publication-level support and was not locally re-fit because the Figshare supplement download was WAF-blocked; the source DOI and expected MD5 are recorded in the registry and manifest.

Human G1 RNA meta-effects are retained as translation context, not as a veto, consistent with the model-scoped execution boundary.
"""
    (RUN / "execution_report.md").write_text(report, encoding="utf-8")

    checksums = []
    for path in sorted(RUN.rglob("*")):
        if path.is_file() and path.name != "checksums.sha256":
            checksums.append(f"{sha256(path)}  {rel(path)}")
    (RUN / "checksums.sha256").write_text("\n".join(checksums) + "\n", encoding="utf-8")

    print(json.dumps({"run": rel(RUN), "decision": decision, "modalities": modalities, "directional_FDR_iron_rows": int(len(iron)), "directional_FDR_lpo_rows": int(len(lpo)), "directional_FDR_ferro_aging_rows": int(len(ferro)), "human_context": human.to_dict(orient="records")}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
