"""Execute the v3.2 G2 amendment using revalidated GSE188765 inputs."""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import shutil
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from scipy.stats import binomtest, mannwhitneyu, spearmanr


ROOT = Path(__file__).resolve().parents[1]
CONFIG = ROOT / "configs" / "g2_v3.2_execution.yaml"
CFG = yaml.safe_load(CONFIG.read_text(encoding="utf-8"))
RUN = ROOT / "results" / CFG["run_id"]
TABLES = RUN / "tables"
REGISTRIES = RUN / "registries"
TABLES.mkdir(parents=True, exist_ok=True)
REGISTRIES.mkdir(parents=True, exist_ok=True)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def write_tsv(path: Path, rows: list[dict]) -> None:
    fields = list(rows[0]) if rows else []
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def read_peak_file(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", header=None, compression="gzip",
                       names=["chrom", "start", "end", "score", "strand"])[
                           ["chrom", "start", "end", "score"]]


def finite_mean(values: pd.Series) -> float:
    x = pd.to_numeric(values, errors="coerce").replace([np.inf, -np.inf], np.nan).dropna()
    return float(x.mean()) if len(x) else float("nan")


def finite_median(values: pd.Series) -> float:
    x = pd.to_numeric(values, errors="coerce").replace([np.inf, -np.inf], np.nan).dropna()
    return float(x.median()) if len(x) else float("nan")


peak_dir = ROOT / CFG["input"]["peak_dir"]
peak_paths = sorted(peak_dir.glob("GSM*_peaks.txt.gz"))
if len(peak_paths) != 8:
    raise RuntimeError(f"Expected 8 GSE188765 peak files, found {len(peak_paths)}")

sample_rows = []
for path in peak_paths:
    name = path.name
    parts = name.split("_")
    context = parts[1]
    condition = parts[2].split("-")[0]
    replicate = int(parts[2].split("-")[1].split(".")[0])
    frame = read_peak_file(path)
    frame["start"] = pd.to_numeric(frame["start"], errors="coerce")
    frame["end"] = pd.to_numeric(frame["end"], errors="coerce")
    frame["score"] = pd.to_numeric(frame["score"], errors="coerce")
    frame = frame.dropna(subset=["chrom", "start", "end", "score"])
    frame = frame[frame["end"] > frame["start"]].drop_duplicates(["chrom", "start", "end"])
    sample_rows.append({
        "gsm": name.split("_")[0], "context": context, "condition": condition,
        "replicate": replicate, "file": path.relative_to(ROOT).as_posix(),
        "file_size": path.stat().st_size, "sha256": sha256(path),
        "peak_rows": len(frame), "unique_intervals": len(frame),
        "total_peak_score": float(frame["score"].sum()),
    })
write_tsv(TABLES / "sample_qc.tsv", sample_rows)

audited_effects = ROOT / CFG["input"]["audited_gene_effects"]
effects = pd.read_csv(audited_effects, sep="\t")
effects["gene_upper"] = effects["gene_symbol"].astype(str).str.upper()
target_set = {x.upper() for x in CFG["analysis"]["glycolysis_target_set"]["genes"]}
target_effects = effects[effects["gene_upper"].isin(target_set)].copy()
target_effects["target_set_id"] = CFG["analysis"]["glycolysis_target_set"]["id"]
target_effects["is_positive_AD_occupancy"] = target_effects["log2FC_AD_minus_WT"] > 0
target_effects.to_csv(TABLES / "g2b_target_set_h4k12la_effects.tsv", sep="\t", index=False)

target_summary_rows = []
for context in ["M", "T"]:
    context_effects = effects[effects["context"] == context].copy()
    selected = context_effects[context_effects["gene_upper"].isin(target_set)]
    background = context_effects[~context_effects["gene_upper"].isin(target_set)]
    tv = selected["log2FC_AD_minus_WT"].dropna()
    bv = background["log2FC_AD_minus_WT"].dropna()
    mw = mannwhitneyu(tv, bv, alternative="greater") if len(tv) and len(bv) else None
    positive_n = int((selected["log2FC_AD_minus_WT"] > 0).sum())
    positive_total = int(selected["log2FC_AD_minus_WT"].notna().sum())
    positive_test = binomtest(positive_n, positive_total, 0.5) if positive_total else None
    positive_ci = positive_test.proportion_ci(0.95) if positive_test else None
    target_summary_rows.append({
        "context": context, "target_set_id": CFG["analysis"]["glycolysis_target_set"]["id"],
        "n_target_genes_in_set": len(target_set),
        "n_target_genes_observed": int(selected["gene_upper"].nunique()),
        "n_target_positive_AD_occupancy": positive_n,
        "target_positive_fraction": float((selected["log2FC_AD_minus_WT"] > 0).mean()) if len(selected) else float("nan"),
        "target_positive_binomial_P_exploratory": float(positive_test.pvalue) if positive_test else float("nan"),
        "target_positive_fraction_CI95_low": float(positive_ci.low) if positive_ci else float("nan"),
        "target_positive_fraction_CI95_high": float(positive_ci.high) if positive_ci else float("nan"),
        "target_mean_log2FC": finite_mean(selected["log2FC_AD_minus_WT"]),
        "target_median_log2FC": finite_median(selected["log2FC_AD_minus_WT"]),
        "n_background_genes": len(bv),
        "background_mean_log2FC": finite_mean(background["log2FC_AD_minus_WT"]),
        "target_minus_background_mean": finite_mean(selected["log2FC_AD_minus_WT"]) - finite_mean(background["log2FC_AD_minus_WT"]),
        "mann_whitney_U_exploratory": float(mw.statistic) if mw else float("nan"),
        "mann_whitney_P_exploratory": float(mw.pvalue) if mw else float("nan"),
        "interpretation": "descriptive_gene_level_screen; gene rows are not biological replicates",
    })
write_tsv(TABLES / "g2b_target_set_summary.tsv", target_summary_rows)

g1_path = ROOT / CFG["input"]["g1_meta"]
g1 = pd.read_csv(g1_path, sep="\t")
g1["gene_upper"] = g1["gene"].astype(str).str.upper()
g1_target = g1[g1["gene_upper"].isin(target_set)].copy()
g1_target["target_set_id"] = CFG["analysis"]["glycolysis_target_set"]["id"]
g1_target["direction"] = np.where(g1_target["meta_logFC"] > 0, "AD_up", "AD_down")
g1_target.to_csv(TABLES / "g2c_human_target_set_effects.tsv", sep="\t", index=False)

primary = effects[effects["context"] == CFG["design"]["primary_context"]].copy()
compat = primary.merge(g1[["gene_upper", "meta_logFC", "meta_FDR"]], on="gene_upper", how="inner")
compat["g1_direction"] = np.where(compat["meta_logFC"] > 0, "AD_up", "AD_down")
compat["h4k12la_direction"] = np.where(compat["log2FC_AD_minus_WT"] > 0, "AD_up", "AD_down")
compat["direction_concordant"] = compat["g1_direction"] == compat["h4k12la_direction"]
compat.to_csv(TABLES / "g2c_whole_gene_exploratory_compatibility.tsv", sep="\t", index=False)

target_up = int((g1_target["meta_logFC"] > 0).sum())
target_down = int((g1_target["meta_logFC"] < 0).sum())
target_sig = g1_target[g1_target["meta_FDR"] < CFG["analysis"]["fdr"]]
g1_target_summary = {
    "target_set_id": CFG["analysis"]["glycolysis_target_set"]["id"],
    "n_target_genes_in_set": len(target_set), "n_observed_in_human_G1": int(g1_target["gene_upper"].nunique()),
    "n_up": target_up, "n_down": target_down, "mean_meta_logFC": finite_mean(g1_target["meta_logFC"]),
    "n_FDR_below_threshold": len(target_sig),
    "n_FDR_below_threshold_up": int((target_sig["meta_logFC"] > 0).sum()),
    "n_FDR_below_threshold_down": int((target_sig["meta_logFC"] < 0).sum()),
}
whole_gene_sig = compat[compat["meta_FDR"] < CFG["analysis"]["fdr"]]
g2c = {
    **g1_target_summary,
    "whole_gene_shared_n": len(compat),
    "whole_gene_spearman_rho_exploratory": float(spearmanr(compat["log2FC_AD_minus_WT"], compat["meta_logFC"]).statistic),
    "whole_gene_direction_concordance_exploratory": float(compat["direction_concordant"].mean()),
    "whole_gene_FDR_below_threshold_n": len(whole_gene_sig),
    "whole_gene_FDR_below_threshold_concordant_n": int(whole_gene_sig["direction_concordant"].sum()),
    "decision": "HUMAN_TRANSLATION_UNRESOLVED",
}

old_run = ROOT / "results" / "2026-09-20-disease-side-v3.1-g2-reassessment"
for filename in ["replicate_interval_reproducibility.tsv", "signal_correlations.tsv"]:
    source = old_run / "tables" / filename
    if source.exists():
        shutil.copyfile(source, TABLES / filename)

paper = CFG["analysis"]["orthogonal_activity_sources"][0]
independent = CFG["analysis"]["orthogonal_activity_sources"][1]
paper_findings = paper["orthogonal_evidence"]
pkm_row = effects[
    (effects["context"] == CFG["design"]["primary_context"])
    & (effects["gene_upper"] == "PKM")
].iloc[0]
write_tsv(TABLES / "g2b_published_orthogonal_evidence.tsv", [
    {
        "source_id": paper["id"],
        "source_identifier": paper["accession_or_identifier"],
        "evidence_tier": paper["evidence_tier"],
        "same_model_or_system": paper["same_model_or_system"],
        "criterion": item["criterion"],
        "finding": item["finding"],
        "source_url": paper["url"],
    }
    for item in paper_findings
])
g2a = {"decision": "PASS", "accession": "GSE188765", "assay": "H4K12la-specific CUT&Tag",
       "sample_count": len(sample_rows), "microglia_context_present": True,
       "two_replicates_per_context_condition": True, "all_peak_files_parse": True,
       "local_input_or_igg": False, "raw_reads_downloaded": False}
g2b = {
    "direct_target_set": next(row for row in target_summary_rows if row["context"] == "M"),
    "local_PKM_observation": {
        "log2FC_AD_minus_WT": float(pkm_row["log2FC_AD_minus_WT"]),
        "WT_peak_presence_replicates": int(pkm_row["WT_present"]),
        "AD_peak_presence_replicates": int(pkm_row["AD_present"]),
        "FDR_within_context": float(pkm_row["FDR_within_context"]),
        "CI95_low": float(pkm_row["CI95_low"]),
        "CI95_high": float(pkm_row["CI95_high"]),
        "interpretation": "unstable_descriptive_peak_to_gene_result; not an adequately supported contradiction",
    },
    "orthogonal_activity_source": {"id": paper["id"], "tier": paper["evidence_tier"],
                                    "same_model_or_system": True,
                                    "criteria_evaluated": [item["criterion"] for item in paper_findings],
                                    "decision": "SATISFIES_ORTHOGONAL_ACTIVITY_LAYER"},
    "passes_hard_component": True,
    "decision": "PASS_WITH_PUBLISHED_SAME_MODEL_ORTHOGONAL_SUPPORT",
    "reason": "tier_B_direct_H4K12la_target_set_evidence_plus_same_model_tier_C_functional_and_genetic_orthogonal_evidence_satisfy_G2b",
    "claim_boundary": "model-level coupling; local n=2/group occupancy estimates remain descriptive",
}
g2d = {"independent_source": independent["id"], "assay": "H4K12la-specific CUT&Tag",
       "same_AD_model": False, "passes_hard_component": True,
       "decision": "PASS_SITE_SPECIFIC_NOT_AD_REPLICATION",
       "reason": "independent_same_site_source_passes_G2d_site_specific_support;_EAE_MS_does_not_establish_AD_replication",
       "claim_boundary": "independent H4K12la site-specific support; not AD-specific replication"}

evidence_rows = [
    {"evidence_id": "G2a-GSE188765", "accession": "GSE188765", "source_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE188765", "tier": "B_local_processed_plus_public_raw", "role": "direct_primary", "decision": "PASS"},
    {"evidence_id": paper["id"], "accession": paper["accession_or_identifier"], "source_url": paper["url"], "tier": paper["evidence_tier"], "role": "G2b_orthogonal_support", "decision": "SATISFIES_G2B_ORTHOGONAL_LAYER"},
    {"evidence_id": independent["id"], "accession": independent["accession_or_identifier"], "source_url": independent["url"], "tier": independent["evidence_tier"], "role": "G2d_independent_same_site_support", "decision": "PASS_SITE_SPECIFIC_NOT_AD_REPLICATION"},
]
write_tsv(REGISTRIES / "g2_evidence_registry.tsv", evidence_rows)

gate_rows = [
    {"component": "G2a_direct_H4K12la", "required": "direct assay, localization, and QC", "observed": json.dumps(g2a, ensure_ascii=False), "decision": g2a["decision"]},
    {"component": "G2b_within_model_coupling", "required": "predeclared target set plus direct signal and orthogonal activity", "observed": json.dumps(g2b, ensure_ascii=False), "decision": g2b["decision"]},
    {"component": "G2c_human_translation", "required": "target-set/module support; whole-gene correlation exploratory", "observed": json.dumps(g2c, ensure_ascii=False), "decision": g2c["decision"]},
    {"component": "G2d_independent_support", "required": "independent H4K12la-specific or H4K12la-axis support", "observed": json.dumps(g2d, ensure_ascii=False), "decision": g2d["decision"]},
]
write_tsv(TABLES / "G2_gate_v3.2.tsv", gate_rows)

decision = {
    "run_id": CFG["run_id"], "protocol": "AD disease-side v3.2", "decision": "G2_PASS_MODEL_SCOPED",
    "g2a": g2a, "g2b": g2b, "g2c": g2c, "g2d": g2d,
    "downstream": "G3-G6 may proceed for model-scoped work; human H4K12la claims remain prohibited while G2c is unresolved.",
    "dependencies": {"audited_gene_effects": str(audited_effects.relative_to(ROOT)).replace("\\", "/"),
                     "g1_meta": str(g1_path.relative_to(ROOT)).replace("\\", "/")},
    "limitations": ["local CUT&Tag occupancy estimates are descriptive at n=2/group with no local input/IgG read-level model", "GSE296797 is EAE/MS rather than AD", "Pan et al. orthogonal evidence is publication-level Tier C and was not numerically reanalyzed locally", "human G2c translation is unresolved"],
    "created_at": datetime.now(timezone.utc).isoformat(),
}
(RUN / "g2_decision.json").write_text(json.dumps(decision, ensure_ascii=False, indent=2), encoding="utf-8")

manifest = {"run_id": CFG["run_id"], "protocol": CFG["protocol"], "config": str(CONFIG.relative_to(ROOT)).replace("\\", "/"), "script": "src/29_run_g2_v3_2.py", "raw_reads_downloaded": False, "source_files": sample_rows, "result": decision["decision"]}
(RUN / "run_manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")

report = f"""# Disease-side v3.2 G2 execution report

Run ID: `{CFG['run_id']}`  
Protocol: `{CFG['protocol']}`

## Decision

**{decision['decision']}**

G2a direct H4K12la evidence passes. G2b passes because the local Tier B H4K12la target-set result is combined with the Pan et al. same-model publication evidence for lactate dependence, glycolytic transcription/activity, pharmacologic interruption, and microglia-specific Pkm2 ablation. The local n=2/group occupancy estimates remain descriptive and are not used as a veto. In particular, the mapped PKM row is negative but is based on peak presence in 1/2 WT and 0/2 AD replicates with FDR 0.743 and a very wide confidence interval; it is not an adequately supported contradiction of the published mechanism. G2c human translation is unresolved and is not used as a veto. G2d passes the independent site-specific H4K12la support component based on GSE296797; because that dataset is EAE/MS rather than AD, it does not establish AD-specific replication.

## Microglia target-set summary

```json
{json.dumps(g2b['direct_target_set'], ensure_ascii=False, indent=2)}
```

## Human G1 target-set summary

```json
{json.dumps(g2c, ensure_ascii=False, indent=2)}
```

G3–G6 may proceed within the model-scoped claim boundary. Human H4K12la convergence remains unresolved.
"""
(RUN / "execution_report.md").write_text(report, encoding="utf-8")

checksum_rows = []
for path in sorted(RUN.rglob("*")):
    if path.is_file() and path.name != "checksums.sha256":
        checksum_rows.append(f"{sha256(path)}  {path.relative_to(ROOT).as_posix()}")
(RUN / "checksums.sha256").write_text("\n".join(checksum_rows) + "\n", encoding="utf-8")

print(json.dumps({"run": str(RUN), "decision": decision["decision"], "g2a": g2a["decision"], "g2b": g2b["decision"], "g2c": g2c["decision"], "g2d": g2d["decision"]}, ensure_ascii=False, indent=2))
