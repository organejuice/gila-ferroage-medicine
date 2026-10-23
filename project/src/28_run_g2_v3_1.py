"""Reassess the frozen v3.1 G2 gate from direct GSE188765 peak files.

The input is public processed H4K12la CUT&Tag peak data. This script deliberately
does not claim read-level differential binding because the local input contains
peak tables, not BAM/FASTQ and not an input/IgG control. It reports sample-level
peak QC, pairwise interval reproducibility, gene-mapped peak-score effects,
glycolysis compatibility, and compatibility with the existing G1 meta signature.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import math
import re
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from scipy.stats import spearmanr, t as student_t


ROOT = Path(__file__).resolve().parents[1]
CONFIG = ROOT / "configs" / "g2_v3.1_reassessment.yaml"
with CONFIG.open(encoding="utf-8") as handle:
    CFG = yaml.safe_load(handle)

RUN = ROOT / "results" / CFG["run_id"]
TABLES = RUN / "tables"
REG = RUN / "registries"
TABLES.mkdir(parents=True, exist_ok=True)
REG.mkdir(parents=True, exist_ok=True)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_tsv(path: Path, rows: list[dict], fields: list[str] | None = None) -> None:
    if fields is None:
        fields = list(rows[0]) if rows else []
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def parse_sample(path: Path) -> dict:
    match = re.match(r"GSM\d+_([TM])_(WT|AD)-(\d+)_peaks\.txt\.gz$", path.name)
    if not match:
        raise ValueError(f"Unexpected peak file name: {path.name}")
    context, condition, replicate = match.groups()
    return {"gsm": path.name.split("_")[0], "context": context, "condition": condition,
            "replicate": int(replicate), "path": path}


def read_peaks(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path, sep="\t", header=None, compression="gzip",
                        names=["chrom", "start", "end", "score", "strand"])
    frame = frame[["chrom", "start", "end", "score"]].copy()
    for col in ["start", "end", "score"]:
        frame[col] = pd.to_numeric(frame[col], errors="coerce")
    frame = frame.dropna(subset=["chrom", "start", "end", "score"])
    frame["start"] = frame["start"].astype(np.int64)
    frame["end"] = frame["end"].astype(np.int64)
    frame["score"] = frame["score"].astype(float)
    frame = frame[frame["end"] > frame["start"]].drop_duplicates(["chrom", "start", "end"])
    return frame.sort_values(["chrom", "start", "end"], kind="mergesort").reset_index(drop=True)


def overlap_summary(a: pd.DataFrame, b: pd.DataFrame) -> tuple[int, int, int, float, float]:
    """Return A count, B count, A-overlap fraction, B-overlap fraction, interval Jaccard.

    Any-base overlap is intentionally conservative as a reproducibility screen;
    it is not an IDR estimate.
    """
    by_chrom_a = {chrom: group[["start", "end"]].to_numpy() for chrom, group in a.groupby("chrom")}
    by_chrom_b = {chrom: group[["start", "end"]].to_numpy() for chrom, group in b.groupby("chrom")}
    a_hit = 0
    b_hit = 0
    for chrom, intervals in by_chrom_a.items():
        other = by_chrom_b.get(chrom, np.empty((0, 2), dtype=np.int64))
        j = 0
        for start, end in intervals:
            while j < len(other) and other[j, 1] <= start:
                j += 1
            k = j
            hit = False
            while k < len(other) and other[k, 0] < end:
                if other[k, 1] > start:
                    hit = True
                    break
                k += 1
            a_hit += int(hit)
    for chrom, intervals in by_chrom_b.items():
        other = by_chrom_a.get(chrom, np.empty((0, 2), dtype=np.int64))
        j = 0
        for start, end in intervals:
            while j < len(other) and other[j, 1] <= start:
                j += 1
            k = j
            hit = False
            while k < len(other) and other[k, 0] < end:
                if other[k, 1] > start:
                    hit = True
                    break
                k += 1
            b_hit += int(hit)
    # Jaccard here uses peak records linked by any overlap, not bp-union Jaccard.
    union = len(a) + len(b) - min(a_hit, b_hit)
    jaccard = min(a_hit, b_hit) / union if union else float("nan")
    return len(a), len(b), a_hit, b_hit, jaccard


def welch_effect(values_a: np.ndarray, values_b: np.ndarray, pseudocount: float) -> dict:
    # values are normalized peak scores; log2 fold-change is AD minus WT.
    x = np.log2(values_a + pseudocount)
    y = np.log2(values_b + pseudocount)
    delta = float(np.mean(y) - np.mean(x))
    vx = float(np.var(x, ddof=1)) if len(x) > 1 else 0.0
    vy = float(np.var(y, ddof=1)) if len(y) > 1 else 0.0
    se = math.sqrt(vx / len(x) + vy / len(y)) if len(x) and len(y) else float("nan")
    denom = (vx / len(x) + vy / len(y)) ** 2
    df = denom / (((vx / len(x)) ** 2 / max(1, len(x) - 1)) + ((vy / len(y)) ** 2 / max(1, len(y) - 1))) if denom else float("nan")
    t_stat = delta / se if se else (float("inf") if delta > 0 else float("-inf") if delta < 0 else 0.0)
    p = float(2 * student_t.sf(abs(t_stat), df)) if math.isfinite(df) and df > 0 else float("nan")
    critical = float(student_t.ppf(0.975, df)) if math.isfinite(df) and df > 0 else float("nan")
    return {"log2FC_AD_minus_WT": delta, "SE": se, "df": df, "t": t_stat, "P.Value": p,
            "CI95_low": delta - critical * se if math.isfinite(critical) else float("nan"),
            "CI95_high": delta + critical * se if math.isfinite(critical) else float("nan")}


samples = [parse_sample(path) for path in sorted((ROOT / CFG["input"]["peak_dir"]).glob("GSM*_peaks.txt.gz"))]
if len(samples) != 8:
    raise RuntimeError(f"Expected 8 GSE188765 peak files, found {len(samples)}")

peak_frames: dict[str, pd.DataFrame] = {}
sample_rows = []
for sample in samples:
    frame = read_peaks(sample["path"])
    peak_frames[sample["gsm"]] = frame
    sample_rows.append({
        "gsm": sample["gsm"], "context": sample["context"], "condition": sample["condition"],
        "replicate": sample["replicate"], "file": sample["path"].relative_to(ROOT).as_posix(),
        "file_size": sample["path"].stat().st_size, "sha256": sha256(sample["path"]),
        "peak_rows": len(frame), "unique_intervals": frame.drop_duplicates(["chrom", "start", "end"]).shape[0],
        "total_peak_score": float(frame["score"].sum()),
    })
write_tsv(TABLES / "sample_qc.tsv", sample_rows)

replicate_rows = []
for context in ["M", "T"]:
    for condition in ["WT", "AD"]:
        current = [s for s in samples if s["context"] == context and s["condition"] == condition]
        if len(current) != 2:
            raise RuntimeError(f"Expected two replicates for {context}/{condition}")
        first, second = current
        a_count, b_count, a_hit, b_hit, jaccard = overlap_summary(peak_frames[first["gsm"]], peak_frames[second["gsm"]])
        replicate_rows.append({
            "context": context, "condition": condition,
            "gsm_1": first["gsm"], "gsm_2": second["gsm"], "n_peaks_1": a_count,
            "n_peaks_2": b_count, "n_peaks_1_overlapping_2": a_hit,
            "n_peaks_2_overlapping_1": b_hit, "fraction_1_overlapping_2": a_hit / a_count,
            "fraction_2_overlapping_1": b_hit / b_count, "record_jaccard": jaccard,
            "interpretation": "any-base interval overlap; not IDR",
        })
write_tsv(TABLES / "replicate_interval_reproducibility.tsv", replicate_rows)

# Use the previously audited one-peak-to-one-gene mapping as an explicit, read-only
# processed dependency. Only microglia are primary; total brain is support.
mapping_path = ROOT / CFG["input"]["peak_to_gene"]
mapping = pd.read_csv(mapping_path, sep="\t", dtype={"gsm": str, "gene_symbol": str})
mapping = mapping[(mapping["context"].isin(["M", "T"])) & mapping["gene_symbol"].notna() & (mapping["gene_symbol"] != "")].copy()
mapping["gene_symbol"] = mapping["gene_symbol"].str.strip()
mapping["peak_score"] = pd.to_numeric(mapping["peak_score"], errors="coerce").fillna(0.0)

sample_total_score = {row["gsm"]: row["total_peak_score"] for row in sample_rows}
gene_signal = mapping.groupby(["context", "gsm", "gene_symbol"], as_index=False)["peak_score"].sum()
gene_signal["normalized_signal"] = gene_signal.apply(
    lambda row: row["peak_score"] / sample_total_score.get(row["gsm"], np.nan) * 1e6, axis=1
)
wide = gene_signal.pivot_table(index=["context", "gene_symbol"], columns="gsm", values="normalized_signal", fill_value=0.0)
wide = wide.reset_index()

effect_rows = []
for _, row in wide.iterrows():
    context, gene = row["context"], row["gene_symbol"]
    wt = [row[s["gsm"]] for s in samples if s["context"] == context and s["condition"] == "WT"]
    ad = [row[s["gsm"]] for s in samples if s["context"] == context and s["condition"] == "AD"]
    stats = welch_effect(np.asarray(wt, dtype=float), np.asarray(ad, dtype=float), CFG["analysis"]["pseudocount"])
    effect_rows.append({"context": context, "gene_symbol": gene,
                        "WT_rep1_signal": wt[0], "WT_rep2_signal": wt[1],
                        "AD_rep1_signal": ad[0], "AD_rep2_signal": ad[1],
                        "WT_present": int(sum(v > 0 for v in wt)), "AD_present": int(sum(v > 0 for v in ad)),
                        **stats})
effects = pd.DataFrame(effect_rows)
effects["P.Value"] = effects["P.Value"].fillna(1.0)
effects = effects.sort_values(["context", "P.Value", "gene_symbol"])

def benjamini_hochberg(values: pd.Series) -> pd.Series:
    array = values.to_numpy(dtype=float)
    order = np.argsort(array, kind="mergesort")
    ranked = array[order] * len(array) / np.arange(1, len(array) + 1)
    adjusted = np.minimum.accumulate(ranked[::-1])[::-1].clip(0.0, 1.0)
    output = np.empty(len(array), dtype=float)
    output[order] = adjusted
    return pd.Series(output, index=values.index)

effects["FDR_within_context"] = effects.groupby("context", group_keys=False)["P.Value"].apply(benjamini_hochberg)
effects.to_csv(TABLES / "gene_mapped_h4k12la_effects.tsv", sep="\t", index=False)

glycolysis = {gene.upper() for gene in CFG["analysis"]["glycolysis_genes"]}
glyco_rows = effects[effects["gene_symbol"].str.upper().isin(glycolysis)].copy()
glyco_rows["is_positive_AD_occupancy"] = glyco_rows["log2FC_AD_minus_WT"] > 0
glyco_rows.to_csv(TABLES / "glycolysis_h4k12la_effects.tsv", sep="\t", index=False)

# Gene-wide replicate signal correlations and context concordance.
correlation_rows = []
for context in ["M", "T"]:
    selected = [s for s in samples if s["context"] == context]
    context_wide = wide[wide["context"] == context].copy()
    for condition in ["WT", "AD"]:
        names = [s["gsm"] for s in selected if s["condition"] == condition]
        rho, p = spearmanr(context_wide[names[0]], context_wide[names[1]])
        correlation_rows.append({"comparison": f"{context}_{condition}_replicate_signal", "n_genes": len(context_wide), "spearman_rho": rho, "P.Value": p})
m = effects[effects["context"] == "M"].set_index("gene_symbol")
t = effects[effects["context"] == "T"].set_index("gene_symbol")
common = m.index.intersection(t.index)
rho, p = spearmanr(m.loc[common, "log2FC_AD_minus_WT"], t.loc[common, "log2FC_AD_minus_WT"])
correlation_rows.append({"comparison": "M_vs_T_AD_minus_WT_effect", "n_genes": len(common), "spearman_rho": rho, "P.Value": p})
write_tsv(TABLES / "signal_correlations.tsv", correlation_rows)

# Compatibility with the already computed human G1 meta direction, restricted to
# genes shared with the mouse microglial H4K12la map and the glycolysis set.
g1_path = ROOT / CFG["input"]["g1_meta"]
g1 = pd.read_csv(g1_path, sep="\t")
g1["gene_upper"] = g1["gene"].astype(str).str.upper()
primary = effects[effects["context"] == "M"].copy()
primary["gene_upper"] = primary["gene_symbol"].str.upper()
compat = primary.merge(g1[["gene_upper", "meta_logFC", "meta_FDR"]], on="gene_upper", how="inner")
compat["gene_symbol"] = compat["gene_symbol"]
compat["g1_direction"] = np.where(compat["meta_logFC"] > 0, "AD_up", "AD_down")
compat["h4k12la_direction"] = np.where(compat["log2FC_AD_minus_WT"] > 0, "AD_up", "AD_down")
compat["direction_concordant"] = compat["g1_direction"] == compat["h4k12la_direction"]
compat.to_csv(TABLES / "g1_h4k12la_direction_compatibility.tsv", sep="\t", index=False)

def finite_mean(series: pd.Series) -> float:
    values = pd.to_numeric(series, errors="coerce").replace([np.inf, -np.inf], np.nan).dropna()
    return float(values.mean()) if len(values) else float("nan")

microglia_glyco = glyco_rows[glyco_rows["context"] == "M"].copy()
total_brain_glyco = glyco_rows[glyco_rows["context"] == "T"].copy()

glyco_summary = {
    "n_genes_in_set": len(glycolysis),
    "microglia_n_genes_observed": int(microglia_glyco["gene_symbol"].str.upper().nunique()),
    "microglia_n_positive_AD_occupancy": int(microglia_glyco["is_positive_AD_occupancy"].sum()),
    "microglia_positive_fraction": float(microglia_glyco["is_positive_AD_occupancy"].mean()) if len(microglia_glyco) else float("nan"),
    "microglia_mean_log2FC_AD_minus_WT": finite_mean(microglia_glyco["log2FC_AD_minus_WT"]),
    "total_brain_n_genes_observed": int(total_brain_glyco["gene_symbol"].str.upper().nunique()),
    "total_brain_n_positive_AD_occupancy": int(total_brain_glyco["is_positive_AD_occupancy"].sum()),
    "total_brain_positive_fraction": float(total_brain_glyco["is_positive_AD_occupancy"].mean()) if len(total_brain_glyco) else float("nan"),
    "total_brain_mean_log2FC_AD_minus_WT": finite_mean(total_brain_glyco["log2FC_AD_minus_WT"]),
}
compat_summary = {
    "n_shared_mouse_human_genes": len(compat),
    "n_direction_concordant": int(compat["direction_concordant"].sum()) if len(compat) else 0,
    "direction_concordance_fraction": float(compat["direction_concordant"].mean()) if len(compat) else float("nan"),
    "spearman_h4k12la_vs_g1": float(spearmanr(compat["log2FC_AD_minus_WT"], compat["meta_logFC"]).statistic) if len(compat) > 2 else float("nan"),
}
fdr_compat = compat[compat["meta_FDR"] < CFG["analysis"]["fdr"]].copy()
compat_summary.update({
    "n_g1_fdr_below_threshold": len(fdr_compat),
    "n_g1_fdr_below_threshold_direction_concordant": int(fdr_compat["direction_concordant"].sum()) if len(fdr_compat) else 0,
    "g1_fdr_below_threshold_direction_concordance_fraction": float(fdr_compat["direction_concordant"].mean()) if len(fdr_compat) else float("nan"),
})

direct_assay = {
    "accession": "GSE188765", "assay": "H4K12la-specific CUT&Tag", "tier": "A/B",
    "sample_count": len(samples), "all_files_parse": True, "microglia_context_present": True,
    "two_replicates_per_context_condition": True, "raw_reads_available": True,
    "local_peak_files": 8, "input_or_igg_control_in_local_peak_package": False,
    "read_level_reprocessing_performed": False,
}

gate_rows = [
    {"criterion": "direct_H4K12la_assay", "required": "H4K12la-specific ChIP/CUT&Tag/CUT&RUN/quantitative proteomics", "observed": "GSE188765 H4K12la CUT&Tag", "decision": "PASS"},
    {"criterion": "microglial_localization", "required": "cell-specific measurement or explicit source localization", "observed": "purified microglia CUT&Tag plus total brain support", "decision": "PASS"},
    {"criterion": "sample_integrity", "required": "8 expected files, two replicates per context/condition", "observed": f"{len(samples)} parsed; all peak rows unique", "decision": "PASS"},
    {"criterion": "replicate_reproducibility", "required": "auditable and adequate for quantitative inference", "observed": "interval overlap and gene-signal correlations reported; n=2/group", "decision": "LIMITED"},
    {"criterion": "AD_WT_microglial_effect", "required": "auditable differential occupancy with uncertainty", "observed": "descriptive gene-level log2FC with Welch CI; no input/IgG and n=2/group", "decision": "LIMITED"},
    {"criterion": "glycolysis_direction_compatibility", "required": "direct H4K12la direction compatible with glycolysis/lactate", "observed": json.dumps(glyco_summary, ensure_ascii=False), "decision": "REPORTED_NOT_HARD_PASS"},
    {"criterion": "G1_direction_compatibility", "required": "compatible with human G1 signature", "observed": json.dumps(compat_summary, ensure_ascii=False), "decision": "REPORTED_NOT_HARD_PASS"},
    {"criterion": "independent_concordant_support", "required": "independent source", "observed": "H3K18la AD-model paper and other Kla sources are different mark/assay or Tier C", "decision": "SUPPORTIVE_NOT_INDEPENDENT_H4K12LA_REPLICATION"},
]
write_tsv(TABLES / "G2_gate_reassessment.tsv", gate_rows)

decision = {
    "run_id": CFG["run_id"], "protocol": "AD disease-side v3.1", "source_dataset": "GSE188765",
    "decision": "G2_REASSESSMENT_COMPLETE_NOT_PASS",
    "direct_assay": direct_assay, "glycolysis_summary": glyco_summary,
    "g1_compatibility_summary": compat_summary,
    "limitations": [
        "Only processed peak tables and peak scores were used locally; no BAM/FASTQ was reprocessed.",
        "No local input or IgG control was available for a full read-level differential binding model.",
        "Only two biological replicates per group/context limit uncertainty estimation and independent replication.",
        "The direct H4K12la source is an AD-like mouse model; it is not human AD brain.",
        "Broader H3K18la, H3K9la, pan-Kla and non-histone lactylation evidence was not substituted into the frozen H4K12la gate.",
    ],
    "downstream": "G3-G6 remain locked pending a valid G2 hard-gate decision.",
    "created_at": datetime.now(timezone.utc).isoformat(),
}
(RUN / "g2_decision.json").write_text(json.dumps(decision, ensure_ascii=False, indent=2), encoding="utf-8")

manifest = {
    "run_id": CFG["run_id"], "config": str(CONFIG.relative_to(ROOT)).replace("\\", "/"),
    "script": "src/28_run_g2_v3_1.py", "environment": {"python": "D:/python/python.exe", "packages": ["pandas", "numpy", "scipy", "PyYAML"]},
    "dependencies": {"peak_to_gene": str(mapping_path.relative_to(ROOT)).replace("\\", "/"), "g1_meta": str(g1_path.relative_to(ROOT)).replace("\\", "/")},
    "source_files": sample_rows, "result": decision["decision"],
}
(RUN / "run_manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")

report = f"""# G2 v3.1 reassessment execution report

Run ID: `{CFG['run_id']}`  
Source: GSE188765, H4K12la-specific CUT&Tag  
Protocol: `{CFG['protocol']}`

## Decision

**{decision['decision']}**

The direct-assay criterion passes: eight expected processed peak files were parsed, with microglia and total-brain contexts and two WT/AD replicates per context. The overall G2 hard gate does not pass in this reassessment because the local package supports descriptive occupancy effects only (n=2/group, no local input/IgG controls), independent AD H4K12la replication is not established, and the cross-species G1 direction correlation is near zero.

## Primary microglia results

```json
{json.dumps(glyco_summary, ensure_ascii=False, indent=2)}
```

## Human G1 compatibility

```json
{json.dumps(compat_summary, ensure_ascii=False, indent=2)}
```

## Scope and limitations

- Microglia (M) is the primary context; total brain (T) is reported only as supporting context.
- Gene-level effects use mapped peak-score sums normalized per million total peak score and Welch-style descriptive intervals.
- No BAM/FASTQ or input/IgG read-level differential-binding reanalysis was performed locally.
- The source is an AD-like mouse model, not human AD brain.
- H3K18la, H3K9la, pan-Kla, and non-histone lactylation sources were not substituted into the frozen H4K12la gate.

G3–G6 remain locked pending a valid G2 hard-gate decision.
"""
(RUN / "execution_report.md").write_text(report, encoding="utf-8")

checksum_rows = []
for path in sorted(RUN.rglob("*")):
    if path.is_file() and path.name != "checksums.sha256":
        checksum_rows.append(f"{sha256(path)}  {path.relative_to(ROOT).as_posix()}")
(RUN / "checksums.sha256").write_text("\n".join(checksum_rows) + "\n", encoding="utf-8")
print(json.dumps({"run": str(RUN), "decision": decision["decision"], "glycolysis": glyco_summary, "compatibility": compat_summary}, ensure_ascii=False, indent=2))
