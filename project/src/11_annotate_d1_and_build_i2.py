#!/usr/bin/env python3
"""Annotate D1 H4K12la peaks and build the descriptive I2 regulatory profile."""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import re
from bisect import bisect_left, bisect_right
from collections import defaultdict
from pathlib import Path

import pandas as pd


PEAK_RE = re.compile(r"(?P<gsm>GSM\d+)_(?P<context>[A-Za-z]+)_(?P<condition>WT|AD)-(?P<replicate>\d+)_peaks\.txt\.gz$")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def load_genes(gtf_path: Path):
    by_chrom = defaultdict(dict)
    with gzip.open(gtf_path, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 13:
                continue
            chrom, name, strand = fields[2], fields[3], fields[4]
            try:
                tx_start, tx_end = int(fields[5]), int(fields[6])
            except ValueError:
                continue
            symbol = fields[12].strip()
            if not symbol or symbol == "none":
                continue
            key = symbol.upper()
            if key not in by_chrom[chrom]:
                by_chrom[chrom][key] = {"symbol": symbol, "start": tx_start, "end": tx_end, "tss": []}
            rec = by_chrom[chrom][key]
            rec["start"] = min(rec["start"], tx_start)
            rec["end"] = max(rec["end"], tx_end)
            rec["tss"].append(tx_start if strand == "+" else tx_end - 1)
    index = {}
    for chrom, records in by_chrom.items():
        rows = list(records.values())
        rows.sort(key=lambda r: (r["start"], r["end"], r["symbol"]))
        for r in rows:
            r["tss_center"] = sorted(r["tss"])[len(r["tss"]) // 2]
        index[chrom] = {"rows": rows, "starts": [r["start"] for r in rows], "tss_rows": sorted(rows, key=lambda r: r["tss_center"]), "tss": sorted(r["tss_center"] for r in rows)}
    return index


def annotate_peak(chrom: str, start: int, end: int, index, max_tss_distance: int):
    if chrom not in index:
        return None
    block = index[chrom]
    center = (start + end - 1) // 2
    right = bisect_right(block["starts"], end)
    overlaps = []
    j = right - 1
    while j >= 0 and block["rows"][j]["end"] > start:
        r = block["rows"][j]
        if r["start"] < end and r["end"] > start:
            overlaps.append(r)
        j -= 1
    if overlaps:
        r = min(overlaps, key=lambda x: (min(abs(center - t) for t in x["tss"]), x["symbol"]))
        return r["symbol"], "gene_body_overlap", min(abs(center - t) for t in r["tss"])
    pos = bisect_left(block["tss"], center)
    candidates = []
    for k in range(max(0, pos - 3), min(len(block["tss_rows"]), pos + 4)):
        r = block["tss_rows"][k]
        d = abs(center - r["tss_center"])
        if d <= max_tss_distance:
            candidates.append((d, r["symbol"], r))
    if not candidates:
        return None
    d, symbol, _ = min(candidates, key=lambda x: (x[0], x[1]))
    return symbol, "nearest_tss", d


def read_peaks(root: Path, annotation_index, max_tss_distance: int):
    rows = []
    sample_gene_counts = defaultdict(int)
    audit = []
    for path in sorted((root / "data" / "raw" / "D1_GSE188765").glob("*_peaks.txt.gz")):
        m = PEAK_RE.match(path.name)
        if not m:
            continue
        labels = m.groupdict()
        seen = set()
        mapped = 0
        for line in gzip.open(path, "rt", encoding="utf-8", errors="replace"):
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 4:
                continue
            chrom, start_s, end_s, score_s = fields[:4]
            try:
                start, end = int(start_s), int(end_s)
                score = float(score_s)
            except ValueError:
                continue
            key = (chrom, start, end)
            if key in seen:
                continue
            seen.add(key)
            ann = annotate_peak(chrom, start, end, annotation_index, max_tss_distance)
            if ann is None:
                gene_symbol, method, distance = "", "unmapped", ""
            else:
                gene_symbol, method, distance = ann
                mapped += 1
                sample_gene_counts[(labels["context"], labels["condition"], labels["gsm"], gene_symbol)] += 1
            rows.append({**labels, "chrom": chrom, "start": start, "end": end, "peak_score": score, "gene_symbol": gene_symbol, "annotation_method": method, "distance_to_tss": distance})
        audit.append({**labels, "unique_peaks": len(seen), "mapped_peaks": mapped, "mapping_rate": mapped / len(seen) if seen else 0.0})
    return rows, sample_gene_counts, audit


def load_orthology(path: Path):
    classes = defaultdict(lambda: {"mouse": set(), "human": set()})
    with path.open("r", encoding="utf-8", errors="replace") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            key = row.get("DB Class Key", "")
            taxon = row.get("NCBI Taxon ID", "")
            symbol = row.get("Symbol", "").strip()
            if not key or not symbol:
                continue
            if taxon == "10090":
                classes[key]["mouse"].add(symbol)
            elif taxon == "9606":
                classes[key]["human"].add(symbol)
    mapping = defaultdict(set)
    for item in classes.values():
        if len(item["mouse"]) == 1 and len(item["human"]) == 1:
            mapping[next(iter(item["mouse"]))].add(next(iter(item["human"])))
    return mapping


def write_outputs(root: Path, rows, sample_gene_counts, audit, orthology):
    tables = root / "results" / "tables"
    tables.mkdir(parents=True, exist_ok=True)
    peak_path = tables / "D1_peak_to_gene.tsv"
    with peak_path.open("w", encoding="utf-8", newline="") as f:
        fields = ["gsm", "context", "condition", "replicate", "chrom", "start", "end", "peak_score", "gene_symbol", "annotation_method", "distance_to_tss"]
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    profile_rows = []
    for context in sorted({r["context"] for r in rows}):
        samples = sorted({r["gsm"] for r in rows if r["context"] == context})
        labels = {r["gsm"]: (r["condition"], r["replicate"]) for r in rows if r["context"] == context}
        genes = sorted({r["gene_symbol"] for r in rows if r["context"] == context and r["gene_symbol"]})
        for gene in genes:
            wt_counts = [sample_gene_counts[(context, "WT", s, gene)] for s in samples if labels[s][0] == "WT"]
            ad_counts = [sample_gene_counts[(context, "AD", s, gene)] for s in samples if labels[s][0] == "AD"]
            wt_present = sum(x > 0 for x in wt_counts) / len(wt_counts)
            ad_present = sum(x > 0 for x in ad_counts) / len(ad_counts)
            delta_presence = ad_present - wt_present
            delta_count = (sum(ad_counts) / len(ad_counts)) - (sum(wt_counts) / len(wt_counts))
            profile_rows.append({"context": context, "mouse_gene_symbol": gene, "wt_present_fraction": wt_present, "ad_present_fraction": ad_present, "delta_presence": delta_presence, "wt_mean_peak_count": sum(wt_counts) / len(wt_counts), "ad_mean_peak_count": sum(ad_counts) / len(ad_counts), "delta_mean_peak_count": delta_count, "direction": "AD_gain" if delta_presence > 0 else ("AD_loss" if delta_presence < 0 else "unchanged"), "human_one_to_one_symbols": ";".join(sorted(orthology.get(gene, set()))), "one_to_one_ortholog": len(orthology.get(gene, set())) == 1})
    profile = pd.DataFrame(profile_rows)
    profile.to_csv(tables / "D1_gene_regulatory_profile.tsv", sep="\t", index=False)
    profile[profile["one_to_one_ortholog"] & (profile["context"] == "M")].to_csv(tables / "I2_regulatory_activity.tsv", sep="\t", index=False)
    primary = profile[(profile["context"] == "M") & profile["one_to_one_ortholog"]].copy()
    primary["module"] = "M2"
    primary["species_mapping"] = "mouse_to_human_one_to_one_MGI"
    primary["signed_activity"] = primary["delta_presence"]
    primary[["module", "mouse_gene_symbol", "human_one_to_one_symbols", "delta_presence", "delta_mean_peak_count", "direction", "species_mapping", "signed_activity"]].to_csv(tables / "M2_H4K12la_regulatory_reference.tsv", sep="\t", index=False)
    pd.DataFrame(audit).to_csv(tables / "D1_annotation_audit.tsv", sep="\t", index=False)
    counts = {"peak_rows": len(rows), "mapped_peak_rows": sum(bool(r["gene_symbol"]) for r in rows), "mapped_peak_rate": sum(bool(r["gene_symbol"]) for r in rows) / len(rows), "primary_M_mouse_genes": int((profile["context"] == "M").sum()), "primary_M_one_to_one_human_genes": int(len(primary)), "primary_M_AD_gain_genes": int((primary["direction"] == "AD_gain").sum()), "primary_M_AD_loss_genes": int((primary["direction"] == "AD_loss").sum())}
    provenance = {"annotation": "UCSC mm10 ncbiRefSeq table", "annotation_version": "2021-04-23", "annotation_file": "data/raw/D1_GSE188765/mm10_ncbiRefSeq_2021-04-23.txt.gz", "annotation_sha256": sha256(root / "data" / "raw" / "D1_GSE188765" / "mm10_ncbiRefSeq_2021-04-23.txt.gz"), "orthology": "MGI HOM_MouseHumanSequence.rpt", "orthology_file": "data/raw/D1_GSE188765/MGI_HOM_MouseHumanSequence.rpt", "orthology_sha256": sha256(root / "data" / "raw" / "D1_GSE188765" / "MGI_HOM_MouseHumanSequence.rpt"), "peak_to_gene_rule": "gene-body overlap first; otherwise nearest TSS within 100000 bp; one primary gene per peak", "primary_context": "M (microglia); T (total brain) retained as support", "results": counts}
    (root / "data" / "interim" / "D1_annotation_provenance.json").write_text(json.dumps(provenance, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(counts, ensure_ascii=False, indent=2))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--max-tss-distance", type=int, default=100000)
    args = parser.parse_args()
    root = args.root.resolve()
    index = load_genes(root / "data" / "raw" / "D1_GSE188765" / "mm10_ncbiRefSeq_2021-04-23.txt.gz")
    rows, counts, audit = read_peaks(root, index, args.max_tss_distance)
    orthology = load_orthology(root / "data" / "raw" / "D1_GSE188765" / "MGI_HOM_MouseHumanSequence.rpt")
    write_outputs(root, rows, counts, audit, orthology)


if __name__ == "__main__":
    main()