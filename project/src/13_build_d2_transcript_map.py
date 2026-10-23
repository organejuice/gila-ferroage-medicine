#!/usr/bin/env python3
"""Extract a reproducible Ensembl transcript-to-gene map for D2."""
from pathlib import Path
import argparse
import csv
import gzip
import json
import re

ATTR_RE = re.compile(r'(?P<key>transcript_id|gene_id|gene_name) "(?P<value>[^"]+)"')

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    args = p.parse_args()
    root = args.root.resolve()
    gtf = root / "data" / "raw" / "D2_GSE227221" / "Homo_sapiens.GRCh38.109.gtf.gz"
    out = root / "data" / "interim" / "D2_transcript_gene_map.tsv"
    rows = {}
    with gzip.open(gtf, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "transcript":
                continue
            attrs = dict((m.group("key"), m.group("value")) for m in ATTR_RE.finditer(fields[8]))
            tx = attrs.get("transcript_id", "").split(".", 1)[0]
            gene = attrs.get("gene_name", "")
            gene_id = attrs.get("gene_id", "").split(".", 1)[0]
            if tx and gene:
                rows.setdefault(tx, {"transcript_id": tx, "gene_id": gene_id, "gene_symbol": gene})
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["transcript_id", "gene_id", "gene_symbol"], delimiter="\t")
        w.writeheader()
        w.writerows(sorted(rows.values(), key=lambda x: x["transcript_id"]))
    audit = {"gtf": str(gtf.relative_to(root)), "release": "Ensembl release 109 / GRCh38", "mapped_transcripts": len(rows), "output": str(out.relative_to(root))}
    (root / "data" / "interim" / "D2_transcript_map_provenance.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
    print(json.dumps(audit, indent=2))

if __name__ == "__main__":
    main()