"""Validate compressed GEO Matrix Market files before pseudobulk use."""

from __future__ import annotations

import csv
import gzip
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
manifest_path = ROOT / "data" / "interim" / "D7_GSE148822" / "D7_sample_manifest.tsv"
rows = list(csv.DictReader(manifest_path.open(encoding="utf-8"), delimiter="\t"))
rows = [row for row in rows if row["include_primary"] == "true"]
results = []
for row in rows:
    path = Path(row["local_matrix"])
    expected = None
    actual = 0
    error = ""
    nrow = ncol = None
    try:
        with gzip.open(path, "rt", encoding="utf-8", errors="strict") as handle:
            for line in handle:
                line = line.strip()
                if not line or line.startswith("%"):
                    continue
                if expected is None:
                    dims = line.split()
                    if len(dims) != 3:
                        raise ValueError(f"bad Matrix Market dimension line: {line}")
                    nrow, ncol, expected = int(dims[0]), int(dims[1]), int(dims[2])
                else:
                    actual += 1
    except Exception as exc:
        error = repr(exc)
    results.append({
        "accession": row["accession"],
        "donor": row["donor"],
        "region": row["region"],
        "file": str(path),
        "size_bytes": path.stat().st_size if path.exists() else 0,
        "expected_nnz": expected,
        "observed_nnz": actual,
        "nrow": nrow,
        "ncol": ncol,
        "ok": (not error) and expected == actual,
        "error": error,
    })

out = ROOT / "data" / "interim" / "D7_GSE148822" / "D7_mtx_integrity_audit.json"
out.write_text(json.dumps(results, ensure_ascii=False, indent=2), encoding="utf-8")
print(json.dumps({"audit": str(out), "n_files": len(results), "invalid": [x for x in results if not x["ok"]]}, ensure_ascii=False, indent=2))
if any(not x["ok"] for x in results):
    raise SystemExit(2)
