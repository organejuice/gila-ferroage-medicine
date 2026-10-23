#!/usr/bin/env python3
"""Audit complete GSE160936 Cell Ranger archives without extracting them."""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import re
import tarfile
from pathlib import Path


def count_gzip_lines(archive: tarfile.TarFile, member: tarfile.TarInfo) -> int:
    stream = archive.extractfile(member)
    if stream is None:
        raise ValueError(f"Could not read archive member {member.name}")
    with gzip.GzipFile(fileobj=stream) as handle:
        return sum(1 for _ in handle)


def matrix_shape(archive: tarfile.TarFile, member: tarfile.TarInfo) -> tuple[int, int, int]:
    stream = archive.extractfile(member)
    if stream is None:
        raise ValueError(f"Could not read archive member {member.name}")
    with gzip.GzipFile(fileobj=stream) as handle:
        for raw in handle:
            line = raw.decode("utf-8", errors="replace").strip()
            if not line or line.startswith("%"):
                continue
            fields = line.split()
            if len(fields) != 3:
                raise ValueError(f"Unexpected Matrix Market dimensions in {member.name}: {line}")
            return tuple(int(value) for value in fields)
    raise ValueError(f"Matrix Market header is missing in {member.name}")


def audit(root: Path) -> dict:
    manifest_path = root / "results" / "tables" / "D3_download_manifest.tsv"
    raw_dir = root / "data" / "raw" / "D3_GSE160936"
    rows = []
    with manifest_path.open("r", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["status"] != "complete":
                continue
            path = root / row["file"]
            with tarfile.open(path, mode="r:gz") as archive:
                members = archive.getmembers()
                feature = next((m for m in members if m.name.endswith("/features.tsv.gz")), None)
                barcode = next((m for m in members if m.name.endswith("/barcodes.tsv.gz")), None)
                matrix = next((m for m in members if m.name.endswith("/matrix.mtx.gz")), None)
                if not all((feature, barcode, matrix)):
                    raise ValueError(f"Archive missing Cell Ranger members: {path.name}")
                n_features = count_gzip_lines(archive, feature)
                n_barcodes = count_gzip_lines(archive, barcode)
                matrix_rows, matrix_cols, nonzero_entries = matrix_shape(archive, matrix)
            match = re.match(r"(GSM\d+)_", path.name)
            rows.append({
                "geo_sample": match.group(1) if match else path.stem,
                "archive": str(path.relative_to(root)),
                "feature_rows": n_features,
                "barcode_rows": n_barcodes,
                "matrix_rows": matrix_rows,
                "matrix_columns": matrix_cols,
                "matrix_nonzero_entries": nonzero_entries,
                "feature_barcode_dimensions_match": matrix_rows == n_features and matrix_cols == n_barcodes,
            })
    out = root / "results" / "tables" / "D3_archive_audit.tsv"
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t") if rows else None
        if writer:
            writer.writeheader()
            writer.writerows(rows)
    result = {"complete_archives_audited": len(rows), "all_dimensions_match": all(row["feature_barcode_dimensions_match"] for row in rows), "table": str(out.relative_to(root))}
    (root / "data" / "interim" / "d3_archive_audit.json").write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    args = parser.parse_args()
    print(json.dumps(audit(args.root.resolve()), ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()