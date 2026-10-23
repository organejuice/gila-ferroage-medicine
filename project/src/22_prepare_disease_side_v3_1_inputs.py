"""Prepare and audit the processed public inputs for disease-side v3.1.

This script only handles auditable metadata, selective GEO downloads, and file
integrity checks. It never edits raw files and never downloads FASTQ/SRA data.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import re
from pathlib import Path
from typing import Dict, Iterable, List


ROOT = Path(__file__).resolve().parents[1]
D7_RAW = ROOT / "data" / "raw" / "D7_GSE148822"
D8_RAW = ROOT / "data" / "raw" / "D8_GSE174367"
D7_INTERIM = ROOT / "data" / "interim" / "D7_GSE148822"
D8_INTERIM = ROOT / "data" / "interim" / "D8_GSE174367"


def parse_soft(path: Path) -> List[Dict[str, str]]:
    """Parse sample-level fields and supplementary URLs from a GEO SOFT file."""
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        text = handle.read()
    records: List[Dict[str, str]] = []
    for block in text.split("^SAMPLE = ")[1:]:
        lines = block.splitlines()
        if not lines:
            continue
        rec: Dict[str, str] = {"accession": lines[0].strip()}
        for line in lines[1:]:
            if not line.startswith("!Sample_"):
                continue
            key, sep, value = line.partition(" = ")
            if not sep:
                continue
            key = key.removeprefix("!Sample_")
            if key.startswith("characteristics_ch1 = "):
                continue
            rec[key] = value.strip()
            if key.startswith("supplementary_file_"):
                rec[f"{key}_url"] = value.strip()
        # GEO stores repeated characteristics as independent lines.
        for line in lines:
            prefix = "!Sample_characteristics_ch1 = "
            if not line.startswith(prefix):
                continue
            key, sep, value = line[len(prefix) :].partition(": ")
            if sep:
                rec[key.strip().lower().replace(" ", "_")] = value.strip()
        records.append(rec)
    return records


def donor_from_title(title: str) -> str:
    match = re.search(r"Donor(\d+)", title, flags=re.IGNORECASE)
    if not match:
        raise ValueError(f"Cannot parse donor from GEO sample title: {title}")
    return f"Donor{int(match.group(1))}"


def d7_region(value: str) -> str:
    value = value.lower()
    if "occipital" in value and "temporal" not in value:
        return "OC"
    if "occipitotemporal" in value:
        return "OTC"
    raise ValueError(f"Unexpected D7 region: {value}")


def d7_primary_decision(rec: Dict[str, str]) -> tuple[bool, str]:
    group = rec.get("sample_group", "")
    donor_number = int(re.search(r"\d+", rec["donor"]).group())
    if group == "AD":
        return True, "primary_AD"
    if group == "CTR" and donor_number not in {6, 7, 8}:
        return True, "primary_strict_control"
    if group == "CTR+":
        return False, "excluded_CTR_plus"
    return False, f"excluded_group_{group or 'missing'}"


def prepare_d7_manifest() -> Path:
    soft = D7_RAW / "GSE148822_family.soft.gz"
    if not soft.exists():
        raise FileNotFoundError(soft)
    rows: List[Dict[str, str]] = []
    for source in parse_soft(soft):
        matrix_url = source.get("supplementary_file_3_url", "")
        if not matrix_url or matrix_url == "NONE":
            continue
        title = source.get("title", "")
        if not title.startswith("SN_"):
            continue
        rec: Dict[str, str] = {
            "cohort": "D7_GSE148822",
            "accession": source["accession"],
            "title": title,
            "donor": donor_from_title(title),
            "region": d7_region(source.get("region", "")),
            "age": source.get("age", ""),
            "sex": source.get("gender", ""),
            "RIN": source.get("rin", ""),
            "diagnosis": source.get("clinical_diagnosis", ""),
            "sample_group": source.get("sample_group", ""),
            "matrix_url": matrix_url.replace("ftp://", "https://"),
            "barcode_url": source.get("supplementary_file_1_url", "").replace("ftp://", "https://"),
            "feature_url": source.get("supplementary_file_2_url", "").replace("ftp://", "https://"),
        }
        include, reason = d7_primary_decision(rec)
        rec["include_primary"] = str(include).lower()
        rec["decision"] = reason
        filename = matrix_url.rsplit("/", 1)[-1]
        rec["matrix_file"] = filename
        rec["barcode_file"] = rec["barcode_url"].rsplit("/", 1)[-1]
        rec["feature_file"] = rec["feature_url"].rsplit("/", 1)[-1]
        rec["local_matrix"] = str(D7_RAW / filename)
        rec["local_barcode"] = str(D7_RAW / rec["barcode_file"])
        rec["local_feature"] = str(D7_RAW / rec["feature_file"])
        rows.append(rec)
    rows.sort(key=lambda row: (row["donor"], row["region"]))
    D7_INTERIM.mkdir(parents=True, exist_ok=True)
    out = D7_INTERIM / "D7_sample_manifest.tsv"
    with out.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=sorted(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    included = [row for row in rows if row["include_primary"] == "true"]
    counts = {}
    for row in included:
        key = (row["region"], row["sample_group"])
        counts[key] = counts.get(key, 0) + 1
    expected = {("OC", "AD"): 10, ("OC", "CTR"): 5, ("OTC", "AD"): 10, ("OTC", "CTR"): 5}
    if counts != expected:
        raise RuntimeError(f"D7 strict primary counts do not match freeze: {counts} != {expected}")
    print(json.dumps({"manifest": str(out), "n_sn_samples": len(rows), "n_primary_samples": len(included), "counts": {f"{k[0]}_{k[1]}": v for k, v in sorted(counts.items())}}, ensure_ascii=False, indent=2))
    return out


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def audit_paths(paths: Iterable[Path]) -> List[Dict[str, object]]:
    result = []
    for path in paths:
        exists = path.exists()
        result.append({
            "file": str(path),
            "exists": exists,
            "size_bytes": path.stat().st_size if exists else 0,
            "sha256": sha256(path) if exists else "",
        })
    return result


def audit_d7() -> None:
    manifest = D7_INTERIM / "D7_sample_manifest.tsv"
    if not manifest.exists():
        prepare_d7_manifest()
    with manifest.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    paths: List[Path] = [D7_RAW / "GSE148822_family.soft.gz"]
    for row in rows:
        if row["include_primary"] == "true":
            paths.extend(Path(row[key]) for key in ("local_matrix", "local_barcode", "local_feature"))
    audit = audit_paths(paths)
    out = D7_INTERIM / "D7_file_audit.json"
    out.write_text(json.dumps(audit, ensure_ascii=False, indent=2), encoding="utf-8")
    missing = [item["file"] for item in audit if not item["exists"]]
    print(json.dumps({"audit": str(out), "n_files": len(audit), "missing": missing}, ensure_ascii=False, indent=2))
    if missing:
        raise SystemExit(2)


def audit_d8() -> None:
    paths = [
        D8_RAW / "GSE174367_snRNA-seq_cell_meta.csv.gz",
        D8_RAW / "GSE174367_snRNA-seq_filtered_feature_bc_matrix.h5",
    ]
    audit = audit_paths(paths)
    out = D8_INTERIM / "D8_file_audit.json"
    D8_INTERIM.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(audit, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps({"audit": str(out), "files": audit}, ensure_ascii=False, indent=2))
    if any(not item["exists"] or item["size_bytes"] < 10000 for item in audit):
        raise SystemExit(2)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("command", choices=["prepare-d7", "audit-d7", "audit-d8"])
    args = parser.parse_args()
    if args.command == "prepare-d7":
        prepare_d7_manifest()
    elif args.command == "audit-d7":
        audit_d7()
    else:
        audit_d8()


if __name__ == "__main__":
    main()
