#!/usr/bin/env python
"""Route G6 candidates to effect-blind validation without a 25% panel veto.

Only donor-count metadata and declared, hashed panel overlap counts are read.
The counts remain subject to independent gene-list review in the SOP.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import platform
import re
import sys
from pathlib import Path

import yaml


CONFIG = "configs/g6_candidate_panel_triage_004.yaml"
HASH_RE = re.compile(r"[0-9a-f]{64}")


def route(size_status: str, missing: list[str], low: list[str], config: dict[str, object]) -> tuple[str, str]:
    """Return candidate next step; legacy overlap never excludes a sized candidate."""
    if size_status == "SIZE_PASS_PENDING_FULL_AUDIT":
        if missing:
            return (str(config["incomplete_coverage_status"]),
                    "OBTAIN_PANEL_GENE_LIST_AND_COMPLETE_BLIND_OVERLAP_AUDIT")
        if low:
            return (str(config["complete_low_coverage_status"]), "RUN_EFFECT_BLIND_MEASUREMENT_VALIDATION")
        return (str(config["complete_adequate_coverage_status"]), "RUN_EFFECT_BLIND_MEASUREMENT_VALIDATION")
    if size_status == "HOLD_DONOR_COUNTS":
        return str(config["donor_hold_status"]), "VERIFY_INDEPENDENT_DONOR_COUNTS"
    return size_status, "AUXILIARY_REVIEW_ONLY"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path, fields: set[str]) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = fields - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{path}: missing columns {sorted(missing)}")
        return list(reader)


def write_once(path: Path, content: str) -> None:
    if path.exists():
        if path.read_text(encoding="utf-8") != content:
            raise FileExistsError(f"Refusing to overwrite different audit output: {path}")
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8", newline="")


def run(root: Path) -> dict[str, object]:
    root = root.resolve()
    config_path = root / CONFIG
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    if config.get("id") != "G6-META-2026-09-27-004":
        raise ValueError("Wrong panel amendment ID")
    for flag in ("below_legacy_coverage_is_hard_exclusion", "full_primary_eligibility_from_triage",
                 "diagnosis_group_effects_allowed"):
        if config.get(flag) is not False:
            raise ValueError(f"{flag} must be false")
    legacy = config["legacy_reference_only"]
    if (int(legacy["G1_genes"]), int(legacy["G1_total"]), float(legacy["G5_module_fraction"])) != (61, 244, 0.25):
        raise ValueError("Legacy comparator must remain 61/244 and 25%, without becoming an exclusion")

    size_path = root / config["size_screen"]
    evidence_path = root / config["panel_overlap_evidence"]
    mapping_path = root / config["G5_orthology_registry"]
    sizes = read_tsv(size_path, {"candidate_id", "size_screen_status", "full_primary_eligibility"})
    evidence = read_tsv(evidence_path, {"candidate_id", "module", "panel_overlap_genes", "source_reference", "source_sha256"})
    mapping = json.loads(mapping_path.read_text(encoding="utf-8"))
    module_sizes = {"G1": 244}
    for entry in mapping["modules"]:
        module_sizes[entry["module"]] = int(entry["n_unique_human_GeneIDs_after_deduplication"])
    if len(module_sizes) != 8:
        raise ValueError("Expected G1 and seven frozen G5 modules")

    candidates = [row["candidate_id"].strip() for row in sizes]
    if len(candidates) != len(set(candidates)) or not all(candidates):
        raise ValueError("Blank or duplicate candidate in donor-count screen")
    if any(row["full_primary_eligibility"] != "NOT_ASSESSED" for row in sizes):
        raise ValueError("Donor-count screen must not assert primary eligibility")

    by_candidate: dict[str, dict[str, int]] = {candidate: {} for candidate in candidates}
    for row in evidence:
        candidate, module = row["candidate_id"].strip(), row["module"].strip()
        if candidate not in by_candidate or module not in module_sizes:
            raise ValueError(f"Unknown candidate/module: {candidate}/{module}")
        if module in by_candidate[candidate]:
            raise ValueError(f"Duplicate panel overlap: {candidate}/{module}")
        try:
            number = int(row["panel_overlap_genes"])
        except ValueError as exc:
            raise ValueError(f"Invalid overlap count: {candidate}/{module}") from exc
        if not 0 <= number <= module_sizes[module]:
            raise ValueError(f"Overlap count outside frozen module size: {candidate}/{module}")
        source = row["source_reference"].strip()
        stated_hash = row["source_sha256"].strip().lower()
        if not source or not HASH_RE.fullmatch(stated_hash):
            raise ValueError(f"Panel source and SHA-256 required: {candidate}/{module}")
        source_path = (root / source).resolve()
        if not source_path.is_relative_to(root) or not source_path.is_file():
            raise ValueError(f"Panel source must be a local file inside the project: {source}")
        if sha256(source_path) != stated_hash:
            raise ValueError(f"Panel source hash mismatch: {source}")
        by_candidate[candidate][module] = number

    rows: list[dict[str, object]] = []
    for row in sizes:
        candidate = row["candidate_id"].strip()
        counts = by_candidate[candidate]
        missing = sorted(set(module_sizes) - set(counts))
        low = []
        if "G1" in counts and counts["G1"] < 61:
            low.append("G1")
        for module, total in module_sizes.items():
            if module != "G1" and module in counts and counts[module] < math.ceil(total * 0.25):
                low.append(module)
        size_status = row["size_screen_status"]
        status, next_step = route(size_status, missing, low, config)
        rows.append({
            "candidate_id": candidate,
            "size_screen_status": size_status,
            "G1_panel_overlap_genes": counts.get("G1", "UNKNOWN"),
            "G1_legacy_reference_genes": 61,
            "modules_with_overlap_below_legacy_25pct": ";".join(low) or "NONE_KNOWN",
            "modules_without_overlap_evidence": ";".join(missing) or "NONE",
            "panel_triage_status": status,
            "next_step": next_step,
            "full_primary_eligibility": "NOT_ASSESSED",
        })

    columns = list(rows[0]) if rows else []
    table = "\t".join(columns) + "\n" + "".join("\t".join(str(row[k]) for k in columns) + "\n" for row in rows)
    status_counts = {status: sum(row["panel_triage_status"] == status for row in rows)
                     for status in sorted({row["panel_triage_status"] for row in rows})}
    audit = {
        "protocol_id": config["id"],
        "scope": "donor_count_and_panel_overlap_triage_only",
        "diagnosis_group_effects_computed": False,
        "full_primary_eligibility_assessed": False,
        "configuration_sha256": sha256(config_path),
        "size_screen_sha256": sha256(size_path),
        "panel_overlap_evidence_sha256": sha256(evidence_path),
        "G5_orthology_registry_sha256": sha256(mapping_path),
        "candidate_count": len(rows),
        "status_counts": status_counts,
        "python_version": platform.python_version(),
        "pyyaml_version": yaml.__version__,
    }
    write_once(root / config["output_table"], table)
    write_once(root / config["output_audit"], json.dumps(audit, ensure_ascii=False, indent=2, sort_keys=True) + "\n")
    return audit


def main() -> int:
    if len(sys.argv) != 2:
        raise SystemExit("Usage: 99_screen_g6_candidate_panel_004.py <project_root>")
    audit = run(Path(sys.argv[1]))
    print(json.dumps(audit["status_counts"], ensure_ascii=False, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
