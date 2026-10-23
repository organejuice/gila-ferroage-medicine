#!/usr/bin/env python
"""Prospective, metadata-only donor-count screen for G6-META-2026-09-27-002.

This does not assess full input eligibility or inspect any expression data.
"""

from __future__ import annotations

import csv
import hashlib
import json
import platform
import sys
from pathlib import Path

import yaml


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def count(value: str, field: str, candidate: str) -> int | None:
    if value.strip().upper() in {"", "UNKNOWN", "NA"}:
        return None
    try:
        result = int(value)
    except ValueError as exc:
        raise ValueError(f"{candidate}: invalid {field}={value!r}") from exc
    if result < 0:
        raise ValueError(f"{candidate}: negative {field}={result}")
    return result


def keep_or_write(path: Path, content: str) -> None:
    if path.exists():
        if path.read_text(encoding="utf-8") != content:
            raise FileExistsError(f"Refusing to replace an existing different result: {path}")
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8", newline="")


def run(root: Path, config_relative: str) -> int:
    root = root.resolve()
    config_path = root / config_relative
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    expected_minima = {
        "G6-META-2026-09-27-002": (8, 8),
        "G6-META-2026-09-27-003": (4, 4),
    }
    if config.get("id") not in expected_minima:
        raise ValueError("Wrong or missing amendment ID")
    if config.get("diagnosis_group_effects_allowed") is not False:
        raise ValueError("Effect analysis must remain prohibited")
    if config.get("full_primary_eligibility_from_size_alone") is not False:
        raise ValueError("A donor-count screen cannot assert full primary eligibility")
    minima = config["minimum_independent_donors"]
    min_ad, min_control = int(minima["AD"]), int(minima["Control"])
    if (min_ad, min_control) != expected_minima[config["id"]]:
        raise ValueError("Donor minima disagree with the numbered amendment")

    registry_path = root / config["input_registry"]
    with registry_path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"candidate_id", "AD_donors", "control_donors", "metadata_screen_status", "primary_source"}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError(f"Missing registry columns: {sorted(required - set(reader.fieldnames or []))}")
        source_rows = list(reader)
    if not source_rows:
        raise ValueError("Empty candidate registry")

    seen: set[str] = set()
    output: list[dict[str, object]] = []
    for row in source_rows:
        candidate = row["candidate_id"].strip()
        if not candidate or candidate in seen:
            raise ValueError(f"Blank or repeated candidate_id: {candidate!r}")
        seen.add(candidate)
        ad = count(row["AD_donors"], "AD_donors", candidate)
        control = count(row["control_donors"], "control_donors", candidate)
        if ad is None or control is None:
            status = config["unknown_count_status"]
        elif ad >= min_ad and control >= min_control:
            status = config["size_pass_status"]
        else:
            status = config["size_fail_status"]
        output.append({
            "candidate_id": candidate,
            "AD_donors_source_level": row["AD_donors"],
            "control_donors_source_level": row["control_donors"],
            "minimum_AD": min_ad,
            "minimum_Control": min_control,
            "size_screen_status": status,
            "full_primary_eligibility": "NOT_ASSESSED",
            "prior_metadata_screen_status": row["metadata_screen_status"],
            "primary_source": row["primary_source"],
        })

    if config.get("rank_by_donor_count"):
        for row in output:
            row["donor_count_priority_rank"] = "NOT_RANKED"
        passing = [row for row in output if row["size_screen_status"] == config["size_pass_status"]]
        passing.sort(key=lambda row: (
            -min(int(row["AD_donors_source_level"]), int(row["control_donors_source_level"])),
            -(int(row["AD_donors_source_level"]) + int(row["control_donors_source_level"])),
            str(row["candidate_id"]),
        ))
        for rank, row in enumerate(passing, start=1):
            row["donor_count_priority_rank"] = rank

    fields = list(output[0])
    lines = ["\t".join(fields)]
    lines.extend("\t".join(str(row[field]) for field in fields) for row in output)
    table_text = "\n".join(lines) + "\n"
    audit = {
        "protocol_id": config["id"],
        "scope": "source_metadata_donor_count_only",
        "diagnosis_group_effects_computed": False,
        "full_primary_eligibility_assessed": False,
        "configuration": str(config_path.relative_to(root)),
        "configuration_sha256": digest(config_path),
        "input_registry": str(registry_path.relative_to(root)),
        "input_registry_sha256": digest(registry_path),
        "candidate_count": len(output),
        "status_counts": {status: sum(row["size_screen_status"] == status for row in output)
                          for status in sorted({row["size_screen_status"] for row in output})},
        "python_version": platform.python_version(),
        "pyyaml_version": yaml.__version__,
    }
    audit_text = json.dumps(audit, ensure_ascii=False, indent=2, sort_keys=True) + "\n"
    keep_or_write(root / config["output_table"], table_text)
    keep_or_write(root / config["output_audit"], audit_text)
    print(json.dumps(audit["status_counts"], ensure_ascii=False, sort_keys=True))
    return 0


def main() -> int:
    if len(sys.argv) != 2:
        raise SystemExit("Usage: 97_screen_g6_candidate_size_002.py <project_root>")
    return run(Path(sys.argv[1]), "configs/g6_candidate_size_screen_002.yaml")


if __name__ == "__main__":
    raise SystemExit(main())
