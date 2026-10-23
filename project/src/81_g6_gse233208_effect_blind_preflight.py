#!/usr/bin/env python3
"""Validate existing GSE233208 author-correction and spatial input audit artifacts.

This gate is metadata/input-only. It never fits a model or computes group effects.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path
from typing import Any


PROTOCOL_ID = "G6-PROT-2026-09-26-003"
AUDIT_ROOT = Path("results/2026-09-24-disease-side-v3.5-g6-input-audit")
OUT_REL = Path("results/2026-09-26-g6-gse233208-visium-input-audit")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def tsv_rows(path: Path) -> list[dict[str, str]]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def kv_audit(path: Path) -> dict[str, str]:
    result: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if "=" in line:
            key, value = line.split("=", 1)
            result[key.strip()] = value.strip()
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path("."))
    args = parser.parse_args()
    root = args.root.resolve()
    audit = root / AUDIT_ROOT
    output = root / OUT_REL
    tables = audit / "tables"
    required = {
        "sample_manifest": audit / "registries/GSE233208_donor_section_manifest.tsv",
        "image_manifest": tables / "GSE233208_spatial_image_coordinate_manifest.tsv.gz",
        "count_audit": tables / "GSE233208_expression_layer_audit.tsv.gz",
        "correction_audit": tables / "GSE233208_metadata_correction_audit.txt",
        "g1_registry": audit / "registries/G1_microglia_signed_sender_244.tsv.gz",
        "g5_registry": audit / "registries/G5_receiver_one_to_one_human_ortholog_registry.tsv.gz",
        "overlap_audit": audit / "GSE174367_GSE233208_overlap_audit.json",
        "reference_gate": root / "configs/g6_v3.6_reference_admissibility.json",
        "raw_object": root / "data/raw/GSE233208/GSE233208_Human_visium_ADDS_seurat_processed.rds.gz",
        "author_table": root / "data/raw/G6_source_metadata_2026-09-23/GSE233208_author_Supplementary_Table1.tsv.gz",
    }
    missing = [name for name, path in required.items() if not path.is_file()]
    if missing:
        raise SystemExit("Missing required audited source/artifact(s): " + ", ".join(missing))

    donors = tsv_rows(required["sample_manifest"])
    images = tsv_rows(required["image_manifest"])
    count_rows = tsv_rows(required["count_audit"])
    audit_kv = kv_audit(required["correction_audit"])
    overlap = json.loads(required["overlap_audit"].read_text(encoding="utf-8"))
    g1 = tsv_rows(required["g1_registry"])
    g5 = tsv_rows(required["g5_registry"])

    eligible = [row for row in donors if row.get("role_status", "").startswith("validation_eligible_")]
    ad = [row for row in eligible if row.get("Diagnosis") == "AD"]
    control = [row for row in eligible if row.get("Diagnosis") == "Control"]
    donor_ids = [row.get("donor_id", "") for row in eligible]
    sample_ids = [row.get("Seq.ID", "") for row in eligible]
    image_by_sample = {row.get("Sample", ""): row for row in images}
    errors: list[str] = []
    checks: dict[str, Any] = {}

    def check(name: str, passed: bool, detail: Any) -> None:
        checks[name] = {"pass": bool(passed), "detail": detail}
        if not passed:
            errors.append(name)

    check("unique_target_donors", len(eligible) == 20 and all(donor_ids) and len(set(donor_ids)) == 20,
          {"n": len(eligible), "unique": len(set(donor_ids))})
    check("frozen_diagnosis_counts", len(ad) == 10 and len(control) == 10,
          {"AD": len(ad), "Control": len(control)})
    check("one_sample_per_target_donor", all(sample_ids) and len(set(sample_ids)) == 20,
          {"n_samples": len(sample_ids), "unique": len(set(sample_ids))})
    check("all_target_images_mapped", all(s in image_by_sample for s in sample_ids),
          {"target_samples": len(sample_ids), "mapped": sum(s in image_by_sample for s in sample_ids)})
    coordinate_ok = all(
        int(image_by_sample[s].get("n_coordinates_unmatched_to_metadata", "-1") or -1) == 0
        and int(image_by_sample[s].get("n_coordinates_matched_to_metadata", "0") or 0) > 0
        and int(image_by_sample[s].get("n_author_samples", "0") or 0) == 1
        for s in sample_ids if s in image_by_sample
    )
    check("target_coordinate_mapping", coordinate_ok, "each target image maps uniquely to metadata spots")
    check("author_sample_spot_counts", audit_kv.get("all_author_Sample_n_spot_counts_match") == "TRUE",
          audit_kv.get("all_author_Sample_n_spot_counts_match"))
    check("integer_spatial_counts", audit_kv.get("expression_noninteger_entries") == "0",
          audit_kv.get("expression_noninteger_entries"))
    check("count_layer_dimensions", audit_kv.get("expression_counts_layer") == "dgCMatrix 36601x115451",
          audit_kv.get("expression_counts_layer"))
    check("frozen_signature_registry", len(g1) == 244 and len(g5) > 0,
          {"G1_rows": len(g1), "G5_rows": len(g5)})
    check("GSE174367_reference_held", overlap.get("candidate_pairs_matching_age_sex_PMI_diagnosis") == 15,
          {"matched_attribute_candidates": overlap.get("candidate_pairs_matching_age_sex_PMI_diagnosis"),
           "target_match_counts": overlap.get("candidate_GSE233208_Diagnosis_counts"),
           "identity_crosswalk": "not available"})

    reference_lock = root / "configs/g6_v3.6_reference_admissibility.json"
    reference_admissible = False
    if reference_lock.is_file():
        try:
            reference_data = json.loads(reference_lock.read_text(encoding="utf-8"))
            donor_gate = reference_data.get("donor_independence", {})
            sha_ok = lambda value: isinstance(value, str) and len(value) == 64 and all(c in "0123456789abcdefABCDEF" for c in value)
            hashed_paths = [
                (reference_data.get("reference_path"), reference_data.get("reference_sha256")),
                (reference_data.get("donor_manifest_path"), reference_data.get("donor_manifest_sha256")),
                (donor_gate.get("identity_crosswalk_path"), donor_gate.get("identity_crosswalk_sha256")),
            ]
            files_match_hashes = all(
                isinstance(path, str) and sha_ok(expected)
                and (root / path).is_file() and sha256(root / path).lower() == expected.lower()
                for path, expected in hashed_paths
            )
            reference_admissible = (
                reference_data.get("status") == "PASS"
                and reference_data.get("effect_blind") is True
                and reference_data.get("group_effects_calculated") is False
                and sha_ok(reference_data.get("reference_sha256"))
                and sha_ok(reference_data.get("donor_manifest_sha256"))
                and donor_gate.get("status") == "PASS"
                and sha_ok(donor_gate.get("identity_crosswalk_sha256"))
                and donor_gate.get("eligible_spatial_donors_checked") == 20
                and donor_gate.get("unresolved_matches") == 0
                and files_match_hashes
            )
        except (json.JSONDecodeError, OSError):
            reference_admissible = False
    blockers = [] if reference_admissible else ["independent_reference_and_donor_crosswalk_pending"]
    metadata_pass = not errors
    status = "READY_FOR_BLIND_REFERENCE_CALIBRATION" if metadata_pass and reference_admissible else "INPUT_METADATA_PASS_REFERENCE_HOLD"

    output.mkdir(parents=True, exist_ok=True)
    manifest = []
    for name, path in required.items():
        manifest.append({"input": str(path.relative_to(root)), "bytes": path.stat().st_size, "sha256": sha256(path)})
    result = {
        "protocol_id": PROTOCOL_ID,
        "status": status,
        "metadata_and_source_checks_pass": metadata_pass,
        "independent_reference_pass": reference_admissible,
        "group_effects_calculated": False,
        "donor_level_effects_calculated": False,
        "checks": checks,
        "errors": errors,
        "blockers": blockers,
        "input_manifest": manifest,
    }
    (output / "GSE233208_effect_blind_input_preflight.json").write_text(
        json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    report = [
        "# GSE233208 effect-blind input preflight",
        "",
        f"- Protocol: `{PROTOCOL_ID}`",
        f"- Status: `{status}`",
        f"- Metadata/source checks passed: `{metadata_pass}`",
        f"- Independent reference passed: `{reference_admissible}`",
        "- Group effects calculated: `false`",
        "",
        "## Checks",
        "",
    ]
    report.extend(f"- {'PASS' if row['pass'] else 'FAIL'} `{name}`: `{row['detail']}`" for name, row in checks.items())
    report.extend(["", "## Blockers", ""])
    report.extend(f"- {item}" for item in (blockers + errors) or ["none"])
    (output / "GSE233208_effect_blind_input_preflight.md").write_text("\n".join(report) + "\n", encoding="utf-8")
    print(f"{status}; group_effects_calculated=false; blockers={','.join(blockers + errors) or 'none'}")
    return 0 if not errors else 2


if __name__ == "__main__":
    raise SystemExit(main())
