"""Audit and fingerprint the frozen G1 region-aware reanalysis outputs."""

from __future__ import annotations

import csv
import hashlib
import json
import math
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def finite(value: str) -> bool:
    try:
        return math.isfinite(float(value))
    except (TypeError, ValueError):
        return False


def git_commit(root: Path) -> str | None:
    try:
        result = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=root,
            check=True,
            capture_output=True,
            text=True,
        )
    except (OSError, subprocess.CalledProcessError):
        return None
    return result.stdout.strip() or None


def main() -> int:
    root = Path(sys.argv[1] if len(sys.argv) > 1 else ".").resolve()
    run_id = sys.argv[2] if len(sys.argv) > 2 else "2026-09-20-disease-side-v3-g1-region-aware"
    out = root / "results" / run_id
    effects = out / "effects" / "D3_microglia_region_gene_effects.tsv"
    feasibility = out / "effects" / "D6_microglia_region_feasibility.tsv"
    meta = out / "meta" / "random_effects_meta_summary.tsv"
    loco = out / "meta" / "leave_one_cohort_out.tsv"
    sensitivity = out / "meta" / "sensitivity_audit.tsv"
    criteria = out / "gates" / "G1_reanalysis_criteria.tsv"
    final_path = out / "final_decision.json"

    required = [effects, feasibility, meta, loco, sensitivity, criteria, final_path]
    errors: list[str] = []
    missing = [str(p.relative_to(root)) for p in required if not p.exists()]
    if missing:
        errors.append("missing required outputs: " + ", ".join(missing))

    effect_rows = read_tsv(effects) if effects.exists() else []
    feasibility_rows = read_tsv(feasibility) if feasibility.exists() else []
    meta_rows = read_tsv(meta) if meta.exists() else []
    loco_rows = read_tsv(loco) if loco.exists() else []
    criteria_rows = read_tsv(criteria) if criteria.exists() else []
    final = json.loads(final_path.read_text(encoding="utf-8")) if final_path.exists() else {}

    regions = sorted({row.get("region", "") for row in effect_rows})
    if len(effect_rows) != 33454:
        errors.append(f"D3 effect row count is {len(effect_rows)}, expected 33454")
    if regions != ["EC", "SSC"]:
        errors.append(f"D3 regions are {regions}, expected ['EC', 'SSC']")
    bad_se = [row.get("gene", "") for row in effect_rows if not finite(row.get("SE", ""))]
    if bad_se:
        errors.append(f"D3 has {len(bad_se)} non-finite SE values")
    if len(feasibility_rows) != 3:
        errors.append(f"D6 feasibility row count is {len(feasibility_rows)}, expected 3")
    if sum(row.get("estimable_for_frozen_G1") == "TRUE" for row in feasibility_rows) != 0:
        errors.append("D6 contains an estimable frozen-G1 region")
    if not meta_rows or meta_rows[0].get("random_effects_meta_status") != "NOT_ESTIMABLE":
        errors.append("random-effects meta status is not NOT_ESTIMABLE")
    if not loco_rows or loco_rows[0].get("status") != "NOT_RUN":
        errors.append("leave-one-cohort-out status is not NOT_RUN")
    if final.get("decision") != "STOPPED_AT_G1":
        errors.append(f"final decision is {final.get('decision')!r}, expected STOPPED_AT_G1")
    if final.get("d6_estimable_regions") != 0:
        errors.append("final decision reports a non-zero D6 estimable region count")
    criteria_map = {row.get("criterion"): row.get("status") for row in criteria_rows}
    for key in ("brain_regions_modeled_separately", "random_effects_meta_BH_FDR_lt_0_05", "leave_one_cohort_out_stable"):
        if criteria_map.get(key) != "FAIL":
            errors.append(f"G1 criterion {key} is not FAIL")

    input_candidates = [
        root / "params.yaml",
        root / "configs" / "disease_side_v3.0_frozen.yaml",
        root / "data" / "interim" / "D3_pseudobulk_counts_microglia.tsv.gz",
        root / "data" / "interim" / "D3_pseudobulk_metadata_microglia.tsv",
        root / "data" / "raw" / "D6_GSE214979" / "GSE214979_cell_metadata.csv.gz",
        root / "results" / "archive" / "2026-09-20-pre-v3-freeze" / "results" / "tables" / "D6_signature_concordance.tsv",
    ]
    inputs = {
        str(path.relative_to(root)): {"sha256": sha256(path), "bytes": path.stat().st_size}
        for path in input_candidates
        if path.exists()
    }
    output_files = [
        path for path in sorted(out.rglob("*"))
        if path.is_file() and path.name not in {"run_manifest.json", "audit.json"}
    ]
    manifest = {
        "run_id": run_id,
        "protocol": "AD disease-side v3.0",
        "status": "PASS" if not errors else "FAIL",
        "created_at": datetime.now(timezone.utc).isoformat(),
        "git_commit": git_commit(root),
        "inputs": inputs,
        "outputs": {
            str(path.relative_to(root)): {"sha256": sha256(path), "bytes": path.stat().st_size}
            for path in output_files
        },
        "checks": {
            "d3_effect_rows": len(effect_rows),
            "d3_regions": regions,
            "d3_nonfinite_se": len(bad_se),
            "d6_feasibility_rows": len(feasibility_rows),
            "d6_estimable_regions": sum(row.get("estimable_for_frozen_G1") == "TRUE" for row in feasibility_rows),
            "meta_status": meta_rows[0].get("random_effects_meta_status") if meta_rows else None,
            "loco_status": loco_rows[0].get("status") if loco_rows else None,
            "final_decision": final.get("decision"),
        },
        "errors": errors,
    }
    out.mkdir(parents=True, exist_ok=True)
    (out / "run_manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    (out / "audit.json").write_text(json.dumps({"status": manifest["status"], "errors": errors}, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": manifest["status"], "errors": errors, "output": str(out)}, ensure_ascii=False))
    return 0 if not errors else 1


if __name__ == "__main__":
    raise SystemExit(main())
