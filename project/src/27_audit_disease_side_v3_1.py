"""Independent post-run audit for the v3.1 disease-side package."""

from __future__ import annotations

import csv
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RUN = ROOT / "results" / "2026-09-20-disease-side-v3.1"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


decision = json.loads((RUN / "final_decision.json").read_text(encoding="utf-8"))
estimability = read_tsv(RUN / "gates" / "primary_cohort_region_estimability.tsv")
meta = read_tsv(RUN / "meta" / "random_effects_meta_summary.tsv")
loo = read_tsv(RUN / "meta" / "leave_one_cohort_out.tsv")
effects = read_tsv(RUN / "effects" / "primary_region_gene_effects.tsv")
d7_audit = json.loads((ROOT / "data" / "interim" / "D7_GSE148822" / "D7_mtx_integrity_audit.json").read_text(encoding="utf-8"))

hits = [row for row in meta if float(row["meta_FDR"]) < 0.05]
strict_hits = [row for row in hits if row["k"] == "5"]
stable_by_gene = {}
for row in loo:
    if row["gene"] in {x["gene"] for x in strict_hits}:
        stable_by_gene.setdefault(row["gene"], []).append(row["direction_stable"] == "TRUE")
strict_loo_stable = all(len(v) == 3 and all(v) for v in stable_by_gene.values()) and len(stable_by_gene) == len(strict_hits)

required = [
    RUN / "final_decision.json", RUN / "execution_report.md", RUN / "run_manifest.json",
    RUN / "registries" / "dataset_registry.tsv", RUN / "registries" / "sample_manifest.tsv",
    RUN / "registries" / "file_manifest.tsv", RUN / "registries" / "evidence_registry.tsv",
    RUN / "registries" / "exclusion_log.tsv", RUN / "registries" / "per_dataset_provenance.json",
    RUN / "registries" / "checksums.sha256", RUN / "gates" / "G2_assessment.tsv",
]
checks = {
    "decision_is_stopped_at_g2": decision.get("decision") == "STOPPED_AT_G2",
    "all_five_primary_region_units_estimable": len(estimability) == 5 and all(row["estimable"] == "TRUE" for row in estimability),
    "d7_all_mtx_files_valid": len(d7_audit) == 30 and all(row["ok"] for row in d7_audit),
    "primary_effects_only_predeclared_cohorts": set(row["cohort"] for row in effects) == {"D3_GSE160936", "D7_GSE148822", "D8_GSE174367"},
    "random_effects_meta_hits_413": len(hits) == 413,
    "strict_all_region_hits_244": len(strict_hits) == 244,
    "strict_hits_loo_stable": strict_loo_stable,
    "required_registry_files_present": all(path.exists() for path in required),
    "historical_v3_archive_present": (ROOT / "results" / "archive" / "2026-09-20-pre-v3.1").exists(),
}
audit = {
    "run_id": decision.get("run_id"),
    "checks": checks,
    "passed": all(checks.values()),
    "observed": {"meta_fdr_hits": len(hits), "strict_all_region_hits": len(strict_hits), "loo_genes_checked": len(stable_by_gene), "primary_effect_rows": len(effects)},
}
(RUN / "audit.json").write_text(json.dumps(audit, ensure_ascii=False, indent=2), encoding="utf-8")
print(json.dumps(audit, ensure_ascii=False, indent=2))
if not audit["passed"]:
    raise SystemExit(2)
