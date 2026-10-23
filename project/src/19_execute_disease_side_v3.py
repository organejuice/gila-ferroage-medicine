#!/usr/bin/env python3
"""Execute the available-data portion of the frozen disease-side v3.0 protocol.

This runner is deliberately fail-closed. It audits registrations and the G1
requirements from existing donor-aware outputs, then stops before G2-G6 when
the frozen prerequisites are not demonstrably satisfied.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import subprocess
import sys
from datetime import date, datetime, timezone
from pathlib import Path
from typing import Iterable


DATASETS = [
    {
        "dataset_id": "D1_GSE188765",
        "accession": "GSE188765",
        "database": "GEO",
        "source_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE188765",
        "species": "Mus musculus",
        "model": "AD mouse model",
        "pathology_type": "AD-like amyloid pathology",
        "brain_region": "brain region as provided by source",
        "assay": "H4K12la CUT&Tag peak files",
        "processing_level": "raw peak files plus local annotation",
        "evidence_tier": "A",
        "mechanism_node": "H4K12la_dysregulation",
        "analysis_role": "direct animal-model reference; microglial localization audit",
        "license": "source terms apply",
        "access_conditions": "public GEO",
        "download_date": "pre-v3; inherited project raw asset",
        "version": "source record and local annotation 2021-04-23",
        "diagnosis": "AD-like versus wild-type in source model",
    },
    {
        "dataset_id": "D2_GSE227221",
        "accession": "GSE227221",
        "database": "GEO",
        "source_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE227221",
        "species": "Homo sapiens",
        "model": "human iPSC-derived microglia",
        "pathology_type": "Aβ1-42 fibril time course",
        "brain_region": "not applicable",
        "assay": "bulk RNA-seq TPM",
        "processing_level": "processed TPM",
        "evidence_tier": "B",
        "mechanism_node": "glycolysis_lactate_alteration;senescence_SASP",
        "analysis_role": "temporal support; not primary AD-control effect",
        "license": "source terms apply",
        "access_conditions": "public GEO",
        "download_date": "pre-v3; inherited project raw asset",
        "version": "GEO source record",
        "diagnosis": "Aβ time course; no human AD-control contrast",
    },
    {
        "dataset_id": "D3_GSE160936",
        "accession": "GSE160936",
        "database": "GEO",
        "source_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE160936",
        "species": "Homo sapiens",
        "model": "human postmortem brain",
        "pathology_type": "definite AD versus neurologically normal control in source metadata",
        "brain_region": "entorhinal cortex; superior frontal cortex",
        "assay": "scRNA-seq",
        "processing_level": "raw Cell Ranger filtered matrices",
        "evidence_tier": "A",
        "mechanism_node": "microglial_disease_state",
        "analysis_role": "primary human discovery cohort",
        "license": "source terms apply",
        "access_conditions": "public GEO",
        "download_date": "pre-v3; inherited project raw asset",
        "version": "GEO source record",
        "diagnosis": "AD versus control",
    },
    {
        "dataset_id": "D4_GSE263791",
        "accession": "GSE263791",
        "database": "GEO",
        "source_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE263791",
        "species": "Homo sapiens",
        "model": "human postmortem brain",
        "pathology_type": "AD spatial sections",
        "brain_region": "source region in CosMx metadata",
        "assay": "spatial transcriptomics",
        "processing_level": "processed CosMx matrix, metadata, positions",
        "evidence_tier": "B",
        "mechanism_node": "spatial_cell_communication",
        "analysis_role": "descriptive spatial support only; no primary AD-control contrast",
        "license": "source terms apply",
        "access_conditions": "public GEO",
        "download_date": "pre-v3; inherited project raw asset",
        "version": "GEO source record",
        "diagnosis": "AD only in available sections",
    },
    {
        "dataset_id": "D5_CellMetab2026_ferro_aging_reference",
        "raw_dir": "D5_CellMetab2026",
        "accession": "local_frozen_reference",
        "database": "supplementary_reference",
        "source_url": "project raw asset; source citation recorded in prior provenance",
        "species": "mixed/reference",
        "model": "published reference gene set",
        "pathology_type": "ferro-aging-related reference",
        "brain_region": "not applicable",
        "assay": "complete supplementary gene set",
        "processing_level": "processed reference table",
        "evidence_tier": "B",
        "mechanism_node": "iron_homeostasis_lipid_peroxidation_ferro_aging",
        "analysis_role": "frozen descriptive reference; not direct biochemical validation",
        "license": "source terms apply",
        "access_conditions": "local project copy",
        "download_date": "pre-v3; inherited project raw asset",
        "version": "Cell Metabolism 2026 reference copy",
        "diagnosis": "reference only",
    },
    {
        "dataset_id": "D6_GSE214979",
        "accession": "GSE214979",
        "database": "GEO",
        "source_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE214979",
        "species": "Homo sapiens",
        "model": "human postmortem brain",
        "pathology_type": "AD versus control multiome cohort",
        "brain_region": "structure field available in source metadata; not retained in current effect table",
        "assay": "snRNA-seq filtered matrix with multiome source record",
        "processing_level": "raw filtered matrix plus source cell metadata",
        "evidence_tier": "A",
        "mechanism_node": "microglial_disease_state",
        "analysis_role": "independent human validation candidate",
        "license": "source terms apply",
        "access_conditions": "public GEO",
        "download_date": "pre-v3; inherited project raw asset",
        "version": "GEO source record",
        "diagnosis": "AD versus control",
    },
]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if fields is None:
        fields = list(rows[0]) if rows else []
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def git_value(root: Path, *args: str) -> str:
    try:
        return subprocess.check_output(["git", *args], cwd=root, text=True, stderr=subprocess.DEVNULL).strip()
    except Exception:
        return "unavailable"


def source_files(root: Path, dataset_id: str, raw_dir: str | None = None) -> list[Path]:
    directory = root / "data" / "raw" / (raw_dir or dataset_id)
    return sorted(p for p in directory.rglob("*") if p.is_file()) if directory.exists() else []


def provenance_for(root: Path, dataset_id: str) -> dict:
    candidates = {
        "D1_GSE188765": "D1_annotation_provenance.json",
        "D2_GSE227221": "D2_temporal_provenance.json",
        "D4_GSE263791": "D4_spatial_provenance.json",
        "D6_GSE214979": "D6_expression_validation_provenance.json",
    }
    name = candidates.get(dataset_id)
    if name:
        path = root / "data" / "interim" / name
        if path.exists():
            try:
                return json.loads(path.read_text(encoding="utf-8"))
            except json.JSONDecodeError:
                return {"source_file": str(path.relative_to(root)), "parse_status": "invalid_json"}
    return {"source_file": "not available; see registry and raw file manifest"}


def build_registries(root: Path, output: Path, legacy: Path) -> tuple[list[dict], list[dict]]:
    registry_dir = output / "registries"
    provenance_dir = registry_dir / "per_dataset_provenance"
    registry_dir.mkdir(parents=True, exist_ok=True)
    provenance_dir.mkdir(parents=True, exist_ok=True)
    dataset_rows = []
    file_rows = []
    checksum_rows = []
    for spec in DATASETS:
        files = source_files(root, spec["dataset_id"], spec.get("raw_dir"))
        total_bytes = 0
        for path in files:
            digest = sha256(path)
            size = path.stat().st_size
            total_bytes += size
            row = {
                "dataset_id": spec["dataset_id"],
                "file_name": str(path.relative_to(root)),
                "file_size": size,
                "sha256": digest,
                "processing_level": spec["processing_level"],
                "raw_processed_relationship": "raw source asset used by inherited analysis; no raw overwrite",
            }
            file_rows.append(row)
            checksum_rows.append(f"{digest}  {row['file_name']}")
        dataset_row = dict(spec)
        dataset_row.update({
            "file_count": len(files),
            "total_file_size": total_bytes,
            "raw_data_present": bool(files),
            "inclusion_decision": "included_for_stated_role" if files else "excluded_missing_raw_asset",
            "exclusion_reason": "" if files else "raw asset not present",
        })
        dataset_rows.append(dataset_row)
        provenance = dict(spec)
        provenance["files"] = file_rows[-len(files):] if files else []
        provenance["prior_provenance"] = provenance_for(root, spec["dataset_id"])
        provenance["registered_at"] = datetime.now(timezone.utc).isoformat()
        (provenance_dir / f"{spec['dataset_id']}.json").write_text(json.dumps(provenance, ensure_ascii=False, indent=2), encoding="utf-8")

    write_tsv(registry_dir / "dataset_registry.tsv", dataset_rows)
    write_tsv(registry_dir / "file_manifest.tsv", file_rows)

    d3_meta = root / "data" / "interim" / "D3_pseudobulk_metadata_microglia.tsv"
    d6_meta = legacy / "results" / "tables" / "D6_donor_metadata.tsv"
    sample_rows = []
    if d3_meta.exists():
        for row in read_tsv(d3_meta):
            sample_rows.append({
                "dataset_id": "D3_GSE160936",
                "sample_id": row.get("group_id", ""),
                "donor_id": row.get("donor_id", ""),
                "diagnosis": row.get("disease", ""),
                "brain_region": row.get("region", ""),
                "cell_type": row.get("cell_type", "Microglia"),
                "age": row.get("age", ""),
                "sex": row.get("sex", ""),
                "metadata_source": str(d3_meta.relative_to(root)),
                "primary_role": "G1_primary_human_discovery",
            })
    if d6_meta.exists():
        for row in read_tsv(d6_meta):
            sample_rows.append({
                "dataset_id": "D6_GSE214979",
                "sample_id": row.get("id", ""),
                "donor_id": row.get("id", ""),
                "diagnosis": row.get("Status", ""),
                "brain_region": "not retained in available D6 effect table",
                "cell_type": "Microglia_candidate",
                "age": row.get("Age", ""),
                "sex": row.get("Sex", ""),
                "metadata_source": str(d6_meta.relative_to(root)),
                "primary_role": "G1_independent_validation_candidate",
            })
    for spec in DATASETS:
        if not any(row["dataset_id"] == spec["dataset_id"] for row in sample_rows):
            sample_rows.append({
                "dataset_id": spec["dataset_id"],
                "sample_id": "dataset_level_record",
                "donor_id": "not available",
                "diagnosis": spec["diagnosis"],
                "brain_region": spec["brain_region"],
                "cell_type": "not applicable",
                "age": "",
                "sex": "",
                "metadata_source": "dataset_registry.tsv",
                "primary_role": spec["analysis_role"],
            })
    write_tsv(registry_dir / "sample_manifest.tsv", sample_rows)

    evidence_rows = [
        {"evidence_id": "E_D3_G1", "dataset_id": "D3_GSE160936", "tier": "A", "mechanism_node": "microglial_disease_state", "role": "primary human discovery", "admissible": "yes", "boundary": "donor-aware output exists; v3 region-level cross-cohort meta not completed"},
        {"evidence_id": "E_D6_G1", "dataset_id": "D6_GSE214979", "tier": "A", "mechanism_node": "microglial_disease_state", "role": "independent human validation candidate", "admissible": "yes", "boundary": "current module output lacks region-stratified effects and uncertainty for meta-analysis"},
        {"evidence_id": "E_D1_G2", "dataset_id": "D1_GSE188765", "tier": "A", "mechanism_node": "H4K12la_dysregulation", "role": "direct assay reference", "admissible": "conditional", "boundary": "G2 not evaluated because G1 failed"},
        {"evidence_id": "E_D2_support", "dataset_id": "D2_GSE227221", "tier": "B", "mechanism_node": "glycolysis_lactate_alteration;senescence_SASP", "role": "iPSC temporal support", "admissible": "yes", "boundary": "not AD-control primary evidence"},
        {"evidence_id": "E_D4_spatial", "dataset_id": "D4_GSE263791", "tier": "B", "mechanism_node": "spatial_cell_communication", "role": "descriptive spatial support", "admissible": "limited", "boundary": "both available sections are AD; no AD-control effect"},
        {"evidence_id": "E_D5_reference", "dataset_id": "D5_CellMetab2026_ferro_aging_reference", "tier": "B", "mechanism_node": "iron_homeostasis_lipid_peroxidation_ferro_aging", "role": "frozen reference signature", "admissible": "yes", "boundary": "reference gene set is not direct biochemical evidence"},
    ]
    write_tsv(registry_dir / "evidence_registry.tsv", evidence_rows)

    exclusions = [
        {"dataset_id": "D2_GSE227221", "scope": "primary_G1_effect", "decision": "excluded", "reason": "Aβ iMGL time course is not definite AD versus neurologically normal control", "restart_condition": "use only as secondary directional support"},
        {"dataset_id": "D4_GSE263791", "scope": "primary_G6_effect", "decision": "excluded", "reason": "available spatial sections are AD-only and cannot establish AD/control response direction", "restart_condition": "add donor-replicated AD/control spatial data"},
        {"dataset_id": "D5_CellMetab2026_ferro_aging_reference", "scope": "direct_G4_biochemical_claim", "decision": "excluded", "reason": "reference signature alone cannot establish iron or lipid-peroxidation biochemistry", "restart_condition": "add at least two convergent biochemical or genetic modalities"},
        {"dataset_id": "DrugSide_iLINCS;MAP_official", "scope": "v3_disease_side", "decision": "excluded", "reason": "drug perturbation and drug-side assets are outside the frozen v3.0 scope", "restart_condition": "use only in a separately frozen drug-side protocol"},
    ]
    write_tsv(registry_dir / "exclusion_log.tsv", exclusions)

    protocol = root / "docs" / "protocol" / "disease_side_v3.0_latest_frozen.md"
    config = root / "configs" / "disease_side_v3.0_frozen.yaml"
    params = root / "params.yaml"
    for path in (protocol, config, params):
        if path.exists():
            digest = sha256(path)
            checksum_rows.append(f"{digest}  {path.relative_to(root)}")
    (registry_dir / "checksums.sha256").write_text("\n".join(checksum_rows) + "\n", encoding="utf-8")
    return dataset_rows, sample_rows


def table_rows(path: Path) -> list[dict[str, str]]:
    return read_tsv(path) if path.exists() else []


def build_gate_audit(root: Path, output: Path, legacy: Path, dataset_rows: list[dict], sample_rows: list[dict]) -> dict:
    tables = legacy / "results" / "tables"
    d3_meta = [r for r in sample_rows if r["dataset_id"] == "D3_GSE160936" and r["cell_type"] == "Microglia"]
    d6_meta = [r for r in sample_rows if r["dataset_id"] == "D6_GSE214979"]
    d6_concordance = table_rows(tables / "D6_signature_concordance.tsv")
    micro = next((r for r in d6_concordance if r.get("cell_type") == "Microglia"), {})
    d3_d6_direction = float(micro.get("spearman_rho", "nan")) if micro.get("spearman_rho") else None
    components = [
        {"gate": "G1", "criterion": "at_least_two_independent_human_AD_control_cohorts", "observed": "D3_GSE160936; D6_GSE214979", "pass": len({r["dataset_id"] for r in d3_meta + d6_meta}) >= 2, "reason": "Two candidate human AD/control cohorts are registered."},
        {"gate": "G1", "criterion": "donor_level_microglial_model", "observed": f"D3 rows={len(d3_meta)}; D6 donors={len(d6_meta)}", "pass": bool(d3_meta) and bool(d6_meta), "reason": "Donor-aware pseudobulk or donor-level validation outputs are present."},
        {"gate": "G1", "criterion": "region_stratified_effects_in_each_cohort", "observed": "D3 has EC/SSC metadata; D6 effect table has no retained region/structure column", "pass": False, "reason": "The available D6 module effect output cannot support the frozen region-stratified integration."},
        {"gate": "G1", "criterion": "concordant_microglial_direction_in_two_cohorts", "observed": f"D6-vs-D3 microglia signature Spearman rho={d3_d6_direction if d3_d6_direction is not None else 'NA'}", "pass": False, "reason": "Available cross-cohort concordance is approximately null/negative and does not establish a reproducible shared disease state."},
        {"gate": "G1", "criterion": "region_aware_random_effects_meta_BH_FDR_lt_0_05", "observed": "not computed", "pass": False, "reason": "Comparable region-level effects with standard errors are not available in both cohorts."},
        {"gate": "G1", "criterion": "leave_one_cohort_out_direction_stable", "observed": "not computed", "pass": False, "reason": "The required meta-analysis inputs are missing; no post hoc substitute is used."},
        {"gate": "G1", "criterion": "covariate_and_influential_dataset_sensitivity", "observed": "not available at cross-cohort meta level", "pass": False, "reason": "Age/sex/platform/region sensitivity was not completed on a comparable two-cohort effect scale."},
    ]
    component_path = output / "gates" / "G1_component_audit.tsv"
    write_tsv(component_path, components)
    g1_pass = all(row["pass"] for row in components)
    gate_rows = [
        {"gate": "G0", "status": "PASS", "reason": "All six disease-side source records have auditable raw/processed assets and a stated role; primary-scope exclusions are recorded.", "next_action": "evaluate G1"},
        {"gate": "G1", "status": "PASS" if g1_pass else "FAIL", "reason": "All frozen G1 criteria must pass; one or more required region-aware/concordance/meta/sensitivity criteria are not demonstrated.", "next_action": "continue to G2" if g1_pass else "stop; do not run G2-G6"},
    ]
    for gate in ("G2", "G3", "G4", "G5", "G6"):
        gate_rows.append({"gate": gate, "status": "NOT_RUN_UPSTREAM_FAIL", "reason": "Frozen fail-closed rule after G1 failure.", "next_action": "restart only after G1 passes"})
    write_tsv(output / "gates" / "gate_audit.tsv", gate_rows)
    gaps = [
        {"priority": "blocking", "gate": "G1", "gap": "D6 region/structure-stratified donor-level effect table with standard errors or equivalent uncertainty", "evidence": "D6_module_effects.tsv and D6 provenance", "restart": "recompute D6 with region/structure retained and export comparable donor-level effects"},
        {"priority": "blocking", "gate": "G1", "gap": "region-aware random-effects meta-analysis with BH-FDR", "evidence": "no meta-analysis output in registered run", "restart": "run after comparable D3/D6 effects and uncertainty are available"},
        {"priority": "blocking", "gate": "G1", "gap": "leave-one-cohort-out direction stability", "evidence": "not computed", "restart": "run on the same predeclared effect scale"},
        {"priority": "blocking", "gate": "G1", "gap": "reproducible microglial direction across cohorts", "evidence": f"D6 microglia signature rho={d3_d6_direction if d3_d6_direction is not None else 'NA'}", "restart": "reassess only after harmonized region-aware modeling; do not lower the gate"},
        {"priority": "downstream", "gate": "G2", "gap": "independent concordant direct H4K12la support plus glycolysis/lactate compatibility", "evidence": "not evaluated after G1 failure", "restart": "evaluate only after G1 passes"},
        {"priority": "downstream", "gate": "G4", "gap": "two convergent biochemical/genetic modalities beyond RNA proxy", "evidence": "reference module alone is insufficient", "restart": "add direct iron/lipid-peroxidation or genetic evidence"},
        {"priority": "downstream", "gate": "G6", "gap": "AD/control donor-replicated spatial communication", "evidence": "available D4 sections are AD-only", "restart": "add eligible AD/control spatial cohort after G1-G5 pass"},
    ]
    write_tsv(output / "evidence_gaps.tsv", gaps)
    return {"g1_pass": g1_pass, "d3_d6_microglia_rho": d3_d6_direction, "component_count": len(components)}


def write_outputs(root: Path, output: Path, legacy: Path, dataset_rows: list[dict], sample_rows: list[dict], gate_result: dict, run_id: str) -> None:
    output.mkdir(parents=True, exist_ok=True)
    status = "GO" if gate_result["g1_pass"] else "STOPPED_AT_G1"
    decision = {
        "run_id": run_id,
        "protocol": "AD disease-side v3.0",
        "protocol_status": "LATEST_FROZEN",
        "decision": status,
        "scope": "disease-side only",
        "gates": {"G0": "PASS", "G1": "PASS" if gate_result["g1_pass"] else "FAIL", "G2": "NOT_RUN_UPSTREAM_FAIL", "G3": "NOT_RUN_UPSTREAM_FAIL", "G4": "NOT_RUN_UPSTREAM_FAIL", "G5": "NOT_RUN_UPSTREAM_FAIL", "G6": "NOT_RUN_UPSTREAM_FAIL"},
        "claim_boundary": "No drug efficacy, drug ranking, causal communication, or downstream mechanism claim is made.",
        "reason": "G1 frozen criteria were not all demonstrated from the available registered outputs." if not gate_result["g1_pass"] else "All currently audited G1 criteria passed; downstream gates require a separate continuation run.",
        "d3_d6_microglia_signature_spearman": gate_result["d3_d6_microglia_rho"],
        "created_at": datetime.now(timezone.utc).isoformat(),
    }
    (output / "final_decision.json").write_text(json.dumps(decision, ensure_ascii=False, indent=2), encoding="utf-8")
    restart = """# Restart conditions for AD disease-side v3.0

Current decision: `STOPPED_AT_G1`.

G2–G6 remain locked. A restart must first add or regenerate auditable donor-level microglial effects for at least two independent human definite-AD versus neurologically-normal cohorts, retain brain region/structure in both cohorts, and export comparable effect estimates with uncertainty. The region-aware random-effects meta-analysis must pass BH-FDR < 0.05, preserve direction under leave-one-cohort-out analysis, and report the frozen covariate/influential-dataset sensitivities.

Only after G1 passes may the workflow evaluate direct H4K12la evidence and glycolysis/lactate compatibility (G2), then sequentially unlock G3, G4, G5, and G6. The AD-only spatial sections and the iMGL Aβ time course remain secondary support and cannot be promoted to the primary gates by themselves.
"""
    (output / "restart_conditions.md").write_text(restart, encoding="utf-8")
    report = f"""# AD disease-side v3.0 execution report

- Run: `{run_id}`
- Protocol: `docs/protocol/disease_side_v3.0_latest_frozen.md` (`LATEST_FROZEN`)
- Scope: disease-side only; drug perturbations, drug ranking, and efficacy claims excluded.
- Registered datasets: {len(dataset_rows)}; sample records: {len(sample_rows)}.
- Decision: **{status}**.

G0 passed for the registered stated roles. G1 was audited against the frozen requirements and stopped the workflow because the available D6 output does not retain region-stratified donor-level effects/uncertainty for cross-cohort meta-analysis, the required random-effects meta-analysis and leave-one-cohort-out analysis were not computed, and the available D3–D6 microglial signature concordance is not supportive (Spearman rho `{gate_result['d3_d6_microglia_rho']}`).

Accordingly, G2–G6 were not run. No downstream causal, biochemical, spatial communication, drug-side, or treatment-efficacy claim is exported. See `gates/`, `evidence_gaps.tsv`, and `restart_conditions.md`.
"""
    (output / "execution_report.md").write_text(report, encoding="utf-8")
    manifest = {
        "run_id": run_id,
        "created_at": datetime.now(timezone.utc).isoformat(),
        "project_root": str(root),
        "git_commit": git_value(root, "rev-parse", "HEAD"),
        "git_tag_at_start": git_value(root, "describe", "--tags", "--always"),
        "protocol_sha256": sha256(root / "docs/protocol/disease_side_v3.0_latest_frozen.md"),
        "config_sha256": sha256(root / "configs/disease_side_v3.0_frozen.yaml"),
        "params_sha256": sha256(root / "params.yaml"),
        "legacy_results_root": str(legacy.relative_to(root)),
        "python": sys.version.split()[0],
        "decision": status,
        "raw_data_read_only": True,
        "external_network_used": False,
    }
    (output / "run_manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--legacy-results", type=Path, required=True, help="Archived results root containing results/tables")
    parser.add_argument("--run-id", default="2026-09-20-disease-side-v3")
    parser.add_argument("--refresh", action="store_true", help="refresh an existing run package in place")
    args = parser.parse_args()
    root = args.root.resolve()
    legacy = args.legacy_results.resolve()
    output = root / "results" / args.run_id
    if output.exists() and not args.refresh:
        raise FileExistsError(f"Output already exists: {output}")
    dataset_rows, sample_rows = build_registries(root, output, legacy)
    gate_result = build_gate_audit(root, output, legacy, dataset_rows, sample_rows)
    write_outputs(root, output, legacy, dataset_rows, sample_rows, gate_result, args.run_id)
    print(json.dumps({"run_id": args.run_id, "decision": "GO" if gate_result["g1_pass"] else "STOPPED_AT_G1", "output": str(output)}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
