"""Create the auditable v3.1 decision package after G1 execution."""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RUN = ROOT / "results" / "2026-09-20-disease-side-v3.1"
REG = RUN / "registries"
GATES = RUN / "gates"
RUN.mkdir(parents=True, exist_ok=True)
REG.mkdir(parents=True, exist_ok=True)
GATES.mkdir(parents=True, exist_ok=True)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if fields is None:
        fields = list(rows[0]) if rows else []
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def numeric(value: str) -> object:
    if value in ("", "NA", "NaN"):
        return ""
    try:
        f = float(value)
        return int(f) if f.is_integer() else f
    except ValueError:
        return value


def build_sample_manifest() -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    d7 = read_tsv(ROOT / "data" / "interim" / "D7_GSE148822" / "D7_sample_manifest.tsv")
    for row in d7:
        rows.append({
            "cohort": row["cohort"], "accession": row["accession"], "donor_id": row["donor"], "region": row["region"],
            "diagnosis": "Control" if row["sample_group"] == "CTR" else row["sample_group"], "sample_group": row["sample_group"],
            "age": numeric(row["age"]), "sex": row["sex"], "RIN": numeric(row["RIN"]), "PMI": "", "batch": "",
            "include_primary": row["include_primary"], "decision": row["decision"], "n_cells_used": "", "source_file": row["matrix_file"],
        })
    d8 = read_tsv(ROOT / "data" / "interim" / "D8_GSE174367" / "D8_pseudobulk_metadata_microglia.tsv")
    for row in d8:
        rows.append({
            "cohort": row["cohort"], "accession": row["donor_id"], "donor_id": row["donor_id"], "region": row["region"],
            "diagnosis": row["disease"], "sample_group": row["disease"], "age": numeric(row["age"]), "sex": row["sex"],
            "RIN": numeric(row["RIN"]), "PMI": numeric(row["PMI"]), "batch": row["batch"], "include_primary": "true",
            "decision": "primary_annotated_MG", "n_cells_used": numeric(row["n_cells_used"]), "source_file": "GSE174367_filtered_feature_bc_matrix.h5",
        })
    d3 = read_tsv(ROOT / "data" / "interim" / "D3_pseudobulk_metadata_microglia.tsv")
    for row in d3:
        rows.append({
            "cohort": "D3_GSE160936", "accession": row.get("group_id", ""), "donor_id": row["donor_id"], "region": row["region"],
            "diagnosis": row["disease"], "sample_group": row["disease"], "age": numeric(row["age"]), "sex": row["sex"],
            "RIN": numeric(row.get("rin", "")), "PMI": "", "batch": "", "include_primary": "true", "decision": "primary_existing_QC",
            "n_cells_used": numeric(row.get("n_cells", "")), "source_file": "D3_existing_processed_pseudobulk",
        })
    return rows


sample_rows = build_sample_manifest()
write_tsv(REG / "sample_manifest.tsv", sample_rows, ["cohort", "accession", "donor_id", "region", "diagnosis", "sample_group", "age", "sex", "RIN", "PMI", "batch", "include_primary", "decision", "n_cells_used", "source_file"])

dataset_rows = [
    {
        "accession": "D3_GSE160936", "database": "NCBI GEO", "source_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE160936", "download_date": "2026-09-20", "version": "GEO public processed matrices",
        "license": "NCBI GEO public data", "access_conditions": "public", "file_name": "D3 existing processed Cell Ranger-derived pseudobulk inputs", "file_size": "see file_manifest", "sha256": "see checksums.sha256", "processing_level": "processed integer pseudobulk counts",
        "raw_processed_relationship": "derived from downloaded processed Cell Ranger matrices; raw source retained read-only", "species": "Homo sapiens", "model": "human AD brain", "pathology_type": "definite AD vs normal control", "brain_region": "EC; SSC", "layer_or_ROI": "not applicable", "age": "see sample_manifest", "sex": "see sample_manifest", "diagnosis": "AD/Control", "donor_count": "12", "sample_count": "24", "cell_count": "see D3 metadata", "assay": "snRNA-seq", "platform": "10x Cell Ranger", "antibody": "not applicable", "genome_build": "GRCh38", "covariates": "age, sex, RIN", "pairing_structure": "paired EC/SSC donor regions", "donor_overlap_audit": "no detectable overlap with D7/D8 labels", "inclusion_decision": "primary G1 discovery", "exclusion_reason": "none", "evidence_tier": "B", "mechanism_node": "microglial disease state", "analysis_role": "discovery",
    },
    {
        "accession": "GSE148822", "database": "NCBI GEO", "source_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE148822", "download_date": "2026-09-20", "version": "GEO public processed Cell Ranger MTX", "license": "NCBI GEO public data", "access_conditions": "public", "file_name": "30 strict sample MTX/feature/barcode triplets", "file_size": "see file_manifest", "sha256": "see checksums.sha256", "processing_level": "processed integer Cell Ranger MTX", "raw_processed_relationship": "processed filtered/feature matrix; FASTQ/SRA not required by v3.1", "species": "Homo sapiens", "model": "human AD brain", "pathology_type": "definite AD vs strict normal control; CTR+ excluded", "brain_region": "OC; OTC", "layer_or_ROI": "not applicable", "age": "see sample_manifest", "sex": "see sample_manifest", "diagnosis": "10 AD + 5 Control per region", "donor_count": "15", "sample_count": "30 primary; 36 series", "cell_count": "source matrix and QC summary", "assay": "snRNA-seq", "platform": "10x Cell Ranger v3.0.1", "antibody": "NeuN-/OLIG2- sorted nuclei", "genome_build": "GRCh38", "covariates": "age, sex, RIN", "pairing_structure": "paired OC/OTC donor regions", "donor_overlap_audit": "no detectable overlap with D3/D8 labels", "inclusion_decision": "primary G1 validation", "exclusion_reason": "CTR+ donors 6-8 excluded", "evidence_tier": "B", "mechanism_node": "microglial disease state", "analysis_role": "primary validation",
    },
    {
        "accession": "GSE174367", "database": "NCBI GEO", "source_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE174367", "download_date": "2026-09-20", "version": "GEO public filtered H5 + cell metadata", "license": "NCBI GEO public data", "access_conditions": "public", "file_name": "filtered_feature_bc_matrix.h5 + cell_meta.csv.gz", "file_size": "see file_manifest", "sha256": "see checksums.sha256", "processing_level": "processed integer Cell Ranger H5", "raw_processed_relationship": "processed H5; FASTQ/SRA not required by v3.1", "species": "Homo sapiens", "model": "human AD brain", "pathology_type": "AD vs Control", "brain_region": "PFC", "layer_or_ROI": "not applicable", "age": "see sample_manifest", "sex": "see sample_manifest", "diagnosis": "11 AD + 7 Control", "donor_count": "18", "sample_count": "18", "cell_count": "4126 MG source-annotated; 298 unlabeled H5 barcodes excluded", "assay": "snRNA-seq", "platform": "10x filtered H5", "antibody": "not applicable", "genome_build": "GRCh38", "covariates": "age, sex, RIN, PMI, batch", "pairing_structure": "one PFC region per donor", "donor_overlap_audit": "no detectable overlap with D3/D7 labels", "inclusion_decision": "primary G1 secondary validation", "exclusion_reason": "Sample-100 excluded only from complete-PMI sensitivity", "evidence_tier": "B", "mechanism_node": "microglial disease state", "analysis_role": "secondary validation",
    },
]
fields = list(dataset_rows[0])
write_tsv(REG / "dataset_registry.tsv", dataset_rows, fields)

evidence_rows = [
    {"evidence_id": "G1-v3.1", "node": "microglial disease state", "source": "D3/GSE160936 + D7/GSE148822 + D8/GSE174367", "tier": "B", "analysis": "donor-region pseudobulk; random-effects meta; LOCO", "decision": "PASS", "claim_boundary": "reproducible human AD microglial association; not causal"},
    {"evidence_id": "G2-direct-H4K12la", "node": "glycolysis/lactate-H4K12la", "source": "registered disease-side evidence package", "tier": "none", "analysis": "accepted direct H4K12la assay audit", "decision": "FAIL", "claim_boundary": "D1 chromatin peak annotation and RNA/module projections are not direct H4K12la measurement"},
    {"evidence_id": "G3-G6", "node": "senescence/SASP; ferro-aging; astrocyte; spatial communication", "source": "v3.1 fail-closed execution", "tier": "not_run", "analysis": "upstream G2 failure", "decision": "STOPPED_UPSTREAM", "claim_boundary": "no downstream gate or causal communication claim"},
]
write_tsv(REG / "evidence_registry.tsv", evidence_rows)

exclusion_rows = [
    {"dataset_or_sample": "D7/GSE148822 Donor6-8", "scope": "primary G1", "reason": "CTR+ excluded by frozen protocol", "action": "retained in metadata; not downloaded for primary matrices"},
    {"dataset_or_sample": "D6/GSE214979", "scope": "G1", "reason": "no region with >=2 AD and >=2 Control donors", "action": "exploratory only; excluded from G1 numerator/denominator"},
    {"dataset_or_sample": "D9/GSE263468-BA9", "scope": "G1", "reason": "predeclared sensitivity only", "action": "not used to rescue or replace primary cohorts"},
    {"dataset_or_sample": "D8 H5 298 barcodes", "scope": "D8 cell input", "reason": "no matching GEO cell metadata", "action": "excluded before MG pseudobulk; auditable count recorded"},
    {"dataset_or_sample": "D8 Sample-100", "scope": "complete-covariate sensitivity", "reason": "PMI missing", "action": "kept in primary model; excluded in complete-PMI sensitivity"},
    {"dataset_or_sample": "GSE264648", "scope": "independent G1", "reason": "detectable donor overlap with D3", "action": "non-independent sensitivity only; not run in hard gate"},
]
write_tsv(REG / "exclusion_log.tsv", exclusion_rows)

g2_rows = [
    {"criterion": "accepted_direct_assay_present", "required": "H4K12la ChIP-seq/CUT&Tag/CUT&RUN/quantitative proteomics", "observed": "0 in registered v3.1 inputs", "pass": "FALSE", "reason": "D1 is peak annotation/reference projection, not H4K12la-specific direct measurement"},
    {"criterion": "independent_concordant_support", "required": "yes", "observed": "not assessable after direct assay failure", "pass": "FALSE", "reason": "G2 fail-closed"},
    {"criterion": "downstream_execution", "required": "only after G2 pass", "observed": "not run", "pass": "FALSE", "reason": "G3-G6 locked by protocol"},
]
write_tsv(GATES / "G2_assessment.tsv", g2_rows)
gate_summary = [
    {"gate": "G0", "decision": "PASS", "reason": "All three primary cohorts have auditable processed integer counts, donor/diagnosis/region metadata, and QC support."},
    {"gate": "G1", "decision": "PASS", "reason": "All five region units estimable; 413 random-effects gene meta FDR<0.05; 244 hits cover all five units; direction stable in all three LOCO analyses."},
    {"gate": "G2", "decision": "FAIL", "reason": "No accepted direct H4K12la assay in registered evidence; transcriptomic/chromatin proxies cannot pass this gate."},
    {"gate": "G3", "decision": "NOT_RUN_UPSTREAM_FAIL", "reason": "G2 failed."},
    {"gate": "G4", "decision": "NOT_RUN_UPSTREAM_FAIL", "reason": "G2 failed."},
    {"gate": "G5", "decision": "NOT_RUN_UPSTREAM_FAIL", "reason": "G2 failed."},
    {"gate": "G6", "decision": "NOT_RUN_UPSTREAM_FAIL", "reason": "G2 failed; no spatial communication claim."},
]
write_tsv(GATES / "gate_summary.tsv", gate_summary)

now = datetime.now(timezone.utc).isoformat()
decision = {
    "run_id": "2026-09-20-disease-side-v3.1",
    "protocol": "AD disease-side v3.1",
    "protocol_status": "LATEST_FROZEN",
    "decision": "STOPPED_AT_G2",
    "gates": {row["gate"]: row["decision"] for row in gate_summary},
    "g1": {"primary_region_units": 5, "estimable": True, "meta_genes_fdr_lt_0_05": 413, "strict_all_5_region_hits": 244, "loo_direction_stable": True, "frozen_module_meta_signal": False},
    "g2": {"direct_h4k12la_assays": 0, "decision": "FAIL"},
    "downstream_lock": "G3-G6 not run after G2 failure",
    "claim_boundary": "Human AD microglial disease association is reproducible across three primary cohorts; no direct H4K12la, downstream mechanism, spatial communication, drug ranking, or efficacy claim.",
    "created_at": now,
}
(RUN / "final_decision.json").write_text(json.dumps(decision, ensure_ascii=False, indent=2), encoding="utf-8")

report = f"""# Disease-side v3.1 execution report

Run: `2026-09-20-disease-side-v3.1`  
Protocol: `LATEST_FROZEN` v3.1  
Execution date: 2026-09-20

## Result

The modified disease-side protocol was executed with processed public count matrices. G0 passed and G1 passed. G1 used five predeclared region units: D3 EC/SSC, D7 OC/OTC, and D8 PFC. Random-effects meta-analysis identified 413 genes at BH-FDR < 0.05; 244 significant genes had effects in all five region units, and all three leave-one-primary-cohort-out analyses preserved their directions. The frozen M1/M3/M5/M6/M7 module composites did not pass BH-FDR, so the positive G1 conclusion is explicitly gene-level rather than a claim that every operational module is replicated.

G2 failed because the registered evidence package contains no accepted direct H4K12la assay. D1 chromatin peak annotation, RNA module scores, and cross-model projections are proxies/reference layers and cannot be promoted to direct H4K12la measurement. Per the frozen fail-closed rule, G3–G6 were not run. The final decision is `STOPPED_AT_G2`.

## Inputs and processing

- D3/GSE160936: existing auditable donor×region microglia pseudobulk input; 6 AD and 6 controls in each EC/SSC region.
- D7/GSE148822: 30 primary processed Cell Ranger MTX triplets; OC/OTC each 10 AD and 5 strict controls; CTR+ donors excluded.
- D8/GSE174367: official filtered H5 plus cell metadata; 11 AD and 7 controls; GEO Cell.Type==MG annotation; 298 H5 barcodes lacking metadata excluded.
- No FASTQ/SRA was downloaded. This is allowed by v3.1 because the processed integer count matrices preserved the required donor, region, diagnosis, and microglia information.
- Cell QC used `nCount>=500`, `nFeature>=200`, and `mt_pct<=20`; D7 additionally required at least two microglia markers and microglia score >= astrocyte score. D8 used the source MG annotation. Every included donor-region exceeded the frozen 20-cell minimum.
- `scDblFinder` was unavailable in the local R 4.3.3 environment, so doublet detection was not claimed; this limitation is recorded in every QC metadata row.

## Reproducibility and audit

The run includes dataset/sample/file/evidence/exclusion registries, per-dataset provenance, checksums, region-specific effects, random-effects meta, cohort-pooled effects, LOCO, module meta, complete-covariate and cell-count sensitivity outputs, G0/G1/G2 gate tables, and machine-readable final decision. DVC inspection remains blocked by the existing shared-cache Windows access-denied condition; raw inputs remain read-only and the limitation is recorded rather than bypassed.

Previous v2/v3 results were retained as historical archives and were not overwritten. The v3.0 STOPPED_AT_G1 packages remain separate from this v3.1 run.
"""
(RUN / "execution_report.md").write_text(report, encoding="utf-8")

provenance = {
    "D3_GSE160936": {"source_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE160936", "role": "discovery", "input_level": "existing processed integer pseudobulk derived from Cell Ranger matrices", "download_date": "pre-existing; reused read-only"},
    "D7_GSE148822": {"source_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE148822", "paper_doi": "10.1007/s00401-021-02263-w", "role": "primary validation", "input_level": "processed Cell Ranger MTX", "selected_samples": 30, "strict_counts": {"OC_AD": 10, "OC_Control": 5, "OTC_AD": 10, "OTC_Control": 5}},
    "D8_GSE174367": {"source_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE174367", "paper_doi": "10.1038/s41588-021-00894-z", "role": "secondary validation", "input_level": "processed Cell Ranger H5 + cell metadata", "samples": 18, "source_MG_cells": 4126, "unmapped_H5_barcodes_excluded": 298},
}
(REG / "per_dataset_provenance.json").write_text(json.dumps(provenance, ensure_ascii=False, indent=2), encoding="utf-8")

manifest = {
    "run_id": "2026-09-20-disease-side-v3.1",
    "protocol": "docs/protocol/disease_side_v3.1_latest_frozen.md",
    "config": "configs/disease_side_v3.1_frozen.yaml",
    "entrypoints": ["src/22_prepare_disease_side_v3_1_inputs.py", "src/23_prepare_d7_d8_pseudobulk.R", "src/24_validate_d7_mtx.py", "src/25_run_g1_v3_1.R", "src/26_finalize_disease_side_v3_1.py"],
    "environment": {"python": "D:/python/python.exe", "Rscript": "D:/R/R-4.3.3/bin/Rscript.exe", "R_packages": ["Seurat", "hdf5r", "rhdf5", "edgeR", "limma", "metafor"], "scDblFinder": "unavailable"},
    "statistical_unit": "donor within region and microglia",
    "meta_method": "DerSimonian-Laird random-effects on region-unit logFC and moderated SE; BH-FDR",
    "result": "STOPPED_AT_G2",
    "old_results": ["results/archive/2026-09-20-pre-v3-freeze", "results/archive/2026-09-20-pre-v3.1"],
}
(RUN / "run_manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")

tracked = [
    ROOT / "data" / "raw" / "D7_GSE148822" / "GSE148822_family.soft.gz",
    ROOT / "data" / "raw" / "D8_GSE174367" / "GSE174367_snRNA-seq_cell_meta.csv.gz",
    ROOT / "data" / "raw" / "D8_GSE174367" / "GSE174367_snRNA-seq_filtered_feature_bc_matrix.h5",
    ROOT / "data" / "interim" / "D3_pseudobulk_counts_microglia.tsv.gz",
    ROOT / "data" / "interim" / "D3_pseudobulk_metadata_microglia.tsv",
    ROOT / "data" / "interim" / "D7_GSE148822" / "D7_sample_manifest.tsv",
    ROOT / "data" / "interim" / "D7_GSE148822" / "D7_mtx_integrity_audit.json",
    ROOT / "data" / "interim" / "D7_GSE148822" / "D7_pseudobulk_counts_microglia.tsv.gz",
    ROOT / "data" / "interim" / "D7_GSE148822" / "D7_pseudobulk_metadata_microglia.tsv",
    ROOT / "data" / "interim" / "D7_GSE148822" / "D7_qc_summary_microglia.tsv",
    ROOT / "data" / "interim" / "D8_GSE174367" / "D8_file_audit.json",
    ROOT / "data" / "interim" / "D8_GSE174367" / "D8_pseudobulk_counts_microglia.tsv.gz",
    ROOT / "data" / "interim" / "D8_GSE174367" / "D8_pseudobulk_metadata_microglia.tsv",
    ROOT / "data" / "interim" / "D8_GSE174367" / "D8_qc_summary_microglia.tsv",
]
tracked.extend(sorted((ROOT / "data" / "raw" / "D7_GSE148822").glob("GSM*.gz")))
tracked.extend(sorted(RUN.rglob("*")))
tracked = [p for p in tracked if p.is_file() and p.name not in {"checksums.sha256", "file_manifest.tsv"}]
checksum_rows = []
file_rows = []
seen = set()
for path in tracked:
    rel = path.relative_to(ROOT).as_posix()
    if rel in seen:
        continue
    seen.add(rel)
    digest = sha256(path)
    size = path.stat().st_size
    checksum_rows.append(f"{digest}  {rel}")
    file_rows.append({"file_name": rel, "file_size": size, "sha256": digest, "processing_level": "raw_public_processed" if "/raw/" in rel else "interim_or_result", "read_only": "true" if "/raw/" in rel else "false"})
write_tsv(REG / "file_manifest.tsv", file_rows, ["file_name", "file_size", "sha256", "processing_level", "read_only"])
(REG / "checksums.sha256").write_text("\n".join(checksum_rows) + "\n", encoding="utf-8")
print(json.dumps({"run": str(RUN), "decision": decision["decision"], "tracked_files": len(file_rows)}, ensure_ascii=False, indent=2))
