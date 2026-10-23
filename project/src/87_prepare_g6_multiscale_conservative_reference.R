#!/usr/bin/env Rscript
suppressPackageStartupMessages({ library(Matrix); library(jsonlite) })
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: 87_prepare_g6_multiscale_conservative_reference.R <project_root>")
root <- args[[1L]]
old <- file.path(root, "data/interim/G6_GSE269906_bin50_tissue/G6_RCTD_reference_cells.rds")
audit <- file.path(root, "results/2026-09-24-disease-side-v3.5-g6-input-audit/GSE269906_GSE174367_overlap_audit.json")
outdir <- file.path(root, "data/interim/G6_GSE269906_multiscale_001/reference")
if (!file.exists(old) || !file.exists(audit)) stop("Missing original reference or overlap audit")
overlap <- jsonlite::fromJSON(audit)
if (overlap$primary_accession != "GSE269906" || overlap$G1_reference_accession != "GSE174367" ||
    overlap$demographic_diagnosis_candidate_pairs < 1L || length(overlap$GSE174367_candidate_donors) < 1L)
  stop("Unexpected or empty donor overlap audit")
ref <- readRDS(old)
if (!identical(colnames(ref$counts), names(ref$donor_id)) ||
    nrow(ref$selected_metadata) != ncol(ref$counts) || length(ref$selected_umi) != ncol(ref$counts))
  stop("Original reference metadata/counts alignment failed")
excluded <- as.character(overlap$GSE174367_candidate_donors)
if (!all(excluded %in% unique(ref$donor_id))) stop("Candidate donor missing from reference")
keep <- !(ref$donor_id %in% excluded)
ref$counts <- ref$counts[, keep, drop = FALSE]
ref$cell_type <- droplevels(ref$cell_type[keep])
ref$donor_id <- ref$donor_id[keep]
ref$selected_metadata <- ref$selected_metadata[keep, , drop = FALSE]
ref$selected_umi <- ref$selected_umi[keep]
if (length(unique(ref$donor_id)) < 8L || length(levels(ref$cell_type)) != 7L ||
    any(table(ref$cell_type) < 25L)) stop("Conservative reference lacks adequate donors/types")
ref$source$multiscale_note <- "Conservative exclusion of all GSE269906 demographic/diagnosis candidate-overlap GSE174367 donors; not a proof of identity crosswalk"
ref$source$excluded_candidate_donors <- excluded
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
outfile <- file.path(outdir, "G6_RCTD_reference_no_candidate_overlap.rds")
saveRDS(ref, outfile, compress = "gzip")
record <- list(pilot_id = "G6-MS-2026-09-26-001", excluded_reference_donors = excluded,
               remaining_donors = length(unique(ref$donor_id)), remaining_cells = ncol(ref$counts),
               cell_type_counts = as.list(table(ref$cell_type)),
               source_md5 = unname(tools::md5sum(old)), audit_md5 = unname(tools::md5sum(audit)),
               output_md5 = unname(tools::md5sum(outfile)),
               limitation = "Conservative metadata exclusion; unresolved donor crosswalk and source-cohort independence remain limitations")
jsonlite::write_json(record, file.path(outdir, "reference_exclusion_audit.json"), auto_unbox = TRUE, pretty = TRUE)
cat("remaining_reference_donors=", record$remaining_donors, "; remaining_cells=", record$remaining_cells, "\n", sep = "")
