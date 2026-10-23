#!/usr/bin/env Rscript

# Effect-blind coverage audit for the frozen G1 sender and G5 receiver genes.

suppressPackageStartupMessages({
  library(Matrix)
})
args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) args[[1L]] else "."
audit_dir <- file.path(root, "results/2026-09-24-disease-side-v3.5-g6-input-audit")
reg_dir <- file.path(audit_dir, "registries")
table_dir <- file.path(audit_dir, "tables")
interim_dir <- file.path(root, "data/interim/G6_GSE269906_bin50_tissue")
g1 <- read.delim(gzfile(file.path(reg_dir, "G1_microglia_signed_sender_244.tsv.gz")), check.names = FALSE, stringsAsFactors = FALSE)
g5 <- read.delim(gzfile(file.path(reg_dir, "G5_receiver_one_to_one_human_ortholog_registry.tsv.gz")), check.names = FALSE, stringsAsFactors = FALSE)
qc <- read.delim(file.path(table_dir, "GSE269906_tissue_bin50_section_qc.tsv"), check.names = FALSE, stringsAsFactors = FALSE)
g1$gene <- toupper(g1$gene)
g5$human_one_to_one_symbol <- toupper(g5$human_one_to_one_symbol)
g5_pass <- g5[g5$G5_module_gate == "PASS" & nzchar(g5$human_one_to_one_symbol), , drop = FALSE]
modules <- sort(unique(g5_pass$module))
signatures <- c(list(G1_signed_sender = unique(g1$gene)),
                setNames(lapply(modules, function(module) unique(g5_pass$human_one_to_one_symbol[g5_pass$module == module])), modules))

summary_rows <- list()
gene_rows <- list()
for (k in seq_len(nrow(qc))) {
  donor <- qc[k, , drop = FALSE]
  rds_path <- file.path(interim_dir, paste0(donor$GSM, "_bin50_tissue_signature_counts.rds"))
  obj <- readRDS(rds_path)
  counts <- obj$counts
  rownames(counts) <- toupper(rownames(counts))
  for (signature in names(signatures)) {
    genes <- signatures[[signature]]
    row_idx <- match(genes, rownames(counts))
    present <- !is.na(row_idx)
    prevalence <- rep(0, length(genes))
    if (any(present)) prevalence[present] <- Matrix::rowMeans(counts[row_idx[present], , drop = FALSE] > 0)
    summary_rows[[length(summary_rows) + 1L]] <- data.frame(
      accession = "GSE269906", GSM = donor$GSM, donor_id = donor$donor_id,
      signature = signature, n_frozen_genes = length(genes),
      n_in_candidate_matrix = sum(present), n_detected_in_any_accepted_bin = sum(prevalence > 0),
      n_detected_in_at_least_1pct_bins = sum(prevalence >= 0.01),
      n_detected_in_at_least_5pct_bins = sum(prevalence >= 0.05),
      n_detected_in_at_least_10pct_bins = sum(prevalence >= 0.10),
      median_gene_prevalence = median(prevalence), n_accepted_bins = ncol(counts),
      stringsAsFactors = FALSE
    )
    gene_rows[[length(gene_rows) + 1L]] <- data.frame(
      accession = "GSE269906", GSM = donor$GSM, donor_id = donor$donor_id,
      signature = signature, gene = genes, in_candidate_matrix = present,
      accepted_bin_prevalence = prevalence,
      stringsAsFactors = FALSE
    )
  }
  message("effect_blind_signature_coverage=", donor$GSM, "; bins=", ncol(counts))
  rm(obj, counts)
  gc(verbose = FALSE)
}
coverage <- do.call(rbind, summary_rows)
gene_coverage <- do.call(rbind, gene_rows)
write.table(coverage, file.path(table_dir, "GSE269906_frozen_signature_coverage_by_donor.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "")
write.table(gene_coverage, file.path(table_dir, "GSE269906_frozen_signature_gene_prevalence.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "")
cat("donors=", nrow(qc), "\n", sep = "")
cat("passing_G5_modules=", paste(modules, collapse = ";"), "\n", sep = "")
