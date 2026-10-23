#!/usr/bin/env Rscript

# Post-gate diagnostic: quantify whether the frozen identity markers are
# detectable in QC-passing CellBin nuclei. This does not relabel any nucleus.
suppressPackageStartupMessages({ library(Matrix); library(data.table) })
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: 95_audit_g6_cellbin_marker_availability.R <project_root>")
root <- args[[1L]]
input_dir <- file.path(root, "data/interim/G6_next_feasibility_001/cellbin")
audit_dir <- file.path(root, "results/2026-09-24-disease-side-v3.5-g6-input-audit")
out_dir <- file.path(root, "results/2026-09-27-g6-next-feasibility/cellbin")
ann_path <- file.path(out_dir, "annotated_nuclei_metadata.tsv")
if (!file.exists(ann_path)) stop("Primary blind annotation output is required")
gsm <- c("GSM8330060", "GSM8330061", "GSM8330062", "GSM8330063", "GSM8330064", "GSM8330065",
         "GSM8330066", "GSM8330067", "GSM8330068", "GSM8330069", "GSM8330070", "GSM8330071")
g1 <- fread(file.path(audit_dir, "registries/G1_microglia_signed_sender_244.tsv.gz"))
g5 <- fread(file.path(audit_dir, "registries/G5_receiver_one_to_one_human_ortholog_registry.tsv.gz"))
state <- unique(toupper(c(g1$gene, g5$human_one_to_one_symbol)))
state <- state[!is.na(state) & nzchar(state)]
sets <- list(
  MG_training = c("CD74", "CD14"),
  MG_holdout_disjoint_from_G1_G5 = setdiff(c("AIF1", "C1QA", "C1QB", "C1QC", "CX3CR1", "LST1", "P2RY12", "TMEM119", "TREM2", "TYROBP"), state),
  ASC_training = c("AQP4", "GLUL"),
  ASC_holdout_disjoint_from_G1_G5 = setdiff(c("ALDH1L1", "GFAP", "GJA1", "S100B", "SLC1A2", "SLC1A3", "SOX9"), state)
)
rows <- list()
for (sample in gsm) {
  x <- readRDS(file.path(input_dir, paste0(sample, "_nucleus_counts.rds")))
  meta <- x$cell_metadata
  keep <- meta$n_genes >= 150L & meta$polygon_valid
  if (any(grepl("diagnos|patholog|braak|cerad|thal|control", tolower(names(meta))))) stop("Forbidden metadata field")
  counts <- x$counts[, keep, drop = FALSE]
  genes_upper <- toupper(rownames(counts))
  for (set_name in names(sets)) for (gene in unique(sets[[set_name]])) {
    ix <- match(gene, genes_upper)
    n_positive <- if (is.na(ix)) NA_integer_ else as.integer(sum(counts[ix, ] > 0))
    rows[[length(rows) + 1L]] <- data.frame(
      GSM = sample, donor_id = unique(meta$donor_id), marker_set = set_name, gene = gene,
      in_gene_catalog = !is.na(ix), n_QC_nuclei = ncol(counts), n_positive_nuclei = n_positive,
      detection_fraction = if (is.na(ix) || ncol(counts) == 0L) NA_real_ else n_positive / ncol(counts),
      meets_1pct_detection = !is.na(ix) && ncol(counts) > 0L && n_positive / ncol(counts) >= 0.01,
      stringsAsFactors = FALSE)
  }
  message("MARKER_AVAILABILITY ", sample, "; QC nuclei=", ncol(counts))
  rm(x, counts, meta); gc(verbose = FALSE)
}
tab <- rbindlist(rows, fill = TRUE)
summary <- tab[, .(genes_in_catalog = sum(in_gene_catalog), genes_total = .N,
  genes_meeting_1pct = sum(meets_1pct_detection, na.rm = TRUE),
  median_gene_detection_fraction = median(detection_fraction, na.rm = TRUE)),
  by = .(marker_set, GSM, donor_id)]
fwrite(tab, file.path(out_dir, "marker_availability_by_donor.tsv"), sep = "\t")
fwrite(summary, file.path(out_dir, "marker_availability_summary.tsv"), sep = "\t")
writeLines(c(
  "# CellBin identity-marker availability diagnostic", "",
  "Diagnosis-blind, post-gate diagnostic only. It reports marker detection in all QC-passing nuclei and does not reassign clusters or relax identity criteria.",
  paste("QC rule: detected genes >=150 and valid CellBin polygon; donors:", length(unique(tab$donor_id))),
  "Training markers are retained for measurement; independent holdout sets exclude G1/G5 overlap.",
  "A marker below 1% or missing from a donor catalog is recorded as measured sparsity/coverage, not as proof that the cell lineage is absent.",
  "No diagnosis-linked measurements or effects were calculated."
), file.path(out_dir, "marker_availability_report.md"))
message("Marker availability audit complete; no labels were changed.")
