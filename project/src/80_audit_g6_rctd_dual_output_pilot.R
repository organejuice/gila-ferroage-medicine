#!/usr/bin/env Rscript

# G6-DIAG-2026-09-26-001: audit only; never reads diagnosis labels.
suppressPackageStartupMessages({ library(Matrix); library(spacexr) })
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: 80_audit_g6_rctd_dual_output_pilot.R <project_root>")
root <- args[[1L]]
if (!dir.exists(root)) stop("Project root does not exist")
gsms <- c("GSM8330063", "GSM8330064")
outdir <- file.path(root, "results/2026-09-26-g6-rctd-dual-output-pilot")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
fitdir <- file.path(root, "data/interim/G6_GSE269906_RCTD_dual_output_pilot")
olddir <- file.path(root, "data/interim/G6_GSE269906_RCTD")
countdir <- file.path(root, "data/interim/G6_GSE269906_bin50_full")
g1path <- file.path(root, "results/2026-09-24-disease-side-v3.5-g6-input-audit/registries/G1_microglia_signed_sender_244.tsv.gz")
paths <- c(file.path(fitdir, paste0(gsms, "_dual_output_pilot.rds")),
           file.path(olddir, paste0(gsms, "_RCTD_weights_effect_blind.rds")),
           file.path(countdir, paste0(gsms, "_bin50_centroid_supported_counts.rds")),
           g1path, file.path(root, "src/80_audit_g6_rctd_dual_output_pilot.R"),
           file.path(root, "configs/g6_rctd_dual_output_diagnostic_001.yaml"))
if (!all(file.exists(paths))) stop("Pilot output or locked audit input missing")
manifest <- data.frame(path = paths, bytes = file.info(paths)$size,
                       md5 = unname(tools::md5sum(paths)), stringsAsFactors = FALSE)
if (anyNA(manifest$md5)) stop("Cannot hash all inputs")
write.csv(manifest, file.path(outdir, "input_manifest.csv"), row.names = FALSE)

g1 <- read.delim(gzfile(g1path), stringsAsFactors = FALSE)
g1genes <- unique(toupper(g1$gene))
if (length(g1genes) != 244L) stop("G1 signature length mismatch")
mg_markers <- c("AIF1", "C1QA", "C1QB", "C1QC", "CX3CR1", "LST1", "P2RY12", "TMEM119", "TREM2", "TYROBP")
asc_markers <- c("ALDH1L1", "GFAP", "GJA1", "GLUL", "S100B", "SLC1A2", "SLC1A3", "SOX9", "CD74", "AQP4")
summary_rows <- list()
gene_rows <- list()
integrity_rows <- list()

marker_score <- function(counts, genes, umi) {
  genes <- intersect(genes, rownames(counts))
  if (length(genes) < 3L) stop("Too few holdout markers")
  mat <- counts[genes, , drop = FALSE]
  cpm <- t(t(mat) / pmax(umi, 1)) * 1e6
  list(score = Matrix::colMeans(log2(as.matrix(cpm) + 1)),
       positive = Matrix::colSums(mat > 0), n_genes = length(genes))
}

for (gsm in gsms) {
  fit <- readRDS(file.path(fitdir, paste0(gsm, "_dual_output_pilot.rds")))
  old <- readRDS(file.path(olddir, paste0(gsm, "_RCTD_weights_effect_blind.rds")))
  spatial <- readRDS(file.path(countdir, paste0(gsm, "_bin50_centroid_supported_counts.rds")))
  ids <- colnames(spatial$counts)
  if (!identical(rownames(fit$all_weights), ids) ||
      !identical(rownames(fit$doublet_weights), ids) ||
      !identical(rownames(fit$mapped_doublet_weights), ids) ||
      !identical(rownames(fit$results_df), ids) ||
      !identical(rownames(old$weights_raw), ids) ||
      !identical(as.character(old$weights$tile_id), ids)) stop("Barcode order mismatch: ", gsm)
  if (!identical(fit$donor_id, as.character(old$weights$donor_id[[1L]]))) stop("Donor ID mismatch: ", gsm)

  class_columns <- c("spot_class", "first_type", "second_type")
  class_match <- vapply(class_columns, function(nm) {
    identical(as.character(fit$results_df[[nm]]), as.character(old$weights[[nm]]))
  }, logical(1))
  full_diff <- max(abs(fit$all_weights - as.matrix(old$weights_raw)))
  full_reproduced <- all(class_match) && is.finite(full_diff) && full_diff <= 1e-6

  all_raw <- fit$all_weights
  doublet_raw <- fit$doublet_weights
  mapped <- fit$mapped_doublet_weights
  class <- as.character(fit$results_df$spot_class)
  valid_class <- class %in% c("singlet", "doublet_certain")
  invalid_class <- class %in% c("reject", "doublet_uncertain")
  if (anyNA(class) || !all(valid_class | invalid_class)) stop("Unrecognized spot class")
  all_sum <- rowSums(all_raw)
  mapped_sum <- rowSums(mapped)
  raw_doublet_sum <- rowSums(doublet_raw)
  finite_nonnegative <- all(is.finite(all_raw)) && all(all_raw >= -1e-8) &&
    all(is.finite(doublet_raw)) && all(doublet_raw >= -1e-8) &&
    all(is.finite(mapped)) && all(mapped >= -1e-8)
  mapped_nonzero_types <- rowSums(mapped > 1e-8)
  n_uncertain_nonzero <- sum(invalid_class & mapped_sum > 1e-8)
  n_singlet_bad_sum <- sum(class == "singlet" & abs(mapped_sum - 1) > 1e-8)
  n_certain_bad_types <- sum(class == "doublet_certain" & mapped_nonzero_types > 2L)
  integrity_rows[[gsm]] <- data.frame(
    GSM = gsm, donor_id = fit$donor_id, n_bins = length(ids),
    n_singlet = sum(class == "singlet"), n_doublet_certain = sum(class == "doublet_certain"),
    n_doublet_uncertain = sum(class == "doublet_uncertain"), n_reject = sum(class == "reject"),
    class_reproduced = all(class_match), full_weight_max_abs_diff = full_diff,
    full_fit_reproduced = full_reproduced, finite_nonnegative = finite_nonnegative,
    n_invalid_class_mapped_nonzero = n_uncertain_nonzero,
    n_singlet_mapped_sum_not_one = n_singlet_bad_sum,
    n_certain_more_than_two_types = n_certain_bad_types,
    raw_doublet_sum_min_valid = min(raw_doublet_sum[valid_class]),
    raw_doublet_sum_median_valid = median(raw_doublet_sum[valid_class]),
    raw_doublet_sum_max_valid = max(raw_doublet_sum[valid_class]),
    mapped_sum_min_valid = min(mapped_sum[valid_class]),
    mapped_sum_median_valid = median(mapped_sum[valid_class]),
    mapped_sum_max_valid = max(mapped_sum[valid_class]), stringsAsFactors = FALSE)

  usable_all <- class != "reject" & all_sum > 0
  all_normalized <- matrix(0, nrow(all_raw), ncol(all_raw), dimnames = dimnames(all_raw))
  all_normalized[usable_all, ] <- all_raw[usable_all, , drop = FALSE] / all_sum[usable_all]
  if (max(abs(all_normalized[, "MG"] - old$weights$MG)) > 1e-6) stop("Original full-weight normalization mismatch")
  umi <- as.numeric(spatial$metadata$total_UMI)
  mgmark <- marker_score(spatial$counts, mg_markers, umi)
  ascmark <- marker_score(spatial$counts, asc_markers, umi)
  gene_present <- intersect(g1genes, rownames(spatial$counts))
  for (estimator in c("all_type_normalized", "doublet_mapped")) {
    weights <- if (estimator == "all_type_normalized") all_normalized else mapped
    usable <- if (estimator == "all_type_normalized") usable_all else valid_class & mapped_sum > 0
    mg <- weights[, "MG"]
    selected <- usable & mg >= 0.5
    n <- sum(selected)
    detection <- setNames(rep(0, length(g1genes)), g1genes)
    if (n > 0L) detection[gene_present] <- Matrix::rowMeans(spatial$counts[gene_present, selected, drop = FALSE] > 0)
    gene_rows[[paste(gsm, estimator)]] <- data.frame(
      GSM = gsm, donor_id = fit$donor_id, estimator = estimator,
      gene = g1genes, detection_fraction = unname(detection),
      meets_1pct = unname(detection) >= 0.01, stringsAsFactors = FALSE)
    rho <- if (sum(usable) >= 3L && sd(mg[usable]) > 0) {
      suppressWarnings(cor(mg[usable], mgmark$score[usable], method = "spearman"))
    } else NA_real_
    summary_rows[[paste(gsm, estimator)]] <- data.frame(
      GSM = gsm, donor_id = fit$donor_id, estimator = estimator,
      n_usable_bins = sum(usable), n_MG_ge_0_5 = n,
      n_G1_genes_detected_at_least_1pct = sum(detection >= 0.01),
      n_MG_bins_with_two_holdout_markers = if (n) sum(mgmark$positive[selected] >= 2) else 0L,
      fraction_MG_bins_with_two_holdout_markers = if (n) mean(mgmark$positive[selected] >= 2) else NA_real_,
      MG_weight_marker_spearman = rho,
      MG_marker_median_selected = if (n) median(mgmark$score[selected]) else NA_real_,
      MG_marker_median_other_usable = if (sum(usable & !selected)) median(mgmark$score[usable & !selected]) else NA_real_,
      ASC_marker_median_selected = if (n) median(ascmark$score[selected]) else NA_real_,
      n_selected_singlet = sum(selected & class == "singlet"),
      n_selected_doublet_certain = sum(selected & class == "doublet_certain"),
      n_selected_below_40 = n < 40L, stringsAsFactors = FALSE)
  }
  message("AUDIT ", gsm, "; max all-weight difference=", format(full_diff, digits = 4),
          "; class match=", all(class_match))
  rm(fit, old, spatial)
  gc(verbose = FALSE)
}
integrity <- do.call(rbind, integrity_rows)
summary <- do.call(rbind, summary_rows)
genes <- do.call(rbind, gene_rows)
write.csv(integrity, file.path(outdir, "fit_reproduction_and_weight_integrity.csv"), row.names = FALSE)
write.csv(summary, file.path(outdir, "donor_estimator_comparison.csv"), row.names = FALSE)
write.csv(genes, file.path(outdir, "G1_gene_detection_by_donor_estimator.csv"), row.names = FALSE)
gate_reproduced <- all(integrity$full_fit_reproduced)
gate_integrity <- all(integrity$finite_nonnegative) &&
  all(integrity$n_invalid_class_mapped_nonzero == 0L) &&
  all(integrity$n_singlet_mapped_sum_not_one == 0L) &&
  all(integrity$n_certain_more_than_two_types == 0L)
status <- if (!gate_reproduced) "PILOT_INVALID_FIT_NOT_REPRODUCED" else
  if (!gate_integrity) "PILOT_INVALID_WEIGHT_INTEGRITY" else
    "PILOT_VALID_DIAGNOSTIC_ONLY"
writeLines(c(
  "# G6 RCTD dual-output pilot: diagnosis-blind result",
  "", paste("Status:", status),
  "Frozen v3.5 G6_NOT_EXECUTABLE remains unchanged. No diagnosis fields or group effects were read.",
  "", "## Numeric comparison", "",
  paste(capture.output(print(summary[, c("GSM", "estimator", "n_MG_ge_0_5",
                                        "n_G1_genes_detected_at_least_1pct",
                                        "MG_weight_marker_spearman", "n_selected_below_40")], row.names = FALSE)), collapse = "\n"),
  "", "## Interpretation limit", "",
  "Two donors cannot establish 10-of-12 G1 coverage, identity specificity, donor-level AD effects, or G6_PASS.",
  paste("R version:", R.version.string),
  paste("spacexr:", as.character(packageVersion("spacexr")))
), file.path(outdir, "pilot_report.md"))
writeLines(capture.output(sessionInfo()), file.path(outdir, "R_sessionInfo.txt"))
message("Audit status: ", status)
if (status != "PILOT_VALID_DIAGNOSTIC_ONLY") quit(status = 2L)
