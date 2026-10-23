#!/usr/bin/env Rscript

# Diagnosis-blind audit for G6-FP-2026-09-27-001, WP1.
suppressPackageStartupMessages(library(Matrix))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: 90_audit_g6_numi_sensitivity.R <project_root>")
root <- args[[1L]]
gsms <- c("GSM8330063", "GSM8330064")
outdir <- file.path(root, "results/2026-09-27-g6-next-feasibility/numi")
fitdir <- file.path(root, "data/interim/G6_next_feasibility_001/numi")
legacy_fitdir <- file.path(root, "data/interim/G6_GSE269906_RCTD_dual_output_pilot")
olddir <- file.path(root, "data/interim/G6_GSE269906_RCTD")
countdir <- file.path(root, "data/interim/G6_GSE269906_bin50_full")
g1path <- file.path(root, "results/2026-09-24-disease-side-v3.5-g6-input-audit/registries/G1_microglia_signed_sender_244.tsv.gz")
g5path <- file.path(root, "results/2026-09-24-disease-side-v3.5-g6-input-audit/registries/G5_receiver_one_to_one_human_ortholog_registry.tsv.gz")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
if (!all(file.exists(c(g1path, g5path)))) stop("Frozen G1/G5 registry missing")
g1 <- read.delim(gzfile(g1path), check.names = FALSE, stringsAsFactors = FALSE)
if (nrow(g1) != 244L || anyDuplicated(g1$gene)) stop("Frozen G1 registry is not 244 unique genes")
g1genes <- unique(toupper(g1$gene))
g5 <- read.delim(gzfile(g5path), check.names = FALSE, stringsAsFactors = FALSE)
state_genes <- unique(toupper(c(g1$gene, g5$human_one_to_one_symbol)))
state_genes <- state_genes[!is.na(state_genes) & nzchar(state_genes)]
mg_markers <- setdiff(c("AIF1", "C1QA", "C1QB", "C1QC", "CX3CR1", "LST1", "P2RY12", "TMEM119", "TREM2", "TYROBP"), state_genes)
asc_markers <- setdiff(c("ALDH1L1", "GFAP", "GJA1", "S100B", "SLC1A2", "SLC1A3", "SOX9"), state_genes)
if (length(mg_markers) < 3L || length(asc_markers) < 3L) stop("Insufficient independent holdout markers")
marker_score <- function(counts, genes, umi) {
  genes <- intersect(genes, rownames(counts))
  if (length(genes) < 3L) stop("Fewer than three measured holdout markers")
  mat <- counts[genes, , drop = FALSE]
  cpm <- t(t(mat) / pmax(umi, 1)) * 1e6
  list(score = Matrix::colMeans(log2(as.matrix(cpm) + 1)), positive = Matrix::colSums(mat > 0), n_genes = length(genes))
}

comparison <- list()
gene_rows <- list()
integrity <- list()
for (gsm in gsms) {
  spatial <- readRDS(file.path(countdir, paste0(gsm, "_bin50_centroid_supported_counts.rds")))
  ids <- colnames(spatial$counts)
  old <- readRDS(file.path(olddir, paste0(gsm, "_RCTD_weights_effect_blind.rds")))
  baseline_path <- file.path(legacy_fitdir, paste0(gsm, "_dual_output_pilot.rds"))
  baseline <- readRDS(baseline_path)
  if (!identical(rownames(baseline$all_weights), ids) || !identical(rownames(baseline$results_df), ids) ||
      !identical(rownames(old$weights_raw), ids) || !identical(as.character(old$weights$tile_id), ids)) stop("Baseline barcode order mismatch: ", gsm)
  old_class_match <- all(vapply(c("spot_class", "first_type", "second_type"), function(nm)
    identical(as.character(baseline$results_df[[nm]]), as.character(old$weights[[nm]])), logical(1)))
  old_diff <- max(abs(baseline$all_weights - as.matrix(old$weights_raw)))
  baseline_ok <- old_class_match && is.finite(old_diff) && old_diff <= 1e-6
  candidate_modes <- list()
  new_full_path <- file.path(fitdir, paste0(gsm, "_full.rds"))
  new_subset_path <- file.path(fitdir, paste0(gsm, "_fit_features.rds"))
  if (file.exists(new_full_path)) candidate_modes$full <- readRDS(new_full_path) else candidate_modes$full <- baseline
  if (!file.exists(new_subset_path)) stop("Missing fitted-feature nUMI sensitivity output: ", new_subset_path)
  candidate_modes$fit_features <- readRDS(new_subset_path)
  fit_features <- candidate_modes$fit_features$metadata$fit_features
  if (is.null(fit_features) || length(fit_features) < 10000L) stop("Cannot verify fitted-feature nUMI definition")
  fit_feature_sum <- as.numeric(Matrix::colSums(spatial$counts[fit_features, , drop = FALSE]))
  for (mode in names(candidate_modes)) {
    fit <- candidate_modes[[mode]]
    if (!identical(rownames(fit$all_weights), ids) || !identical(rownames(fit$doublet_weights), ids) ||
        !identical(rownames(fit$mapped_doublet_weights), ids) || !identical(rownames(fit$results_df), ids)) stop("Sensitivity output barcode mismatch: ", gsm, "/", mode)
    if (!identical(fit$donor_id, as.character(old$weights$donor_id[[1L]]))) stop("Donor ID mismatch: ", gsm)
    if (mode == "full") {
      fit_reproduced <- all(vapply(c("spot_class", "first_type", "second_type"), function(nm)
        identical(as.character(fit$results_df[[nm]]), as.character(old$weights[[nm]])), logical(1))) &&
        max(abs(fit$all_weights - as.matrix(old$weights_raw))) <= 1e-6
    } else fit_reproduced <- NA
    class <- as.character(fit$results_df$spot_class)
    valid_class <- class %in% c("singlet", "doublet_certain")
    if (anyNA(class) || !all(valid_class | class %in% c("reject", "doublet_uncertain"))) stop("Invalid RCTD class labels")
    all_raw <- fit$all_weights
    all_sum <- rowSums(all_raw)
    all_normalized <- matrix(0, nrow(all_raw), ncol(all_raw), dimnames = dimnames(all_raw))
    usable_all <- class != "reject" & all_sum > 0
    all_normalized[usable_all, ] <- all_raw[usable_all, , drop = FALSE] / all_sum[usable_all]
    mapped <- fit$mapped_doublet_weights
    mapped_sum <- rowSums(mapped)
    estimators <- list(all_type_normalized = list(w = all_normalized, usable = usable_all),
                       doublet_mapped = list(w = mapped, usable = valid_class & mapped_sum > 0))
    umi <- if (!is.null(fit$nUMI)) as.numeric(fit$nUMI) else if (mode == "full") as.numeric(spatial$metadata$total_UMI) else numeric()
    if (length(umi) != length(ids) || any(!is.finite(umi)) || any(umi <= 0)) stop("Invalid saved nUMI vector")
    mgmark <- marker_score(spatial$counts, mg_markers, umi)
    ascmark <- marker_score(spatial$counts, asc_markers, umi)
    present_g1 <- intersect(g1genes, rownames(spatial$counts))
    for (estimator in names(estimators)) {
      w <- estimators[[estimator]]$w
      usable <- estimators[[estimator]]$usable
      mg <- w[, "MG"]
      selected <- usable & mg >= 0.5
      nsel <- sum(selected)
      detection <- setNames(rep(0, length(g1genes)), g1genes)
      if (nsel) detection[present_g1] <- Matrix::rowMeans(spatial$counts[present_g1, selected, drop = FALSE] > 0)
      gene_rows[[paste(gsm, mode, estimator)]] <- data.frame(GSM = gsm, donor_id = fit$donor_id, nUMI_mode = mode,
        estimator = estimator, gene = g1genes, detection_fraction = unname(detection), meets_1pct = unname(detection) >= 0.01)
      rho <- if (sum(usable) >= 3L && sd(mg[usable]) > 0) suppressWarnings(cor(mg[usable], mgmark$score[usable], method = "spearman")) else NA_real_
      selected_median <- if (nsel) median(mgmark$score[selected]) else NA_real_
      other <- usable & !selected
      other_median <- if (sum(other)) median(mgmark$score[other]) else NA_real_
      comparison[[paste(gsm, mode, estimator)]] <- data.frame(
        GSM = gsm, donor_id = fit$donor_id, nUMI_mode = mode, estimator = estimator,
        n_usable_bins = sum(usable), n_MG_ge_0_5 = nsel,
        n_G1_genes_detected_at_least_1pct = sum(detection >= 0.01),
        n_MG_bins_with_two_holdout_markers = if (nsel) sum(mgmark$positive[selected] >= 2L) else 0L,
        MG_holdout_marker_rho = rho, MG_marker_median_selected = selected_median,
        MG_marker_median_other = other_median,
        MG_selected_marker_gt_other = is.finite(selected_median) && is.finite(other_median) && selected_median > other_median,
        baseline_reproduced = if (mode == "full") fit_reproduced else baseline_ok,
        baseline_full_weight_max_abs_diff = if (mode == "full") max(abs(fit$all_weights - as.matrix(old$weights_raw))) else old_diff,
        reference_status = "REFERENCE_IDENTITY_UNRESOLVED_DIAGNOSTIC_ONLY", stringsAsFactors = FALSE)
    }
    integrity[[paste(gsm, mode)]] <- data.frame(GSM = gsm, donor_id = fit$donor_id, nUMI_mode = mode,
      baseline_reproduced = if (mode == "full") fit_reproduced else baseline_ok,
      full_weight_max_abs_diff_vs_locked = if (mode == "full") max(abs(fit$all_weights - as.matrix(old$weights_raw))) else old_diff,
      class_order_exact = if (mode == "full") fit_reproduced else baseline_ok,
      nUMI_matches_mode_definition = if (mode == "full") max(abs(umi - spatial$metadata$total_UMI)) < 1e-8 else
        max(abs(umi - fit_feature_sum)) < 1e-8,
      nUMI_vs_fit_matrix_sum_max_abs_diff = max(abs(umi - fit_feature_sum)),
      stringsAsFactors = FALSE)
  }
  rm(spatial, old, baseline)
  gc(verbose = FALSE)
}
summary <- do.call(rbind, comparison)
integrity_table <- do.call(rbind, integrity)
genes <- do.call(rbind, gene_rows)
write.table(summary, file.path(outdir, "donor_estimator_comparison.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")
write.table(integrity_table, file.path(outdir, "fit_integrity.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")
write.table(genes, file.path(outdir, "G1_gene_detection_by_donor_mode_estimator.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")
summary$gate_pass <- summary$n_MG_ge_0_5 >= 40L & summary$n_G1_genes_detected_at_least_1pct >= 61L &
  is.finite(summary$MG_holdout_marker_rho) & summary$MG_holdout_marker_rho > 0 & summary$MG_selected_marker_gt_other
write.table(summary[, c("GSM", "donor_id", "nUMI_mode", "estimator", "n_MG_ge_0_5", "n_G1_genes_detected_at_least_1pct",
                        "MG_holdout_marker_rho", "MG_marker_median_selected", "MG_marker_median_other", "MG_selected_marker_gt_other",
                        "baseline_reproduced", "gate_pass", "reference_status")],
            file.path(outdir, "mechanical_gate_by_donor.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")
eligible <- any(vapply(split(summary, interaction(summary$nUMI_mode, summary$estimator, drop = TRUE)),
                       function(z) nrow(z) == 2L && all(z$gate_pass) && all(z$baseline_reproduced), logical(1)))
status <- if (any(!integrity_table$baseline_reproduced | !integrity_table$nUMI_matches_mode_definition)) "PILOT_INVALID" else
  if (eligible) "RECOMMEND_FULL_COHORT_METHOD_SENSITIVITY_REVIEW" else "STOP_NO_FULL_COHORT_REFIT"
writeLines(c(
  "# G6 nUMI sensitivity: diagnosis-blind feasibility result", "",
  paste("Status:", status),
  "No diagnosis labels or AD-control effects were read or computed.",
  "The full-library baseline is the pre-existing dual-output pilot, re-audited against the locked legacy fit (exact class/order and max full-weight difference <= 1e-6); the new fit_features mode was fitted under the same locked model settings.",
  paste("Fitted-feature outputs:", paste(file.exists(file.path(fitdir, paste0(gsms, "_fit_features.rds"))), collapse = ", ")),
  "Reference identity/independence remains unresolved; this result is diagnostic only.", "", "## Donor by mode and estimator", "",
  paste(capture.output(print(summary[, c("GSM", "nUMI_mode", "estimator", "n_MG_ge_0_5", "n_G1_genes_detected_at_least_1pct",
                                             "MG_holdout_marker_rho", "MG_selected_marker_gt_other", "gate_pass")], row.names = FALSE)), collapse = "\n"),
  "", "## Interpretation limit", "",
  "Passing only recommends review of a locked full-cohort model sensitivity run. It does not authorize diagnosis effects, a new G6 protocol, or G6_PASS.",
  paste("R version:", R.version.string), paste("Matrix:", as.character(packageVersion("Matrix")))
), file.path(outdir, "report.md"))
writeLines(capture.output(sessionInfo()), file.path(outdir, "R_sessionInfo.txt"))
message("WP1 status: ", status)
if (status == "PILOT_INVALID") quit(status = 2L)
