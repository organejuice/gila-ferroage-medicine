#!/usr/bin/env Rscript

# Diagnosis-blind audit of the existing G6 RCTD fit. Never loads diagnosis.
suppressPackageStartupMessages(library(Matrix))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: 78_diagnose_g6_mg_support_effect_blind.R <project_root> <output_dir>")
root <- args[[1L]]
if (!dir.exists(root)) stop("Project root does not exist")
outdir <- args[[2L]]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
fitdir <- file.path(root, "data/interim/G6_GSE269906_RCTD")
countdir <- file.path(root, "data/interim/G6_GSE269906_bin50_full")
auditdir <- file.path(root, "results/2026-09-24-disease-side-v3.5-g6-input-audit/registries")
gsm <- sub("_RCTD_weights_effect_blind.rds$", "", list.files(fitdir, pattern = "^GSM[0-9]+_RCTD_weights_effect_blind[.]rds$"))
if (length(gsm) != 12L || anyDuplicated(gsm)) stop("Expected exactly 12 distinct RCTD donor files")
input_paths <- c(file.path(fitdir, paste0(gsm, "_RCTD_weights_effect_blind.rds")),
                 file.path(countdir, paste0(gsm, "_bin50_centroid_supported_counts.rds")),
                 file.path(auditdir, "G1_microglia_signed_sender_244.tsv.gz"),
                 file.path(auditdir, "G5_receiver_one_to_one_human_ortholog_registry.tsv.gz"))
if (!all(file.exists(input_paths))) stop("Missing diagnostic input")
input_manifest <- data.frame(path = input_paths, bytes = file.info(input_paths)$size,
                             md5 = unname(tools::md5sum(input_paths)), stringsAsFactors = FALSE)
if (anyNA(input_manifest$md5)) stop("Could not hash every input")
write.csv(input_manifest, file.path(outdir, "input_manifest.csv"), row.names = FALSE)
g1 <- read.delim(gzfile(file.path(auditdir, "G1_microglia_signed_sender_244.tsv.gz")), stringsAsFactors = FALSE)
g1genes <- unique(toupper(g1$gene))
if (length(g1genes) != 244L) stop("Unexpected G1 signature length")
g5 <- read.delim(gzfile(file.path(auditdir, "G5_receiver_one_to_one_human_ortholog_registry.tsv.gz")), stringsAsFactors = FALSE)
g5genes <- unique(toupper(g5$human_one_to_one_symbol))
g5genes <- g5genes[!is.na(g5genes) & nzchar(g5genes)]
markers <- c("AIF1", "C1QA", "C1QB", "C1QC", "CX3CR1", "LST1", "P2RY12", "TMEM119", "TREM2", "TYROBP")
asc_markers <- c("ALDH1L1", "GFAP", "GJA1", "GLUL", "S100B", "SLC1A2", "SLC1A3", "SOX9", "CD74", "AQP4")
excluded <- unique(c(g1genes, g5genes, markers, asc_markers))
thresholds <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.7)
donor_rows <- vector("list", length(gsm))
support_rows <- vector("list", length(gsm) * length(thresholds))
gene_rows <- vector("list", length(gsm) * length(thresholds))
marker_rows <- vector("list", length(gsm) * 2L)
marker_gene_rows <- vector("list", length(gsm) * 2L)
class_rows <- vector("list", length(gsm))
class_gene_rows <- vector("list", length(gsm))
package_rows <- vector("list", length(gsm))
package_gene_rows <- vector("list", length(gsm))
row_i <- 0L
marker_i <- 0L

for (k in seq_along(gsm)) {
  fit <- readRDS(file.path(fitdir, paste0(gsm[[k]], "_RCTD_weights_effect_blind.rds")))
  spatial <- readRDS(file.path(countdir, paste0(gsm[[k]], "_bin50_centroid_supported_counts.rds")))
  w <- fit$weights
  if (!identical(as.character(w$tile_id), colnames(spatial$counts)) ||
      !identical(rownames(fit$weights_raw), colnames(spatial$counts))) {
    stop("Bin order mismatch for ", gsm[[k]])
  }
  usable <- w$spot_class != "reject" & w$raw_weight_sum > 0
  if (!all(is.finite(w$MG[usable])) || !all(is.finite(w$total_UMI))) stop("Nonfinite input")
  if (max(abs(rowSums(fit$weights_raw[usable, , drop = FALSE]) - w$raw_weight_sum[usable])) > 1e-8 ||
      max(abs(fit$weights_raw[usable, "MG"] / w$raw_weight_sum[usable] - w$MG[usable])) > 1e-8) {
    stop("Raw/normalized MG weight mismatch for ", gsm[[k]])
  }
  mg <- w$MG[usable]
  x <- spatial$counts[, usable, drop = FALSE]
  umi <- w$total_UMI[usable]
  frac <- Matrix::colSums(x) / umi
  retained <- setdiff(rownames(x), excluded)
  retained_frac <- Matrix::colSums(x[retained, , drop = FALSE]) / umi
  marker_genes <- intersect(markers, rownames(x))
  if (length(marker_genes) < 3L) stop("Too few MG markers for ", gsm[[k]])
  marker_counts <- x[marker_genes, , drop = FALSE]
  marker_n <- Matrix::colSums(marker_counts > 0)
  marker_cpm <- t(t(marker_counts) / umi) * 1e6
  marker_score <- Matrix::colMeans(log2(as.matrix(marker_cpm) + 1))
  top_marker <- marker_score >= as.numeric(quantile(marker_score, 0.99, names = FALSE))
  top_mg <- mg >= as.numeric(quantile(mg, 0.99, names = FALSE))
  donor_rows[[k]] <- data.frame(
    GSM = gsm[[k]], donor_id = as.character(w$donor_id[[1L]]),
    n_bins = nrow(w), n_usable = sum(usable), n_rejected = sum(!usable),
    n_markers = length(marker_genes),
    median_umi = median(umi), median_detected_genes = median(w$n_features[usable]),
    median_common_gene_umi_fraction = median(frac),
    median_fit_input_umi_fraction = median(retained_frac),
    n_fit_input_features = length(retained),
    median_raw_weight_sum = median(w$raw_weight_sum[usable]),
    max_MG_weight = max(mg), q99_MG_weight = as.numeric(quantile(mg, 0.99)),
    median_MG_weight = median(mg),
    fraction_MG_weight_positive = mean(mg > 0),
    fraction_two_or_more_MG_markers = mean(marker_n >= 2),
    fraction_three_or_more_MG_markers = mean(marker_n >= 3),
    median_MG_weight_top_marker_1pct = median(mg[top_marker]),
    max_MG_weight_top_marker_1pct = max(mg[top_marker]),
    median_marker_score_top_MG_1pct = median(marker_score[top_mg]),
    median_marker_score_all_bins = median(marker_score),
    n_first_type_MG = sum(w$first_type[usable] == "MG", na.rm = TRUE),
    n_second_type_MG = sum(w$second_type[usable] == "MG", na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  g1present <- intersect(g1genes, rownames(x))
  confident_mg <- w$first_type[usable] == "MG" |
    (w$spot_class[usable] == "doublet_certain" & w$second_type[usable] == "MG")
  confident_asc <- w$first_type[usable] == "ASC" |
    (w$spot_class[usable] == "doublet_certain" & w$second_type[usable] == "ASC")
  if (anyNA(confident_mg) || anyNA(confident_asc)) stop("Missing doublet classifications")
  n_confident_mg <- sum(confident_mg)
  class_detection <- setNames(rep(NA_real_, length(g1genes)), g1genes)
  if (n_confident_mg > 0L) {
    class_detection[g1present] <- Matrix::rowMeans(x[g1present, confident_mg, drop = FALSE] > 0)
    class_detection[setdiff(g1genes, g1present)] <- 0
  }
  class_rows[[k]] <- data.frame(
    GSM = gsm[[k]], donor_id = as.character(w$donor_id[[1L]]),
    n_confident_MG_by_doublet_class = n_confident_mg,
    n_confident_ASC_by_doublet_class = sum(confident_asc),
    n_G1_genes_detected_1pct_in_confident_MG = if (n_confident_mg) sum(class_detection >= 0.01) else NA_integer_,
    n_confident_MG_with_full_weight_ge_0_5 = sum(confident_mg & mg >= 0.5),
    n_confident_MG_with_two_or_more_holdout_markers = sum(confident_mg & marker_n >= 2),
    stringsAsFactors = FALSE
  )
  class_gene_rows[[k]] <- data.frame(GSM = gsm[[k]], gene = g1genes,
                                     detection_fraction = as.numeric(class_detection))
  mg_singlet <- w$spot_class[usable] == "singlet" & w$first_type[usable] == "MG"
  mg_certain_doublet <- w$spot_class[usable] == "doublet_certain" &
    (w$first_type[usable] == "MG" | w$second_type[usable] == "MG")
  package_mg <- mg_singlet | mg_certain_doublet
  package_detection <- setNames(rep(NA_real_, length(g1genes)), g1genes)
  if (sum(package_mg) > 0L) {
    package_detection[g1present] <- Matrix::rowMeans(x[g1present, package_mg, drop = FALSE] > 0)
    package_detection[setdiff(g1genes, g1present)] <- 0
  }
  package_rows[[k]] <- data.frame(
    GSM = gsm[[k]], donor_id = as.character(w$donor_id[[1L]]),
    n_MG_singlet = sum(mg_singlet),
    n_MG_certain_doublet = sum(mg_certain_doublet),
    n_MG_bin_nonzero_under_package_get_doublet_weights = sum(package_mg),
    n_G1_genes_detected_1pct = if (sum(package_mg)) sum(package_detection >= 0.01) else NA_integer_,
    stringsAsFactors = FALSE
  )
  package_gene_rows[[k]] <- data.frame(GSM = gsm[[k]], gene = g1genes,
                                       detection_fraction = as.numeric(package_detection))
  for (t in thresholds) {
    row_i <- row_i + 1L
    selected <- if (t == 0) rep(TRUE, length(mg)) else mg >= t
    n <- sum(selected)
    detection <- setNames(rep(NA_real_, length(g1genes)), g1genes)
    if (n > 0L) {
      detection[g1present] <- Matrix::rowMeans(x[g1present, selected, drop = FALSE] > 0)
      detection[setdiff(g1genes, g1present)] <- 0
    }
    support_rows[[row_i]] <- data.frame(
      GSM = gsm[[k]], donor_id = as.character(w$donor_id[[1L]]),
      MG_weight_threshold = t, n_bins = n,
      n_G1_genes_detected_1pct = if (n) sum(detection >= 0.01) else NA_integer_,
      n_bins_two_or_more_MG_markers = if (n) sum(marker_n[selected] >= 2) else 0L,
      stringsAsFactors = FALSE
    )
    gene_rows[[row_i]] <- data.frame(GSM = gsm[[k]], MG_weight_threshold = t,
                                    gene = g1genes, detection_fraction = as.numeric(detection),
                                    stringsAsFactors = FALSE)
  }
  for (m in c(2L, 3L)) {
    marker_i <- marker_i + 1L
    selected <- marker_n >= m
    n <- sum(selected)
    detection <- setNames(rep(NA_real_, length(g1genes)), g1genes)
    if (n > 0L) {
      detection[g1present] <- Matrix::rowMeans(x[g1present, selected, drop = FALSE] > 0)
      detection[setdiff(g1genes, g1present)] <- 0
    }
    marker_rows[[marker_i]] <- data.frame(
      GSM = gsm[[k]], donor_id = as.character(w$donor_id[[1L]]),
      minimum_detected_MG_markers = m, n_bins = n,
      n_G1_genes_detected_1pct = if (n) sum(detection >= 0.01) else NA_integer_,
      max_MG_weight = if (n) max(mg[selected]) else NA_real_,
      median_MG_weight = if (n) median(mg[selected]) else NA_real_,
      stringsAsFactors = FALSE
    )
    marker_gene_rows[[marker_i]] <- data.frame(GSM = gsm[[k]],
                                              minimum_detected_MG_markers = m,
                                              gene = g1genes,
                                              detection_fraction = as.numeric(detection))
  }
  message("Audited ", gsm[[k]], ": usable=", sum(usable), ", max MG=", round(max(mg), 4))
  rm(fit, spatial, x, marker_counts, marker_cpm)
  gc(verbose = FALSE)
}
donor <- do.call(rbind, donor_rows)
support <- do.call(rbind, support_rows)
genes <- do.call(rbind, gene_rows)
marker_support <- do.call(rbind, marker_rows)
marker_genes <- do.call(rbind, marker_gene_rows)
class_support <- do.call(rbind, class_rows)
class_genes <- do.call(rbind, class_gene_rows)
package_support <- do.call(rbind, package_rows)
package_genes <- do.call(rbind, package_gene_rows)
overall <- do.call(rbind, lapply(thresholds, function(t) {
  s <- support[support$MG_weight_threshold == t, , drop = FALSE]
  gg <- genes[genes$MG_weight_threshold == t, , drop = FALSE]
  per_gene <- aggregate(detection_fraction ~ gene, gg, function(z) sum(z >= 0.01, na.rm = TRUE))
  data.frame(MG_weight_threshold = t, n_donors_with_bins = sum(s$n_bins > 0),
             n_bins_total = sum(s$n_bins), n_G1_genes_supported_in_10_donors = sum(per_gene$detection_fraction >= 10),
             n_G1_genes_supported_in_all_12 = sum(per_gene$detection_fraction == 12))
}))
write.csv(donor, file.path(outdir, "G6_MG_donor_diagnostic.csv"), row.names = FALSE)
write.csv(support, file.path(outdir, "G6_MG_threshold_support.csv"), row.names = FALSE)
write.csv(overall, file.path(outdir, "G6_G1_threshold_feasibility.csv"), row.names = FALSE)
write.csv(genes, file.path(outdir, "G6_G1_gene_detection_by_donor_threshold.csv"), row.names = FALSE)
write.csv(marker_support, file.path(outdir, "G6_MG_marker_positive_support.csv"), row.names = FALSE)
write.csv(marker_genes, file.path(outdir, "G6_G1_gene_detection_by_marker_support.csv"), row.names = FALSE)
write.csv(class_support, file.path(outdir, "G6_MG_doublet_class_support.csv"), row.names = FALSE)
write.csv(class_genes, file.path(outdir, "G6_G1_gene_detection_by_doublet_class.csv"), row.names = FALSE)
write.csv(package_support, file.path(outdir, "G6_MG_package_doublet_nonzero_support.csv"), row.names = FALSE)
write.csv(package_genes, file.path(outdir, "G6_G1_gene_detection_by_package_doublet_nonzero.csv"), row.names = FALSE)
writeLines(c("G6 diagnosis-blind RCTD support diagnostic", "No diagnosis labels or group effects loaded.",
             paste("R version:", R.version.string),
             paste("Script MD5:", unname(tools::md5sum(file.path(root, "src/78_diagnose_g6_mg_support_effect_blind.R"))))),
           file.path(outdir, "provenance.txt"))
