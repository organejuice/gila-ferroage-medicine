#!/usr/bin/env Rscript

# Locked G6 primary analysis. The script first runs effect-blind identity and
# signature-coverage gates. It reads donor diagnoses only after those gates pass.
suppressPackageStartupMessages({
  library(Matrix)
})
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: 73_run_g6_primary_analysis.R <project_root>")
root <- args[[1L]]
run_id <- "2026-09-24-disease-side-v3.5-g6-primary"
lock_id <- "G6-EXEC-2026-09-24-001"
audit <- file.path(root, "results/2026-09-24-disease-side-v3.5-g6-input-audit")
input_dir <- file.path(root, "data/interim/G6_GSE269906_bin50_full")
weight_dir <- file.path(root, "data/interim/G6_GSE269906_RCTD")
mask_dir <- file.path(root, "data/interim/G6_GSE269906_bin50_masks")
outdir <- file.path(root, "results/2026-09-24-disease-side-v3.5-g6-primary")
dir.create(file.path(outdir, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "figures"), recursive = TRUE, showWarnings = FALSE)

read_tsv_gz <- function(path) read.delim(gzfile(path), check.names = FALSE,
                                          stringsAsFactors = FALSE)
g1 <- read_tsv_gz(file.path(audit, "registries/G1_microglia_signed_sender_244.tsv.gz"))
g5 <- read_tsv_gz(file.path(audit, "registries/G5_receiver_one_to_one_human_ortholog_registry.tsv.gz"))
donor_manifest_path <- file.path(audit, "registries/GSE269906_donor_section_manifest.tsv")
donor_manifest_columns <- names(read.delim(donor_manifest_path, nrows = 0L, check.names = FALSE))
if (!all(c("GSM", "donor_id") %in% donor_manifest_columns)) stop("Primary donor manifest lacks effect-blind identity fields")
donor_manifest_classes <- ifelse(donor_manifest_columns %in% c("GSM", "donor_id"), NA_character_, "NULL")
donors_blind <- read.delim(donor_manifest_path, check.names = FALSE, stringsAsFactors = FALSE,
                           colClasses = donor_manifest_classes)
reference_overlap <- jsonlite::fromJSON(file.path(audit, "GSE269906_GSE174367_overlap_audit.json"), simplifyVector = FALSE)
reference_overlap_summary <- list(
  candidate_pairs = reference_overlap$demographic_diagnosis_candidate_pairs,
  exact_crosswalk_found = reference_overlap$donor_identity_crosswalk_found,
  pairs_matching_braak_and_plaque = reference_overlap$pairs_also_matching_braak_and_plaque_stage,
  interpretation = "possible reference subject overlap remains unresolved; locked reference retained"
)
if (nrow(g1) != 244L || anyDuplicated(g1$gene) || nrow(donors_blind) != 12L ||
    anyDuplicated(donors_blind$GSM)) stop("Frozen signatures or primary donor manifest failed integrity checks")
modules <- sort(unique(g5$module))
module_info <- unique(g5[, c("module", "G5_primary_status", "G5_module_gate",
                              "G5_module_logFC", "G5_module_FDR_BH_7",
                              "G5_module_direction_hypothesis")])
if (nrow(module_info) != length(modules) || length(modules) != 7L) stop("Expected seven uniquely defined frozen G5 modules")
module_info$expected_sign <- ifelse(module_info$G5_module_direction_hypothesis == "positive_model", 1,
                                    ifelse(module_info$G5_module_direction_hypothesis == "negative_model", -1, NA_real_))
g5 <- merge(g5, module_info[, c("module", "expected_sign")], by = "module", sort = FALSE)
g1$gene <- toupper(g1$gene)
g5$human_one_to_one_symbol <- toupper(g5$human_one_to_one_symbol)
g5 <- g5[!is.na(g5$human_one_to_one_symbol) & nzchar(g5$human_one_to_one_symbol), , drop = FALSE]
g1_sign <- setNames(ifelse(g1$direction == "positive_AD", 1, -1), g1$gene)
module_genes <- lapply(modules, function(m) unique(g5$human_one_to_one_symbol[g5$module == m]))
names(module_genes) <- modules
all_state_genes <- unique(c(names(g1_sign), unlist(module_genes, use.names = FALSE)))

weight_path <- function(gsm) file.path(weight_dir, paste0(gsm, "_RCTD_weights_effect_blind.rds"))
input_path <- function(gsm) file.path(input_dir, paste0(gsm, "_bin50_centroid_supported_counts.rds"))
get_aligned_weights <- function(gsm, metadata) {
  object <- readRDS(weight_path(gsm))
  w <- object$weights
  raw <- object$weights_raw
  required <- c("tile_id", "MG", "ASC", "spot_class", "weight_sum_normalized", "raw_weight_sum")
  if (!all(required %in% names(w)) || is.null(raw) || ncol(raw) != 7L) stop("Missing RCTD raw/normalized weight fields for ", gsm)
  idx <- match(metadata$tile_id, w$tile_id)
  if (anyNA(idx) || anyDuplicated(w$tile_id) || nrow(w) != nrow(metadata)) stop("RCTD bins do not map one-to-one for ", gsm)
  w <- w[idx, , drop = FALSE]
  raw_idx <- match(as.character(metadata$tile_id), rownames(raw))
  if (anyNA(raw_idx) || anyDuplicated(rownames(raw)) || nrow(raw) != nrow(metadata)) stop("Raw RCTD weights do not map one-to-one for ", gsm)
  raw <- raw[raw_idx, , drop = FALSE]
  if (any(!is.finite(raw)) || any(raw < -1e-8)) stop("Raw RCTD weights contain non-finite or negative values for ", gsm)
  raw_sum <- rowSums(raw)
  if (any(!is.finite(raw_sum)) || any(abs(raw_sum - w$raw_weight_sum) > 1e-8)) stop("Raw RCTD weight sums disagree with saved bin metadata for ", gsm)
  if (any(as.character(w$spot_class) == "reject" & (w$MG != 0 | w$ASC != 0))) stop("Rejected RCTD bins retain usable weights for ", gsm)
  if (any(as.character(w$spot_class) != "reject" & raw_sum <= 0)) stop("RCTD usable class and raw-weight support disagree for ", gsm)
  if (any(as.character(w$spot_class) != "reject" & abs(w$weight_sum_normalized - 1) > 1e-7)) stop("Normalized RCTD compositions fail for ", gsm)
  if (any(!is.finite(w$MG)) || any(!is.finite(w$ASC)) || any(w$MG < 0 | w$ASC < 0)) stop("Invalid RCTD lineage weights for ", gsm)
  sum_qc_path <- file.path(weight_dir, paste0(gsm, "_RCTD_weight_sum_qc.csv"))
  if (!file.exists(sum_qc_path)) stop("Missing RCTD raw-weight-sum QC record for ", gsm)
  sum_qc <- read.csv(sum_qc_path, stringsAsFactors = FALSE)
  if (nrow(sum_qc) != 1L || sum_qc$n_bins != nrow(metadata) ||
      sum_qc$n_rejected != sum(as.character(w$spot_class) == "reject") ||
      sum_qc$n_nonfinite != 0L || sum_qc$n_negative != 0L ||
      sum_qc$normalized_sum_min_usable < 1 - 1e-7 || sum_qc$normalized_sum_max_usable > 1 + 1e-7) {
    stop("Saved RCTD weight-sum QC table failed reconciliation for ", gsm)
  }
  list(object = object, weights = w, sum_qc = sum_qc)
}

# Effect-blind audit: all donor identity outputs, reference-fit success,
# disjoint marker count/correlation, signature feature coverage, and bin mapping.
marker_rows <- vector("list", nrow(donors_blind))
coverage_rows <- list()
donor_blind_qc <- vector("list", nrow(donors_blind))
weight_sum_qc_rows <- vector("list", nrow(donors_blind))
for (i in seq_len(nrow(donors_blind))) {
  gsm <- donors_blind$GSM[i]
  if (!file.exists(weight_path(gsm)) || !file.exists(input_path(gsm))) stop("Effect-blind RCTD is incomplete; missing donor result for ", gsm)
  spatial <- readRDS(input_path(gsm))
  aligned <- get_aligned_weights(gsm, spatial$metadata)
  w <- aligned$weights
  weight_sum_qc_rows[[i]] <- aligned$sum_qc
  n_bins <- nrow(w)
  usable <- as.character(w$spot_class) != "reject"
  marker_file <- file.path(weight_dir, paste0(gsm, "_RCTD_holdout_marker_qc.csv"))
  if (!file.exists(marker_file)) stop("Missing disjoint-marker QC for ", gsm)
  marker <- read.csv(marker_file, stringsAsFactors = FALSE)
  if (nrow(marker) != 2L || anyDuplicated(marker$lineage)) stop("Holdout-marker QC table malformed for ", gsm)
  marker_rows[[i]] <- marker
  donor_blind_qc[[i]] <- data.frame(
    GSM = gsm, donor_id = donors_blind$donor_id[i], n_bins = n_bins,
    n_rejected = sum(!usable), rejection_fraction = mean(!usable),
    n_MG_weight_ge_0_5 = sum(usable & w$MG >= 0.5),
    n_ASC_weight_ge_0_5 = sum(usable & w$ASC >= 0.5),
    n_doublet_converged = sum(usable & as.logical(w$conv_doublet), na.rm = TRUE),
    n_usable = sum(usable), stringsAsFactors = FALSE
  )
  if (n_bins < 5000L || sum(usable) == 0L || anyNA(w$conv_doublet[usable])) stop("Primary bin/fit eligibility failed for ", gsm)
  mask <- readRDS(file.path(mask_dir, paste0(gsm, "_bin50_CellBin_mask.rds")))
  if (!all(spatial$metadata$tile_id %in% mask$polygon_mask_tile_id)) stop("Count/RCTD bin outside frozen polygon mask for ", gsm)
  for (lineage in c("MG", "ASC")) {
    col <- lineage
    relevant <- which(usable & w[[col]] >= 0.5)
    genes <- if (lineage == "MG") names(g1_sign) else unique(unlist(module_genes, use.names = FALSE))
    available <- intersect(genes, rownames(spatial$counts))
    det_n <- if (length(available) && length(relevant)) Matrix::rowSums(spatial$counts[available, relevant, drop = FALSE] > 0) else numeric()
    detected <- setNames(rep(FALSE, length(genes)), genes)
    detection_fraction <- setNames(if (length(relevant)) rep(0, length(genes)) else rep(NA_real_, length(genes)), genes)
    if (length(available) && length(relevant)) {
      detection_fraction[available] <- as.numeric(det_n) / length(relevant)
      detected[available] <- det_n / length(relevant) >= 0.01
    }
    coverage_rows[[length(coverage_rows) + 1L]] <- data.frame(
      GSM = gsm, donor_id = donors_blind$donor_id[i], lineage = lineage,
      gene = genes, feature_present = genes %in% rownames(spatial$counts),
      n_weight_supported_bins = length(relevant), detection_fraction = as.numeric(detection_fraction[genes]),
      passes_1pct = as.logical(detected[genes]), stringsAsFactors = FALSE
    )
  }
  rm(spatial, aligned, w, mask); gc(verbose = FALSE)
  message("G6 blind preflight donor complete: ", gsm)
}
marker_qc <- do.call(rbind, marker_rows)
donor_qc <- do.call(rbind, donor_blind_qc)
coverage <- do.call(rbind, coverage_rows)
weight_sum_qc <- do.call(rbind, weight_sum_qc_rows)
write.csv(marker_qc, file.path(outdir, "tables/G6_RCTD_holdout_marker_qc_all_donors.csv"), row.names = FALSE)
write.csv(donor_qc, file.path(outdir, "tables/G6_RCTD_donor_fit_qc.csv"), row.names = FALSE)
write.csv(coverage, file.path(outdir, "tables/G6_state_gene_coverage_by_donor_effect_blind.csv"), row.names = FALSE)
write.csv(weight_sum_qc, file.path(outdir, "tables/G6_RCTD_weight_sum_qc_all_donors.csv"), row.names = FALSE)

marker_summary <- do.call(rbind, lapply(c("MG", "ASC"), function(lineage) {
  x <- marker_qc[marker_qc$lineage == lineage, , drop = FALSE]
  data.frame(lineage = lineage, donors = nrow(x), donors_positive_rho = sum(x$spearman_rho > 0, na.rm = TRUE),
             minimum_disjoint_markers = min(x$n_holdout_markers, na.rm = TRUE),
             median_rho = median(x$spearman_rho, na.rm = TRUE),
             pass = sum(x$spearman_rho > 0, na.rm = TRUE) >= 10L && min(x$n_holdout_markers, na.rm = TRUE) >= 3L)
}))
g1_coverage <- coverage[coverage$lineage == "MG" & coverage$gene %in% names(g1_sign), , drop = FALSE]
g1_supported <- tapply(g1_coverage$passes_1pct, g1_coverage$gene, sum) >= 10L
g1_coverage_pass <- sum(g1_supported, na.rm = TRUE) >= ceiling(0.25 * length(g1_sign))
module_coverage_rows <- lapply(modules, function(m) {
  genes <- module_genes[[m]]
  x <- coverage[coverage$lineage == "ASC" & coverage$gene %in% genes, , drop = FALSE]
  supported <- tapply(x$passes_1pct, x$gene, sum) >= 10L
  needed <- ceiling(0.25 * length(genes))
  data.frame(module = m, n_frozen_one_to_one_genes = length(genes), n_genes_supported_in_ge10_donors = sum(supported, na.rm = TRUE),
             minimum_genes_required = needed, G5_primary_status = module_info$G5_primary_status[match(m, module_info$module)],
             G5_module_gate = module_info$G5_module_gate[match(m, module_info$module)],
             pass_coverage = sum(supported, na.rm = TRUE) >= needed,
             G6_state_estimable = sum(supported, na.rm = TRUE) >= needed &&
               module_info$G5_primary_status[match(m, module_info$module)] == "estimable",
             stringsAsFactors = FALSE)
})
module_coverage <- do.call(rbind, module_coverage_rows)
primary_weight_support_pass <- all(donor_qc$n_MG_weight_ge_0_5 > 0L & donor_qc$n_ASC_weight_ge_0_5 > 0L)
failed_weight_support <- donor_qc[donor_qc$n_MG_weight_ge_0_5 == 0L | donor_qc$n_ASC_weight_ge_0_5 == 0L,
                                  c("GSM", "donor_id", "n_MG_weight_ge_0_5", "n_ASC_weight_ge_0_5"), drop = FALSE]
write.csv(marker_summary, file.path(outdir, "tables/G6_RCTD_holdout_marker_gate.csv"), row.names = FALSE)
write.csv(data.frame(G1_signature_genes = length(g1_sign), G1_supported_genes = sum(g1_supported, na.rm = TRUE),
                     G1_minimum_supported_genes = ceiling(0.25 * length(g1_sign)), G1_coverage_pass = g1_coverage_pass),
          file.path(outdir, "tables/G6_G1_coverage_gate.csv"), row.names = FALSE)
write.csv(module_coverage, file.path(outdir, "tables/G6_G5_module_coverage_gate.csv"), row.names = FALSE)
preflight_pass <- all(marker_summary$pass) && g1_coverage_pass && any(module_coverage$G6_state_estimable) && primary_weight_support_pass
preflight_blockers <- character()
if (!all(marker_summary$pass)) preflight_blockers <- c(preflight_blockers, "Disjoint holdout-marker acceptance criteria failed.")
if (!g1_coverage_pass) preflight_blockers <- c(preflight_blockers, paste0("Frozen G1 coverage failed: ", sum(g1_supported, na.rm = TRUE), " of ", length(g1_sign), " genes met the locked >=1% detection in >=10 of 12 donors rule (", ceiling(0.25 * length(g1_sign)), " required)."))
if (!any(module_coverage$G6_state_estimable)) preflight_blockers <- c(preflight_blockers, "No passing G5 module was estimable in G6 under the locked gene-coverage rule.")
if (!primary_weight_support_pass) preflight_blockers <- c(preflight_blockers, paste0(nrow(failed_weight_support), " of 12 donors lack at least one non-rejected bin with MG and ASC weights >=0.5."))
preflight <- list(run_id = run_id, execution_lock_id = lock_id, effects_calculated = FALSE,
                  protocol_amendment_id = "G6-PROT-2026-09-25-001",
                  execution_addendum_ids = c("G6-EXEC-2026-09-25-002", "G6-EXEC-2026-09-25-003", "G6-EXEC-2026-09-25-004", "G6-EXEC-2026-09-25-005", "G6-EXEC-2026-09-25-006", "G6-EXEC-2026-09-25-007", "G6-EXEC-2026-09-25-008", "G6-EXEC-2026-09-25-009", "G6-EXEC-2026-09-25-010", "G6-EXEC-2026-09-25-011"),
                  all_12_RCTD_outputs_present = TRUE, marker_gate = marker_summary,
                  G1_coverage_pass = g1_coverage_pass, G1_supported_genes = sum(g1_supported, na.rm = TRUE),
                  G5_coverage = module_coverage, reference_overlap_audit = reference_overlap_summary,
                  primary_weight_support_pass = primary_weight_support_pass,
                  donors_failing_primary_weight_support = failed_weight_support,
                  preflight_blockers = as.list(preflight_blockers), preflight_pass = preflight_pass)
jsonlite::write_json(preflight, file.path(outdir, "G6_effect_blind_preflight.json"), auto_unbox = TRUE, pretty = TRUE, na = "null")
if (!preflight_pass) {
  terminal_decision <- list(
    run_id = run_id, execution_lock_id = lock_id,
    protocol = "disease_side_v3.5_latest_frozen",
    protocol_amendment_id = "G6-PROT-2026-09-25-001",
    execution_addendum_ids = preflight$execution_addendum_ids,
    primary_input = "GSE269906 / BA10 / 6 AD and 6 author-labelled Normal Aging controls",
    preflight_pass = FALSE, disease_effects_calculated = FALSE,
    decision = "G6_NOT_EXECUTABLE",
    decision_note = paste(c(preflight_blockers,
                            "GSE233208 has 0 eligible independent AD donors after fail-closed overlap exclusion.",
                            "GSE220442 GSM-to-patient mapping remains unresolved and cannot serve as donor-level validation."), collapse = " "),
    blockers = as.list(c(preflight_blockers,
                         "GSE233208 has 0 eligible independent AD donors after fail-closed overlap exclusion.",
                         "GSE220442 GSM-to-patient mapping remains unresolved and cannot serve as donor-level validation.")),
    independent_validation = list(accession = "GSE233208", eligible_AD_donors = 0L, eligible_control_donors = 5L, validation_pass = FALSE),
    cross_region_support = list(accession = "GSE220442", donor_mapping_resolved = FALSE, validation_pass = FALSE)
  )
  jsonlite::write_json(terminal_decision, file.path(outdir, "G6_decision.json"), auto_unbox = TRUE, pretty = TRUE, na = "null")
  stop("Effect-blind G6 gate failed. No G6 group effects were calculated.")
}

# Diagnosis labels enter only after the locked effect-blind gates pass.
donors <- read.delim(donor_manifest_path, check.names = FALSE, stringsAsFactors = FALSE)
donors$age_midpoint <- as.numeric(sub("^([0-9]+)-([0-9]+)$", "\\1", donors$age_band)) / 2 +
  as.numeric(sub("^([0-9]+)-([0-9]+)$", "\\2", donors$age_band)) / 2
if (any(!is.finite(donors$age_midpoint)) || !all(table(donors$source_diagnosis) == 6L)) stop("Diagnosis/age donor model is not the frozen 6-versus-6 design")
donors$diagnosis <- factor(donors$source_diagnosis, levels = c("Control", "AD"))

# Diagnosis-blind gene-wise logCPM moments pooled over all eligible bins as
# specified in the lock. A missing feature is NA (not an observed zero).
gene_mean_sum <- setNames(numeric(length(all_state_genes)), all_state_genes)
gene_second_sum <- gene_mean_sum
gene_observation_n <- gene_mean_sum
for (i in seq_len(nrow(donors))) {
  gsm <- donors$GSM[i]
  spatial <- readRDS(input_path(gsm))
  genes <- intersect(all_state_genes, rownames(spatial$counts))
  mat <- as.matrix(spatial$counts[genes, , drop = FALSE])
  scale <- 1e6 / pmax(spatial$metadata$total_UMI, 1)
  expr <- log2(sweep(mat, 2L, scale, `*`) + 1)
  gene_mean_sum[genes] <- gene_mean_sum[genes] + rowSums(expr)
  gene_second_sum[genes] <- gene_second_sum[genes] + rowSums(expr * expr)
  gene_observation_n[genes] <- gene_observation_n[genes] + ncol(expr)
  rm(spatial, mat, expr); gc(verbose = FALSE)
  message("G6 diagnosis-blind gene moments complete: ", gsm)
}
gene_mean <- gene_mean_sum / pmax(gene_observation_n, 1)
gene_variance <- gene_second_sum / pmax(gene_observation_n, 1) - gene_mean^2
gene_sd <- sqrt(pmax(gene_variance, 0))
gene_sd[!is.finite(gene_sd) | gene_sd == 0] <- 1

score_rows <- list()
abundance_rows <- list()
spatial_score_objects <- vector("list", nrow(donors))
for (i in seq_len(nrow(donors))) {
  gsm <- donors$GSM[i]
  spatial <- readRDS(input_path(gsm))
  aligned <- get_aligned_weights(gsm, spatial$metadata)
  w <- aligned$weights
  mat <- as.matrix(spatial$counts[intersect(all_state_genes, rownames(spatial$counts)), , drop = FALSE])
  genes <- rownames(mat)
  scale <- 1e6 / pmax(spatial$metadata$total_UMI, 1)
  expr <- log2(sweep(mat, 2L, scale, `*`) + 1)
  z <- sweep(expr, 1L, gene_mean[genes], `-`)
  z <- sweep(z, 1L, gene_sd[genes], `/`)
  z[!is.finite(z)] <- 0
  g1_here <- intersect(names(g1_sign), genes)
  sender_score <- if (length(g1_here)) colMeans(z[g1_here, , drop = FALSE] * g1_sign[g1_here]) else rep(NA_real_, ncol(z))
  module_scores <- matrix(NA_real_, nrow = ncol(z), ncol = length(modules),
                          dimnames = list(NULL, modules))
  for (m in modules) {
    here <- intersect(module_genes[[m]], genes)
    if (length(here) && is.finite(module_info$expected_sign[match(m, module_info$module)])) {
      module_scores[, m] <- colMeans(z[here, , drop = FALSE]) * module_info$expected_sign[match(m, module_info$module)]
    }
  }
  usable <- as.character(w$spot_class) != "reject"
  state_means <- c(G1_sender = weighted.mean(sender_score[usable], w$MG[usable]),
                   setNames(vapply(modules, function(m) weighted.mean(module_scores[usable, m], w$ASC[usable], na.rm = TRUE), numeric(1)), modules))
  score_rows[[i]] <- data.frame(GSM = gsm, donor_id = donors$donor_id[i], state = names(state_means),
                                weighted_mean_score = as.numeric(state_means), stringsAsFactors = FALSE)
  abundance_rows[[i]] <- data.frame(
    GSM = gsm, donor_id = donors$donor_id[i],
    MG_mean_normalized_RCTD_weight = mean(w$MG[usable]),
    ASC_mean_normalized_RCTD_weight = mean(w$ASC[usable]),
    n_nonrejected_bins = sum(usable), stringsAsFactors = FALSE
  )
  spatial_score_objects[[i]] <- list(
    GSM = gsm, metadata = spatial$metadata, weights = w,
    sender_score = sender_score, module_scores = module_scores,
    MG_fraction = mean(w$MG[usable]), ASC_fraction = mean(w$ASC[usable])
  )
  rm(spatial, aligned, mat, expr, z); gc(verbose = FALSE)
  message("G6 state scoring complete: ", gsm)
}
state_donor <- do.call(rbind, score_rows)
state_donor <- merge(state_donor, donors[, c("GSM", "donor_id", "source_diagnosis", "age_midpoint", "NIA_AA_2012_ADNC_level")],
                     by = c("GSM", "donor_id"), all.x = TRUE, sort = FALSE)
abundance_donor <- merge(do.call(rbind, abundance_rows),
                         donors[, c("GSM", "donor_id", "source_diagnosis", "age_midpoint", "NIA_AA_2012_ADNC_level")],
                         by = c("GSM", "donor_id"), all.x = TRUE, sort = FALSE)

fit_contrast <- function(d) {
  d <- d[is.finite(d$weighted_mean_score) & is.finite(d$age_midpoint) & !is.na(d$source_diagnosis), , drop = FALSE]
  ad <- d$weighted_mean_score[d$source_diagnosis == "AD"]
  ctrl <- d$weighted_mean_score[d$source_diagnosis == "Control"]
  if (length(ad) < 2L || length(ctrl) < 2L) {
    return(data.frame(n_AD = length(ad), n_Control = length(ctrl), AD_minus_Control = NA_real_,
                      Welch_CI_low = NA_real_, Welch_CI_high = NA_real_, Welch_p = NA_real_,
                      age_adjusted_effect = NA_real_, age_adjusted_CI_low = NA_real_, age_adjusted_CI_high = NA_real_,
                      age_adjusted_p = NA_real_))
  }
  tt <- tryCatch(t.test(ad, ctrl), error = function(e) NULL)
  fit <- lm(weighted_mean_score ~ diagnosis + age_midpoint,
            data = transform(d, diagnosis = factor(source_diagnosis, levels = c("Control", "AD"))))
  cf <- summary(fit)$coefficients
  term <- "diagnosisAD"
  ci <- tryCatch(confint(fit, term, level = 0.95), error = function(e) c(NA_real_, NA_real_))
  data.frame(n_AD = length(ad), n_Control = length(ctrl), AD_minus_Control = mean(ad) - mean(ctrl),
             Welch_CI_low = if (is.null(tt)) NA_real_ else unname(tt$conf.int[1]),
             Welch_CI_high = if (is.null(tt)) NA_real_ else unname(tt$conf.int[2]),
             Welch_p = if (is.null(tt)) NA_real_ else tt$p.value,
             age_adjusted_effect = cf[term, "Estimate"], age_adjusted_CI_low = ci[1], age_adjusted_CI_high = ci[2],
             age_adjusted_p = cf[term, "Pr(>|t|)"], stringsAsFactors = FALSE)
}
abundance_metrics <- c("MG_mean_normalized_RCTD_weight", "ASC_mean_normalized_RCTD_weight")
abundance_effects <- do.call(rbind, lapply(abundance_metrics, function(metric) {
  d <- abundance_donor
  d$weighted_mean_score <- d[[metric]]
  cbind(cell_type = if (startsWith(metric, "MG_")) "MG" else "ASC", fit_contrast(d))
}))
abundance_effects$FDR_BH_2 <- p.adjust(abundance_effects$age_adjusted_p, method = "BH")
write.csv(abundance_donor, file.path(outdir, "tables/G6_cell_type_abundance_by_donor.csv"), row.names = FALSE)
write.csv(abundance_effects, file.path(outdir, "tables/G6_cell_type_abundance_group_effects.csv"), row.names = FALSE)
state_effects <- do.call(rbind, lapply(split(state_donor, state_donor$state), function(d) cbind(state = d$state[1], fit_contrast(d))))
state_effects$FDR_BH_8 <- p.adjust(state_effects$age_adjusted_p, method = "BH")
state_effects$direction_concordant <- state_effects$age_adjusted_effect > 0
state_effects$state_estimable <- state_effects$n_AD == 6L & state_effects$n_Control == 6L & is.finite(state_effects$age_adjusted_p)
write.csv(state_donor, file.path(outdir, "tables/G6_state_scores_by_donor.csv"), row.names = FALSE)
write.csv(state_effects, file.path(outdir, "tables/G6_state_group_effects.csv"), row.names = FALSE)

# Grid geometry and mask-constrained edges are constructed independently for
# each section; only 25 um bin centers and CellBin polygons define adjacency.
offsets_for_radius <- function(radius_bins) {
  out <- expand.grid(dx = -radius_bins:radius_bins, dy = -radius_bins:radius_bins)
  out <- out[out$dx^2 + out$dy^2 > 0 & out$dx^2 + out$dy^2 <= radius_bins^2, , drop = FALSE]
  out[order(out$dx, out$dy), , drop = FALSE]
}
bresenham <- function(dx, dy) {
  x <- 0L; y <- 0L; x1 <- as.integer(dx); y1 <- as.integer(dy)
  ax <- abs(x1); ay <- abs(y1); sx <- if (x1 >= 0L) 1L else -1L; sy <- if (y1 >= 0L) 1L else -1L
  err <- ax - ay; cells <- matrix(c(x, y), ncol = 2L)
  while (x != x1 || y != y1) {
    e2 <- 2L * err
    if (e2 > -ay) { err <- err - ay; x <- x + sx }
    if (e2 < ax) { err <- err + ax; y <- y + sy }
    cells <- rbind(cells, c(x, y))
  }
  cells
}
make_sender_hit_map <- function(sender_idx, mask, nx, ny, radius_bins) {
  hit <- logical(length(mask))
  sx <- (sender_idx - 1L) %% nx
  sy <- (sender_idx - 1L) %/% nx
  offsets <- offsets_for_radius(radius_bins)
  for (j in seq_len(nrow(offsets))) {
    dx <- offsets$dx[j]; dy <- offsets$dy[j]
    tx <- sx - dx; ty <- sy - dy
    inside <- tx >= 0L & tx < nx & ty >= 0L & ty < ny
    if (!any(inside)) next
    si <- sender_idx[inside]; ti <- as.integer(tx[inside] + nx * ty[inside]) + 1L
    line <- bresenham(dx, dy)
    ok <- mask[ti] & mask[si]
    if (nrow(line) > 2L) for (q in 2L:(nrow(line) - 1L)) {
      px <- tx[inside] + line[q, 1L]
      py <- ty[inside] + line[q, 2L]
      ok <- ok & mask[as.integer(px + nx * py) + 1L]
    }
    hit[ti[ok]] <- TRUE
  }
  hit
}
candidate_shifts <- function(receiver_idx, mask, nx, ny, seed, nperm = 999L) {
  grid <- expand.grid(dx = -40L:40L, dy = -40L:40L)
  distance <- sqrt(grid$dx^2 + grid$dy^2)
  grid <- grid[distance >= 10 & distance <= 40, , drop = FALSE]
  rx <- (receiver_idx - 1L) %% nx; ry <- (receiver_idx - 1L) %/% nx
  legal <- logical(nrow(grid))
  for (j in seq_len(nrow(grid))) {
    tx <- rx + grid$dx[j]; ty <- ry + grid$dy[j]
    in_bounds <- tx >= 0L & tx < nx & ty >= 0L & ty < ny
    ti <- as.integer(tx[in_bounds] + nx * ty[in_bounds]) + 1L
    overlap <- sum(mask[ti]) / length(receiver_idx)
    legal[j] <- overlap >= 0.8
  }
  grid <- grid[legal, , drop = FALSE]
  if (nrow(grid) < nperm) stop("Fewer than 999 unique mask-constrained 250-1000 um translations")
  set.seed(seed)
  grid[sample(seq_len(nrow(grid)), nperm, replace = FALSE), , drop = FALSE]
}
translate_idx <- function(idx, dx, dy, nx, ny) {
  x <- (idx - 1L) %% nx + dx
  y <- (idx - 1L) %/% nx + dy
  inside <- x >= 0L & x < nx & y >= 0L & y < ny
  out <- rep(NA_integer_, length(idx))
  out[inside] <- as.integer(x[inside] + nx * y[inside]) + 1L
  out
}
translation_neighbor_fraction <- function(hit, translated_idx) {
  if (!length(translated_idx) || any(!is.na(translated_idx) &
      (translated_idx < 1L | translated_idx > length(hit)))) {
    stop("Translated receiver map has invalid support")
  }
  in_bounds <- !is.na(translated_idx)
  sum(hit[translated_idx[in_bounds]]) / length(translated_idx)
}

# Deterministic boundary microfixture: the off-grid receiver slot remains in
# the four-bin denominator; only one of the remaining slots has a sender hit.
translation_fixture <- translation_neighbor_fraction(
  hit = c(FALSE, TRUE, FALSE, FALSE, FALSE, FALSE),
  translated_idx = c(2L, 4L, NA_integer_, 6L)
)
if (!identical(translation_fixture, 0.25)) {
  stop("Fixed-denominator spatial translation boundary check failed")
}

spatial_rows <- list()
for (i in seq_len(nrow(donors))) {
  obj <- spatial_score_objects[[i]]
  gsm <- obj$GSM
  mask_obj <- readRDS(file.path(mask_dir, paste0(gsm, "_bin50_CellBin_mask.rds")))
  nx <- mask_obj$grid$nx; ny <- mask_obj$grid$ny
  mask <- logical(nx * ny); mask[mask_obj$polygon_mask_tile_id] <- TRUE
  w <- obj$weights
  meta <- obj$metadata
  tile <- as.integer(meta$tile_id)
  if (anyDuplicated(tile) || any(!mask[tile])) stop("Spatial bin mapping or mask failed for ", gsm)
  for (weight_threshold in c(0.4, 0.5, 0.7)) {
    mg_ok <- as.character(w$spot_class) != "reject" & w$MG >= weight_threshold
    mg_cut <- quantile(obj$sender_score[mg_ok], 0.75, names = FALSE, type = 7)
    sender_idx <- tile[mg_ok & obj$sender_score >= mg_cut]
    if (!length(sender_idx)) stop("No high-state MG senders in ", gsm, " at weight threshold ", weight_threshold)
    receiver_map_idx <- tile[as.character(w$spot_class) != "reject" & w$ASC >= weight_threshold]
    shifts <- candidate_shifts(receiver_map_idx, mask, nx, ny,
                               seed = 42000L + i * 100L + as.integer(weight_threshold * 10))
    module_maps <- lapply(modules, function(m) {
      score <- obj$module_scores[, m]
      eligible <- as.character(w$spot_class) != "reject" & w$ASC >= weight_threshold & is.finite(score)
      is_estimable <- module_coverage$G6_state_estimable[match(m, module_coverage$module)] && any(eligible)
      if (!is_estimable) return(list(estimable = FALSE, n_high = 0L, observed_idx = NULL, shifted = NULL))
      cut <- quantile(score[eligible], 0.75, names = FALSE, type = 7)
      receiver_high_idx <- tile[eligible & score >= cut]
      if (!length(receiver_high_idx)) stop("No high-state ASC receiver bins in ", gsm, " / ", m)
      shifted <- lapply(seq_len(nrow(shifts)), function(p) {
        translate_idx(receiver_high_idx, shifts$dx[p], shifts$dy[p], nx, ny)
      })
      if (any(lengths(shifted) != length(receiver_high_idx))) {
        stop("A translation did not retain every high receiver-bin slot")
      }
      list(estimable = TRUE, n_high = length(receiver_high_idx),
           observed_idx = receiver_high_idx, shifted = shifted)
    })
    names(module_maps) <- modules
    for (radius_um in c(25, 50, 75)) {
      radius_bins <- as.integer(radius_um / 25)
      hit <- make_sender_hit_map(sender_idx, mask, nx, ny, radius_bins)
      for (m in modules) {
        receiver <- module_maps[[m]]
        if (!receiver$estimable) {
          spatial_rows[[length(spatial_rows) + 1L]] <- data.frame(
            GSM = gsm, donor_id = donors$donor_id[i], module = m, radius_um = radius_um,
            cell_weight_threshold = weight_threshold, n_high_MG_sender = length(sender_idx),
            n_high_ASC_receiver = 0L, n_ASC_map_bins_for_translation = length(receiver_map_idx),
            n_legal_unique_translations = nrow(shifts), observed_fraction_with_high_MG_neighbor = NA_real_,
            null_mean_fraction = NA_real_, null_sd = NA_real_, spatial_excess = NA_real_,
            translation_tail_p = NA_real_, stringsAsFactors = FALSE
          )
          next
        }
        observed <- mean(hit[receiver$observed_idx])
        null <- vapply(receiver$shifted, function(shifted) {
          translation_neighbor_fraction(hit, shifted)
        }, numeric(1))
        spatial_rows[[length(spatial_rows) + 1L]] <- data.frame(
          GSM = gsm, donor_id = donors$donor_id[i], module = m, radius_um = radius_um,
          cell_weight_threshold = weight_threshold, n_high_MG_sender = length(sender_idx),
          n_high_ASC_receiver = receiver$n_high, n_ASC_map_bins_for_translation = length(receiver_map_idx),
          n_legal_unique_translations = nrow(shifts), observed_fraction_with_high_MG_neighbor = observed,
          null_mean_fraction = mean(null), null_sd = sd(null), spatial_excess = observed - mean(null),
          translation_tail_p = (1 + sum(null >= observed)) / (length(null) + 1),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  message("G6 spatial map complete: ", gsm)
}
spatial_donor <- do.call(rbind, spatial_rows)
spatial_donor <- merge(spatial_donor, donors[, c("GSM", "source_diagnosis", "age_midpoint", "NIA_AA_2012_ADNC_level")],
                       by = "GSM", all.x = TRUE, sort = FALSE)
spatial_effects <- do.call(rbind, lapply(split(spatial_donor, interaction(spatial_donor$cell_weight_threshold, spatial_donor$radius_um, spatial_donor$module, drop = TRUE)), function(d) {
  contrast <- fit_contrast(transform(d, weighted_mean_score = spatial_excess))
  cbind(module = d$module[1], radius_um = d$radius_um[1], cell_weight_threshold = d$cell_weight_threshold[1], contrast)
}))
spatial_effects$FDR_BH_7_within_radius <- NA_real_
for (group in unique(interaction(spatial_effects$cell_weight_threshold, spatial_effects$radius_um, drop = TRUE))) {
  idx <- which(interaction(spatial_effects$cell_weight_threshold, spatial_effects$radius_um, drop = TRUE) == group)
  spatial_effects$FDR_BH_7_within_radius[idx] <- p.adjust(spatial_effects$age_adjusted_p[idx], method = "BH")
}
write.csv(spatial_donor, file.path(outdir, "tables/G6_neighborhood_excess_by_donor.csv"), row.names = FALSE)
write.csv(spatial_effects, file.path(outdir, "tables/G6_neighborhood_group_effects.csv"), row.names = FALSE)

primary_state <- state_effects[state_effects$state == "G1_sender", , drop = FALSE]
primary_spatial <- spatial_effects[spatial_effects$radius_um == 25 & spatial_effects$cell_weight_threshold == 0.5, , drop = FALSE]
pass_modules <- module_info$module[module_info$G5_module_gate == "PASS" & module_info$G5_primary_status == "estimable"]
primary_spatial$state_FDR <- state_effects$FDR_BH_8[match(primary_spatial$module, state_effects$state)]
primary_spatial$state_direction_concordant <- state_effects$direction_concordant[match(primary_spatial$module, state_effects$state)]
primary_spatial$state_estimable <- module_coverage$G6_state_estimable[match(primary_spatial$module, module_coverage$module)]
primary_spatial$G5_v3_4_pass <- primary_spatial$module %in% pass_modules
primary_spatial$spatial_gate_pass <- primary_spatial$n_AD == 6L & primary_spatial$n_Control == 6L &
  primary_spatial$G5_v3_4_pass & primary_spatial$state_estimable &
  primary_spatial$state_direction_concordant & primary_spatial$state_FDR < 0.05 &
  primary_spatial$age_adjusted_effect > 0 & primary_spatial$FDR_BH_7_within_radius < 0.05
validation <- list(GSE233208_independent_AD_n = 0L, GSE233208_independent_control_n = 5L,
                   GSE233208_eligible = FALSE, GSE220442_GSM_patient_mapping_resolved = FALSE,
                   validation_pass = FALSE)
g6_pass <- isTRUE(primary_state$state_estimable[1]) &&
  isTRUE(primary_state$direction_concordant[1]) && primary_state$FDR_BH_8[1] < 0.05 &&
  any(primary_spatial$spatial_gate_pass) && validation$validation_pass

# Frozen pathology-only control sensitivity and donor influence checks.
not_adnc_state <- state_donor[state_donor$source_diagnosis == "AD" |
                                state_donor$NIA_AA_2012_ADNC_level == "Not_AD_neuropathologic_change", , drop = FALSE]
not_adnc_state_effects <- do.call(rbind, lapply(split(not_adnc_state, not_adnc_state$state), function(d) cbind(state = d$state[1], fit_contrast(d))))
not_adnc_spatial <- spatial_donor[spatial_donor$radius_um == 25 & spatial_donor$cell_weight_threshold == 0.5 &
                                    (spatial_donor$source_diagnosis == "AD" |
                                       spatial_donor$NIA_AA_2012_ADNC_level == "Not_AD_neuropathologic_change"), , drop = FALSE]
not_adnc_spatial_effects <- do.call(rbind, lapply(split(not_adnc_spatial, not_adnc_spatial$module), function(d) {
  cbind(module = d$module[1], fit_contrast(transform(d, weighted_mean_score = spatial_excess)))
}))
loo_rows <- list()
for (state in unique(state_donor$state)) {
  d <- state_donor[state_donor$state == state, , drop = FALSE]
  for (drop_gsm in d$GSM) {
    res <- fit_contrast(d[d$GSM != drop_gsm, , drop = FALSE])
    loo_rows[[length(loo_rows) + 1L]] <- cbind(analysis = "state", metric = state, dropped_GSM = drop_gsm,
                                               effect = res$age_adjusted_effect, p = res$age_adjusted_p,
                                               n_AD = res$n_AD, n_Control = res$n_Control)
  }
}
for (m in modules) {
  d <- spatial_donor[spatial_donor$module == m & spatial_donor$radius_um == 25 & spatial_donor$cell_weight_threshold == 0.5, , drop = FALSE]
  for (drop_gsm in d$GSM) {
    res <- fit_contrast(transform(d[d$GSM != drop_gsm, , drop = FALSE], weighted_mean_score = spatial_excess))
    loo_rows[[length(loo_rows) + 1L]] <- cbind(analysis = "neighborhood", metric = m, dropped_GSM = drop_gsm,
                                               effect = res$age_adjusted_effect, p = res$age_adjusted_p,
                                               n_AD = res$n_AD, n_Control = res$n_Control)
  }
}
loo_effects <- do.call(rbind, loo_rows)
write.csv(not_adnc_state_effects, file.path(outdir, "tables/G6_state_effects_Not_ADNC_only_controls.csv"), row.names = FALSE)
write.csv(not_adnc_spatial_effects, file.path(outdir, "tables/G6_neighborhood_effects_Not_ADNC_only_controls.csv"), row.names = FALSE)
write.csv(loo_effects, file.path(outdir, "tables/G6_leave_one_donor_out.csv"), row.names = FALSE)
write.csv(data.frame(sensitivity = "leave_one_reference_donor_out_RCTD_refit", performed = FALSE,
                     reason = "Requires 18 x 12 additional section-wide RCTD refits after the primary reference fit; one-core runtime exceeded feasible local execution window.",
                     stringsAsFactors = FALSE),
          file.path(outdir, "tables/G6_reference_donor_sensitivity_status.csv"), row.names = FALSE)
decision <- list(
  run_id = run_id, execution_lock_id = lock_id, protocol = "disease_side_v3.5_latest_frozen",
  protocol_amendment_id = "G6-PROT-2026-09-25-001",
  execution_addendum_ids = c("G6-EXEC-2026-09-25-002", "G6-EXEC-2026-09-25-003", "G6-EXEC-2026-09-25-004", "G6-EXEC-2026-09-25-005", "G6-EXEC-2026-09-25-006", "G6-EXEC-2026-09-25-007", "G6-EXEC-2026-09-25-008", "G6-EXEC-2026-09-25-009", "G6-EXEC-2026-09-25-010", "G6-EXEC-2026-09-25-011"),
  primary_input = "GSE269906 / BA10 / 6 AD and 6 author-labelled Normal Aging controls",
  preflight_pass = preflight_pass, disease_effects_calculated = TRUE,
  state_primary_radius_um = 25, state_donor_table = "tables/G6_state_scores_by_donor.csv",
  state_effect_table = "tables/G6_state_group_effects.csv",
  neighborhood_donor_table = "tables/G6_neighborhood_excess_by_donor.csv",
  neighborhood_effect_table = "tables/G6_neighborhood_group_effects.csv",
  primary_validations = validation,
  primary_reference_identity_audit = reference_overlap_summary,
  decision = if (g6_pass) "G6_PASS" else "G6_NOT_PASS",
  decision_note = if (!validation$validation_pass) "Independent AD validation is ineligible under the frozen fail-closed overlap rule; primary signals are reported as supportive only." else "Primary gate criteria evaluated under the frozen protocol.",
  G1_primary_effect = as.list(primary_state[1, , drop = FALSE]),
  primary_receiver_neighborhood = primary_spatial,
  control_pathology_sensitivity_required = TRUE,
  layer_ROI_histology_covariates_available = FALSE,
  interpretation = "Bin-level deconvolution weights and expression-state scores are inferred, not segmented cell identities; spatial association is not causal."
)
jsonlite::write_json(decision, file.path(outdir, "G6_decision.json"), auto_unbox = TRUE, pretty = TRUE, na = "null")

# Compact donor-level plots with individual donors displayed; no spot-level
# pseudo-replication is used in the inferential test.
pdf(file.path(outdir, "figures/G6_primary_donor_level_effects.pdf"), width = 14, height = 9)
par(mfrow = c(2, 4), mar = c(7, 4, 3, 1))
for (state in c("G1_sender", modules)) {
  d <- state_donor[state_donor$state == state, , drop = FALSE]
  d <- d[is.finite(d$weighted_mean_score), , drop = FALSE]
  if (!nrow(d)) {
    plot.new(); title(main = paste(state, "not estimable"))
    next
  }
  d$source_diagnosis <- factor(d$source_diagnosis, levels = c("Control", "AD"))
  boxplot(weighted_mean_score ~ source_diagnosis, data = d, outline = FALSE,
          xlab = "Diagnosis", ylab = "Donor-level weighted state score", main = state,
          col = c("grey85", "#80B1D3"), names = c("Normal Aging", "AD"))
  stripchart(weighted_mean_score ~ source_diagnosis, data = d, vertical = TRUE,
             method = "jitter", add = TRUE, pch = 19, col = "#333333")
}
dev.off()
pdf(file.path(outdir, "figures/G6_primary_neighborhood_donor_effects.pdf"), width = 14, height = 9)
par(mfrow = c(2, 4), mar = c(7, 4, 3, 1))
for (m in modules) {
  d <- spatial_donor[spatial_donor$module == m & spatial_donor$radius_um == 25 &
                       spatial_donor$cell_weight_threshold == 0.5 & is.finite(spatial_donor$spatial_excess), , drop = FALSE]
  if (!nrow(d)) {
    plot.new(); title(main = paste(m, "not estimable"))
    next
  }
  d$source_diagnosis <- factor(d$source_diagnosis, levels = c("Control", "AD"))
  boxplot(spatial_excess ~ source_diagnosis, data = d, outline = FALSE,
          xlab = "Diagnosis", ylab = "Observed minus translated-null fraction", main = m,
          col = c("grey85", "#80B1D3"), names = c("Normal Aging", "AD"))
  abline(h = 0, lty = 2, col = "grey35")
  stripchart(spatial_excess ~ source_diagnosis, data = d, vertical = TRUE,
             method = "jitter", add = TRUE, pch = 19, col = "#333333")
}
dev.off()
message("G6 primary analysis completed; decision=", decision$decision)
