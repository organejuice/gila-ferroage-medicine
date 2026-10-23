#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(Matrix))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: 86_summarize_g6_multiscale_blind.R <project_root>")
root <- args[[1L]]
base <- file.path(root, "data/interim/G6_GSE269906_multiscale_001")
audit <- file.path(root, "results/2026-09-24-disease-side-v3.5-g6-input-audit/registries")
outdir <- file.path(root, "results/2026-09-26-g6-multiscale-blind-pilot")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
manifest <- read.delim(file.path(audit, "GSE269906_tissue_CNGB_file_manifest.tsv"))
g1 <- read.delim(gzfile(file.path(audit, "G1_microglia_signed_sender_244.tsv.gz")))
g5 <- read.delim(gzfile(file.path(audit, "G5_receiver_one_to_one_human_ortholog_registry.tsv.gz")))
if (nrow(manifest) != 12L || nrow(g1) != 244L) stop("Frozen manifest or G1 registry changed")
rows <- list()
g5_rows_all <- list()
for (width in c(50L, 110L)) for (gsm in manifest$GSM) {
  prefix <- paste0(gsm, "_bin", width)
  sigpath <- file.path(base, paste0("bin", width), "signature", paste0(prefix, "_tissue_signature_counts.rds"))
  fitpath <- file.path(base, paste0("bin", width), "rctd_full_mode", paste0(prefix, "_rctd_full.rds"))
  if (!file.exists(sigpath) || !file.exists(fitpath)) stop("Missing completed blind input/fit: ", prefix)
  sig <- readRDS(sigpath); fit <- readRDS(fitpath)
  matched <- match(fit$metadata$tile_id, sig$metadata$tile_id)
  if (anyNA(matched) || !identical(sig$metadata$donor_id[matched], fit$metadata$donor_id))
    stop("Fit/signature bin mapping mismatch: ", prefix)
  mg <- fit$weights[as.character(fit$metadata$tile_id), "MG"]
  asc <- fit$weights[as.character(fit$metadata$tile_id), "ASC"]
  support <- is.finite(mg) & mg >= 0.5
  asc_support <- is.finite(asc) & asc >= 0.5
  gene_rows <- match(toupper(g1$gene), rownames(sig$counts))
  detected <- integer(nrow(g1))
  if (any(support)) {
    valid <- which(!is.na(gene_rows))
    detected[valid] <- as.integer(Matrix::rowMeans(sig$counts[gene_rows[valid], matched[support], drop = FALSE] > 0) >= 0.01)
  }
  module_cov <- list()
  for (module in sort(unique(g5$module))) {
    genes <- unique(toupper(g5$human_one_to_one_symbol[g5$module == module]))
    genes <- genes[!is.na(genes) & nzchar(genes)]
    gene_row <- match(genes, rownames(sig$counts))
    present <- which(!is.na(gene_row))
    prevalence <- rep(NA_real_, length(genes))
    if (length(present) && any(asc_support)) {
      prevalence[present] <- Matrix::rowMeans(sig$counts[gene_row[present], matched[asc_support], drop = FALSE] > 0)
    }
    detected_n <- sum(prevalence >= 0.01, na.rm = TRUE)
    min_detected <- ceiling(0.25 * length(genes))
    module_cov[[module]] <- detected_n
    g5_rows_all[[length(g5_rows_all) + 1L]] <- data.frame(
      GSM = gsm, donor_id = sig$metadata$donor_id[[1L]], width_DNB = width,
      module = module, n_module_genes = length(genes), n_detected_ge_1pct = detected_n,
      minimum_25pct = min_detected, pass_donor_module = detected_n >= min_detected,
      n_ASC_ge_0_5 = sum(asc_support), stringsAsFactors = FALSE)
  }
  rows[[length(rows) + 1L]] <- data.frame(
    GSM = gsm, donor_id = sig$metadata$donor_id[[1L]], width_DNB = width,
    width_um = width * 0.5, n_bins = nrow(fit$metadata),
    median_features = median(fit$metadata$n_features), median_UMI = median(fit$metadata$total_UMI),
    n_MG_ge_0_5_diagnostic = sum(support), G1_genes_ge_1pct_in_MG_bins = sum(detected),
    n_ASC_ge_0_5_diagnostic = sum(asc_support),
    n_disjoint_MG_holdout_markers = fit$qc$n_disjoint_MG_holdout_markers,
    n_disjoint_ASC_holdout_markers = fit$qc$n_disjoint_ASC_holdout_markers,
    MG_marker_rho = fit$qc$MG_marker_rho, ASC_marker_rho = fit$qc$ASC_marker_rho,
    source_sha256 = sig$construction$source_sha256, fit_signature = fit$signature)
}
tab <- do.call(rbind, rows)
module_tab <- do.call(rbind, g5_rows_all)
if (nrow(tab) != 24L || anyDuplicated(tab[, c("GSM", "width_DNB")])) stop("Expected 12 donor sections at two widths")
write.table(tab, file.path(outdir, "effect_blind_scale_comparison.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(module_tab, file.path(outdir, "effect_blind_g5_module_coverage.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
gates <- list()
for (width in c(50L, 110L)) {
  d <- tab[tab$width_DNB == width, , drop = FALSE]
  m <- module_tab[module_tab$width_DNB == width, , drop = FALSE]
  g1_donors <- sum(d$G1_genes_ge_1pct_in_MG_bins >= 61L)
  g5_by_module <- aggregate(pass_donor_module ~ module, m, sum)
  g5_all <- nrow(g5_by_module) == length(unique(g5$module)) && all(g5_by_module$pass_donor_module >= 10L)
  marker_mg <- sum(is.finite(d$MG_marker_rho) & d$MG_marker_rho > 0)
  marker_asc <- sum(is.finite(d$ASC_marker_rho) & d$ASC_marker_rho > 0)
  support_pass <- all(d$n_MG_ge_0_5_diagnostic >= 20L & d$n_ASC_ge_0_5_diagnostic >= 20L)
  g1_pass <- g1_donors >= 10L
  marker_pass <- marker_mg >= 10L && marker_asc >= 10L &&
    all(d$n_disjoint_MG_holdout_markers >= 3L) && all(d$n_disjoint_ASC_holdout_markers >= 3L)
  passed <- support_pass && g1_pass && g5_all && marker_pass
  gates[[as.character(width)]] <- data.frame(
    width_DNB = width, donor_rows = nrow(d), support_gate = support_pass,
    MG_donors_ge_20_bins = sum(d$n_MG_ge_0_5_diagnostic >= 20L),
    ASC_donors_ge_20_bins = sum(d$n_ASC_ge_0_5_diagnostic >= 20L),
    MG_marker_gate = marker_mg >= 10L && all(d$n_disjoint_MG_holdout_markers >= 3L),
    ASC_marker_gate = marker_asc >= 10L && all(d$n_disjoint_ASC_holdout_markers >= 3L),
    G1_coverage_gate = g1_pass, G1_donors_ge_61_genes = g1_donors,
    G5_module_gate = g5_all, G5_modules_passing = sum(g5_by_module$pass_donor_module >= 10L),
    G5_failed_modules = paste(g5_by_module$module[g5_by_module$pass_donor_module < 10L], collapse = ";"),
    candidate_scale_gate = passed, decision = if (passed) "ELIGIBLE_FOR_NEXT_PROTOCOL_REVIEW" else "NOT_ELIGIBLE",
    stringsAsFactors = FALSE)
}
gate_tab <- do.call(rbind, gates)
write.table(gate_tab, file.path(outdir, "effect_blind_scale_gates.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
readme_lines <- c(
  "# G6 GSE269906 multiscale diagnosis-blind result",
  "",
  "Decision follows gates locked in G6-MS-2026-09-26-002 before fitting. Eligibility means entry to a new protocol review only; it does not authorize group effects or imply G6_PASS.",
  "No diagnosis labels or AD-control effects were used.",
  "",
  "| DNB width | Physical width | MG donors with >=20 bins | ASC donors with >=20 bins | G1 donors with >=61 genes | G5 modules passing | Failed gate | Decision |",
  "|---:|---:|---:|---:|---:|---:|---|---|"
)
for (i in seq_len(nrow(gate_tab))) {
  g <- gate_tab[i, , drop = FALSE]
  failed <- c(
    if (!g$support_gate) "MG/ASC support" else NULL,
    if (!g$MG_marker_gate) "MG markers" else NULL,
    if (!g$ASC_marker_gate) "ASC markers" else NULL,
    if (!g$G1_coverage_gate) "G1 coverage" else NULL,
    if (!g$G5_module_gate) paste0("G5: ", g$G5_failed_modules) else NULL
  )
  readme_lines <- c(readme_lines, sprintf(
    "| %d | %d um | %d/12 | %d/12 | %d/12 | %d/7 | %s | %s |",
    g$width_DNB, g$width_DNB * 0.5, g$MG_donors_ge_20_bins,
    g$ASC_donors_ge_20_bins, g$G1_donors_ge_61_genes, g$G5_modules_passing,
    paste(failed, collapse = ", "), g$decision))
}
readme_lines <- c(readme_lines, "",
  "The 55 um scale is a mixed spatial unit and would need a separately frozen neighborhood definition even if eligible.",
  "The source-data `nUMI` value is the whole-library depth while RCTD is fitted to a prespecified gene subset; spacexr emitted the expected mismatch warning. Quantification with subset-only depth was not part of this locked pilot and remains a method sensitivity to resolve before a confirmatory analysis.")
writeLines(readme_lines, file.path(outdir, "readme.md"))
cat("effect_blind_rows=", nrow(tab), "; scales_eligible_for_review=", sum(gate_tab$candidate_scale_gate), "\n", sep = "")
