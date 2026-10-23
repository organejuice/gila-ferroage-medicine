#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, width = 200)
suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
})

root <- "."
run_id <- "2026-09-20-disease-side-v3.1"
out_root <- file.path("results", run_id)
dir.create(file.path(out_root, "effects"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_root, "meta"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_root, "gates"), recursive = TRUE, showWarnings = FALSE)

modules <- list(
  M1_Glycolysis = c("HK1", "HK2", "PFKFB3", "PFKFB4", "ALDOA", "ALDOC", "GAPDH", "PGK1", "PGAM1", "ENO1", "ENO2", "PKM", "LDHA", "PDK1", "SLC2A1", "HIF1A", "PFKP"),
  M3_SASP_Senescence = c("CDKN1A", "CDKN2A", "IL6", "IL1B", "TNF", "CCL2", "CCL3", "CCL4", "CSF1", "CXCL10", "MMP9", "MMP14", "SERPINE1", "IGFBP3", "IGFBP7", "NFKB1", "RELA", "TGFB1"),
  M5_Iron_Uptake = c("TFRC", "SLC11A2", "STEAP3", "DMT1"),
  M5_Iron_Storage = c("FTH1", "FTL"),
  M5_Iron_Export = c("SLC40A1", "CP"),
  M5_Ferritinophagy = c("NCOA4", "ATG3", "ATG5", "ATG7", "MAP1LC3B", "SQSTM1"),
  M6_LPO_ACSL4 = c("ACSL4", "ALOX15", "ALOX5", "ALOX12", "NOX4", "GPX4", "LPCAT3", "SAT1"),
  M7_FerroAging_Author = character()
)
ref_file <- file.path("data", "raw", "D5_CellMetab2026", "Table_S2_ferro_aging_reference.xlsx")
ref_tsv <- file.path("results", "archive", "2026-09-20-pre-v3-freeze", "results", "tables", "D5_ferroaging_reference.tsv")
if (file.exists(ref_tsv)) {
  ref <- tryCatch(read.delim(ref_tsv, stringsAsFactors = FALSE), error = function(e) NULL)
  if (!is.null(ref) && "gene_symbol" %in% names(ref)) modules$M7_FerroAging_Author <- unique(as.character(ref$gene_symbol))
}

read_pb <- function(count_file, meta_file, cohort) {
  x <- read.delim(gzfile(count_file), check.names = FALSE, stringsAsFactors = FALSE)
  genes <- make.unique(as.character(x[[1]]))
  counts <- as.matrix(x[, -1, drop = FALSE])
  storage.mode(counts) <- "numeric"
  rownames(counts) <- genes
  meta <- read.delim(meta_file, check.names = FALSE, stringsAsFactors = FALSE)
  meta <- meta[match(colnames(counts), meta$group_id), , drop = FALSE]
  if (anyNA(meta$group_id)) stop("metadata/count columns do not match for ", cohort)
  if (!"n_cells_used" %in% names(meta) && "n_cells" %in% names(meta)) meta$n_cells_used <- meta$n_cells
  if (!"RIN" %in% names(meta) && "rin" %in% names(meta)) meta$RIN <- meta$rin
  meta$cohort <- cohort
  meta$disease <- factor(meta$disease, levels = c("Control", "AD"))
  if ("sex" %in% names(meta)) meta$sex <- factor(meta$sex)
  if ("batch" %in% names(meta)) meta$batch <- factor(meta$batch)
  list(counts = counts, meta = meta)
}

inputs <- list(
  D3_GSE160936 = list(count = file.path("data", "interim", "D3_pseudobulk_counts_microglia.tsv.gz"), meta = file.path("data", "interim", "D3_pseudobulk_metadata_microglia.tsv")),
  D7_GSE148822 = list(count = file.path("data", "interim", "D7_GSE148822", "D7_pseudobulk_counts_microglia.tsv.gz"), meta = file.path("data", "interim", "D7_GSE148822", "D7_pseudobulk_metadata_microglia.tsv")),
  D8_GSE174367 = list(count = file.path("data", "interim", "D8_GSE174367", "D8_pseudobulk_counts_microglia.tsv.gz"), meta = file.path("data", "interim", "D8_GSE174367", "D8_pseudobulk_metadata_microglia.tsv"))
)

build_design <- function(meta, cohort, include_cell_count = FALSE) {
  keep <- rep(TRUE, nrow(meta))
  vars <- c("age", "sex", "RIN")
  # D8 has an auditable library batch; D3 and D7 do not.
  if (cohort == "D8_GSE174367" && sum(!is.na(meta$batch)) == nrow(meta) && nlevels(droplevels(meta$batch)) > 1) vars <- c(vars, "batch")
  if (include_cell_count) vars <- c(vars, "log_cells")
  if (include_cell_count) meta$log_cells <- log1p(meta$n_cells_used)
  for (v in vars) if (v %in% names(meta)) keep <- keep & !is.na(meta[[v]])
  meta <- meta[keep, , drop = FALSE]
  meta$sex <- droplevels(factor(meta$sex))
  meta$disease <- droplevels(factor(meta$disease, levels = c("Control", "AD")))
  if ("batch" %in% vars) meta$batch <- droplevels(factor(meta$batch))
  if (nlevels(meta$disease) != 2) stop("both disease groups are required")
  rhs <- c("disease", vars)
  form <- as.formula(paste("~", paste(rhs, collapse = " + ")))
  list(meta = meta, design = model.matrix(form, data = meta), variables = rhs)
}

fit_one_region <- function(cohort, region, obj, include_cell_count = FALSE, complete_only = FALSE) {
  m <- obj$meta
  keep <- m$region == region & m$n_cells_used >= 20
  m <- m[keep, , drop = FALSE]
  if (complete_only && cohort == "D8_GSE174367") m <- m[!is.na(m$PMI), , drop = FALSE]
  if (nrow(m) < 4 || length(unique(m$disease)) < 2) stop("non-estimable region: ", cohort, " / ", region)
  d <- build_design(m, cohort, include_cell_count = include_cell_count)
  m <- d$meta
  counts <- obj$counts[, match(m$group_id, colnames(obj$counts)), drop = FALSE]
  dge <- DGEList(counts = counts)
  keep_gene <- filterByExpr(dge, design = d$design)
  dge <- dge[keep_gene, , keep.lib.sizes = FALSE]
  dge <- calcNormFactors(dge)
  v <- voom(dge, d$design, plot = FALSE)
  fit <- eBayes(lmFit(v, d$design))
  coef <- which(colnames(d$design) == "diseaseAD")
  if (length(coef) != 1) stop("disease coefficient missing for ", cohort, " / ", region)
  top <- topTable(fit, coef = coef, number = Inf, sort.by = "none", adjust.method = "BH")
  se <- sqrt(fit$s2.post) * fit$stdev.unscaled[, coef]
  gene_eff <- data.frame(
    cohort = cohort, region = region, cell_type = "Microglia", gene = rownames(top),
    n_groups = nrow(m), n_AD_donors = sum(m$disease == "AD"), n_Control_donors = sum(m$disease == "Control"),
    logFC = as.numeric(top$logFC), SE = as.numeric(se[rownames(top)]), t = as.numeric(top$t),
    P.Value = as.numeric(top$P.Value), FDR_within_region = as.numeric(top$adj.P.Val),
    model = paste0("limma_voom_donor_region_", paste(d$variables[-1], collapse = "_")),
    stringsAsFactors = FALSE
  )
  module_eff <- do.call(rbind, lapply(names(modules), function(module) {
    genes <- intersect(modules[[module]], rownames(v$E))
    if (!length(genes)) return(data.frame())
    score <- colMeans(v$E[genes, , drop = FALSE])
    mf <- eBayes(lmFit(matrix(score, nrow = 1, dimnames = list(module, colnames(v$E))), d$design))
    data.frame(cohort = cohort, region = region, module = module, n_present = length(genes), n_total = length(modules[[module]]),
      coverage = length(genes) / max(1, length(modules[[module]])), n_groups = nrow(m), n_AD_donors = sum(m$disease == "AD"),
      n_Control_donors = sum(m$disease == "Control"), logFC = as.numeric(mf$coefficients[1, coef]),
      SE = as.numeric(sqrt(mf$s2.post) * mf$stdev.unscaled[1, coef]), t = as.numeric(mf$t[1, coef]),
      P.Value = as.numeric(mf$p.value[1, coef]), model = paste0("limma_voom_module_", paste(d$variables[-1], collapse = "_")),
      stringsAsFactors = FALSE)
  }))
  list(gene = gene_eff, module = module_eff, meta = m)
}

run_all_regions <- function(include_cell_count = FALSE, complete_only = FALSE) {
  effects <- list(); module_effects <- list(); used_meta <- list(); k <- 0L
  for (cohort in names(inputs)) {
    obj <- read_pb(inputs[[cohort]]$count, inputs[[cohort]]$meta, cohort)
    for (region in sort(unique(obj$meta$region))) {
      k <- k + 1L
      ans <- fit_one_region(cohort, region, obj, include_cell_count, complete_only)
      effects[[k]] <- ans$gene
      module_effects[[k]] <- ans$module
      used_meta[[paste(cohort, region, sep = "__")]] <- ans$meta
    }
  }
  list(gene = do.call(rbind, effects), module = do.call(rbind, module_effects), used_meta = used_meta)
}

random_one <- function(y, se) {
  ok <- is.finite(y) & is.finite(se) & se > 0
  y <- y[ok]; se <- se[ok]; k <- length(y)
  if (k < 1) return(c(k = 0, beta = NA, se = NA, z = NA, p = NA, tau2 = NA, I2 = NA, Q_p = NA))
  w <- 1 / se^2
  fixed <- sum(w * y) / sum(w)
  if (k == 1) return(c(k = 1, beta = y, se = se, z = y / se, p = 2 * pnorm(-abs(y / se)), tau2 = 0, I2 = 0, Q_p = NA))
  Q <- sum(w * (y - fixed)^2)
  cterm <- sum(w) - sum(w^2) / sum(w)
  tau2 <- max(0, (Q - (k - 1)) / max(cterm, .Machine$double.eps))
  wr <- 1 / (se^2 + tau2)
  beta <- sum(wr * y) / sum(wr)
  se_beta <- sqrt(1 / sum(wr))
  z <- beta / se_beta
  p <- 2 * pnorm(-abs(z))
  I2 <- max(0, (Q - (k - 1)) / Q)
  c(k = k, beta = beta, se = se_beta, z = z, p = p, tau2 = tau2, I2 = I2, Q_p = pchisq(Q, df = k - 1, lower.tail = FALSE))
}

meta_by_gene <- function(effects, exclude_cohort = NULL) {
  if (is.null(exclude_cohort)) x <- effects else x <- effects[effects$cohort != exclude_cohort, , drop = FALSE]
  split_x <- split(x, x$gene)
  out <- lapply(split_x, function(z) {
    r <- random_one(z$logFC, z$SE)
    data.frame(gene = z$gene[1], k = r[["k"]], meta_logFC = r[["beta"]], meta_SE = r[["se"]], meta_z = r[["z"]], meta_P = r[["p"]], tau2 = r[["tau2"]], I2 = r[["I2"]], Q_p = r[["Q_p"]], stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, out); rownames(out) <- NULL
  if (is.null(dim(out))) stop("meta_by_gene produced a non-tabular object")
  out$meta_FDR <- p.adjust(out$meta_P, method = "BH")
  out
}

fixed_pool <- function(z) {
  r <- random_one(z$logFC, z$SE)
  data.frame(beta = r[["beta"]], se = r[["se"]], k = r[["k"]], direction = ifelse(r[["beta"]] > 0, "AD_up", ifelse(r[["beta"]] < 0, "AD_down", "zero")), stringsAsFactors = FALSE)
}

pool_cohorts <- function(effects) {
  genes <- sort(unique(effects$gene)); cohorts <- names(inputs)
  rows <- list(); j <- 0L
  for (g in genes) for (c in cohorts) {
    z <- effects[effects$gene == g & effects$cohort == c, , drop = FALSE]
    if (!nrow(z)) next
    q <- fixed_pool(z); j <- j + 1L
    rows[[j]] <- data.frame(gene = g, cohort = c, cohort_logFC = q$beta, cohort_SE = q$se, cohort_k_regions = q$k, cohort_direction = q$direction, stringsAsFactors = FALSE)
  }
  do.call(rbind, rows)
}

loo_table <- function(effects, pooled, pooled_meta) {
  rows <- list(); j <- 0L
  for (ex in names(inputs)) {
    x <- meta_by_gene(effects, ex)
    names(x)[names(x) == "meta_logFC"] <- "loo_logFC"
    names(x)[names(x) == "meta_SE"] <- "loo_SE"
    names(x)[names(x) == "meta_P"] <- "loo_P"
    names(x)[names(x) == "meta_FDR"] <- "loo_FDR"
    names(x)[names(x) == "k"] <- "loo_k"
    x$excluded_cohort <- ex
    base <- pooled_meta[, c("gene", "meta_logFC"), drop = FALSE]
    x <- merge(x, base, by = "gene", all.x = TRUE)
    x$direction_stable <- sign(x$loo_logFC) == sign(x$meta_logFC) & x$loo_logFC != 0 & x$meta_logFC != 0
    j <- j + 1L; rows[[j]] <- x
  }
  do.call(rbind, rows)
}

primary <- run_all_regions()
write.table(primary$gene, file.path(out_root, "effects", "primary_region_gene_effects.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(primary$module, file.path(out_root, "effects", "primary_region_module_effects.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
meta_gene <- meta_by_gene(primary$gene)
write.table(meta_gene, file.path(out_root, "meta", "random_effects_meta_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cohort_pooled <- pool_cohorts(primary$gene)
write.table(cohort_pooled, file.path(out_root, "meta", "cohort_pooled_effects.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
loo <- loo_table(primary$gene, cohort_pooled, meta_gene)
write.table(loo, file.path(out_root, "meta", "leave_one_cohort_out.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

meta_module <- do.call(rbind, lapply(split(primary$module, primary$module$module), function(z) {
  r <- random_one(z$logFC, z$SE)
  data.frame(module = z$module[1], k = r[["k"]], meta_logFC = r[["beta"]], meta_SE = r[["se"]], meta_z = r[["z"]], meta_P = r[["p"]], tau2 = r[["tau2"]], I2 = r[["I2"]], Q_p = r[["Q_p"]], n_present_min = min(z$n_present), coverage_min = min(z$coverage), stringsAsFactors = FALSE)
}))
meta_module$meta_FDR <- p.adjust(meta_module$meta_P, method = "BH")
write.table(meta_module, file.path(out_root, "meta", "random_effects_module_meta.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

complete <- run_all_regions(complete_only = TRUE)
cell_sens <- run_all_regions(include_cell_count = TRUE)
complete_meta <- meta_by_gene(complete$gene)
cell_meta <- meta_by_gene(cell_sens$gene)
write.table(complete$gene, file.path(out_root, "effects", "complete_covariate_region_gene_effects.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(cell_sens$gene, file.path(out_root, "effects", "cell_count_adjusted_region_gene_effects.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(complete_meta, file.path(out_root, "meta", "complete_covariate_meta.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(cell_meta, file.path(out_root, "meta", "cell_count_adjusted_meta.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

estimability <- do.call(rbind, lapply(names(inputs), function(cohort) {
  obj <- read_pb(inputs[[cohort]]$count, inputs[[cohort]]$meta, cohort)
  do.call(rbind, lapply(sort(unique(obj$meta$region)), function(region) {
    m <- obj$meta[obj$meta$region == region & obj$meta$n_cells_used >= 20, , drop = FALSE]
    data.frame(cohort = cohort, region = region, n_donors = nrow(m), n_AD_donors = sum(m$disease == "AD"), n_Control_donors = sum(m$disease == "Control"), min_cells_used = ifelse(nrow(m), min(m$n_cells_used), NA), estimable = nrow(m) >= 4 && sum(m$disease == "AD") >= 2 && sum(m$disease == "Control") >= 2, stringsAsFactors = FALSE)
  }))
}))
write.table(estimability, file.path(out_root, "gates", "primary_cohort_region_estimability.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

gate_genes <- subset(meta_gene, meta_FDR < 0.05)
gate_genes$concordant_cohorts <- NA_integer_
for (i in seq_len(nrow(gate_genes))) {
  q <- cohort_pooled[cohort_pooled$gene == gate_genes$gene[i], , drop = FALSE]
  gate_genes$concordant_cohorts[i] <- max(sum(q$cohort_direction == "AD_up"), sum(q$cohort_direction == "AD_down"))
}
gene_loo <- loo[loo$gene %in% gate_genes$gene, , drop = FALSE]
loo_stable <- if (nrow(gene_loo)) all(gene_loo$direction_stable) else FALSE
gene_signal <- nrow(gate_genes) > 0 && any(gate_genes$concordant_cohorts >= 2) && loo_stable
module_signal <- any(meta_module$meta_FDR < 0.05)

sensitivity <- data.frame(
  analysis = c("sex_adjusted_primary", "age_adjusted_primary", "RIN_adjusted_primary", "platform_handled_by_within_cohort_meta", "brain_region_stratified", "cell_proportion_adjusted", "incomplete_covariate_donor_exclusion", "influential_dataset_leave_one_out"),
  status = c("PASS_AS_MODEL_TERM", "PASS_AS_MODEL_TERM", "PASS_AS_MODEL_TERM", "PASS_BY_DESIGN", "PASS_BY_DESIGN", ifelse(any(cell_meta$meta_FDR < 0.05), "SIGNAL_RETAINED", "NO_SIGNAL"), ifelse(any(complete_meta$meta_FDR < 0.05), "SIGNAL_RETAINED", "NO_SIGNAL"), ifelse(loo_stable, "DIRECTION_STABLE", "DIRECTION_CHANGED")),
  detail = c("sex included where complete", "age included where complete", "RIN included where complete", "no cross-platform pooling; cohort-specific effects then meta", "EC/SSC, OC/OTC, and PFC analyzed separately", "added log1p(n_cells_used)", "D8 Sample-100 with missing PMI excluded in complete-covariate sensitivity", "all three primary cohorts excluded one at a time"),
  stringsAsFactors = FALSE
)
write.table(sensitivity, file.path(out_root, "meta", "sensitivity_audit.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

criteria <- data.frame(
  criterion = c("all_primary_cohorts_estimable", "random_effects_gene_FDR_lt_0.05", "concordant_direction_in_at_least_2_cohorts", "every_leave_one_cohort_out_direction_stable", "frozen_module_meta_FDR_lt_0.05"),
  observed = c(all(estimability$estimable), nrow(gate_genes), ifelse(nrow(gate_genes), max(gate_genes$concordant_cohorts), 0), loo_stable, module_signal),
  pass = c(all(estimability$estimable), nrow(gate_genes) > 0, ifelse(nrow(gate_genes), any(gate_genes$concordant_cohorts >= 2), FALSE), loo_stable, module_signal),
  notes = c("D3/D7/D8 region strata require >=2 AD and >=2 controls plus >=20 cells per donor-region", "BH-FDR<0.05 in random-effects gene meta", "cohort-level region-pooled direction", "LOO compares direction to pooled estimate", "frozen M1/M3/M5/M6/M7 module set"),
  stringsAsFactors = FALSE
)
write.table(criteria, file.path(out_root, "gates", "G1_criteria.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(gate_genes, file.path(out_root, "gates", "G1_gene_signal.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

cat(jsonlite::toJSON(list(run_id = run_id, n_primary_region_units = nrow(estimability), all_primary_cohorts_estimable = all(estimability$estimable), n_meta_genes = nrow(meta_gene), n_meta_genes_fdr_lt_0_05 = nrow(gate_genes), gene_signal = gene_signal, module_signal = module_signal, g1_signal = gene_signal || module_signal), auto_unbox = TRUE, pretty = TRUE), file = file.path(out_root, "gates", "g1_summary.json"))
cat("G1 v3.1 statistical execution completed\n")
