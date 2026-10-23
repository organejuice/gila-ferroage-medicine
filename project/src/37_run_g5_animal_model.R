#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, width = 220)
suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
  library(yaml)
})

root <- "."
cfg <- yaml::read_yaml(file.path(root, "configs", "g5_model_scoped_v3.2_execution.yaml"))
run_dir <- file.path(root, "results", cfg$run_id)
if (dir.exists(run_dir)) stop("Refusing to overwrite existing run directory: ", run_dir)
dir.create(file.path(run_dir, "tables"), recursive = TRUE)
dir.create(file.path(run_dir, "gates"), recursive = TRUE)

modules <- list(
  M3_SASP_Senescence = c("CDKN1A", "CDKN2A", "IL6", "IL1B", "TNF", "CCL2", "CCL3", "CCL4", "CSF1", "CXCL10", "MMP9", "MMP14", "SERPINE1", "IGFBP3", "IGFBP7", "NFKB1", "RELA", "TGFB1"),
  M5_Iron_Uptake = c("TFRC", "SLC11A2", "STEAP3", "DMT1"),
  M5_Iron_Storage = c("FTH1", "FTL"),
  M5_Iron_Export = c("SLC40A1", "CP"),
  M5_Ferritinophagy = c("NCOA4", "ATG3", "ATG5", "ATG7", "MAP1LC3B", "SQSTM1"),
  M6_LPO_ACSL4 = c("ACSL4", "ALOX15", "ALOX5", "ALOX12", "NOX4", "GPX4", "LPCAT3", "SAT1"),
  M7_FerroAging_Author = unique(as.character(read.delim(file.path(root, "results", "archive", "2026-09-20-pre-v3-freeze", "results", "tables", "D5_ferroaging_reference.tsv"), sep = "\t", stringsAsFactors = FALSE)$gene_symbol))
)

counts_file <- file.path(root, cfg$inputs$counts)
meta_file <- file.path(root, cfg$inputs$metadata)
count_df <- read.delim(gzfile(counts_file), check.names = FALSE, stringsAsFactors = FALSE)
rownames(count_df) <- make.unique(toupper(as.character(count_df[[1]])))
counts <- as.matrix(count_df[, -1, drop = FALSE])
storage.mode(counts) <- "numeric"
meta <- read.delim(meta_file, check.names = FALSE, stringsAsFactors = FALSE)
meta <- meta[match(colnames(counts), meta$group_id), , drop = FALSE]
if (anyNA(meta$group_id)) stop("metadata/count mismatch")
meta$disease <- factor(meta$disease, levels = c("WT_control", "AD_model"))
meta$region <- factor(meta$region)
meta$age_months <- as.numeric(meta$age_months)

fit_stratum <- function(m, stratum_type, stratum_name) {
  m <- m[m$n_cells_used >= cfg$analysis$minimum_cells_per_animal, , drop = FALSE]
  m <- m[complete.cases(m[, c("disease", "age_months")]), , drop = FALSE]
  if (nrow(m) < 4 || length(unique(m$disease)) < 2) return(NULL)
  vars <- character(0)
  if (stratum_type == "region" && length(unique(m$age_months)) > 1) {
    m$age_z <- as.numeric(scale(m$age_months)); vars <- "age_z"
  }
  formula_text <- paste(c("~ disease", vars), collapse = " + ")
  design <- model.matrix(as.formula(formula_text), data = m)
  if (qr(design)$rank < ncol(design) || !"diseaseAD_model" %in% colnames(design)) return(NULL)
  sub_counts <- counts[, match(m$group_id, colnames(counts)), drop = FALSE]
  dge <- DGEList(counts = sub_counts)
  keep <- filterByExpr(dge, design = design)
  dge <- calcNormFactors(dge[keep, , keep.lib.sizes = FALSE])
  v <- voom(dge, design, plot = FALSE)
  coef <- which(colnames(design) == "diseaseAD_model")
  rows <- lapply(names(modules), function(module) {
    genes <- intersect(toupper(modules[[module]]), rownames(v$E))
    if (length(genes) < 2) return(NULL)
    score <- colMeans(v$E[genes, , drop = FALSE])
    fit <- eBayes(lmFit(matrix(score, nrow = 1, dimnames = list(module, colnames(v$E))), design))
    data.frame(
      dataset = "GSE143758", model = "5xFAD_vs_WT", region = as.character(stratum_name),
      stratum_type = stratum_type, stratum = as.character(stratum_name), module = module,
      n_present = length(genes), n_total = length(modules[[module]]), coverage = length(genes) / length(modules[[module]]),
      n_animals = nrow(m), n_AD_model = sum(m$disease == "AD_model"), n_WT_control = sum(m$disease == "WT_control"),
      logFC = as.numeric(fit$coefficients[1, coef]),
      SE = as.numeric(sqrt(fit$s2.post) * fit$stdev.unscaled[1, coef]),
      t = as.numeric(fit$t[1, coef]), P.Value = as.numeric(fit$p.value[1, coef]),
      model_formula = formula_text, stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

out <- list(); idx <- 0L
for (region in levels(meta$region)) {
  ans <- fit_stratum(meta[meta$region == region, , drop = FALSE], "region", region)
  if (!is.null(ans)) { idx <- idx + 1L; out[[idx]] <- ans }
}
effects <- do.call(rbind, out)
effects$FDR_within_stratum <- ave(effects$P.Value, effects$stratum_type, effects$stratum, FUN = function(x) p.adjust(x, method = "BH"))
effects$direction <- ifelse(effects$logFC > 0, "AD_model_up", "AD_model_down")
write.table(effects, file.path(run_dir, "tables", "GSE143758_animal_region_module_effects.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(meta, file.path(run_dir, "tables", "GSE143758_animal_metadata_used.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(region = levels(meta$region), n_animals = as.integer(table(meta$region)), n_AD_model = as.integer(tapply(meta$disease == "AD_model", meta$region, sum)), n_WT_control = as.integer(tapply(meta$disease == "WT_control", meta$region, sum))), file.path(run_dir, "tables", "GSE143758_animal_estimability.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("wrote %d animal region-module rows\n", nrow(effects)))
