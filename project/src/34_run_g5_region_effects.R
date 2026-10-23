#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, width = 200)
suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
  library(yaml)
})

root <- "."
cfg <- yaml::read_yaml(file.path(root, "configs", "g5_v3.2_execution.yaml"))
run_dir <- file.path(root, "results", cfg$run_id)
if (dir.exists(run_dir)) stop("Refusing to overwrite existing run directory: ", run_dir)
dir.create(file.path(run_dir, "tables"), recursive = TRUE)
dir.create(file.path(run_dir, "gates"), recursive = TRUE)
dir.create(file.path(run_dir, "registries"), recursive = TRUE)

modules <- list(
  M3_SASP_Senescence = c("CDKN1A", "CDKN2A", "IL6", "IL1B", "TNF", "CCL2", "CCL3", "CCL4", "CSF1", "CXCL10", "MMP9", "MMP14", "SERPINE1", "IGFBP3", "IGFBP7", "NFKB1", "RELA", "TGFB1"),
  M5_Iron_Uptake = c("TFRC", "SLC11A2", "STEAP3", "DMT1"),
  M5_Iron_Storage = c("FTH1", "FTL"),
  M5_Iron_Export = c("SLC40A1", "CP"),
  M5_Ferritinophagy = c("NCOA4", "ATG3", "ATG5", "ATG7", "MAP1LC3B", "SQSTM1"),
  M6_LPO_ACSL4 = c("ACSL4", "ALOX15", "ALOX5", "ALOX12", "NOX4", "GPX4", "LPCAT3", "SAT1"),
  M7_FerroAging_Author = unique(as.character(read.delim(file.path(root, cfg$inputs$d5_reference), sep = "\t", stringsAsFactors = FALSE)$gene_symbol))
)

read_pb <- function(count_file, meta_file, cohort) {
  x <- read.delim(gzfile(count_file), check.names = FALSE, stringsAsFactors = FALSE)
  rownames(x) <- make.unique(toupper(as.character(x[[1]])))
  counts <- as.matrix(x[, -1, drop = FALSE]); storage.mode(counts) <- "numeric"
  meta <- read.delim(meta_file, check.names = FALSE, stringsAsFactors = FALSE)
  meta <- meta[match(colnames(counts), meta$group_id), , drop = FALSE]
  if (anyNA(meta$group_id)) stop("metadata/count mismatch for ", cohort)
  if (!"n_cells_used" %in% names(meta) && "n_cells" %in% names(meta)) meta$n_cells_used <- meta$n_cells
  if (!"RIN" %in% names(meta) && "rin" %in% names(meta)) meta$RIN <- meta$rin
  meta$cohort <- cohort
  meta$disease <- factor(meta$disease, levels = c("Control", "AD"))
  meta$sex <- factor(meta$sex)
  if ("batch" %in% names(meta)) meta$batch <- factor(meta$batch)
  list(counts = counts, meta = meta)
}

fit_region <- function(cohort, region, obj) {
  m <- obj$meta[obj$meta$region == region & obj$meta$n_cells_used >= cfg$analysis$minimum_cells_per_donor_region, , drop = FALSE]
  if (nrow(m) < 4 || length(unique(m$disease)) < 2) return(NULL)
  vars <- intersect(c("age", "sex", "RIN"), names(m))
  m <- m[complete.cases(m[, vars, drop = FALSE]), , drop = FALSE]
  if (nrow(m) < 4 || length(unique(m$disease)) < 2) return(NULL)
  if ("age" %in% vars && sd(m$age) > 0) { m$age_z <- as.numeric(scale(m$age)); vars[vars == "age"] <- "age_z" }
  if ("RIN" %in% vars && sd(m$RIN) > 0) { m$RIN_z <- as.numeric(scale(m$RIN)); vars[vars == "RIN"] <- "RIN_z" }
  if (cohort == "D8_GSE174367" && "batch" %in% names(m) && all(!is.na(m$batch)) && nlevels(droplevels(m$batch)) > 1) vars <- c(vars, "batch")
  design <- model.matrix(as.formula(paste("~", paste(c("disease", vars), collapse = " + "))), data = m)
  if (qr(design)$rank < ncol(design) || !"diseaseAD" %in% colnames(design)) return(NULL)
  counts <- obj$counts[, match(m$group_id, colnames(obj$counts)), drop = FALSE]
  dge <- DGEList(counts = counts); keep <- filterByExpr(dge, design = design)
  dge <- calcNormFactors(dge[keep, , keep.lib.sizes = FALSE]); v <- voom(dge, design, plot = FALSE)
  coef <- which(colnames(design) == "diseaseAD")
  rows <- lapply(names(modules), function(module) {
    genes <- intersect(toupper(modules[[module]]), rownames(v$E))
    if (length(genes) < 2) return(NULL)
    score <- colMeans(v$E[genes, , drop = FALSE])
    mf <- eBayes(lmFit(matrix(score, nrow = 1, dimnames = list(module, colnames(v$E))), design))
    data.frame(cohort = cohort, region = region, module = module, n_present = length(genes), n_total = length(modules[[module]]), coverage = length(genes) / length(modules[[module]]), n_groups = nrow(m), n_AD_donors = sum(m$disease == "AD"), n_Control_donors = sum(m$disease == "Control"), logFC = as.numeric(mf$coefficients[1, coef]), SE = as.numeric(sqrt(mf$s2.post) * mf$stdev.unscaled[1, coef]), t = as.numeric(mf$t[1, coef]), P.Value = as.numeric(mf$p.value[1, coef]), model = paste0("limma_voom_astrocyte_region_", paste(c("disease", vars), collapse = "_")), stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

inputs <- list(
  D3_GSE160936 = read_pb(file.path(root, cfg$inputs$d3_counts), file.path(root, cfg$inputs$d3_metadata), "D3_GSE160936"),
  D7_GSE148822 = read_pb(file.path(root, cfg$inputs$d7_counts), file.path(root, cfg$inputs$d7_metadata), "D7_GSE148822"),
  D8_GSE174367 = read_pb(file.path(root, cfg$inputs$d8_counts), file.path(root, cfg$inputs$d8_metadata), "D8_GSE174367")
)

out <- list(); idx <- 0L
for (cohort in names(inputs)) for (region in sort(unique(inputs[[cohort]]$meta$region))) {
  ans <- fit_region(cohort, region, inputs[[cohort]])
  if (!is.null(ans)) { idx <- idx + 1L; out[[idx]] <- ans }
}
if (!length(out)) stop("No estimable astrocyte region effects")
effects <- do.call(rbind, out)
effects$FDR_within_stratum <- ave(effects$P.Value, effects$cohort, effects$region, FUN = function(x) p.adjust(x, method = "BH"))
effects$direction <- ifelse(effects$logFC > 0, "AD_up", "AD_down")
write.table(effects, file.path(run_dir, "tables", "astrocyte_region_module_effects.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("wrote %d region-module rows\n", nrow(effects)))
