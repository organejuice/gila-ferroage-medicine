options(stringsAsFactors = FALSE, width = 220)
suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
  library(yaml)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})

root <- "."
cfg <- yaml::read_yaml(file.path(root, "configs", "g5_model_scoped_v3.2_cell_rescue_execution.yaml"))
run_dir <- file.path(root, "results", cfg$run_id)
dir.create(file.path(run_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(run_dir, "gates"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(run_dir, "registries"), recursive = TRUE, showWarnings = FALSE)

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
raw_meta <- read.csv(gzfile(meta_file), stringsAsFactors = FALSE, check.names = FALSE)
names(raw_meta)[1] <- "Sample"
if (!"Sample" %in% names(raw_meta)) stop("GSE319380 metadata lacks Sample column after BOM normalization")
raw_counts <- read.csv(gzfile(counts_file), stringsAsFactors = FALSE, check.names = FALSE, row.names = 1)
raw_counts <- as.matrix(raw_counts[, raw_meta$Sample, drop = FALSE])
storage.mode(raw_counts) <- "numeric"
ens <- sub("\\..*$", "", rownames(raw_counts))
symbols <- AnnotationDbi::mapIds(org.Hs.eg.db, keys = ens, keytype = "ENSEMBL", column = "SYMBOL", multiVals = "first")
mapping <- data.frame(ensembl = ens, symbol = unname(symbols), stringsAsFactors = FALSE)
mapping <- mapping[!is.na(mapping$symbol) & nzchar(mapping$symbol), , drop = FALSE]
mapped_counts <- rowsum(raw_counts[match(mapping$ensembl, rownames(raw_counts)), , drop = FALSE], group = mapping$symbol, reorder = FALSE)
rownames(mapped_counts) <- toupper(rownames(mapped_counts))

meta <- raw_meta[match(colnames(raw_counts), raw_meta$Sample), , drop = FALSE]
meta$disease <- ifelse(meta$Type == "sAD", "sAD", "control")
meta$Line <- factor(meta$Line)
meta$Stressor <- factor(meta$Stressor, levels = c("CTR", "AB", "H2O2"))
meta$unit_id <- paste(meta$Line, meta$Stressor, sep = "__")
unit_ids <- unique(meta$unit_id)
unit_counts <- do.call(cbind, lapply(unit_ids, function(u) rowSums(mapped_counts[, meta$unit_id == u, drop = FALSE])))
colnames(unit_counts) <- unit_ids
unit_meta <- meta[match(unit_ids, meta$unit_id), c("unit_id", "Line", "Stressor", "disease"), drop = FALSE]
rownames(unit_meta) <- unit_meta$unit_id
unit_meta$Line <- droplevels(unit_meta$Line)
unit_meta$Stressor <- droplevels(unit_meta$Stressor)

design <- model.matrix(~ disease + Stressor, data = unit_meta)
dge <- DGEList(counts = unit_counts)
keep <- filterByExpr(dge, design = design)
dge <- calcNormFactors(dge[keep, , keep.lib.sizes = FALSE])
v0 <- voom(dge, design, plot = FALSE)
corfit <- tryCatch(duplicateCorrelation(v0, design = design, block = unit_meta$Line), error = function(e) list(consensus = 0, error = conditionMessage(e)))
rho <- as.numeric(corfit$consensus)
if (!is.finite(rho)) rho <- 0
v <- voom(dge, design, plot = FALSE, block = unit_meta$Line, correlation = rho)
coef_idx <- which(colnames(design) == "diseasesAD")
if (length(coef_idx) != 1) stop("Disease coefficient is not estimable")

fit_score <- function(module, genes) {
  present <- intersect(toupper(genes), rownames(v$E))
  if (length(present) < 2) return(NULL)
  score <- matrix(colMeans(v$E[present, , drop = FALSE]), nrow = 1, dimnames = list(module, colnames(v$E)))
  fit <- lmFit(score, design, block = unit_meta$Line, correlation = rho)
  beta <- as.numeric(fit$coefficients[1, coef_idx])
  se <- as.numeric(sqrt(fit$sigma[1]^2) * fit$stdev.unscaled[1, coef_idx])
  tval <- beta / se
  df <- as.numeric(fit$df.residual[1])
  pval <- 2 * pt(-abs(tval), df = df)
  data.frame(
    dataset = "GSE319380", model = "sAD_iPSC_astrocyte_vs_control", route = "cell_model",
    contrast = "disease_adjusted_for_stressor", module = module,
    n_present = length(present), n_total = length(genes), coverage = length(present) / length(genes),
    n_lines = length(unique(unit_meta$Line)), n_units = nrow(unit_meta),
    n_sAD_lines = sum(unique(unit_meta[, c("Line", "disease")])$disease == "sAD"),
    n_control_lines = sum(unique(unit_meta[, c("Line", "disease")])$disease == "control"),
    logFC = beta, SE = se, t = tval, P.Value = pval, model_formula = "~ disease + Stressor; block=Line",
    stringsAsFactors = FALSE
  )
}

effects <- do.call(rbind, lapply(names(modules), function(nm) fit_score(nm, modules[[nm]])))
effects$FDR_within_contrast <- p.adjust(effects$P.Value, method = "BH")
effects$direction <- ifelse(effects$logFC > 0, "sAD_up", "sAD_down")
write.table(effects, file.path(run_dir, "tables", "GSE319380_cell_module_effects.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(unit_meta, file.path(run_dir, "tables", "GSE319380_cell_units_used.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(mapping, file.path(run_dir, "registries", "GSE319380_ensembl_symbol_mapping.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(metric = c("raw_genes", "mapped_genes", "tested_genes", "units", "lines", "duplicate_correlation"), value = c(nrow(raw_counts), nrow(mapped_counts), nrow(dge), nrow(unit_meta), length(unique(unit_meta$Line)), rho)), file.path(run_dir, "tables", "GSE319380_cell_qc_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(coefficient = colnames(design), estimable = qr(design)$rank == ncol(design)), file.path(run_dir, "tables", "GSE319380_cell_design.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("wrote %d cell module rows; duplicateCorrelation=%.4f\n", nrow(effects), rho))
