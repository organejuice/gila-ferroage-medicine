suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
})

parse_root <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) >= 1 && nzchar(args[[1]])) return(args[[1]])
  getwd()
}
root <- gsub("\\\\", "/", parse_root())
tables <- file.path(root, "results", "tables")
out_path <- file.path(tables, "D3_module_effects_all.tsv")

pb <- read.delim(file.path(tables, "D3_pseudobulk_counts.tsv"), row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
meta <- read.delim(file.path(tables, "D3_pseudobulk_metadata.tsv"), check.names = FALSE, stringsAsFactors = FALSE)
ref <- read.delim(file.path(tables, "D5_ferroaging_reference.tsv"), check.names = FALSE, stringsAsFactors = FALSE)
i2 <- read.delim(file.path(tables, "I2_regulatory_activity.tsv"), check.names = FALSE, stringsAsFactors = FALSE)
stopifnot(all(meta$group_id %in% colnames(pb)))

module_sets <- list(
  M1 = c("HK1", "HK2", "PFKFB3", "PFKFB4", "ALDOA", "ALDOC", "GAPDH", "PGK1", "PGAM1", "ENO1", "ENO2", "PKM", "LDHA", "PDK1", "SLC2A1", "HIF1A", "PFKP"),
  M3 = c("CDKN1A", "CDKN2A", "IL6", "IL1B", "TNF", "CCL2", "CCL3", "CCL4", "CSF1", "CXCL10", "MMP9", "MMP14", "SERPINE1", "IGFBP3", "IGFBP7", "NFKB1", "RELA", "TGFB1"),
  M5_Iron_Uptake = c("TFRC", "SLC11A2", "STEAP3", "DMT1"),
  M5_Iron_Storage = c("FTH1", "FTL"),
  M5_Iron_Export = c("SLC40A1", "CP"),
  M5_Ferritinophagy = c("NCOA4", "ATG3", "ATG5", "ATG7", "MAP1LC3B", "SQSTM1"),
  M6 = c("ACSL4", "ALOX15", "ALOX5", "ALOX12", "NOX4", "GPX4", "LPCAT3", "SAT1"),
  M7 = unique(ref$gene_symbol[ref$module == "M7"])
)

counts <- as.matrix(pb)
storage.mode(counts) <- "numeric"
logcpm <- edgeR::cpm(edgeR::DGEList(counts = round(counts)), log = TRUE, prior.count = 2)
gene_upper <- toupper(rownames(logcpm))
rownames(logcpm) <- gene_upper
logcpm <- logcpm[!duplicated(rownames(logcpm)), , drop = FALSE]

# M2 is a signed, one-to-one human orthology projection of the mouse microglia D1 peak profile.
i2 <- i2[i2$context == "M" & is.finite(i2$delta_presence) & i2$delta_presence != 0 & !is.na(i2$human_one_to_one_symbols) & nzchar(i2$human_one_to_one_symbols), , drop = FALSE]
m2_weights <- tapply(i2$delta_presence, toupper(i2$human_one_to_one_symbols), sum)
module_sets$M2 <- names(m2_weights)

module_rows <- list()
row_id <- 0
for (ct in c("Microglia", "Astrocyte")) {
  idx <- which(meta$cell_type == ct)
  m <- meta[idx, , drop = FALSE]
  m$disease <- factor(m$disease, levels = c("Control", "AD"))
  m$region <- factor(m$region)
  m$sex <- factor(m$sex)
  m$age_c <- as.numeric(m$age) - mean(as.numeric(m$age), na.rm = TRUE)
  m$rin_c <- as.numeric(m$rin) - mean(as.numeric(m$rin), na.rm = TRUE)
  design <- model.matrix(~ disease + region + age_c + sex + rin_c, data = m)
  for (module in names(module_sets)) {
    requested <- toupper(unique(module_sets[[module]]))
    present <- intersect(requested, rownames(logcpm))
    if (length(present) < 2) next
    if (module == "M2") {
      weights <- as.numeric(m2_weights[present])
      names(weights) <- present
      score <- colSums(sweep(logcpm[present, m$group_id, drop = FALSE], 1, weights, "*")) / sum(abs(weights))
      source_type <- "D1_mouse_M_one_to_one_human_signed_projection"
    } else {
      score <- colMeans(logcpm[present, m$group_id, drop = FALSE], na.rm = TRUE)
      source_type <- ifelse(module == "M7", "D5_author_reference", "project_operational_gene_set")
    }
    score <- as.numeric(score)
    score_mat <- matrix(score, nrow = 1, dimnames = list(module, m$group_id))
    corfit <- tryCatch(limma::duplicateCorrelation(score_mat, design, block = m$donor_id), error = function(e) NULL)
    if (!is.null(corfit) && is.finite(corfit$consensus)) {
      fit <- limma::lmFit(score_mat, design, block = m$donor_id, correlation = corfit$consensus)
    } else {
      fit <- limma::lmFit(score_mat, design)
    }
    fit <- limma::eBayes(fit)
    row_id <- row_id + 1
    module_rows[[row_id]] <- data.frame(
      module = module,
      cell_type = ct,
      source_type = source_type,
      n_total = length(requested),
      n_present = length(present),
      coverage = length(present) / length(requested),
      n_groups = nrow(m),
      n_AD_groups = sum(m$disease == "AD"),
      n_Control_groups = sum(m$disease == "Control"),
      logFC = unname(fit$coefficients[1, "diseaseAD"]),
      t = unname(fit$t[1, "diseaseAD"]),
      P.Value = unname(fit$p.value[1, "diseaseAD"]),
      stringsAsFactors = FALSE
    )
  }
}
res <- do.call(rbind, module_rows)
res$FDR <- p.adjust(res$P.Value, method = "BH")
res$direction <- ifelse(res$logFC > 0, "AD_up", "AD_down")
write.table(res, out_path, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("D3 all-module effects written: %d rows\n", nrow(res)))