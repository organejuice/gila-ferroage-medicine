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

pb_path <- file.path(root, "results", "tables", "D3_pseudobulk_counts.tsv")
meta_path <- file.path(root, "results", "tables", "D3_pseudobulk_metadata.tsv")
ref_path <- file.path(root, "results", "tables", "D5_ferroaging_reference.tsv")
out_path <- file.path(root, "results", "tables", "D3_module_effects.tsv")

pb <- read.delim(pb_path, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
meta <- read.delim(meta_path, check.names = FALSE, stringsAsFactors = FALSE)
ref <- read.delim(ref_path, check.names = FALSE, stringsAsFactors = FALSE)
stopifnot(all(meta$group_id %in% colnames(pb)))

counts <- as.matrix(pb)
storage.mode(counts) <- "numeric"
logcpm <- edgeR::cpm(edgeR::DGEList(counts = round(counts)), log = TRUE, prior.count = 2)

module_rows <- list()
for (ct in c("Microglia", "Astrocyte")) {
  genes <- unique(ref$gene_symbol[ref$module == "M7"])
  genes <- intersect(genes, rownames(logcpm))
  idx <- which(meta$cell_type == ct)
  m <- meta[idx, , drop = FALSE]
  m$disease <- factor(m$disease, levels = c("Control", "AD"))
  m$region <- factor(m$region)
  m$sex <- factor(m$sex)
  m$age_c <- as.numeric(m$age) - mean(as.numeric(m$age), na.rm = TRUE)
  m$rin_c <- as.numeric(m$rin) - mean(as.numeric(m$rin), na.rm = TRUE)
  score <- colMeans(logcpm[genes, m$group_id, drop = FALSE], na.rm = TRUE)
  score_mat <- matrix(as.numeric(score), nrow = 1, dimnames = list("M7_score", m$group_id))
  design <- model.matrix(~ disease + region + age_c + sex + rin_c, data = m)
  block <- m$donor_id
  corfit <- tryCatch(limma::duplicateCorrelation(score_mat, design, block = block), error = function(e) NULL)
  if (!is.null(corfit) && is.finite(corfit$consensus)) {
    fit <- limma::lmFit(score_mat, design, block = block, correlation = corfit$consensus)
  } else {
    fit <- limma::lmFit(score_mat, design)
  }
  fit <- limma::eBayes(fit)
  coef_name <- "diseaseAD"
  ad_score <- mean(score[m$disease == "AD"])
  control_score <- mean(score[m$disease == "Control"])
  module_rows[[ct]] <- data.frame(
    module = "M7",
    cell_type = ct,
    n_genes_reference = length(unique(ref$gene_symbol[ref$module == "M7"])),
    n_genes_present = length(genes),
    n_groups = nrow(m),
    n_AD_groups = sum(m$disease == "AD"),
    n_Control_groups = sum(m$disease == "Control"),
    AD_mean_logCPM = ad_score,
    Control_mean_logCPM = control_score,
    logFC = unname(fit$coefficients[1, coef_name]),
    t = unname(fit$t[1, coef_name]),
    P.Value = unname(fit$p.value[1, coef_name]),
    stringsAsFactors = FALSE
  )
}
res <- do.call(rbind, module_rows)
res$FDR <- p.adjust(res$P.Value, method = "BH")
res$direction <- ifelse(res$logFC > 0, "AD_up", "AD_down")
write.table(res, out_path, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("D3 M7 module effects written: %s\\n", out_path))