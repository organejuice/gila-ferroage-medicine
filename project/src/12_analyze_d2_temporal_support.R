suppressPackageStartupMessages({
  library(limma)
  library(edgeR)
})

parse_root <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) >= 1 && nzchar(args[[1]])) return(args[[1]])
  getwd()
}
root <- gsub("\\\\", "/", parse_root())
raw_path <- file.path(root, "data", "raw", "D2_GSE227221", "GSE227221_iMGLs.bulk.tpm.txt.gz")
meta_path <- file.path(root, "results", "tables", "D2_sample_metadata.tsv")
tables <- file.path(root, "results", "tables")
interim <- file.path(root, "data", "interim")
dir.create(tables, recursive = TRUE, showWarnings = FALSE)
dir.create(interim, recursive = TRUE, showWarnings = FALSE)

raw <- read.delim(gzfile(raw_path), check.names = FALSE, row.names = 1, stringsAsFactors = FALSE)
meta <- read.delim(meta_path, check.names = FALSE, stringsAsFactors = FALSE)
stopifnot(all(meta$sample_id_in_expression_table %in% colnames(raw)))
raw <- raw[, meta$sample_id_in_expression_table, drop = FALSE]
raw_ids <- sub("\\..*$", "", rownames(raw))
map_path <- file.path(root, "data", "interim", "D2_transcript_gene_map.tsv")
map <- read.delim(map_path, stringsAsFactors = FALSE, check.names = FALSE)
map <- map[nzchar(map$transcript_id) & nzchar(map$gene_symbol), c("transcript_id", "gene_symbol"), drop = FALSE]
map <- map[!duplicated(map$transcript_id), , drop = FALSE]
tx_to_gene <- setNames(map$gene_symbol, map$transcript_id)
keep <- raw_ids %in% names(tx_to_gene)
raw_keep <- raw[keep, , drop = FALSE]
gene_names <- unname(tx_to_gene[raw_ids[keep]])
gene_tpm <- rowsum(as.matrix(raw_keep), group = gene_names, reorder = FALSE)
storage.mode(gene_tpm) <- "numeric"
write.table(gene_tpm, file.path(tables, "D2_gene_tpm.tsv"), sep = "\t", quote = FALSE, col.names = NA)

map_audit <- data.frame(
  transcript_rows = nrow(raw),
  transcript_rows_mapped = sum(keep),
  transcript_mapping_rate = mean(keep),
  unique_gene_symbols = nrow(gene_tpm),
  samples = ncol(raw),
  stringsAsFactors = FALSE
)
write.table(map_audit, file.path(tables, "D2_transcript_gene_audit.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

meta$time_factor <- factor(paste0("h", meta$time_hours), levels = c("h0", "h4", "h24", "h96"))
design <- model.matrix(~ 0 + time_factor, data = meta)
colnames(design) <- sub("time_factor", "", colnames(design))
fit <- eBayes(limma::lmFit(log2(gene_tpm + 1), design))
contrast_names <- c(h4_vs_h0 = "h4-h0", h24_vs_h0 = "h24-h0", h96_vs_h0 = "h96-h0")
all_signatures <- list()
for (nm in names(contrast_names)) {
  contrast <- limma::makeContrasts(contrasts = contrast_names[[nm]], levels = design)
  fc <- limma::contrasts.fit(fit, contrast)
  fc <- limma::eBayes(fc)
  tt <- limma::topTable(fc, number = Inf, sort.by = "none")
  out <- data.frame(gene = rownames(tt), logFC = tt$logFC, t = tt$t, P.Value = tt$P.Value, FDR = tt$adj.P.Val, B = tt$B, comparison = nm, direction = ifelse(tt$logFC > 0, "A_beta_up", "A_beta_down"), n_samples = nrow(meta), n_replicates_per_time = 3, stringsAsFactors = FALSE, row.names = NULL)
  out <- out[order(out$FDR, -abs(out$t)), ]
  write.table(out, file.path(tables, paste0("D2_temporal_signature_", nm, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
  all_signatures[[nm]] <- out
}

module_sets <- list(
  M1_Glycolysis = c("HK1", "HK2", "PFKFB3", "PFKFB4", "ALDOA", "ALDOC", "GAPDH", "PGK1", "PGAM1", "ENO1", "ENO2", "PKM", "LDHA", "PDK1", "SLC2A1", "HIF1A", "PFKP"),
  M3_SASP_Senescence = c("CDKN1A", "CDKN2A", "IL6", "IL1B", "TNF", "CCL2", "CCL3", "CCL4", "CSF1", "CXCL10", "MMP9", "MMP14", "SERPINE1", "IGFBP3", "IGFBP7", "NFKB1", "RELA", "TGFB1"),
  M5_Iron_Uptake = c("TFRC", "SLC11A2", "STEAP3", "DMT1"),
  M5_Iron_Storage = c("FTH1", "FTL"),
  M5_Iron_Export = c("SLC40A1", "CP"),
  M5_Ferritinophagy = c("NCOA4", "ATG3", "ATG5", "ATG7", "MAP1LC3B", "SQSTM1"),
  M6_LPO_ACSL4 = c("ACSL4", "ALOX15", "ALOX5", "ALOX12", "NOX4", "GPX4", "LPCAT3", "SAT1"),
  M7_FerroAging = unique(read.delim(file.path(tables, "D5_ferroaging_reference.tsv"), stringsAsFactors = FALSE)$gene_symbol)
)
m2 <- read.delim(file.path(tables, "I2_regulatory_activity.tsv"), stringsAsFactors = FALSE)
m2 <- m2[is.finite(m2$delta_presence) & m2$delta_presence != 0 & nzchar(m2$human_one_to_one_symbols), , drop = FALSE]
module_rows <- list()
log_expr <- log2(gene_tpm + 1)
zscore <- t(apply(log_expr, 1, function(x) {
  s <- sd(x)
  if (!is.finite(s) || s == 0) return(rep(0, length(x)))
  (x - mean(x)) / s
}))
for (module in names(module_sets)) {
  genes <- intersect(module_sets[[module]], rownames(zscore))
  if (module == "M2_H4K12la") next
  if (length(genes) >= 2) {
    score <- colMeans(zscore[genes, , drop = FALSE])
    module_rows[[module]] <- data.frame(module = module, n_total = length(module_sets[[module]]), n_present = length(genes), coverage = length(genes) / length(module_sets[[module]]), score = score, stringsAsFactors = FALSE)
  }
}
if (nrow(m2) > 0) {
  m2_genes <- unique(unlist(strsplit(m2$human_one_to_one_symbols, ";", fixed = TRUE)))
  m2_genes <- intersect(m2_genes, rownames(zscore))
  weights <- setNames(m2$delta_presence[match(m2_genes, m2$human_one_to_one_symbols)], m2_genes)
  weights[!is.finite(weights)] <- 0
  if (length(m2_genes) >= 2 && sum(abs(weights)) > 0) {
    score <- colSums(sweep(zscore[m2_genes, , drop = FALSE], 1, weights[m2_genes], "*")) / sum(abs(weights[m2_genes]))
    module_rows[["M2_H4K12la"]] <- data.frame(module = "M2_H4K12la", n_total = nrow(m2), n_present = length(m2_genes), coverage = length(m2_genes) / max(1, nrow(m2)), score = score, stringsAsFactors = FALSE)
  }
}
module_effects <- list()
for (module in names(module_rows)) {
  row <- module_rows[[module]]
  score_mat <- matrix(row$score, nrow = 1, dimnames = list(module, meta$sample_id_in_expression_table))
  fit_m <- eBayes(limma::lmFit(score_mat, design))
  for (nm in names(contrast_names)) {
    contrast <- limma::makeContrasts(contrasts = contrast_names[[nm]], levels = design)
    cf <- eBayes(limma::contrasts.fit(fit_m, contrast))
    module_effects[[length(module_effects) + 1]] <- data.frame(module = module, comparison = nm, n_total = row$n_total[[1]], n_present = row$n_present[[1]], coverage = row$coverage[[1]], logFC = cf$coefficients[1, 1], t = cf$t[1, 1], P.Value = cf$p.value[1, 1], stringsAsFactors = FALSE)
  }
}
module_effects <- do.call(rbind, module_effects)
module_effects$FDR <- p.adjust(module_effects$P.Value, method = "BH")
module_effects$direction <- ifelse(module_effects$logFC > 0, "A_beta_up", "A_beta_down")
write.table(module_effects, file.path(tables, "D2_module_effects.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

provenance <- list(dataset = "D2_GSE227221", role = "iMGL A_beta_1_42 fibril time-course support", input_feature = "versioned Ensembl transcript TPM", transcript_to_gene = "Ensembl release 109 GRCh38 GTF; version suffix stripped; first transcript record retained", differential_model = "limma on log2(TPM+1); timepoint contrasts vs h0; 3 replicates per timepoint", module_score = "per-gene z-scored log2(TPM+1), donor-free time-course summary; M2 weighted by D1 delta_presence", boundaries = c("TPM is not raw count and is not used in count-based DESeq2/edgeR", "D2 does not replace human postmortem AD-Control discovery"), mapped_transcripts = sum(keep), total_transcripts = nrow(raw), mapped_genes = nrow(gene_tpm))
jsonlite::write_json(provenance, file.path(interim, "D2_temporal_provenance.json"), auto_unbox = TRUE, pretty = TRUE)
cat(sprintf("D2 temporal support completed: %d genes, %d samples\\n", nrow(gene_tpm), ncol(gene_tpm)))