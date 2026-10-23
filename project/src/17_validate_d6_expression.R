suppressPackageStartupMessages({
  library(rhdf5)
  library(Matrix)
  library(edgeR)
  library(limma)
  library(jsonlite)
})

parse_root <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) >= 1 && nzchar(args[[1]])) return(args[[1]])
  getwd()
}
root <- gsub("\\\\", "/", parse_root())
tables <- file.path(root, "results", "tables")
interim <- file.path(root, "data", "interim")
h5_path <- Sys.getenv("D6_H5_PATH", unset = file.path(root, "data", "raw", "D6_GSE214979", "GSE214979_filtered_feature_bc_matrix.h5"))
meta_path <- file.path(root, "data", "raw", "D6_GSE214979", "GSE214979_cell_metadata.csv.gz")
expected_bytes <- 1369492123
actual_bytes <- file.info(h5_path)$size
if (!is.finite(actual_bytes) || actual_bytes != expected_bytes) stop(sprintf("D6 H5 incomplete: %s/%s bytes", actual_bytes, expected_bytes))
con <- file(h5_path, "rb"); magic <- readBin(con, "raw", 8); close(con)
if (!identical(as.integer(magic), as.integer(charToRaw("\x89HDF\r\n\x1a\n")))) stop("D6 H5 header validation failed")

meta <- read.csv(gzfile(meta_path), check.names = FALSE, stringsAsFactors = FALSE)
rownames(meta) <- meta[[1]]
meta[[1]] <- NULL
meta$cell_type <- ifelse(meta$predicted.id == "Microglia", "Microglia", ifelse(meta$predicted.id == "Astrocytes", "Astrocyte", NA))
meta <- meta[!is.na(meta$cell_type), , drop = FALSE]
meta$Status <- factor(meta$Status, levels = c("Ctrl", "AD"))
meta$Sex <- factor(meta$Sex)
meta$structure <- factor(meta$structure)
meta$Age_c <- as.numeric(meta$Age) - mean(as.numeric(meta$Age), na.rm = TRUE)
meta$PMI_c <- as.numeric(meta$PMI..h.) - mean(as.numeric(meta$PMI..h.), na.rm = TRUE)
meta$RIN_c <- as.numeric(meta$RIN) - mean(as.numeric(meta$RIN), na.rm = TRUE)

read_10x_gene_expression <- function(path, target_barcodes, chunk_cols = 256L) {
  barcodes <- as.character(rhdf5::h5read(path, "matrix/barcodes"))
  feature_type <- as.character(rhdf5::h5read(path, "matrix/features/feature_type"))
  feature_names <- make.unique(as.character(rhdf5::h5read(path, "matrix/features/name")))
  rna_idx <- which(feature_type == "Gene Expression")
  if (length(rna_idx) < 1000) stop("D6 H5 has too few Gene Expression features")
  indptr <- as.numeric(rhdf5::h5read(path, "matrix/indptr"))
  shape <- as.numeric(rhdf5::h5read(path, "matrix/shape"))
  target_idx <- match(target_barcodes, barcodes)
  target_idx <- sort(target_idx[is.finite(target_idx)])
  if (length(target_idx) < 1000) stop("Too few D6 metadata-expression barcode matches")
  chunks <- list()
  chunk_barcodes <- character()
  chunk_id <- 0L
  for (c_start in seq.int(1L, length(barcodes), by = chunk_cols)) {
    c_end <- min(length(barcodes), c_start + chunk_cols - 1L)
    local_cols <- seq.int(c_start, c_end)
    target_local <- which(local_cols %in% target_idx)
    if (!length(target_local)) next
    entry_start <- indptr[c_start] + 1
    entry_end <- indptr[c_end + 1]
    if (entry_end < entry_start) {
      chunk <- Matrix::Matrix(0, nrow = length(rna_idx), ncol = length(target_local), sparse = TRUE)
    } else {
      entry_idx <- as.integer(entry_start):as.integer(entry_end)
      data_chunk <- as.numeric(rhdf5::h5read(path, "matrix/data", index = list(entry_idx)))
      index_chunk <- as.integer(rhdf5::h5read(path, "matrix/indices", index = list(entry_idx)))
      column_ids <- rep(seq_len(length(local_cols)), diff(indptr[c_start:(c_end + 1)]))
      keep_entry <- (index_chunk + 1L) %in% rna_idx & column_ids %in% target_local
      if (any(keep_entry)) {
        chunk <- Matrix::sparseMatrix(
          i = match(index_chunk[keep_entry] + 1L, rna_idx),
          j = match(column_ids[keep_entry], target_local),
          x = data_chunk[keep_entry],
          dims = c(length(rna_idx), length(target_local)),
          dimnames = list(feature_names[rna_idx], barcodes[local_cols[target_local]])
        )
      } else {
        chunk <- Matrix::Matrix(0, nrow = length(rna_idx), ncol = length(target_local), sparse = TRUE,
                                dimnames = list(feature_names[rna_idx], barcodes[local_cols[target_local]]))
      }
    }
    chunk_id <- chunk_id + 1L
    chunks[[chunk_id]] <- chunk
    chunk_barcodes <- c(chunk_barcodes, barcodes[local_cols[target_local]])
  }
  if (!length(chunks)) stop("No target cells were read from D6 H5")
  x <- do.call(cbind, chunks)
  colnames(x) <- chunk_barcodes
  x
}

x <- read_10x_gene_expression(h5_path, rownames(meta))
keep_cells <- intersect(colnames(x), rownames(meta))
if (length(keep_cells) < 1000) stop("Too few D6 metadata-expression barcode matches")
x <- x[, keep_cells, drop = FALSE]
meta <- meta[keep_cells, , drop = FALSE]
meta$group_id <- paste(meta$id, meta$structure, meta$cell_type, sep = "__")

# Use donor x structure x cell type units; retain groups with at least 20 cells and record this threshold.
group_n <- table(meta$group_id)
keep_groups <- names(group_n[group_n >= 20])
keep <- meta$group_id %in% keep_groups
x <- x[, keep, drop = FALSE]
meta <- meta[keep, , drop = FALSE]
group_n <- table(meta$group_id)
agg <- Matrix::sparse.model.matrix(~ 0 + factor(meta$group_id))
colnames(agg) <- sub("^factor\\(meta\\$group_id\\)", "", colnames(agg))
pb <- x %*% agg
meta_group <- meta[match(colnames(pb), meta$group_id), , drop = FALSE]
meta_group <- meta_group[!duplicated(meta_group$group_id), , drop = FALSE]
rownames(meta_group) <- meta_group$group_id
meta_group$n_cells <- as.integer(group_n[meta_group$group_id])
meta_group$Status <- factor(meta_group$Status, levels = c("Ctrl", "AD"))
meta_group$structure <- factor(meta_group$structure)
meta_group$Sex <- factor(meta_group$Sex)

logcpm <- edgeR::cpm(edgeR::DGEList(counts = round(as.matrix(pb))), log = TRUE, prior.count = 2)
rownames(logcpm) <- toupper(rownames(logcpm))
logcpm <- logcpm[!duplicated(rownames(logcpm)), , drop = FALSE]

m1 <- c("HK1", "HK2", "PFKFB3", "PFKFB4", "ALDOA", "ALDOC", "GAPDH", "PGK1", "PGAM1", "ENO1", "ENO2", "PKM", "LDHA", "PDK1", "SLC2A1", "HIF1A", "PFKP")
m3 <- c("CDKN1A", "CDKN2A", "IL6", "IL1B", "TNF", "CCL2", "CCL3", "CCL4", "CSF1", "CXCL10", "MMP9", "MMP14", "SERPINE1", "IGFBP3", "IGFBP7", "NFKB1", "RELA", "TGFB1")
m5 <- c("TFRC", "SLC11A2", "STEAP3", "DMT1", "FTH1", "FTL", "SLC40A1", "CP", "NCOA4", "ATG3", "ATG5", "ATG7", "MAP1LC3B", "SQSTM1")
m6 <- c("ACSL4", "ALOX15", "ALOX5", "ALOX12", "NOX4", "GPX4", "LPCAT3", "SAT1")
m7 <- read.delim(file.path(tables, "D5_ferroaging_reference.tsv"), stringsAsFactors = FALSE)$gene_symbol
m7 <- unique(toupper(m7))
i2 <- read.delim(file.path(tables, "I2_regulatory_activity.tsv"), stringsAsFactors = FALSE)
i2 <- i2[i2$context == "M" & is.finite(i2$delta_presence) & i2$delta_presence != 0 & !is.na(i2$human_one_to_one_symbols) & nzchar(i2$human_one_to_one_symbols), , drop = FALSE]
m2_weights <- tapply(i2$delta_presence, toupper(i2$human_one_to_one_symbols), sum)
modules <- list(M1 = m1, M2 = names(m2_weights), M3 = m3, M5 = m5, M6 = m6, M7 = m7)
module_rows <- list()
row_id <- 0
for (ct in c("Microglia", "Astrocyte")) {
  idx <- which(meta_group$cell_type == ct)
  m <- meta_group[idx, , drop = FALSE]
  design <- model.matrix(~ Status + structure + Age_c + Sex + PMI_c + RIN_c, data = m)
  for (module in names(modules)) {
    requested <- unique(toupper(modules[[module]]))
    present <- intersect(requested, rownames(logcpm))
    if (length(present) < 2) next
    if (module == "M2") {
      w <- as.numeric(m2_weights[present]); names(w) <- present
      score <- colSums(sweep(logcpm[present, m$group_id, drop = FALSE], 1, w, "*")) / sum(abs(w))
    } else {
      score <- colMeans(logcpm[present, m$group_id, drop = FALSE], na.rm = TRUE)
    }
    score_mat <- matrix(as.numeric(score), nrow = 1, dimnames = list(module, m$group_id))
    corfit <- tryCatch(limma::duplicateCorrelation(score_mat, design, block = m$id), error = function(e) NULL)
    if (!is.null(corfit) && is.finite(corfit$consensus)) {
      fit <- limma::lmFit(score_mat, design, block = m$id, correlation = corfit$consensus)
    } else {
      fit <- limma::lmFit(score_mat, design)
    }
    fit <- limma::eBayes(fit)
    row_id <- row_id + 1
    module_rows[[row_id]] <- data.frame(module = module, cell_type = ct, n_total = length(requested), n_present = length(present), coverage = length(present)/length(requested), n_groups = nrow(m), n_AD_groups = sum(m$Status == "AD"), n_Control_groups = sum(m$Status == "Ctrl"), logFC = unname(fit$coefficients[1, "StatusAD"]), t = unname(fit$t[1, "StatusAD"]), P.Value = unname(fit$p.value[1, "StatusAD"]), stringsAsFactors = FALSE)
  }
}
module_effects <- do.call(rbind, module_rows)
module_effects$FDR <- p.adjust(module_effects$P.Value, method = "BH")
module_effects$direction <- ifelse(module_effects$logFC > 0, "AD_up", "AD_down")
write.table(module_effects, file.path(tables, "D6_module_effects.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# Independent direction/concordance check against the frozen D3 ranked tables.
counts <- round(as.matrix(pb))
dge <- edgeR::DGEList(counts = counts)
dge <- edgeR::calcNormFactors(dge)
vo <- limma::voom(dge, design = model.matrix(~ Status + structure + Age_c + Sex + PMI_c + RIN_c, data = meta_group), plot = FALSE)
concordance <- list()
for (ct in c("Microglia", "Astrocyte")) {
  idx <- which(meta_group$cell_type == ct)
  m <- meta_group[idx, , drop = FALSE]
  design <- model.matrix(~ Status + structure + Age_c + Sex + PMI_c + RIN_c, data = m)
  v <- vo$E[, m$group_id, drop = FALSE]
  corfit <- tryCatch(limma::duplicateCorrelation(v, design, block = m$id), error = function(e) NULL)
  if (!is.null(corfit) && is.finite(corfit$consensus)) fit <- limma::lmFit(v, design, block = m$id, correlation = corfit$consensus) else fit <- limma::lmFit(v, design)
  fit <- limma::eBayes(fit)
  d3 <- read.delim(file.path(tables, paste0("DiseaseSignature_", tolower(ct), ".tsv")), stringsAsFactors = FALSE)
  d3$gene_upper <- toupper(d3$gene)
  d6 <- data.frame(gene_upper = toupper(rownames(fit$coefficients)), d6_logFC = fit$coefficients[, "StatusAD"], stringsAsFactors = FALSE)
  merged <- merge(d3[, c("gene_upper", "logFC")], d6, by = "gene_upper")
  concordance[[ct]] <- data.frame(cell_type = ct, n_shared_genes = nrow(merged), spearman_rho = cor(merged$logFC, merged$d6_logFC, method = "spearman"), direction_concordance = mean(sign(merged$logFC) == sign(merged$d6_logFC)), stringsAsFactors = FALSE)
}
concordance <- do.call(rbind, concordance)
write.table(concordance, file.path(tables, "D6_signature_concordance.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
prov <- list(dataset = "D6_GSE214979", h5_expected_bytes = expected_bytes, h5_actual_bytes = actual_bytes, matched_cells = length(keep_cells), analyzed_cells = ncol(x), min_cells_per_donor_region_celltype = 20, primary_unit = "donor x structure x glial cell type pseudobulk", covariates = "Status + structure + Age + Sex + PMI + RIN", outputs = c("D6_module_effects.tsv", "D6_signature_concordance.tsv"), boundaries = c("independent validation, not cross-dataset causal modeling", "module/RNA signatures do not measure H4K12la, lactate, iron or lipid flux directly", "duplicateCorrelation blocks are all size 1 after donor x structure x cell type aggregation; intrablock correlation is therefore set to zero"))
jsonlite::write_json(prov, file.path(interim, "D6_expression_validation_provenance.json"), auto_unbox = TRUE, pretty = TRUE)
cat(sprintf("D6 validation completed: %d matched cells, %d pseudobulk groups\n", length(keep_cells), ncol(pb)))