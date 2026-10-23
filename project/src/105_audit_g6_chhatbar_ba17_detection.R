#!/usr/bin/env Rscript
# Per-sample BA17 measurement inventory; no between-diagnosis effect computed.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop("Usage: Rscript src/105_audit_g6_chhatbar_ba17_detection.R VISIUM_RDS G1_TSV_GZ G5_TSV_GZ OUT_TSV")
}
object <- readRDS(args[[1]])
if (!inherits(object, "Seurat")) stop("Expected Seurat object")
metadata <- object@meta.data
selected <- metadata$Region == "Occipital_cortex" &
  metadata$Clinical_Brain_Diagnosis_Short %in% c("AD", "Control")
metadata <- metadata[selected, , drop = FALSE]
if (length(unique(metadata$Subject_ID)) != 10L) stop("Unexpected BA17 sample count")
counts <- object@assays$RNA@counts
if (nrow(counts) != 18085L || ncol(counts) != 96409L) stop("Unexpected count matrix")
g1 <- unique(read.delim(gzfile(args[[2]]), check.names = FALSE)$gene)
g5 <- read.delim(gzfile(args[[3]]), check.names = FALSE)
accepted <- startsWith(g5$one_to_one_mapping_status, "accepted") &
  !is.na(g5$human_one_to_one_symbol) & g5$human_one_to_one_symbol != ""
modules <- c(list(G1_Microglia = g1),
             lapply(split(g5$human_one_to_one_symbol[accepted], g5$module[accepted]), unique))
results <- list()
for (sample in sort(unique(as.character(metadata$Subject_ID)))) {
  for (module in names(modules)) {
    target_type <- if (module == "G1_Microglia") "Microglia" else "Astrocytes"
    ids <- rownames(metadata)[as.character(metadata$Subject_ID) == sample &
                             as.character(metadata$Final_celltype) == target_type]
    genes <- unique(modules[[module]])
    present <- intersect(genes, rownames(counts))
    detected <- if (length(ids) && length(present)) {
      sum(Matrix::rowSums(counts[present, ids, drop = FALSE] > 0) > 0)
    } else 0L
    results[[length(results) + 1L]] <- data.frame(
      sample_id = sample,
      diagnosis = as.character(metadata$Clinical_Brain_Diagnosis_Short[match(sample, metadata$Subject_ID)]),
      region = "Occipital_cortex_BA17",
      module = module,
      target_celltype = target_type,
      target_cells = length(ids),
      frozen_genes = length(genes),
      object_features = length(present),
      genes_detected_at_least_one_target_cell = detected,
      stringsAsFactors = FALSE
    )
  }
}
dir.create(dirname(args[[4]]), recursive = TRUE, showWarnings = FALSE)
write.table(do.call(rbind, results), args[[4]], sep = "\t", row.names = FALSE, quote = FALSE)
