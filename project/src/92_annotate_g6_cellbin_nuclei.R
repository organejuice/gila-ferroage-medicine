#!/usr/bin/env Rscript

# Exploratory cell identity reconstruction for G6-FP-2026-09-27-001, WP2.
# This is not an exact rerun of the authors' Scanpy 1.9.3 workflow: the pinned
# author repository/release could not be retrieved or verified in this session.
suppressPackageStartupMessages({
  library(Matrix); library(data.table); library(Seurat); library(harmony)
})
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: 92_annotate_g6_cellbin_nuclei.R <project_root>")
root <- args[[1L]]
input_dir <- file.path(root, "data/interim/G6_next_feasibility_001/cellbin")
audit_dir <- file.path(root, "results/2026-09-24-disease-side-v3.5-g6-input-audit")
out_dir <- file.path(root, "results/2026-09-27-g6-next-feasibility/cellbin")
if (length(list.files(out_dir, pattern = "annotated_nuclei_metadata\\.tsv$"))) stop("Existing annotation result; refusing overwrite")
gsm <- c("GSM8330060", "GSM8330061", "GSM8330062", "GSM8330063", "GSM8330064", "GSM8330065",
         "GSM8330066", "GSM8330067", "GSM8330068", "GSM8330069", "GSM8330070", "GSM8330071")
paths <- file.path(input_dir, paste0(gsm, "_nucleus_counts.rds"))
if (!all(file.exists(paths))) stop("All 12 direct CellBin nucleus matrices are required")
g1 <- fread(file.path(audit_dir, "registries/G1_microglia_signed_sender_244.tsv.gz"))
g5 <- fread(file.path(audit_dir, "registries/G5_receiver_one_to_one_human_ortholog_registry.tsv.gz"))
if (nrow(g1) != 244L || uniqueN(g1$gene) != 244L) stop("G1 registry integrity check failed")
state_genes <- unique(toupper(c(g1$gene, g5$human_one_to_one_symbol)))
state_genes <- state_genes[!is.na(state_genes) & nzchar(state_genes)]
identity_markers_raw <- list(
  MG = c("CD74", "CD14"), ASC = c("AQP4", "GLUL"),
  End = c("CLDN5", "IFITM2"), Ex = c("NRGN", "SYT1"), Inh = c("GAD1", "GAD2"),
  Oli = c("MBP", "PLP1"), OPC = "VCAN"
)
# Training markers may also be members of the frozen state registry. They are
# retained for annotation but excluded from variable-gene clustering below.
identity_markers <- identity_markers_raw
mg_holdout <- setdiff(c("AIF1", "C1QA", "C1QB", "C1QC", "CX3CR1", "LST1", "P2RY12", "TMEM119", "TREM2", "TYROBP"), state_genes)
asc_holdout <- setdiff(c("ALDH1L1", "GFAP", "GJA1", "S100B", "SLC1A2", "SLC1A3", "SOX9"), state_genes)
if (length(mg_holdout) < 3L || length(asc_holdout) < 3L) stop("Insufficient disjoint holdout markers")
for (nm in names(identity_markers)) identity_markers[[nm]] <- unique(identity_markers[[nm]])
excluded_features <- unique(c(state_genes, unlist(identity_markers, use.names = FALSE), mg_holdout, asc_holdout))

matrix_parts <- vector("list", length(paths))
metadata_parts <- vector("list", length(paths))
qc_rows <- vector("list", length(paths))
for (i in seq_along(paths)) {
  x <- readRDS(paths[[i]])
  if (!identical(x$study_id, "G6-FP-2026-09-27-001") ||
      any(grepl("diagnos|patholog|braak|cerad|thal|control", tolower(names(x$cell_metadata))))) {
    stop("Invalid or diagnosis-contaminated CellBin object: ", gsm[[i]])
  }
  meta <- x$cell_metadata
  if (!identical(colnames(x$counts), rownames(meta))) stop("Count/cell metadata mismatch: ", gsm[[i]])
  keep <- meta$n_genes >= 150L & meta$polygon_valid
  qc_rows[[i]] <- data.frame(GSM = gsm[[i]], donor_id = unique(meta$donor_id),
    n_cells_pre_filter = nrow(meta), n_cells_ge150_and_polygon_valid = sum(keep),
    median_n_genes_pre = median(meta$n_genes), median_n_genes_post = if (any(keep)) median(meta$n_genes[keep]) else NA_real_,
    median_n_umi_post = if (any(keep)) median(meta$n_umi[keep]) else NA_real_,
    stringsAsFactors = FALSE)
  matrix_parts[[i]] <- x$counts[, keep, drop = FALSE]
  metadata_parts[[i]] <- meta[keep, , drop = FALSE]
  rm(x)
  gc(verbose = FALSE)
}
qc <- rbindlist(qc_rows, fill = TRUE)
if (nrow(qc) != 12L || uniqueN(qc$donor_id) != 12L) stop("CellBin sample crosswalk is not 12/12")
if (sum(qc$n_cells_ge150_and_polygon_valid) < 100L) stop("Too few nuclei survive the predeclared <150 gene filter to cluster")
all_meta <- rbindlist(metadata_parts, fill = TRUE)
gene_universe <- sort(unique(unlist(lapply(matrix_parts, rownames), use.names = FALSE)))
if (!length(gene_universe) || anyNA(gene_universe) || any(!nzchar(gene_universe))) stop("CellBin gene symbols are missing or invalid")
align_gene_universe <- function(m) {
  if (is.null(rownames(m)) || anyDuplicated(rownames(m))) stop("A donor matrix lacks unique gene symbols")
  source_i <- m@i + 1L
  source_j <- rep.int(seq_len(ncol(m)), diff(m@p))
  sparseMatrix(i = match(rownames(m)[source_i], gene_universe), j = source_j, x = m@x,
    dims = c(length(gene_universe), ncol(m)), dimnames = list(gene_universe, colnames(m)), giveCsparse = TRUE)
}
aligned_parts <- lapply(matrix_parts, align_gene_universe)
all_counts <- do.call(cbind, aligned_parts)
rm(aligned_parts, gene_universe)
expected_cells <- paste(all_meta$GSM, sprintf("cell%08d", all_meta$cell_id), sep = "::")
if (!identical(colnames(all_counts), expected_cells) || anyDuplicated(colnames(all_counts))) stop("Merged count matrix/cell map mismatch")
if (any(Matrix::colSums(all_counts > 0) != all_meta$n_genes) || any(Matrix::colSums(all_counts) != all_meta$n_umi)) stop("Merged counts fail per-cell QC")
rm(matrix_parts, metadata_parts)
gc(verbose = FALSE)

# Diagnosis-blind exploratory analogue of the paper's pooled, donor-harmonized clustering.
set.seed(42L)
meta_df <- as.data.frame(all_meta)
rownames(meta_df) <- expected_cells
obj <- CreateSeuratObject(counts = all_counts, meta.data = meta_df,
                          project = "G6_CellBin_blind", min.cells = 0, min.features = 0)
rm(all_counts, all_meta, meta_df, expected_cells)
obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
allowed_features <- rownames(obj)[!toupper(rownames(obj)) %in% excluded_features]
if (length(allowed_features) < 2000L) stop("Fewer than 2,000 non-signature clustering genes remain")
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000L,
                            features = allowed_features, verbose = FALSE)
cluster_features <- intersect(VariableFeatures(obj), allowed_features)
if (length(cluster_features) < 500L) stop("Fewer than 500 permitted variable genes remain")
npcs <- min(50L, length(cluster_features) - 1L, ncol(obj) - 1L)
# Seurat ScaleData materializes a dense 2,000 x ~400,000 matrix (>6 GiB).
# Use irlba's sparse centering/scaling operator, which computes the same
# feature-centered, feature-scaled PCA without constructing that dense matrix.
if (!requireNamespace("irlba", quietly = TRUE)) stop("The locked R environment lacks irlba")
lognorm <- GetAssayData(obj, assay = "RNA", layer = "data")[cluster_features, , drop = FALSE]
pca_input <- Matrix::t(lognorm)
feature_center <- Matrix::colMeans(pca_input)
feature_second_moment <- Matrix::colMeans(pca_input ^ 2)
feature_scale <- sqrt(pmax(feature_second_moment - feature_center ^ 2, 0))
feature_scale[!is.finite(feature_scale) | feature_scale == 0] <- 1
set.seed(42L)
pca <- irlba::irlba(pca_input, nv = npcs, nu = npcs, center = feature_center,
                    scale = feature_scale, fastpath = TRUE, tol = 1e-5)
embeddings <- pca$u %*% diag(pca$d)
rownames(embeddings) <- rownames(pca_input)
colnames(embeddings) <- paste0("PC_", seq_len(npcs))
loadings <- pca$v
rownames(loadings) <- colnames(pca_input)
colnames(loadings) <- colnames(embeddings)
obj[["pca"]] <- CreateDimReducObject(embeddings = embeddings, loadings = loadings,
                                     key = "PC_", assay = "RNA")
rm(lognorm, pca_input, pca, embeddings, loadings, feature_center, feature_second_moment, feature_scale)
gc(verbose = FALSE)
obj <- RunHarmony(obj, group.by.vars = "GSM", reduction.use = "pca", dims.use = seq_len(npcs),
                  reduction.save = "harmony", project.dim = FALSE, verbose = FALSE)
obj <- FindNeighbors(obj, reduction = "harmony", dims = seq_len(npcs), k.param = 15L, verbose = FALSE)
obj <- FindClusters(obj, algorithm = 4L, resolution = 1.0, random.seed = 42L, verbose = FALSE)
cluster_id <- as.character(Idents(obj))
names(cluster_id) <- colnames(obj)
cluster_ids <- sort(unique(cluster_id))
norm <- GetAssayData(obj, assay = "RNA", layer = "data")
raw <- GetAssayData(obj, assay = "RNA", layer = "counts")
gene_upper <- toupper(rownames(norm))
measured_identity <- lapply(identity_markers, function(g) {
  idx <- match(g, gene_upper)
  unique(rownames(norm)[idx[!is.na(idx)]])
})
if (length(measured_identity$MG) < 2L || length(measured_identity$ASC) < 2L) stop("Training marker pairs are not measurable")

# Cluster label only when the required training marker(s) are co-detected and the
# corresponding marker-module score is the unique maximum; all other clusters stay unknown.
cluster_label_rows <- vector("list", length(cluster_ids))
cluster_scores <- matrix(NA_real_, nrow = length(cluster_ids), ncol = length(identity_markers),
                         dimnames = list(cluster_ids, names(identity_markers)))
for (k in seq_along(cluster_ids)) {
  id <- cluster_ids[[k]]
  cells <- names(cluster_id)[cluster_id == id]
  ix <- match(cells, colnames(norm))
  for (lineage in names(identity_markers)) {
    genes <- measured_identity[[lineage]]
    if (!length(genes)) next
    cluster_scores[k, lineage] <- mean(Matrix::rowMeans(norm[genes, ix, drop = FALSE]))
    pct <- Matrix::rowMeans(raw[genes, ix, drop = FALSE] > 0)
  }
  required <- list(MG = c("CD74", "CD14"), ASC = c("AQP4", "GLUL"), End = c("CLDN5", "IFITM2"),
                   Ex = c("NRGN", "SYT1"), Inh = c("GAD1", "GAD2"), Oli = c("MBP", "PLP1"), OPC = "VCAN")
  supported <- names(required)[vapply(names(required), function(lineage) {
    gs <- intersect(required[[lineage]], measured_identity[[lineage]])
    if (!length(gs)) return(FALSE)
    all(Matrix::rowMeans(raw[gs, ix, drop = FALSE] > 0) > 0)
  }, logical(1))]
  ord <- order(cluster_scores[k, ], decreasing = TRUE, na.last = NA)
  best <- if (length(ord)) names(identity_markers)[ord[[1L]]] else NA_character_
  second <- if (length(ord) > 1L) cluster_scores[k, ord[[2L]]] else NA_real_
  label <- if (length(supported) && best %in% supported && is.finite(second) && cluster_scores[k, best] > second) best else "unknown"
  cluster_label_rows[[k]] <- data.frame(cluster = id, n_cells = length(cells), label = label,
    top_marker_lineage = best, top_score = if (length(ord)) cluster_scores[k, best] else NA_real_,
    second_score = second, training_marker_supported = label != "unknown", stringsAsFactors = FALSE)
}
cluster_table <- rbindlist(cluster_label_rows, fill = TRUE)
label_map <- setNames(cluster_table$label, cluster_table$cluster)
cell_label <- unname(label_map[cluster_id])
names(cell_label) <- names(cluster_id)

score_holdout <- function(genes, cells) {
  rows <- match(toupper(genes), gene_upper)
  rows <- unique(rows[!is.na(rows)])
  if (length(rows) < 3L) return(rep(NA_real_, length(cells)))
  idx <- match(cells, colnames(raw))
  umi <- as.numeric(Matrix::colSums(raw[, idx, drop = FALSE]))
  cpm <- t(t(raw[rows, idx, drop = FALSE]) / pmax(umi, 1)) * 1e6
  as.numeric(Matrix::colMeans(log2(as.matrix(cpm) + 1)))
}
cell_tables <- vector("list", 12L)
marker_rows <- list(); g1_rows <- list(); g5_rows <- list()
for (i in seq_len(nrow(qc))) {
  donor <- qc$donor_id[[i]]; sample <- qc$GSM[[i]]
  cells <- rownames(obj[[]])[obj[[]]$donor_id == donor]
  labs <- cell_label[cells]
  mg_score <- score_holdout(mg_holdout, cells)
  asc_score <- score_holdout(asc_holdout, cells)
  for (lineage in c("MG", "ASC")) {
    target <- labs == lineage
    scores <- if (lineage == "MG") mg_score else asc_score
    target_median <- if (any(target)) median(scores[target], na.rm = TRUE) else NA_real_
    other_median <- if (any(!target)) median(scores[!target], na.rm = TRUE) else NA_real_
    holdout_count <- length(unique(na.omit(match(toupper(if (lineage == "MG") mg_holdout else asc_holdout), gene_upper))))
    marker_rows[[paste(sample, lineage)]] <- data.frame(GSM = sample, donor_id = donor, lineage = lineage,
      n_target_nuclei = sum(target), n_other_nuclei = sum(!target), n_holdout_markers = holdout_count,
      median_holdout_score_target = target_median,
      median_holdout_score_other = other_median,
      target_marker_score_higher = is.finite(target_median) && is.finite(other_median) && target_median > other_median,
      stringsAsFactors = FALSE)
  }
  mg_cells <- cells[labs == "MG"]
  g1_symbols <- toupper(g1$gene)
  g1_detect <- setNames(rep(0, length(g1_symbols)), g1_symbols)
  g1_rows_in_matrix <- match(g1_symbols, gene_upper)
  if (length(mg_cells)) {
    use <- which(!is.na(g1_rows_in_matrix))
    if (length(use)) g1_detect[use] <- Matrix::rowMeans(raw[g1_rows_in_matrix[use], match(mg_cells, colnames(raw)), drop = FALSE] > 0)
  }
  g1_rows[[sample]] <- data.frame(GSM = sample, donor_id = donor, gene = names(g1_detect),
    detection_fraction_in_MG = unname(g1_detect), meets_1pct = unname(g1_detect) >= 0.01, stringsAsFactors = FALSE)
  asc_cells <- cells[labs == "ASC"]
  for (module in sort(unique(g5$module))) {
    genes <- unique(toupper(g5$human_one_to_one_symbol[g5$module == module]))
    genes <- genes[!is.na(genes) & nzchar(genes)]
    gene_rows <- match(genes, gene_upper)
    present <- which(!is.na(gene_rows))
    fractions <- rep(0, length(genes))
    if (length(asc_cells) && length(present)) fractions[present] <- Matrix::rowMeans(raw[gene_rows[present], match(asc_cells, colnames(raw)), drop = FALSE] > 0)
    n_detected <- sum(fractions >= 0.01)
    required_n <- ceiling(length(genes) * 0.25)
    g5_rows[[paste(sample, module)]] <- data.frame(GSM = sample, donor_id = donor, module = module,
      n_module_genes = length(genes), n_detected_at_1pct = n_detected, required_at_25pct = required_n,
      module_pass = length(asc_cells) > 0L && n_detected >= required_n, stringsAsFactors = FALSE)
  }
  cell_tables[[i]] <- data.frame(cell_id = cells, GSM = sample, donor_id = donor, x_um = obj[[]][cells, "x_um"],
    y_um = obj[[]][cells, "y_um"], n_genes = obj[[]][cells, "n_genes"], n_umi = obj[[]][cells, "n_umi"],
    cluster = cluster_id[cells], label = labs, polygon_valid = obj[[]][cells, "polygon_valid"], stringsAsFactors = FALSE)
  message("ANNOTATED ", sample, "; filtered_nuclei=", length(cells), "; MG=", sum(labs == "MG"), "; ASC=", sum(labs == "ASC"))
  gc(verbose = FALSE)
}
marker_table <- rbindlist(marker_rows, fill = TRUE)
g1_table <- rbindlist(g1_rows, fill = TRUE)
g5_table <- rbindlist(g5_rows, fill = TRUE)
cell_table <- rbindlist(cell_tables, fill = TRUE)
marker_gate <- marker_table[, .(donors_with_direction = sum(target_marker_score_higher), donors_total = .N), by = lineage]
g1_donor_gate <- g1_table[, .(n_G1_genes_at_1pct = sum(meets_1pct)), by = .(GSM, donor_id)]
g5_module_gate <- g5_table[, .(donors_meeting_module_coverage = sum(module_pass), donors_total = .N), by = module]

write.table(qc, file.path(out_dir, "cellbin_nucleus_qc_after_150_genes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")
write.table(cluster_table, file.path(out_dir, "exploratory_cluster_marker_annotation.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")
write.table(marker_table, file.path(out_dir, "independent_marker_validation_by_donor.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")
write.table(marker_gate, file.path(out_dir, "identity_marker_gate.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")
write.table(g1_table, file.path(out_dir, "G1_detection_by_donor.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")
write.table(g1_donor_gate, file.path(out_dir, "G1_coverage_gate_by_donor.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")
write.table(g5_table, file.path(out_dir, "G5_module_detection_by_donor.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")
write.table(g5_module_gate, file.path(out_dir, "G5_coverage_gate_by_module.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")
write.table(cell_table, file.path(out_dir, "annotated_nuclei_metadata.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")
write.table(data.frame(feature = cluster_features), file.path(out_dir, "identity_clustering_feature_list.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "Seurat_exploratory_sessionInfo.txt"))

identity_pass <- all(marker_gate$donors_with_direction >= 10L & marker_gate$donors_total == 12L)
g1_pass <- sum(g1_donor_gate$n_G1_genes_at_1pct >= 61L) >= 10L
g5_pass <- all(g5_module_gate$donors_meeting_module_coverage >= 10L & g5_module_gate$donors_total == 12L)
status <- if (!identity_pass) "IDENTITY_NOT_CALIBRATED" else if (!g1_pass || !g5_pass) "CELLBIN_NOT_ELIGIBLE_COVERAGE" else "IDENTITY_AND_COVERAGE_ELIGIBLE_FOR_NEIGHBOR_CHECK"
writeLines(c(
  "# G6 CellBin exploratory annotation result", "",
  paste("Status:", status),
  "Author Scanpy code/version could not be verified; this is a separately marked exploratory Seurat/Harmony reconstruction, not an exact author-pipeline replication.",
  "No diagnosis or pathology columns were loaded. G1/G5 and all training/holdout identity markers were excluded from clustering features; training markers labeled clusters and disjoint holdout markers validated labels.",
  paste("Filtered nuclei:", nrow(cell_table), "from 12 donors; 150-gene threshold and valid polygons applied."),
  paste("Identity gate:", paste(capture.output(print(marker_gate)), collapse = " | ")),
  paste("G1 donor gate:", paste(capture.output(print(g1_donor_gate)), collapse = " | ")),
  paste("G5 module gate:", paste(capture.output(print(g5_module_gate)), collapse = " | ")),
  "This analysis is diagnosis-blind and cannot produce a G6_PASS. Neighborhood eligibility is assessed separately and requires an auditable tissue-hole mask."
), file.path(out_dir, "exploratory_annotation_report.md"))
message("WP2 annotation status: ", status, "; identity=", identity_pass, "; G1=", g1_pass, "; G5=", g5_pass)
