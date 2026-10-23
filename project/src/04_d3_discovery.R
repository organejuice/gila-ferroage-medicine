#!/usr/bin/env Rscript
# D3 GSE160936: sample-aware QC, glial annotation, pseudobulk, and signed DE.
# The workflow follows the local single-nucleus QC and Seurat/Harmony module,
# while deriving QC cutoffs from sample distributions rather than importing
# project-specific hard thresholds.

suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(Matrix)
  library(edgeR)
  library(limma)
  library(yaml)
})

parse_root <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  hit <- grep("^--root=", args, value = TRUE)
  if (length(hit)) sub("^--root=", "", hit[[1]]) else getwd()
}

root <- gsub("\\\\", "/", parse_root())
config <- yaml::read_yaml(file.path(root, "configs", "d3_params.yaml"))$d3
raw_dir <- file.path(root, "data", "raw", "D3_GSE160936")
interim_dir <- file.path(root, "data", "interim")
processed_dir <- file.path(root, "data", "processed")
tables_dir <- file.path(root, "results", "tables")
figures_dir <- file.path(root, "results", "figures")
for (d in c(interim_dir, processed_dir, tables_dir, figures_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

metadata <- read.delim(file.path(interim_dir, "D3_sample_metadata.tsv"), stringsAsFactors = FALSE, check.names = FALSE)
stopifnot(nrow(metadata) == 24, all(metadata$archive_exists))

extract_one <- function(row) {
  sample_dir <- file.path(interim_dir, "D3_matrices", row$geo_sample)
  matrix_dir <- file.path(sample_dir, "filtered_feature_bc_matrix")
  expected <- file.path(matrix_dir, c("features.tsv.gz", "barcodes.tsv.gz", "matrix.mtx.gz"))
  if (!all(file.exists(expected))) {
    dir.create(sample_dir, recursive = TRUE, showWarnings = FALSE)
    utils::untar(file.path(root, row$archive), exdir = sample_dir)
    found <- list.files(sample_dir, pattern = "(features|barcodes|matrix)\\.tsv\\.gz$", recursive = TRUE, full.names = TRUE)
    if (length(found) != 3) stop("Could not locate all MEX members for ", row$geo_sample)
    dir.create(matrix_dir, recursive = TRUE, showWarnings = FALSE)
    for (f in found) file.copy(f, file.path(matrix_dir, basename(f)), overwrite = TRUE)
  }
  if (!all(file.exists(expected))) stop("MEX extraction incomplete for ", row$geo_sample)
  matrix_dir
}

message("Extracting/verifying 24 Cell Ranger matrices")
matrix_dirs <- vapply(seq_len(nrow(metadata)), function(i) extract_one(metadata[i, ]), character(1))

feature_reference <- NULL
read_one <- function(i) {
  gsm <- metadata$geo_sample[[i]]
  mat <- Seurat::Read10X(data.dir = matrix_dirs[[i]], gene.column = 2)
  if (is.list(mat)) {
    if (!"Gene Expression" %in% names(mat)) stop("No Gene Expression matrix in ", gsm)
    mat <- mat[["Gene Expression"]]
  }
  if (is.null(feature_reference)) feature_reference <<- rownames(mat)
  if (!identical(rownames(mat), feature_reference)) stop("Feature order differs in ", gsm)
  colnames(mat) <- paste(gsm, colnames(mat), sep = "_")
  row <- metadata[i, ]
  cell_meta <- data.frame(
    sample_id = rep(gsm, ncol(mat)),
    donor_id = rep(row$donor_id, ncol(mat)),
    disease = rep(row$disease, ncol(mat)),
    region = rep(row$region, ncol(mat)),
    age = rep(as.numeric(row$age), ncol(mat)),
    sex = rep(row$sex, ncol(mat)),
    rin = rep(as.numeric(row$rin), ncol(mat)),
    braak_tangle_stage = rep(as.numeric(row$braak_tangle_stage), ncol(mat)),
    ptau_percent_positive_cells = rep(as.numeric(row$ptau_percent_positive_cells), ncol(mat)),
    amyloid_percent_area = rep(as.numeric(row$amyloid_percent_area), ncol(mat)),
    sample_title = rep(row$sample_title, ncol(mat)),
    row.names = colnames(mat),
    stringsAsFactors = FALSE
  )
  list(counts = mat, meta = cell_meta)
}

message("Reading sparse matrices")
# Read and cbind once; repeated Seurat merge copies the 59k-feature matrix many times.
loaded <- lapply(seq_len(nrow(metadata)), read_one)
counts_list <- lapply(loaded, `[[`, "counts")
combined_counts <- do.call(cbind, counts_list)
cell_meta <- do.call(rbind, lapply(loaded, `[[`, "meta"))
stopifnot(identical(colnames(combined_counts), rownames(cell_meta)))
combined <- Seurat::CreateSeuratObject(counts = combined_counts, meta.data = cell_meta, min.cells = 0, min.features = 0, project = "GSE160936")
DefaultAssay(combined) <- "RNA"
combined[["percent.mt"]] <- Seurat::PercentageFeatureSet(combined, pattern = "^MT-")
combined[["percent.ribo"]] <- Seurat::PercentageFeatureSet(combined, pattern = "^(RPS|RPL)")

md0 <- combined[[]]
qc_summary <- do.call(rbind, lapply(split(md0, md0$sample_id), function(df) {
  data.frame(
    sample_id = df$sample_id[[1]],
    n_cells_before = nrow(df),
    median_features = median(df$nFeature_RNA),
    q01_features = as.numeric(quantile(df$nFeature_RNA, 0.01)),
    q99_features = as.numeric(quantile(df$nFeature_RNA, 0.99)),
    median_umis = median(df$nCount_RNA),
    median_mt_percent = median(df$percent.mt),
    q99_mt_percent = as.numeric(quantile(df$percent.mt, 0.99)),
    stringsAsFactors = FALSE
  )
}))
qc_thresholds <- do.call(rbind, lapply(split(md0, md0$sample_id), function(df) {
  data.frame(
    sample_id = df$sample_id[[1]],
    min_features = as.numeric(config$min_features),
    max_features = as.numeric(quantile(df$nFeature_RNA, config$max_feature_quantile, names = FALSE)),
    max_mt_percent = as.numeric(quantile(df$percent.mt, config$max_mt_quantile, names = FALSE)),
    stringsAsFactors = FALSE
  )
}))
write.table(qc_summary, file.path(tables_dir, "D3_qc_summary_before_filter.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(qc_thresholds, file.path(tables_dir, "D3_qc_thresholds.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

max_features <- setNames(qc_thresholds$max_features, qc_thresholds$sample_id)
max_mt <- setNames(qc_thresholds$max_mt_percent, qc_thresholds$sample_id)
keep <- md0$nFeature_RNA >= as.numeric(config$min_features) &
  md0$nFeature_RNA <= max_features[md0$sample_id] &
  md0$percent.mt <= max_mt[md0$sample_id]
qc_audit <- data.frame(
  sample_id = names(table(md0$sample_id)),
  n_cells_before = as.integer(table(md0$sample_id)),
  n_cells_after = as.integer(table(md0$sample_id)[names(table(md0$sample_id))]),
  stringsAsFactors = FALSE
)
after_counts <- table(md0$sample_id[keep])
qc_audit$n_cells_after <- as.integer(after_counts[qc_audit$sample_id])
qc_audit$n_cells_after[is.na(qc_audit$n_cells_after)] <- 0
qc_audit$retention_fraction <- qc_audit$n_cells_after / qc_audit$n_cells_before
qc_audit$doublet_filter <- "not_available_in_public_filtered_matrix"
write.table(qc_audit, file.path(tables_dir, "D3_qc_audit.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

filtered <- subset(combined, cells = rownames(md0)[keep])
filtered <- Seurat::NormalizeData(filtered, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
filtered <- Seurat::FindVariableFeatures(filtered, selection.method = "vst", nfeatures = as.numeric(config$nfeatures), verbose = FALSE)
filtered <- Seurat::ScaleData(filtered, features = Seurat::VariableFeatures(filtered), vars.to.regress = "percent.mt", verbose = FALSE)
filtered <- Seurat::RunPCA(filtered, features = Seurat::VariableFeatures(filtered), npcs = as.numeric(config$npcs), verbose = FALSE)
filtered <- harmony::RunHarmony(filtered, group.by.vars = "sample_id", reduction.use = "pca", dims.use = seq_len(as.numeric(config$npcs)), verbose = FALSE)
filtered <- Seurat::FindNeighbors(filtered, reduction = "harmony", dims = seq_len(as.numeric(config$npcs)), k.param = as.numeric(config$neighbor_k), verbose = FALSE)
filtered <- Seurat::FindClusters(filtered, resolution = as.numeric(config$cluster_resolution), algorithm = 1, verbose = FALSE)
filtered <- Seurat::RunUMAP(filtered, reduction = "harmony", dims = seq_len(as.numeric(config$npcs)), verbose = FALSE)

marker_sets <- list(
  Microglia = c("AIF1", "CX3CR1", "P2RY12", "TMEM119", "C1QA", "C1QB", "C1QC", "TYROBP", "LST1", "SPI1"),
  Astrocyte = c("AQP4", "GFAP", "SLC1A2", "SLC1A3", "ALDH1L1", "SOX9", "GJA1", "CLU"),
  Neuron = c("SNAP25", "RBFOX3", "SYT1", "STMN2", "GRIN1"),
  Oligodendrocyte = c("MBP", "PLP1", "MOG", "MOBP", "OLIG1", "OLIG2"),
  OPC = c("PDGFRA", "CSPG4", "VCAN"),
  Endothelial = c("CLDN5", "KDR", "EMCN", "FLT1"),
  Pericyte = c("RGS5", "PDGFRB", "CSPG4", "MCAM"),
  Inhibitory = c("GAD1", "GAD2", "SLC6A1"),
  Excitatory = c("SLC17A7", "CAMK2A", "TBR1")
)
marker_sets <- lapply(marker_sets, function(x) intersect(x, rownames(filtered)))
marker_sets <- marker_sets[lengths(marker_sets) >= 3]
for (nm in names(marker_sets)) {
  filtered <- Seurat::AddModuleScore(filtered, features = list(marker_sets[[nm]]), name = paste0(nm, "Score"), assay = "RNA", search = FALSE)
}
score_cols <- paste0(names(marker_sets), "Score1")
md <- filtered[[]]
cluster_ids <- as.character(Idents(filtered))
cluster_levels <- sort(unique(cluster_ids))
cluster_scores <- do.call(rbind, lapply(cluster_levels, function(cl) {
  idx <- which(cluster_ids == cl)
  vals <- vapply(score_cols, function(col) mean(md[idx, col], na.rm = TRUE), numeric(1))
  data.frame(cluster = cl, t(vals), stringsAsFactors = FALSE, check.names = FALSE)
}))
label_from_scores <- function(x) {
  vals <- as.numeric(x[score_cols])
  names(vals) <- names(x[score_cols])
  ord <- order(vals, decreasing = TRUE)
  best <- sub("Score1$", "", names(vals)[ord[[1]]])
  margin <- vals[ord[[1]]] - vals[ord[[2]]]
  if (!is.finite(margin) || margin < 0.05) "Unresolved" else best
}
cluster_scores$cell_type <- apply(cluster_scores, 1, label_from_scores)
cluster_to_type <- setNames(cluster_scores$cell_type, cluster_scores$cluster)
filtered$cell_type <- unname(cluster_to_type[cluster_ids])
write.table(cluster_scores, file.path(tables_dir, "D3_cluster_marker_scores.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

cell_summary <- as.data.frame(table(filtered$sample_id, filtered$disease, filtered$region, filtered$cell_type), stringsAsFactors = FALSE)
colnames(cell_summary) <- c("sample_id", "disease", "region", "cell_type", "n_cells")
write.table(cell_summary, file.path(tables_dir, "D3_celltype_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

if (requireNamespace("ggplot2", quietly = TRUE)) {
  p <- Seurat::DimPlot(filtered, reduction = "umap", group.by = "cell_type", raster = TRUE) + ggplot2::ggtitle("GSE160936 D3 cell types")
  ggplot2::ggsave(file.path(figures_dir, "D3_umap_celltype.png"), p, width = 10, height = 7, dpi = 200)
}

# Pseudobulk is the inferential unit for disease contrast; cells are not treated as independent replicates.
counts <- tryCatch(Seurat::GetAssayData(filtered, assay = "RNA", layer = "counts"), error = function(e) Seurat::GetAssayData(filtered, assay = "RNA", slot = "counts"))
md <- filtered[[]]
primary_types <- as.character(config$primary_cell_types)
keep_cells <- md$cell_type %in% primary_types
if (sum(keep_cells) == 0) stop("No primary glial cells were assigned; inspect D3_cluster_marker_scores.tsv")
md_glia <- md[keep_cells, , drop = FALSE]
counts_glia <- counts[, rownames(md_glia), drop = FALSE]
group_key <- paste(md_glia$donor_id, md_glia$region, md_glia$cell_type, sep = "__")
group_sizes <- table(group_key)
valid_groups <- names(group_sizes)[group_sizes >= as.numeric(config$min_cells_per_pseudobulk)]
if (length(valid_groups) < 2) stop("Too few glial pseudobulk groups after cell-count filter")

pb <- do.call(cbind, lapply(valid_groups, function(g) Matrix::rowSums(counts_glia[, which(group_key == g), drop = FALSE])))
pb <- Matrix::Matrix(pb, sparse = TRUE)
colnames(pb) <- valid_groups
pb_meta <- do.call(rbind, lapply(valid_groups, function(g) {
  idx <- which(group_key == g)[[1]]
  data.frame(
    group_id = g,
    donor_id = md_glia$donor_id[[idx]],
    disease = md_glia$disease[[idx]],
    region = md_glia$region[[idx]],
    cell_type = md_glia$cell_type[[idx]],
    age = md_glia$age[[idx]],
    sex = md_glia$sex[[idx]],
    rin = md_glia$rin[[idx]],
    ptau_percent_positive_cells = md_glia$ptau_percent_positive_cells[[idx]],
    amyloid_percent_area = md_glia$amyloid_percent_area[[idx]],
    n_cells = as.integer(group_sizes[[g]]),
    stringsAsFactors = FALSE
  )
}))
rownames(pb_meta) <- pb_meta$group_id
write.table(pb_meta, file.path(tables_dir, "D3_pseudobulk_metadata.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(as.matrix(pb), file.path(tables_dir, "D3_pseudobulk_counts.tsv"), sep = "\t", quote = FALSE, col.names = NA)

run_de <- function(cell_type) {
  sel <- which(pb_meta$cell_type == cell_type)
  meta_ct <- pb_meta[sel, , drop = FALSE]
  if (length(unique(meta_ct$disease)) < 2 || min(table(meta_ct$disease)) < as.numeric(config$min_replicate_groups)) {
    return(data.frame())
  }
  dat <- pb[, sel, drop = FALSE]
  meta_ct$disease <- relevel(factor(meta_ct$disease), ref = "Control")
  meta_ct$region <- factor(meta_ct$region)
  meta_ct$sex <- factor(meta_ct$sex)
  meta_ct$age_c <- as.numeric(scale(meta_ct$age))
  meta_ct$rin_c <- as.numeric(scale(meta_ct$rin))
  design <- model.matrix(~ disease + region + age_c + sex + rin_c, data = meta_ct)
  keep_genes <- Matrix::rowSums(dat >= as.numeric(config$min_count_genes)) >= as.numeric(config$min_replicate_groups)
  dat <- dat[keep_genes, , drop = FALSE]
  y <- edgeR::DGEList(counts = dat, samples = meta_ct)
  y <- edgeR::calcNormFactors(y)
  v <- limma::voom(y, design, plot = FALSE)
  corfit <- tryCatch(limma::duplicateCorrelation(v, design = design, block = meta_ct$donor_id), error = function(e) NULL)
  if (is.null(corfit) || !is.finite(corfit$consensus)) {
    fit <- limma::lmFit(v, design)
  } else {
    fit <- limma::lmFit(v, design, block = meta_ct$donor_id, correlation = corfit$consensus)
  }
  disease_coef <- grep("^diseaseAD$", colnames(design))
  if (length(disease_coef) != 1) stop("Could not identify AD coefficient")
  fit <- limma::eBayes(fit)
  tt <- limma::topTable(fit, coef = disease_coef, number = Inf, sort.by = "none")
  out <- data.frame(
    gene = rownames(tt),
    logFC = tt$logFC,
    t = tt$t,
    P.Value = tt$P.Value,
    FDR = tt$adj.P.Val,
    B = tt$B,
    cell_type = cell_type,
    n_pseudobulk_groups = nrow(meta_ct),
    n_AD_groups = sum(meta_ct$disease == "AD"),
    n_Control_groups = sum(meta_ct$disease == "Control"),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  out$direction <- ifelse(out$logFC > 0, "AD_up", "AD_down")
  out[order(out$FDR, -abs(out$t)), ]
}

signatures <- lapply(primary_types, run_de)
names(signatures) <- primary_types
for (cell_type in primary_types) {
  out <- signatures[[cell_type]]
  if (nrow(out) == 0) next
  file_name <- ifelse(cell_type == "Microglia", "DiseaseSignature_microglia.tsv", "DiseaseSignature_astrocyte.tsv")
  write.table(out, file.path(tables_dir, file_name), sep = "\t", quote = FALSE, row.names = FALSE)
}

summary <- list(
  dataset = "D3_GSE160936",
  samples = nrow(metadata),
  donors = length(unique(metadata$donor_id)),
  cells_before_qc = ncol(combined),
  cells_after_qc = ncol(filtered),
  qc_retention = ncol(filtered) / ncol(combined),
  clusters = length(unique(cluster_ids)),
  marker_sets = lapply(marker_sets, length),
  cell_type_counts = as.list(table(filtered$cell_type)),
  pseudobulk_groups = nrow(pb_meta),
  pseudobulk_groups_by_type = as.list(table(pb_meta$cell_type)),
  de_outputs = names(signatures)[vapply(signatures, nrow, integer(1)) > 0],
  qc_note = "Sample-wise empirical upper quantiles; min_features configured; no doublet caller available in public filtered matrix.",
  disease_effect = "AD minus Control on donor-region-celltype pseudobulk; region, age, sex, RIN covariates; donor block via duplicateCorrelation when estimable."
)
jsonlite::write_json(summary, file.path(interim_dir, "D3_discovery_summary.json"), auto_unbox = TRUE, pretty = TRUE)
saveRDS(filtered, file.path(processed_dir, "D3_GSE160936_glia_seurat.rds"), compress = TRUE)
message("D3 discovery completed: ", ncol(filtered), " cells after QC; ", nrow(pb_meta), " pseudobulk groups")