#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(Matrix)
  library(hdf5r)
})

root <- "."
d7_raw <- file.path("data", "raw", "D7_GSE148822")
d7_interim <- file.path("data", "interim", "D7_GSE148822")
d8_raw <- file.path("data", "raw", "D8_GSE174367")
d8_interim <- file.path("data", "interim", "D8_GSE174367")
dir.create(d7_interim, recursive = TRUE, showWarnings = FALSE)
dir.create(d8_interim, recursive = TRUE, showWarnings = FALSE)

micro_markers <- c("AIF1", "CX3CR1", "P2RY12", "TMEM119", "C1QA", "C1QB", "C1QC", "TYROBP", "LST1", "CSF1R", "SPI1")
astro_markers <- c("ALDH1L1", "AQP4", "GFAP", "SLC1A2", "SLC1A3", "SOX9", "GJA1")

read_features <- function(path) {
  x <- read.delim(gzfile(path), header = FALSE, sep = "\t", quote = "", comment.char = "", check.names = FALSE)
  gene <- if (ncol(x) >= 2) as.character(x[[2]]) else as.character(x[[1]])
  first <- as.character(x[[1]])
  gene[is.na(gene) | gene == ""] <- first[is.na(gene) | gene == ""]
  gene <- sub("\\.[0-9]+$", "", gene)
  make.unique(gene)
}

marker_rows <- function(genes, markers) {
  g <- toupper(sub("\\.[0-9]+$", "", genes))
  which(g %in% toupper(markers))
}

safe_median <- function(x) if (length(x)) median(x, na.rm = TRUE) else NA_real_

write_pseudobulk <- function(count_list, meta, out_counts, out_meta, out_qc) {
  if (!length(count_list)) stop("No pseudobulk samples were generated")
  common_genes <- Reduce(intersect, lapply(count_list, names))
  all_genes <- unique(unlist(lapply(count_list, names), use.names = FALSE))
  # Preserve the first-seen feature order while retaining all genes.
  mat <- matrix(0, nrow = length(all_genes), ncol = length(count_list), dimnames = list(all_genes, names(count_list)))
  for (j in seq_along(count_list)) mat[match(names(count_list[[j]]), all_genes), j] <- as.numeric(count_list[[j]])
  con <- gzfile(out_counts, open = "wt")
  write.table(data.frame(gene = rownames(mat), mat, check.names = FALSE), con, sep = "\t", quote = FALSE, row.names = FALSE)
  close(con)
  write.table(meta, out_meta, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
  write.table(meta, out_qc, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
  invisible(list(n_genes = nrow(mat), n_samples = ncol(mat), n_common_genes = length(common_genes)))
}

prepare_d7 <- function() {
  manifest <- read.delim(file.path(d7_interim, "D7_sample_manifest.tsv"), sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
  manifest <- manifest[manifest$include_primary == "true", , drop = FALSE]
  count_list <- list()
  meta_rows <- list()
  qc_rows <- list()
  for (k in seq_len(nrow(manifest))) {
    row <- manifest[k, ]
    message(sprintf("D7 %d/%d: %s %s", k, nrow(manifest), row$donor, row$region))
    feat <- read_features(row$local_feature)
    mat <- as(Matrix::readMM(gzfile(row$local_matrix)), "dgCMatrix")
    if (nrow(mat) != length(feat)) stop("D7 feature/matrix dimension mismatch: ", row$accession)
    rownames(mat) <- feat
    total <- Matrix::colSums(mat)
    nfeature <- Matrix::colSums(mat > 0)
    mt <- which(grepl("^(MT-|MT\\.)", toupper(feat)))
    mt_pct <- if (length(mt)) 100 * Matrix::colSums(mat[mt, , drop = FALSE]) / pmax(total, 1) else rep(0, ncol(mat))
    mi <- marker_rows(feat, micro_markers)
    ai <- marker_rows(feat, astro_markers)
    if (!length(mi)) stop("No microglia marker genes found in D7 features: ", row$accession)
    micro_detect <- Matrix::colSums(mat[mi, , drop = FALSE] > 0)
    astro_detect <- if (length(ai)) Matrix::colSums(mat[ai, , drop = FALSE] > 0) else rep(0, ncol(mat))
    micro_score <- Matrix::colSums(mat[mi, , drop = FALSE]) / pmax(total, 1) * 10000
    astro_score <- if (length(ai)) Matrix::colSums(mat[ai, , drop = FALSE]) / pmax(total, 1) * 10000 else rep(0, ncol(mat))
    qc <- total >= 500 & nfeature >= 200 & mt_pct <= 20
    selected <- qc & micro_detect >= 2 & micro_score >= astro_score
    group_id <- paste(row$donor, row$region, "Microglia", sep = "__")
    n_used <- sum(selected)
    if (n_used < 20) warning("D7 donor-region below frozen 20-cell threshold: ", group_id)
    count_list[[group_id]] <- Matrix::rowSums(mat[, selected, drop = FALSE])
    meta_rows[[group_id]] <- data.frame(
      group_id = group_id, cohort = "D7_GSE148822", dataset_id = "GSE148822", donor_id = row$donor,
      disease = ifelse(row$sample_group == "AD", "AD", "Control"), region = row$region,
      cell_type = "Microglia", age = as.numeric(row$age), sex = row$sex, RIN = as.numeric(row$RIN),
      PMI = NA_real_, batch = NA_character_, n_cells_total = ncol(mat), n_cells_qc = sum(qc),
      n_cells_microglia = n_used, n_cells_used = n_used, median_nCount = safe_median(total[selected]),
      median_nFeature = safe_median(nfeature[selected]), median_mt_pct = safe_median(mt_pct[selected]),
      annotation_source = "marker_based_from_GEO_filtered_Cell_Ranger_MTX",
      qc_rule = "GEO_filtered; nCount>=500; nFeature>=200; mt_pct<=20; >=2 microglia markers; micro_score>=astro_score",
      doublet_detection = "not_run_scDblFinder_unavailable", stringsAsFactors = FALSE
    )
    qc_rows[[group_id]] <- data.frame(accession = row$accession, title = row$title, group_id = group_id,
      cohort = "D7_GSE148822", donor_id = row$donor, region = row$region, diagnosis = row$sample_group,
      n_cells_total = ncol(mat), n_cells_qc = sum(qc), n_cells_microglia = n_used,
      median_nCount_all = safe_median(total), median_nFeature_all = safe_median(nfeature), median_mt_pct_all = safe_median(mt_pct),
      min_nCount = min(total), max_nCount = max(total), min_nFeature = min(nfeature), max_nFeature = max(nfeature),
      marker_rule = "micro_detect>=2 and micro_score>=astro_score", stringsAsFactors = FALSE)
    rm(mat)
    gc(verbose = FALSE)
  }
  meta <- do.call(rbind, meta_rows)
  rownames(meta) <- NULL
  qc <- do.call(rbind, qc_rows)
  rownames(qc) <- NULL
  out <- write_pseudobulk(count_list, meta,
    file.path(d7_interim, "D7_pseudobulk_counts_microglia.tsv.gz"),
    file.path(d7_interim, "D7_pseudobulk_metadata_microglia.tsv"),
    file.path(d7_interim, "D7_qc_summary_microglia.tsv"))
  print(out)
}

prepare_d8 <- function() {
  meta <- read.csv(gzfile(file.path(d8_raw, "GSE174367_snRNA-seq_cell_meta.csv.gz")), check.names = FALSE, stringsAsFactors = FALSE)
  # rhdf5 on this Windows installation cannot open non-ASCII absolute paths;
  # use the project-relative path while keeping all outputs under the project.
  h5 <- file.path("data", "raw", "D8_GSE174367", "GSE174367_snRNA-seq_filtered_feature_bc_matrix.h5")
  hf <- H5File$new(h5, mode = "r")
  barcodes <- as.character(hf[["matrix/barcodes"]][])
  genes <- as.character(hf[["matrix/features/name"]][])
  shape <- as.integer(hf[["matrix/shape"]][])
  if (length(barcodes) != shape[2]) stop("D8 barcode count does not match H5 shape")
  order_meta <- match(barcodes, meta$Barcode)
  unmatched <- is.na(order_meta)
  if (any(unmatched)) message("D8 excluding ", sum(unmatched), " H5 barcodes without GEO cell metadata")
  valid_barcode <- !unmatched
  barcodes <- barcodes[valid_barcode]
  meta <- meta[order_meta[valid_barcode], , drop = FALSE]
  if (anyDuplicated(meta$Barcode)) stop("D8 metadata barcodes are duplicated")
  ind <- as.integer(hf[["matrix/indices"]][])+1L
  # 10x HDF5 stores CSC column pointers as 0-based cumulative offsets;
  # Matrix::sparseMatrix expects the same p convention and 1-based row indices.
  ptr <- as.integer(hf[["matrix/indptr"]][])
  dat <- as.numeric(hf[["matrix/data"]][])
  hf$close_all()
  if (length(ind) != length(dat) || length(ptr) != shape[2] + 1L) stop("D8 H5 sparse arrays are inconsistent")
  mat <- sparseMatrix(i = ind, p = ptr, x = dat, dims = shape, repr = "C")
  if (any(unmatched)) mat <- mat[, valid_barcode, drop = FALSE]
  rownames(mat) <- make.unique(sub("\\.[0-9]+$", "", genes))
  mg <- which(meta$Cell.Type == "MG")
  if (!length(mg)) stop("D8 source annotation has no MG cells")
  mgmat <- mat[, mg, drop = FALSE]
  total <- Matrix::colSums(mgmat)
  nfeature <- Matrix::colSums(mgmat > 0)
  mt <- which(grepl("^(MT-|MT\\.)", toupper(rownames(mgmat))))
  mt_pct <- if (length(mt)) 100 * Matrix::colSums(mgmat[mt, , drop = FALSE]) / pmax(total, 1) else rep(0, ncol(mgmat))
  qc <- total >= 500 & nfeature >= 200 & mt_pct <= 20
  count_list <- list(); meta_rows <- list(); qc_rows <- list()
  for (s in unique(meta$SampleID[mg])) {
    local <- which(meta$SampleID[mg] == s)
    selected <- local[qc[local]]
    group_id <- paste(s, "PFC", "Microglia", sep = "__")
    count_list[[group_id]] <- Matrix::rowSums(mgmat[, selected, drop = FALSE])
    source <- meta[mg[local[1]], , drop = FALSE]
    disease <- ifelse(source$Diagnosis == "AD", "AD", "Control")
    meta_rows[[group_id]] <- data.frame(
      group_id = group_id, cohort = "D8_GSE174367", dataset_id = "GSE174367", donor_id = s, disease = disease,
      region = "PFC", cell_type = "Microglia", age = as.numeric(source$Age), sex = source$Sex,
      RIN = as.numeric(source$RIN), PMI = as.numeric(source$PMI), batch = source$Batch,
      n_cells_total = length(local), n_cells_qc = sum(qc[local]), n_cells_microglia = length(local),
      n_cells_used = length(selected), median_nCount = safe_median(total[selected]),
      median_nFeature = safe_median(nfeature[selected]), median_mt_pct = safe_median(mt_pct[selected]),
      annotation_source = "GEO_cell_meta_Cell.Type_MG",
      qc_rule = "GEO_filtered; Cell.Type==MG; nCount>=500; nFeature>=200; mt_pct<=20",
      doublet_detection = "not_run_scDblFinder_unavailable", stringsAsFactors = FALSE
    )
    qc_rows[[group_id]] <- data.frame(sample_id = s, group_id = group_id, cohort = "D8_GSE174367",
      diagnosis = source$Diagnosis, n_cells_total = length(local), n_cells_qc = sum(qc[local]),
      n_cells_microglia = length(local), n_cells_used = length(selected), median_nCount_all = safe_median(total[local]),
      median_nFeature_all = safe_median(nfeature[local]), median_mt_pct_all = safe_median(mt_pct[local]),
      min_nCount = min(total[local]), max_nCount = max(total[local]), min_nFeature = min(nfeature[local]), max_nFeature = max(nfeature[local]),
      stringsAsFactors = FALSE)
    if (length(selected) < 20) warning("D8 donor below frozen 20-cell threshold: ", group_id)
  }
  meta_out <- do.call(rbind, meta_rows); rownames(meta_out) <- NULL
  qc_out <- do.call(rbind, qc_rows); rownames(qc_out) <- NULL
  out <- write_pseudobulk(count_list, meta_out,
    file.path(d8_interim, "D8_pseudobulk_counts_microglia.tsv.gz"),
    file.path(d8_interim, "D8_pseudobulk_metadata_microglia.tsv"),
    file.path(d8_interim, "D8_qc_summary_microglia.tsv"))
  print(out)
  write.table(data.frame(barcode = barcodes, SampleID = meta$SampleID, Diagnosis = meta$Diagnosis, Cell.Type = meta$Cell.Type),
    file.path(d8_interim, "D8_barcode_metadata_mapping.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  rm(mat, mgmat, ind, ptr, dat)
  gc(verbose = FALSE)
}

args <- commandArgs(trailingOnly = TRUE)
if (!length(args) || "all" %in% args || "d7" %in% args) prepare_d7()
if (!length(args) || "all" %in% args || "d8" %in% args) prepare_d8()
