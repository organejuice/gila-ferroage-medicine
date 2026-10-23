#!/usr/bin/env Rscript

# G6-FP-2026-09-27-001, WP2. Builds auditable per-cell matrices from the
# direct CellBin GEF cellExp layer; does not use cellTypeID or diagnosis.
suppressPackageStartupMessages({ library(rhdf5); library(Matrix); library(data.table); library(digest) })
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: 91_build_g6_cellbin_nucleus_counts.R <project_root>")
root <- args[[1L]]
raw_cell_dir <- file.path(root, "data/raw/GSE269906_cellbin")
raw_tissue_dir <- file.path(root, "data/raw/GSE269906_tissue")
audit_dir <- file.path(root, "results/2026-09-24-disease-side-v3.5-g6-input-audit")
cell_manifest_path <- file.path(audit_dir, "registries/GSE269906_file_acquisition_manifest.tsv")
tissue_qc_path <- file.path(audit_dir, "tables/GSE269906_tissue_bin50_section_qc.tsv")
cell_struct_path <- file.path(audit_dir, "tables/GSE269906_cellbin_structural_gene_coverage.tsv.gz")
if (!all(file.exists(c(cell_manifest_path, tissue_qc_path, cell_struct_path)))) stop("Prior locked input audit is incomplete")
out_root <- file.path(root, "data/interim/G6_next_feasibility_001/cellbin")
result_dir <- file.path(root, "results/2026-09-27-g6-next-feasibility/cellbin")
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
if (length(list.files(out_root, pattern = "_nucleus_counts\\.rds$"))) stop("Existing CellBin outputs found; refusing overwrite")

# Whitelist columns at read time; disease/pathology metadata never enters this package.
acq <- fread(cell_manifest_path, select = c("GSM", "donor_id", "file", "expected_bytes", "actual_bytes", "status", "sha256"))
acq <- acq[grepl("\\.cellbin\\.gef$", file)]
tissue_qc <- fread(tissue_qc_path, select = c("GSM", "donor_id", "source_file", "source_sha256", "source_bytes",
                                              "total_UMI", "n_bin1_expression_records", "coordinate_resolution_nm", "minX", "minY"))
cell_struct <- fread(cell_struct_path, select = c("GSM", "donor_id", "file", "n_cellbin_rows", "n_unique_gene_symbols",
                                                  "GEF_resolution_attribute", "offsetX", "offsetY"))
if (nrow(acq) != 12L || uniqueN(acq$GSM) != 12L || nrow(tissue_qc) != 12L || uniqueN(tissue_qc$GSM) != 12L) stop("Expected 12 uniquely mapped paired donors")
if (!all(acq$status == "complete") || any(acq$actual_bytes != acq$expected_bytes)) stop("CellBin acquisition manifest is incomplete")
setkey(acq, GSM, donor_id)
setkey(tissue_qc, GSM, donor_id)
setkey(cell_struct, GSM, donor_id)
cell_files <- list.files(raw_cell_dir, pattern = "\\.cellbin\\.gef$", full.names = TRUE)
tissue_files <- list.files(raw_tissue_dir, pattern = "\\.tissue\\.gef$", full.names = TRUE)
if (length(cell_files) != 12L || length(tissue_files) != 12L) stop("Expected 12 raw GEF files per modality")
parse_pair <- function(path, suffix) {
  base <- basename(path)
  m <- regexec("^(GSM[0-9]+)_([^.]+)\\.(cellbin|tissue)\\.gef$", base, perl = TRUE)
  z <- regmatches(base, m)[[1L]]
  if (length(z) != 4L || z[[4L]] != suffix) stop("Unrecognized raw GEF name: ", base)
  list(GSM = z[[2L]], donor_id = z[[3L]], file = base)
}
cell_info <- rbindlist(lapply(cell_files, parse_pair, suffix = "cellbin"))
tissue_info <- rbindlist(lapply(tissue_files, parse_pair, suffix = "tissue"))
if (uniqueN(cell_info, by = c("GSM", "donor_id")) != 12L || uniqueN(tissue_info, by = c("GSM", "donor_id")) != 12L) stop("Duplicate GSM/donor file crosswalk")
if (!setequal(paste(cell_info$GSM, cell_info$donor_id), paste(tissue_info$GSM, tissue_info$donor_id))) stop("Tissue/CellBin donor crosswalk mismatch")
if (!all(paste(cell_info$GSM, cell_info$donor_id) %in% paste(acq$GSM, acq$donor_id))) stop("Raw files do not match the locked donor manifest")

qc_rows <- vector("list", 12L)
input_rows <- vector("list", 24L)
for (k in seq_len(nrow(cell_info))) {
  ci <- cell_info[k]
  ti <- tissue_info[GSM == ci$GSM & donor_id == ci$donor_id]
  ac <- acq[GSM == ci$GSM & donor_id == ci$donor_id]
  tq <- tissue_qc[GSM == ci$GSM & donor_id == ci$donor_id]
  cs <- cell_struct[GSM == ci$GSM & donor_id == ci$donor_id]
  if (nrow(ti) != 1L || nrow(ac) != 1L || nrow(tq) != 1L || nrow(cs) != 1L) stop("Non-unique crosswalk for ", ci$GSM)
  cell_path <- file.path(raw_cell_dir, ci$file)
  tissue_path <- file.path(raw_tissue_dir, ti$file)
  cell_bytes <- as.numeric(file.info(cell_path)$size)
  tissue_bytes <- as.numeric(file.info(tissue_path)$size)
  cell_sha <- digest(cell_path, algo = "sha256", file = TRUE, serialize = FALSE)
  tissue_sha <- digest(tissue_path, algo = "sha256", file = TRUE, serialize = FALSE)
  if (cell_bytes != ac$expected_bytes || cell_bytes != ac$actual_bytes || cell_sha != ac$sha256) stop("CellBin raw file hash/size mismatch: ", ci$file)
  if (tissue_bytes != tq$source_bytes || tissue_sha != tq$source_sha256) stop("Tissue raw file hash/size mismatch: ", ti$file)
  input_rows[[2L * k - 1L]] <- data.frame(GSM = ci$GSM, donor_id = ci$donor_id, modality = "CellBin",
    file = ci$file, bytes = cell_bytes, sha256 = cell_sha, stringsAsFactors = FALSE)
  input_rows[[2L * k]] <- data.frame(GSM = ti$GSM, donor_id = ti$donor_id, modality = "tissue_bin1",
    file = ti$file, bytes = tissue_bytes, sha256 = tissue_sha, stringsAsFactors = FALSE)

  tmp <- tempfile(pattern = "g6_cellbin_", tmpdir = Sys.getenv("TEMP"), fileext = ".gef")
  if (grepl("[^ -~]", tmp) || !file.copy(cell_path, tmp, overwrite = FALSE)) stop("Cannot stage CellBin GEF at ASCII temp path")
  root_attrs <- h5readAttributes(tmp, "/")
  paths <- h5ls(tmp, recursive = TRUE)
  required <- c("cell", "cellBorder", "cellExp", "gene", "geneExp")
  if (!all(required %in% paths$name[paths$group == "/cellBin"])) stop("Direct CellBin expression/geometry datasets missing: ", ci$file)
  if (!identical(as.character(root_attrs$omics), "Transcriptomics") || as.numeric(root_attrs$resolution)[1L] != tq$coordinate_resolution_nm ||
      as.numeric(root_attrs$resolution)[1L] != 500 || as.numeric(root_attrs$offsetX)[1L] != tq$minX ||
      as.numeric(root_attrs$offsetY)[1L] != tq$minY) stop("CellBin/tissue coordinate transform is not exact: ", ci$file)
  cell <- h5read(tmp, "/cellBin/cell")
  expression <- h5read(tmp, "/cellBin/cellExp")
  gene <- h5read(tmp, "/cellBin/gene")
  gene_exp <- h5read(tmp, "/cellBin/geneExp")
  borders <- h5read(tmp, "/cellBin/cellBorder")
  n_cell <- nrow(cell); n_gene <- nrow(gene)
  if (n_cell != cs$n_cellbin_rows || nrow(expression) != sum(cell$geneCount) ||
      nrow(gene_exp) != sum(gene$cellCount) || n_gene != cs$n_unique_gene_symbols) stop("Direct CellBin matrix dimensions disagree with audited catalog")
  if (anyDuplicated(cell$id) || !identical(as.integer(cell$id), seq.int(0L, n_cell - 1L))) stop("Cell IDs are not unique ordered zero-based keys")
  if (anyDuplicated(gene$geneName) || anyNA(gene$geneName) || any(!nzchar(gene$geneName))) stop("CellBin gene catalog has ambiguous symbols")
  expected_cell_offset <- c(0, head(cumsum(cell$geneCount), -1L))
  expected_gene_offset <- numeric(n_gene)
  active_genes <- which(gene$cellCount > 0L)
  expected_gene_offset[active_genes] <- c(0, head(cumsum(gene$cellCount[active_genes]), -1L))
  if (any(cell$offset != expected_cell_offset) || any(gene$offset != expected_gene_offset)) stop("Cell/ gene offsets are not contiguous")
  if (any(expression$geneID < 0L | expression$geneID >= n_gene) || any(expression$count <= 0L) ||
      any(gene_exp$cellID < 0L | gene_exp$cellID >= n_cell) || any(gene_exp$count <= 0L)) stop("CellBin index/count values out of range")
  cell_ix <- rep.int(seq_len(n_cell), as.integer(cell$geneCount))
  if (length(cell_ix) != nrow(expression)) stop("cellExp records do not map one-to-one to cell offsets")
  cell_key <- as.double(cell_ix - 1L) * n_gene + as.double(expression$geneID)
  if (anyDuplicated(cell_key)) stop("Duplicate cell-by-gene records in direct cellExp layer")
  cell_totals <- rowsum(as.numeric(expression$count), cell_ix, reorder = TRUE)
  cell_sum <- numeric(n_cell)
  cell_sum[as.integer(rownames(cell_totals))] <- cell_totals[, 1L]
  gene_totals <- rowsum(as.numeric(expression$count), as.integer(expression$geneID) + 1L, reorder = TRUE)
  gene_sum <- numeric(n_gene)
  gene_sum[as.integer(rownames(gene_totals))] <- gene_totals[, 1L]
  gene_n <- tabulate(expression$geneID + 1L, nbins = n_gene)
  gene_catalog_modulus_ok <- all(as.numeric(gene$expCount) == (gene_sum %% 65536))
  if (!identical(as.numeric(cell_sum), as.numeric(cell$expCount)) ||
      !identical(as.integer(gene_n), as.integer(gene$cellCount)) || !gene_catalog_modulus_ok) {
    stop("CellExp marginals disagree with CellBin cell/gene indices or encoded gene totals")
  }
  # CellExp and geneExp are transposed sparse maps. Compare every cell-gene pair
  # exactly; GEF's gene.expCount summary is UInt16 and wraps after 65,535.
  active_genes <- which(gene$cellCount > 0L)
  gene_ix <- rep.int(active_genes - 1L, as.integer(gene$cellCount[active_genes]))
  if (length(gene_ix) != nrow(gene_exp)) stop("geneExp records do not map one-to-one to gene offsets")
  cell_key_gene <- as.double(gene_exp$cellID) * n_gene + as.double(gene_ix)
  if (anyDuplicated(cell_key_gene)) stop("Duplicate cell-by-gene records in transposed geneExp layer")
  order_cell <- order(cell_key)
  order_gene <- order(cell_key_gene)
  if (!identical(cell_key[order_cell], cell_key_gene[order_gene]) ||
      any(as.integer(expression$count)[order_cell] != as.integer(gene_exp$count)[order_gene])) {
    stop("cellExp and geneExp differ at one or more cell-by-gene records")
  }
  cell_totals_inverse <- rowsum(as.numeric(gene_exp$count), as.integer(gene_exp$cellID) + 1L, reorder = TRUE)
  cell_umi_inverse <- numeric(n_cell)
  cell_umi_inverse[as.integer(rownames(cell_totals_inverse))] <- cell_totals_inverse[, 1L]
  cell_n_inverse <- tabulate(gene_exp$cellID + 1L, nbins = n_cell)
  if (!identical(as.numeric(cell_umi_inverse), as.numeric(cell$expCount)) ||
      !identical(as.integer(cell_n_inverse), as.integer(cell$geneCount))) stop("GeneExp inverse cell marginals disagree with cell table")
  valid_vertex <- borders[1L, , ] != 32767 & borders[2L, , ] != 32767
  valid_polygon <- colSums(valid_vertex) >= 3L & cell$area > 0
  if (length(valid_polygon) != n_cell || mean(valid_polygon) < 0.99) stop("CellBin segmentation geometry is invalid")
  assigned_umi <- sum(expression$count)
  raw_umi <- as.numeric(tq$total_UMI)
  if (assigned_umi != sum(cell$expCount) || assigned_umi > raw_umi) stop("CellBin UMI accounting does not balance")
  cell_names <- paste(ci$GSM, sprintf("cell%08d", cell$id), sep = "::")
  counts <- sparseMatrix(i = as.integer(expression$geneID) + 1L, j = cell_ix,
                         x = as.integer(expression$count), dims = c(n_gene, n_cell),
                         dimnames = list(gene$geneName, cell_names), giveCsparse = TRUE)
  if (!identical(rownames(counts), gene$geneName) || !identical(colnames(counts), cell_names)) {
    stop("Materialized matrix lost its gene or nucleus identifiers: ", ci$GSM)
  }
  meta <- data.frame(cell_id = as.integer(cell$id), GSM = ci$GSM, donor_id = ci$donor_id,
    x_dnb = as.numeric(cell$x), y_dnb = as.numeric(cell$y),
    x_um = as.numeric(cell$x) * 0.5, y_um = as.numeric(cell$y) * 0.5,
    n_genes = as.integer(cell$geneCount), n_umi = as.numeric(cell$expCount),
    dnb_count = as.integer(cell$dnbCount), area_dnb2 = as.numeric(cell$area),
    polygon_valid = as.logical(valid_polygon), row.names = cell_names, stringsAsFactors = FALSE)
  if (sum(counts) != assigned_umi || any(Matrix::colSums(counts) != meta$n_umi) ||
      any(Matrix::colSums(counts > 0) != meta$n_genes)) stop("Materialized sparse matrix fails cell-level UMI/geneCount checks")
  out <- list(study_id = "G6-FP-2026-09-27-001", GSM = ci$GSM, donor_id = ci$donor_id,
              counts = counts, cell_metadata = meta,
              provenance = list(cellbin_file = ci$file, tissue_file = ti$file,
                cellbin_sha256 = cell_sha, tissue_sha256 = tissue_sha,
                cellbin_resolution_nm = as.numeric(root_attrs$resolution)[1L],
                cellbin_offsetX = as.numeric(root_attrs$offsetX)[1L], cellbin_offsetY = as.numeric(root_attrs$offsetY)[1L],
                tissue_minX = as.numeric(tq$minX), tissue_minY = as.numeric(tq$minY),
                assignment_layer = "/cellBin/cellExp; direct cellID/geneID mapping",
                cell_type_annotation = "cellTypeID excluded; all default-zero labels ignored",
                author_pipeline_status = "AUTHOR_PIPELINE_NOT_REPRODUCIBLE"))
  out_path <- file.path(out_root, paste0(ci$GSM, "_nucleus_counts.rds"))
  tmp_out <- paste0(out_path, ".tmp")
  saveRDS(out, tmp_out, compress = "gzip")
  if (!file.rename(tmp_out, out_path)) stop("Cannot atomically write nucleus matrix")
  qc_rows[[k]] <- data.frame(GSM = ci$GSM, donor_id = ci$donor_id, n_cells_total = n_cell,
    n_cells_ge150_genes = sum(cell$geneCount >= 150L), fraction_cells_ge150_genes = mean(cell$geneCount >= 150L),
    median_geneCount = median(cell$geneCount), median_expCount = median(cell$expCount),
    n_valid_polygons = sum(valid_polygon), fraction_valid_polygons = mean(valid_polygon),
    n_cellExp_records = nrow(expression), n_geneExp_records = nrow(gene_exp),
    cellExp_UMI = assigned_umi, tissue_UMI_from_locked_QC = raw_umi,
    non_CellBin_residual_UMI = raw_umi - assigned_umi,
    fraction_tissue_UMI_in_CellBin = assigned_umi / raw_umi,
    gene_catalog_expCount_uint16_overflow_genes = sum(gene_sum > 65535L),
    gene_catalog_expCount_overflow_gene_symbols = paste(gene$geneName[gene_sum > 65535L], collapse = ";"),
    cell_count_sum_matches_metadata = TRUE, gene_count_sum_matches_catalog = TRUE,
    cell_ID_gene_ID_crosswalk_unique = TRUE, polygon_coordinate_match = TRUE,
    histology_validation = "UNAVAILABLE_IN_PRIOR_AUDIT", author_pipeline_status = "AUTHOR_PIPELINE_NOT_REPRODUCIBLE",
    stringsAsFactors = FALSE)
  message("CELLBIN_COUNTS ", ci$GSM, "; cells=", n_cell, "; >=150 genes=", sum(cell$geneCount >= 150L),
          "; assigned fraction=", signif(assigned_umi / raw_umi, 4))
  unlink(tmp)
  rm(cell, expression, gene, gene_exp, borders, counts, meta, out, cell_ix, pair_key, cell_sum, gene_sum,
     gene_n, cell_umi_inverse, cell_n_inverse)
  gc(verbose = FALSE)
}
qc <- rbindlist(qc_rows, fill = TRUE)
if (nrow(qc) != 12L || uniqueN(qc$donor_id) != 12L) stop("Not all 12 donor matrices completed")
fwrite(qc, file.path(result_dir, "cellbin_nucleus_counts_qc.tsv"), sep = "\t")
fwrite(rbindlist(input_rows), file.path(result_dir, "input_sha256_manifest.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(result_dir, "R_sessionInfo.txt"))
message("WP2 direct count build complete; donors=12; nuclei after >=150 gene screen=", sum(qc$n_cells_ge150_genes))
