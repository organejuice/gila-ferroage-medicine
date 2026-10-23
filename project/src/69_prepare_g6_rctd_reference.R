#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Matrix)
  library(rhdf5)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) {
  stop("Usage: 69_prepare_g6_rctd_reference.R <project_root> <staged_10x_h5> <metadata_csv_gz> <output_dir> <chunk_cells>")
}
project_root <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
h5_path <- normalizePath(args[[2]], winslash = "/", mustWork = TRUE)
metadata_path <- normalizePath(args[[3]], winslash = "/", mustWork = TRUE)
output_dir <- args[[4]]
chunk_cells <- as.integer(args[[5]])
if (chunk_cells < 100L) stop("chunk_cells must be >=100")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

meta <- read.csv(gzfile(metadata_path), stringsAsFactors = FALSE, check.names = FALSE)
required <- c("Barcode", "SampleID", "Cell.Type")
if (!all(required %in% names(meta))) stop("Reference metadata are missing required fields")
barcodes <- as.character(h5read(h5_path, "/matrix/barcodes"))
genes <- as.character(h5read(h5_path, "/matrix/features/name"))
genes <- toupper(genes)
gene_levels <- unique(genes)
gene_row <- match(genes, gene_levels)
indptr <- as.numeric(h5read(h5_path, "/matrix/indptr"))
n_cells <- length(barcodes)
if (length(indptr) != n_cells + 1L || length(genes) < 1000L) {
  stop("Unexpected 10x HDF5 dimensions")
}
meta_i <- match(barcodes, meta$Barcode)
unmapped_barcodes <- sum(is.na(meta_i))
if (unmapped_barcodes / length(barcodes) > 0.01) {
  stop("More than 1% of matrix barcodes are absent from source metadata")
}
meta <- meta[meta_i, , drop = FALSE]
if (anyDuplicated(meta$Barcode[!is.na(meta$Barcode)])) {
  stop("Duplicate reference barcode metadata")
}

# First streamed pass computes UMI totals without materializing the 177M-entry
# matrix. It deliberately never reads or uses the Diagnosis column.
umi <- numeric(n_cells)
chunk_starts <- seq.int(1L, n_cells, by = chunk_cells)
for (start_cell in chunk_starts) {
  end_cell <- min(n_cells, start_cell + chunk_cells - 1L)
  ptr <- indptr[start_cell:(end_cell + 1L)]
  n_nz <- ptr[length(ptr)] - ptr[[1L]]
  if (n_nz <= 0) next
  values <- as.numeric(h5read(h5_path, "/matrix/data", start = ptr[[1L]] + 1,
                              count = n_nz))
  local_cell <- rep.int(seq_len(end_cell - start_cell + 1L), diff(ptr))
  umi[start_cell:end_cell] <- 0
  cell_totals <- rowsum(values, local_cell, reorder = FALSE)
  umi[start_cell - 1L + as.integer(rownames(cell_totals))] <- cell_totals[, 1]
  if (start_cell %% (chunk_cells * 10L) == 1L || start_cell == 1L) {
    message("UMI scan through reference nucleus ", end_cell, "/", n_cells)
  }
}

cell_types <- c("ASC", "EX", "INH", "MG", "ODC", "OPC", "PER.END")
if (!all(cell_types %in% unique(meta$Cell.Type))) stop("Reference cell-type ontology changed")
eligible <- which(umi >= 100 & meta$Cell.Type %in% cell_types &
                    !is.na(meta$SampleID) & nzchar(meta$SampleID))
set.seed(42L)
selected <- integer()
balance_rows <- list()
groups <- unique(meta[eligible, c("SampleID", "Cell.Type"), drop = FALSE])
groups <- groups[order(groups$SampleID, groups$Cell.Type), , drop = FALSE]
for (k in seq_len(nrow(groups))) {
  donor <- groups$SampleID[[k]]
  type <- groups$Cell.Type[[k]]
  ix <- eligible[meta$SampleID[eligible] == donor & meta$Cell.Type[eligible] == type]
  take <- min(100L, length(ix))
  chosen <- if (take == length(ix)) ix else sample(ix, take, replace = FALSE)
  selected <- c(selected, chosen)
  balance_rows[[length(balance_rows) + 1L]] <- data.frame(
    donor_id = donor, cell_type = type, available_umi_ge_100 = length(ix),
    selected = take, median_selected_umi = median(umi[chosen]),
    stringsAsFactors = FALSE
  )
}
selected <- sort(unique(selected))
if (length(selected) < 5000L) stop("Too few donor-balanced nuclei selected")
selected_col <- integer(n_cells)
selected_col[selected] <- seq_along(selected)

# Second streamed pass reads only selected columns into sparse triplets.
i_parts <- vector("list", length(chunk_starts))
j_parts <- vector("list", length(chunk_starts))
x_parts <- vector("list", length(chunk_starts))
for (q in seq_along(chunk_starts)) {
  start_cell <- chunk_starts[[q]]
  end_cell <- min(n_cells, start_cell + chunk_cells - 1L)
  ptr <- indptr[start_cell:(end_cell + 1L)]
  n_nz <- ptr[length(ptr)] - ptr[[1L]]
  if (n_nz <= 0) next
  values <- as.integer(h5read(h5_path, "/matrix/data", start = ptr[[1L]] + 1,
                              count = n_nz))
  row0 <- as.integer(h5read(h5_path, "/matrix/indices", start = ptr[[1L]] + 1,
                            count = n_nz))
  local_cell <- rep.int(seq_len(end_cell - start_cell + 1L), diff(ptr))
  global_cell <- start_cell + local_cell - 1L
  keep <- selected_col[global_cell] > 0L
  if (any(keep)) {
    i_parts[[q]] <- gene_row[row0[keep] + 1L]
    j_parts[[q]] <- selected_col[global_cell[keep]]
    x_parts[[q]] <- values[keep]
  }
  if (q %% 10L == 1L || q == length(chunk_starts)) {
    message("Sparse extraction chunk ", q, "/", length(chunk_starts))
  }
}
ii <- unlist(i_parts, use.names = FALSE)
jj <- unlist(j_parts, use.names = FALSE)
xx <- unlist(x_parts, use.names = FALSE)
selected_meta <- meta[selected, c("Barcode", "SampleID", "Cell.Type"), drop = FALSE]
counts <- sparseMatrix(i = ii, j = jj, x = xx,
                       dims = c(length(gene_levels), length(selected)),
                       dimnames = list(gene_levels, selected_meta$Barcode),
                       repr = "C")
rm(i_parts, j_parts, x_parts, ii, jj, xx)
gc()
cell_type <- factor(selected_meta$Cell.Type, levels = cell_types)
names(cell_type) <- selected_meta$Barcode
donor_id <- selected_meta$SampleID
names(donor_id) <- selected_meta$Barcode

out <- list(
  counts = counts,
  cell_type = cell_type,
  donor_id = donor_id,
  selected_metadata = selected_meta,
  selected_umi = as.integer(umi[selected]),
  source = list(
    accession = "GSE174367", platform = "human cortical single-nucleus RNA-seq",
    matrix_file_sha256 = system2("certutil", c("-hashfile", h5_path, "SHA256"),
                                 stdout = TRUE, stderr = FALSE),
    metadata_file = metadata_path, selected_per_donor_type_cap = 100L,
    min_umi = 100L, seed = 42L, diagnosis_used_for_selection = FALSE,
    matrix_barcodes_absent_from_metadata = unmapped_barcodes
  )
)
saveRDS(out, file.path(output_dir, "G6_RCTD_reference_cells.rds"), compress = "gzip")
write.csv(do.call(rbind, balance_rows),
          file.path(output_dir, "G6_RCTD_reference_balance.csv"), row.names = FALSE)
write.csv(data.frame(n_nuclei = ncol(counts), n_genes = nrow(counts),
                     n_nonzero = length(counts@x), n_donors = length(unique(donor_id)),
                     n_cell_types = nlevels(droplevels(cell_type)),
                     min_umi = min(umi[selected]), median_umi = median(umi[selected]),
                     max_umi = max(umi[selected])),
          file.path(output_dir, "G6_RCTD_reference_qc.csv"), row.names = FALSE)
message("Prepared ", ncol(counts), " reference nuclei; ", nrow(counts),
        " genes; ", length(counts@x), " nonzeros")
