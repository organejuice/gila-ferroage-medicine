#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(rhdf5)
  library(Matrix)
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop("Usage: 84_extract_g6_multiscale_rctd_inputs.R <project_root> <ascii_staged_tissue_gef> <GSM> <50|110>")
}
root <- args[[1]]
gef_path <- normalizePath(args[[2]], winslash = "/", mustWork = TRUE)
if (grepl("[^ -~]", gef_path)) stop("HDF5 input must be staged at an ASCII path")
gsm <- args[[3]]
bin_width <- as.integer(args[[4]])
if (is.na(bin_width) || !bin_width %in% c(50L, 110L)) stop("Allowed widths: 50, 110 DNB")
audit <- file.path(root, "results/2026-09-24-disease-side-v3.5-g6-input-audit")
interim <- file.path(root, "data/interim/G6_GSE269906_multiscale_001", paste0("bin", bin_width), "signature")
outdir <- file.path(root, "data/interim/G6_GSE269906_multiscale_001", paste0("bin", bin_width), "full")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

manifest <- read.delim(file.path(audit, "registries/GSE269906_tissue_CNGB_file_manifest.tsv"),
                       check.names = FALSE, stringsAsFactors = FALSE)
donors <- read.delim(file.path(audit, "registries/GSE269906_donor_section_manifest.tsv"),
                     check.names = FALSE, stringsAsFactors = FALSE)
file_row <- manifest[manifest$GSM == gsm, , drop = FALSE]
donor_row <- donors[donors$GSM == gsm, , drop = FALSE]
if (nrow(file_row) != 1L || nrow(donor_row) != 1L) stop("GSM must map to exactly one audited donor and tissue GEF")
if (file.info(gef_path)$size != file_row$GEO_expected_bytes) stop("Staged GEF byte size differs from the audited manifest")
if (!requireNamespace("digest", quietly = TRUE)) stop("The digest package is required for source verification")
staged_sha <- digest::digest(file = gef_path, algo = "sha256", serialize = FALSE)
if (!identical(tolower(staged_sha), tolower(file_row$sha256))) stop("Staged GEF SHA-256 differs from the audited source manifest")

gene_rds <- readRDS(file.path(interim, paste0(gsm, "_bin", bin_width, "_tissue_signature_counts.rds")))
eligible_meta <- gene_rds$metadata[gene_rds$metadata$n_segmented_cell_centers > 0L, , drop = FALSE]
eligible_tiles <- as.integer(eligible_meta$tile_id)
if (!length(eligible_tiles) || anyDuplicated(eligible_tiles)) stop("No unique centroid-supported bins")
if (any(eligible_meta$n_features < 400L)) stop("Selected-width metadata violate the pilot >=400 detected-gene screen")

ref <- readRDS(file.path(root, "data/interim/G6_GSE269906_multiscale_001/reference/G6_RCTD_reference_no_candidate_overlap.rds"))
ref_genes <- rownames(ref$counts)
g1 <- read.delim(gzfile(file.path(audit, "registries/G1_microglia_signed_sender_244.tsv.gz")),
                 check.names = FALSE, stringsAsFactors = FALSE)
g5 <- read.delim(gzfile(file.path(audit, "registries/G5_receiver_one_to_one_human_ortholog_registry.tsv.gz")),
                 check.names = FALSE, stringsAsFactors = FALSE)
state_genes <- unique(toupper(c(g1$gene, g5$human_one_to_one_symbol)))
state_genes <- state_genes[!is.na(state_genes) & nzchar(state_genes)]
holdout_marker_candidates <- c("AIF1", "C1QA", "C1QB", "C1QC", "CX3CR1", "LST1", "P2RY12", "TMEM119", "TREM2", "TYROBP",
                               "ALDH1L1", "GFAP", "GJA1", "GLUL", "S100B", "SLC1A2", "SLC1A3", "SOX9", "CD74", "AQP4")
holdout_markers <- setdiff(holdout_marker_candidates, state_genes)
feature_universe <- ref_genes
if (length(feature_universe) < 10000L) stop("Unexpectedly small reference feature universe")
gene_map <- unique(toupper(as.character(h5read(gef_path, "/geneExp/bin1/gene")$gene)))
output_genes <- intersect(gene_map, feature_universe)
if (length(output_genes) < 10000L) stop("Unexpectedly small common GEF/reference gene universe")

expr_path <- "/geneExp/bin1/expression"
gene_path <- "/geneExp/bin1/gene"
attrs <- h5readAttributes(gef_path, expr_path)
min_x <- as.numeric(attrs$minX)[1L]; max_x <- as.numeric(attrs$maxX)[1L]
min_y <- as.numeric(attrs$minY)[1L]; max_y <- as.numeric(attrs$maxY)[1L]
resolution_nm <- as.numeric(attrs$resolution)[1L]
if (resolution_nm != 500) stop("Expected audited 500 nm DNB pitch")
origin_x <- floor((min_x - 1) / bin_width) * bin_width + 1
origin_y <- floor((min_y - 1) / bin_width) * bin_width + 1
nx <- floor((max_x - 1) / bin_width) - floor((min_x - 1) / bin_width) + 1
ny <- floor((max_y - 1) / bin_width) - floor((min_y - 1) / bin_width) + 1
gene <- h5read(gef_path, gene_path)
gene_names <- toupper(as.character(gene$gene))
gene_offsets <- as.double(gene$offset)
gene_counts <- as.double(gene$count)
if (any(gene_offsets != c(0, head(cumsum(gene_counts), -1L)))) stop("GEF gene offsets are not contiguous")
selected_gene_ids <- which(gene_names %in% output_genes)
output_levels <- intersect(output_genes, unique(gene_names[selected_gene_ids]))
output_row <- match(gene_names, output_levels)
candidate_by_id <- which(!is.na(output_row))

accepted_column <- integer(max(eligible_tiles))
accepted_column[eligible_tiles] <- seq_along(eligible_tiles)
target_chunk_records <- 2000000L
active <- which(gene_counts > 0)
batches <- list(); current <- integer(); current_n <- 0
for (gene_id in active) {
  if (length(current) && current_n + gene_counts[gene_id] > target_chunk_records) {
    batches[[length(batches) + 1L]] <- current; current <- integer(); current_n <- 0
  }
  current <- c(current, gene_id); current_n <- current_n + gene_counts[gene_id]
}
if (length(current)) batches[[length(batches) + 1L]] <- current

ii_parts <- vector("list", length(batches)); jj_parts <- vector("list", length(batches)); xx_parts <- vector("list", length(batches))
for (b in seq_along(batches)) {
  ids <- batches[[b]]
  if (!any(ids %in% candidate_by_id)) next
  lo <- gene_offsets[min(ids)]
  hi <- gene_offsets[max(ids)] + gene_counts[max(ids)] - 1
  n <- as.integer(hi - lo + 1)
  expression <- h5read(gef_path, expr_path, start = lo + 1, count = n)
  expression_gene_id <- findInterval(seq.int(lo, hi), gene_offsets)
  use_gene <- expression_gene_id %in% candidate_by_id
  if (!any(use_gene)) next
  x <- as.numeric(expression$x[use_gene]); y <- as.numeric(expression$y[use_gene])
  count <- as.integer(expression$count[use_gene])
  ids_for_records <- expression_gene_id[use_gene]
  tile_id <- as.integer((x - origin_x) %/% bin_width + nx * ((y - origin_y) %/% bin_width)) + 1L
  col <- accepted_column[tile_id]
  keep <- !is.na(col) & col > 0L & count > 0L
  if (!any(keep)) next
  pairs <- data.table(gene_row = output_row[ids_for_records[keep]],
                      bin_col = col[keep], count = count[keep])[
                        , .(count = sum(count)), by = .(gene_row, bin_col)]
  ii_parts[[b]] <- pairs$gene_row
  jj_parts[[b]] <- pairs$bin_col
  xx_parts[[b]] <- pairs$count
  if (b %% 10L == 0L || b == length(batches)) message("Extracted eligible GEF batches ", b, "/", length(batches), " for ", gsm)
}
ii <- unlist(ii_parts, use.names = FALSE)
jj <- unlist(jj_parts, use.names = FALSE)
xx <- unlist(xx_parts, use.names = FALSE)
counts <- sparseMatrix(i = ii, j = jj, x = xx,
                       dims = c(length(output_levels), length(eligible_tiles)),
                       dimnames = list(output_levels, as.character(eligible_tiles)), repr = "C")
rm(ii_parts, jj_parts, xx_parts, ii, jj, xx, gene_rds, ref); gc()
if (ncol(counts) != nrow(eligible_meta) || length(counts@x) == 0L) stop("Constructed spatial count matrix is empty or misaligned")
if (any(Matrix::colSums(counts) > eligible_meta$total_UMI)) stop("Submatrix UMI count exceeds its bin's audited total")
eligible_meta$donor_id <- donor_row$donor_id
eligible_meta$GSM <- gsm
out <- list(counts = counts, metadata = eligible_meta,
            construction = list(
              run_id = "2026-09-26-g6-multiscale-blind-pilot",
              execution_lock_id = "G6-MS-2026-09-26-001",
              source_file = file_row$donor_file,
              source_sha256 = staged_sha,
              n_centroid_supported_bins = ncol(counts),
              n_common_genes = nrow(counts),
              n_G1_sender_genes_in_matrix = sum(rownames(counts) %in% toupper(g1$gene)),
              n_G5_ortholog_genes_in_matrix = sum(rownames(counts) %in% toupper(g5$human_one_to_one_symbol)),
              cell_identity = "bin-level RCTD-inferred weights; all GEF expression retained in the count submatrix",
              tile_size_um = bin_width * 0.5,
              min_features_per_bin = 400L,
              valid_polygon_overlap = TRUE,
              at_least_one_segmented_cell_centroid = TRUE,
              inference_features_excluded = "G1 signed sender, all G5 human one-to-one orthologs, and holdout markers not overlapping the frozen signatures",
              holdout_marker_candidates = holdout_marker_candidates,
              holdout_markers_disjoint_from_frozen_states = holdout_markers
            ))
saveRDS(out, file.path(outdir, paste0(gsm, "_bin", bin_width, "_centroid_supported_counts.rds")), compress = "gzip")
write.csv(data.frame(GSM = gsm, donor_id = donor_row$donor_id,
                     n_bins = ncol(counts), n_genes = nrow(counts),
                     n_nonzero = length(counts@x), median_features = median(eligible_meta$n_features),
                     median_umi = median(eligible_meta$total_UMI),
                     file_bytes = file.info(file.path(outdir, paste0(gsm, "_bin", bin_width, "_centroid_supported_counts.rds")))$size),
          file.path(outdir, paste0(gsm, "_bin", bin_width, "_centroid_supported_qc.csv")), row.names = FALSE)
message("Saved ", gsm, ": ", ncol(counts), " bins x ", nrow(counts), " genes; ", length(counts@x), " nonzeros")

