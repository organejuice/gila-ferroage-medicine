#!/usr/bin/env Rscript

# Effect-blind construction and QC of 50x50 Stereo-seq pseudo-spots from the
# deposited cellBin objects. This is a derived interim input, not an AD effect
# analysis and not the authors' original DNB-level bin50 matrix.

suppressPackageStartupMessages({
  library(rhdf5)
  library(Matrix)
})

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) args[[1L]] else "."
out_dir <- file.path(root, "results/2026-09-24-disease-side-v3.5-g6-input-audit")
reg_dir <- file.path(out_dir, "registries")
table_dir <- file.path(out_dir, "tables")
figure_dir <- file.path(out_dir, "figures")
interim_dir <- file.path(root, "data/interim/G6_GSE269906_bin50")
manifest <- read.delim(file.path(reg_dir, "GSE269906_file_acquisition_manifest.tsv"),
                       check.names = FALSE, stringsAsFactors = FALSE)
manifest <- manifest[manifest$file_type == "cellbin" & manifest$status == "complete", , drop = FALSE]
donors <- read.delim(file.path(reg_dir, "GSE269906_donor_section_manifest.tsv"),
                     check.names = FALSE, stringsAsFactors = FALSE)
g1 <- read.delim(gzfile(file.path(reg_dir, "G1_microglia_signed_sender_244.tsv.gz")),
                 check.names = FALSE, stringsAsFactors = FALSE)
g5 <- read.delim(gzfile(file.path(reg_dir, "G5_receiver_one_to_one_human_ortholog_registry.tsv.gz")),
                 check.names = FALSE, stringsAsFactors = FALSE)
dir.create(interim_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

bin_width <- 50L
min_features <- 400L
resolution_nm <- 500
if (nrow(manifest) != 12L || nrow(donors) != 12L) stop("Expected exactly 12 mapped GSE269906 donors/sections.")
if (length(unique(vapply(manifest$file, function(f) file.info(file.path(root, "data/raw/GSE269906_cellbin", f))$size, numeric(1)))) < 1L) {
  stop("No complete cellBin GEF files are available.")
}
if (nrow(g1) != 244L || length(unique(g1$gene)) != 244L) stop("The frozen G1 sender signature is not the audited 244-gene set.")

write_gz <- function(x, path) {
  con <- gzfile(path, open = "wt", compression = 6)
  on.exit(close(con))
  write.table(x, con, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
}

qc_rows <- list()
coverage_rows <- list()
plot_rows <- list()

pdf(file.path(figure_dir, "GSE269906_bin50_effect_blind_qc.pdf"), width = 12, height = 8.5, onefile = TRUE)
par(mfrow = c(3, 4), mar = c(2.5, 2.5, 3, 1), oma = c(0, 0, 2, 0))
for (k in seq_len(nrow(manifest))) {
  row <- manifest[k, , drop = FALSE]
  source <- file.path(root, "data/raw/GSE269906_cellbin", row$file)
  if (!file.exists(source) || file.info(source)$size != row$expected_bytes) stop("Raw file size differs from the frozen acquisition manifest: ", row$file)
  h5 <- tempfile(pattern = "g6_cellbin_qc_", fileext = ".gef", tmpdir = tempdir())
  on.exit(unlink(h5), add = TRUE)
  if (!file.copy(source, h5, overwrite = TRUE)) stop("Could not stage GEF at an ASCII temporary path: ", row$file)

  attrs <- h5readAttributes(h5, "/")
  res <- as.numeric(attrs$resolution)
  if (length(res) != 1L || res != resolution_nm) stop("Unexpected physical resolution for ", row$file, ": ", paste(res, collapse = ","))
  cells <- h5read(h5, "/cellBin/cell")
  cell_exp <- h5read(h5, "/cellBin/cellExp")
  genes <- h5read(h5, "/cellBin/gene")
  if (!all(c("x", "y", "offset", "geneCount", "expCount", "area") %in% names(cells))) stop("Required cell-bin fields absent in ", row$file)
  if (!all(c("geneID", "count") %in% names(cell_exp)) || !all(c("geneName", "offset", "cellCount") %in% names(genes))) stop("Unexpected GEF expression structure in ", row$file)

  expected_offsets <- c(0, head(cumsum(as.double(cells$geneCount)), -1L))
  if (any(cells$offset != expected_offsets) || sum(as.double(cells$geneCount)) != nrow(cell_exp)) {
    stop("Cell-bin expression offsets/counts are not contiguous and complete in ", row$file)
  }
  if (sum(as.double(cell_exp$count)) != sum(as.double(cells$expCount))) stop("Cell-bin expression sums do not reproduce cell expCount in ", row$file)
  if (anyNA(cells$x) || anyNA(cells$y) || anyNA(cells$area) || any(cells$area <= 0) ||
      anyNA(cell_exp$geneID) || any(cell_exp$geneID < 0) || any(cell_exp$geneID >= nrow(genes)) ||
      anyNA(cell_exp$count) || any(cell_exp$count < 0) || any(abs(cell_exp$count - round(cell_exp$count)) > 1e-8)) {
    stop("Invalid coordinates, segmentation area, gene IDs, or expression counts in ", row$file)
  }

  gene_symbols <- toupper(as.character(genes$geneName))
  if (anyNA(gene_symbols) || any(!nzchar(gene_symbols))) stop("Missing gene symbols in ", row$file)
  unique_symbols <- unique(gene_symbols)
  gene_row <- match(gene_symbols, unique_symbols)
  tile_x <- floor(as.numeric(cells$x) / bin_width)
  tile_y <- floor(as.numeric(cells$y) / bin_width)
  tile_key_per_cell <- paste(tile_x, tile_y, sep = ":")
  tile_keys <- unique(tile_key_per_cell)
  tile_index_per_cell <- match(tile_key_per_cell, tile_keys)
  tile_index_per_expression <- rep.int(tile_index_per_cell, as.integer(cells$geneCount))
  if (length(tile_index_per_expression) != nrow(cell_exp)) stop("Cell-to-expression indexing failed in ", row$file)
  expression_gene_row <- gene_row[as.integer(cell_exp$geneID) + 1L]
  if (anyNA(expression_gene_row)) stop("Expression gene IDs do not map to the GEF gene table in ", row$file)

  counts <- sparseMatrix(
    i = expression_gene_row,
    j = tile_index_per_expression,
    x = as.numeric(cell_exp$count),
    dims = c(length(unique_symbols), length(tile_keys)),
    dimnames = list(unique_symbols, tile_keys),
    giveCsparse = TRUE
  )
  tile_xy <- do.call(rbind, strsplit(tile_keys, ":", fixed = TRUE))
  spot_meta <- data.frame(
    tile_id = tile_keys,
    grid_x = as.integer(tile_xy[, 1L]),
    grid_y = as.integer(tile_xy[, 2L]),
    x_um = (as.integer(tile_xy[, 1L]) * bin_width + bin_width / 2) * res / 1000,
    y_um = (as.integer(tile_xy[, 2L]) * bin_width + bin_width / 2) * res / 1000,
    n_segment_units = tabulate(tile_index_per_cell, nbins = length(tile_keys)),
    n_features = as.integer(Matrix::colSums(counts > 0)),
    total_UMI = as.numeric(Matrix::colSums(counts)),
    qc_min400_features = Matrix::colSums(counts > 0) >= min_features,
    stringsAsFactors = FALSE
  )
  if (any(spot_meta$total_UMI < 0) || any(spot_meta$n_features > nrow(counts))) stop("Invalid aggregated pseudo-spot QC metrics.")
  accepted <- which(spot_meta$qc_min400_features)
  counts_qc <- counts[, accepted, drop = FALSE]
  spot_meta$donor_id <- row$donor_id
  spot_meta$GSM <- row$GSM

  saveRDS(list(counts = counts_qc, metadata = spot_meta[accepted, , drop = FALSE],
               all_spot_metadata = spot_meta, gene_symbols = unique_symbols,
               construction = list(platform = "GSE269906 Stereo-seq cellBin GEF", bin_width_pixels = bin_width,
                                   pixel_resolution_nm = res, tile_width_um = bin_width * res / 1000,
                                   minimum_detected_features = min_features,
                                   expression_source = "cellBin/cellExp assigned to cell-bin centroids",
                                   not_authors_original_bin50 = TRUE)),
          file.path(interim_dir, paste0(row$GSM, "_bin50_cellbin_aggregated.rds")), compress = "gzip")

  qc_rows[[k]] <- data.frame(
    accession = "GSE269906", donor_id = row$donor_id, GSM = row$GSM,
    n_segment_units = nrow(cells), n_all_bin50_tiles = nrow(spot_meta),
    n_tiles_ge400_features = length(accepted), fraction_tiles_ge400_features = length(accepted) / nrow(spot_meta),
    median_features_all_tiles = median(spot_meta$n_features), median_features_qc_tiles = if (length(accepted)) median(spot_meta$n_features[accepted]) else NA_real_,
    median_UMI_all_tiles = median(spot_meta$total_UMI), median_UMI_qc_tiles = if (length(accepted)) median(spot_meta$total_UMI[accepted]) else NA_real_,
    median_segment_units_per_tile = median(spot_meta$n_segment_units),
    median_segment_area_px2 = median(as.numeric(cells$area)),
    coordinate_resolution_nm = res, pseudo_spot_width_um = bin_width * res / 1000,
    all_cell_ids_unique = !anyDuplicated(cells$id), cellbin_type_annotation = paste(sort(unique(cells$cellTypeID)), collapse = ";"),
    polygon_mask_available = "/cellBin/cellBorder present; centroid-occupied 25um tiles used as conservative mask",
    histology_image_available = FALSE, anatomical_layer_ROI_available = FALSE,
    stringsAsFactors = FALSE
  )

  symbols_in_qc <- unique_symbols[Matrix::rowSums(counts_qc > 0) > 0]
  coverage_rows[[k]] <- rbind(
    data.frame(accession = "GSE269906", donor_id = row$donor_id, GSM = row$GSM,
               signature = "G1_signed_sender", n_frozen_human_genes = length(unique(g1$gene)),
               n_genes_detected_in_qc_tiles = sum(toupper(g1$gene) %in% symbols_in_qc),
               stringsAsFactors = FALSE),
    do.call(rbind, lapply(sort(unique(g5$module[g5$G5_module_gate == "PASS"])), function(module) {
      module_genes <- unique(toupper(g5$human_one_to_one_symbol[g5$module == module & g5$G5_module_gate == "PASS" & nzchar(g5$human_one_to_one_symbol)]))
      data.frame(accession = "GSE269906", donor_id = row$donor_id, GSM = row$GSM,
                 signature = module, n_frozen_human_genes = length(module_genes),
                 n_genes_detected_in_qc_tiles = sum(module_genes %in% symbols_in_qc),
                 stringsAsFactors = FALSE)
    }))
  )

  qc_colors <- hcl.colors(100, "YlOrRd", rev = FALSE)
  color_index <- pmax(1L, pmin(100L, as.integer(1 + 99 * (log1p(spot_meta$n_features) - min(log1p(spot_meta$n_features))) /
                                                        max(1e-9, diff(range(log1p(spot_meta$n_features)))))))
  plot(spot_meta$x_um, spot_meta$y_um, pch = 15, cex = 0.45, col = qc_colors[color_index],
       asp = 1, axes = FALSE, xlab = "", ylab = "", main = row$GSM)
  box()
  plot_rows[[k]] <- data.frame(GSM = row$GSM, donor_id = row$donor_id,
                               x_um = spot_meta$x_um, y_um = spot_meta$y_um,
                               n_features = spot_meta$n_features,
                               qc_min400_features = spot_meta$qc_min400_features)
  unlink(h5)
  message("effect_blind_bin50_qc_complete=", row$GSM, "; tiles_ge400=", length(accepted), "/", nrow(spot_meta))
}
dev.off()

qc <- do.call(rbind, qc_rows)
coverage <- do.call(rbind, coverage_rows)
write.table(qc, file.path(table_dir, "GSE269906_bin50_section_qc.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")
write.table(coverage, file.path(table_dir, "GSE269906_state_gene_coverage_by_section.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")
write_gz(do.call(rbind, plot_rows), file.path(table_dir, "GSE269906_bin50_spatial_qc_coordinates.tsv.gz"))
cat("complete_sections=", nrow(qc), "\n", sep = "")
cat("bin50_width_um=", bin_width * resolution_nm / 1000, "\n", sep = "")
cat("min_detected_features=", min_features, "\n", sep = "")
cat("qc_tiles_total=", sum(qc$n_tiles_ge400_features), "\n", sep = "")
cat("g1_detected_gene_coverage_range=", paste(range(coverage$n_genes_detected_in_qc_tiles[coverage$signature == "G1_signed_sender"] / coverage$n_frozen_human_genes[coverage$signature == "G1_signed_sender"]), collapse = ".."), "\n", sep = "")
