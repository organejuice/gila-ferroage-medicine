#!/usr/bin/env Rscript

# Effect-blind reconstruction of section tissue masks from CellBin cellBorder
# polygons, rasterized to the globally anchored 50x50 DNB (25 um) grid.

suppressPackageStartupMessages({
  library(rhdf5)
  library(sp)
})

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) args[[1L]] else "."
stage_dir <- if (length(args) >= 2L) args[[2L]] else Sys.getenv("G6_ASCII_HDF5_STAGE_DIR", unset = tempdir())
audit_dir <- file.path(root, "results/2026-09-24-disease-side-v3.5-g6-input-audit")
reg_dir <- file.path(audit_dir, "registries")
table_dir <- file.path(audit_dir, "tables")
figure_dir <- file.path(audit_dir, "figures")
mask_dir <- file.path(root, "data/interim/G6_GSE269906_bin50_masks")
donors <- read.delim(file.path(reg_dir, "GSE269906_donor_section_manifest.tsv"), check.names = FALSE, stringsAsFactors = FALSE)
qc <- read.delim(file.path(table_dir, "GSE269906_tissue_bin50_section_qc.tsv"), check.names = FALSE, stringsAsFactors = FALSE)
if (nrow(donors) != 12L || nrow(qc) != 12L || anyDuplicated(donors$GSM) || anyDuplicated(qc$GSM)) {
  stop("Expected 12 unique primary donors and 12 effect-blind tissue QC records.")
}
dir.create(stage_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(mask_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

polygon_intersects_tile <- function(px, py, left, right, top, bottom) {
  if (any(px >= left & px <= right & py >= top & py <= bottom)) return(TRUE)
  corners <- sp::point.in.polygon(c(left, right, right, left), c(top, top, bottom, bottom), px, py)
  if (any(corners > 0L)) return(TRUE)
  n <- length(px)
  i <- seq_len(n)
  j <- c(seq.int(2L, n), 1L)
  dx <- px[j] - px[i]
  dy <- py[j] - py[i]
  for (x_edge in c(left, right)) {
    crosses <- ((px[i] <= x_edge & px[j] >= x_edge) | (px[j] <= x_edge & px[i] >= x_edge)) & dx != 0
    if (any(crosses)) {
      y_cross <- py[i[crosses]] + (x_edge - px[i[crosses]]) / dx[crosses] * dy[crosses]
      if (any(y_cross >= top & y_cross <= bottom)) return(TRUE)
    }
  }
  for (y_edge in c(top, bottom)) {
    crosses <- ((py[i] <= y_edge & py[j] >= y_edge) | (py[j] <= y_edge & py[i] >= y_edge)) & dy != 0
    if (any(crosses)) {
      x_cross <- px[i[crosses]] + (y_edge - py[i[crosses]]) / dy[crosses] * dx[crosses]
      if (any(x_cross >= left & x_cross <= right)) return(TRUE)
    }
  }
  FALSE
}

summary_rows <- vector("list", nrow(donors))
pdf(file.path(figure_dir, "GSE269906_bin50_CellBin_polygon_mask_overlay.pdf"), width = 12, height = 8.5, onefile = TRUE)
par(mfrow = c(3, 4), mar = c(2.2, 2.2, 3, 0.8), oma = c(0, 0, 2, 0))

for (k in seq_len(nrow(donors))) {
  donor <- donors[k, , drop = FALSE]
  section_qc <- qc[match(donor$GSM, qc$GSM), , drop = FALSE]
  cellbin_path <- file.path(root, "data/raw/GSE269906_cellbin", donor$cellbin_file)
  stage_path <- file.path(stage_dir, basename(cellbin_path))
  if (!file.exists(cellbin_path) || file.info(cellbin_path)$size != donor$cellbin_expected_bytes) {
    stop("Missing or non-exact-size CellBin input for ", donor$GSM)
  }
  if (file.exists(stage_path)) unlink(stage_path)
  if (!file.copy(cellbin_path, stage_path, overwrite = FALSE, copy.mode = FALSE)) {
    stop("Could not stage CellBin input at an HDF5-compatible ASCII path: ", donor$GSM)
  }
  root_attrs <- h5readAttributes(stage_path, "/")
  cell <- h5read(stage_path, "/cellBin/cell")
  border <- h5read(stage_path, "/cellBin/cellBorder")
  if (!identical(as.character(root_attrs$bin_type), "CellBin") ||
      as.numeric(root_attrs$resolution)[1L] != section_qc$coordinate_resolution_nm ||
      !all(c("x", "y", "area", "cellTypeID") %in% names(cell)) ||
      length(dim(border)) != 3L || dim(border)[1L] != 2L || dim(border)[3L] != nrow(cell)) {
    stop("CellBin polygon geometry or coordinate scale is inconsistent for ", donor$GSM)
  }

  bin_width <- 50L
  origin_x <- floor((section_qc$minX - 1) / bin_width) * bin_width + 1
  origin_y <- floor((section_qc$minY - 1) / bin_width) * bin_width + 1
  nx <- floor((section_qc$maxX - 1) / bin_width) - floor((section_qc$minX - 1) / bin_width) + 1L
  ny <- floor((section_qc$maxY - 1) / bin_width) - floor((section_qc$minY - 1) / bin_width) + 1L
  n_tiles <- as.integer(nx * ny)

  border_x <- border[1L, , ]
  border_y <- border[2L, , ]
  valid_vertex <- border_x != 32767 & border_y != 32767
  n_vertices <- colSums(valid_vertex)
  valid_polygon <- n_vertices >= 3L & cell$area > 0
  if (!any(valid_polygon)) stop("No valid CellBin segmentation polygons for ", donor$GSM)

  polygon_mask <- logical(n_tiles)
  for (cell_id in which(valid_polygon)) {
    valid <- valid_vertex[, cell_id]
    px <- cell$x[cell_id] + border_x[valid, cell_id]
    py <- cell$y[cell_id] + border_y[valid, cell_id]
    x0 <- max(0L, as.integer(floor((min(px) - (origin_x - 0.5)) / bin_width)))
    x1 <- min(nx - 1L, as.integer(floor((max(px) - (origin_x - 0.5)) / bin_width)))
    y0 <- max(0L, as.integer(floor((min(py) - (origin_y - 0.5)) / bin_width)))
    y1 <- min(ny - 1L, as.integer(floor((max(py) - (origin_y - 0.5)) / bin_width)))
    if (x0 > x1 || y0 > y1) next
    for (tile_x in x0:x1) for (tile_y in y0:y1) {
      left <- origin_x - 0.5 + tile_x * bin_width
      right <- left + bin_width
      top <- origin_y - 0.5 + tile_y * bin_width
      bottom <- top + bin_width
      if (polygon_intersects_tile(px, py, left, right, top, bottom)) {
        polygon_mask[tile_x + nx * tile_y + 1L] <- TRUE
      }
    }
  }

  mask_tiles <- which(polygon_mask)
  count_path <- file.path(root, "data/interim/G6_GSE269906_bin50_tissue", paste0(donor$GSM, "_bin50_tissue_signature_counts.rds"))
  counts <- readRDS(count_path)
  accepted_tiles <- counts$metadata$tile_id
  if (!all(accepted_tiles %in% mask_tiles)) stop("A retained bin falls outside its CellBin polygon mask: ", donor$GSM)

  mask_object <- list(
    accession = "GSE269906",
    GSM = donor$GSM,
    donor_id = donor$donor_id,
    grid = list(bin_width_DNB = bin_width, bin_width_um = bin_width * section_qc$coordinate_resolution_nm / 1000,
                origin_x = origin_x, origin_y = origin_y, nx = nx, ny = ny),
    polygon_mask_tile_id = mask_tiles,
    accepted_bin_tile_id = accepted_tiles,
    segmentation = list(n_units = nrow(cell), n_valid_polygons = sum(valid_polygon),
                       n_invalid_polygons = sum(!valid_polygon), cellTypeID_all_default_zero = all(cell$cellTypeID == 0L),
                       polygons_from = "CellBin/cellBorder; valid vertex count >=3 and cell area >0"),
    image_available = FALSE,
    layer_or_ROI_available = FALSE
  )
  saveRDS(mask_object, file.path(mask_dir, paste0(donor$GSM, "_bin50_CellBin_mask.rds")), compress = "gzip")

  tile_x <- (mask_tiles - 1L) %% nx
  tile_y <- (mask_tiles - 1L) %/% nx
  accepted <- match(accepted_tiles, mask_tiles)
  accepted_x <- (accepted_tiles - 1L) %% nx
  accepted_y <- (accepted_tiles - 1L) %/% nx
  plot(origin_x + (tile_x + 0.5) * bin_width, origin_y + (tile_y + 0.5) * bin_width,
       pch = 15, cex = 0.24, col = "grey85", asp = 1, axes = FALSE,
       xlab = "", ylab = "", main = paste(donor$GSM, donor$donor_id))
  points(origin_x + (accepted_x + 0.5) * bin_width, origin_y + (accepted_y + 0.5) * bin_width,
         pch = 15, cex = 0.24, col = "#0072B2")
  box()

  summary_rows[[k]] <- data.frame(
    accession = "GSE269906", GSM = donor$GSM, donor_id = donor$donor_id,
    n_CellBin_units = nrow(cell), n_valid_polygons = sum(valid_polygon), n_invalid_polygons = sum(!valid_polygon),
    fraction_valid_polygons = mean(valid_polygon), n_polygon_mask_bin50 = length(mask_tiles),
    n_accepted_bin50_ge400_and_masked = length(accepted_tiles),
    n_accepted_bins_outside_mask = sum(!accepted_tiles %in% mask_tiles),
    fraction_CellBin_cellTypeID_default_zero = mean(cell$cellTypeID == 0L),
    coordinate_resolution_nm = section_qc$coordinate_resolution_nm, bin_width_um = bin_width * section_qc$coordinate_resolution_nm / 1000,
    histology_image_available = FALSE, layer_or_ROI_available = FALSE,
    stringsAsFactors = FALSE
  )
  unlink(stage_path)
  rm(cell, border, counts, polygon_mask)
  gc(verbose = FALSE)
  message("effect_blind_polygon_mask_complete=", donor$GSM, "; mask_tiles=", length(mask_tiles), "; accepted_bins=", length(accepted_tiles))
}
dev.off()

mask_summary <- do.call(rbind, summary_rows)
write.table(mask_summary, file.path(table_dir, "GSE269906_bin50_CellBin_polygon_mask_qc.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "")
if (sum(mask_summary$n_accepted_bins_outside_mask) != 0L) stop("Retained bin/mask consistency audit failed.")
cat("sections=", nrow(mask_summary), "\n", sep = "")
cat("valid_polygons_total=", sum(mask_summary$n_valid_polygons), "\n", sep = "")
cat("polygon_mask_tiles_total=", sum(mask_summary$n_polygon_mask_bin50), "\n", sep = "")
cat("accepted_bins_total=", sum(mask_summary$n_accepted_bin50_ge400_and_masked), "\n", sep = "")
