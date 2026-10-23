#!/usr/bin/env Rscript

# Effect-blind processing of the official GSE269906 tissue GEFs. Reconstructs
# bin50 or bin110 pseudo-spots from bin1 expression, applies the pilot's
# common >=400 detected-feature screen for QC, and keeps frozen G1/G5 genes for
# subsequent locked state scoring. This script never compares diagnoses.

suppressPackageStartupMessages({
  library(rhdf5)
  library(Matrix)
  library(data.table)
  library(sp)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L || length(args) > 3L) stop("Usage: 83_build_g6_multiscale_blind.R <project_root> <50|110> [ascii_stage_dir]")
root <- args[[1L]]
bin_width <- as.integer(args[[2L]])
if (is.na(bin_width) || !bin_width %in% c(50L, 110L)) stop("Allowed widths: 50, 110 DNB")
stage_dir <- if (length(args) >= 3L) args[[3L]] else Sys.getenv("G6_ASCII_HDF5_STAGE_DIR", unset = tempdir())
raw_root <- normalizePath(file.path(root, "data/raw"), winslash = "/", mustWork = TRUE)
stage_resolved <- normalizePath(stage_dir, winslash = "/", mustWork = TRUE)
if (startsWith(tolower(stage_resolved), paste0(tolower(raw_root), "/")) ||
    identical(tolower(stage_resolved), tolower(raw_root)) || grepl("[^ -~]", stage_resolved))
  stop("ASCII stage directory must exist outside read-only raw data and use ASCII characters")
audit_dir <- file.path(root, "results/2026-09-24-disease-side-v3.5-g6-input-audit")
pilot_id <- "G6-MS-2026-09-26-001"
pilot_dir <- file.path(root, "results/2026-09-26-g6-multiscale-blind-pilot", paste0("bin", bin_width))
reg_dir <- file.path(audit_dir, "registries")
table_dir <- file.path(pilot_dir, "tables")
figure_dir <- file.path(pilot_dir, "figures")
input_dir <- file.path(root, "data/raw/GSE269906_tissue")
interim_dir <- file.path(root, "data/interim/G6_GSE269906_multiscale_001", paste0("bin", bin_width), "signature")
manifest_path <- file.path(reg_dir, "GSE269906_tissue_CNGB_file_manifest.tsv")
donors_path <- file.path(reg_dir, "GSE269906_donor_section_manifest.tsv")
g1_path <- file.path(reg_dir, "G1_microglia_signed_sender_244.tsv.gz")
g5_path <- file.path(reg_dir, "G5_receiver_one_to_one_human_ortholog_registry.tsv.gz")

minimum_features <- 400L
target_chunk_records <- 2000000L
if (!file.exists(manifest_path) || !file.exists(donors_path) || !file.exists(g1_path) || !file.exists(g5_path)) {
  stop("Missing one or more audited manifests or frozen signature registries.")
}
manifest <- read.delim(manifest_path, check.names = FALSE, stringsAsFactors = FALSE)
donors <- read.delim(donors_path, check.names = FALSE, stringsAsFactors = FALSE)[, c("GSM", "donor_id")]
g1 <- read.delim(gzfile(g1_path), check.names = FALSE, stringsAsFactors = FALSE)
g5 <- read.delim(gzfile(g5_path), check.names = FALSE, stringsAsFactors = FALSE)
if (nrow(manifest) != 12L || nrow(donors) != 12L || anyDuplicated(manifest$GSM) || anyDuplicated(donors$GSM)) {
  stop("Expected 12 uniquely mapped primary donors and tissue GEF files.")
}
if (nrow(g1) != 244L || anyDuplicated(g1$gene)) stop("Frozen G1 sender set does not contain 244 unique genes.")
if (any(!manifest$status %in% c("complete_exact_length", "existing_exact_reassembled")) ||
    any(manifest$GEO_expected_bytes != manifest$CNGB_actual_bytes)) {
  stop("At least one tissue GEF is incomplete or differs from the GEO-listed byte size.")
}

dir.create(interim_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(stage_dir, recursive = TRUE, showWarnings = FALSE)

g5_modules <- sort(unique(g5$module))
state_symbols <- unique(toupper(c(g1$gene, g5$human_one_to_one_symbol)))
lineage_markers_raw <- list(
  MG_lineage_markers = c("AIF1", "C1QA", "C1QB", "C1QC", "CX3CR1", "LST1", "P2RY12", "TMEM119", "TREM2", "TYROBP"),
  ASC_lineage_markers = c("ALDH1L1", "GFAP", "GJA1", "GLUL", "S100B", "SLC1A2", "SLC1A3", "SOX9")
)
# Keep lineage-abundance proxies disjoint from the frozen sender/receiver genes.
lineage_markers <- lapply(lineage_markers_raw, setdiff, y = state_symbols)
identity_markers <- list(MG_paper_marker_CD74 = "CD74", ASC_paper_marker_AQP4 = "AQP4")
candidate_symbols <- unique(toupper(c(
  g1$gene,
  g5$human_one_to_one_symbol[nzchar(g5$human_one_to_one_symbol)],
  unlist(lineage_markers, use.names = FALSE),
  unlist(identity_markers, use.names = FALSE)
)))
candidate_symbols <- candidate_symbols[!is.na(candidate_symbols) & nzchar(candidate_symbols)]

build_batches <- function(gene_counts, target) {
  active <- which(gene_counts > 0)
  batches <- list()
  current <- integer()
  current_n <- 0
  for (gene_id in active) {
    if (length(current) && current_n + gene_counts[gene_id] > target) {
      batches[[length(batches) + 1L]] <- current
      current <- integer()
      current_n <- 0
    }
    current <- c(current, gene_id)
    current_n <- current_n + gene_counts[gene_id]
  }
  if (length(current)) batches[[length(batches) + 1L]] <- current
  batches
}

polygon_intersects_tile <- function(px, py, left, right, top, bottom) {
  if (any(px >= left & px <= right & py >= top & py <= bottom)) return(TRUE)
  corners <- sp::point.in.polygon(
    c(left, right, right, left), c(top, top, bottom, bottom), px, py
  )
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

qc_rows <- vector("list", nrow(manifest))
coverage_rows <- list()
pdf(file.path(figure_dir, paste0("GSE269906_tissue_bin", bin_width, "_effect_blind_qc.pdf")), width = 12, height = 8.5, onefile = TRUE)
on.exit(dev.off(), add = TRUE)
par(mfrow = c(3, 4), mar = c(2.5, 2.5, 3, 1), oma = c(0, 0, 2, 0))

for (k in seq_len(nrow(manifest))) {
  file_row <- manifest[k, , drop = FALSE]
  donor_row <- donors[match(file_row$GSM, donors$GSM), , drop = FALSE]
  raw_input_path <- file.path(input_dir, file_row$donor_file)
  if (!file.exists(raw_input_path) || file.info(raw_input_path)$size != file_row$GEO_expected_bytes) {
    stop("Missing tissue GEF or exact-size check failed: ", file_row$donor_file)
  }
  if (!requireNamespace("digest", quietly = TRUE)) stop("digest required")
  observed_sha <- digest::digest(file = raw_input_path, algo = "sha256", serialize = FALSE)
  if (!identical(tolower(observed_sha), tolower(file_row$sha256))) stop("Raw GEF SHA-256 mismatch: ", file_row$GSM)
  input_path <- file.path(stage_dir, file_row$donor_file)
  if (file.exists(input_path)) unlink(input_path)
  if (!file.copy(raw_input_path, input_path, overwrite = FALSE, copy.mode = FALSE) ||
      file.info(input_path)$size != file_row$GEO_expected_bytes) {
    stop("Could not stage tissue GEF at an HDF5-compatible ASCII path: ", file_row$donor_file)
  }

  root_attrs <- h5readAttributes(input_path, "/")
  expression_path <- "/geneExp/bin1/expression"
  gene_path <- "/geneExp/bin1/gene"
  paths <- h5ls(input_path, recursive = TRUE)$name
  if (!all(c("expression", "gene") %in% paths)) stop("Tissue GEF is missing bin1 expression/gene datasets: ", file_row$donor_file)
  exp_attrs <- h5readAttributes(input_path, expression_path)
  min_x <- as.numeric(exp_attrs$minX)[1L]
  max_x <- as.numeric(exp_attrs$maxX)[1L]
  min_y <- as.numeric(exp_attrs$minY)[1L]
  max_y <- as.numeric(exp_attrs$maxY)[1L]
  resolution_nm <- as.numeric(exp_attrs$resolution)[1L]
  max_exp <- as.numeric(exp_attrs$maxExp)[1L]
  if (!identical(as.character(root_attrs$omics), "Transcriptomics") ||
      any(!is.finite(c(min_x, max_x, min_y, max_y, resolution_nm, max_exp))) ||
      resolution_nm != 500 || max_x < min_x || max_y < min_y) {
    stop("Unexpected root/dataset attributes or physical scale: ", file_row$donor_file)
  }

  genes <- h5read(input_path, gene_path)
  gene_names <- toupper(as.character(genes$gene))
  gene_offsets <- as.double(genes$offset)
  gene_counts <- as.double(genes$count)
  expected_offsets <- c(0, head(cumsum(gene_counts), -1L))
  if (!all(c("gene", "offset", "count") %in% names(genes)) || anyNA(gene_names) || any(!nzchar(gene_names)) ||
      any(gene_offsets != expected_offsets) || any(gene_counts < 0) || sum(gene_counts) <= 0) {
    stop("Invalid gene catalog or non-contiguous offsets: ", file_row$donor_file)
  }

  # Anchor selected-width DNB tiles to the global 1-based DNB coordinate grid.
  origin_x <- floor((min_x - 1) / bin_width) * bin_width + 1
  origin_y <- floor((min_y - 1) / bin_width) * bin_width + 1
  nx <- floor((max_x - 1) / bin_width) - floor((min_x - 1) / bin_width) + 1
  ny <- floor((max_y - 1) / bin_width) - floor((min_y - 1) / bin_width) + 1
  n_tiles <- as.integer(nx * ny)
  n_features <- integer(n_tiles)
  total_umi <- numeric(n_tiles)
  selected_gene_ids <- which(gene_names %in% candidate_symbols)
  selected_row <- match(gene_names[selected_gene_ids], candidate_symbols)
  batches <- build_batches(gene_counts, target_chunk_records)
  sig_i <- vector("list", length(batches))
  sig_j <- vector("list", length(batches))
  sig_x <- vector("list", length(batches))
  n_records <- sum(gene_counts)
  t0 <- proc.time()[3L]

  for (b in seq_along(batches)) {
    ids <- batches[[b]]
    lo <- gene_offsets[min(ids)]
    hi <- gene_offsets[max(ids)] + gene_counts[max(ids)] - 1
    n <- as.integer(hi - lo + 1)
    expression <- h5read(input_path, expression_path, start = lo + 1, count = n)
    x <- as.numeric(expression$x)
    y <- as.numeric(expression$y)
    count <- as.integer(expression$count)
    if (length(count) != n || anyNA(x) || anyNA(y) || anyNA(count) ||
        any(x < min_x | x > max_x | y < min_y | y > max_y) || any(count <= 0) || any(count > max_exp)) {
      stop("Invalid coordinates or expression counts in ", file_row$donor_file, " batch ", b)
    }
    expression_gene_id <- findInterval(seq.int(lo, hi), gene_offsets)
    tile_id <- as.integer((x - origin_x) %/% bin_width + nx * ((y - origin_y) %/% bin_width)) + 1L
    pairs <- data.table(gene_id = expression_gene_id, tile_id = tile_id, umi = count)[,
      .(umi = sum(umi)), by = .(gene_id, tile_id)]
    n_features <- n_features + tabulate(pairs$tile_id, nbins = n_tiles)
    tile_umi <- pairs[, .(umi = sum(umi)), by = tile_id]
    total_umi[tile_umi$tile_id] <- total_umi[tile_umi$tile_id] + tile_umi$umi

    keep <- pairs$gene_id %in% selected_gene_ids
    if (any(keep)) {
      kept <- pairs[keep]
      gene_row <- match(gene_names[kept$gene_id], candidate_symbols)
      sig_i[[b]] <- gene_row
      sig_j[[b]] <- kept$tile_id
      sig_x[[b]] <- kept$umi
    }
    if (b %% 10L == 0L) message("effect_blind_bin_batch=", b, "/", length(batches), "; GSM=", file_row$GSM)
  }
  if (sum(gene_counts) != n_records || sum(total_umi) <= 0) stop("Expression record or UMI total check failed: ", file_row$donor_file)

  accepted_gene_qc <- which(n_features >= minimum_features)
  occupied <- which(n_features > 0)

  # CellBin polygons are the available section mask; actual polygon/tile
  # intersections provide cell-associated support. cellTypeID is not used.
  raw_cellbin_path <- file.path(root, "data/raw/GSE269906_cellbin", sub("tissue\\.gef$", "cellbin.gef", file_row$donor_file))
  if (!file.exists(raw_cellbin_path)) stop("Missing paired CellBin GEF: ", raw_cellbin_path)
  cellbin_path <- file.path(stage_dir, basename(raw_cellbin_path))
  if (file.exists(cellbin_path)) unlink(cellbin_path)
  if (!file.copy(raw_cellbin_path, cellbin_path, overwrite = FALSE, copy.mode = FALSE)) {
    stop("Could not stage paired CellBin GEF at an HDF5-compatible ASCII path: ", basename(raw_cellbin_path))
  }
  cell_table <- h5read(cellbin_path, "/cellBin/cell")
  borders <- h5read(cellbin_path, "/cellBin/cellBorder")
  if (!identical(as.character(h5readAttributes(cellbin_path, "/")$omics), "Transcriptomics") ||
      as.numeric(h5readAttributes(cellbin_path, "/")$resolution)[1L] != resolution_nm ||
      !all(c("x", "y", "cellTypeID") %in% names(cell_table)) || length(dim(borders)) != 3L || dim(borders)[1L] != 2L) {
    stop("Paired CellBin geometry or coordinate scale does not match tissue GEF: ", basename(cellbin_path))
  }
  valid_vertex <- borders[1L, , ] != 32767 & borders[2L, , ] != 32767
  n_valid_vertex <- colSums(valid_vertex)
  valid_polygon <- n_valid_vertex >= 3L & cell_table$area > 0
  if (length(n_valid_vertex) != nrow(cell_table) || !any(valid_polygon)) {
    stop("CellBin contains no usable segmentation polygons: ", basename(cellbin_path))
  }
  cell_tile_x <- as.integer((cell_table$x - origin_x) %/% bin_width)
  cell_tile_y <- as.integer((cell_table$y - origin_y) %/% bin_width)
  center_in_bounds <- cell_tile_x >= 0L & cell_tile_x < nx & cell_tile_y >= 0L & cell_tile_y < ny
  center_tile <- cell_tile_x[center_in_bounds] + nx * cell_tile_y[center_in_bounds] + 1L
  segmented_cell_count <- tabulate(center_tile, nbins = n_tiles)
  polygon_supported <- logical(n_tiles)
  border_x <- borders[1L, , ]
  border_y <- borders[2L, , ]
  for (cell_id in which(valid_polygon)) {
    valid <- valid_vertex[, cell_id]
    px <- cell_table$x[cell_id] + border_x[valid, cell_id]
    py <- cell_table$y[cell_id] + border_y[valid, cell_id]
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
        tile_id <- tile_x + nx * tile_y + 1L
        polygon_supported[tile_id] <- TRUE
      }
    }
  }
  accepted_centroid <- which(n_features >= minimum_features & segmented_cell_count > 0L)
  accepted <- which(n_features >= minimum_features & polygon_supported)
  if (!length(accepted)) stop("No bins pass the exploratory >=400 feature screen and verified CellBin support: ", file_row$donor_file)
  sig_i <- unlist(sig_i, use.names = FALSE)
  sig_j <- unlist(sig_j, use.names = FALSE)
  sig_x <- unlist(sig_x, use.names = FALSE)
  accepted_column <- match(sig_j, accepted)
  use_sig <- !is.na(accepted_column)
  sig_counts <- sparseMatrix(
    i = sig_i[use_sig], j = accepted_column[use_sig], x = sig_x[use_sig],
    dims = c(length(candidate_symbols), length(accepted)),
    dimnames = list(candidate_symbols, as.character(accepted)), giveCsparse = TRUE
  )
  tile_x <- (accepted - 1L) %% nx
  tile_y <- (accepted - 1L) %/% nx
  metadata <- data.frame(
    tile_id = accepted,
    grid_x = tile_x,
    grid_y = tile_y,
    x_um = (origin_x + tile_x * bin_width + (bin_width - 1) / 2) * resolution_nm / 1000,
    y_um = (origin_y + tile_y * bin_width + (bin_width - 1) / 2) * resolution_nm / 1000,
    n_features = n_features[accepted],
    total_UMI = total_umi[accepted],
    n_segmented_cell_centers = segmented_cell_count[accepted],
    donor_id = donor_row$donor_id,
    GSM = file_row$GSM,
    stringsAsFactors = FALSE
  )

  saveRDS(list(
    counts = sig_counts,
    metadata = metadata,
    candidate_gene_symbols = candidate_symbols,
    construction = list(
      source_accession = "GSE269906",
      source_file = file_row$donor_file,
      source_sha256 = file_row$sha256,
      expression_layer = expression_path,
      platform_unit = "DNB/bin1",
      bin_width_DNB = bin_width,
      bin_width_um = bin_width * resolution_nm / 1000,
      coordinate_resolution_nm = resolution_nm,
      minimum_detected_features = minimum_features,
      global_DNB_grid_origin = c(x = 1L, y = 1L),
      cell_identity = "mixed bin; retained only if a verified CellBin segmentation polygon intersects its selected-width DNB tile; CellBin cellTypeID is not used as an annotation",
      mask_definition = "cell-associated bin support from exact valid CellBin polygon/tile intersection; occupied tissue-GEF bin mask additionally used for spatial QC",
      author_scale_reconstruction = "selected-width DNB grid reconstructed from tissue GEF bin1; original matrix not deposited"
    )
  ), file.path(interim_dir, paste0(file_row$GSM, "_bin", bin_width, "_tissue_signature_counts.rds")), compress = "gzip")

  qc_rows[[k]] <- data.frame(
    accession = "GSE269906",
    donor_id = donor_row$donor_id,
    GSM = file_row$GSM,

    source_file = file_row$donor_file,
    source_sha256 = file_row$sha256,
    source_bytes = file.info(input_path)$size,
    n_gene_catalog_rows = nrow(genes),
    n_bin1_expression_records = n_records,
    total_UMI = sum(total_umi),
    n_occupied_bin = length(occupied),
    n_bin_ge400_features_unmasked = length(accepted_gene_qc),
    n_bin_ge400_and_cellbin_supported = length(accepted),
    n_bin_ge400_and_centroid_supported = length(accepted_centroid),
    fraction_occupied_bin_ge400 = length(accepted_gene_qc) / length(occupied),
    fraction_ge400_with_CellBin_polygon_support = length(accepted) / length(accepted_gene_qc),
    fraction_ge400_with_cell_centroid_support = length(accepted_centroid) / length(accepted_gene_qc),
    median_features_occupied = median(n_features[occupied]),
    median_features_ge400 = median(n_features[accepted]),
    median_UMI_ge400 = median(total_umi[accepted]),
    n_cellbin_segmentation_units = nrow(cell_table),
    n_CellBin_polygon_mask_tiles = sum(polygon_supported),
    n_invalid_CellBin_polygons = sum(!valid_polygon),
    fraction_valid_CellBin_polygons = mean(valid_polygon),
    fraction_CellBin_cellTypeID_default_zero = mean(cell_table$cellTypeID == 0L),
    fraction_CellBin_centers_in_tissue_GEF_bounds = mean(center_in_bounds),
    coordinate_resolution_nm = resolution_nm,
    pseudo_spot_width_um = bin_width * resolution_nm / 1000,
    minX = min_x, maxX = max_x, minY = min_y, maxY = max_y,
    histology_image_available = FALSE,
    mask_source = "paired CellBin/cellBorder segmentation; support requires a valid polygon/tile intersection",
    layer_or_ROI_available = FALSE,
    elapsed_seconds = round(proc.time()[3L] - t0, 1),
    stringsAsFactors = FALSE
  )

  coverage_rows[[k]] <- do.call(rbind, lapply(c("G1_signed_sender", g5_modules, names(lineage_markers), names(identity_markers)), function(signature) {
    signature_genes <- if (signature == "G1_signed_sender") {
      unique(toupper(g1$gene))
    } else if (signature %in% names(lineage_markers)) {
      unique(toupper(lineage_markers[[signature]]))
    } else if (signature %in% names(identity_markers)) {
      identity_markers[[signature]]
    } else {
      unique(toupper(g5$human_one_to_one_symbol[g5$module == signature & nzchar(g5$human_one_to_one_symbol)]))
    }
    signature_genes <- signature_genes[!is.na(signature_genes) & nzchar(signature_genes)]
    gene_rows <- match(signature_genes, candidate_symbols)
    usable_rows <- gene_rows[!is.na(gene_rows)]
    prevalence <- if (length(usable_rows)) Matrix::rowMeans(sig_counts[usable_rows, , drop = FALSE] > 0) else numeric()
    data.frame(
      accession = "GSE269906", donor_id = donor_row$donor_id, GSM = file_row$GSM,
      signature = signature, n_frozen_human_genes = length(signature_genes),
      n_in_GEF_gene_catalog = length(usable_rows),
      n_detected_in_any_ge400_bin = sum(prevalence > 0),
      n_detected_in_at_least_5pct_ge400_bins = sum(prevalence >= 0.05),
      median_bin_prevalence = if (length(prevalence)) median(prevalence) else NA_real_,
      stringsAsFactors = FALSE
    )
  }))

  color_index <- pmax(1L, pmin(100L, as.integer(1 + 99 *
    (log1p(n_features[occupied]) - min(log1p(n_features[occupied]))) /
    max(1e-9, diff(range(log1p(n_features[occupied])))))))
  plot(
    origin_x + ((occupied - 1L) %% nx + 0.5) * bin_width,
    origin_y + ((occupied - 1L) %/% nx + 0.5) * bin_width,
    pch = 15, cex = 0.38, col = hcl.colors(100, "YlOrRd")[color_index],
    asp = 1, axes = FALSE, xlab = "", ylab = "", main = file_row$GSM
  )
  box()
  message("effect_blind_bin_qc_complete=", file_row$GSM,
          "; ge400=", length(accepted), "/", length(occupied),
          "; median_features=", median(n_features[accepted]))
  unlink(c(input_path, cellbin_path))
  rm(genes, sig_counts, n_features, total_umi, expression, pairs, cell_table, borders)
  gc(verbose = FALSE)
}
dev.off()

qc <- do.call(rbind, qc_rows)
coverage <- do.call(rbind, coverage_rows)
write.table(qc, file.path(table_dir, paste0("GSE269906_tissue_bin", bin_width, "_section_qc.tsv")), sep = "\t", quote = FALSE, row.names = FALSE, na = "")
write.table(coverage, file.path(table_dir, paste0("GSE269906_tissue_bin", bin_width, "_state_gene_coverage_by_section.tsv")), sep = "\t", quote = FALSE, row.names = FALSE, na = "")
cat("complete_sections=", nrow(qc), "\n", sep = "")
cat("bin_width_um=", bin_width * 500 / 1000, "\n", sep = "")
cat("min_detected_features=", minimum_features, "\n", sep = "")
cat("qc_tiles_total_unmasked=", sum(qc$n_bin_ge400_features_unmasked), "\n", sep = "")
cat("qc_tiles_total_CellBin_supported=", sum(qc$n_bin_ge400_and_cellbin_supported), "\n", sep = "")
cat("g1_catalog_coverage_range=", paste(range(coverage$n_in_GEF_gene_catalog[coverage$signature == "G1_signed_sender"] / coverage$n_frozen_human_genes[coverage$signature == "G1_signed_sender"]), collapse = ".."), "\n", sep = "")



