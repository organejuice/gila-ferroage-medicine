#!/usr/bin/env Rscript

# Effect-blind follow-up audit for the public Chhatbar Visium HD object.
# Reports per-donor measurement/QC and geometry linkage; never compares diagnoses.
suppressPackageStartupMessages({ library(Seurat); library(Matrix) })
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop("Usage: Rscript src/107_audit_g6_chhatbar_visium_hd_ba17_qc.R VISIUM_RDS G1_TSV_GZ G5_TSV_GZ OUT_DIR")
}
input <- args[[1L]]
g1_path <- args[[2L]]
g5_path <- args[[3L]]
out_dir <- args[[4L]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

object <- readRDS(input)
if (!inherits(object, "Seurat")) stop("Expected a Seurat object")
metadata <- object@meta.data
needed <- c("Subject_ID", "Sample_id", "orig.ident", "Region",
            "Tissue_Region", "Clinical_Brain_Diagnosis_Short", "Final_celltype",
            "nCount_RNA", "nFeature_RNA", "percent.mt")
if (!all(needed %in% names(metadata))) stop("Required metadata fields are missing")
selected <- metadata$Region == "Occipital_cortex" &
  metadata$Clinical_Brain_Diagnosis_Short %in% c("AD", "Control")
metadata <- metadata[selected, , drop = FALSE]
if (length(unique(metadata$Subject_ID)) != 10L) stop("Unexpected BA17 sample-label count")
counts <- object@assays$RNA@counts
if (nrow(counts) != 18085L || ncol(counts) != 96409L) stop("Unexpected count matrix")
if (!all(rownames(metadata) %in% colnames(counts))) stop("Metadata/count barcodes do not match")

write_tsv <- function(x, name) {
  write.table(x, file.path(out_dir, name), sep = "\t", row.names = FALSE,
              quote = FALSE, na = "NA")
}
safe_q <- function(x, prob) {
  if (!length(x)) return(rep(NA_real_, length(prob)))
  as.numeric(stats::quantile(x, probs = prob, names = FALSE, type = 7))
}

# Sample-level QC distributions are descriptive; no generic single-cell cutoff is applied.
donors <- sort(unique(as.character(metadata$Subject_ID)))
qc_rows <- list()
type_rows <- list()
for (donor in donors) {
  md <- metadata[as.character(metadata$Subject_ID) == donor, , drop = FALSE]
  diagnosis <- unique(as.character(md$Clinical_Brain_Diagnosis_Short))
  tma <- unique(as.character(md$orig.ident))
  if (length(diagnosis) != 1L || length(tma) != 1L) stop("Donor label maps to multiple diagnosis/TMA values")
  classes <- c("All_nuclei", "Microglia", "Astrocytes")
  for (cell_class in classes) {
    z <- if (cell_class == "All_nuclei") md else md[as.character(md$Final_celltype) == cell_class, , drop = FALSE]
    qc_rows[[length(qc_rows) + 1L]] <- data.frame(
      sample_id = donor, diagnosis_label = diagnosis, TMA = tma,
      final_celltype = cell_class, n_nuclei = nrow(z),
      nCount_median = safe_q(z$nCount_RNA, .5),
      nCount_q25 = safe_q(z$nCount_RNA, .25), nCount_q75 = safe_q(z$nCount_RNA, .75),
      nFeature_median = safe_q(z$nFeature_RNA, .5),
      nFeature_q25 = safe_q(z$nFeature_RNA, .25), nFeature_q75 = safe_q(z$nFeature_RNA, .75),
      percent_mt_median = safe_q(z$percent.mt, .5),
      percent_mt_q25 = safe_q(z$percent.mt, .25), percent_mt_q75 = safe_q(z$percent.mt, .75),
      n_below_author_5_feature_cutoff = sum(z$nFeature_RNA < 5L),
      n_above_author_75pct_mt_cutoff = sum(z$percent.mt > 75),
      stringsAsFactors = FALSE)
  }
  tab <- sort(table(as.character(md$Final_celltype)), decreasing = TRUE)
  type_rows[[length(type_rows) + 1L]] <- data.frame(
    sample_id = donor, diagnosis_label = diagnosis, TMA = tma,
    final_celltype = names(tab), n_nuclei = as.integer(tab),
    stringsAsFactors = FALSE)
}
write_tsv(do.call(rbind, qc_rows), "visium_hd_ba17_donor_qc.tsv")
write_tsv(do.call(rbind, type_rows), "visium_hd_ba17_celltype_counts.tsv")

# Frozen-state detection at a descriptive 1% per-donor level. This is not a pass line
# for this 10-donor candidate; the historical 10/12 threshold is not transferred.
g1 <- unique(read.delim(gzfile(g1_path), check.names = FALSE)$gene)
g5 <- read.delim(gzfile(g5_path), check.names = FALSE, stringsAsFactors = FALSE)
accepted <- startsWith(g5$one_to_one_mapping_status, "accepted") &
  !is.na(g5$human_one_to_one_symbol) & g5$human_one_to_one_symbol != ""
g5 <- g5[accepted, , drop = FALSE]
modules <- c(list(G1_Microglia = unique(g1)),
             lapply(split(g5$human_one_to_one_symbol, g5$module), unique))
gene_upper <- toupper(rownames(counts))
detection_rows <- list()
for (donor in donors) {
  md <- metadata[as.character(metadata$Subject_ID) == donor, , drop = FALSE]
  for (module in names(modules)) {
    target_type <- if (module == "G1_Microglia") "Microglia" else "Astrocytes"
    cells <- rownames(md)[as.character(md$Final_celltype) == target_type]
    genes <- unique(modules[[module]])
    idx <- match(toupper(genes), gene_upper)
    present <- which(!is.na(idx))
    rate <- if (length(cells) && length(present)) {
      Matrix::rowMeans(counts[idx[present], cells, drop = FALSE] > 0)
    } else numeric()
    detection_rows[[length(detection_rows) + 1L]] <- data.frame(
      sample_id = donor,
      diagnosis_label = as.character(md$Clinical_Brain_Diagnosis_Short[1L]),
      TMA = as.character(md$orig.ident[1L]), target_celltype = target_type,
      module = module, target_cells = length(cells), frozen_genes = length(genes),
      object_features = length(present), genes_at_or_above_1pct = sum(rate >= .01),
      max_gene_detection_fraction = if (length(rate)) max(rate) else NA_real_,
      median_gene_detection_fraction = if (length(rate)) median(rate) else NA_real_,
      stringsAsFactors = FALSE)
  }
}
write_tsv(do.call(rbind, detection_rows), "visium_hd_ba17_detection_rate.tsv")

# Independent identity markers from the project's pre-existing blind SOP. State-gene
# overlaps are removed; this reports label/marker concordance without relabeling nuclei.
state_genes <- toupper(unique(c(g1, g5$human_one_to_one_symbol)))
marker_sets <- list(
  Microglia = setdiff(c("AIF1", "C1QA", "C1QB", "C1QC", "CX3CR1", "LST1",
                        "P2RY12", "TMEM119", "TREM2", "TYROBP"), state_genes),
  Astrocytes = setdiff(c("ALDH1L1", "GFAP", "GJA1", "S100B", "SLC1A2", "SLC1A3", "SOX9"), state_genes)
)
marker_rows <- list()
marker_summary_rows <- list()
for (donor in donors) {
  md <- metadata[as.character(metadata$Subject_ID) == donor, , drop = FALSE]
  all_cells <- rownames(md)
  for (target_type in names(marker_sets)) {
    target_cells <- rownames(md)[as.character(md$Final_celltype) == target_type]
    other_cells <- setdiff(all_cells, target_cells)
    genes <- marker_sets[[target_type]]
    ix <- match(genes, gene_upper)
    detected <- which(!is.na(ix))
    if (length(detected)) {
      target_rate <- if (length(target_cells)) Matrix::rowMeans(counts[ix[detected], target_cells, drop = FALSE] > 0) else rep(NA_real_, length(detected))
      other_rate <- if (length(other_cells)) Matrix::rowMeans(counts[ix[detected], other_cells, drop = FALSE] > 0) else rep(NA_real_, length(detected))
      target_score <- if (length(target_cells)) Matrix::colMeans(counts[ix[detected], target_cells, drop = FALSE] > 0) else numeric()
      other_score <- if (length(other_cells)) Matrix::colMeans(counts[ix[detected], other_cells, drop = FALSE] > 0) else numeric()
    } else {
      target_rate <- other_rate <- numeric()
      target_score <- other_score <- numeric()
    }
    for (j in seq_along(detected)) {
      marker_rows[[length(marker_rows) + 1L]] <- data.frame(
        sample_id = donor, diagnosis_label = as.character(md$Clinical_Brain_Diagnosis_Short[1L]),
        TMA = as.character(md$orig.ident[1L]), target_celltype = target_type,
        marker = genes[detected[j]], target_cells = length(target_cells), other_cells = length(other_cells),
        target_detection_fraction = target_rate[j], other_detection_fraction = other_rate[j],
        stringsAsFactors = FALSE)
    }
    target_median <- if (length(target_score)) stats::median(target_score) else NA_real_
    other_median <- if (length(other_score)) stats::median(other_score) else NA_real_
    marker_summary_rows[[length(marker_summary_rows) + 1L]] <- data.frame(
      sample_id = donor, diagnosis_label = as.character(md$Clinical_Brain_Diagnosis_Short[1L]),
      TMA = as.character(md$orig.ident[1L]), target_celltype = target_type,
      markers_after_state_overlap_removal = length(genes), markers_in_object = length(detected),
      target_cells = length(target_cells), other_cells = length(other_cells),
      target_median_fraction_markers_detected = target_median,
      other_median_fraction_markers_detected = other_median,
      target_fraction_with_any_marker = if (length(target_score)) mean(target_score > 0) else NA_real_,
      other_fraction_with_any_marker = if (length(other_score)) mean(other_score > 0) else NA_real_,
      target_fraction_with_2plus_markers = if (length(target_score)) mean(target_score >= 2 / max(1, length(detected))) else NA_real_,
      other_fraction_with_2plus_markers = if (length(other_score)) mean(other_score >= 2 / max(1, length(detected))) else NA_real_,
      target_marker_score_higher = if (is.na(target_median) || is.na(other_median)) NA else target_median > other_median,
      stringsAsFactors = FALSE)
  }
}
write_tsv(if (length(marker_rows)) do.call(rbind, marker_rows) else data.frame(),
          "visium_hd_ba17_identity_marker_detection.tsv")
write_tsv(do.call(rbind, marker_summary_rows), "visium_hd_ba17_identity_marker_summary.tsv")

# Verify coordinate/polygon cell linkage within the public Seurat object. These are
# nuclear geometries; they are not a tissue mask and coordinate units are not inferred.
fov_rows <- list()
fov_cell_rows <- list()
for (fov_name in names(object@images)) {
  fov <- object@images[[fov_name]]
  boundaries <- slot(fov, "boundaries")
  centroid_ids <- character()
  centroid_n <- 0L
  centroid_coords <- matrix(numeric(), ncol = 2L)
  if ("centroids" %in% names(boundaries)) {
    centroids <- boundaries[["centroids"]]
    centroid_coords <- slot(centroids, "coords")
    centroid_n <- nrow(centroid_coords)
    centroid_ids <- as.character(slot(centroids, "cells"))
    if (length(centroid_ids) != centroid_n) {
      if (length(rownames(centroid_coords)) == centroid_n) centroid_ids <- rownames(centroid_coords)
      else stop(paste("Centroid IDs do not match coordinate rows in", fov_name))
    }
    md_idx <- match(centroid_ids, rownames(metadata))
    keep <- !is.na(md_idx)
    if (any(keep)) fov_cell_rows[[length(fov_cell_rows) + 1L]] <- data.frame(
      fov = fov_name, cell_id = centroid_ids[keep],
      sample_id = as.character(metadata$Subject_ID[md_idx[keep]]),
      TMA = as.character(metadata$orig.ident[md_idx[keep]]),
      x = as.numeric(centroid_coords[keep, 1L]), y = as.numeric(centroid_coords[keep, 2L]),
      stringsAsFactors = FALSE)
  }
  seg_n <- 0L
  seg_ids <- character()
  if ("segmentation" %in% names(boundaries)) {
    seg <- boundaries[["segmentation"]]
    # SeuratObject stores this Segmentation as named polygon geometry; the
    # advertised sf.data slot is not present in the serialized class instance.
    seg_ids <- names(seg)
    if (!length(seg_ids) && "polygons" %in% slotNames(seg)) {
      polygons <- slot(seg, "polygons")
      seg_ids <- names(polygons)
      if (!length(seg_ids)) seg_ids <- vapply(polygons, function(p) as.character(p@ID), character(1))
    }
    seg_n <- length(seg_ids)
  }
  centroid_match <- match(centroid_ids, rownames(metadata))
  seg_match <- match(seg_ids, rownames(metadata))
  selected_centroid <- !is.na(centroid_match)
  selected_seg <- !is.na(seg_match)
  fov_rows[[length(fov_rows) + 1L]] <- data.frame(
    fov = fov_name, centroid_rows = centroid_n, centroid_cell_ids = length(centroid_ids),
    centroid_ids_matching_ba17_metadata = sum(selected_centroid),
    centroid_ba17_donor_labels = length(unique(as.character(metadata$Subject_ID[centroid_match[selected_centroid]]))),
    segmentation_rows = seg_n, segmentation_cell_ids = length(seg_ids),
    segmentation_ids_matching_ba17_metadata = sum(selected_seg),
    segmentation_ba17_donor_labels = length(unique(as.character(metadata$Subject_ID[seg_match[selected_seg]]))),
    centroid_segmentation_id_intersection = length(intersect(centroid_ids, seg_ids)),
    has_raster_image = any(c("image", "raster", "data") %in% slotNames(fov)),
    has_tissue_mask = FALSE, coordinate_unit_confirmed = FALSE,
    stringsAsFactors = FALSE)
}
write_tsv(do.call(rbind, fov_rows), "visium_hd_fov_geometry_linkage.tsv")
if (length(fov_cell_rows)) {
  fov_cells <- do.call(rbind, fov_cell_rows)
  write_tsv(fov_cells, "visium_hd_fov_ba17_cell_linkage.tsv")
  repeated <- names(which(table(fov_cells$cell_id) > 1L))
  duplicate_rows <- list()
  if (length(repeated)) for (cell_id in repeated) {
    z <- fov_cells[fov_cells$cell_id == cell_id, , drop = FALSE]
    duplicate_rows[[length(duplicate_rows) + 1L]] <- data.frame(
      cell_id = cell_id, sample_id = z$sample_id[1L], TMA = z$TMA[1L],
      fov_count = length(unique(z$fov)), fovs = paste(sort(unique(z$fov)), collapse = ";"),
      x_min = min(z$x), x_max = max(z$x), y_min = min(z$y), y_max = max(z$y),
      coordinate_range_x = diff(range(z$x)), coordinate_range_y = diff(range(z$y)),
      stringsAsFactors = FALSE)
  }
  if (length(duplicate_rows)) {
    dup <- do.call(rbind, duplicate_rows)
    write_tsv(dup, "visium_hd_fov_duplicate_cell_geometry.tsv")
    summary <- do.call(rbind, lapply(split(dup, dup$sample_id), function(z) data.frame(
      sample_id = z$sample_id[1L], TMA = z$TMA[1L], duplicate_cell_ids = nrow(z),
      repeated_with_same_coordinates = sum(z$coordinate_range_x < 1e-8 & z$coordinate_range_y < 1e-8),
      repeated_with_different_coordinates = sum(z$coordinate_range_x >= 1e-8 | z$coordinate_range_y >= 1e-8),
      fov_pairs = paste(sort(unique(z$fovs)), collapse = ","), stringsAsFactors = FALSE)))
    write_tsv(summary, "visium_hd_fov_duplicate_summary.tsv")
  } else {
    write_tsv(data.frame(sample_id = character(), TMA = character(), duplicate_cell_ids = integer(),
                         repeated_with_same_coordinates = integer(), repeated_with_different_coordinates = integer(),
                         fov_pairs = character()), "visium_hd_fov_duplicate_summary.tsv")
  }
}
message("Visium HD BA17 audit complete; measurements are donor-level and effect-blind.")
