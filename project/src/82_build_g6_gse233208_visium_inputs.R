#!/usr/bin/env Rscript

# Build donor-scoped Visium inputs only after the independent-reference gate.
# Diagnosis/group labels are intentionally removed from the saved inputs.
suppressPackageStartupMessages(library(Matrix))
args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1L) args[[1L]] else "."
if (!identical(root, ".")) tryCatch(setwd(root), error = function(e) stop("Could not switch to the project root: ", root))
root <- "."
protocol <- file.path(root, "docs/protocol/disease_side_v3.6_G6_GSE233208.md")
if (!file.exists(protocol)) stop("Project protocol file is missing under root=", root, ": ", protocol)

audit_dir <- file.path(root, "results/2026-09-24-disease-side-v3.5-g6-input-audit/tables")
preflight <- file.path(root, "results/2026-09-26-g6-gse233208-visium-input-audit/GSE233208_effect_blind_input_preflight.json")
reference_lock <- file.path(root, "configs/g6_v3.6_reference_admissibility.json")
object_path <- file.path(root, "data/interim/GSE233208/GSE233208_Human_visium_ADDS_seurat_processed.rds")
manifest_path <- file.path(root, "results/2026-09-24-disease-side-v3.5-g6-input-audit/registries/GSE233208_donor_section_manifest.tsv")
image_path <- file.path(audit_dir, "GSE233208_spatial_image_coordinate_manifest.tsv.gz")
out_dir <- file.path(root, "data/interim/G6_GSE233208_Visium_inputs")
if (!all(file.exists(c(preflight, reference_lock, object_path, manifest_path, image_path)))) {
  stop("Missing a required input gate, author-corrected object, or audit manifest.")
}
pre <- jsonlite::fromJSON(preflight, simplifyVector = FALSE)
ref_lock <- jsonlite::fromJSON(reference_lock, simplifyVector = FALSE)
valid_sha <- function(x) is.character(x) && length(x) == 1L && grepl("^[A-Fa-f0-9]{64}$", x)
donor_gate <- ref_lock$donor_independence
if (!identical(pre$status, "READY_FOR_BLIND_REFERENCE_CALIBRATION") ||
    !identical(ref_lock$status, "PASS") || !identical(ref_lock$effect_blind, TRUE) ||
    !identical(ref_lock$group_effects_calculated, FALSE) || !valid_sha(ref_lock$reference_sha256) ||
    !valid_sha(ref_lock$donor_manifest_sha256) || !identical(donor_gate$status, "PASS") ||
    !valid_sha(donor_gate$identity_crosswalk_sha256) ||
    !identical(as.integer(donor_gate$eligible_spatial_donors_checked), 20L) ||
    !identical(as.integer(donor_gate$unresolved_matches), 0L)) {
  stop("Independent reference/donor crosswalk is not locked PASS; no RCTD inputs built.")
}

donors <- read.delim(gzfile(manifest_path), check.names = FALSE, stringsAsFactors = FALSE)
images <- read.delim(gzfile(image_path), check.names = FALSE, stringsAsFactors = FALSE)
keep <- startsWith(donors$role_status, "validation_eligible_") & donors$Diagnosis %in% c("AD", "Control")
donors <- donors[keep, , drop = FALSE]
if (nrow(donors) != 20L || length(unique(donors$donor_id)) != 20L ||
    sum(donors$Diagnosis == "AD") != 10L || sum(donors$Diagnosis == "Control") != 10L) {
  stop("Author-corrected donor manifest no longer matches the 10+10 lock.")
}
obj <- readRDS(object_path)
if (!inherits(obj, "Seurat") || !"Spatial" %in% names(obj@assays)) stop("Expected corrected Seurat object with Spatial assay.")
md <- obj@meta.data
if (!"Sample" %in% names(md)) stop("Spatial sample key missing from Seurat metadata.")
counts <- methods::slot(obj@assays[["Spatial"]], "counts")
if (!inherits(counts, "sparseMatrix") || ncol(counts) != nrow(md) || any(counts@x < 0) ||
    any(abs(counts@x - round(counts@x)) > 1e-8)) stop("Spatial counts are not aligned nonnegative integers.")
if (!"images" %in% slotNames(obj) || length(obj@images) != 39L) stop("Expected the audited 39 image objects.")

ref_hash <- tolower(as.character(ref_lock$reference_sha256))
if (!grepl("^[0-9a-f]{64}$", ref_hash)) stop("Reference SHA-256 is malformed.")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
saved <- character()
for (i in seq_len(nrow(donors))) {
  row <- donors[i, , drop = FALSE]
  sample_id <- as.character(row$Seq.ID)
  donor_id <- as.character(row$donor_id)
  spots <- which(as.character(md$Sample) == sample_id)
  if (!length(spots) || anyDuplicated(colnames(counts)[spots])) stop("No spots or duplicate barcodes for ", sample_id)
  image_row <- images[images$Sample == sample_id, , drop = FALSE]
  if (nrow(image_row) != 1L || image_row$n_author_samples[[1L]] != 1L ||
      image_row$n_coordinates_unmatched_to_metadata[[1L]] != 0L) stop("Image mapping not unique for ", sample_id)
  image_obj <- obj@images[[as.character(image_row$image_name[[1L]])]]
  coord <- image_obj@coordinates
  ix <- match(colnames(counts)[spots], rownames(coord))
  if (anyNA(ix)) stop("Spatial coordinates do not cover all spots for ", sample_id)
  coord <- coord[ix, , drop = FALSE]
  x_name <- intersect(c("imagecol", "col", "pxl_col_in_fullres"), colnames(coord))
  y_name <- intersect(c("imagerow", "row", "pxl_row_in_fullres"), colnames(coord))
  if (!length(x_name) || !length(y_name)) stop("Recognized Visium image coordinate columns absent for ", sample_id)
  tissue_name <- intersect(c("tissue", "in_tissue"), colnames(coord))
  umi <- Matrix::colSums(counts[, spots, drop = FALSE])
  spot_meta <- data.frame(
    spot_id = colnames(counts)[spots], donor_id = donor_id, sample_id = sample_id,
    x_px = as.numeric(coord[[x_name[[1L]]]]), y_px = as.numeric(coord[[y_name[[1L]]]]),
    tissue = if (length(tissue_name)) as.logical(coord[[tissue_name[[1L]]]]) else NA,
    n_UMI = as.numeric(umi), stringsAsFactors = FALSE
  )
  if (any(!is.finite(spot_meta$x_px)) || any(!is.finite(spot_meta$y_px)) || any(!is.finite(spot_meta$n_UMI))) {
    stop("Non-finite coordinates or UMI values in ", sample_id)
  }
  names <- rownames(counts)
  if (anyDuplicated(names)) stop("Duplicated gene symbols require a separately locked resolution rule.")
  donor_counts <- counts[, spots, drop = FALSE]
  # Do not persist author Diagnosis, age, sex, clinical data, or pathology.
  input <- list(counts = donor_counts, metadata = spot_meta,
                image_name = as.character(image_row$image_name[[1L]]),
                reference_sha256 = ref_hash,
                source = list(accession = "GSE233208", protocol_id = "G6-PROT-2026-09-26-003",
                              author_sample_id = sample_id, donor_id = donor_id,
                              coordinates = "Visium image pixels; physical conversion remains execution-lock controlled"))
  path <- file.path(out_dir, paste0(gsub("[^A-Za-z0-9_-]", "_", donor_id), "_visium_counts.rds"))
  if (file.exists(path)) stop("Refusing to overwrite existing processed input: ", path)
  tmp <- paste0(path, ".tmp")
  saveRDS(input, tmp, compress = "gzip")
  if (!file.rename(tmp, path)) stop("Could not atomically finalize ", path)
  saved <- c(saved, path)
  message("Prepared blinded donor input ", i, "/", nrow(donors), " spots=", nrow(spot_meta))
}
if (length(saved) != 20L) stop("Did not produce exactly 20 donor-scoped inputs.")
write.csv(data.frame(path = saved, bytes = file.info(saved)$size,
                     md5 = unname(tools::md5sum(saved))),
          file.path(out_dir, "input_manifest.csv"), row.names = FALSE)
message("Completed 20 donor-scoped inputs; no group effects calculated.")
