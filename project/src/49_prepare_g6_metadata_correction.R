#!/usr/bin/env Rscript

# Reproduce the author-documented GSE233208 patient-metadata correction.
# The source Seurat object is read only and is never overwritten.

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1L) args[[1L]] else "."
if (!file.exists(file.path(root, "docs/protocol/disease_side_v3.6_G6_GSE233208.md"))) {
  stop("Run this script from the project root, or pass the project root as its first argument.")
}
raw_gz <- file.path(root, "data/raw/GSE233208/GSE233208_Human_visium_ADDS_seurat_processed.rds.gz")
interim_dir <- file.path(root, "data/interim/GSE233208")
interim_rds <- file.path(interim_dir, "GSE233208_Human_visium_ADDS_seurat_processed.rds")
author_table <- file.path(root, "data/raw/G6_source_metadata_2026-09-23/GSE233208_author_Supplementary_Table1.tsv.gz")
out_dir <- file.path(root, "results/2026-09-26-g6-gse233208-visium-input-audit/tables")

dir.create(interim_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(interim_rds)) {
  if (!file.exists(raw_gz)) stop("Missing immutable GEO source object: ", raw_gz)
  con_in <- gzfile(raw_gz, open = "rb")
  con_out <- file(interim_rds, open = "wb")
  repeat {
    block <- readBin(con_in, what = "raw", n = 8L * 1024L * 1024L)
    if (!length(block)) break
    writeBin(block, con_out)
  }
  close(con_out)
  close(con_in)
}

obj <- readRDS(interim_rds)
if (!inherits(obj, "Seurat")) stop("Expected a Seurat object; received: ", paste(class(obj), collapse = ","))
md <- obj@meta.data
if (!"Sample" %in% names(md)) stop("Object has no Sample key required by the author correction.")
auth <- read.delim(gzfile(author_table), check.names = FALSE, stringsAsFactors = FALSE)
required <- c("Seq.ID", "Case.Year", "Case.Num", "Diagnosis")
if (!all(required %in% names(auth))) stop("Author table lacks required fields: ", paste(setdiff(required, names(auth)), collapse = ","))
if (anyNA(auth$Seq.ID) || anyDuplicated(auth$Seq.ID)) stop("Author Seq.ID values are missing or duplicated.")

spot_sample <- as.character(md$Sample)
if (anyNA(spot_sample) || any(!nzchar(spot_sample))) stop("Object Sample key contains missing or empty values.")
if (!setequal(unique(spot_sample), auth$Seq.ID)) {
  stop("Object Sample IDs and author Seq.ID values are not a complete one-to-one set.")
}
idx <- match(spot_sample, auth$Seq.ID)
if (anyNA(idx)) stop("Unmatched object Sample IDs remain after matching to author Seq.ID.")

donor_key <- paste(auth$Case.Year, auth$Case.Num, sep = "_")
if (anyNA(donor_key) || anyDuplicated(donor_key)) stop("Author Case.Year/Case.Num donor keys are missing or duplicated.")
dx_counts <- table(factor(auth$Diagnosis, levels = c("Control", "AD", "earlyAD", "AD_DS")))
if (unname(dx_counts[["Control"]]) != 10L || unname(dx_counts[["AD"]]) != 10L) {
  stop("Expected 10 Control and 10 late sporadic AD cases in the author table.")
}

# Preserve the uncorrected object metadata as an audit artifact, then remove
# patient/section fields from the object and replace them with author values.
original <- data.frame(spot_barcode = rownames(md), md, check.names = FALSE)
conflicting_patient_fields <- c(
  "orig.ident", "SampleID", "Case.Year", "Case.Num", "Diagnosis", "Sex", "Age",
  "PMI", "APoE", "Clinical.Syndrome", "Tangle.Stage", "Plaque.Stage", "RIN",
  "Slide", "Sectioning", "LibraryDay1", "LibraryDay2", "CaptureArea"
)
drop_fields <- intersect(names(md), conflicting_patient_fields)
corrected <- md[, setdiff(names(md), drop_fields), drop = FALSE]
corrected$spot_barcode <- rownames(md)
corrected$Seq.ID <- auth$Seq.ID[idx]
for (field in setdiff(names(auth), "Seq.ID")) corrected[[field]] <- auth[[field]][idx]
corrected <- corrected[, c("spot_barcode", setdiff(names(corrected), "spot_barcode")), drop = FALSE]

spots_by_sample <- table(factor(spot_sample, levels = auth$Seq.ID))
sample_manifest <- auth
sample_manifest$n_spots_object <- as.integer(spots_by_sample)
if ("n_Visium_spots" %in% names(auth)) {
  sample_manifest$spot_count_delta_object_minus_author <- sample_manifest$n_spots_object - as.integer(auth$n_Visium_spots)
}

image_manifest <- lapply(names(obj@images), function(image_name) {
  image <- obj@images[[image_name]]
  coordinates <- image@coordinates
  coord_idx <- match(rownames(coordinates), rownames(md))
  matched_samples <- unique(as.character(md$Sample[coord_idx[!is.na(coord_idx)]]))
  author_idx <- match(matched_samples, auth$Seq.ID)
  scales <- unlist(image@scale.factors)
  image_dims <- dim(image@image)
  data.frame(
    accession = "GSE233208",
    image_name = image_name,
    n_image_coordinates = nrow(coordinates),
    n_coordinates_matched_to_metadata = sum(!is.na(coord_idx)),
    n_coordinates_unmatched_to_metadata = sum(is.na(coord_idx)),
    n_author_samples = length(matched_samples),
    Sample = if (length(matched_samples) == 1L) matched_samples else paste(matched_samples, collapse = ";"),
    Case.Year = if (length(author_idx) == 1L) auth$Case.Year[author_idx] else NA,
    Case.Num = if (length(author_idx) == 1L) auth$Case.Num[author_idx] else NA,
    Diagnosis = if (length(author_idx) == 1L) auth$Diagnosis[author_idx] else NA,
    n_spots_object_for_sample = if (length(matched_samples) == 1L) as.integer(spots_by_sample[matched_samples]) else NA,
    coordinate_columns = paste(colnames(coordinates), collapse = ";"),
    tissue_flag_present = "tissue" %in% colnames(coordinates),
    image_class = paste(class(image), collapse = ";"),
    image_dimensions = paste(image_dims, collapse = "x"),
    spot_scale_factor = if ("spot" %in% names(scales)) unname(scales[["spot"]]) else NA,
    fiducial_scale_factor = if ("fiducial" %in% names(scales)) unname(scales[["fiducial"]]) else NA,
    hires_scale_factor = if ("hires" %in% names(scales)) unname(scales[["hires"]]) else NA,
    lowres_scale_factor = if ("lowres" %in% names(scales)) unname(scales[["lowres"]]) else NA,
    stringsAsFactors = FALSE
  )
})
image_manifest <- if (length(image_manifest)) do.call(rbind, image_manifest) else data.frame()
if (nrow(image_manifest) != nrow(auth) || any(image_manifest$n_author_samples != 1L) ||
    any(image_manifest$n_coordinates_unmatched_to_metadata != 0L) ||
    !setequal(image_manifest$Sample, auth$Seq.ID)) {
  stop("Spatial image-to-spot-to-author-donor mapping is incomplete or ambiguous.")
}

assay_names <- names(obj@assays)
if (!"Spatial" %in% assay_names || !"counts" %in% methods::slotNames(obj@assays[["Spatial"]])) {
  stop("The expected Spatial assay count layer is missing.")
}
counts <- methods::slot(obj@assays[["Spatial"]], "counts")
if (is.null(dim(counts)) || !identical(as.integer(dim(counts)), as.integer(dim(obj)))) {
  stop("Spatial count dimensions do not match the Seurat object dimensions.")
}
counts_summary <- data.frame(
  accession = "GSE233208",
  assay = "Spatial",
  layer = "counts",
  matrix_class = paste(class(counts), collapse = ";"),
  n_features = nrow(counts),
  n_spots = ncol(counts),
  nonzero_entries = length(counts@x),
  minimum_nonzero = min(counts@x),
  maximum_nonzero = max(counts@x),
  noninteger_nonzero_entries = sum(abs(counts@x - round(counts@x)) > 1e-8),
  feature_id_example = paste(head(rownames(counts), 5L), collapse = ";"),
  stringsAsFactors = FALSE
)

write_gz_tsv <- function(x, path) {
  con <- gzfile(path, open = "wt", compression = 6)
  on.exit(close(con))
  write.table(x, con, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
}

write_gz_tsv(original, file.path(out_dir, "GSE233208_original_object_metadata.tsv.gz"))
write_gz_tsv(corrected, file.path(out_dir, "GSE233208_author_corrected_spot_metadata.tsv.gz"))
write_gz_tsv(sample_manifest, file.path(out_dir, "GSE233208_author_corrected_sample_donor_manifest.tsv.gz"))
write_gz_tsv(image_manifest, file.path(out_dir, "GSE233208_spatial_image_coordinate_manifest.tsv.gz"))
write_gz_tsv(counts_summary, file.path(out_dir, "GSE233208_expression_layer_audit.tsv.gz"))
write_gz_tsv(data.frame(object_patient_field_removed = drop_fields), file.path(out_dir, "GSE233208_removed_conflicting_object_fields.tsv.gz"))

audit <- c(
  paste0("object_class=", paste(class(obj), collapse = ";")),
  paste0("object_dimensions_features_x_spots=", paste(dim(obj), collapse = "x")),
  paste0("object_metadata_rows=", nrow(md)),
  paste0("object_unique_Sample=", length(unique(spot_sample))),
  paste0("author_rows=", nrow(auth)),
  paste0("author_unique_Seq.ID=", length(unique(auth$Seq.ID))),
  paste0("author_unique_Case.Year_Case.Num=", length(unique(donor_key))),
  paste0("matched_spot_rows=", sum(!is.na(idx))),
  paste0("unmatched_spot_rows=", sum(is.na(idx))),
  paste0("unmatched_author_samples=", sum(!auth$Seq.ID %in% unique(spot_sample))),
  paste0("late_sporadic_AD_donors=", unname(dx_counts[["AD"]])),
  paste0("Control_donors=", unname(dx_counts[["Control"]])),
  paste0("earlyAD_excluded_donors=", unname(dx_counts[["earlyAD"]])),
  paste0("AD_DS_excluded_donors=", unname(dx_counts[["AD_DS"]])),
  paste0("spatial_images=", nrow(image_manifest)),
  paste0("images_one_to_one_to_author_Samples=", all(image_manifest$n_author_samples == 1L) && setequal(image_manifest$Sample, auth$Seq.ID)),
  paste0("image_coordinates_matched_to_metadata=", sum(image_manifest$n_coordinates_matched_to_metadata)),
  paste0("image_coordinates_unmatched_to_metadata=", sum(image_manifest$n_coordinates_unmatched_to_metadata)),
  paste0("expression_counts_layer=", counts_summary$matrix_class, " ", counts_summary$n_features, "x", counts_summary$n_spots),
  paste0("expression_nonzero_entries=", counts_summary$nonzero_entries),
  paste0("expression_noninteger_entries=", counts_summary$noninteger_nonzero_entries),
  paste0("object_fields_removed=", paste(drop_fields, collapse = ",")),
  paste0("all_author_Sample_n_spot_counts_match=", if ("spot_count_delta_object_minus_author" %in% names(sample_manifest)) all(sample_manifest$spot_count_delta_object_minus_author == 0) else NA)
)
writeLines(audit, file.path(out_dir, "GSE233208_metadata_correction_audit.txt"), useBytes = TRUE)
cat(paste(audit, collapse = "\n"), "\n")
