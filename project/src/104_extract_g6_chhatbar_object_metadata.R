#!/usr/bin/env Rscript
# Read a public Chhatbar Seurat object and export only source-level audit fields.
# No expression values or donor-level effects are computed.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("Usage: Rscript src/104_extract_g6_chhatbar_object_metadata.R INPUT_RDS OUTPUT_DIR PLATFORM")
}
input <- args[[1]]
output_dir <- args[[2]]
platform <- args[[3]]
if (!platform %in% c("cosmx", "visium_hd")) stop("Unknown platform")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
object <- readRDS(input)
if (!inherits(object, "Seurat")) stop("Expected a Seurat object")
metadata <- object@meta.data
write.table(data.frame(gene = rownames(object)),
            file.path(output_dir, paste0(platform, "_features.tsv")),
            sep = "\t", row.names = FALSE, quote = FALSE)

fields <- data.frame(
  field = names(metadata),
  class = vapply(metadata, function(x) paste(class(x), collapse = "/"), character(1)),
  missing_n = vapply(metadata, function(x) sum(is.na(x)), integer(1)),
  unique_n = vapply(metadata, function(x) length(unique(x)), integer(1)),
  stringsAsFactors = FALSE
)
write.table(fields, file.path(output_dir, paste0(platform, "_metadata_fields.tsv")),
            sep = "\t", row.names = FALSE, quote = FALSE)

keys <- if (platform == "cosmx") {
  c("Sample", "Pathology")
} else {
  c("Subject_ID", "Sample_id", "orig.ident", "Region", "Tissue_Region",
    "Tissue_Type", "Reported_Sex", "Age", "Clinical_Brain_Diagnosis", "Clinical_Brain_Diagnosis_Short",
    "Disease_region", "Sample_id_disease")
}
keys <- intersect(keys, names(metadata))
if (!length(keys)) stop("No expected grouping keys found")
sample_rows <- metadata[, keys, drop = FALSE]
sample_rows[] <- lapply(sample_rows, function(x) {
  value <- as.character(x)
  value[is.na(value)] <- ""
  value
})
sample_rows$audited_cells <- 1L
summary <- aggregate(sample_rows$audited_cells, sample_rows[, keys, drop = FALSE], sum)
names(summary)[ncol(summary)] <- "cell_count"
if ("Sample.FOV" %in% names(metadata)) {
  fov <- aggregate(as.character(metadata$Sample.FOV), sample_rows[, keys, drop = FALSE],
                   function(x) length(unique(x)))
  names(fov)[ncol(fov)] <- "fov_count"
  summary <- merge(summary, fov, by = keys, all.x = TRUE, sort = FALSE)
}
write.table(summary, file.path(output_dir, paste0(platform, "_sample_summary.tsv")),
            sep = "\t", row.names = FALSE, quote = FALSE)

inventory <- data.frame(
  platform = platform,
  class = paste(class(object), collapse = "/"),
  cells = ncol(object),
  features = nrow(object),
  assays = paste(names(object@assays), collapse = ";"),
  images = length(object@images),
  metadata_fields = ncol(metadata),
  counts_class = class(object@assays[[1]]@counts)[1],
  counts_nonzero = if (inherits(object@assays[[1]]@counts, "sparseMatrix"))
    length(object@assays[[1]]@counts@x) else sum(object@assays[[1]]@counts != 0)
)
write.table(inventory, file.path(output_dir, paste0(platform, "_inventory.tsv")),
            sep = "\t", row.names = FALSE, quote = FALSE)

image_rows <- lapply(names(object@images), function(nm) {
  im <- object@images[[nm]]
  slots <- slotNames(im)
  bounds <- if ("boundaries" %in% slots) names(slot(im, "boundaries")) else character()
  data.frame(image_name = nm, image_class = class(im)[1],
             boundary_types = paste(bounds, collapse = ";"),
             has_raster_slot = any(c("image", "raster", "data") %in% slots),
             stringsAsFactors = FALSE)
})
write.table(do.call(rbind, image_rows),
            file.path(output_dir, paste0(platform, "_spatial_images.tsv")),
            sep = "\t", row.names = FALSE, quote = FALSE)
