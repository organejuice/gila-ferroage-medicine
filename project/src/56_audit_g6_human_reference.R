#!/usr/bin/env Rscript

# Effect-blind audit of the local human snRNA-seq reference proposed for
# cell-type inference in the G6 Visium cohort.

suppressPackageStartupMessages(library(rhdf5))
args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) args[[1L]] else "."
out_dir <- file.path(root, "results/2026-09-24-disease-side-v3.5-g6-input-audit")
raw_h5 <- file.path(root, "data/raw/D8_GSE174367/GSE174367_snRNA-seq_filtered_feature_bc_matrix.h5")
raw_meta <- file.path(root, "data/raw/D8_GSE174367/GSE174367_snRNA-seq_cell_meta.csv.gz")
if (!file.exists(raw_h5) || !file.exists(raw_meta)) stop("Missing local D8/GSE174367 reference inputs.")

meta <- read.csv(gzfile(raw_meta), stringsAsFactors = FALSE, check.names = FALSE)
required <- c("Barcode", "SampleID", "Diagnosis", "Cell.Type", "Age", "Sex")
if (!all(required %in% names(meta))) stop("D8 cell metadata lacks required cell-type/donor fields.")
if (anyNA(meta$Barcode) || anyDuplicated(meta$Barcode)) stop("D8 metadata barcodes are missing or duplicated.")

temp_h5 <- tempfile(pattern = "g6_human_reference_", fileext = ".h5", tmpdir = tempdir())
if (!file.copy(raw_h5, temp_h5, overwrite = TRUE)) stop("Cannot stage the reference H5 at a temporary ASCII path.")
barcodes <- as.character(h5read(temp_h5, "/matrix/barcodes"))
feature_names <- as.character(h5read(temp_h5, "/matrix/features/name"))
feature_ids <- as.character(h5read(temp_h5, "/matrix/features/id"))
feature_types <- as.character(h5read(temp_h5, "/matrix/features/feature_type"))
feature_genomes <- as.character(h5read(temp_h5, "/matrix/features/genome"))
shape <- as.integer(h5read(temp_h5, "/matrix/shape"))
counts_data <- h5read(temp_h5, "/matrix/data")

matched_meta <- match(barcodes, meta$Barcode)
matched_matrix <- match(meta$Barcode, barcodes)
if (anyNA(matched_matrix) || anyDuplicated(meta$Barcode)) {
  stop("A source-labelled D8 metadata barcode is absent from or duplicated in the expression matrix.")
}
unlabeled_matrix_barcodes <- barcodes[is.na(matched_meta)]
if (anyDuplicated(feature_ids) || anyNA(feature_ids) || anyNA(feature_names)) stop("D8 reference feature IDs are duplicated or feature identifiers are missing.")
if (length(shape) != 2L || shape[1L] != length(feature_names) || shape[2L] != length(barcodes)) {
  stop("D8 10x matrix shape does not match its features and barcodes.")
}
if (any(counts_data < 0) || any(abs(counts_data - round(counts_data)) > 1e-8)) {
  stop("D8 reference count layer is not nonnegative integer counts.")
}
if (anyNA(meta$Cell.Type) || any(!nzchar(meta$Cell.Type))) stop("D8 source cell-type labels are missing.")

donor_groups <- split(meta, meta$SampleID)
sample_manifest <- do.call(rbind, lapply(names(donor_groups), function(sample) {
  current <- donor_groups[[sample]]
  unique_value <- function(column) {
    values <- unique(as.character(current[[column]]))
    values <- values[!is.na(values) & nzchar(values)]
    if (length(values) == 1L) values else if (!length(values)) NA_character_ else paste(sort(values), collapse = ";")
  }
  data.frame(
    accession = "GSE174367",
    donor_id = sample,
    Diagnosis = unique_value("Diagnosis"),
    Age = unique_value("Age"),
    Sex = unique_value("Sex"),
    n_nuclei = nrow(current),
    n_cell_types = length(unique(current$Cell.Type)),
    Cell.Type_labels = paste(sort(unique(current$Cell.Type)), collapse = ";"),
    n_MG = sum(current$Cell.Type == "MG"),
    n_ASC = sum(current$Cell.Type == "ASC"),
    n_MG_ASC = sum(current$Cell.Type %in% c("MG", "ASC")),
    stringsAsFactors = FALSE
  )
}))

cell_type_counts <- as.data.frame(table(meta$Cell.Type), stringsAsFactors = FALSE)
names(cell_type_counts) <- c("Cell.Type", "n_nuclei")
feature_catalog <- data.frame(
  accession = "GSE174367",
  feature_id = feature_ids,
  feature_name = feature_names,
  feature_type = feature_types,
  genome = feature_genomes,
  stringsAsFactors = FALSE
)
expression_audit <- data.frame(
  accession = "GSE174367",
  assay = "10x snRNA-seq filtered feature-barcode matrix",
  n_features = shape[1L],
  n_barcodes = shape[2L],
  n_nonzero_entries = length(counts_data),
  nonnegative_integer_counts = TRUE,
  metadata_barcodes_unique = !anyDuplicated(meta$Barcode),
  labelled_metadata_barcodes_all_present_in_matrix = TRUE,
  n_unlabelled_matrix_barcodes_excluded = length(unlabeled_matrix_barcodes),
  n_donors = length(donor_groups),
  n_cell_type_labels = length(unique(meta$Cell.Type)),
  n_duplicate_feature_names = sum(duplicated(feature_names)),
  metadata_cell_type_missing = sum(is.na(meta$Cell.Type) | !nzchar(meta$Cell.Type)),
  stringsAsFactors = FALSE
)

write_gz <- function(x, path) {
  con <- gzfile(path, "wt", compression = 6)
  on.exit(close(con))
  write.table(x, con, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
}
write_gz(sample_manifest, file.path(out_dir, "registries/GSE174367_human_reference_donor_manifest.tsv.gz"))
write_gz(cell_type_counts, file.path(out_dir, "tables/GSE174367_human_reference_celltype_counts.tsv.gz"))
write_gz(feature_catalog, file.path(out_dir, "registries/GSE174367_human_reference_feature_catalog.tsv.gz"))
write_gz(expression_audit, file.path(out_dir, "tables/GSE174367_human_reference_expression_audit.tsv.gz"))
write_gz(data.frame(barcode = unlabeled_matrix_barcodes), file.path(out_dir, "registries/GSE174367_unlabeled_matrix_barcodes_excluded.tsv.gz"))
unlink(temp_h5)

cat("reference_accession=GSE174367\n")
cat("matrix_dimensions_features_x_nuclei=", paste(shape, collapse = "x"), "\n", sep = "")
cat("metadata_nuclei=", nrow(meta), "\n", sep = "")
cat("unlabeled_matrix_barcodes_excluded=", length(unlabeled_matrix_barcodes), "\n", sep = "")
cat("donors=", length(donor_groups), "\n", sep = "")
cat("unique_cell_types=", length(unique(meta$Cell.Type)), "\n", sep = "")
cat("MG_nuclei=", sum(meta$Cell.Type == "MG"), "\n", sep = "")
cat("ASC_nuclei=", sum(meta$Cell.Type == "ASC"), "\n", sep = "")
cat("barcode_match_n=", sum(!is.na(matched_meta)), "\n", sep = "")
cat("reference_counts_nonnegative_integer=TRUE\n")
