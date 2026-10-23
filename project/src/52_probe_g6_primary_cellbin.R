#!/usr/bin/env Rscript

# Effect-blind structural and gene-coverage audit of complete GSE269906
# cell-bin GEFs. Each source file is copied to a temporary ASCII path because
# the installed rhdf5 build cannot open this project's Unicode path directly.

suppressPackageStartupMessages(library(rhdf5))
args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) args[[1L]] else "."
out_dir <- file.path(root, "results/2026-09-24-disease-side-v3.5-g6-input-audit")
manifest_path <- file.path(out_dir, "registries/GSE269906_file_acquisition_manifest.tsv")
g1_path <- file.path(out_dir, "registries/G1_microglia_signed_sender_244.tsv.gz")
if (!file.exists(manifest_path) || !file.exists(g1_path)) stop("Run scripts 50 and 51 before the GEF audit.")

files <- read.delim(manifest_path, check.names = FALSE, stringsAsFactors = FALSE)
files <- files[files$file_type == "cellbin" & files$status == "complete", , drop = FALSE]
donors <- read.delim(file.path(out_dir, "registries/GSE269906_donor_section_manifest.tsv"), check.names = FALSE, stringsAsFactors = FALSE)
files$source_diagnosis <- donors$source_diagnosis[match(files$GSM, donors$GSM)]
if (anyNA(files$source_diagnosis)) stop("A complete GEF file lacks an author-verified donor diagnosis mapping.")
g1 <- read.delim(gzfile(g1_path), check.names = FALSE, stringsAsFactors = FALSE)
g1_symbols <- toupper(as.character(g1$gene))
if (length(g1_symbols) != 244L) stop("G1 frozen signature does not contain 244 genes.")

audit_one <- function(row) {
  source <- file.path(root, "data/raw/GSE269906_cellbin", row$file)
  if (!file.exists(source) || file.info(source)$size != row$expected_bytes) stop("Complete-file manifest mismatch: ", row$file)
  temp_file <- tempfile(pattern = "g6_gef_", fileext = ".gef", tmpdir = tempdir())
  on.exit(unlink(temp_file), add = TRUE)
  if (!file.copy(source, temp_file, overwrite = TRUE)) stop("Unable to create temporary ASCII-path GEF copy.")
  listing <- h5ls(temp_file, recursive = TRUE, all = TRUE)
  root_attrs <- h5readAttributes(temp_file, "/")
  cell_attrs <- h5readAttributes(temp_file, "/cellBin")
  cell <- h5read(temp_file, "/cellBin/cell")
  if (is.matrix(cell)) cell <- as.data.frame(cell)
  if (!is.data.frame(cell) && !is.list(cell)) stop("Unexpected /cellBin/cell format in ", row$file)
  cell_cols <- names(cell)
  if (is.null(cell_cols)) stop("Cannot resolve cell fields in ", row$file)
  field <- function(name) {
    if (!name %in% cell_cols) return(NULL)
    as.numeric(cell[[name]])
  }
  gene_field_candidates <- c("geneCount", "gene_count")
  exp_field_candidates <- c("expCount", "MIDCount", "exp_count")
  area_field_candidates <- c("area", "Area")
  pick_field <- function(candidates) candidates[candidates %in% cell_cols][1L]
  gene_field <- pick_field(gene_field_candidates)
  exp_field <- pick_field(exp_field_candidates)
  area_field <- pick_field(area_field_candidates)
  type_values <- field("cellTypeID")
  cluster_values <- field("clusterID")
  border <- listing[listing$group == "/cellBin" & listing$name == "cellBorder", , drop = FALSE]
  gene_rows <- listing[listing$group == "/cellBin" & listing$name == "gene", , drop = FALSE]
  gene_table <- h5read(temp_file, "/cellBin/gene")
  if (is.matrix(gene_table)) gene_table <- as.data.frame(gene_table)
  gene_cols <- names(gene_table)
  if (is.null(gene_cols)) gene_cols <- colnames(gene_table)
  possible_names <- c("gene", "geneName", "gene_name", "name")
  name_col <- possible_names[possible_names %in% gene_cols][1L]
  if (is.na(name_col) || !length(name_col)) stop("Unable to locate gene names in /cellBin/gene for ", row$file)
  symbols <- toupper(as.character(gene_table[[name_col]]))
  symbols <- symbols[!is.na(symbols) & nzchar(symbols)]
  unique_symbols <- unique(symbols)
  matched <- intersect(unique_symbols, g1_symbols)
  gene_data <- data.frame(
    accession = "GSE269906", donor_id = row$donor_id, GSM = row$GSM,
    gene_symbol_upper = unique_symbols,
    G1_signed_sender_present = unique_symbols %in% g1_symbols,
    stringsAsFactors = FALSE
  )
  cell_types <- if (is.null(type_values)) "field_absent" else paste(sort(unique(type_values)), collapse = ";")
  clusters <- if (is.null(cluster_values)) "field_absent" else paste(sort(unique(cluster_values)), collapse = ";")
  attr_value <- function(x, key) {
    if (!key %in% names(x)) return(NA_character_)
    sep <- if (key == "geftool_ver") "." else ";"
    paste(as.character(x[[key]]), collapse = sep)
  }
  summary_values <- list(
    accession = "GSE269906", donor_id = row$donor_id, GSM = row$GSM,
    diagnosis = row$source_diagnosis, file = row$file,
    n_cellbin_rows = length(field("id")), cell_fields = paste(cell_cols, collapse = ";"),
    median_geneCount = if (!is.null(gene_field)) median(as.numeric(cell[[gene_field]]), na.rm = TRUE) else NA_real_,
    median_expCount = if (!is.null(exp_field)) median(as.numeric(cell[[exp_field]]), na.rm = TRUE) else NA_real_,
    median_area = if (!is.null(area_field)) median(as.numeric(cell[[area_field]]), na.rm = TRUE) else NA_real_,
    cellTypeID_values = cell_types, clusterID_values = clusters,
    cellBorder_dimensions = if (nrow(border)) paste(border$dim, collapse = "x") else NA_character_,
    gene_table_rows = if (nrow(gene_rows)) as.character(gene_rows$dim[1L]) else NA_character_,
    n_unique_gene_symbols = length(unique_symbols),
    n_G1_signed_sender_genes_present = length(matched),
    G1_signed_sender_coverage = length(matched) / length(g1_symbols),
    GEF_resolution_attribute = attr_value(root_attrs, "resolution"),
    GEF_omics_attribute = attr_value(root_attrs, "omics"),
    GEF_version_attribute = attr_value(root_attrs, "version"),
    GEF_tool_version_attribute = attr_value(root_attrs, "geftool_ver"),
    offsetX = attr_value(root_attrs, "offsetX"), offsetY = attr_value(root_attrs, "offsetY"),
    cellBin_attributes = paste(names(cell_attrs), collapse = ";"),
    segmentation_representation = if (nrow(border)) "cellBorder polygon dataset present" else "cellBorder dataset absent",
    nucleus_assignment = "not established by audited GEF cell fields",
    cell_type_annotation = if (identical(cell_types, "0")) "unassigned_default_ID_0" else "IDs_present_labels_require_source_audit"
  )
  bad_lengths <- lengths(summary_values) != 1L
  if (any(bad_lengths)) {
    message("Non-scalar summary fields for ", row$file, ": ", paste(names(summary_values)[bad_lengths], lengths(summary_values)[bad_lengths], collapse = ", "))
    print(root_attrs)
    print(cell_attrs)
    print(border)
    print(gene_rows)
    stop("Structural summary is not scalar per donor.")
  }
  summary <- as.data.frame(summary_values, stringsAsFactors = FALSE)
  list(summary = summary, genes = gene_data)
}

results <- lapply(seq_len(nrow(files)), function(i) audit_one(files[i, , drop = FALSE]))
summary <- do.call(rbind, lapply(results, `[[`, "summary"))
genes <- do.call(rbind, lapply(results, `[[`, "genes"))

write_gz <- function(x, file) {
  con <- gzfile(file, open = "wt", compression = 6)
  on.exit(close(con))
  write.table(x, con, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
}
write_gz(summary, file.path(out_dir, "tables/GSE269906_cellbin_structural_gene_coverage.tsv.gz"))
write_gz(genes, file.path(out_dir, "tables/GSE269906_cellbin_gene_catalog.tsv.gz"))
cat("complete_cellbin_files_audited=", nrow(summary), "\n", sep = "")
cat("G1_signed_sender_coverage_range=", paste(range(summary$G1_signed_sender_coverage), collapse = ".."), "\n", sep = "")
cat("cellbin_rows_total=", sum(summary$n_cellbin_rows), "\n", sep = "")
cat("files_with_unassigned_default_cellTypeID=", sum(summary$cell_type_annotation == "unassigned_default_ID_0"), "\n", sep = "")
