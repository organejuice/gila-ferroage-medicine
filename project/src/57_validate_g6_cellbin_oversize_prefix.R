#!/usr/bin/env Rscript

# Effect-blind validation of oversized GSE269906 GEO transfers. The script
# preserves the downloaded source artifact, copies only the GEO-listed prefix
# into data/interim, and reads every core CellBin dataset from each candidate.

suppressPackageStartupMessages(library(rhdf5))
args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) args[[1L]] else "."
filelist_path <- file.path(root, "data/raw/G6_source_metadata_2026-09-23/GSE269906_filelist.txt")
source_dir <- file.path(root, "data/raw/GSE269906_cellbin")
reacquire_dir <- file.path(root, "data/interim/GSE269906_cellbin_reacquire")
candidate_dir <- file.path(root, "data/interim/GSE269906_cellbin_prefix_validation")
report_path <- file.path(root, "results/2026-09-24-disease-side-v3.5-g6-input-audit/tables/GSE269906_oversize_prefix_validation.tsv")
dir.create(candidate_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(report_path), recursive = TRUE, showWarnings = FALSE)

listing <- read.delim(filelist_path, check.names = FALSE, comment.char = "", stringsAsFactors = FALSE)
names(listing)[1L] <- "kind"
listing <- listing[listing$kind == "File" & grepl("\\.cellbin\\.gef$", listing$Name), , drop = FALSE]
listing$expected_bytes <- as.numeric(gsub(",", "", listing$Size, fixed = TRUE))

source_files <- list.files(reacquire_dir, pattern = "\\.reassembled$", full.names = TRUE)
reacquired <- length(source_files) > 0L
if (!reacquired) {
  source_files <- list.files(source_dir, pattern = "\\.cellbin\\.gef$", full.names = TRUE)
}
rows <- list()

copy_prefix <- function(source, target, n_bytes) {
  input <- file(source, open = "rb")
  output <- file(target, open = "wb")
  on.exit(close(input), add = TRUE)
  on.exit(close(output), add = TRUE)
  remaining <- as.double(n_bytes)
  while (remaining > 0) {
    n <- as.integer(min(8 * 1024 * 1024, remaining))
    chunk <- readBin(input, what = "raw", n = n)
    if (length(chunk) != n) stop("Unexpected EOF while prefix-copying ", basename(source))
    writeBin(chunk, output)
    remaining <- remaining - length(chunk)
  }
}

for (source in source_files) {
  filename <- if (reacquired) sub("\\.reassembled$", "", basename(source)) else basename(source)
  expected <- listing$expected_bytes[match(filename, listing$Name)]
  actual <- file.info(source)$size
  if (is.na(expected) || (!reacquired && actual <= expected)) next
  if (reacquired && actual != expected) stop("Reacquired file size does not match the official GEO listing: ", filename)
  candidate <- file.path(candidate_dir, filename)
  copy_prefix(source, candidate, expected)
  hdf5_candidate <- tempfile(pattern = "g6_prefix_", fileext = ".gef", tmpdir = tempdir())
  if (!file.copy(candidate, hdf5_candidate, overwrite = TRUE)) stop("Could not stage ASCII-path validation copy for ", filename)

  core_paths <- c("/cellBin/cell", "/cellBin/cellBorder", "/cellBin/geneExp", "/cellBin/gene")
  checks <- setNames(rep(FALSE, length(core_paths)), core_paths)
  dimensions <- setNames(rep("", length(core_paths)), core_paths)
  core_values <- setNames(vector("list", length(core_paths)), core_paths)
  error_text <- ""
  tryCatch({
    for (dataset in core_paths) {
      value <- h5read(hdf5_candidate, dataset)
      core_values[[dataset]] <- value
      dims <- dim(value)
      dimensions[[dataset]] <- if (is.null(dims)) paste0("length=", length(value)) else paste(dims, collapse = "x")
      checks[[dataset]] <- length(value) > 0
    }
    cell <- core_values[["/cellBin/cell"]]
    if (is.matrix(cell)) cell <- as.data.frame(cell)
    if (!all(c("id", "x", "y", "geneCount", "expCount") %in% names(cell))) {
      stop("Cell table lacks expected identifiers, coordinates, or count fields.")
    }
    if (anyNA(cell$id) || anyNA(cell$x) || anyNA(cell$y) || any(cell$geneCount < 0) || any(cell$expCount < 0)) {
      stop("Cell table contains missing coordinates/IDs or negative count values.")
    }
  }, error = function(e) {
    error_text <<- conditionMessage(e)
  })
  unlink(hdf5_candidate)

  rows[[length(rows) + 1L]] <- data.frame(
    file = filename,
    expected_bytes = expected,
    original_bytes = if (reacquired && file.exists(file.path(source_dir, filename))) file.info(file.path(source_dir, filename))$size else actual,
    reacquired_bytes = actual,
    extra_trailing_bytes = if (reacquired) 0 else actual - expected,
    candidate_bytes = file.info(candidate)$size,
    core_dataset_checks = paste(names(checks), checks, sep = "=", collapse = ";"),
    core_dataset_dimensions = paste(names(dimensions), dimensions, sep = "=", collapse = ";"),
    hdf5_core_data_valid = all(checks) && !nzchar(error_text),
    error = error_text,
    validation_source = normalizePath(source, winslash = "/", mustWork = TRUE),
    candidate_path = normalizePath(candidate, winslash = "/", mustWork = TRUE),
    stringsAsFactors = FALSE
  )
}

if (!length(rows)) stop("No oversized cellbin files were found for validation.")
report <- do.call(rbind, rows)
write.table(report, report_path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
writeLines(capture.output(print(report, row.names = FALSE)), con = stdout())
if (!all(report$hdf5_core_data_valid)) quit(status = 2L)
