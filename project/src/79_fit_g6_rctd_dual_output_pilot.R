#!/usr/bin/env Rscript

# G6-DIAG-2026-09-26-001: diagnosis-blind two-donor RCTD output audit.
suppressPackageStartupMessages({ library(Matrix); library(spacexr) })
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: 79_fit_g6_rctd_dual_output_pilot.R <project_root> <GSM>")
root <- args[[1L]]
if (!dir.exists(root)) stop("Project root does not exist")
gsm <- args[[2L]]
if (!gsm %in% c("GSM8330063", "GSM8330064")) stop("GSM is outside the locked pilot")
if (as.character(packageVersion("spacexr")) != "2.2.1") stop("spacexr version mismatch")

interim <- file.path(root, "data/interim/G6_GSE269906_bin50_tissue")
auditdir <- file.path(root, "results/2026-09-24-disease-side-v3.5-g6-input-audit/registries")
input_path <- file.path(root, "data/interim/G6_GSE269906_bin50_full",
                        paste0(gsm, "_bin50_centroid_supported_counts.rds"))
ref_path <- file.path(interim, "G6_RCTD_reference_cells.rds")
g1_path <- file.path(auditdir, "G1_microglia_signed_sender_244.tsv.gz")
g5_path <- file.path(auditdir, "G5_receiver_one_to_one_human_ortholog_registry.tsv.gz")
old_path <- file.path(root, "data/interim/G6_GSE269906_RCTD",
                      paste0(gsm, "_RCTD_weights_effect_blind.rds"))
config_path <- file.path(root, "configs/g6_rctd_dual_output_diagnostic_001.yaml")
script_path <- file.path(root, "src/79_fit_g6_rctd_dual_output_pilot.R")
paths <- c(input_path, ref_path, g1_path, g5_path, old_path, config_path, script_path)
if (!all(file.exists(paths))) stop("A locked input or config is missing: ", paste(paths[!file.exists(paths)], collapse = "; "))
path_md5 <- unname(tools::md5sum(paths))
if (anyNA(path_md5)) stop("Cannot hash all locked inputs")
names(path_md5) <- paths
outdir <- file.path(root, "data/interim/G6_GSE269906_RCTD_dual_output_pilot")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
out_path <- file.path(outdir, paste0(gsm, "_dual_output_pilot.rds"))
if (file.exists(out_path)) stop("Existing final pilot output; refusing overwrite: ", out_path)

reference <- readRDS(ref_path)
spatial <- readRDS(input_path)
g1 <- read.delim(gzfile(g1_path), check.names = FALSE, stringsAsFactors = FALSE)
g5 <- read.delim(gzfile(g5_path), check.names = FALSE, stringsAsFactors = FALSE)
state_genes <- unique(toupper(c(g1$gene, g5$human_one_to_one_symbol)))
state_genes <- state_genes[!is.na(state_genes) & nzchar(state_genes)]
holdout_candidates <- c("AIF1", "C1QA", "C1QB", "C1QC", "CX3CR1", "LST1", "P2RY12", "TMEM119", "TREM2", "TYROBP",
                        "ALDH1L1", "GFAP", "GJA1", "GLUL", "S100B", "SLC1A2", "SLC1A3", "SOX9", "CD74", "AQP4")
holdout_markers <- setdiff(holdout_candidates, state_genes)
features <- intersect(rownames(spatial$counts), rownames(reference$counts))
features <- setdiff(features, c(state_genes, holdout_markers))
if (length(features) < 10000L) stop("Fewer than 10,000 fit features")
spatial_counts <- spatial$counts[features, , drop = FALSE]
reference_counts <- reference$counts[features, , drop = FALSE]
cell_type <- droplevels(reference$cell_type)
names(cell_type) <- colnames(reference_counts)
if (length(levels(cell_type)) != 7L || any(table(cell_type) < 25L)) stop("Reference ontology/QC mismatch")
coords <- data.frame(x = spatial$metadata$x_um, y = spatial$metadata$y_um,
                     row.names = colnames(spatial_counts))
if (any(!is.finite(as.matrix(coords))) || anyDuplicated(rownames(coords))) stop("Invalid coordinates")
spatial_umi <- setNames(as.numeric(spatial$metadata$total_UMI), colnames(spatial_counts))
reference_umi <- setNames(as.numeric(reference$selected_umi), colnames(reference_counts))
puck <- SpatialRNA(coords = coords, counts = spatial_counts, nUMI = spatial_umi, require_int = TRUE)
ref <- Reference(counts = reference_counts, cell_types = cell_type, nUMI = reference_umi,
                 require_int = TRUE, n_max_cells = ncol(reference_counts), min_UMI = 100)
rctd <- create.RCTD(
  spatialRNA = puck, reference = ref, max_cores = 1, test_mode = FALSE,
  gene_cutoff = 0.000125, fc_cutoff = 0.5, gene_cutoff_reg = 0.0002,
  fc_cutoff_reg = 0.75, UMI_min = 100, UMI_max = 2e7, counts_MIN = 10,
  UMI_min_sigma = 300, CONFIDENCE_THRESHOLD = 5, DOUBLET_THRESHOLD = 20,
  keep_reference = FALSE
)
message("FIT start ", gsm, "; bins=", ncol(spatial_counts), "; input_features=", length(features))
t0 <- proc.time()[3L]
set.seed(42L)
rctd <- spacexr:::fitBulk(rctd)
rctd <- spacexr:::choose_sigma_c(rctd)
model_features <- rctd@internal_vars$gene_list_reg
message("Regression features=", length(model_features))

chunk_size <- 5000L
n_spots <- ncol(rctd@spatialRNA@counts)
barcodes <- colnames(rctd@spatialRNA@counts)
cell_types <- rctd@cell_type_info$renorm[[2L]]
all_weights <- matrix(NA_real_, n_spots, length(cell_types), dimnames = list(barcodes, cell_types))
doublet_weights <- matrix(NA_real_, n_spots, 2L,
                          dimnames = list(barcodes, c("first_type", "second_type")))
parts <- vector("list", ceiling(n_spots / chunk_size))
signature <- paste("G6-DIAG-2026-09-26-001", gsm, paste(path_md5, collapse = ":"),
                   as.character(packageVersion("spacexr")), n_spots, length(features), sep = "|")
checkpoint_dir <- file.path(outdir, "chunk_checkpoints", gsm)
dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
starts <- seq.int(1L, n_spots, by = chunk_size)
for (i in seq_along(starts)) {
  ix <- starts[[i]]:min(n_spots, starts[[i]] + chunk_size - 1L)
  ids <- barcodes[ix]
  checkpoint_path <- file.path(checkpoint_dir, sprintf("chunk_%03d.rds", i))
  if (file.exists(checkpoint_path)) {
    chunk <- readRDS(checkpoint_path)
    if (!identical(chunk$signature, signature) || !identical(chunk$barcodes, ids)) {
      stop("Checkpoint provenance or barcode mismatch: ", checkpoint_path)
    }
    message("RESTORE ", gsm, " chunk ", i, "/", length(starts))
  } else {
    puck_full <- rctd@spatialRNA
    puck_chunk <- SpatialRNA(coords = puck_full@coords[ix, , drop = FALSE],
                             counts = puck_full@counts[, ix, drop = FALSE],
                             nUMI = puck_full@nUMI[ix], require_int = TRUE)
    fit <- rctd
    fit@spatialRNA <- puck_chunk
    fit <- spacexr:::fitPixels(fit, doublet_mode = "doublet")
    chunk <- list(signature = signature, barcodes = ids,
                  all_weights = as.matrix(fit@results$weights),
                  doublet_weights = as.matrix(fit@results$weights_doublet),
                  results_df = fit@results$results_df)
    tmp <- paste0(checkpoint_path, ".tmp")
    saveRDS(chunk, tmp, compress = "gzip")
    if (!file.rename(tmp, checkpoint_path)) stop("Cannot atomically save chunk")
    rm(fit, puck_chunk)
    message("FIT ", gsm, " chunk ", i, "/", length(starts))
  }
  if (!identical(rownames(chunk$all_weights), ids) ||
      !identical(rownames(chunk$doublet_weights), ids) ||
      !identical(rownames(chunk$results_df), ids)) stop("Chunk barcode mismatch")
  all_weights[ix, ] <- chunk$all_weights
  doublet_weights[ix, ] <- chunk$doublet_weights
  parts[[i]] <- chunk$results_df
  rm(chunk)
  gc(verbose = FALSE)
}
results_df <- do.call(rbind, parts)
if (!identical(rownames(results_df), barcodes) || anyNA(all_weights)) stop("Incomplete pilot fit")
rctd@results <- list(results_df = results_df, weights_doublet = doublet_weights)
mapped_doublet <- spacexr::get_doublet_weights(rctd)
if (!identical(rownames(mapped_doublet), barcodes) ||
    !identical(colnames(mapped_doublet), cell_types)) stop("Doublet mapping dimensions mismatch")
out <- list(
  GSM = gsm, donor_id = as.character(spatial$metadata$donor_id[[1L]]),
  all_weights = all_weights, doublet_weights = doublet_weights,
  mapped_doublet_weights = mapped_doublet,
  results_df = results_df,
  metadata = list(id = "G6-DIAG-2026-09-26-001", input_md5 = path_md5,
                  source_sha256 = spatial$construction$source_sha256,
                  R_version = R.version.string, spacexr_version = as.character(packageVersion("spacexr")),
                  seed = 42L, n_bins = n_spots, n_input_features = length(features),
                  n_regression_features = length(model_features), chunk_size = chunk_size,
                  elapsed_seconds = proc.time()[3L] - t0, parameters = rctd@config)
)
tmp <- paste0(out_path, ".tmp")
saveRDS(out, tmp, compress = "gzip")
if (!file.rename(tmp, out_path)) stop("Cannot atomically save final pilot output")
message("FIT complete ", gsm, "; elapsed_s=", round(out$metadata$elapsed_seconds, 1))
