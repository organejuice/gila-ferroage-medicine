#!/usr/bin/env Rscript

# G6-FP-2026-09-27-001, WP1. Diagnosis-blind RCTD nUMI sensitivity.
suppressPackageStartupMessages({ library(Matrix); library(spacexr) })
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("Usage: 89_fit_g6_numi_sensitivity.R <project_root> <GSM> <full|fit_features>")
root <- args[[1L]]
gsm <- args[[2L]]
mode <- args[[3L]]
if (!gsm %in% c("GSM8330063", "GSM8330064")) stop("GSM is outside the locked donor pair")
if (!mode %in% c("full", "fit_features")) stop("Unknown nUMI mode")
if (as.character(packageVersion("spacexr")) != "2.2.1") stop("spacexr version mismatch")

interim <- file.path(root, "data/interim/G6_GSE269906_bin50_tissue")
auditdir <- file.path(root, "results/2026-09-24-disease-side-v3.5-g6-input-audit/registries")
input_path <- file.path(root, "data/interim/G6_GSE269906_bin50_full", paste0(gsm, "_bin50_centroid_supported_counts.rds"))
ref_path <- file.path(interim, "G6_RCTD_reference_cells.rds")
g1_path <- file.path(auditdir, "G1_microglia_signed_sender_244.tsv.gz")
g5_path <- file.path(auditdir, "G5_receiver_one_to_one_human_ortholog_registry.tsv.gz")
config_path <- file.path(root, "configs/g6_next_blind_feasibility_001.yaml")
script_path <- file.path(root, "src/89_fit_g6_numi_sensitivity.R")
paths <- c(input_path, ref_path, g1_path, g5_path, config_path, script_path)
if (!all(file.exists(paths))) stop("A locked input or config is missing: ", paste(paths[!file.exists(paths)], collapse = "; "))
input_md5 <- unname(tools::md5sum(paths))
names(input_md5) <- paths

outdir <- file.path(root, "data/interim/G6_next_feasibility_001/numi")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
out_path <- file.path(outdir, paste0(gsm, "_", mode, ".rds"))
if (file.exists(out_path)) stop("Existing final output; refusing overwrite: ", out_path)

reference <- readRDS(ref_path)
spatial <- readRDS(input_path)
metadata_names <- tolower(names(spatial$metadata))
forbidden <- grepl("diagnos|patholog|\bad\b|control|braak|cerad|thal", metadata_names)
if (any(forbidden)) stop("Diagnosis/pathology metadata entered WP1 input: ", paste(names(spatial$metadata)[forbidden], collapse = ","))
g1 <- read.delim(gzfile(g1_path), check.names = FALSE, stringsAsFactors = FALSE)
g5 <- read.delim(gzfile(g5_path), check.names = FALSE, stringsAsFactors = FALSE)
state_genes <- unique(toupper(c(g1$gene, g5$human_one_to_one_symbol)))
state_genes <- state_genes[!is.na(state_genes) & nzchar(state_genes)]
holdout_candidates <- c("AIF1", "C1QA", "C1QB", "C1QC", "CX3CR1", "LST1", "P2RY12", "TMEM119", "TREM2", "TYROBP",
                        "ALDH1L1", "GFAP", "GJA1", "GLUL", "S100B", "SLC1A2", "SLC1A3", "SOX9", "CD74", "AQP4")
holdout_markers <- setdiff(holdout_candidates, state_genes)
features <- intersect(rownames(spatial$counts), rownames(reference$counts))
features <- setdiff(features, c(state_genes, holdout_markers))
if (length(features) < 10000L) stop("Fewer than 10,000 eligible fit features")
spatial_counts <- spatial$counts[features, , drop = FALSE]
reference_counts <- reference$counts[features, , drop = FALSE]
cell_type <- droplevels(reference$cell_type)
names(cell_type) <- colnames(reference_counts)
if (length(levels(cell_type)) != 7L || any(table(cell_type) < 25L)) stop("Reference ontology/QC mismatch")
coords <- data.frame(x = spatial$metadata$x_um, y = spatial$metadata$y_um, row.names = colnames(spatial_counts))
if (any(!is.finite(as.matrix(coords))) || anyDuplicated(rownames(coords))) stop("Invalid spatial coordinates")
n_umi <- if (mode == "full") {
  setNames(as.numeric(spatial$metadata$total_UMI), colnames(spatial_counts))
} else {
  setNames(as.numeric(Matrix::colSums(spatial_counts)), colnames(spatial_counts))
}
if (any(!is.finite(n_umi)) || any(n_umi <= 0)) stop("Invalid nUMI vector")
puck <- SpatialRNA(coords = coords, counts = spatial_counts, nUMI = n_umi, require_int = TRUE)
reference_umi <- setNames(as.numeric(reference$selected_umi), colnames(reference_counts))
ref <- Reference(counts = reference_counts, cell_types = cell_type, nUMI = reference_umi,
                 require_int = TRUE, n_max_cells = ncol(reference_counts), min_UMI = 100)
rctd <- create.RCTD(
  spatialRNA = puck, reference = ref, max_cores = 1, test_mode = FALSE,
  gene_cutoff = 0.000125, fc_cutoff = 0.5, gene_cutoff_reg = 0.0002,
  fc_cutoff_reg = 0.75, UMI_min = 100, UMI_max = 2e7, counts_MIN = 10,
  UMI_min_sigma = 300, CONFIDENCE_THRESHOLD = 5, DOUBLET_THRESHOLD = 20,
  keep_reference = FALSE
)
message("FIT start ", gsm, "; mode=", mode, "; bins=", ncol(spatial_counts), "; features=", length(features))
t0 <- proc.time()[3L]
set.seed(42L)
rctd <- spacexr:::fitBulk(rctd)
rctd <- spacexr:::choose_sigma_c(rctd)
model_features <- rctd@internal_vars$gene_list_reg
chunk_size <- 5000L
n_spots <- ncol(rctd@spatialRNA@counts)
barcodes <- colnames(rctd@spatialRNA@counts)
cell_types <- rctd@cell_type_info$renorm[[2L]]
all_weights <- matrix(NA_real_, n_spots, length(cell_types), dimnames = list(barcodes, cell_types))
doublet_weights <- matrix(NA_real_, n_spots, 2L, dimnames = list(barcodes, c("first_type", "second_type")))
parts <- vector("list", ceiling(n_spots / chunk_size))
signature <- paste("G6-FP-2026-09-27-001", gsm, mode, paste(input_md5, collapse = ":"),
                   as.character(packageVersion("spacexr")), n_spots, length(features),
                   paste(format(n_umi, scientific = FALSE), collapse = ","), sep = "|")
checkpoint_dir <- file.path(outdir, "chunk_checkpoints", gsm, mode)
dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
starts <- seq.int(1L, n_spots, by = chunk_size)
for (i in seq_along(starts)) {
  ix <- starts[[i]]:min(n_spots, starts[[i]] + chunk_size - 1L)
  ids <- barcodes[ix]
  checkpoint_path <- file.path(checkpoint_dir, sprintf("chunk_%03d.rds", i))
  if (file.exists(checkpoint_path)) {
    chunk <- readRDS(checkpoint_path)
    if (!identical(chunk$signature, signature) || !identical(chunk$barcodes, ids)) stop("Checkpoint signature/barcode mismatch")
    message("RESTORE ", gsm, " ", mode, " chunk ", i, "/", length(starts))
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
    message("FIT ", gsm, " ", mode, " chunk ", i, "/", length(starts))
  }
  if (!identical(rownames(chunk$all_weights), ids) || !identical(rownames(chunk$doublet_weights), ids) ||
      !identical(rownames(chunk$results_df), ids)) stop("Chunk barcode order mismatch")
  all_weights[ix, ] <- chunk$all_weights
  doublet_weights[ix, ] <- chunk$doublet_weights
  parts[[i]] <- chunk$results_df
  rm(chunk)
  gc(verbose = FALSE)
}
results_df <- do.call(rbind, parts)
if (!identical(rownames(results_df), barcodes) || anyNA(all_weights)) stop("Incomplete fit")
rctd@results <- list(results_df = results_df, weights_doublet = doublet_weights)
mapped_doublet <- spacexr::get_doublet_weights(rctd)
if (!identical(rownames(mapped_doublet), barcodes) || !identical(colnames(mapped_doublet), cell_types)) stop("Doublet map mismatch")
out <- list(
  study_id = "G6-FP-2026-09-27-001", GSM = gsm, mode = mode,
  donor_id = as.character(spatial$metadata$donor_id[[1L]]),
  all_weights = all_weights, doublet_weights = doublet_weights,
  mapped_doublet_weights = mapped_doublet, results_df = results_df,
  nUMI = n_umi,
  metadata = list(input_md5 = input_md5, fit_features = features,
                  source_sha256 = spatial$construction$source_sha256,
                  R_version = R.version.string, spacexr_version = as.character(packageVersion("spacexr")),
                  seed = 42L, n_bins = n_spots, n_input_features = length(features),
                  n_regression_features = length(model_features), chunk_size = chunk_size,
                  elapsed_seconds = proc.time()[3L] - t0, parameters = rctd@config,
                  reference_status = "REFERENCE_IDENTITY_UNRESOLVED_DIAGNOSTIC_ONLY")
)
tmp <- paste0(out_path, ".tmp")
saveRDS(out, tmp, compress = "gzip")
if (!file.rename(tmp, out_path)) stop("Cannot atomically save final output")
writeLines(capture.output(sessionInfo()), file.path(outdir, paste0(gsm, "_", mode, "_sessionInfo.txt")))
message("FIT complete ", gsm, "; mode=", mode, "; elapsed_s=", round(out$metadata$elapsed_seconds, 1))
