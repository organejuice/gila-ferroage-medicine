#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Matrix)
  library(spacexr)
})
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: 72_fit_g6_rctd_effect_blind.R <project_root> <GSM>")
root <- args[[1]]
gsm <- args[[2]]
run_id <- "2026-09-24-disease-side-v3.5-g6-primary"
lock_id <- "G6-EXEC-2026-09-24-001"
interim <- file.path(root, "data/interim/G6_GSE269906_bin50_tissue")
input_path <- file.path(root, "data/interim/G6_GSE269906_bin50_full",
                        paste0(gsm, "_bin50_centroid_supported_counts.rds"))
outdir <- file.path(root, "data/interim/G6_GSE269906_RCTD")
auditdir <- file.path(root, "results/2026-09-24-disease-side-v3.5-g6-input-audit")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(input_path)) stop("Missing full bin50 RCTD input: ", input_path)

reference <- readRDS(file.path(interim, "G6_RCTD_reference_cells.rds"))
spatial <- readRDS(input_path)
g1 <- read.delim(gzfile(file.path(auditdir, "registries/G1_microglia_signed_sender_244.tsv.gz")),
                 check.names = FALSE, stringsAsFactors = FALSE)
g5 <- read.delim(gzfile(file.path(auditdir, "registries/G5_receiver_one_to_one_human_ortholog_registry.tsv.gz")),
                 check.names = FALSE, stringsAsFactors = FALSE)
state_genes <- unique(toupper(c(g1$gene, g5$human_one_to_one_symbol)))
state_genes <- state_genes[!is.na(state_genes) & nzchar(state_genes)]
holdout_candidates <- c("AIF1", "C1QA", "C1QB", "C1QC", "CX3CR1", "LST1", "P2RY12", "TMEM119", "TREM2", "TYROBP",
                        "ALDH1L1", "GFAP", "GJA1", "GLUL", "S100B", "SLC1A2", "SLC1A3", "SOX9", "CD74", "AQP4")
holdout_markers <- setdiff(holdout_candidates, state_genes)
features <- intersect(rownames(spatial$counts), rownames(reference$counts))
features <- setdiff(features, c(state_genes, holdout_markers))
if (length(features) < 10000L) stop("Fewer than 10,000 RCTD features remain after locked leakage exclusion")

spatial_counts <- spatial$counts[features, , drop = FALSE]
reference_counts <- reference$counts[features, , drop = FALSE]
cell_type <- droplevels(reference$cell_type)
names(cell_type) <- colnames(reference_counts)
if (length(levels(cell_type)) != 7L || any(table(cell_type) < 25L)) {
  stop("Reference does not retain all 7 locked broad cell types with at least 25 nuclei each")
}

coords <- data.frame(x = spatial$metadata$x_um, y = spatial$metadata$y_um,
                     row.names = colnames(spatial_counts))
if (any(!is.finite(as.matrix(coords))) || anyDuplicated(rownames(coords))) stop("Invalid or duplicated spatial coordinates")
spatial_umi <- setNames(as.numeric(spatial$metadata$total_UMI), colnames(spatial_counts))
reference_umi <- setNames(as.numeric(reference$selected_umi), colnames(reference_counts))
puck <- SpatialRNA(coords = coords, counts = spatial_counts, nUMI = spatial_umi,
                   require_int = TRUE)
ref <- Reference(counts = reference_counts, cell_types = cell_type, nUMI = reference_umi,
                 require_int = TRUE, n_max_cells = ncol(reference_counts), min_UMI = 100)

message("RCTD start: ", gsm, "; bins=", ncol(spatial_counts), "; features=", nrow(spatial_counts),
        "; reference nuclei=", ncol(reference_counts), "; reference features=", nrow(reference_counts))
t0 <- proc.time()[3L]
rctd <- create.RCTD(
  spatialRNA = puck,
  reference = ref,
  max_cores = 1,
  test_mode = FALSE,
  gene_cutoff = 0.000125,
  fc_cutoff = 0.5,
  gene_cutoff_reg = 0.0002,
  fc_cutoff_reg = 0.75,
  UMI_min = 100,
  UMI_max = 2e7,
  counts_MIN = 10,
  UMI_min_sigma = 300,
  CONFIDENCE_THRESHOLD = 5,
  DOUBLET_THRESHOLD = 20,
  keep_reference = FALSE
)
model_features <- rctd@internal_vars$gene_list_reg
message("RCTD regression feature count: ", length(model_features))
set.seed(42L)
# fitBulk and sigma selection are section-wide. fitPixels evaluates each bin
# independently, so run it in fixed-size chunks to avoid the package's dense
# all-bin temporary matrix while preserving the same per-bin RCTD algorithm.
rctd <- spacexr:::fitBulk(rctd)
rctd <- spacexr:::choose_sigma_c(rctd)
chunk_size <- 5000L
n_spots <- ncol(rctd@spatialRNA@counts)
source_sha <- spatial$construction$source_sha256
input_md5 <- unname(tools::md5sum(input_path))
reference_path <- file.path(interim, "G6_RCTD_reference_cells.rds")
reference_md5 <- unname(tools::md5sum(reference_path))
script_md5 <- unname(tools::md5sum(file.path(root, "src/72_fit_g6_rctd_effect_blind.R")))
config_paths <- file.path(root, c(
  "configs/g6_v3.5_execution.yaml",
  "configs/g6_v3.5_protocol_amendment_001.yaml",
  "configs/g6_v3.5_execution_addendum_002.yaml",
  "configs/g6_v3.5_execution_addendum_003.yaml",
  "configs/g6_v3.5_execution_addendum_004.yaml",
  "configs/g6_v3.5_execution_addendum_006.yaml"
))
config_md5_values <- unname(tools::md5sum(config_paths))
config_md5 <- paste(config_md5_values, collapse = ":")
config_key <- paste(substr(config_md5_values, 1L, 8L), collapse = "_")
if (anyNA(c(source_sha, input_md5, reference_md5, script_md5, config_md5))) stop("Could not construct a stable RCTD checkpoint signature")
checkpoint_signature <- paste(
  lock_id, gsm, source_sha, input_md5, reference_md5, script_md5, config_md5,
  nrow(spatial_counts), ncol(spatial_counts), chunk_size, as.character(packageVersion("spacexr")),
  sep = "|"
)
checkpoint_key <- paste0(gsm, "_", substr(input_md5, 1L, 8L), "_", substr(reference_md5, 1L, 8L),
                         "_", substr(script_md5, 1L, 8L), "_", config_key)
checkpoint_dir <- file.path(outdir, "chunk_checkpoints", checkpoint_key)
dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
weights <- matrix(0, nrow = n_spots, ncol = length(rctd@cell_type_info$renorm[[2]]),
                  dimnames = list(colnames(rctd@spatialRNA@counts), rctd@cell_type_info$renorm[[2]]))
results_parts <- vector("list", ceiling(n_spots / chunk_size))
chunk_starts <- seq.int(1L, n_spots, by = chunk_size)
for (chunk_i in seq_along(chunk_starts)) {
  ix <- chunk_starts[[chunk_i]]:min(n_spots, chunk_starts[[chunk_i]] + chunk_size - 1L)
  checkpoint_path <- file.path(checkpoint_dir, sprintf("chunk_%03d.rds", chunk_i))
  expected_ids <- colnames(rctd@spatialRNA@counts)[ix]
  if (file.exists(checkpoint_path)) {
    checkpoint <- readRDS(checkpoint_path)
    if (!identical(checkpoint$signature, checkpoint_signature) ||
        !identical(checkpoint$bin_ids, expected_ids)) {
      stop("RCTD checkpoint signature or bin order mismatch: ", checkpoint_path)
    }
    chunk_weights <- checkpoint$weights
    chunk_results <- checkpoint$results_df
    message("RCTD effect-blind fitPixels restored checkpoint ", chunk_i, "/", length(chunk_starts), " for ", gsm)
  } else {
    puck_full <- rctd@spatialRNA
    puck_chunk <- SpatialRNA(
      coords = puck_full@coords[ix, , drop = FALSE],
      counts = puck_full@counts[, ix, drop = FALSE],
      nUMI = puck_full@nUMI[ix], require_int = TRUE
    )
    chunk_rctd <- rctd
    chunk_rctd@spatialRNA <- puck_chunk
    chunk_rctd <- spacexr:::fitPixels(chunk_rctd, doublet_mode = "doublet")
    chunk_weights <- as.matrix(chunk_rctd@results$weights)
    chunk_results <- chunk_rctd@results$results_df
    checkpoint_tmp <- paste0(checkpoint_path, ".tmp")
    saveRDS(list(signature = checkpoint_signature, bin_ids = expected_ids,
                 weights = chunk_weights, results_df = chunk_results), checkpoint_tmp, compress = "gzip")
    if (!file.rename(checkpoint_tmp, checkpoint_path)) stop("Could not atomically save RCTD chunk checkpoint: ", checkpoint_path)
    rm(chunk_rctd, puck_chunk)
  }
  if (!all(rownames(chunk_weights) %in% rownames(weights)) || ncol(chunk_weights) != ncol(weights)) {
    stop("Chunked RCTD returned unexpected weights")
  }
  weights[rownames(chunk_weights), ] <- chunk_weights
  results_parts[[chunk_i]] <- chunk_results
  message("RCTD effect-blind fitPixels chunk ", chunk_i, "/", length(chunk_starts),
          " completed for ", gsm)
  rm(chunk_weights, chunk_results)
  gc(verbose = FALSE)
}
results_df <- do.call(rbind, results_parts)
results_df <- results_df[rownames(weights), , drop = FALSE]
if (!identical(rownames(weights), colnames(spatial_counts)) || ncol(weights) != 7L) {
  stop("RCTD returned unexpected weight dimensions or barcodes")
}
if (!identical(rownames(results_df), rownames(weights))) stop("RCTD classification table does not align with weights")
if (any(!is.finite(weights)) || any(weights < -1e-8)) stop("RCTD returned non-finite or negative cell-type weights")
weight_sum_raw <- rowSums(weights)
if (any(!is.finite(weight_sum_raw)) || any(weight_sum_raw < 0)) stop("RCTD returned invalid total weight")
usable <- results_df$spot_class != "reject" & weight_sum_raw > 0
if (!any(usable)) stop("RCTD rejected every spatial bin")
weights_normalized <- matrix(0, nrow = nrow(weights), ncol = ncol(weights),
                             dimnames = dimnames(weights))
weights_normalized[usable, ] <- weights[usable, , drop = FALSE] / weight_sum_raw[usable]
weight_sum_normalized <- rowSums(weights_normalized)
weight_sum_summary <- data.frame(
  GSM = gsm, donor_id = spatial$metadata$donor_id[[1L]],
  n_bins = nrow(weights), n_nonfinite = sum(!is.finite(weights)),
  n_negative = sum(weights < -1e-8),
  n_rejected = sum(results_df$spot_class == "reject"),
  n_nonrejected_zero_weight = sum(results_df$spot_class != "reject" & weight_sum_raw == 0),
  n_nonrejected_weight_sum_outside_098_102 = sum(usable & (weight_sum_raw < 0.98 | weight_sum_raw > 1.02)),
  raw_weight_sum_min = min(weight_sum_raw), raw_weight_sum_median = median(weight_sum_raw),
  raw_weight_sum_max = max(weight_sum_raw),
  normalized_sum_min_usable = min(weight_sum_normalized[usable]),
  normalized_sum_max_usable = max(weight_sum_normalized[usable]),
  stringsAsFactors = FALSE
)
write.csv(weight_sum_summary, file.path(outdir, paste0(gsm, "_RCTD_weight_sum_qc.csv")), row.names = FALSE)

# Holdout markers are excluded from the fit and used only for identity QC.
marker_rows <- list()
marker_sets <- list(
  MG = intersect(c("AIF1", "C1QA", "C1QB", "C1QC", "CX3CR1", "LST1", "P2RY12", "TMEM119", "TREM2", "TYROBP"), holdout_markers),
  ASC = intersect(c("ALDH1L1", "GFAP", "GJA1", "GLUL", "S100B", "SLC1A2", "SLC1A3", "SOX9", "CD74", "AQP4"), holdout_markers)
)
for (lineage in names(marker_sets)) {
  genes <- intersect(marker_sets[[lineage]], rownames(spatial$counts))
  if (length(genes) >= 3L) {
    marker_logcpm <- log2(t(t(as.matrix(spatial$counts[genes, , drop = FALSE])) /
                               pmax(spatial$metadata$total_UMI, 1)) * 1e6 + 1)
    marker_score <- colMeans(marker_logcpm)
    target_col <- if (lineage == "MG") "MG" else "ASC"
    marker_rows[[lineage]] <- data.frame(
      GSM = gsm, donor_id = spatial$metadata$donor_id[[1L]], lineage = lineage,
      n_holdout_markers = length(genes), markers = paste(genes, collapse = ";"),
      spearman_rho = suppressWarnings(cor(marker_score[usable], weights_normalized[usable, target_col], method = "spearman")),
      n_bins = sum(usable), stringsAsFactors = FALSE
    )
  } else {
    marker_rows[[lineage]] <- data.frame(
      GSM = gsm, donor_id = spatial$metadata$donor_id[[1L]], lineage = lineage,
      n_holdout_markers = length(genes), markers = paste(genes, collapse = ";"),
      spearman_rho = NA_real_, n_bins = ncol(weights), stringsAsFactors = FALSE
    )
  }
}
marker_qc <- do.call(rbind, marker_rows)

out <- data.frame(
  spatial$metadata,
  weights_normalized,
  raw_weight_sum = weight_sum_raw,
  spot_class = as.character(results_df$spot_class),
  first_type = as.character(results_df$first_type),
  second_type = as.character(results_df$second_type),
  first_class = results_df$first_class,
  second_class = results_df$second_class,
  conv_all = results_df$conv_all,
  conv_doublet = results_df$conv_doublet,
  weight_sum_normalized = weight_sum_normalized,
  row.names = NULL,
  check.names = FALSE
)
saveRDS(list(
  weights = out,
  weights_raw = weights,
  marker_qc = marker_qc,
  weight_sum_qc = weight_sum_summary,
  metadata = list(
    run_id = run_id, execution_lock_id = lock_id, accession = "GSE269906",
    execution_addendum_ids = c("G6-EXEC-2026-09-25-002", "G6-EXEC-2026-09-25-003", "G6-EXEC-2026-09-25-004", "G6-EXEC-2026-09-25-006"),
    GSM = gsm, donor_id = spatial$metadata$donor_id[[1L]],
    reference = "GSE174367 human cortical snRNA-seq",
    reference_nuclei = ncol(reference_counts),
    reference_donors = length(unique(reference$donor_id)),
    reference_cell_type_counts = as.list(table(cell_type)),
    n_model_features_reg = length(model_features),
    n_spatial_bins = ncol(spatial_counts),
    n_RCTD_rejected = sum(results_df$spot_class == "reject"),
    n_RCTD_nonrejected_zero_weight = sum(results_df$spot_class != "reject" & weight_sum_raw == 0),
    rctd_version = as.character(packageVersion("spacexr")),
    rctd_mode = "doublet", max_cores = 1L,
    fitpixels_chunk_size_bins = chunk_size,
    chunk_checkpointing = TRUE,
    checkpoint_signature = checkpoint_signature,
    checkpoint_directory = checkpoint_dir,
    parameters = rctd@config,
    elapsed_seconds = proc.time()[3L] - t0,
    R_version = R.version.string
  )
), file.path(outdir, paste0(gsm, "_RCTD_weights_effect_blind.rds")), compress = "gzip")
write.csv(marker_qc, file.path(outdir, paste0(gsm, "_RCTD_holdout_marker_qc.csv")), row.names = FALSE)
message("RCTD complete: ", gsm, "; rejected=", sum(results_df$spot_class == "reject"),
        "; elapsed_s=", round(proc.time()[3L] - t0, 1))
