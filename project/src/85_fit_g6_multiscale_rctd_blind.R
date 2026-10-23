#!/usr/bin/env Rscript
# Diagnostic-blind full-mode RCTD for the bin50/bin110 scale experiment.
suppressPackageStartupMessages({ library(Matrix); library(spacexr) })
if (as.character(packageVersion("spacexr")) != "2.2.1") stop("Expected audited spacexr 2.2.1")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("Usage: 85_fit_g6_multiscale_rctd_blind.R <root> <GSM> <50|110>")
root <- args[[1L]]; gsm <- args[[2L]]; width <- as.integer(args[[3L]])
if (is.na(width) || !width %in% c(50L, 110L) || !grepl("^GSM[0-9]+$", gsm)) stop("Invalid GSM or width")
base <- file.path(root, "data/interim/G6_GSE269906_multiscale_001", paste0("bin", width))
input <- file.path(base, "full", paste0(gsm, "_bin", width, "_centroid_supported_counts.rds"))
reference_path <- file.path(root, "data/interim/G6_GSE269906_multiscale_001/reference/G6_RCTD_reference_no_candidate_overlap.rds")
outdir <- file.path(base, "rctd_full_mode")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(input) || !file.exists(reference_path)) stop("Missing multiscale counts or audited reference")
spatial <- readRDS(input); reference <- readRDS(reference_path)
audit <- file.path(root, "results/2026-09-24-disease-side-v3.5-g6-input-audit/registries")
g1 <- read.delim(gzfile(file.path(audit, "G1_microglia_signed_sender_244.tsv.gz")))
g5 <- read.delim(gzfile(file.path(audit, "G5_receiver_one_to_one_human_ortholog_registry.tsv.gz")))
states <- unique(toupper(c(g1$gene, g5$human_one_to_one_symbol)))
states <- states[!is.na(states) & nzchar(states)]
markers <- c("AIF1", "C1QA", "C1QB", "C1QC", "CX3CR1", "LST1", "P2RY12", "TMEM119", "TREM2", "TYROBP",
             "ALDH1L1", "GFAP", "GJA1", "GLUL", "S100B", "SLC1A2", "SLC1A3", "SOX9", "CD74", "AQP4")
holdout <- setdiff(markers, states)
features <- setdiff(intersect(rownames(spatial$counts), rownames(reference$counts)), c(states, holdout))
if (length(features) < 10000L) stop("Insufficient independent reference features")
counts <- spatial$counts[features, , drop = FALSE]
refcounts <- reference$counts[features, , drop = FALSE]
celltype <- droplevels(reference$cell_type); names(celltype) <- colnames(refcounts)
if (length(levels(celltype)) != 7L || any(table(celltype) < 25L)) stop("Reference cell-type audit failed")
coords <- data.frame(x = spatial$metadata$x_um, y = spatial$metadata$y_um,
                     row.names = colnames(counts))
if (any(!is.finite(as.matrix(coords))) || anyDuplicated(rownames(coords))) stop("Invalid coordinates")
puck <- SpatialRNA(coords, counts, nUMI = setNames(as.numeric(spatial$metadata$total_UMI), colnames(counts)), require_int = TRUE)
ref <- Reference(refcounts, cell_types = celltype,
                 nUMI = setNames(as.numeric(reference$selected_umi), colnames(refcounts)),
                 require_int = TRUE, n_max_cells = ncol(refcounts), min_UMI = 100)
rctd <- create.RCTD(puck, ref, max_cores = 1, test_mode = FALSE,
                    gene_cutoff = 0.000125, fc_cutoff = 0.5,
                    gene_cutoff_reg = 0.0002, fc_cutoff_reg = 0.75,
                    UMI_min = 100, UMI_max = 2e7, counts_MIN = 10,
                    UMI_min_sigma = 300, CONFIDENCE_THRESHOLD = 5,
                    DOUBLET_THRESHOLD = 20, keep_reference = FALSE)
set.seed(42L)
rctd <- spacexr:::fitBulk(rctd)
rctd <- spacexr:::choose_sigma_c(rctd)
input_md5 <- unname(tools::md5sum(input)); ref_md5 <- unname(tools::md5sum(reference_path))
script_md5 <- unname(tools::md5sum(file.path(root, "src/85_fit_g6_multiscale_rctd_blind.R")))
signature <- paste("G6-MS-2026-09-26-001", gsm, width, input_md5, ref_md5, script_md5,
                   as.character(packageVersion("spacexr")), sep = "|")
checkpoint_dir <- file.path(outdir, paste0(gsm, "_", substr(input_md5, 1, 8), "_", substr(script_md5, 1, 8)))
dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
starts <- seq.int(1L, ncol(counts), by = 2000L)
parts <- vector("list", length(starts))
for (i in seq_along(starts)) {
  ix <- starts[[i]]:min(ncol(counts), starts[[i]] + 1999L)
  ids <- colnames(counts)[ix]
  path <- file.path(checkpoint_dir, sprintf("chunk_%04d.rds", i))
  if (file.exists(path)) {
    part <- readRDS(path)
    if (!identical(part$signature, signature) || !identical(part$ids, ids)) stop("Checkpoint mismatch")
  } else {
    chunk <- rctd
    chunk@spatialRNA <- SpatialRNA(rctd@spatialRNA@coords[ix, , drop = FALSE],
                                  rctd@spatialRNA@counts[, ix, drop = FALSE],
                                  nUMI = rctd@spatialRNA@nUMI[ix], require_int = TRUE)
    chunk <- spacexr:::fitPixels(chunk, doublet_mode = "full")
    part <- list(signature = signature, ids = ids, weights = as.matrix(chunk@results$weights),
                 results_df = chunk@results$results_df)
    tmp <- paste0(path, ".tmp")
    saveRDS(part, tmp, compress = "gzip")
    if (!file.rename(tmp, path)) stop("Could not save checkpoint")
    rm(chunk)
  }
  if (!identical(rownames(part$weights), ids)) stop("RCTD bin order mismatch")
  parts[[i]] <- part
  message(gsm, " bin", width, " RCTD chunk ", i, "/", length(starts))
}
weights <- do.call(rbind, lapply(parts, `[[`, "weights"))
if (!identical(rownames(weights), colnames(counts)) || any(!is.finite(weights)) || any(weights < -1e-8)) stop("Invalid RCTD weights")
totals <- rowSums(weights)
valid <- totals > 0
weights[valid, ] <- weights[valid, , drop = FALSE] / totals[valid]
weights[!valid, ] <- 0
if (!all(c("MG", "ASC") %in% colnames(weights))) stop("MG/ASC missing from reference fit")
marker_score <- function(genes) {
  genes <- intersect(setdiff(genes, states), rownames(spatial$counts))
  if (length(genes) < 3L) return(rep(NA_real_, nrow(weights)))
  m <- as.matrix(spatial$counts[genes, , drop = FALSE])
  colMeans(log2(t(t(m) / pmax(spatial$metadata$total_UMI, 1)) * 1e6 + 1))
}
mg_score <- marker_score(markers[1:10]); asc_score <- marker_score(markers[11:20])
safe_cor <- function(x, y) if (sum(is.finite(x) & is.finite(y)) < 10L || sd(x, na.rm = TRUE) == 0 || sd(y, na.rm = TRUE) == 0) NA_real_ else suppressWarnings(cor(x, y, method = "spearman", use = "complete.obs"))
qc <- data.frame(GSM = gsm, donor_id = spatial$metadata$donor_id[[1L]], width_DNB = width,
                 n_bins = nrow(weights), n_valid_weights = sum(valid),
                 median_MG_weight = median(weights[valid, "MG"]),
                 n_MG_ge_0_5_diagnostic = sum(valid & weights[, "MG"] >= 0.5),
                 n_ASC_ge_0_5_diagnostic = sum(valid & weights[, "ASC"] >= 0.5),
                 n_disjoint_MG_holdout_markers = sum(markers[1:10] %in% rownames(spatial$counts)),
                 n_disjoint_ASC_holdout_markers = sum(markers[11:20] %in% rownames(spatial$counts)),
                 MG_marker_rho = safe_cor(mg_score[valid], weights[valid, "MG"]),
                 ASC_marker_rho = safe_cor(asc_score[valid], weights[valid, "ASC"]),
                 input_md5 = input_md5, reference_md5 = ref_md5, script_md5 = script_md5)
saveRDS(list(weights = weights, metadata = spatial$metadata, qc = qc, method = "RCTD full mode; exploratory scale comparison; diagnosis blinded",
             signature = signature), file.path(outdir, paste0(gsm, "_bin", width, "_rctd_full.rds")), compress = "gzip")
write.csv(qc, file.path(outdir, paste0(gsm, "_bin", width, "_rctd_full_qc.csv")), row.names = FALSE)
