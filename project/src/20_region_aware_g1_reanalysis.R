#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
root_arg <- if (length(args) >= 1) args[[1]] else "."
root <- if (root_arg %in% c("", ".")) "." else normalizePath(root_arg, winslash = "/", mustWork = TRUE)
run_id <- if (length(args) >= 2) args[[2]] else "2026-09-20-disease-side-v3-g1-region-aware"
out <- file.path(root, "results", run_id)
dir.create(out, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out, "effects"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out, "gates"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out, "meta"), recursive = TRUE, showWarnings = FALSE)

write_tsv <- function(x, path) {
  write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}

safe_num <- function(x) {
  x <- as.numeric(x)
  x[!is.finite(x)] <- NA_real_
  x
}

fit_region <- function(counts, meta, cohort, region) {
  m <- meta[meta$region == region, , drop = FALSE]
  m$disease <- factor(m$disease, levels = c("Control", "AD"))
  m$sex <- factor(m$sex)
  m$age_c <- as.numeric(m$age) - mean(as.numeric(m$age), na.rm = TRUE)
  m$rin_c <- as.numeric(m$rin) - mean(as.numeric(m$rin), na.rm = TRUE)
  n_ad <- sum(m$disease == "AD")
  n_ctrl <- sum(m$disease == "Control")
  if (n_ad < 2 || n_ctrl < 2) return(NULL)

  dat <- counts[, m$group_id, drop = FALSE]
  keep <- Matrix::rowSums(dat >= 10) >= 2
  dat <- dat[keep, , drop = FALSE]
  design <- model.matrix(~ disease + age_c + sex + rin_c, data = m)
  y <- DGEList(counts = round(as.matrix(dat)), samples = m)
  y <- calcNormFactors(y)
  v <- voom(y, design, plot = FALSE)
  corfit <- tryCatch(duplicateCorrelation(v, design = design, block = m$donor_id), error = function(e) NULL)
  if (!is.null(corfit) && is.finite(corfit$consensus)) {
    fit <- lmFit(v, design, block = m$donor_id, correlation = corfit$consensus)
  } else {
    fit <- lmFit(v, design)
  }
  coef_name <- grep("^diseaseAD$", colnames(design), value = TRUE)
  if (length(coef_name) != 1) stop(sprintf("AD coefficient not identifiable for %s", region))
  fit <- eBayes(fit)
  se <- sqrt(fit$s2.post) * fit$stdev.unscaled[, coef_name]
  data.frame(
    cohort = cohort,
    region = region,
    cell_type = "Microglia",
    gene = rownames(fit$coefficients),
    n_groups = nrow(m),
    n_AD_donors = n_ad,
    n_Control_donors = n_ctrl,
    logFC = safe_num(fit$coefficients[, coef_name]),
    SE = safe_num(se),
    t = safe_num(fit$t[, coef_name]),
    P.Value = safe_num(fit$p.value[, coef_name]),
    FDR_within_region = safe_num(p.adjust(fit$p.value[, coef_name], method = "BH")),
    model = "limma_voom_donor_region_pseudobulk_age_sex_RIN",
    stringsAsFactors = FALSE
  )
}

# D3: existing raw-derived donor x region x microglia pseudobulk.
d3_counts_path <- file.path(root, "data", "interim", "D3_pseudobulk_counts_microglia.tsv.gz")
d3_meta_path <- file.path(root, "data", "interim", "D3_pseudobulk_metadata_microglia.tsv")
d3_counts <- read.delim(gzfile(d3_counts_path), check.names = FALSE, stringsAsFactors = FALSE)
d3_genes <- d3_counts[[1]]
d3_counts[[1]] <- NULL
rownames(d3_counts) <- make.unique(as.character(d3_genes))
d3_meta <- read.delim(d3_meta_path, check.names = FALSE, stringsAsFactors = FALSE)
stopifnot(all(d3_meta$group_id %in% colnames(d3_counts)))
d3_counts <- d3_counts[, d3_meta$group_id, drop = FALSE]
rownames(d3_counts) <- toupper(rownames(d3_counts))
d3_counts <- d3_counts[!duplicated(rownames(d3_counts)), , drop = FALSE]
d3_results <- do.call(rbind, lapply(sort(unique(d3_meta$region)), function(r) fit_region(d3_counts, d3_meta, "D3_GSE160936", r)))
write_tsv(d3_results, file.path(out, "effects", "D3_microglia_region_gene_effects.tsv"))

# D6: audit whether a region/structure-specific AD-Control donor contrast is estimable.
d6_meta_path <- file.path(root, "data", "raw", "D6_GSE214979", "GSE214979_cell_metadata.csv.gz")
d6 <- read.csv(gzfile(d6_meta_path), check.names = FALSE, stringsAsFactors = FALSE)
d6 <- d6[d6$predicted.id == "Microglia", , drop = FALSE]
d6$donor_id <- d6$id
group_counts <- as.data.frame(table(d6$Status, d6$structure, d6$donor_id), stringsAsFactors = FALSE)
colnames(group_counts) <- c("Status", "region", "donor_id", "n_cells")
group_counts <- group_counts[group_counts$n_cells > 0, , drop = FALSE]
keep_groups <- group_counts[group_counts$n_cells >= 20, , drop = FALSE]
regions <- sort(unique(d6$structure))
d6_region_rows <- lapply(regions, function(r) {
  m <- keep_groups[keep_groups$region == r, , drop = FALSE]
  ad <- unique(m$donor_id[m$Status == "AD"])
  ctrl <- unique(m$donor_id[m$Status == "Ctrl"])
  status_present <- length(ad) > 0 && length(ctrl) > 0
  stable_replication <- length(ad) >= 2 && length(ctrl) >= 2
  reason <- if (!status_present) {
    "Only one disease status is present after the >=20-cell donor-region filter"
  } else if (!stable_replication) {
    "Both statuses are present, but at least one status has fewer than two donor replicates"
  } else {
    "Estimable candidate"
  }
  data.frame(
    cohort = "D6_GSE214979",
    region = r,
    cell_type = "Microglia",
    n_AD_donors = length(ad),
    n_Control_donors = length(ctrl),
    n_AD_cells = sum(d6$Status == "AD" & d6$structure == r),
    n_Control_cells = sum(d6$Status == "Ctrl" & d6$structure == r),
    status_present = status_present,
    stable_replication_for_region_model = stable_replication,
    estimable_for_frozen_G1 = stable_replication,
    reason = reason,
    model = "donor_region_feasibility_audit_before_expression_fit",
    stringsAsFactors = FALSE
  )
})
d6_region <- do.call(rbind, d6_region_rows)
write_tsv(d6_region, file.path(out, "effects", "D6_microglia_region_feasibility.tsv"))

# Explicitly record non-estimable D6 effects; no cross-region substitution is allowed.
d6_effects <- d6_region[, c("cohort", "region", "cell_type", "n_AD_donors", "n_Control_donors", "estimable_for_frozen_G1", "reason")]
d6_effects$gene <- NA_character_
d6_effects$logFC <- NA_real_
d6_effects$SE <- NA_real_
d6_effects$P.Value <- NA_real_
d6_effects$FDR_within_region <- NA_real_
d6_effects <- d6_effects[, c("cohort", "region", "cell_type", "gene", "n_AD_donors", "n_Control_donors", "logFC", "SE", "P.Value", "FDR_within_region", "estimable_for_frozen_G1", "reason")]
write_tsv(d6_effects, file.path(out, "effects", "D6_microglia_region_gene_effects.tsv"))

# Frozen meta gate: require at least two independent cohorts with comparable region-level effects.
meta_summary <- data.frame(
  endpoint = "microglial_gene_state",
  candidate_cohorts = 2,
  cohorts_with_region_effects = 1,
  regions_with_two_statuses_and_stable_replication = sum(d6_region$estimable_for_frozen_G1),
  random_effects_meta_status = "NOT_ESTIMABLE",
  BH_FDR_status = "NOT_ESTIMABLE",
  reason = "D6 has no structure with >=2 AD and >=2 Control microglial donor replicates after the predeclared >=20-cell group filter; a one-cohort or cross-structure substitute would violate the frozen region-aware design.",
  stringsAsFactors = FALSE
)
write_tsv(meta_summary, file.path(out, "meta", "random_effects_meta_summary.tsv"))

loco <- data.frame(
  analysis = "leave_one_cohort_out",
  status = "NOT_RUN",
  reason = "No valid two-cohort region-aware meta input exists; removing either cohort would leave zero or one estimable cohort.",
  stringsAsFactors = FALSE
)
write_tsv(loco, file.path(out, "meta", "leave_one_cohort_out.tsv"))

sensitivity <- data.frame(
  sensitivity = c("sex", "age", "APOE", "platform", "brain_region", "cell_proportion", "influential_dataset"),
  status = c("D3_region_model_included", "D3_region_model_included", "NOT_AVAILABLE_CROSS_COHORT", "NOT_ESTIMABLE_CROSS_COHORT", "D3_COMPLETED_D6_NOT_ESTIMABLE", "NOT_AVAILABLE_COMPARABLE", "NOT_RUN_NO_META"),
  reason = c("D3 region-specific primary model includes sex", "D3 region-specific primary model includes centered age", "APOE is not available in a comparable D3/D6 region-level effect table", "D3 and D6 platforms/cohort effects cannot be compared after D6 regional non-estimability", "D3 EC/SSC effects fitted; D6 has no stable within-structure contrast", "No comparable donor-level cell-proportion adjustment table was produced", "Influential-cohort analysis requires valid two-cohort meta inputs"),
  stringsAsFactors = FALSE
)
write_tsv(sensitivity, file.path(out, "meta", "sensitivity_audit.tsv"))

# Existing full-cohort concordance is retained as a separate descriptive check.
old_concordance <- file.path(root, "results", "archive", "2026-09-20-pre-v3-freeze", "results", "tables", "D6_signature_concordance.tsv")
rho <- NA_real_
if (file.exists(old_concordance)) {
  cc <- read.delim(old_concordance, stringsAsFactors = FALSE)
  rho <- cc$spearman_rho[cc$cell_type == "Microglia"][1]
}

criteria <- data.frame(
  criterion = c("two_independent_human_cohorts", "donor_level_microglial_model", "brain_regions_modeled_separately", "concordant_direction_in_two_cohorts", "random_effects_meta_BH_FDR_lt_0_05", "leave_one_cohort_out_stable", "heterogeneity_and_covariate_sensitivity"),
  status = c("PASS", "PASS", "FAIL", "FAIL", "FAIL", "FAIL", "FAIL"),
  observed = c("D3 and D6 registered", "D3 donor-region pseudobulk; D6 donor-structure metadata", "D3 estimable; D6 no stable region contrast", sprintf("Existing descriptive D3-D6 Microglia Spearman rho=%s", rho), "NOT_ESTIMABLE", "NOT_RUN", "D3-only partial sensitivity; cross-cohort sensitivity not estimable"),
  reason = c("Two candidate cohorts are available", "Cells are not treated as independent donors", "D6 status is confounded with structure after donor/cell filter", "Existing genome-wide concordance is near zero and not a region-aware replication test", "No two-cohort region-level effect set", "No valid meta input", "No valid two-cohort meta input"),
  stringsAsFactors = FALSE
)
write_tsv(criteria, file.path(out, "gates", "G1_reanalysis_criteria.tsv"))

final <- list(
  run_id = run_id,
  protocol = "AD disease-side v3.0",
  decision = "STOPPED_AT_G1",
  d3_region_effects = nrow(d3_results),
  d3_regions = sort(unique(d3_results$region)),
  d6_region_rows = nrow(d6_region),
  d6_estimable_regions = sum(d6_region$estimable_for_frozen_G1),
  d6_microglia_full_cohort_spearman_rho = rho,
  meta_status = "NOT_ESTIMABLE",
  loco_status = "NOT_RUN",
  downstream = list(G2 = "NOT_RUN_UPSTREAM_FAIL", G3 = "NOT_RUN_UPSTREAM_FAIL", G4 = "NOT_RUN_UPSTREAM_FAIL", G5 = "NOT_RUN_UPSTREAM_FAIL", G6 = "NOT_RUN_UPSTREAM_FAIL"),
  claim_boundary = "No downstream mechanism, spatial communication, drug ranking, or efficacy claim is made.",
  created_at = format(Sys.time(), tz = "UTC", usetz = TRUE)
)
write_json(final, file.path(out, "final_decision.json"), auto_unbox = TRUE, pretty = TRUE)

report <- c(
  paste0("# G1 region-aware reanalysis: ", run_id),
  "",
  "The D3 region-level donor model was recomputed from the existing raw-derived pseudobulk counts and exported with logFC, SE, t, P.Value, and within-region FDR.",
  "",
  "D6 was audited at the donor-by-structure level before expression fitting. After the existing >=20-cell donor-structure filter, BA46 contains AD donors but only one Control donor, while BA9 and DLPFC contain Control donors without an AD counterpart. Therefore no D6 structure supports a stable within-structure AD-Control estimate with the required donor replication.",
  "",
  "Because a valid two-cohort region-aware input set does not exist, random-effects meta-analysis, leave-one-cohort-out analysis, and cross-cohort sensitivity analysis are reported as not estimable/not run. The frozen G1 gate remains failed and G2-G6 remain locked.",
  "",
  "No raw data were modified."
)
writeLines(report, file.path(out, "execution_report.md"), useBytes = TRUE)
cat(jsonlite::toJSON(final, auto_unbox = TRUE, pretty = TRUE), "\n")
