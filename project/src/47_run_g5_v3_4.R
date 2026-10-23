#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, width = 220)
suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
  library(yaml)
  library(jsonlite)
})

root_override <- Sys.getenv("G5_ROOT_ASCII", unset = "")
root <- if (nzchar(root_override)) root_override else normalizePath(".", winslash = "/", mustWork = TRUE)
cfg <- yaml::read_yaml(file.path(root, "configs", "g5_v3.4_execution.yaml"))
run_dir <- file.path(root, cfg$outputs$result_dir)
dir.create(file.path(run_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(run_dir, "gates"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(run_dir, "registries"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(run_dir, "logs"), recursive = TRUE, showWarnings = FALSE)

audit_path <- file.path(run_dir, cfg$outputs$audit_json)
if (!file.exists(audit_path)) stop("G0 input audit is missing: ", audit_path)
audit <- jsonlite::fromJSON(audit_path)
if (!isTRUE(audit$primary_input_pass)) stop("G0 v3.4 input gate did not pass; refusing to fit the primary model")

write_tsv <- function(x, path) {
  write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}

as_numeric_matrix <- function(x) {
  m <- as.matrix(x)
  storage.mode(m) <- "numeric"
  m
}

extract_symbol <- function(ids, fallback = NULL) {
  ids <- as.character(ids)
  out <- toupper(ids)
  for (i in seq_along(ids)) {
    pieces <- strsplit(ids[[i]], "_", fixed = TRUE)[[1]]
    if (length(pieces) >= 2 && grepl("^ENSMUSG", toupper(pieces[[1]]))) out[[i]] <- toupper(pieces[[2]])
    if (!is.null(fallback) && toupper(pieces[[1]]) %in% names(fallback)) out[[i]] <- toupper(fallback[[toupper(pieces[[1]])]])
  }
  out[out %in% c("", "NA", "N/A", "NULL")] <- NA_character_
  out
}

sum_duplicate_symbols <- function(mat, symbols) {
  ok <- !is.na(symbols) & nzchar(symbols)
  mat <- mat[ok, , drop = FALSE]
  symbols <- symbols[ok]
  if (!nrow(mat)) stop("No gene symbols remain after identifier normalization")
  out <- rowsum(mat, group = symbols, reorder = FALSE)
  rownames(out) <- toupper(rownames(out))
  out
}

read_primary_csv <- function(path) {
  frame <- read.csv(gzfile(path), check.names = FALSE, stringsAsFactors = FALSE)
  id_col <- match("gene_id", names(frame))
  if (is.na(id_col)) stop("gene_id column missing from ", path)
  sample_names <- setdiff(names(frame), names(frame)[id_col])
  list(ids = as.character(frame[[id_col]]), samples = sample_names, values = as_numeric_matrix(frame[sample_names]))
}

read_count_text <- function(path) {
  con <- gzfile(path, open = "rt")
  lines <- readLines(con, warn = FALSE)
  close(con)
  tokens <- strsplit(lines, "\t", fixed = TRUE)
  header <- gsub("^\"|\"$", "", tokens[[1]])
  if (length(header) >= 1L && tolower(header[[1]]) == "gene_id") header <- header[-1L]
  if (!length(header) || any(!nzchar(header))) stop("Invalid count matrix header in ", path)
  width <- length(header) + 1L
  rows <- lapply(tokens[-1], function(x) {
    if (length(x) < width) x <- c(x, rep(NA_character_, width - length(x)))
    x[seq_len(width)]
  })
  ids <- vapply(rows, function(x) gsub("^\"|\"$", "", x[[1]]), character(1))
  values <- do.call(rbind, lapply(rows, function(x) as.numeric(gsub("^\"|\"$", "", x[-1]))))
  colnames(values) <- header
  list(ids = ids, samples = colnames(values), values = values)
}

parse_soft_samples <- function(path) {
  lines <- readLines(gzfile(path), warn = FALSE, encoding = "UTF-8")
  starts <- grep("^\\^SAMPLE", lines)
  if (!length(starts)) stop("No sample records in ", path)
  ends <- c(starts[-1] - 1L, length(lines))
  rows <- vector("list", length(starts))
  get_one <- function(block, key) {
    hit <- grep(paste0("^!Sample_", key, " = "), block, value = TRUE)
    if (!length(hit)) return(NA_character_)
    sub(paste0("^!Sample_", key, " = "), "", hit[[1]])
  }
  get_many <- function(block, key) {
    hit <- grep(paste0("^!Sample_", key, " = "), block, value = TRUE)
    sub(paste0("^!Sample_", key, " = "), "", hit)
  }
  for (i in seq_along(starts)) {
    block <- lines[starts[[i]]:ends[[i]]]
    chars <- get_many(block, "characteristics_ch1")
    rows[[i]] <- data.frame(
      geo_accession = get_one(block, "geo_accession"),
      title = get_one(block, "title"),
      source_name = get_one(block, "source_name_ch1"),
      characteristics = paste(chars, collapse = " | "),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

primary_meta <- function(samples) {
  pattern <- cfg$inputs$primary$sample_id_pattern
  hits <- regexec(pattern, samples)
  parsed <- regmatches(samples, hits)
  rows <- lapply(seq_along(samples), function(i) {
    x <- parsed[[i]]
    if (length(x) != 6L) return(data.frame(sample_id = samples[[i]], metadata_parse_pass = FALSE))
    source_group <- x[[2]]
    numeric_id <- x[[3]]
    time_label <- x[[4]]
    sex <- x[[5]]
    batch <- paste0("batch", x[[6]])
    time_values <- c(`12am` = 0, `2am` = 2, `4am` = 4, `6am` = 6, `8am` = 8, `10am` = 10,
                     `12pm` = 12, `2pm` = 14, `4pm` = 16, `6pm` = 18, `8pm` = 20, `10pm` = 22)
    hour <- unname(time_values[[time_label]])
    data.frame(
      sample_id = samples[[i]], metadata_parse_pass = TRUE, source_group = source_group,
      disease = ifelse(source_group == "APP", "APP_PS1", "WT"), animal_id = paste0(source_group, numeric_id),
      clock_label = time_label, clock_hour = hour, sin_time = sin(2 * pi * hour / 24),
      cos_time = cos(2 * pi * hour / 24), time_bin = ifelse(hour < 12, "AM", "PM"),
      sex = sex, batch = batch, stringsAsFactors = FALSE
    )
  })
  meta <- do.call(rbind, rows)
  meta$disease <- factor(meta$disease, levels = c("WT", "APP_PS1"))
  meta$sex <- factor(meta$sex)
  meta$batch <- factor(meta$batch)
  meta
}

soft_meta_to_astro <- function(soft, count_samples, dataset) {
  idx <- grepl("astro", paste(soft$title, soft$source_name, soft$characteristics), ignore.case = TRUE)
  soft <- soft[idx, , drop = FALSE]
  if (nrow(soft) != length(count_samples)) stop(dataset, ": astrocyte SOFT records do not match count columns")
  soft$sample_id <- count_samples
  soft$disease <- ifelse(grepl("APP", paste(soft$title, soft$source_name, soft$characteristics), ignore.case = TRUE), "APP_PS1", "WT")
  soft$sex <- sub(".*Sex: *([^|]+).*", "\\1", soft$characteristics)
  soft$sex[!grepl("Sex:", soft$characteristics)] <- NA_character_
  soft$disease <- factor(soft$disease, levels = c("WT", "APP_PS1"))
  soft$sex <- factor(trimws(soft$sex))
  soft
}

read_m7 <- function() {
  path <- file.path(root, cfg$inputs$m7_reference)
  ref <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  unique(toupper(as.character(ref$gene_symbol[ref$module == "M7"])))
}

modules <- cfg$analysis$modules
if ("M7_FerroAging_Author" %in% names(modules)) {
  m7_value <- modules[["M7_FerroAging_Author"]]
  if (length(m7_value) == 1L && identical(as.character(m7_value), "reference:m7_reference")) modules[["M7_FerroAging_Author"]] <- read_m7()
} else {
  modules[["M7_FerroAging_Author"]] <- read_m7()
}
modules <- lapply(modules, function(x) unique(toupper(as.character(x))))
if (length(modules) != 7L) stop("Expected exactly seven frozen G5 modules, found ", length(modules))

fit_score_table <- function(scores, meta, formula_text, dataset, analysis, coefficient = "diseaseAPP_PS1",
                            block = NULL, note = "", normalization = "primary_voom") {
  meta <- meta[match(colnames(scores), meta$sample_id), , drop = FALSE]
  if (anyNA(meta$sample_id)) stop("metadata/count mismatch for ", dataset, " / ", analysis)
  design <- model.matrix(as.formula(formula_text), data = meta)
  if (qr(design)$rank < ncol(design) || !(coefficient %in% colnames(design))) {
    stop("rank deficient or missing coefficient in ", dataset, " / ", analysis, " / ", formula_text)
  }
  if (is.null(block)) {
    fit <- eBayes(lmFit(scores, design))
    cor_value <- NA_real_
  } else {
    cor_fit <- duplicateCorrelation(scores, design, block = block)
    fit <- eBayes(lmFit(scores, design, block = block, correlation = cor_fit$consensus))
    cor_value <- unname(cor_fit$consensus)
  }
  coef_idx <- match(coefficient, colnames(design))
  rows <- lapply(seq_len(nrow(scores)), function(i) {
    se <- sqrt(fit$s2.post[[i]]) * fit$stdev.unscaled[i, coef_idx]
    beta <- fit$coefficients[i, coef_idx]
    df <- fit$df.total[[i]]
    data.frame(
      dataset = dataset, analysis = analysis, module = rownames(scores)[[i]],
      formula = formula_text, coefficient = coefficient,
      n_present = NA_integer_, n_total = NA_integer_, coverage = NA_real_,
      n_units = nrow(meta), n_case = sum(meta$disease == "APP_PS1", na.rm = TRUE),
      n_control = sum(meta$disease == "WT", na.rm = TRUE),
      logFC = beta, SE = se, t = fit$t[i, coef_idx], P.Value = fit$p.value[i, coef_idx],
      df_total = df, CI_low = beta - qt(0.975, df) * se, CI_high = beta + qt(0.975, df) * se,
      status = "estimable", normalization = normalization, block_correlation = cor_value,
      note = note, stringsAsFactors = FALSE
    )
  })
  list(effects = do.call(rbind, rows), design = design, fit = fit, meta = meta)
}

fit_counts <- function(counts, meta, formula_text, dataset, analysis, coefficient = "diseaseAPP_PS1",
                       block = NULL, note = "", normalization = "edgeR_TMM_then_limma_voom") {
  meta <- meta[match(colnames(counts), meta$sample_id), , drop = FALSE]
  if (anyNA(meta$sample_id)) stop("metadata/count mismatch for ", dataset, " / ", analysis)
  design <- model.matrix(as.formula(formula_text), data = meta)
  if (qr(design)$rank < ncol(design) || !(coefficient %in% colnames(design))) {
    stop("rank deficient or missing coefficient in ", dataset, " / ", analysis, " / ", formula_text)
  }
  dge <- DGEList(counts = counts)
  keep <- filterByExpr(dge, design = design)
  if (sum(keep) < 2L) stop("fewer than two genes remain after filterByExpr in ", dataset, " / ", analysis)
  dge <- calcNormFactors(dge[keep, , keep.lib.sizes = FALSE], method = "TMM")
  v <- voom(dge, design, plot = FALSE)
  block_correlation <- NA_real_
  if (is.null(block)) {
    fit <- eBayes(lmFit(v, design))
  } else {
    cor_fit <- duplicateCorrelation(v$E, design, block = block)
    block_correlation <- unname(cor_fit$consensus)
    v <- voom(dge, design, plot = FALSE, block = block, correlation = block_correlation)
    fit <- eBayes(lmFit(v, design, block = block, correlation = block_correlation))
  }
  coef_idx <- match(coefficient, colnames(design))
  rows <- vector("list", length(modules))
  scores <- matrix(NA_real_, nrow = length(modules), ncol = ncol(v$E), dimnames = list(names(modules), colnames(v$E)))
  for (i in seq_along(modules)) {
    module <- names(modules)[[i]]
    genes <- intersect(modules[[i]], rownames(v$E))
    scores[i, ] <- if (length(genes) >= 2L) colMeans(v$E[genes, , drop = FALSE]) else NA_real_
    if (length(genes) < 2L) {
      rows[[i]] <- data.frame(
        dataset = dataset, analysis = analysis, module = module, formula = formula_text,
        coefficient = coefficient, n_present = length(genes), n_total = length(modules[[i]]),
        coverage = length(genes) / length(modules[[i]]), n_units = nrow(meta),
        n_case = sum(meta$disease == "APP_PS1", na.rm = TRUE), n_control = sum(meta$disease == "WT", na.rm = TRUE),
        logFC = NA_real_, SE = NA_real_, t = NA_real_, P.Value = NA_real_, df_total = NA_real_,
        CI_low = NA_real_, CI_high = NA_real_, status = "non_estimable", normalization = normalization,
        block_correlation = block_correlation, note = paste(note, "fewer than two present genes"), stringsAsFactors = FALSE
      )
    } else {
      score_fit <- if (is.null(block)) eBayes(lmFit(matrix(scores[i, ], nrow = 1, dimnames = list(module, colnames(v$E))), design)) else eBayes(lmFit(matrix(scores[i, ], nrow = 1, dimnames = list(module, colnames(v$E))), design, block = block, correlation = block_correlation))
      beta <- score_fit$coefficients[1, coef_idx]
      se <- sqrt(score_fit$s2.post[[1]]) * score_fit$stdev.unscaled[1, coef_idx]
      df <- score_fit$df.total[[1]]
      rows[[i]] <- data.frame(
        dataset = dataset, analysis = analysis, module = module, formula = formula_text,
        coefficient = coefficient, n_present = length(genes), n_total = length(modules[[i]]),
        coverage = length(genes) / length(modules[[i]]), n_units = nrow(meta),
        n_case = sum(meta$disease == "APP_PS1", na.rm = TRUE), n_control = sum(meta$disease == "WT", na.rm = TRUE),
        logFC = beta, SE = se, t = score_fit$t[1, coef_idx], P.Value = score_fit$p.value[1, coef_idx],
        df_total = df, CI_low = beta - qt(0.975, df) * se, CI_high = beta + qt(0.975, df) * se,
        status = "estimable", normalization = normalization, block_correlation = block_correlation,
        note = note, stringsAsFactors = FALSE
      )
    }
  }
  list(effects = do.call(rbind, rows), scores = scores, voom = v, design = design, meta = meta, keep = keep)
}

primary_path_app <- file.path(root, cfg$inputs$primary$app_counts)
primary_path_wt <- file.path(root, cfg$inputs$primary$wt_counts)
app <- read_primary_csv(primary_path_app)
wt <- read_primary_csv(primary_path_wt)
primary_map <- setNames(extract_symbol(app$ids), toupper(sub("_.*$", "", app$ids)))
app_symbols <- extract_symbol(app$ids, primary_map)
wt_symbols <- extract_symbol(wt$ids, primary_map)
app_counts <- sum_duplicate_symbols(app$values, app_symbols)
wt_counts <- sum_duplicate_symbols(wt$values, wt_symbols)
common <- intersect(rownames(app_counts), rownames(wt_counts))
primary_counts <- cbind(app_counts[common, , drop = FALSE], wt_counts[common, , drop = FALSE])
meta_primary <- primary_meta(colnames(primary_counts))
meta_primary$sample_id <- as.character(meta_primary$sample_id)
meta_primary$animal_id <- as.character(meta_primary$animal_id)
write_tsv(meta_primary, file.path(run_dir, "tables/GSE261698_model_metadata.tsv"))
write_tsv(data.frame(raw_id = app$ids, gene_symbol = app_symbols, stringsAsFactors = FALSE), file.path(run_dir, "registries/GSE261698_gene_mapping.tsv"))

primary_formula <- "~ disease + sin_time + cos_time + batch + sex"
primary <- fit_counts(primary_counts, meta_primary, primary_formula, "GSE261698", "primary", note = "frozen v3.4 primary animal model")
primary$effects$FDR_BH_7 <- p.adjust(primary$effects$P.Value, method = "BH", n = 7L)
write_tsv(primary$effects, file.path(run_dir, "tables/GSE261698_primary_module_effects.tsv"))
write_tsv(data.frame(gene = rownames(primary$voom$E), filtered_in = TRUE, stringsAsFactors = FALSE), file.path(run_dir, "tables/GSE261698_filter_genes.tsv"))
write_tsv(as.data.frame(primary$scores), file.path(run_dir, "tables/GSE261698_voom_module_scores.tsv"))

run_sensitivity <- function(label, subset_idx, formula_text, note = "") {
  m <- meta_primary[subset_idx, , drop = FALSE]
  cts <- primary_counts[, subset_idx, drop = FALSE]
  ans <- tryCatch(fit_counts(cts, m, formula_text, "GSE261698", label, note = note)$effects,
                  error = function(e) data.frame(dataset = "GSE261698", analysis = label, module = names(modules), formula = formula_text,
                    coefficient = "diseaseAPP_PS1", n_present = NA, n_total = vapply(modules, length, integer(1)), coverage = NA,
                    n_units = nrow(m), n_case = sum(m$disease == "APP_PS1", na.rm = TRUE), n_control = sum(m$disease == "WT", na.rm = TRUE),
                    logFC = NA, SE = NA, t = NA, P.Value = NA, df_total = NA, CI_low = NA, CI_high = NA, status = "not_estimable",
                    normalization = "edgeR_TMM_then_limma_voom", block_correlation = NA, note = paste(note, conditionMessage(e)), stringsAsFactors = FALSE))
  ans$FDR_BH_7 <- p.adjust(ans$P.Value, method = "BH", n = 7L)
  ans
}

sens <- list()
sens[[length(sens) + 1L]] <- run_sensitivity("AM", meta_primary$time_bin == "AM", primary_formula, "AM stratum")
sens[[length(sens) + 1L]] <- run_sensitivity("PM", meta_primary$time_bin == "PM", primary_formula, "PM stratum")
for (tp in sort(unique(meta_primary$clock_hour))) {
  sens[[length(sens) + 1L]] <- run_sensitivity(paste0("leave_one_time_point_out_", tp), meta_primary$clock_hour != tp, primary_formula, paste0("excluded clock hour ", tp))
}
for (batch in levels(meta_primary$batch)) {
  keep <- meta_primary$batch != batch
  formula_batch <- "~ disease + sin_time + cos_time + sex"
  sens[[length(sens) + 1L]] <- run_sensitivity(paste0("leave_one_batch_out_", batch), keep, formula_batch, paste0("excluded ", batch, "; batch term omitted because it is constant"))
}
for (sex in levels(meta_primary$sex)) {
  keep <- meta_primary$sex == sex
  formula_sex <- "~ disease + sin_time + cos_time + batch"
  sens[[length(sens) + 1L]] <- run_sensitivity(paste0("sex_stratified_", sex), keep, formula_sex, paste0("sex stratum ", sex, "; sex term omitted because it is constant"))
}
sens_df <- do.call(rbind, sens)
write_tsv(sens_df, file.path(run_dir, "tables/GSE261698_primary_sensitivities.tsv"))

message("stage: disease_by_circular_time")
interaction_rows <- lapply(seq_len(nrow(primary$scores)), function(i) {
  module <- rownames(primary$scores)[[i]]
  if (any(!is.finite(primary$scores[i, ]))) {
    return(data.frame(module = module, joint_F = NA_real_, joint_P = NA_real_, status = "not_estimable", stringsAsFactors = FALSE))
  }
  meta <- meta_primary
  design <- model.matrix(~ disease * sin_time + disease * cos_time + batch + sex, data = meta)
  fit <- eBayes(lmFit(matrix(primary$scores[i, ], nrow = 1, dimnames = list(module, colnames(primary$scores))), design))
  idx <- grep("^diseaseAPP_PS1:(sin_time|cos_time)$", colnames(design))
  if (length(idx) != 2L) return(data.frame(module = module, joint_F = NA, joint_P = NA, status = "not_estimable", stringsAsFactors = FALSE))
  b <- as.numeric(fit$coefficients[1, idx])
  covb <- fit$cov.coefficients[idx, idx, drop = FALSE] * fit$s2.post[[1]]
  fval <- as.numeric(sum(b * solve(covb, b)) / length(idx))
  data.frame(module = module, joint_F = fval, joint_P = pf(fval, length(idx), fit$df.total[[1]], lower.tail = FALSE), status = "estimable", stringsAsFactors = FALSE)
})
write_tsv(do.call(rbind, interaction_rows), file.path(run_dir, "tables/GSE261698_disease_by_circular_time.tsv"))

message("stage: rounded_matrix_sensitivity")
round_counts <- round(primary_counts)
rounded <- tryCatch(fit_counts(round_counts, meta_primary, primary_formula, "GSE261698", "rounded_matrix_sensitivity", note = "non-primary rounded matrix; cannot change the gate", normalization = "rounded_counts_TMM_then_limma_voom")$effects,
                    error = function(e) data.frame(module = names(modules), status = "not_estimable", note = conditionMessage(e), stringsAsFactors = FALSE))
if ("P.Value" %in% names(rounded)) rounded$FDR_BH_7 <- p.adjust(rounded$P.Value, method = "BH", n = 7L)
write_tsv(rounded, file.path(run_dir, "tables/GSE261698_rounded_matrix_sensitivity.tsv"))

message("stage: leave_one_animal_out")
influence_rows <- list(); k <- 0L
for (animal in unique(meta_primary$animal_id)) {
  keep <- meta_primary$animal_id != animal
  ans <- tryCatch(fit_score_table(primary$scores[, keep, drop = FALSE], meta_primary[keep, , drop = FALSE], primary_formula,
                                  "GSE261698", "leave_one_animal_out", normalization = "primary_voom_reused", note = paste0("excluded animal ", animal))$effects,
                  error = function(e) data.frame(module = names(modules), status = "not_estimable", note = conditionMessage(e), stringsAsFactors = FALSE))
  ans$excluded_animal <- animal
  k <- k + 1L; influence_rows[[k]] <- ans
}
write_tsv(do.call(rbind, influence_rows), file.path(run_dir, "tables/GSE261698_leave_one_animal_out.tsv"))

message("stage: GSE291239")
g291 <- read_count_text(file.path(root, cfg$inputs$validations$GSE291239$counts))
soft291 <- parse_soft_samples(file.path(root, cfg$inputs$validations$GSE291239$family_soft))
meta291 <- soft_meta_to_astro(soft291, g291$samples, "GSE291239")
g291_counts <- sum_duplicate_symbols(g291$values, extract_symbol(g291$ids, primary_map))
formula291 <- if (length(unique(na.omit(meta291$sex))) > 1L && qr(model.matrix(~ disease + sex, meta291))$rank == 3L) "~ disease + sex" else "~ disease"
val291 <- fit_counts(g291_counts, meta291, formula291, "GSE291239", "independent_animal_validation", note = "independent AppNL-G-F astrocyte animal validation")$effects
val291$FDR_BH_7 <- p.adjust(val291$P.Value, method = "BH", n = 7L)
write_tsv(meta291, file.path(run_dir, "tables/GSE291239_metadata_used.tsv"))
write_tsv(val291, file.path(run_dir, "tables/GSE291239_module_effects.tsv"))

read_series_matrix <- function(path) {
  lines <- readLines(gzfile(path), warn = FALSE, encoding = "UTF-8")
  start <- grep("^!series_matrix_table_begin", lines)
  end <- grep("^!series_matrix_table_end", lines)
  if (length(start) != 1L || length(end) != 1L) stop("series matrix table markers missing")
  tab <- read.delim(text = paste(lines[(start + 1L):(end - 1L)], collapse = "\n"), header = TRUE, sep = "\t", quote = "\"", check.names = FALSE, comment.char = "", stringsAsFactors = FALSE)
  ids <- as.character(tab[[1]])
  mat <- as_numeric_matrix(tab[, -1, drop = FALSE])
  colnames(mat) <- gsub("^\"|\"$", "", names(tab)[-1])
  list(ids = ids, samples = colnames(mat), values = mat)
}

parse_platform_table <- function(path) {
  lines <- readLines(gzfile(path), warn = FALSE, encoding = "UTF-8")
  start <- grep("^!platform_table_begin", lines)
  end <- grep("^!platform_table_end", lines)
  if (length(start) != 1L || length(end) != 1L) stop("platform table markers missing")
  read.delim(text = paste(lines[(start + 1L):(end - 1L)], collapse = "\n"), header = TRUE, sep = "\t", quote = "", fill = TRUE, check.names = FALSE, comment.char = "", stringsAsFactors = FALSE)
}

message("stage: GSE74614")
g746 <- read_series_matrix(file.path(root, cfg$inputs$validations$GSE74614$series_matrix))
soft746 <- parse_soft_samples(file.path(root, cfg$inputs$validations$GSE74614$family_soft))
soft746$sample_id <- soft746$geo_accession
desc746 <- paste(soft746$title, soft746$source_name, soft746$characteristics)
sel746 <- grepl("astrocytes", desc746, ignore.case = TRUE) & grepl("GFAP\\+", desc746) & grepl("Vim\\+", desc746)
sel746_meta <- soft746[sel746, , drop = FALSE]
sel746_meta$disease <- factor(ifelse(grepl("APPswe", paste(sel746_meta$title, sel746_meta$source_name, sel746_meta$characteristics), ignore.case = TRUE), "APP_PS1", "WT"), levels = c("WT", "APP_PS1"))
sel746_meta$sample_id <- as.character(sel746_meta$sample_id)
if (nrow(sel746_meta) != 8L) stop("GSE74614 GFAP+/+ Vim+/+ astrocyte selection yielded ", nrow(sel746_meta), " samples, expected 8")
mat746 <- g746$values[, sel746_meta$sample_id, drop = FALSE]
platform746 <- parse_platform_table(file.path(root, cfg$inputs$validations$GSE74614$family_soft))
platform746$ID <- as.character(platform746$ID)
platform746$GENE_SYMBOL <- toupper(trimws(as.character(platform746$GENE_SYMBOL)))
probe_idx <- match(g746$ids, platform746$ID)
gene_symbols746 <- platform746$GENE_SYMBOL[probe_idx]
micro_scores <- matrix(NA_real_, nrow = length(modules), ncol = ncol(mat746), dimnames = list(names(modules), colnames(mat746)))
micro_rows <- list()
for (i in seq_along(modules)) {
  module <- names(modules)[[i]]
  genes <- intersect(modules[[i]], unique(gene_symbols746))
  probe_rows <- which(gene_symbols746 %in% genes)
  gene_means <- do.call(rbind, lapply(genes, function(g) colMeans(mat746[probe_rows[gene_symbols746[probe_rows] == g], , drop = FALSE], na.rm = TRUE)))
  if (length(genes) >= 2L) micro_scores[i, ] <- colMeans(gene_means) else micro_scores[i, ] <- NA_real_
  micro_rows[[i]] <- data.frame(dataset = "GSE74614", module = module, n_present = length(genes), n_total = length(modules[[i]]), coverage = length(genes) / length(modules[[i]]), n_probes = length(probe_rows), stringsAsFactors = FALSE)
}
micro <- fit_score_table(micro_scores, sel746_meta, "~ disease", "GSE74614", "orthogonal_microarray_validation", normalization = "source_processed_microarray_values", note = "GFAP+/+ Vim+/+ cortical astrocytes only")$effects
micro <- merge(micro, do.call(rbind, micro_rows), by = c("dataset", "module"), all.x = TRUE, suffixes = c("", "_map"))
micro$FDR_BH_7 <- p.adjust(micro$P.Value, method = "BH", n = 7L)
write_tsv(sel746_meta, file.path(run_dir, "tables/GSE74614_metadata_used.tsv"))
write_tsv(data.frame(probe_id = g746$ids, gene_symbol = gene_symbols746, stringsAsFactors = FALSE), file.path(run_dir, "registries/GSE74614_probe_gene_map.tsv"))
write_tsv(micro, file.path(run_dir, "tables/GSE74614_module_effects.tsv"))

message("stage: GSE290101")
g290 <- read_count_text(file.path(root, cfg$inputs$validations$GSE290101$counts))
soft290 <- parse_soft_samples(file.path(root, cfg$inputs$validations$GSE290101$family_soft))
if (nrow(soft290) != length(g290$samples)) stop("GSE290101 SOFT records do not match count columns")
meta290 <- soft290
meta290$sample_id <- g290$samples
meta290$disease <- factor(ifelse(grepl("APP/PS1", paste(meta290$title, meta290$characteristics), ignore.case = TRUE), "APP_PS1", "WT"), levels = c("WT", "APP_PS1"))
meta290$subpopulation <- factor(sub(".*subpopulation: *([^|]+).*", "\\1", meta290$characteristics))
meta290$animal_id <- sub("^((AD|WT)[0-9]+)-.*", "\\1", meta290$title)
g290_counts <- sum_duplicate_symbols(g290$values, extract_symbol(g290$ids, primary_map))
val290 <- fit_counts(g290_counts, meta290, "~ disease + subpopulation", "GSE290101", "animal_aware_repeated_measure_sensitivity", block = meta290$animal_id, note = "animal-blocked duplicateCorrelation across reporter-defined astrocyte subpopulations")$effects
val290$FDR_BH_7 <- p.adjust(val290$P.Value, method = "BH", n = 7L)
write_tsv(meta290, file.path(run_dir, "tables/GSE290101_metadata_used.tsv"))
write_tsv(val290, file.path(run_dir, "tables/GSE290101_module_effects.tsv"))

cat(jsonlite::toJSON(list(
  run_id = cfg$run_id, primary_rows = nrow(primary$effects), sensitivity_rows = nrow(sens_df),
  GSE291239_rows = nrow(val291), GSE74614_rows = nrow(micro), GSE290101_rows = nrow(val290),
  primary_modules = names(modules), primary_gene_count_after_symbol_sum = nrow(primary_counts)
), auto_unbox = TRUE, pretty = TRUE), "\n")
