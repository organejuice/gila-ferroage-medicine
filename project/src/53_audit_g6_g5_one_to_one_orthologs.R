#!/usr/bin/env Rscript

# Resolve the frozen v3.4 G5 receiver module symbols through mouse gene IDs
# and a versioned Ensembl one-to-one orthology release. No G6 expression is read.

suppressPackageStartupMessages(library(AnnotationDbi))
suppressPackageStartupMessages(library(org.Mm.eg.db))
suppressPackageStartupMessages(library(biomaRt))
args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) args[[1L]] else "."
reg_dir <- file.path(root, "results/2026-09-24-disease-side-v3.5-g6-input-audit/registries")
input_path <- file.path(reg_dir, "G5_frozen_receiver_module_gene_registry.tsv.gz")
model_map_path <- file.path(root, "results/2026-09-22-disease-side-v3.4-g5/registries/GSE261698_gene_mapping.tsv")
if (!file.exists(input_path)) stop("Run src/51_freeze_g6_upstream_signatures.py first.")

modules <- read.delim(gzfile(input_path), check.names = FALSE, stringsAsFactors = FALSE)
model_map <- read.delim(model_map_path, check.names = FALSE, stringsAsFactors = FALSE)
model_map$mouse_ensembl_id <- sub("_.*$", "", model_map$raw_id)

# Link frozen module symbols to observed mouse Ensembl features by the upstream
# v3.4 mapping. Use the installed mouse OrgDb alias map only when no exact match
# exists; ambiguity is retained rather than silently selecting one gene.
mouse_options <- lapply(modules$gene_symbol_in_frozen_G5_list, function(symbol) {
  exact <- unique(model_map$mouse_ensembl_id[toupper(model_map$gene_symbol) == toupper(symbol)])
  exact <- exact[!is.na(exact) & nzchar(exact)]
  if (length(exact)) return(data.frame(mouse_ensembl_id = exact, resolution = "GSE261698 exact normalized symbol", stringsAsFactors = FALSE))
  alias_rows <- tryCatch(
    AnnotationDbi::select(org.Mm.eg.db, keys = symbol, keytype = "ALIAS", columns = c("SYMBOL", "ENSEMBL")),
    error = function(e) data.frame(ALIAS = character(), SYMBOL = character(), ENSEMBL = character())
  )
  alias_rows <- alias_rows[!is.na(alias_rows$ENSEMBL) & nzchar(alias_rows$ENSEMBL), , drop = FALSE]
  if (!nrow(alias_rows)) return(data.frame(mouse_ensembl_id = "", resolution = "no mouse ID mapping", stringsAsFactors = FALSE))
  data.frame(mouse_ensembl_id = unique(alias_rows$ENSEMBL), resolution = "org.Mm.eg.db ALIAS", stringsAsFactors = FALSE)
})
mouse_map <- do.call(rbind, lapply(seq_len(nrow(modules)), function(i) {
  opt <- mouse_options[[i]]
  data.frame(
    input_row = i,
    module = modules$module[i],
    gene_symbol_in_frozen_G5_list = modules$gene_symbol_in_frozen_G5_list[i],
    mouse_ensembl_id = opt$mouse_ensembl_id,
    mouse_symbol_resolution = opt$resolution,
    n_mouse_ids_for_input = length(opt$mouse_ensembl_id),
    mouse_mapping_status = if (!length(opt$mouse_ensembl_id) || !nzchar(opt$mouse_ensembl_id[1L])) "no_mouse_gene_mapping" else if (length(opt$mouse_ensembl_id) == 1L) "unique_mouse_gene_mapping" else "ambiguous_mouse_gene_mapping",
    stringsAsFactors = FALSE
  )
}))

eligible_ids <- unique(mouse_map$mouse_ensembl_id[mouse_map$mouse_mapping_status == "unique_mouse_gene_mapping" & nzchar(mouse_map$mouse_ensembl_id)])
attempt <- list(
  date = as.character(Sys.Date()),
  orthology_source = "Ensembl BioMart archived release 113; mouse-to-human gene orthology attributes",
  biomaRt_version = as.character(packageVersion("biomaRt")),
  AnnotationDbi_version = as.character(packageVersion("AnnotationDbi")),
  org.Mm.eg.db_version = as.character(packageVersion("org.Mm.eg.db")),
  n_frozen_module_gene_rows = nrow(modules),
  n_mouse_ensembl_queries = length(eligible_ids),
  status = "pending"
)

mart <- tryCatch(
  {
    Sys.setenv(BIOMART_CACHE = file.path(tempdir(), "biomaRt_cache"))
    Sys.unsetenv(c("ALL_PROXY", "all_proxy", "HTTP_PROXY", "http_proxy", "HTTPS_PROXY", "https_proxy"))
    options(timeout = 30)
    biomaRt::useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl", version = 113)
  },
  error = function(e) {
    attempt$status <<- "network_or_archive_error"
    attempt$error <<- conditionMessage(e)
    jsonlite::write_json(attempt, file.path(reg_dir, "G6_orthology_query_status.json"), auto_unbox = TRUE, pretty = TRUE)
    stop("Unable to connect to Ensembl release 113: ", conditionMessage(e))
  }
)
available <- biomaRt::listAttributes(mart)$name
wanted <- c(
  "ensembl_gene_id", "external_gene_name", "hsapiens_homolog_ensembl_gene",
  "hsapiens_homolog_associated_gene_name", "hsapiens_homolog_orthology_type",
  "hsapiens_homolog_orthology_confidence", "hsapiens_homolog_perc_id"
)
if (!all(wanted %in% available)) stop("Ensembl release 113 lacks expected orthology attributes: ", paste(setdiff(wanted, available), collapse = ","))
bm <- if (length(eligible_ids)) biomaRt::getBM(
  attributes = wanted,
  filters = "ensembl_gene_id",
  values = eligible_ids,
  mart = mart
) else data.frame()

one2one <- bm[
  bm$hsapiens_homolog_orthology_type == "ortholog_one2one" &
    bm$hsapiens_homolog_orthology_confidence == 1 &
    !is.na(bm$hsapiens_homolog_ensembl_gene) & nzchar(bm$hsapiens_homolog_ensembl_gene),
  , drop = FALSE
]
one2one_n <- table(one2one$ensembl_gene_id)
accepted_ids <- names(one2one_n)[one2one_n == 1L]
one2one <- one2one[one2one$ensembl_gene_id %in% accepted_ids, , drop = FALSE]

joined <- merge(mouse_map, bm, by.x = "mouse_ensembl_id", by.y = "ensembl_gene_id", all.x = TRUE, sort = FALSE)
joined <- joined[order(joined$input_row), , drop = FALSE]
joined$one_to_one_human_ortholog_status <- ifelse(
  joined$mouse_mapping_status == "ambiguous_mouse_gene_mapping", "ambiguous_mouse_gene_mapping",
  ifelse(joined$mouse_mapping_status == "no_mouse_gene_mapping", "no_mouse_gene_mapping",
    ifelse(joined$mouse_ensembl_id %in% accepted_ids, "accepted_Ensembl113_one_to_one_confidence1", "not_one_to_one_or_missing"))
)
joined$human_one_to_one_ensembl_id <- ifelse(joined$one_to_one_human_ortholog_status == "accepted_Ensembl113_one_to_one_confidence1", joined$hsapiens_homolog_ensembl_gene, "")
joined$human_one_to_one_gene_symbol <- ifelse(joined$one_to_one_human_ortholog_status == "accepted_Ensembl113_one_to_one_confidence1", joined$hsapiens_homolog_associated_gene_name, "")
joined$G5_module_gate <- modules$G5_module_gate[joined$input_row]
joined$G5_module_direction_hypothesis <- modules$G5_module_direction_hypothesis[joined$input_row]

out_file <- file.path(reg_dir, "G5_receiver_one_to_one_human_orthologs_Ensembl113.tsv.gz")
con <- gzfile(out_file, "wt", compression = 6)
write.table(joined, con, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
close(con)
attempt$status <- "completed"
attempt$n_mapped_one_to_one_rows <- sum(joined$one_to_one_human_ortholog_status == "accepted_Ensembl113_one_to_one_confidence1")
attempt$n_unique_human_genes <- length(unique(joined$human_one_to_one_ensembl_id[nzchar(joined$human_one_to_one_ensembl_id)]))
jsonlite::write_json(attempt, file.path(reg_dir, "G6_orthology_query_status.json"), auto_unbox = TRUE, pretty = TRUE)
cat("orthology_release=Ensembl113\n")
cat("frozen_module_gene_rows=", nrow(modules), "\n", sep = "")
cat("accepted_one_to_one_rows=", attempt$n_mapped_one_to_one_rows, "\n", sep = "")
cat("unique_human_genes=", attempt$n_unique_human_genes, "\n", sep = "")
