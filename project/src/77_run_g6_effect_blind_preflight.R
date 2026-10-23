#!/usr/bin/env Rscript

# Execute the locked G6 primary script only through its effect-blind gate.
# The gate expression is located by its fail-closed message so source edits
# cannot silently move diagnosis-aware work into this preflight entry point.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: 77_run_g6_effect_blind_preflight.R <project_root>")
root <- args[[1L]]
analysis_path <- file.path(root, "src/73_run_g6_primary_analysis.R")
expressions <- parse(analysis_path, keep.source = FALSE)
sentinel <- "Effect-blind G6 gate failed. No G6 group effects were calculated."
gate_idx <- which(vapply(expressions, function(expr) {
  grepl(sentinel, paste(deparse(expr), collapse = " "), fixed = TRUE)
}, logical(1)))
if (length(gate_idx) != 1L) stop("Expected exactly one locked effect-blind stop sentinel in the G6 analysis script")

execution_env <- new.env(parent = globalenv())
for (i in seq_len(gate_idx)) eval(expressions[[i]], envir = execution_env)

preflight_path <- file.path(root, "results/2026-09-24-disease-side-v3.5-g6-primary/G6_effect_blind_preflight.json")
if (!file.exists(preflight_path)) stop("G6 effect-blind preflight did not write its required JSON")
preflight <- jsonlite::fromJSON(preflight_path, simplifyVector = FALSE)
if (!identical(preflight$effects_calculated, FALSE)) stop("Effect-blind preflight unexpectedly reports calculated effects")
if (!isTRUE(preflight$preflight_pass)) stop("Effect-blind G6 input gate failed; no group effects were calculated")
message("G6 effect-blind preflight passed; diagnosis-aware expressions were not evaluated.")
