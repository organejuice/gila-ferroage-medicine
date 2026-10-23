#!/usr/bin/env Rscript

# Install the open-source RCTD implementation and missing dependencies into a
# project-local R library; do not alter the user's global R installation.

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) args[[1L]] else "."
lib <- if (length(args) >= 2L) args[[2L]] else file.path(root, "external/R-library/4.3.3/G6")
package_source <- if (length(args) >= 3L) args[[3L]] else ""
source_commit <- if (length(args) >= 4L) args[[4L]] else "9f5dc33c8060f946c6072a138b70e189636e1435"
source_archive_sha256 <- if (length(args) >= 5L) args[[5L]] else "ff525a671e1bfd5e2abace76dc6915ab27a3e3ae15e3d73c4e75ec49a8ee70e0"
dir.create(lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(lib, .libPaths()))
options(repos = c(CRAN = "https://cloud.r-project.org"), timeout = 600)

required_cran <- c("pals", "quadprog", "fields", "locfdr", "remotes")
missing_cran <- required_cran[!vapply(required_cran, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_cran)) {
  install.packages(missing_cran, lib = lib, dependencies = NA)
}
if (!requireNamespace("spacexr", quietly = TRUE)) {
  if (nzchar(package_source)) {
    remotes::install_local(package_source, lib = lib, dependencies = FALSE, upgrade = "never")
  } else {
    remotes::install_github("dmcable/spacexr", lib = lib, build_vignettes = FALSE,
                            dependencies = FALSE, upgrade = "never")
  }
}
if (!requireNamespace("spacexr", quietly = TRUE)) stop("spacexr did not install into the project-local library.")

desc <- packageDescription("spacexr", lib.loc = lib)
out_dir <- file.path(root, "results/2026-09-24-disease-side-v3.5-g6-input-audit")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
manifest <- c(
  paste0("R_version\t", R.version.string),
  paste0("local_library\t", normalizePath(lib, winslash = "/", mustWork = TRUE)),
  paste0("package\tspacexr"),
  paste0("version\t", desc$Version),
  paste0("remote_type\t", ifelse(is.null(desc$RemoteType), ifelse(nzchar(package_source), "github_archive", "github"), desc$RemoteType)),
  paste0("remote_repo\t", ifelse(is.null(desc$RemoteRepo), "dmcable/spacexr", desc$RemoteRepo)),
  paste0("remote_sha\t", ifelse(is.null(desc$RemoteSha), source_commit, desc$RemoteSha)),
  paste0("source_archive_sha256\t", source_archive_sha256),
  paste0("source_url\thttps://codeload.github.com/dmcable/spacexr/zip/", source_commit),
  paste0("license\t", desc$License)
)
writeLines(manifest, file.path(out_dir, "G6_RCTD_software_manifest.tsv"), useBytes = TRUE)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "G6_RCTD_R_sessionInfo.txt"), useBytes = TRUE)
cat(paste(manifest, collapse = "\n"), "\n", sep = "")
