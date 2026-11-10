#!/usr/bin/env Rscript

# ============================================================================
# install_r_packages.R -- install the R packages HEAP needs
# ----------------------------------------------------------------------------
# Installs the packages in config/r_packages.tsv: CRAN packages with
# install.packages, Bioconductor packages with BiocManager, GitHub packages with
# remotes. Current releases are installed, except for the versions pinned below
# (PINNED), which are installed at exactly that version.
#
# Target library: $HEAP_RLIB if set; else the shared library on the original
# cluster if it exists; else your default R library (as install.packages uses).
#
# Usage (from the repository root, with R 4.4.x):
#   Rscript scripts/setup/install_r_packages.R [--only-missing] [package ...]
#
#   --only-missing  install only packages that are missing, or whose installed
#                   version differs from a pin (default when no packages are named)
#   package ...     install only these
#
# Packages built from source on Linux need C, C++ and Fortran compilers.
# ============================================================================

local({
  cand <- c(Sys.getenv("HEAP_PATHS_FILE", unset = ""),
            file.path(getwd(), "workflow", "00_paths.R"),
            "/n/groups/patel/shakson_ukb/HEAP/workflow/00_paths.R")
  hit <- cand[nzchar(cand) & file.exists(cand)][1]
  if (!is.na(hit)) source(hit)
})

args        <- commandArgs(trailingOnly = TRUE)
only_missing <- ("--only-missing" %in% args) ||
  !length(setdiff(args, c("--only-missing")))
explicit    <- setdiff(args, "--only-missing")

# Versions installed exactly rather than at the current release:
#   ggplot2 4.x breaks the ggtree -> enrichplot -> clusterProfiler stack;
#   qs has been archived on CRAN, so it is installed from the CRAN archive.
PINNED <- c(ggplot2 = "3.5.2", qs = "0.27.3")

# Target library: HEAP_RLIB, else the original cluster's shared library (same
# path 00_paths.R prepends to .libPaths()), else the user's default library.
cluster_lib <- file.path("/n/groups/patel/IGLOO", "Rlib", as.character(getRversion()))
rlib <- Sys.getenv("HEAP_RLIB", unset = "")
if (!nzchar(rlib)) {
  rlib <- if (dir.exists(cluster_lib) && file.access(cluster_lib, 2) == 0) cluster_lib else {
    user_lib <- Sys.getenv("R_LIBS_USER")
    writable <- .libPaths()[file.access(.libPaths(), 2) == 0]
    if (length(writable)) writable[1] else user_lib
  }
}
dir.create(rlib, recursive = TRUE, showWarnings = FALSE)
if (file.access(rlib, 2) != 0) stop("Cannot write to the R library ", rlib,
                                    "; set HEAP_RLIB to a writable directory.")
.libPaths(c(rlib, .libPaths()))
message("Installing into R library: ", rlib)

# --- read manifest ----------------------------------------------------------
manifest_path <- if (exists("heap_config")) heap_config("r_packages.tsv") else
  file.path(getwd(), "config", "r_packages.tsv")
man <- read.delim(manifest_path, stringsAsFactors = FALSE)
if (length(explicit)) man <- man[man$package %in% explicit, , drop = FALSE]

resolvable <- function(p) {
  if (!requireNamespace(p, quietly = TRUE)) return(FALSE)
  if (p %in% names(PINNED)) return(as.character(packageVersion(p)) == PINNED[[p]])
  TRUE
}
if (only_missing) man <- man[!vapply(man$package, resolvable, logical(1)), , drop = FALSE]

if (!nrow(man)) { message("Nothing to install (all resolvable)."); quit(status = 0) }
message("To install (", nrow(man), "): ", paste(man$package, collapse = ", "))

pinned <- intersect(man$package, names(PINNED))
cran   <- setdiff(man$package[man$repo == "CRAN"], pinned)
bioc   <- man$package[man$repo == "Bioconductor"]
github <- man[grepl("^GitHub:", man$repo), ]

ok <- character(0); failed <- character(0)
try_install <- function(p, fn) {
  res <- tryCatch({ fn(); if (requireNamespace(p, quietly = TRUE)) "ok" else "no" },
                  error = function(e) conditionMessage(e))
  if (identical(res, "ok")) { ok <<- c(ok, p); message("  [ok]   ", p) }
  else { failed <<- c(failed, p); message("  [FAIL] ", p, " : ", res) }
}

if (length(cran)) {
  repos <- "https://cloud.r-project.org"
  for (p in cran) try_install(p, function() install.packages(p, lib = rlib, repos = repos))
}
if (length(pinned)) {
  if (!requireNamespace("remotes", quietly = TRUE))
    install.packages("remotes", lib = rlib, repos = "https://cloud.r-project.org")
  for (p in pinned) try_install(p, function()
    remotes::install_version(p, version = PINNED[[p]], lib = rlib, upgrade = "never",
                             repos = "https://cloud.r-project.org"))
}
if (length(bioc)) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager", lib = rlib, repos = "https://cloud.r-project.org")
  for (p in bioc) try_install(p, function()
    BiocManager::install(p, lib = rlib, update = FALSE, ask = FALSE))
}
if (nrow(github)) {
  if (!requireNamespace("remotes", quietly = TRUE))
    install.packages("remotes", lib = rlib, repos = "https://cloud.r-project.org")
  for (i in seq_len(nrow(github))) {
    repo <- sub("^GitHub:", "", github$repo[i]); p <- github$package[i]
    try_install(p, function() remotes::install_github(repo, lib = rlib, upgrade = "never"))
  }
}

message("\nInstalled: ", length(ok), " | Failed: ", length(failed))
if (length(failed)) message("Failed: ", paste(failed, collapse = ", "))
# On the original cluster, keep the shared library group-readable.
if (identical(rlib, cluster_lib)) try(system(paste("chmod -R g+rX", shQuote(rlib))), silent = TRUE)
