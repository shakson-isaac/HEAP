#!/usr/bin/env Rscript

# ============================================================================
# install_r_packages.R — build/refresh the shared HEAP R library
# ----------------------------------------------------------------------------
# Installs the HEAP package set (config/r_packages.tsv) into the SHARED group
# R library so any hpc_patel member can run the HEAP R scripts without
# installing anything themselves. Run by the library maintainer.
#
# Target library: $HEAP_RLIB, else /n/groups/patel/IGLOO/Rlib/<R-version>
# (the same path 00_paths.R prepends to .libPaths()).
#
# Usage (maintainer):
#   module load gcc/14.2.0 R/4.4.2
#   HEAP_PATHS_FILE=/n/groups/patel/shakson_ukb/HEAP/workflow/00_paths.R \
#     Rscript scripts/setup/install_r_packages.R [--only-missing] [package ...]
#
#   --only-missing  install just packages not already resolvable (default if no
#                   explicit packages are given)
#   package ...     install only these (overrides the manifest)
#
# Notes:
#   - Bioconductor packages are installed via BiocManager.
#   - GitHub packages (TwoSampleMR, genetics.binaRies) via remotes.
#   - Compiled packages must match the loaded R module (gcc/14.2.0 R/4.4.2);
#     everyone running HEAP should use that module so the shared lib is ABI-safe.
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

# Target shared library (matches 00_paths.R logic)
rlib <- Sys.getenv("HEAP_RLIB",
  unset = file.path("/n/groups/patel/IGLOO", "Rlib", as.character(getRversion())))
dir.create(rlib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(rlib, .libPaths()))
message("Shared HEAP R library: ", rlib)

# --- read manifest ----------------------------------------------------------
manifest_path <- if (exists("heap_config")) heap_config("r_packages.tsv") else
  "/n/groups/patel/shakson_ukb/HEAP/config/r_packages.tsv"
man <- read.delim(manifest_path, stringsAsFactors = FALSE)
if (length(explicit)) man <- man[man$package %in% explicit, , drop = FALSE]

resolvable <- function(p) requireNamespace(p, quietly = TRUE)
if (only_missing) man <- man[!vapply(man$package, resolvable, logical(1)), , drop = FALSE]

if (!nrow(man)) { message("Nothing to install (all resolvable)."); quit(status = 0) }
message("To install (", nrow(man), "): ", paste(man$package, collapse = ", "))

cran   <- man$package[man$repo == "CRAN"]
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
# keep the shared lib group-readable for the team
try(system(paste("chmod -R g+rX", shQuote(rlib))), silent = TRUE)
