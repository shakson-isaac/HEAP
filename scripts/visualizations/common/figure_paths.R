#!/usr/bin/env Rscript

# ============================================================================
# figure_paths.R — canonical IGLOO figure + result paths for HEAP visualization
# ----------------------------------------------------------------------------
# This is the single source of truth for *where* visualization scripts read
# module outputs and *where* they write figures. It builds on the module path
# system in workflow/00_paths.R; it never hard-codes scratch or UK_Biobank
# locations.
#
# Module outputs are read from the canonical IGLOO project output
# (heap_project_output(), i.e. /n/groups/patel/IGLOO/UKB/HEAP/output) with a
# transparent fall-back to the local HEAP/output staging directory while the
# IGLOO migration of individual modules is still in progress.
#
# Figures are written under the canonical IGLOO figure tree:
#   /n/groups/patel/IGLOO/UKB/HEAP/figures/
#     main/        publication main-text figure panels
#     supplement/  supplementary figure panels
#     exploratory/ scratch / diagnostic figures
#     website/     website-ready exported assets
#     data/        underlying plotted-data tables (one per panel)
#     logs/        figure build logs
#
# Source this file from any visualization script:
#   source(".../scripts/visualizations/common/figure_paths.R")
# ============================================================================

# --- Locate and source the module path system (workflow/00_paths.R) ---------
local({
  if (exists("heap_project_output", mode = "function")) return(invisible())
  candidates <- c(
    Sys.getenv("HEAP_PATHS_FILE", unset = ""),
    file.path(getwd(), "workflow", "00_paths.R"),
    file.path(getwd(), "..", "workflow", "00_paths.R"),
    file.path(getwd(), "..", "..", "workflow", "00_paths.R"),
    file.path(getwd(), "..", "..", "..", "workflow", "00_paths.R"),
    file.path(getwd(), "..", "..", "..", "..", "workflow", "00_paths.R"),
    "/n/groups/patel/shakson_ukb/HEAP/workflow/00_paths.R"
  )
  candidates <- candidates[nzchar(candidates)]
  hit <- candidates[file.exists(candidates)][1]
  if (is.na(hit))
    stop("figure_paths.R: could not locate workflow/00_paths.R. ",
         "Set HEAP_PATHS_FILE=/path/to/HEAP/workflow/00_paths.R")
  source(hit)
})

# ---------------------------------------------------------------------------
# Figure output tree (IGLOO-rooted, canonical)
# ---------------------------------------------------------------------------
heap_figures_root <- function(...) heap_project_root("figures", ...)

# Allowed figure categories (also the immediate sub-directory names).
HEAP_FIGURE_CATEGORIES <- c("main", "supplement", "exploratory",
                            "website", "data", "logs")

#' Directory for a figure category, created on demand.
#' @param category one of HEAP_FIGURE_CATEGORIES
heap_figure_dir <- function(category = "exploratory", ...) {
  category <- match.arg(category, HEAP_FIGURE_CATEGORIES)
  d <- heap_figures_root(category, ...)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}

#' Full path to a figure file inside a category directory.
#' @param category one of HEAP_FIGURE_CATEGORIES
#' @param filename file name, e.g. "figure2_r2_decomposition.pdf"
#' @param subdir   optional nested sub-directory (e.g. a figure_id)
heap_figure_path <- function(filename, category = "exploratory", subdir = NULL) {
  base <- if (is.null(subdir)) heap_figure_dir(category)
          else heap_figure_dir(category, subdir)
  file.path(base, filename)
}

#' Path for a website-ready figure data table (TSV/JSON for the HEAP website).
heap_figure_data_path <- function(filename, subdir = NULL) {
  heap_figure_path(filename, category = "data", subdir = subdir)
}

#' Path for a website-ready exported asset (PNG/SVG/JSON consumed by the site).
heap_figure_website_path <- function(filename, subdir = NULL) {
  heap_figure_path(filename, category = "website", subdir = subdir)
}

#' Path for a figure build log.
heap_figure_log_path <- function(filename) {
  heap_figure_path(filename, category = "logs")
}

# ---------------------------------------------------------------------------
# Module-output resolver
#
# Modules are mid-migration: some write to IGLOO (heap_project_output, e.g.
# Module 2, Module 5), some still write to local HEAP/output (heap_output, e.g.
# Module 1, Module 6) during pilots. Visualization code should treat IGLOO as
# canonical but transparently fall back to local staging so figures can be
# built from whichever location currently holds the run.
#
# Override search order entirely with HEAP_RESULTS_ROOT to pin a single root.
# ---------------------------------------------------------------------------

#' Resolve a module-output sub-directory to an existing canonical location.
#'
#' Search order:
#'   1. $HEAP_RESULTS_ROOT/<subdir>          (explicit pin, if set)
#'   2. heap_project_output(<subdir>)        (IGLOO canonical)
#'   3. heap_output(<subdir>)                (local HEAP/output staging)
#'
#' @param subdir module output sub-directory, e.g. "module2/Type3"
#' @param must_exist if TRUE (default) stop with an informative message when no
#'   candidate exists; if FALSE return the canonical IGLOO path regardless.
#' @return absolute path to the resolved directory
heap_resolve_output <- function(subdir, must_exist = TRUE) {
  pin <- Sys.getenv("HEAP_RESULTS_ROOT", unset = "")
  candidates <- c(
    if (nzchar(pin)) file.path(pin, subdir) else character(0),
    heap_project_output(subdir),
    heap_output(subdir)
  )
  hit <- candidates[dir.exists(candidates)][1]
  if (!is.na(hit)) return(normalizePath(hit, mustWork = FALSE))
  if (!must_exist) return(heap_project_output(subdir))
  stop(
    "Missing module output: '", subdir, "'.\n",
    "Looked in:\n  - ", paste(candidates, collapse = "\n  - "), "\n",
    "Run the upstream module/experiment that produces this output, or set ",
    "HEAP_RESULTS_ROOT to the directory that holds it.",
    call. = FALSE
  )
}

#' Convenience: TRUE if a module output sub-directory resolves to something real.
heap_output_exists <- function(subdir) {
  res <- tryCatch(heap_resolve_output(subdir, must_exist = TRUE),
                  error = function(e) NULL)
  !is.null(res)
}
