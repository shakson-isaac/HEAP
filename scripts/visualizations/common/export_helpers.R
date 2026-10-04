#!/usr/bin/env Rscript

# ============================================================================
# export_helpers.R — figure + figure-data export for HEAP visualization
# ----------------------------------------------------------------------------
# Standardizes how figures and their underlying plotted-data tables are written
# to the canonical IGLOO figure tree, so every panel ships with reproducible
# data and (optionally) a website-ready asset.
#
# Principle: never save a plot without also saving the data behind it.
# ============================================================================

local({
  if (exists("heap_figure_path", mode = "function")) return(invisible())
  cand <- c(file.path(getwd(), "scripts", "visualizations", "common", "figure_paths.R"),
            "/n/groups/patel/shakson_ukb/HEAP/scripts/visualizations/common/figure_paths.R")
  hit <- cand[file.exists(cand)][1]
  if (is.na(hit)) stop("export_helpers.R: cannot find figure_paths.R")
  source(hit)
})

# Registry-driven module subdir lookup (heap_figure_module_subdir). Figure scripts
# source export_helpers but not figure_registry, so source it here (guarded) — else
# the per-module nesting silently won't happen.
local({
  if (exists("heap_figure_module_subdir", mode = "function")) return(invisible())
  cand <- c(file.path(getwd(), "scripts", "visualizations", "common", "figure_registry.R"),
            "/n/groups/patel/shakson_ukb/HEAP/scripts/visualizations/common/figure_registry.R")
  hit <- cand[file.exists(cand)][1]
  if (!is.na(hit)) source(hit)
})

# Resolve the output sub-directory for a figure_id, tolerating a missing/failed
# registry (falls back to NULL = flat category dir, so emission never breaks).
.heap_emit_subdir <- function(figure_id, subdir) {
  if (!is.null(subdir)) return(subdir)
  if (!exists("heap_figure_module_subdir", mode = "function")) return(NULL)
  tryCatch(heap_figure_module_subdir(figure_id), error = function(e) NULL)
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(data.table)
})

#' Save a ggplot (or base/grob via a draw function) to a figure category.
#'
#' Writes <figure_id>.<ext> for each requested format into the category dir.
#'
#' @param plot a ggplot object (or NULL if `draw` is supplied)
#' @param figure_id stable id / file stem, e.g. "figure2_r2_decomposition"
#' @param category one of main/supplement/exploratory/website
#' @param formats one or more of "pdf","png","svg"
#' @param width,height inches
#' @param dpi raster dpi
#' @param draw optional function(){...} that draws to the current device
#'   (for non-ggplot figures); used when `plot` is NULL
#' @return invisibly, the vector of written file paths
heap_save_figure <- function(plot, figure_id, category = "exploratory",
                             formats = c("pdf", "png"),
                             width = 7, height = 5, dpi = 300, draw = NULL,
                             subdir = NULL) {
  formats <- match.arg(formats, c("pdf", "png", "svg"), several.ok = TRUE)
  subdir  <- .heap_emit_subdir(figure_id, subdir)
  written <- character(0)
  for (fmt in formats) {
    path <- heap_figure_path(paste0(figure_id, ".", fmt), category = category,
                             subdir = subdir)
    if (!is.null(plot) && inherits(plot, "ggplot")) {
      ggsave(path, plot = plot, width = width, height = height, dpi = dpi,
             device = fmt, limitsize = FALSE)
    } else if (!is.null(draw)) {
      dev <- switch(fmt,
        pdf = grDevices::pdf(path, width = width, height = height),
        png = grDevices::png(path, width = width, height = height,
                             units = "in", res = dpi),
        svg = grDevices::svg(path, width = width, height = height))
      draw(); grDevices::dev.off()
    } else {
      stop("heap_save_figure: provide a ggplot `plot` or a `draw` function.")
    }
    written <- c(written, path)
  }
  loc <- if (is.null(subdir)) category else file.path(category, subdir)
  message("  [figure] ", figure_id, " -> ", loc, "/ (",
          paste(formats, collapse = ","), ")")
  invisible(written)
}

#' Save the data table behind a figure panel (always, for reproducibility).
#'
#' @param data data.frame/data.table of plotted values
#' @param figure_id stable id; file written to figures/data/<figure_id>.tsv
#' @return invisibly the written path
heap_save_figure_data <- function(data, figure_id, subdir = NULL) {
  subdir <- .heap_emit_subdir(figure_id, subdir)
  path <- heap_figure_data_path(paste0(figure_id, ".tsv"), subdir = subdir)
  fwrite(as.data.table(data), path, sep = "\t")
  loc <- if (is.null(subdir)) "data" else file.path("data", subdir)
  message("  [data]   ", figure_id, " -> ", loc, "/", basename(path))
  invisible(path)
}

#' Export a website-ready asset (JSON for charts, or a copied PNG/SVG).
#'
#' @param data optional data.frame -> written as JSON for the site to consume
#' @param figure_id stable id
#' @param copy_from optional existing image path to copy into figures/website/
#' @return invisibly the written path(s)
heap_export_website <- function(figure_id, data = NULL, copy_from = NULL) {
  out <- character(0)
  if (!is.null(data)) {
    if (!requireNamespace("jsonlite", quietly = TRUE))
      stop("heap_export_website: install 'jsonlite' to export website JSON.")
    p <- heap_figure_website_path(paste0(figure_id, ".json"))
    writeLines(jsonlite::toJSON(data, dataframe = "rows", na = "null",
                                auto_unbox = TRUE, pretty = TRUE), p)
    out <- c(out, p)
  }
  if (!is.null(copy_from) && file.exists(copy_from)) {
    p <- heap_figure_website_path(basename(copy_from))
    file.copy(copy_from, p, overwrite = TRUE)
    out <- c(out, p)
  }
  if (length(out)) message("  [website] ", figure_id, " -> website/")
  invisible(out)
}

#' One-call convenience: save figure + its data (+ optional website export).
heap_emit_figure <- function(plot, figure_id, data, category = "main",
                             formats = c("pdf", "png"),
                             width = 7, height = 5, website = FALSE,
                             subdir = NULL, ...) {
  subdir <- .heap_emit_subdir(figure_id, subdir)
  fig <- heap_save_figure(plot, figure_id, category = category, subdir = subdir,
                          formats = formats, width = width, height = height, ...)
  dat <- heap_save_figure_data(data, figure_id, subdir = subdir)
  web <- if (website) heap_export_website(figure_id, data = data) else NULL
  invisible(list(figure = fig, data = dat, website = web))
}
