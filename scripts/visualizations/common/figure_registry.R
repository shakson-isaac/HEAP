#!/usr/bin/env Rscript

# ============================================================================
# figure_registry.R — read + query the HEAP figure registry
# ----------------------------------------------------------------------------
# The registry (config/figures/figure_registry.tsv) is the catalog of every
# figure the repository can build: id, type, generating script, required module
# outputs, loader, target output directory, website-export flag, and status.
#
# This accessor lets the build workflow and docs query the registry without
# re-parsing the TSV by hand.
# ============================================================================

local({
  if (exists("heap_config", mode = "function")) return(invisible())
  cand <- c(file.path(getwd(), "scripts", "visualizations", "common", "figure_paths.R"),
            "/n/groups/patel/shakson_ukb/HEAP/scripts/visualizations/common/figure_paths.R")
  hit <- cand[file.exists(cand)][1]
  if (!is.na(hit)) source(hit)
})

suppressPackageStartupMessages({ library(data.table) })

heap_figure_registry_path <- function() heap_config("figures", "figure_registry.tsv")

#' Load the figure registry as a data.table.
heap_figure_registry <- function() {
  f <- heap_figure_registry_path()
  if (!file.exists(f))
    stop("Figure registry not found: ", f,
         "\nExpected config/figures/figure_registry.tsv", call. = FALSE)
  fread(f, sep = "\t", fill = TRUE)
}

#' Look up a single registry row by figure_id (exact match).
heap_figure_lookup <- function(figure_id) {
  reg <- heap_figure_registry()
  row <- reg[reg$figure_id == figure_id, ]
  if (nrow(row) == 0L)
    stop("Unknown figure_id '", figure_id, "'. Known ids:\n  ",
         paste(reg$figure_id, collapse = "\n  "), call. = FALSE)
  row
}

#' Filter the registry by output_dir (main/supplement/exploratory/website),
#' status, figure_type, or website_export flag.
heap_figures_where <- function(output_dir = NULL, status = NULL,
                               figure_type = NULL, website_export = NULL) {
  reg <- heap_figure_registry()
  if (!is.null(output_dir))     reg <- reg[reg$output_dir %in% output_dir, ]
  if (!is.null(status))         reg <- reg[reg$status %in% status, ]
  if (!is.null(figure_type))    reg <- reg[reg$figure_type %in% figure_type, ]
  if (!is.null(website_export)) reg <- reg[reg$website_export %in% website_export, ]
  reg
}

#' Print a compact human-readable summary of the registry.
heap_figure_registry_summary <- function() {
  reg <- heap_figure_registry()
  cat("HEAP figure registry:", nrow(reg), "figures\n")
  cat("By status:\n"); print(table(reg$status))
  cat("\nBy output_dir:\n"); print(table(reg$output_dir))
  cat("\nWebsite-export figures:", sum(reg$website_export == "yes"), "\n")
  invisible(reg)
}

# ===========================================================================
# Module-subdir mapping (registry-driven figure output organization)
# ---------------------------------------------------------------------------
# Each figure is grouped on disk under its producing module, so the IGLOO figure
# tree is figures/<category>/<module>/<figure_id>.{pdf,png} instead of a flat
# dump. The grouping is DERIVED from the registry (required_modules), so adding a
# figure needs no path edits — only a registry row + a generating script.
# ===========================================================================

# canonical producing-module -> output sub-directory name. The enrichment module
# (module4_enrichment) only ever appears as the 2nd token of a module2 figure, so
# its manuscript panels (Fig2E/Fig3H/FigS8/FigS9) live under module2.
HEAP_MODULE_SUBDIRS <- c(
  module1                 = "module1",
  module2                 = "module2",
  module3                 = "module3",
  module4_enrichment      = "module2",
  module5                 = "module5",
  module6                 = "module6",
  gwas_regenie            = "gwas",
  ldsc                    = "gwas",
  population_architecture = "population_architecture"
)

#' Canonical output sub-directory for a figure_id, from the registry.
#'
#' Uses the FIRST token of the row's `required_modules` (";"-separated). The
#' first-token rule is correct for every multi-module registry row
#' (module2;module4_enrichment -> module2, module1;population_architecture ->
#' module1, module2;module5 -> module2). Returns:
#'   * the mapped module subdir (e.g. "module1");
#'   * "shared_qc" when the figure has no module input ("none"/"none(loader)")
#'     or an unmapped module token;
#'   * NULL when the figure_id is not in the registry, so callers fall back to
#'     the flat category directory (covers ad-hoc panel ids like
#'     fig_mr_shared_paths_01).
#' @param figure_id registry figure_id
heap_figure_module_subdir <- function(figure_id) {
  reg <- heap_figure_registry()
  # use which() on a plain vector — NOT reg[reg$figure_id == figure_id, ] — to
  # avoid data.table's i-scoping, where the `figure_id` argument would resolve to
  # the same-named column and match every row.
  idx <- which(reg$figure_id == figure_id)
  if (length(idx) == 0L) return(NULL)
  mods  <- trimws(strsplit(as.character(reg$required_modules[idx[1]]), ";", fixed = TRUE)[[1]])
  first <- mods[nzchar(mods)][1]
  if (is.na(first) || first %in% c("none", "none(loader)")) return("shared_qc")
  sub <- HEAP_MODULE_SUBDIRS[first]   # single-bracket: unknown key -> NA (not an error)
  if (is.na(sub)) "shared_qc" else unname(sub)
}

# ===========================================================================
# Manuscript-panel shaping helpers (shared by build_figures.R + build_legends.R)
# ---------------------------------------------------------------------------
# `manuscript_ref` ties a figure_id to one or more manuscript panels (e.g.
# "Fig2A", "Fig5C;Fig5D;Fig5E", or "extra"). These pure registry-shaping helpers
# expand and order panels; they live here so the figure-build and legend-build
# drivers share one copy.
# ===========================================================================

#' Sortable key for a manuscript panel ref: main before supplement, then figure
#' number, then panel letter (e.g. "Fig2A" -> "0-02-A", "FigS1B" -> "1-01-B").
heap_ms_sort_key <- function(ref) {
  supp <- grepl("^FigS", ref)
  num  <- suppressWarnings(as.integer(sub("^FigS?([0-9]+).*$", "\\1", ref)))
  pan  <- sub("^FigS?[0-9]+", "", ref)
  sprintf("%d-%02d-%s", as.integer(supp), ifelse(is.na(num), 99L, num), pan)
}

#' Expand a registry subset into one row per manuscript panel (splitting
#' ";"-listed refs), carrying figure_id/script/status/output_dir/figure_class.
heap_ms_expand <- function(d) {
  rbindlist(lapply(seq_len(nrow(d)), function(i) {
    panels <- trimws(strsplit(d$manuscript_ref[i], ";", fixed = TRUE)[[1]])
    data.table(panel = panels, figure_id = d$figure_id[i], script = d$script[i],
               status = d$status[i], output_dir = d$output_dir[i],
               figure_class = d$figure_class[i])
  }))
}
