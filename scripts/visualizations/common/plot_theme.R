#!/usr/bin/env Rscript

# ============================================================================
# plot_theme.R — shared ggplot2 theme + palettes for HEAP figures
# ----------------------------------------------------------------------------
# Centralizes the look-and-feel so every figure script produces consistent,
# publication-ready output. Source after loading ggplot2.
# ============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
})

#' Base HEAP theme — clean, manuscript-ready.
theme_heap <- function(base_size = 11, base_family = "") {
  theme_bw(base_size = base_size, base_family = base_family) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),   # clean white default; add geom_hline/vline only for meaningful 0/boundary lines
      panel.border     = element_rect(colour = "grey40", linewidth = 0.4),
      axis.title       = element_text(face = "bold"),
      plot.title       = element_text(face = "bold", hjust = 0, size = rel(1.05)),
      plot.subtitle    = element_text(colour = "grey30", size = rel(0.82)),
      plot.title.position   = "plot",
      plot.caption.position = "plot",
      strip.background = element_rect(fill = "grey95", colour = "grey40", linewidth = 0.3),
      strip.text       = element_text(face = "bold"),
      legend.key       = element_blank(),
      legend.position  = "right"
    )
}

#' Wrap a long plot subtitle/caption so it does not run off the device edge.
#' ggplot2's element_text does not wrap; long one-line subtitles get clipped at
#' the right edge when the figure is saved at a fixed width. Insert line breaks
#' at word boundaries (base strwrap; no extra dependency). NULL/expression/empty
#' are returned unchanged.
heap_sub <- function(x, width = 100) {
  if (is.null(x) || !is.character(x) || length(x) != 1 || !nzchar(x)) return(x)
  paste(strwrap(x, width = width), collapse = "\n")
}

# ---------------------------------------------------------------------------
# Canonical palettes
# ---------------------------------------------------------------------------

#' Variance-component palette (genetic vs exposomic vs interaction vs covars).
HEAP_PAL_COMPONENT <- c(
  Covars  = "#9E9E9E",
  Genetic = "#1B6CA8",  # blue
  G       = "#1B6CA8",
  PGS     = "#1B6CA8",
  Exposome= "#2E9E48",  # green (exposome = environment; was orange #E07B39)
  PXS     = "#2E9E48",
  E       = "#2E9E48",
  GxE     = "#7B3FA0",  # purple
  Interaction = "#7B3FA0",
  Residual= "#D9D9D9"
)

# ---------------------------------------------------------------------------
# CANONICAL exposure-category palette (colour-blind-safe distinct; palette "C").
# SINGLE SOURCE OF TRUTH for category colour anywhere in the paper. Keys are the
# exact fine `Category` values used in the Module 2 output; the documented copy
# lives in config/figures/exposure_category_palette.tsv (generated from here).
# Use the scale_*_exposure() / heap_category_factor() helpers below rather than
# referencing the vector directly, so colour + legend order stay consistent.
# ---------------------------------------------------------------------------

# fine categories in canonical (broad-group) order = legend order
HEAP_ECAT_LEVELS <- c(
  "Alcohol", "Smoking", "Vitamins",                                  # Substance Use
  "Diet_Weekly", "Exercise_Freq", "Exercise_MET",                    # Diet & Activity
  "Internet_Usage", "Sleep", "Sexual_Factors", "Sun_Exposure",       # Behaviour & Screen
  "Residential_Air_Pollution", "Residential_Noise_Pollution",        # Environmental: Physical
  "Deprivation_Indices")                                             # Environmental: Socioeconomic

HEAP_ECAT_COLORS <- c(
  Alcohol                     = "#D55E00",
  Smoking                     = "#444444",
  Vitamins                    = "#E69F00",
  Diet_Weekly                 = "#009E73",
  Exercise_Freq               = "#117733",
  Exercise_MET                = "#44AA99",
  Internet_Usage              = "#CC79A7",
  Sleep                       = "#0072B2",
  Sexual_Factors              = "#AA4499",
  Sun_Exposure                = "#DDCC77",
  Residential_Air_Pollution   = "#882255",
  Residential_Noise_Pollution = "#999933",
  Deprivation_Indices         = "#56B4E9",
  Other                       = "#BDBDBD"
)

# broad-group anchor colours (keyed by BOTH the broad_group code and the
# broad_group_label from config/exposure_sets/analysis_exposure_category_groups.tsv,
# so heap_broad_category(label=TRUE/FALSE) output maps either way).
HEAP_BROAD_COLORS <- c(
  Lifestyle_Substance_Use                 = "#D55E00",
  Lifestyle_Diet_Activity                 = "#117733",
  Lifestyle_Behaviour                     = "#0072B2",
  Environment_Physical                    = "#882255",
  Environment_Socioeconomic               = "#56B4E9",
  "Lifestyle: Substance Use"              = "#D55E00",
  "Lifestyle: Diet and Activity"          = "#117733",
  "Lifestyle: Behavior and Screen Time"  = "#0072B2",
  "Environmental: Physical Exposures"     = "#882255",
  "Environmental: Socioeconomic"          = "#56B4E9"
)

# Backwards-compatible alias (older scripts referenced HEAP_PAL_ECAT).
HEAP_PAL_ECAT <- HEAP_ECAT_COLORS

#' Canonical named colour vector for fine exposure categories.
heap_category_colors <- function() HEAP_ECAT_COLORS

#' Canonical named colour vector for broad exposure groups (code or label keys).
heap_broad_colors <- function() HEAP_BROAD_COLORS

#' Coerce an exposure-category vector to a factor in the canonical legend order.
#' Categories outside the known set are appended after the canonical ones.
heap_category_factor <- function(x) {
  x <- as.character(x)
  extra <- setdiff(unique(x[!is.na(x)]), HEAP_ECAT_LEVELS)
  factor(x, levels = c(HEAP_ECAT_LEVELS, sort(extra)))
}

#' DISPLAY name for a fine exposure category: the canonical code is kept for
#' colour/factor mapping, but underscores become spaces for readability on plots
#' (e.g. "Exercise_Freq" -> "Exercise Freq"). Used as the default `labels` of the
#' exposure scales (legends) and applied to category AXES via scale_*_discrete.
heap_category_pretty <- function(x) gsub("_", " ", as.character(x))

#' ggplot scale: colour by fine exposure category using the canonical palette.
#' Legend labels are prettified (spaces, not underscores) by default.
#' @param drop FALSE keeps every category in the legend (consistent across panels)
#' @param na.value colour for unmapped/NA categories (default grey)
scale_colour_exposure <- function(name = "Exposure category", drop = TRUE,
                                  na.value = "grey80", labels = heap_category_pretty, ...) {
  ggplot2::scale_colour_manual(values = HEAP_ECAT_COLORS, breaks = HEAP_ECAT_LEVELS,
                               name = name, drop = drop, na.value = na.value,
                               labels = labels, ...)
}
scale_color_exposure <- scale_colour_exposure
scale_fill_exposure <- function(name = "Exposure category", drop = TRUE,
                                na.value = "grey80", labels = heap_category_pretty, ...) {
  ggplot2::scale_fill_manual(values = HEAP_ECAT_COLORS, breaks = HEAP_ECAT_LEVELS,
                             name = name, drop = drop, na.value = na.value,
                             labels = labels, ...)
}

#' ggplot scale: colour/fill by broad exposure group using the canonical anchors.
scale_colour_broad <- function(name = "Exposure group", drop = TRUE,
                               na.value = "grey80", ...)
  ggplot2::scale_colour_manual(values = HEAP_BROAD_COLORS, name = name,
                               drop = drop, na.value = na.value, ...)
scale_color_broad <- scale_colour_broad
scale_fill_broad <- function(name = "Exposure group", drop = TRUE,
                             na.value = "grey80", ...)
  ggplot2::scale_fill_manual(values = HEAP_BROAD_COLORS, name = name,
                             drop = drop, na.value = na.value, ...)

# ---------------------------------------------------------------------------
# Canonical GxE genetic-component palette (cis vs trans vs joint-only). cis =
# the protein's own (local) genetic score interacting with the exposure; trans
# = a distal genetic score; joint-only = significant on the joint F-test but
# neither cis nor trans alone. Used wherever GxE is split by component.
# ---------------------------------------------------------------------------
HEAP_GXE_COLORS <- c(cis = "#1B6CA8", trans = "#D55E00", `joint-only` = "#9E9E9E")
HEAP_GXE_LEVELS <- c("cis", "trans", "joint-only")
scale_colour_gxe <- function(name = "GxE component", drop = TRUE, ...)
  ggplot2::scale_colour_manual(values = HEAP_GXE_COLORS, breaks = HEAP_GXE_LEVELS,
                               name = name, drop = drop, ...)
scale_color_gxe <- scale_colour_gxe
scale_fill_gxe <- function(name = "GxE component", drop = TRUE, ...)
  ggplot2::scale_fill_manual(values = HEAP_GXE_COLORS, breaks = HEAP_GXE_LEVELS,
                             name = name, drop = drop, ...)

#' Significance-star helper (shared across many legacy scripts).
heap_pstars <- function(p) {
  ifelse(is.na(p), "",
    ifelse(p < 0.001, "***",
      ifelse(p < 0.01, "**",
        ifelse(p < 0.05, "*", ""))))
}

#' Null-coalescing operator used widely in the legacy plotting code.
`%||%` <- function(a, b) if (!is.null(a)) a else b

# --- rasterize dense point layers -------------------------------------------
# Manhattan/Miami plots embed millions of points; as pure vector the PDF runs to
# 15-25 MB and a viewer stalls rendering the page. heap_rasterize() rasterizes ONLY
# the point layers and leaves text, axes and reference lines as vector, so labels
# stay sharp and selectable while the file drops ~100x. Apply to any plot whose
# point count is in the 10^5-10^7 range (exwas/gxe Miami, GWAS Manhattan, QQ).
# ggrastr is not in the module R -- it lives in a project lib.
local({
  lib <- "/n/groups/patel/shakson_ukb/Rlib/rastr"
  if (dir.exists(lib) && !lib %in% .libPaths()) .libPaths(c(lib, .libPaths()))
})

heap_rasterize <- function(p, dpi = 300) {
  if (!requireNamespace("ggrastr", quietly = TRUE)) {
    warning("ggrastr unavailable; emitting full-vector plot (large file)", call. = FALSE)
    return(p)
  }
  ggrastr::rasterise(p, layers = "Point", dpi = dpi)
}
