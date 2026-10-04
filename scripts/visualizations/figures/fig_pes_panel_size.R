#!/usr/bin/env Rscript

# ============================================================================
# fig_pes_panel_size.R  [figure_id: fig_pes_panel_size]
# ----------------------------------------------------------------------------
# Module 6 (PES) supplement: predictive accuracy vs PES sparsity. Each point is
# one exposure: x = number of proteins the LASSO retained in the proteome-only
# PES, y = cross-validated R2. Colored by category. Shows accuracy is not simply
# a function of panel size — some exposures reach high R2 with compact panels.
#
# Reads ONLY precomputed tables (SelectedProteins + TrainOverallMetrics); the
# only transform is a join on exposure_id.
#
# Run directly:
#   module load gcc/14.2.0 R/4.4.2
#   HEAP_PATHS_FILE=/n/groups/patel/shakson_ukb/HEAP/workflow/00_paths.R \
#     Rscript scripts/visualizations/figures/fig_pes_panel_size.R base
# Or:
#   Rscript scripts/visualizations/build_figures.R --figure fig_pes_panel_size
# ============================================================================

local({
  cand <- c(file.path(getwd(), "scripts", "visualizations", "common"),
            "/n/groups/patel/shakson_ukb/HEAP/scripts/visualizations/common")
  common <- cand[dir.exists(cand)][1]
  if (is.na(common)) stop("fig_pes_panel_size.R: cannot locate common/ helpers")
  for (f in c("figure_paths.R", "load_heap_results.R", "plot_theme.R",
              "label_helpers.R", "export_helpers.R"))
    source(file.path(common, f))
})
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })

a <- commandArgs(trailingOnly = TRUE)
a <- a[!startsWith(a, "--")]
a <- a[!a %in% c("fig_pes_panel_size", "all_main", "all_supplement", "all", "website")]
covarType <- if (length(a) >= 1) a[1] else "base"
figure_id <- Sys.getenv("HEAP_FIGURE_ID", unset = "fig_pes_panel_size")

# --- load panel sizes + CV accuracy; proteome-only PES ----------------------
L  <- load_module6_pes_longitudinal(covarType, c("selected", "overall"))
sel <- L$selected[model == "prot_only", .(exposure_id, n_proteins = n_selected_proteins)]
ov  <- L$overall[model == "prot_only" & is.finite(r2),
                 .(exposure_id, category, r2, correlation)]
d <- merge(ov, sel, by = "exposure_id")
d <- d[is.finite(n_proteins) & n_proteins > 0]
if (!nrow(d)) stop("No joined panel-size/accuracy rows for covarType=", covarType)

d[, exposure_label := heap_exposure_label(exposure_id)]
d[, category := heap_category_factor(category)]
# label the most predictive exposures only (keep the panel readable)
d[, lab := ifelse(rank(-r2) <= 8L, exposure_label, "")]

# --- plot: scatter, log-x panel size ----------------------------------------
p <- ggplot(d, aes(n_proteins, r2, colour = category)) +
  geom_point(size = 2.4, alpha = 0.9) +
  ggrepel::geom_text_repel(aes(label = lab), size = 2.7, colour = "grey20",
                           max.overlaps = 20, min.segment.length = 0,
                           seed = 1, show.legend = FALSE) +
  scale_colour_exposure(drop = TRUE) +
  scale_x_log10() +
  labs(title = "Predictive accuracy vs PES panel size",
       subtitle = NULL,
       x = "Proteins retained in PES (log scale)",
       y = expression("Cross-validated out-of-fold "*R^2)) +
  theme_heap()

heap_emit_figure(p, figure_id,
                 data = d[, .(exposure_id, exposure_label, category, n_proteins, r2, correlation)],
                 category = "supplement", formats = c("pdf", "png"),
                 width = 7.6, height = 5.4, website = FALSE)

message("fig_pes_panel_size: done (", covarType, "; ", nrow(d), " exposures).")
