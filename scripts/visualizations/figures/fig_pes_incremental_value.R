#!/usr/bin/env Rscript

# ============================================================================
# fig_pes_incremental_value.R  [figure_id: fig_pes_incremental_value]
# ----------------------------------------------------------------------------
# Module 6 (PES): incremental predictive value of the proteome OVER covariates.
# Dumbbell of cross-validated R2 from the covariate-only baseline (cov_only) to
# the proteome+covariate model (prot_plus_cov) for each exposure, ranked by the
# gain. Shows that plasma proteins carry exposure information beyond demographics.
#
# Reads ONLY the precomputed CV metric table (TrainOverallMetrics); the only
# transform is a wide reshape across the `model` levels - no statistics here.
#
# Run directly:
#   module load gcc/14.2.0 R/4.4.2
#   HEAP_PATHS_FILE=/n/groups/patel/shakson_ukb/HEAP/workflow/00_paths.R \
#     Rscript scripts/visualizations/figures/fig_pes_incremental_value.R base
# Or:
#   Rscript scripts/visualizations/build_figures.R --figure fig_pes_incremental_value
# ============================================================================

local({
  cand <- c(file.path(getwd(), "scripts", "visualizations", "common"),
            "/n/groups/patel/shakson_ukb/HEAP/scripts/visualizations/common")
  common <- cand[dir.exists(cand)][1]
  if (is.na(common)) stop("fig_pes_incremental_value.R: cannot locate common/ helpers")
  for (f in c("figure_paths.R", "load_heap_results.R", "plot_theme.R",
              "label_helpers.R", "export_helpers.R"))
    source(file.path(common, f))
})
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })

a <- commandArgs(trailingOnly = TRUE)
a <- a[!startsWith(a, "--")]
a <- a[!a %in% c("fig_pes_incremental_value", "all_main", "all_supplement", "all", "website")]
covarType <- if (length(a) >= 1) a[1] else "base"
top_n     <- if (length(a) >= 2) as.integer(a[2]) else 30L
figure_id <- Sys.getenv("HEAP_FIGURE_ID", unset = "fig_pes_incremental_value")

PAL <- c(`Covariate baseline` = "#9E9E9E", `Proteome + covariates` = "#0072B2")

# --- load + reshape wide across models --------------------------------------
# Performance metric per exposure: R2 for continuous, AUC for binary one-hot
# exposures (R2 undefined). Facet by metric so ALL exposures are represented.
ov <- load_module6_pes_longitudinal(covarType, "overall")
ov[, perf := fifelse(exposure_type == "binary",
                     suppressWarnings(as.numeric(auc)),
                     suppressWarnings(as.numeric(r2)))]
w  <- dcast(ov, exposure_id + category + exposure_type ~ model, value.var = "perf")
need <- c("cov_only", "prot_plus_cov")
if (!all(need %in% names(w)))
  stop("TrainOverallMetrics missing model(s): ", paste(setdiff(need, names(w)), collapse = ", "))
w <- w[is.finite(cov_only) & is.finite(prot_plus_cov)]
w[, gain := prot_plus_cov - cov_only]
w[, metric := fifelse(exposure_type == "binary",
                      "Binary exposures - CV AUC",
                      "Continuous exposures - CV R²")]
n_cont <- w[exposure_type == "continuous", .N]
n_bin  <- w[exposure_type == "binary", .N]

w[, exposure_label := heap_exposure_label(exposure_id)]
# rank within each facet by gain (facets hold disjoint exposures, so a global gain
# order yields the correct per-facet ranking under scales = "free_y")
setorder(w, gain)
w[, ykey := factor(exposure_id, levels = unique(exposure_id))]
ylabs <- setNames(w$exposure_label, w$exposure_id)
w[, metric := factor(metric, levels = c("Continuous exposures - CV R²",
                                        "Binary exposures - CV AUC"))]

# long form for the two endpoint points (carry metric for faceting)
pts <- melt(w[, .(ykey, metric, `Covariate baseline` = cov_only,
                  `Proteome + covariates` = prot_plus_cov)],
            id.vars = c("ykey", "metric"), variable.name = "model", value.name = "perf")

# --- plot: dumbbell, faceted by metric --------------------------------------
p <- ggplot() +
  geom_segment(data = w,
               aes(x = cov_only, xend = prot_plus_cov, y = ykey, yend = ykey),
               colour = "grey70", linewidth = 0.5) +
  geom_point(data = pts, aes(perf, ykey, colour = model), size = 0.9) +
  scale_colour_manual(values = PAL, name = NULL) +
  scale_y_discrete(labels = ylabs) +
  facet_wrap(~ metric, scales = "free", ncol = 2) +
  labs(title = NULL, subtitle = NULL,
       x = "Cross-validated out-of-fold performance  (R² for continuous | AUC for binary)", y = NULL) +
  theme_heap(base_size = 8) +
  theme(panel.grid.major.y = element_blank(), legend.position = "top",
        axis.text.y = element_text(size = 4.4),
        axis.text.x = element_text(size = 5.5),
        axis.title = element_text(size = 7),
        strip.text = element_text(size = 7),
        legend.text = element_text(size = 6.5),
        panel.spacing.x = grid::unit(6, "pt"))

# COMPREHENSIVE REFERENCE figure: all exposures kept (split into a continuous and
# a binary facet, side by side, so the row count is max(n_cont, n_bin) not their
# sum). Authored at 6.5in = the supplement's \textwidth so LaTeX places it at
# scale 1.0 -- it used to be emitted 11in wide into main/, which meant the
# manuscript scaled it down by ~40% and the exposure labels came out unreadable.
# Height carries the rows instead of width: ~0.09in/row keeps 5.5pt labels legible
# while staying inside the supplement's page box.
n_rows_max <- max(n_cont, n_bin)
fig_h <- min(8.6, max(6.5, 0.09 * n_rows_max + 1.4))
heap_emit_figure(p, figure_id,
                 data = w[, .(exposure_id, exposure_label, exposure_type, category,
                              metric, cov_only, prot_plus_cov, gain)],
                 category = "supplement", subdir = "module6",
                 formats = c("pdf", "png"),
                 width = 6.5, height = fig_h, website = TRUE)

message("fig_pes_incremental_value: done (", covarType, "; ", n_cont, " continuous + ", n_bin, " binary).")
