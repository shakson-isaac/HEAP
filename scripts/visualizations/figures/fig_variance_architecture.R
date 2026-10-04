#!/usr/bin/env Rscript

# ============================================================================
# fig_variance_architecture.R  [figure_id: fig_variance_architecture]
# ----------------------------------------------------------------------------
# CONSOLIDATED Module-1 variance-partition supplement figure. Replaces five
# figures that each showed one slice of the same partition:
#   fig_variance_raincloud   (ridgeline of Total/Cov/G/E/GxE)
#   fig_genetic_subblocks    (cis vs trans violin)
#   fig_covariate_vs_pgs_r2  (covariate vs genetic violin -- same plot, other cols)
#   fig_r2_metrics           (G/E/GxE jitter -- duplicate of the ridgeline)
#   fig_total_r2_density     (whole-model density -- one row of the ridgeline)
#
# Three panels, three genuinely different questions:
#   a  How big is each component, per protein?   -> violin+box, log10, nonzero
#      only, with the zero fraction stated (40-63% of proteins take exactly 0
#      for the genetic/exposomic blocks, so a median over ALL proteins is 0 for
#      cis/trans and hides everything).
#   b  How FAR does each component reach?        -> % of proteins above an R2
#      threshold. Zeros are handled natively, and this is the panel that
#      actually justifies the ordering claim: the components separate in the
#      tail, not at the median.
#   c  How much do we explain at all?            -> whole-model R2 density.
#
# Authored at 6.5in wide = the supplement's \textwidth (469.76pt) so LaTeX
# places it at scale 1.0 (WYSIWYG; no downscaled text).
#
# Input : module1_predictive_r2_score_partition/<exp>/<covarType>/<method>/
#           predictive_r2_coarse_*  +  predictive_r2_genetic_subblocks_*
# Output: figures/supplement/module1/fig_variance_architecture.{pdf,png} + data tsv
#
# Run:
#   HEAP_PATHS_FILE=.../workflow/00_paths.R \
#     Rscript scripts/visualizations/figures/fig_variance_architecture.R [covarType] [method]
# ============================================================================

local({
  cand <- c(file.path(getwd(), "scripts", "visualizations", "common"),
            "/n/groups/patel/shakson_ukb/HEAP/scripts/visualizations/common")
  common <- cand[dir.exists(cand)][1]
  if (is.na(common)) stop("cannot locate common/ helpers")
  for (f in c("figure_paths", "load_heap_results", "plot_theme",
              "label_helpers", "export_helpers"))
    source(file.path(common, paste0(f, ".R")))
})
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork); library(scales)
})

a <- commandArgs(trailingOnly = TRUE); a <- a[!startsWith(a, "--")]
a <- a[!a %in% c("fig_variance_architecture", "all_main", "all_supplement", "all", "website")]
covarType <- if (length(a) >= 1) a[1] else "base"
method    <- if (length(a) >= 2) a[2] else "lasso"
figure_id <- Sys.getenv("HEAP_FIGURE_ID", unset = "fig_variance_architecture")

BS <- 8   # base_size: at 6.5in wide + scale 1.0 this lands text at ~7-8 pt

# ---------------------------------------------------------------- data ------
coarse <- load_module1_predictive_r2(covarType = covarType, method = method, level = "coarse")
mt <- coarse[get("method") == "score_model_total", .(r2 = mean(r2)), by = .(omic, block)]
ud <- coarse[get("method") == "score_unique_drop", .(r2 = mean(r2)), by = .(omic, block)]

sub <- load_module1_predictive_r2(covarType = covarType, method = method,
                                  level = "genetic_subblocks")
gsu <- sub[get("method") == "score_unique_drop" & block %in% c("Gcis", "Gtrans"),
           .(r2 = mean(r2)), by = .(omic, block)]

total <- mt[block == "C+G+E+GxE", .(omic, r2)]

# NB: use score_unique_drop "Covars" -- NOT score_model_total "C". The latter is
# the covariate-ONLY model R2 and is not the same quantity as the unique G/E/GxE
# blocks; mixing them on one axis (as the old ridgeline did) is apples-to-oranges.
# In practice the two barely differ (median 0.0190 vs 0.0198) because covariates
# are near-orthogonal to G/E/GxE -- but the axis is now one quantity throughout.
d <- rbindlist(list(
  ud[block == "Covars",   .(omic, comp = "Covariates", r2)],
  ud[block == "G",        .(omic, comp = "Genetic",    r2)],
  gsu[block == "Gcis",    .(omic, comp = "cis",        r2)],
  gsu[block == "Gtrans",  .(omic, comp = "trans",      r2)],
  ud[block == "E",        .(omic, comp = "Exposomic",  r2)],
  ud[block == "GxE",      .(omic, comp = "GxE",        r2)]))

LEV <- c("Covariates", "Genetic", "cis", "trans", "Exposomic", "GxE")
d[, comp := factor(comp, levels = LEV)]
NP <- uniqueN(d$omic)

# cis/trans are SUBDIVISIONS of Genetic -> lighter blue, and flagged as such
BLU <- HEAP_PAL_COMPONENT[["Genetic"]]
PAL <- c("Covariates" = HEAP_PAL_COMPONENT[["Covars"]],
         "Genetic"    = BLU,
         "cis"        = "#5C93BE",
         "trans"      = "#9CC0DA",
         "Exposomic"  = HEAP_PAL_COMPONENT[["Exposome"]],
         "GxE"        = HEAP_PAL_COMPONENT[["GxE"]])

stats <- d[, .(pct0    = 100 * mean(r2 <= 0),
               med_nz  = median(r2[r2 > 0]),
               med_all = median(r2),
               p95_nz  = quantile(r2[r2 > 0], .95),
               max     = max(r2)), by = comp]

# ------------------------------------------------- a: per-component violins --
# nonzero only + log10: the only rendering in which cis/trans are not a flat
# line at zero. The zero fraction is printed rather than silently dropped.
FLOOR <- 1e-4                            # low enough that clamping is negligible
nz <- d[r2 > 0]
nz[, r2p := pmax(r2, FLOOR)]            # clamp (not drop) the sub-floor tail
n_clamped <- nz[r2 < FLOOR, .N]

pa <- ggplot(nz, aes(comp, r2p, fill = comp)) +
  geom_violin(width = .88, alpha = .45, colour = "grey30",
              linewidth = .25, scale = "width") +
  geom_boxplot(width = .15, outlier.shape = NA, colour = "grey15", linewidth = .3) +
  geom_point(data = stats, aes(comp, p95_nz), inherit.aes = FALSE,
             shape = 18, size = 1.7, colour = "grey10") +
  scale_y_log10(breaks = c(.0001, .001, .01, .1, .7),
                labels = label_number(drop0trailing = TRUE),
                limits = c(FLOOR, 3.2)) +
  scale_fill_manual(values = PAL, guide = "none") +
  coord_cartesian(clip = "off") +
  labs(x = NULL, y = expression("Unique predictive"~R^2~"(drop-one)")) +
  theme_heap(base_size = BS) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(size = BS - 0.5),
        plot.margin = margin(10, 5, 2, 4))

# ------------------------------------------------------------- b: reach ------
grid  <- 10^seq(log10(0.002), log10(0.7), length.out = 300)
reach <- d[, .(x = grid, frac = sapply(grid, function(t) 100 * mean(r2 >= t))), by = comp]
LTY <- c("Covariates" = "solid", "Genetic" = "solid", "cis" = "dashed",
         "trans" = "dotted", "Exposomic" = "solid", "GxE" = "solid")
LWD <- c("Covariates" = .85, "Genetic" = .85, "cis" = .4,
         "trans" = .4, "Exposomic" = .85, "GxE" = .85)

pb <- ggplot(reach, aes(x, frac, colour = comp, linetype = comp, linewidth = comp)) +
  geom_vline(xintercept = 0.01, linetype = "dashed", colour = "grey70", linewidth = .25) +
  geom_line() +
  annotate("text", x = 0.0107, y = 79, label = "R² = 0.01", hjust = 0,
           size = 2.2, colour = "grey45") +
  scale_x_log10(breaks = c(.002, .01, .05, .2, .5),
                labels = label_number(drop0trailing = TRUE)) +
  scale_colour_manual(values = PAL, name = NULL) +
  scale_linetype_manual(values = LTY, name = NULL) +
  scale_linewidth_manual(values = LWD, name = NULL) +
  labs(x = expression("Component"~R^2~"threshold"),
       y = "% of proteins above threshold") +
  theme_heap(base_size = BS) +
  theme(panel.grid = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(0.99, 0.98),
        legend.justification = c(1, 1),
        legend.background = element_rect(fill = alpha("white", 0.8), colour = NA),
        legend.key = element_blank(),
        legend.key.width = unit(12, "pt"),
        legend.key.height = unit(8, "pt"),
        legend.text = element_text(size = BS - 1.5),
        plot.margin = margin(4, 4, 2, 4))

# ----------------------------------------------------- c: whole-model R2 -----
med_tot <- median(total$r2)
pc <- ggplot(total, aes(r2)) +
  geom_density(fill = "#7FA8C9", colour = "grey20", alpha = .75, linewidth = .35) +
  geom_vline(xintercept = med_tot, linetype = "dashed", colour = "grey25", linewidth = .35) +
  annotate("text", x = med_tot + .012, y = 5.2,
           label = sprintf("median %.3f", med_tot), hjust = 0, size = 2.3, colour = "grey25") +
  scale_x_continuous(limits = c(-0.03, NA)) +
  labs(x = expression("Whole-model test"~R^2), y = "Density") +
  theme_heap(base_size = BS) +
  theme(panel.grid = element_blank(), plot.margin = margin(10, 4, 2, 4))

# ------------------------------------------------------------- assemble ------
p <- (pa / (pb | pc)) +
  plot_layout(heights = c(0.95, 1)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 9),
        plot.tag.position = c(0, 1))

# emitted data = the exact plotted values (long form) + the per-component summary
out_data <- rbindlist(list(
  d[, .(panel = "a_b", comp = as.character(comp), omic, r2)],
  total[, .(panel = "c", comp = "Total", omic, r2)]))

# subdir = "module1" so this lands beside its Module-1 siblings in
# figures/supplement/module1/ -- the manifest's source path must match the disk
# path or the manuscript silently keeps rendering a stale figure.
heap_emit_figure(p, figure_id, data = out_data, category = "supplement",
                 subdir = "module1",
                 formats = c("pdf", "png"), width = 6.5, height = 7.0, website = TRUE)

message(sprintf("fig_variance_architecture: done (%s/%s; %d proteins; %d nonzero values clamped to floor %.4f).",
                covarType, method, NP, n_clamped, FLOOR))
print(stats)
