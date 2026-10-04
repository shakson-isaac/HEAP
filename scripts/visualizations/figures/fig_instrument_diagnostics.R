#!/usr/bin/env Rscript

# ============================================================================
# fig_instrument_diagnostics.R  [figure_id: fig_instrument_diagnostics]
# ----------------------------------------------------------------------------
# CONSOLIDATED MR-instrument diagnostics for the exposure GWAS. Merges:
#   fig_gwas_h2_loci   (SNP h2 vs # independent lead loci)
#   fig_gwas_summary   (# loci vs genomic inflation lambda_GC)
#   fig_ldsc_intercept (lambda_GC vs LDSC intercept = confounding component)
#   fig_ldsc_h2        (per-exposure SNP heritability ranking -- absorbed into a)
#
# These were four scatter plots over the SAME four per-exposure quantities
# (h2, #loci, lambda_GC, LDSC intercept), shown two at a time. Here they are one
# figure, read left to right as the actual instrument-selection logic:
#   a  is there signal?      more heritable exposures map more independent loci
#   b  how much signal?      loci count vs genomic inflation
#   c  is the signal REAL?   LDSC splits inflation into polygenicity vs
#                            confounding -- the intercept stays near 1, so the
#                            inflation is polygenic signal, not structure
#
# Supports the results sentence: "Instrument strength, SNP heritability, genomic
# calibration, and genetic-correlation structure varied and informed which
# exposures were suitable for MR" (genetic-correlation structure = fig_ldsc_rg,
# which stays its own figure).
#
# Input : output/gwas/gwas_locus_summary.tsv    via load_gwas_locus_summary()
#         output/gwas/ldsc/ldsc_h2_summary.tsv  via load_ldsc_h2()
# Output: figures/supplement/gwas/fig_instrument_diagnostics.{pdf,png} + data tsv
#
# Authored at 6.5in = the supplement's \textwidth (scale 1.0).
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
  library(data.table); library(ggplot2); library(ggrepel); library(patchwork)
})

figure_id <- Sys.getenv("HEAP_FIGURE_ID", unset = "fig_instrument_diagnostics")
BS <- 7.5

# ---------------------------------------------------------------- data ------
loci <- load_gwas_locus_summary()[, .(exposure, n_lead, n_gwsig, lambda_gc,
                                      top_locus, top_p)]
h2   <- load_ldsc_h2()[, .(exposure, category, h2, h2_se, h2_z,
                           intercept, intercept_se, mean_chi2, ratio)]
dt   <- merge(loci, h2, by = "exposure")
dt   <- dt[is.finite(h2) & is.finite(n_lead)]
if (!nrow(dt)) stop("No exposures with both GWAS loci and LDSC h2.")
dt[, Category := heap_category_factor(category)]
dt[, label := heap_exposure_label(exposure)]

rho_a <- suppressWarnings(cor(dt$h2, dt$n_lead, method = "spearman", use = "complete.obs"))
rho_b <- suppressWarnings(cor(dt$n_lead, dt$lambda_gc, method = "spearman", use = "complete.obs"))
med_int <- median(dt$intercept, na.rm = TRUE)
med_rat <- median(dt$ratio[is.finite(dt$ratio)], na.rm = TRUE)
message(sprintf("instruments: %d exposures | rho(h2,loci)=%.2f rho(loci,lambda)=%.2f | median LDSC intercept=%.3f (ratio %.0f%%)",
                nrow(dt), rho_a, rho_b, med_int, 100 * med_rat))

# NB: annotate() at x=-Inf/y=Inf is silently DROPPED on a sqrt scale (sqrt(-Inf)
# = NaN), so anchor the annotation at real data coordinates instead.
ann <- function(txt, x, y) annotate("text", x = x, y = y, label = txt,
                                    hjust = 0, vjust = 1, size = 2.1, colour = "grey25")

# ---- a: heritability -> loci ------------------------------------------------
lab_a <- dt[n_lead >= 40 | h2 >= 0.09][order(-n_lead)]
pa <- ggplot(dt, aes(h2, n_lead, colour = Category)) +
  geom_point(size = 1.2, alpha = .85) +
  geom_text_repel(data = lab_a, aes(label = label), size = 1.8, max.overlaps = Inf,
                  min.segment.length = 0, box.padding = .45, force = 4,
                  segment.size = .18, segment.colour = "grey70", seed = 1,
                  show.legend = FALSE) +
  ann(sprintf("Spearman rho = %.2f", rho_a), x = min(dt$h2), y = max(dt$n_lead)) +
  scale_y_sqrt(breaks = c(0, 1, 5, 10, 25, 50, 100, 200),
               expand = expansion(mult = c(.02, .10))) +
  scale_x_continuous(expand = expansion(mult = c(.02, .06))) +
  scale_colour_exposure(drop = TRUE) +
  labs(x = expression(SNP~heritability~italic(h)[SNP]^2~"(LDSC, observed scale)"),
       y = "Independent lead loci\n(sqrt scale)") +
  theme_heap(base_size = BS) +
  theme(panel.grid = element_blank(), plot.margin = margin(10, 4, 2, 2))

# ---- b: loci -> genomic inflation -------------------------------------------
lab_b <- dt[n_lead >= 60 | lambda_gc >= 1.45][order(-n_lead)]
pb <- ggplot(dt, aes(n_lead, lambda_gc, colour = Category)) +
  geom_point(size = 1.2, alpha = .85) +
  geom_text_repel(data = lab_b, aes(label = label), size = 1.8, max.overlaps = Inf,
                  min.segment.length = 0, box.padding = .45, force = 4,
                  segment.size = .18, segment.colour = "grey70", seed = 2,
                  show.legend = FALSE) +
  ann(sprintf("Spearman rho = %.2f", rho_b), x = 0, y = max(dt$lambda_gc)) +
  scale_x_sqrt(breaks = c(0, 1, 5, 10, 25, 50, 100, 200),
               expand = expansion(mult = c(.02, .08))) +
  scale_colour_exposure(drop = TRUE) +
  labs(x = "Independent genome-wide-significant lead loci (sqrt scale)",
       y = expression("Genomic inflation "*lambda[GC])) +
  theme_heap(base_size = BS) +
  theme(panel.grid = element_blank(), plot.margin = margin(10, 4, 2, 2))

# ---- c: inflation = polygenicity or confounding? ----------------------------
dc <- dt[is.finite(intercept) & is.finite(lambda_gc)]
lab_c <- dc[lambda_gc >= 1.5 | ratio >= 0.6][order(-lambda_gc)][seq_len(min(10, .N))]
pc <- ggplot(dc, aes(lambda_gc, intercept)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              colour = "#C0392B", linewidth = .35) +
  geom_hline(yintercept = 1, colour = "grey55", linewidth = .3) +
  geom_linerange(aes(ymin = intercept - intercept_se, ymax = intercept + intercept_se),
                 colour = "grey80", linewidth = .2) +
  geom_point(aes(fill = ratio), shape = 21, size = 1.5, stroke = .15, colour = "grey30") +
  geom_text_repel(data = lab_c, aes(label = label), size = 1.8, max.overlaps = Inf,
                  min.segment.length = 0, box.padding = .45,
                  segment.size = .18, segment.colour = "grey70", seed = 3) +
  scale_fill_viridis_c(option = "plasma", direction = -1, limits = c(0, 1),
                       name = "Confounding\nratio") +
  labs(x = expression("Genomic inflation "*lambda[GC]*" (median "*chi^2*")"),
       y = "LDSC intercept\n(confounding component)") +
  theme_heap(base_size = BS) +
  theme(panel.grid = element_blank(), plot.margin = margin(10, 4, 2, 2))

# ------------------------------------------------------------- assemble ------
# guides="collect" gathers the shared exposure-category legend (a, b) and the
# confounding colorbar (c) into ONE right-hand column, so all three panels keep
# the same plotting width.
p <- (pa / pb / pc) + plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(legend.position = "right",
        legend.key.size = unit(7, "pt"),
        legend.text  = element_text(size = BS - 2),
        legend.title = element_text(size = BS - 1),
        plot.tag = element_text(face = "bold", size = 9),
        plot.tag.position = c(0, 1))

out <- dt[order(-n_lead),
          .(exposure, label = as.character(label), category,
            n_lead, n_gwsig, h2, h2_se, h2_z, lambda_gc, mean_chi2,
            intercept, intercept_se, ratio, top_locus, top_p)]

heap_emit_figure(p, figure_id, data = out, category = "supplement", subdir = "gwas",
                 formats = c("pdf", "png"), width = 6.5, height = 8.0, website = TRUE)

message("fig_instrument_diagnostics: done.")
