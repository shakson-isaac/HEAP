#!/usr/bin/env Rscript

# ============================================================================
# fig_pes_specificity.R  [figure_id: fig_pes_specificity]
# ----------------------------------------------------------------------------
# CONSOLIDATED (2026-07-11). Merges fig_pes_score_correlation + fig_pes_shared_proteins,
# which were cited by a catch-all sentence ("additional analyses ... are shown in")
# that made no claim at all. Together they answer a question the paper needs
# answered: does each PES read ITS OWN exposure, or do they all read one common
# health axis?
#
# They read their own exposure. That is the claim this figure makes.
#
#   a  PES scores are only weakly correlated with one another: median |r| = 0.17
#      across the 26,732 off-diagonal exposure pairs, and 57% of pairs fall below
#      |r| = 0.2. A common health axis would put this mass near |r| = 1.
#   b  and they are built from largely different proteins: the median protein
#      appears in just 3 of the 164 PES panels, 30% appear in exactly ONE, and NO
#      protein is used by even half the exposures.
#
# NB the predecessor fig_pes_shared_proteins plotted only the TOP 25 proteins
# (TOPN <- 25L) and titled itself "A shared proteomic backbone runs across exposure
# signatures". Ranking by sharedness and then showing only the most-shared tail
# guarantees that conclusion. The FULL distribution says the opposite, and is what
# is plotted here.
#
# Input : module6_pes_longitudinal/pes_score_correlation/pes_score_correlation_long_<cov>.tsv
#         module6_pes_longitudinal/deployable_lasso_weights/<cov>/**/*_lasso_k50_weights.txt
# Output: figures/supplement/module6/fig_pes_specificity.{pdf,png} + data tsv
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
  library(data.table); library(ggplot2); library(patchwork); library(scales)
})

a <- commandArgs(trailingOnly = TRUE); a <- a[!startsWith(a, "--")]
a <- a[!a %in% c("fig_pes_specificity", "all_main", "all_supplement", "all", "website")]
covarType <- if (length(a) >= 1) a[1] else "base"
figure_id <- Sys.getenv("HEAP_FIGURE_ID", unset = "fig_pes_specificity")
BS <- 7.5

# ---------------------------------------- a: are the SCORES redundant? -------
cor_dir <- heap_project_output("module6_pes_longitudinal", "pes_score_correlation")
ordf <- file.path(cor_dir, paste0("pes_score_correlation_order_", covarType, ".tsv"))
lngf <- file.path(cor_dir, paste0("pes_score_correlation_long_",  covarType, ".tsv"))
ord  <- fread(ordf)$exposure_id
long <- fread(lngf)[is.finite(r)]
idx  <- setNames(seq_along(ord), ord)
long[, `:=`(xi = idx[exposure_i], yi = idx[exposure_j])]

off <- long[exposure_i != exposure_j]
med_r  <- median(abs(off$r))
pct_lo <- 100 * mean(abs(off$r) < 0.2)
message(sprintf("PES score correlation: %s off-diagonal pairs | median |r| = %.2f | %.0f%% below 0.2",
                comma(nrow(off)), med_r, pct_lo))

pa <- ggplot(long, aes(xi, yi, fill = r)) +
  geom_raster() +
  scale_fill_gradient2(low = "#3a86ff", mid = "white", high = "#d62828",
                       midpoint = 0, limits = c(-1, 1), name = "r") +
  coord_equal(expand = FALSE) +
  labs(x = sprintf("%d exposures (ordered by cluster)", length(ord)), y = NULL) +
  theme_heap(base_size = BS) +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(), axis.ticks = element_blank(),
        legend.position = "right", legend.key.width = unit(5, "pt"),
        legend.key.height = unit(20, "pt"),
        legend.text = element_text(size = BS - 2),
        legend.title = element_text(size = BS - 1),
        plot.margin = margin(10, 4, 2, 2))

# marginal density of |r| makes the "weakly correlated" claim readable
pa2 <- ggplot(off, aes(abs(r))) +
  geom_histogram(bins = 40, fill = "#7FB3D5", colour = NA) +
  geom_vline(xintercept = med_r, linetype = "dashed", colour = "#C0392B", linewidth = .35) +
  annotate("text", x = med_r, y = Inf, hjust = -0.1, vjust = 1.6, size = 2.0, colour = "#C0392B",
           label = sprintf("median |r| = %.2f", med_r)) +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, .5, 1)) +
  scale_y_continuous(labels = comma) +
  labs(x = "|r| between two PES scores", y = "Exposure pairs") +
  theme_heap(base_size = BS) +
  theme(panel.grid = element_blank(), plot.margin = margin(10, 4, 2, 2))

# ------------------------------- b: are the PANELS built from the same proteins? --
wdir <- file.path(heap_project_output("module6_pes_longitudinal", "deployable_lasso_weights"),
                  covarType)
fs <- list.files(wdir, pattern = "_lasso_k50_weights\\.txt$", recursive = TRUE, full.names = TRUE)
if (!length(fs)) stop("no k50 weight files under ", wdir)
n_exp <- length(fs)
cnt <- rbindlist(lapply(fs, function(f) fread(f, select = "protein")), fill = TRUE)[
  , .(n_exposures = .N), by = protein][order(-n_exposures)]
cnt[, pct := 100 * n_exposures / n_exp]

med_n   <- median(cnt$n_exposures)
pct_one <- 100 * mean(cnt$n_exposures == 1)
n_half  <- cnt[pct >= 50, .N]
message(sprintf("PES panels: %d proteins over %d exposures | median in %d panels | %.0f%% in exactly 1 | %d in >=50%% of exposures",
                nrow(cnt), n_exp, med_n, pct_one, n_half))

pb <- ggplot(cnt, aes(n_exposures)) +
  geom_histogram(binwidth = 1, fill = "#82C09A", colour = NA) +
  geom_vline(xintercept = med_n, linetype = "dashed", colour = "#C0392B", linewidth = .35) +
  annotate("text", x = med_n, y = Inf, hjust = -0.15, vjust = 1.6, size = 2.0, colour = "#C0392B",
           label = sprintf("median = %d of %d panels", med_n, n_exp)) +
  scale_x_continuous(breaks = pretty_breaks(6)) +
  scale_y_continuous(labels = comma) +
  labs(x = sprintf("Number of the %d PES panels a protein appears in", n_exp),
       y = "Proteins") +
  theme_heap(base_size = BS) +
  theme(panel.grid = element_blank(), plot.margin = margin(10, 4, 2, 2))

p <- ((pa | pa2) / pb) +
  plot_layout(heights = c(1.15, 1)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 9),
        plot.tag.position = c(0, 1))

out <- rbindlist(list(
  off[, .(panel = "a", key = paste(exposure_i, exposure_j, sep = " | "), value = r)],
  cnt[, .(panel = "b", key = protein, value = as.numeric(n_exposures))]), use.names = TRUE)

heap_emit_figure(p, figure_id, data = out, category = "supplement", subdir = "module6",
                 formats = c("pdf", "png"), width = 6.5, height = 5.2, website = TRUE)

message("fig_pes_specificity: done.")
