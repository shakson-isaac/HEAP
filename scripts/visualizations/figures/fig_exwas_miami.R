#!/usr/bin/env Rscript

# ============================================================================
# fig_exwas_miami.R  [figure_id: fig_exwas_miami]   == Supplement Fig S6 ==
# ----------------------------------------------------------------------------
# "Entire Proteome ExWAS" — signed Miami/Manhattan plot of every exposure-term x
# protein association, faithful to the legacy Module2 miami plot
# (Visualizations/Module2/HEAPassoc_main.R). Demoted from the manuscript main
# Fig 3A to a supplement panel; the descriptive summary (fig_assoc_summary) is
# the main-text Figure 3.
#
#   x = exposure terms, ordered by exposure Category (ticks hidden)
#   y = -log10(p_train + 1e-300) * sign(beta_train)   [SIGNED significance]
#   grey      = not replicated across the 80/20 train/test split
#   colored  = REPLICATED (Bonferroni-significant in BOTH train AND test),
#               colored by exposure Category
#   blue dashed lines = +/- Bonferroni threshold (manuscript p < 7e-8)
#
# Because only sign(beta) enters the y-axis, the ordered-factor polynomial
# contrast terms (e.g. deprivation indices) whose raw coefficients explode are
# harmless here — and they are grey anyway, since they do not replicate.
#
# Input : module2/<experiment>/<covarType>/univar_assoc_*.rds (statE, train+test)
#         via load_module2_replicated()
# Output: figures/main/fig_exwas_miami.{pdf,png} + figures/data/...tsv
#
# Run:
#   HEAP_PATHS_FILE=.../workflow/00_paths.R \
#     Rscript scripts/visualizations/figures/fig_exwas_miami.R [covarType] [experiment]
# Defaults to the base spec (M2_base_main/base). The manuscript uses the maximal
# specification; rerun with that experiment once it has finished.
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
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })

a <- commandArgs(trailingOnly = TRUE); a <- a[!startsWith(a, "--")]
a <- a[!a %in% c("fig_exwas_miami", "all_main", "all_supplement", "all", "website")]
covarType  <- if (length(a) >= 1) a[1] else "base"
experiment <- if (length(a) >= 2) a[2] else "M2_base_main"
figure_id <- Sys.getenv("HEAP_FIGURE_ID", unset = "fig_exwas_miami")

# --- merged train/test statE with replication flag --------------------------
m <- load_module2_replicated(covarType = covarType, experiment = experiment)
thr <- attr(m, "pval_thresh")
m <- m[is.finite(p_train) & is.finite(beta_train)]

# Order x by Category (canonical broad-group order), then exposure term within.
m[, Category := heap_category_factor(Category)]
setorder(m, Category, ID)
m[, idx := as.integer(factor(ID, levels = unique(ID)))]
m[, signed := -log10(p_train + 1e-300) * sign(beta_train)]

# per-category x positions: block midpoint (label) + right boundary (separator)
xax <- m[, .(mid = mean(range(idx)), right = max(idx) + 0.5), by = Category]
setorder(xax, mid)

n_rep <- sum(m$replicated)
message(sprintf("ExWAS: %d term-protein points | replicated=%d | %d exposure terms | thr=%.2g",
                nrow(m), n_rep, max(m$idx), thr))

# --- plot: grey base, colored replicated, +/- Bonferroni lines --------------
# category blocks labeled along x and separated by faint rules, colored by the
# canonical exposure-category palette (common/plot_theme.R::scale_colour_exposure).
p <- ggplot(m, aes(idx, signed)) +
  geom_vline(xintercept = head(xax$right, -1), colour = "grey90", linewidth = 0.2) +
  geom_point(data = m[replicated == FALSE], colour = "grey80", size = 0.5, alpha = 0.5) +
  geom_point(data = m[replicated == TRUE], aes(colour = Category), size = 0.7, alpha = 0.85) +
  geom_hline(yintercept =  -log10(thr + 1e-300), linetype = "dashed", colour = "blue") +
  geom_hline(yintercept = log10(thr + 1e-300), linetype = "dashed", colour = "blue") +
  scale_colour_exposure(drop = FALSE) +
  scale_x_continuous(breaks = xax$mid, labels = heap_category_pretty(xax$Category),
                     expand = expansion(mult = 0.01)) +
  labs(title = "Entire-proteome ExWAS",
       subtitle = NULL,
       x = NULL,
       y = expression(-log[10]~"(P) "%*%" sign"~(beta))) +
  theme_heap() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank(), legend.position = "right") +
  guides(colour = guide_legend(override.aes = list(size = 2)))

out <- m[, .(exposure_id = ID, protein = omicID, Category, idx,
             beta_train, p_train, signed, replicated)]
# ~10^6 points: rasterize the point clouds, keep text/axes/threshold lines vector
p <- heap_rasterize(p)
heap_emit_figure(p, figure_id, data = out, category = "supplement",
                 formats = c("pdf", "png"), width = 12, height = 3.4, website = TRUE)
message("fig_exwas_miami: done (", experiment, "/", covarType, ").")
