#!/usr/bin/env Rscript

# ============================================================================
# fig_ldsc_rg.R  [figure_id: fig_ldsc_rg]
# ----------------------------------------------------------------------------
# Genetic correlation (rg) between the exposures, from bivariate LD Score
# Regression on the REGENIE step-2 exposure GWAS. A symmetric heatmap of the
# pairwise rg among the well-powered exposures (those with a reliable SNP h2),
# ordered by hierarchical clustering so genetically-related exposures sit in
# blocks. This is the panel that shows which exposures share common-variant
# genetic architecture (e.g. smoking/alcohol behaviours, diet groups).
#
#   tile fill = rg (genetic correlation), diverging blue-white-red, [-1, 1]
#   *         = pair survives FDR (BH q < 0.05) for rg != 0
#   diagonal  = 1 (self), drawn grey
#
# rg is noisy for low-h2 traits, so the exposure set is restricted upstream
# (slurm/ldsc/ldsc_rg_exposures.txt = well-powered exposures) before rg is run;
# this plotter just draws whatever pairs are in the collected table.
#
# Input : output/gwas/ldsc/ldsc_rg_summary.tsv via load_ldsc_rg()
#         (built by the LDSC rg stage + scripts/ldsc/collect_ldsc_rg.R)
# Output: figures/supplement/gwas/fig_ldsc_rg.{pdf,png} + figures/data/...tsv
#
# Run:
#   HEAP_PATHS_FILE=.../workflow/00_paths.R \
#     Rscript scripts/visualizations/figures/fig_ldsc_rg.R
# ============================================================================

local({
  cand <- c(file.path(getwd(), "scripts", "visualizations", "common"),
            "/n/groups/patel/shakson_ukb/HEAP/scripts/visualizations/common")
  common <- cand[dir.exists(cand)][1]
  if (is.na(common)) stop("fig_ldsc_rg.R: cannot locate common/ helpers")
  for (f in c("figure_paths", "load_heap_results", "plot_theme",
              "label_helpers", "export_helpers"))
    source(file.path(common, paste0(f, ".R")))
})
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })

figure_id <- Sys.getenv("HEAP_FIGURE_ID", unset = "fig_ldsc_rg")

edges <- load_ldsc_rg()
edges <- edges[is.finite(rg)]
if (!nrow(edges)) stop("No LDSC rg estimates available — run the LDSC rg stage first.")

# Clamp out-of-range rg (LDSC can return |rg| slightly > 1 for noisy pairs).
edges[, rg_c := pmin(pmax(rg, -1), 1)]
# BH FDR over the off-diagonal pairs (display-level multiple-testing flag).
edges[, q := p.adjust(p, method = "BH")]
edges[, sig := is.finite(q) & q < 0.05]

exps <- sort(unique(c(edges$p1, edges$p2)))
# tiles are keyed on the UNIQUE exposure id (short labels can collide, which would
# duplicate factor levels); labels are mapped onto the axes via scale_*_discrete.
lab  <- heap_exposure_label(exps); names(lab) <- exps

# --- cluster order: hclust on 1 - rg (NA/unestimated -> 0) -------------------
m <- matrix(0, length(exps), length(exps), dimnames = list(exps, exps))
diag(m) <- 1
for (i in seq_len(nrow(edges))) { a <- edges$p1[i]; b <- edges$p2[i]
  m[a, b] <- edges$rg_c[i]; m[b, a] <- edges$rg_c[i] }
ord <- tryCatch(hclust(as.dist(1 - m), method = "average")$order,
                error = function(e) seq_along(exps))
lev <- exps[ord]   # exposure ids in cluster order

# --- long table for both triangles + diagonal (keyed on exposure id) --------
sym <- rbind(
  edges[, .(a = p1, b = p2, rg_c, sig)],
  edges[, .(a = p2, b = p1, rg_c, sig)],
  data.table(a = exps, b = exps, rg_c = 1, sig = FALSE))
sym[, `:=`(a = factor(a, levels = lev),
           b = factor(b, levels = lev))]

n <- length(exps)
n_sig <- sum(edges$sig, na.rm = TRUE)
p <- ggplot(sym, aes(a, b, fill = rg_c)) +
  geom_tile(colour = "grey92", linewidth = 0.15) +
  geom_point(data = sym[sig == TRUE], shape = 8, size = 0.7, colour = "grey15",
             show.legend = FALSE) +
  scale_fill_gradient2(name = expression(italic(r)[g]), low = "#2166AC",
                       mid = "white", high = "#B2182B", midpoint = 0,
                       limits = c(-1, 1), na.value = "grey88") +
  scale_x_discrete(labels = lab) +
  scale_y_discrete(labels = lab) +
  coord_equal() +
  labs(title = "Genetic correlation between exposures (LDSC)",
       subtitle = NULL,
       x = NULL, y = NULL) +
  theme_heap() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 6),
        axis.text.y = element_text(size = 6),
        panel.grid = element_blank())

out <- edges[order(-abs(rg)),
             .(p1, p2, label1 = lab[p1], label2 = lab[p2],
               rg, se, z, p, q, sig)]
sz <- max(7, 0.20 * n + 2)
heap_emit_figure(p, figure_id, data = out, category = "supplement",
                 formats = c("pdf", "png"), width = sz, height = sz, website = TRUE)
message(sprintf("fig_ldsc_rg: done (%d exposures, %d pairs, %d FDR<0.05).",
                n, nrow(edges), n_sig))
