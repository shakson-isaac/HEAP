#!/usr/bin/env Rscript

# ============================================================================
# fig_greml_vs_r2.R  [figure_id: fig_greml_vs_r2]
# ----------------------------------------------------------------------------
# Cross-method concordance of the Module-1 G/E/GxE decomposition: per-protein
# GREML variance component (multi-kernel REML, V/Vp) vs HEAP unique drop-one
# predictive R2 (out-of-fold). Three facets (G/E/GxE), y=x ceiling line, lm
# trend; Spearman rho in each strip. Supersedes the legacy fig_h2_vs_r2
# (which used external Sun et al. h2 + pre-refactor data) with the internal GREML.
#
# SINGLE PANEL: component-level concordance only. The former panel b (a COUNT of
# significant GSEA gene sets per ranking) is retired here; the actual per-component
# pathway enrichments live in the dedicated figure fig_varcomp_pathways.
#
# Input : population_architecture/<covar>/grm_cutoff_<cut>/
#           concordance_greml_vs_heap_{r2,stats}.tsv  (built by
#           scripts/support/module1_greml_vs_r2.R -- NO stats computed here)
# Output: figures/supplement/module1/fig_greml_vs_r2.{pdf,png} + data tsv
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
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork) })

covarType <- Sys.getenv("PA_COVAR",      "base")
grm_cut   <- Sys.getenv("PA_GRM_CUTOFF", "0p025")
figure_id <- Sys.getenv("HEAP_FIGURE_ID", unset = "fig_greml_vs_r2")
CELL <- nzchar(Sys.getenv("HEAP_CELL"))   # composite rigor cell = the EXPOSOMIC concordance only
BS  <- if (CELL) 7 else 11
B_W <- 3.4; B_H <- 2.6            # single exposomic concordance panel (G & GxE -> supplement)

pa_dir <- file.path(heap_project_output("population_architecture"), covarType,
                    paste0("grm_cutoff_", grm_cut))
tab   <- fread(file.path(pa_dir, "concordance_greml_vs_heap_r2.tsv"))
stats <- fread(file.path(pa_dir, "concordance_greml_vs_heap_stats.tsv"))

lev <- c("G", "E", "GxE")
labmap <- setNames(
  if (CELL) sprintf("%s  rho=%.2f", c("Genetic","Exposomic","GxE")[match(stats$component, lev)], stats$spearman)
  else sprintf("%s  (rho = %.2f, n = %s)",
               c("Genetic","Exposomic","GxE")[match(stats$component, lev)],
               stats$spearman, formatC(stats$n, big.mark = ",")),
  stats$component)
tab[, component := factor(component, levels = lev)]
tab[, facet := factor(labmap[as.character(component)], levels = labmap[lev])]
tab[, r2_plot := pmax(0, r2)]   # floor OOF negatives for the ceiling view (raw kept in tsv)

pal <- c(G = HEAP_PAL_COMPONENT[["Genetic"]], E = HEAP_PAL_COMPONENT[["Exposome"]],
         GxE = HEAP_PAL_COMPONENT[["GxE"]])

# ---- panel a: per-protein concordance (GREML V/Vp vs HEAP predictive R2) ----
p_scatter <- ggplot(tab, aes(greml, r2_plot)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey55") +
  geom_point(aes(colour = component), alpha = 0.30, size = 0.7, stroke = 0) +
  geom_smooth(method = "lm", se = FALSE, colour = "grey20", linewidth = 0.55) +
  scale_colour_manual(values = pal, guide = "none") +
  facet_wrap(~ facet, scales = "free") +
  # encoding key ("dashed = y=x ceiling") belongs in the caption, not baked on the plot
  labs(title = if (CELL) "Estimates replicate by GREML" else "GREML vs HEAP variance components",
       x = if (CELL) "GREML V/Vp" else "GREML variance component (V/Vp)",
       y = if (CELL) expression("HEAP unique"~R^2) else expression("HEAP unique predictive"~R^2)) +
  theme_heap(base_size = BS) +
  theme(panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold", size = rel(if (CELL) 0.8 else 1)),
        plot.title = element_text(face = "bold", size = if (CELL) 8 else rel(1.0), hjust = 0.5),
        axis.title = element_text(face = if (CELL) "plain" else "bold"))

if (CELL) {
  # main-figure rigor = the EXPOSOMIC concordance only (G & GxE go to the supplement)
  e <- tab[component == "E"]; rhoE <- stats[component == "E", spearman]
  pe <- ggplot(e, aes(greml, r2_plot)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey55") +
    geom_point(colour = HEAP_PAL_COMPONENT[["Exposome"]], alpha = 0.35, size = 0.8, stroke = 0) +
    geom_smooth(method = "lm", se = FALSE, colour = "grey20", linewidth = 0.6) +
    annotate("text", x = 0, y = max(e$r2_plot, na.rm = TRUE) * 0.96, hjust = 0,
             label = sprintf("Spearman rho = %.2f", rhoE), size = 2.6, colour = "grey25") +
    labs(title = "Exposomic estimates replicate (GREML)",
         x = "GREML variance component (V/Vp)", y = expression("HEAP unique"~R^2)) +
    theme_heap(base_size = BS) +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold", size = 8, hjust = 0.5),
          plot.title.position = "panel", axis.title = element_text(face = "plain"),
          plot.margin = margin(3, 4, 2, 3))
  FIGDIR_CELL <- file.path(Sys.getenv("HEAP_PROJECT_ROOT", "/n/groups/patel/IGLOO/UKB/HEAP"), "figures", "exploratory/module1")
  dir.create(FIGDIR_CELL, recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path(FIGDIR_CELL, paste0(figure_id, "_cell.png")), pe, width = B_W, height = B_H, dpi = 400, bg = "white")
  ggsave(file.path(FIGDIR_CELL, paste0(figure_id, "_cell.pdf")), pe, width = B_W, height = B_H, bg = "white")
  message("fig_greml_vs_r2 CELL (exposomic concordance) done"); quit(save = "no")
}

# ---- single-panel figure: per-protein component-level concordance only -------
# (the former pathway-concordance bar is now its own figure, fig_varcomp_pathways)
# ONE short centred title, matching the other supplementary figures. The banner
# title plus the bold left-aligned "Per-protein concordance" line under it read
# as two competing headings (author, 2026-09-01).
p <- p_scatter

heap_emit_figure(p, figure_id, data = tab[, .(protein, component, greml, greml_se, r2, r2_plot)],
                 category = "supplement", formats = c("pdf", "png"),
                 width = 11, height = 3.7, website = FALSE)
message("fig_greml_vs_r2: done (", covarType, "/grm_", grm_cut,
        "; G rho=", round(stats[component=="G"]$spearman,2),
        ", E rho=", round(stats[component=="E"]$spearman,2),
        ", GxE rho=", round(stats[component=="GxE"]$spearman,2), ").")
