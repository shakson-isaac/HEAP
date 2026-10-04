#!/usr/bin/env Rscript

# ============================================================================
# fig_pathways_by_component.R  [figure_id: fig_pathways_by_component]
# ----------------------------------------------------------------------------
# CONSOLIDATED per-component pathway enrichment. Replaces two figures that were
# the same lollipop-GSEA plot ranked two different ways:
#   fig_heap_partition_pathways  (proteins ranked by HEAP unique predictive R2)
#   fig_varcomp_pathways         (proteins ranked by the GREML variance component)
#
# Showing both rankings side by side in one 3x2 facet grid is a STRONGER result
# than either alone, because the GxE row is empty under BOTH rankings: GxE
# enriches zero gene sets no matter how you order the proteome. That is the
# "showed no enriched pathways" clause of the Module-1 GxE claim, and it is a
# negative that only reads as decisive when the two rankings agree.
#
# Input : docs/manuscript_stats/module1_enrichment/gsea_all_significant.tsv
#           (HEAP ranking; component = Exposomic/Genetic; 53 + 2 sets, GxE absent)
#         docs/manuscript_stats/module1_enrichment_greml/gsea_all_significant.tsv
#           (GREML ranking; metric == "greml"; component = E/G; 29 + 34, GxE absent)
# Output: figures/supplement/module1/fig_pathways_by_component.{pdf,png} + data tsv
#
# Authored at 6.5in = the supplement's \textwidth, so LaTeX places it at scale 1.0.
#
# Run:
#   HEAP_PATHS_FILE=.../workflow/00_paths.R \
#     Rscript scripts/visualizations/figures/fig_pathways_by_component.R
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
  library(data.table); library(ggplot2); library(scales)
})

figure_id <- Sys.getenv("HEAP_FIGURE_ID", unset = "fig_pathways_by_component")
TOP <- as.integer(Sys.getenv("HEAP_TOP_SETS", "10"))   # top gene sets per facet
BS  <- 8

ROOT <- Sys.getenv("HEAP_ROOT", "/n/groups/patel/shakson_ukb/HEAP")
f_heap  <- file.path(ROOT, "docs/manuscript_stats/module1_enrichment/gsea_all_significant.tsv")
f_greml <- file.path(ROOT, "docs/manuscript_stats/module1_enrichment_greml/gsea_all_significant.tsv")

# ---------------------------------------------------------------- data ------
R_HEAP  <- "HEAP R²"        # short strip labels -- long ones get clipped by the
R_GREML <- "GREML V/Vp"     # narrow facet panels (y-axis eats ~45% of the width)

H <- fread(f_heap)
H <- H[, .(ranking = R_HEAP, component, Description, setSize, NES, p.adjust)]

G <- fread(f_greml)
G <- G[metric == "greml"]
G[, component := c(G = "Genetic", E = "Exposomic", GxE = "GxE")[component]]
G <- G[, .(ranking = R_GREML, component, Description, setSize, NES, p.adjust)]

d <- rbindlist(list(H, G))
d <- d[!is.na(component) & component %in% c("Genetic", "Exposomic", "GxE")]

# de-duplicate identical gene-set names within a facet (GO/Reactome overlap)
setorder(d, ranking, component, p.adjust)
d <- unique(d, by = c("ranking", "component", "Description"))

# top N per facet by adjusted p
d[, rk := seq_len(.N), by = .(ranking, component)]
d <- d[rk <= TOP]

# GxE has ZERO enriched sets under either ranking -- that is the point of the
# figure, so its row must be PRESENT and visibly empty, not silently dropped.
LEVC <- c("Genetic", "Exposomic", "GxE")
LEVR <- c(R_HEAP, R_GREML)
d[, component := factor(component, levels = LEVC)]
d[, ranking   := factor(ranking,   levels = LEVR)]

# keep the long GO/Reactome names from eating the panel
d[, Description := ifelse(nchar(Description) > 44,
                          paste0(substr(Description, 1, 42), "…"), Description)]

empty <- CJ(ranking = factor(LEVR, levels = LEVR),
            component = factor(LEVC, levels = LEVC))
have  <- unique(d[, .(ranking, component)])
empty <- empty[!have, on = .(ranking, component)]
empty[, lab := "no enriched gene sets"]
XMID <- mean(range(d$NES, na.rm = TRUE))   # centre the annotation in the panel

n_gxe <- d[component == "GxE", .N]
stopifnot(n_gxe == 0)   # if GxE ever enriches something, this figure's claim changes

# order gene sets within each facet by NES
d[, ylab := factor(Description, levels = unique(Description[order(NES)])), by = .(ranking, component)]
d[, negl := -log10(p.adjust)]

PAL_C <- c(Genetic = HEAP_PAL_COMPONENT[["Genetic"]],
           Exposomic = HEAP_PAL_COMPONENT[["Exposome"]],
           GxE = HEAP_PAL_COMPONENT[["GxE"]])

p <- ggplot(d, aes(x = NES, y = ylab)) +
  geom_vline(xintercept = 0, colour = "grey65", linewidth = .3) +
  geom_segment(aes(x = 0, xend = NES, y = ylab, yend = ylab, colour = component),
               linewidth = .35, show.legend = FALSE) +
  geom_point(aes(size = setSize, fill = negl), shape = 21, colour = "grey25", stroke = .25) +
  geom_text(data = empty, aes(x = XMID, y = 1, label = lab), inherit.aes = FALSE,
            size = 2.2, colour = "grey45", fontface = "italic", hjust = 0.5) +
  scale_colour_manual(values = PAL_C, guide = "none") +
  scale_fill_distiller(palette = "YlOrRd", direction = 1,
                       name = expression(-log[10]~p[adj])) +
  scale_size_continuous(range = c(1.2, 4.2), name = "Gene set size") +
  facet_grid(component ~ ranking, scales = "free_y", space = "free_y") +
  labs(x = "Normalized enrichment score (NES)", y = NULL) +
  theme_heap(base_size = BS) +
  theme(panel.grid = element_blank(),
        strip.text.y = element_text(face = "bold", size = BS),
        strip.text.x = element_text(face = "bold", size = BS),
        axis.text.y  = element_text(size = BS - 2),
        legend.position = "right",
        legend.key.height = unit(10, "pt"),
        legend.text = element_text(size = BS - 2),
        legend.title = element_text(size = BS - 1))

heap_emit_figure(p, figure_id, data = d, category = "supplement", subdir = "module1",
                 formats = c("pdf", "png"), width = 6.5, height = 6.2, website = TRUE)

message(sprintf("fig_pathways_by_component: done (%d sets shown; GxE = %d sets under BOTH rankings).",
                nrow(d), n_gxe))
print(d[, .N, by = .(ranking, component)])
