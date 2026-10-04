#!/usr/bin/env Rscript

# ============================================================================
# fig_tissue_enrichment.R  [figure_id: fig_tissue_enrichment]
# ----------------------------------------------------------------------------
# Tissue-specificity GSEA of Module 2 exposure->protein associations: for each
# significant exposure, its protein t-value ranking was tested against GTEx
# tissue gene sets (clusterProfiler::GSEA). This figure summarises the flattened
# result table as a dotplot/heatmap of normalised enrichment (NES) across
# exposures x tissues, showing the significant (p.adjust < 0.05) enrichments.
#
#   x          = exposure (cID)
#   y          = GTEx tissue (ID = GSEA term key)
#   colour     = NES               (diverging blue-white-red)
#   size       = -log10(p.adjust)  (enrichment significance)
#
# This is a PLOTTING-ONLY script. The GSEA itself lives in
# scripts/module4_enrichment/ (03_run_gsea.R -> 04_enrichment_tables.R); here we
# only read the flattened CSV and draw it.
#
# Input  : heap_project_output("module4_enrichment","tissue_enrichment.csv")
#            columns: ID, setSize, NES, p.adjust, cID   (cID = exposure id)
#            (override the directory for testing via HEAP_ENRICH_DIR)
# Output : figures/main/fig_tissue_enrichment.{pdf,png} + figures/data/...tsv
#
# Run:
#   HEAP_PATHS_FILE=.../workflow/00_paths.R \
#     Rscript scripts/visualizations/figures/fig_tissue_enrichment.R
#
# Supersedes (legacy): scripts/visualizations/Visualizations/Module2/HEAPassoc_pathwayviz.R
#   (ComplexHeatmap of TE_df_dir NES; tissue half of the combined F2 heatmap).
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
a <- a[!a %in% c("fig_tissue_enrichment", "all_main", "all_supplement", "all", "website")]
figure_id <- Sys.getenv("HEAP_FIGURE_ID", unset = "fig_tissue_enrichment")
top_n <- suppressWarnings(as.integer(Sys.getenv("HEAP_ENRICH_TOPN", unset = "")))  # optional cap on tissues

# --- locate the flattened Module 4 tissue-enrichment table ------------------
# Canonical: heap_project_output("module4_enrichment","tissue_enrichment.csv").
# HEAP_ENRICH_DIR overrides the directory (used for synthetic-data testing).
enrich_dir <- Sys.getenv("HEAP_ENRICH_DIR", unset = "")
tissue_csv <- if (nzchar(enrich_dir)) {
  file.path(enrich_dir, "tissue_enrichment.csv")
} else {
  heap_project_output("module4_enrichment", "tissue_enrichment.csv")
}

if (!file.exists(tissue_csv))
  stop("[fig_tissue_enrichment] missing tissue enrichment table:\n  ", tissue_csv,
       "\nRun the Module 4 enrichment pipeline first:\n",
       "  scripts/module4_enrichment/03_run_gsea.R  (GSEA -> HEAPgsea.qs)\n",
       "  scripts/module4_enrichment/04_enrichment_tables.R  (flatten -> tissue_enrichment.csv)",
       call. = FALSE)

te <- fread(tissue_csv)

# verified schema from 04_enrichment_tables.R: ID, setSize, NES, p.adjust, cID
req <- c("ID", "NES", "p.adjust", "cID")
miss <- setdiff(req, names(te))
if (length(miss))
  stop("[fig_tissue_enrichment] ", tissue_csv, " missing column(s): ",
       paste(miss, collapse = ", "),
       "\nExpected the 04_enrichment_tables.R schema (ID, setSize, NES, p.adjust, cID).",
       call. = FALSE)

# --- prepare plotted data: significant tissue enrichments -------------------
te <- te[is.finite(NES) & is.finite(p.adjust)]
te[, sig := p.adjust < 0.05]
plot_dt <- te[sig == TRUE]
if (!nrow(plot_dt))
  stop("[fig_tissue_enrichment] no significant (p.adjust < 0.05) tissue enrichments in\n  ",
       tissue_csv, call. = FALSE)

# canonical short exposure labels; keep tissue ID verbatim (GTEx term key)
plot_dt[, exposure := heap_exposure_label(cID)]
plot_dt[, tissue   := as.character(ID)]
plot_dt[, neglog10p := -log10(pmax(p.adjust, .Machine$double.xmin))]

# optionally keep only the most-significant tissues (across all exposures)
if (!is.na(top_n) && top_n > 0L) {
  keep_tissue <- plot_dt[, .(score = max(neglog10p)), by = tissue][
    order(-score)][seq_len(min(top_n, .N)), tissue]
  plot_dt <- plot_dt[tissue %in% keep_tissue]
}

# map exposure -> category, GROUP exposures by category on y (and colour labels),
# tissues ordered by total signal along x. NES sign (direction) is meaningful per exposure.
.ae <- tryCatch(heap_exposure_table(), error = function(e) NULL)
plot_dt[, cat := if (!is.null(.ae)) as.character(setNames(.ae$category, .ae$variable)[cID]) else "Other"]
plot_dt[, cat := heap_category_factor(cat)]
# exposures on y: grouped by category, then by signal; reverse so the first
# category sits at the TOP of the panel (ggplot draws y bottom-up).
ex_ord <- plot_dt[, .(s = sum(neglog10p)), by = .(exposure, cat)][order(cat, -s)]
plot_dt[, exposure := factor(exposure, levels = rev(ex_ord$exposure))]
ti_order <- plot_dt[, .(s = sum(neglog10p)), by = tissue][order(-s), tissue]
# clean GTEx tissue display labels, preserving the signal-based ordering (x-axis)
plot_dt[, tissue   := factor(heap_pretty_tissue(tissue),
                             levels = unique(heap_pretty_tissue(ti_order)))]
# y-axis label colours, ordered to match the reversed exposure factor levels
ycols <- rev(HEAP_ECAT_COLORS[as.character(ex_ord$cat)]); ycols[is.na(ycols)] <- "grey30"
n_ex <- nlevels(plot_dt$exposure); n_ti <- nlevels(plot_dt$tissue)
nes_lim <- max(abs(plot_dt$NES), na.rm = TRUE)

# --- plot: diverging NES dotplot/heatmap (full supplement) ------------------
# EXPOSURES on the long (y) axis, one row each; GTEx tissues run along x. This is
# a COMPREHENSIVE REFERENCE figure (every significant enrichment kept), so the
# canvas is sized generously (height scales with the ~122 exposure rows, width
# with the ~51 tissue columns) so axis labels are legible when the PDF is zoomed.
fig_h <- max(9, 0.16 * n_ex + 1.6)   # ~0.16 in/exposure row
fig_w <- max(9, 0.16 * n_ti + 3.4)   # ~0.16 in/tissue col + legend/label gutter
p <- ggplot(plot_dt, aes(tissue, exposure)) +
  geom_point(aes(colour = NES, size = neglog10p)) +
  scale_colour_gradient2(low = HEAP_PAL_COMPONENT[["Genetic"]], mid = "white",
                         high = "#B2182B", midpoint = 0,
                         limits = c(-nes_lim, nes_lim), name = "NES") +
  scale_size_continuous(range = c(0.5, 3.0), name = expression(-log[10](p[adj]))) +
  labs(title = "Tissue-specificity enrichment of exposure-associated proteins",
       subtitle = NULL,
       x = "GTEx tissue", y = NULL) +
  theme_heap(base_size = 8) +
  theme(plot.title    = element_text(size = 11),
        plot.subtitle = element_text(size = 8),
        axis.title.x  = element_text(size = 9),
        axis.text.x   = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
        axis.text.y   = element_text(size = 6.5, colour = ycols),
        axis.ticks    = element_line(linewidth = 0.2),
        panel.grid.major = element_line(linewidth = 0.15),
        legend.position = "right", legend.key.size = unit(0.4, "cm"),
        legend.title  = element_text(size = 8), legend.text = element_text(size = 7),
        plot.margin   = margin(3, 3, 3, 3))

heap_emit_figure(p, figure_id, data = plot_dt, category = "supplement",
                 formats = c("pdf", "png"),
                 width = fig_w, height = fig_h,
                 website = FALSE)
message("fig_tissue_enrichment: done (", nrow(plot_dt), " significant tissue enrichments).")
