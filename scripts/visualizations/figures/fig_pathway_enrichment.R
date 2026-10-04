#!/usr/bin/env Rscript

# ============================================================================
# fig_pathway_enrichment.R  [figure_id: fig_pathway_enrichment]
# ----------------------------------------------------------------------------
# Pathway (Reactome) GSEA of Module 2 exposure->protein associations: for each
# significant exposure, its protein t-value ranking was tested against Reactome
# pathways (ReactomePA::gsePathway). This figure summarises the flattened result
# table as a dotplot/heatmap of normalised enrichment (NES) across exposures x
# the top-N most-significant pathways.
#
#   x          = exposure (cID)
#   y          = Reactome pathway (Description, wrapped for long names)
#   colour     = NES               (diverging blue-white-red)
#   size       = -log10(p.adjust)  (enrichment significance)
#
# This is a PLOTTING-ONLY script. The GSEA itself lives in
# scripts/module4_enrichment/ (03_run_gsea.R -> 04_enrichment_tables.R); here we
# only read the flattened CSV and draw it.
#
# Input  : heap_project_output("module4_enrichment","pathway_enrichment.csv")
#            columns: Description, setSize, NES, p.adjust, cID  (cID = exposure id)
#            (override the directory for testing via HEAP_ENRICH_DIR)
# Output : figures/main/fig_pathway_enrichment.{pdf,png} + figures/data/...tsv
#
# Run:
#   HEAP_PATHS_FILE=.../workflow/00_paths.R \
#     Rscript scripts/visualizations/figures/fig_pathway_enrichment.R [top_n]
#
# Supersedes (legacy): scripts/visualizations/Visualizations/Module2/HEAPassoc_pathwayviz.R
#   (ComplexHeatmap of Path_dir NES; pathway half of the combined F2 heatmap).
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
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(stringr) })

a <- commandArgs(trailingOnly = TRUE); a <- a[!startsWith(a, "--")]
a <- a[!a %in% c("fig_pathway_enrichment", "all_main", "all_supplement", "all", "website")]
figure_id <- Sys.getenv("HEAP_FIGURE_ID", unset = "fig_pathway_enrichment")
# top-N pathways by significance (1st positional arg, or HEAP_ENRICH_TOPN, default 25)
top_n <- {
  v <- if (length(a) >= 1) a[1] else Sys.getenv("HEAP_ENRICH_TOPN", unset = "25")
  v <- suppressWarnings(as.integer(v)); if (is.na(v) || v <= 0L) 25L else v
}

# --- locate the flattened Module 4 pathway-enrichment table -----------------
# Canonical: heap_project_output("module4_enrichment","pathway_enrichment.csv").
# HEAP_ENRICH_DIR overrides the directory (used for synthetic-data testing).
enrich_dir <- Sys.getenv("HEAP_ENRICH_DIR", unset = "")
pathway_csv <- if (nzchar(enrich_dir)) {
  file.path(enrich_dir, "pathway_enrichment.csv")
} else {
  heap_project_output("module4_enrichment", "pathway_enrichment.csv")
}

if (!file.exists(pathway_csv))
  stop("[fig_pathway_enrichment] missing pathway enrichment table:\n  ", pathway_csv,
       "\nRun the Module 4 enrichment pipeline first:\n",
       "  scripts/module4_enrichment/03_run_gsea.R  (GSEA -> HEAPgsea.qs)\n",
       "  scripts/module4_enrichment/04_enrichment_tables.R  (flatten -> pathway_enrichment.csv)",
       call. = FALSE)

pe <- fread(pathway_csv)

# verified schema from 04_enrichment_tables.R: Description, setSize, NES, p.adjust, cID
req <- c("Description", "NES", "p.adjust", "cID")
miss <- setdiff(req, names(pe))
if (length(miss))
  stop("[fig_pathway_enrichment] ", pathway_csv, " missing column(s): ",
       paste(miss, collapse = ", "),
       "\nExpected the 04_enrichment_tables.R schema (Description, setSize, NES, p.adjust, cID).",
       call. = FALSE)

# --- prepare plotted data: top-N significant pathway enrichments ------------
pe <- pe[is.finite(NES) & is.finite(p.adjust)]
pe[, sig := p.adjust < 0.05]
sig_dt <- pe[sig == TRUE]
if (!nrow(sig_dt))
  stop("[fig_pathway_enrichment] no significant (p.adjust < 0.05) pathway enrichments in\n  ",
       pathway_csv, call. = FALSE)

sig_dt[, exposure  := heap_exposure_label(cID)]   # canonical short labels
sig_dt[, neglog10p := -log10(pmax(p.adjust, .Machine$double.xmin))]

# keep only the top-N pathways by peak significance (across all exposures)
keep_path <- sig_dt[, .(score = max(neglog10p)), by = Description][
  order(-score)][seq_len(min(top_n, .N)), Description]
plot_dt <- sig_dt[Description %in% keep_path]

# truncate over-long Reactome pathway names so the rotated x-axis ticks stay on a
# single readable line (full names are preserved in the figure data TSV)
plot_dt[, pathway := str_trunc(as.character(Description), 58)]

# map exposure -> category, GROUP exposures by category on y (colour labels),
# pathways ordered by total signal along x. NES sign is meaningful per individual exposure.
.ae <- tryCatch(heap_exposure_table(), error = function(e) NULL)
plot_dt[, cat := if (!is.null(.ae)) as.character(setNames(.ae$category, .ae$variable)[cID]) else "Other"]
plot_dt[, cat := heap_category_factor(cat)]
# exposures on y: grouped by category, then signal; reverse so first category is at TOP
ex_ord <- plot_dt[, .(s = sum(neglog10p)), by = .(exposure, cat)][order(cat, -s)]
plot_dt[, exposure := factor(exposure, levels = rev(ex_ord$exposure))]
# order by total signal; keep label levels unique even if two names share a
# truncated 45-char prefix (rare) by ordering on the unique Description.
pa_order <- plot_dt[, .(s = sum(neglog10p)), by = .(Description, pathway)][order(-s)]
plot_dt[, pathway  := factor(pathway, levels = unique(pa_order$pathway))]
ycols <- rev(HEAP_ECAT_COLORS[as.character(ex_ord$cat)]); ycols[is.na(ycols)] <- "grey30"
n_ex <- nlevels(plot_dt$exposure); n_pa <- nlevels(plot_dt$pathway)
nes_lim <- max(abs(plot_dt$NES), na.rm = TRUE)

# --- plot: diverging NES dotplot/heatmap (full supplement) ------------------
# EXPOSURES on the long (y) axis, one row each; Reactome pathways run along x.
# COMPREHENSIVE REFERENCE figure: height scales with the ~105 exposure rows and
# width with the top-N pathway columns so the y labels (exposures) and rotated
# x labels (truncated pathway names) are legible when the PDF is zoomed.
fig_h <- max(9, 0.16 * n_ex + 1.6)   # ~0.16 in/exposure row
fig_w <- max(8, 0.30 * n_pa + 4.0)   # ~0.30 in/pathway col + label/legend gutter
p <- ggplot(plot_dt, aes(pathway, exposure)) +
  geom_point(aes(colour = NES, size = neglog10p)) +
  scale_colour_gradient2(low = HEAP_PAL_COMPONENT[["Genetic"]], mid = "white",
                         high = "#B2182B", midpoint = 0,
                         limits = c(-nes_lim, nes_lim), name = "NES") +
  scale_size_continuous(range = c(0.5, 3.0), name = expression(-log[10](p[adj]))) +
  labs(title = "Reactome pathway enrichment of exposure-associated proteins",
       subtitle = NULL,
       x = "Reactome pathway", y = NULL) +
  theme_heap(base_size = 8) +
  theme(plot.title    = element_text(size = 11),
        plot.subtitle = element_text(size = 8),
        axis.title.x  = element_text(size = 9),
        axis.text.x   = element_text(angle = 45, hjust = 1, vjust = 1, size = 7),
        axis.text.y   = element_text(size = 6.5, colour = ycols),
        axis.ticks    = element_line(linewidth = 0.2),
        panel.grid.major = element_line(linewidth = 0.15),
        legend.position = "right", legend.key.size = unit(0.4, "cm"),
        legend.title  = element_text(size = 8), legend.text = element_text(size = 7),
        plot.margin   = margin(3, 3, 3, 3))

heap_emit_figure(p, figure_id, data = plot_dt, category = "supplement",
                 formats = c("pdf", "png"), website = FALSE,
                 width = fig_w, height = fig_h)
message("fig_pathway_enrichment: done (top ", length(unique(plot_dt$pathway)),
        " pathways, ", nrow(plot_dt), " significant enrichments).")
