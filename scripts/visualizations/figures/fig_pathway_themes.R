#!/usr/bin/env Rscript
# ============================================================================
# fig_pathway_themes.R  [figure_id: fig_pathway_themes]
# ----------------------------------------------------------------------------
# Exposure-associated proteins enrich a BROAD space of pathways, organised into
# biological THEMES. Dotplot of every significant category-level Reactome
# enrichment (NES>0, FDR<0.05), x = exposure category, y = pathway grouped into
# themes (Immune/Inflammation, Growth-factor/Hormone, ECM/Structural,
# Lipid/Metabolism, Protein processing, Infection). Shows breadth + structure —
# inflammation is one theme among several, not the whole story.
#
# Input : module4_enrichment/pathway_enrichment_category.csv
# Run   : HEAP_PATHS_FILE=.../00_paths.R Rscript .../fig_pathway_themes.R
# ============================================================================

local({
  cand <- c(file.path(getwd(), "scripts", "visualizations", "common"),
            "/n/groups/patel/shakson_ukb/HEAP/scripts/visualizations/common")
  common <- cand[dir.exists(cand)][1]; if (is.na(common)) stop("cannot locate common/")
  for (f in c("figure_paths","plot_theme","label_helpers","export_helpers")) source(file.path(common, paste0(f, ".R")))
})
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
figure_id <- Sys.getenv("HEAP_FIGURE_ID", unset = "fig_pathway_themes")

# --- pathway -> theme map (curated over the enriched Reactome set) ----------
THEME <- list(
  "Immune & inflammation" = c("Neutrophil degranulation","Immunoregulatory interactions between a Lymphoid and a non-Lymphoid cell",
    "Interleukin-10 signaling","TNFs bind their physiological receptors","Innate Immune System",
    "TNFR2 non-canonical NF-kB pathway","Adaptive Immune System","Complement cascade",
    "Initial triggering of complement","Chemokine receptors bind chemokines"),
  "Growth factor & hormone" = c("Regulation of Insulin-like Growth Factor (IGF) transport and uptake by Insulin-like Growth Factor Binding Proteins (IGFBPs)",
    "Signaling by TGFB family members","Peptide ligand-binding receptors","Peptide hormone metabolism",
    "EPH-ephrin mediated repulsion of cells"),
  "Extracellular matrix" = c("Extracellular matrix organization","Integrin cell surface interactions","Collagen degradation",
    "Degradation of the extracellular matrix","Elastic fibre formation","Activation of Matrix Metalloproteinases",
    "Glycosaminoglycan metabolism","Keratinization","Formation of the cornified envelope"),
  "Lipid & metabolism" = c("Plasma lipoprotein remodeling","Plasma lipoprotein assembly, remodeling, and clearance",
    "Adipogenesis","Transcriptional regulation of white adipocyte differentiation","Diseases of metabolism","Diseases of glycosylation"),
  "Protein processing" = c("Post-translational protein phosphorylation","Drug ADME"),
  "Infection" = c("Attachment and Entry","Early SARS-CoV-2 Infection Events"))
theme_dt <- rbindlist(lapply(names(THEME), function(t) data.table(Description = THEME[[t]], theme = t)))
THEME_LEVELS <- names(THEME)
THEME_COLS <- setNames(c("#CC3344","#117733","#6699CC","#E69F00","#999999","#AA4499"), THEME_LEVELS)

# short display labels for the long Reactome names
RELAB <- c(
  "Regulation of Insulin-like Growth Factor (IGF) transport and uptake by Insulin-like Growth Factor Binding Proteins (IGFBPs)" = "IGF transport & uptake (IGFBPs)",
  "Immunoregulatory interactions between a Lymphoid and a non-Lymphoid cell" = "Immunoregulatory lymphoid interactions",
  "Plasma lipoprotein assembly, remodeling, and clearance" = "Plasma lipoprotein assembly/clearance",
  "Transcriptional regulation of white adipocyte differentiation" = "White adipocyte differentiation",
  "TNFs bind their physiological receptors" = "TNF receptor binding",
  "TNFR2 non-canonical NF-kB pathway" = "TNFR2 / NF-kB signaling",
  "Early SARS-CoV-2 Infection Events" = "SARS-CoV-2 early infection",
  "Activation of Matrix Metalloproteinases" = "MMP activation")
relab <- function(x) heap_americanize(ifelse(x %in% names(RELAB), RELAB[x], x))

ED <- heap_project_output("module4_enrichment")
pa <- fread(file.path(ED,"pathway_enrichment_category.csv"))[is.finite(NES) & NES > 0 & p.adjust < 0.05]
pa[, neglogp := -log10(pmax(p.adjust, .Machine$double.xmin))]
pa[, cat := heap_category_factor(cID)]
pa <- merge(pa, theme_dt, by = "Description", all.x = TRUE)
unmapped <- unique(pa[is.na(theme), Description]); if (length(unmapped)) message("UNMAPPED: ", paste(unmapped, collapse=" | "))
pa <- pa[!is.na(theme)]
pa[, theme := factor(theme, levels = THEME_LEVELS)]
pa[, plab := relab(Description)]
# order pathways within theme by total breadth; order exposure categories by total signal
ford <- pa[, .(s = sum(neglogp), nc = uniqueN(cID)), by = .(theme, plab)][order(theme, nc, s)]
pa[, plab := factor(plab, levels = ford$plab)]
catord <- pa[, .(s = .N), by = cat][order(s), cat]; pa[, cat := factor(as.character(cat), levels = as.character(catord))]
bg <- data.table(theme = factor(THEME_LEVELS, levels = THEME_LEVELS))   # one row/theme -> colours each facet panel
message(sprintf("pathway themes: %d pathways in %d themes x %d exposure categories", uniqueN(pa$plab), uniqueN(pa$theme), uniqueN(pa$cat)))

p <- ggplot(pa, aes(cat, plab)) +
  geom_rect(data = bg, aes(fill = theme), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,
            alpha = 0.13, inherit.aes = FALSE) +
  scale_fill_manual(values = THEME_COLS, guide = "none") +
  geom_point(aes(colour = NES, size = neglogp)) +
  facet_grid(theme ~ ., scales = "free_y", space = "free_y", switch = "y",
             labeller = labeller(theme = label_wrap_gen(14))) +
  scale_colour_gradient(low = "#FEE0D2", high = "#A50F15", name = "NES") +
  scale_size_continuous(range = c(1.4, 5), name = expression(-log[10](p[adj]))) +
  scale_x_discrete(labels = heap_category_pretty) +
  labs(title = "Pathway themes of exposure-associated proteins",
       subtitle = NULL,
       x = NULL, y = NULL) +
  theme_heap() +
  theme(axis.text.x = element_text(angle = 40, hjust = 1, size = 8),
        axis.text.y = element_text(size = 6.8),
        strip.text.y.left = element_text(angle = 0, face = "bold", size = 7.5),
        strip.placement = "outside", panel.spacing.y = unit(2.5, "pt"),
        legend.position = "right", plot.subtitle = element_text(colour = "grey35", size = 9))

heap_emit_figure(p, figure_id, data = pa[, .(category = cID, pathway = Description, theme, NES, p.adjust)],
                 category = "supplement", formats = c("pdf","png"), width = 9, height = 7, website = TRUE)
message("fig_pathway_themes: done.")
