#!/usr/bin/env Rscript
# ============================================================================
# fig_tissue_themes.R  [figure_id: fig_tissue_themes]
# ----------------------------------------------------------------------------
# Companion to fig_pathway_themes: exposure-associated proteins are enriched in a
# BROAD set of GTEx tissues, grouped into ORGAN SYSTEMS. Dotplot of every
# significant category-level tissue enrichment (NES>0, FDR<0.05), x = exposure
# category, y = tissue grouped by organ system (colored panel bands). Shows the
# breadth of exposure→tissue connections with anatomical structure.
#
# Input : module4_enrichment/tissue_enrichment_category.csv
# Run   : HEAP_PATHS_FILE=.../00_paths.R Rscript .../fig_tissue_themes.R
# ============================================================================

local({
  cand <- c(file.path(getwd(), "scripts", "visualizations", "common"),
            "/n/groups/patel/shakson_ukb/HEAP/scripts/visualizations/common")
  common <- cand[dir.exists(cand)][1]; if (is.na(common)) stop("cannot locate common/")
  for (f in c("figure_paths","plot_theme","label_helpers","export_helpers")) source(file.path(common, paste0(f, ".R")))
})
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
figure_id <- Sys.getenv("HEAP_FIGURE_ID", unset = "fig_tissue_themes")

# --- tissue -> organ-system map (over the enriched GTEx set) ----------------
ORGAN <- list(
  "Respiratory" = c("lung"),
  "Hepatic & digestive" = c("liver","stomach","small_intestine_terminal_ileum","colon_transverse",
                            "minor_salivary_gland","esophagus_mucosa"),
  "Adipose" = c("adipose_subcutaneous","adipose_visceral_omentum"),
  "Renal" = c("kidney_cortex","kidney_medulla"),
  "Vascular" = c("artery_coronary","artery_aorta","artery_tibial"),
  "Immune / blood" = c("spleen","whole_blood"),
  "Endocrine" = c("thyroid"),
  "Reproductive" = c("breast_mammary_tissue","vagina","cervix_endocervix"),
  "Neural" = c("brain_hippocampus","brain_putamen_basal_ganglia","brain_amygdala","nerve_tibial"),
  "Skin & connective" = c("skin_not_sun_exposed_suprapubic","skin_sun_exposed_lower_leg","cells_cultured_fibroblasts"))
organ_dt <- rbindlist(lapply(names(ORGAN), function(o) data.table(ID = ORGAN[[o]], organ = o)))
ORGAN_LEVELS <- names(ORGAN)
ORGAN_COLS <- setNames(c("#88CCEE","#117733","#E69F00","#44AA99","#882255","#CC3344",
                         "#DDCC77","#AA4499","#332288","#999933"), ORGAN_LEVELS)

ED <- heap_project_output("module4_enrichment")
ti <- fread(file.path(ED,"tissue_enrichment_category.csv"))[is.finite(NES) & NES > 0 & p.adjust < 0.05]
ti[, neglogp := -log10(pmax(p.adjust, .Machine$double.xmin))]
ti[, cat := heap_category_factor(cID)]
ti <- merge(ti, organ_dt, by = "ID", all.x = TRUE)
unmapped <- unique(ti[is.na(organ), ID]); if (length(unmapped)) message("UNMAPPED: ", paste(unmapped, collapse=" | "))
ti <- ti[!is.na(organ)]
ti[, organ := factor(organ, levels = ORGAN_LEVELS)]
ti[, tlab := heap_pretty_tissue(ID)]
ford <- ti[, .(s = sum(neglogp), nc = uniqueN(cID)), by = .(organ, tlab)][order(organ, nc, s)]
ti[, tlab := factor(tlab, levels = ford$tlab)]
catord <- ti[, .(s = .N), by = cat][order(s), cat]; ti[, cat := factor(as.character(cat), levels = as.character(catord))]
bg <- data.table(organ = factor(ORGAN_LEVELS, levels = ORGAN_LEVELS))
message(sprintf("tissue themes: %d tissues in %d organ systems x %d exposure categories", uniqueN(ti$tlab), uniqueN(ti$organ), uniqueN(ti$cat)))

p <- ggplot(ti, aes(cat, tlab)) +
  geom_rect(data = bg, aes(fill = organ), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,
            alpha = 0.13, inherit.aes = FALSE) +
  scale_fill_manual(values = ORGAN_COLS, guide = "none") +
  geom_point(aes(colour = NES, size = neglogp)) +
  facet_grid(organ ~ ., scales = "free_y", space = "free_y", switch = "y",
             labeller = labeller(organ = label_wrap_gen(12))) +
  scale_colour_gradient(low = "#FEE0D2", high = "#A50F15", name = "NES") +
  scale_size_continuous(range = c(1.4, 5), name = expression(-log[10](p[adj]))) +
  scale_x_discrete(labels = heap_category_pretty) +
  labs(title = "Exposure-associated proteins span a broad set of tissues (organ systems)",
       subtitle = NULL,
       x = NULL, y = NULL) +
  theme_heap() +
  theme(axis.text.x = element_text(angle = 40, hjust = 1, size = 8),
        axis.text.y = element_text(size = 6.8),
        strip.text.y.left = element_text(angle = 0, face = "bold", size = 7.5),
        strip.placement = "outside", panel.spacing.y = unit(2.5, "pt"),
        legend.position = "right", plot.subtitle = element_text(colour = "grey35", size = 9))

heap_emit_figure(p, figure_id, data = ti[, .(category = cID, tissue = ID, organ, NES, p.adjust)],
                 category = "supplement", formats = c("pdf","png"), width = 9, height = 6.5, website = TRUE)
message("fig_tissue_themes: done.")
