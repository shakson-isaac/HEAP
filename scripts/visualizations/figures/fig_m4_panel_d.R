#!/usr/bin/env Rscript
# ============================================================================
# fig_m4_panel_d.R  -- Fig6 panel d: MR-causal proteins (tiered + coloc)
# ----------------------------------------------------------------------------
# Two contrasting exposures (Strenuous sports; Processed meat) vs GLP1 STEP2 ->
# T2D, faceted. x = HEAP exposure->protein beta, y = trial protein shift. Points
# colored by the protein->T2D MR CONFIDENCE TIER (new Module-5 tiered+coloc
# output, mr_pd_tiered.tsv): Tier1 (robust) / Tier2 / Suggestive / not-causal;
# a ring = cis-pQTL colocalized (PP.H4>=0.8, gold). Tier2+ proteins labeled.
# ICAM1 (Tier1 cis-coloc) is the gold T2D edge; FABP4 is only Suggestive.
# Drug annotation removed (see support/druggability/drug_target_candidates.tsv).
#
# Inputs: support/intervention_compare/{intervention_scatter_mr,mr_pd_tiered,intervention_correlations}.tsv
# Run   : HEAP_PATHS_FILE=.../00_paths.R [HEAP_CELL=1] Rscript .../fig_m4_panel_d.R
# ============================================================================
local({
  cand <- c(file.path(getwd(), "scripts", "visualizations", "common"),
            "/n/groups/patel/shakson_ukb/HEAP/scripts/visualizations/common")
  cm <- cand[dir.exists(cand)][1]
  for (f in c("figure_paths.R","load_heap_results.R","plot_theme.R",
              "label_helpers.R","export_helpers.R")) source(file.path(cm, f))
})
suppressPackageStartupMessages({
  library(data.table); library(ggplot2)
  has_repel <- requireNamespace("ggrepel", quietly = TRUE)
})

CELL <- nzchar(Sys.getenv("HEAP_CELL"))
BS   <- if (CELL) 7 else 10
TTL  <- if (CELL) 9 else 13
LBL  <- if (CELL) 1.9 else 2.7
B_W  <- 2.85; B_H <- 3.75

INT  <- Sys.getenv("HEAP_PANELD_INTERVENTION", unset = "GLP1_1")   # STEP1 covers all 6 Tier1 causal proteins (STEP2 only 3)
DZ   <- Sys.getenv("HEAP_PANELD_DISEASE", unset = "finngen_R12_T2D")
EXPS <- data.table(
  exposure_id = c(Sys.getenv("HEAP_PANELD_EXP1", "types_of_physical_activity_in_last_4_weeks_f6164_0_0.multi_Strenuous_sports"),
                  Sys.getenv("HEAP_PANELD_EXP2", "processed_meat_intake_f1349_0_04")),
  elab        = c(Sys.getenv("HEAP_PANELD_LAB1", "Strenuous sports"),
                  Sys.getenv("HEAP_PANELD_LAB2", "Processed meat")))

EFF  <- c(HERITAGE = "HERITAGE_effect", GLP1_1 = "GLP1_effect1", GLP1_2 = "GLP1_effect2")
ILAB <- c(HERITAGE = "HERITAGE (exercise)", GLP1_1 = "GLP1 STEP1", GLP1_2 = "GLP1 STEP2")
# disease palette for the causal intermediates (the SAME 4-disease hub as panel d)
DZ_COL   <- c(`Type-2 diabetes` = "#D55E00", Obesity = "#0072B2", `Lipid disorder` = "#009E73", Hypertension = "#CC79A7")
DZ_SHORT <- c(`Type-2 diabetes` = "T2D", Obesity = "Obesity", `Lipid disorder` = "Lipids", Hypertension = "HTN")

FIGDIR <- file.path(Sys.getenv("HEAP_PROJECT_ROOT", "/n/groups/patel/IGLOO/UKB/HEAP"), "figures", "exploratory/module4")
dir.create(FIGDIR, recursive = TRUE, showWarnings = FALSE)
.read <- function(f) fread(file.path(heap_project_output("support","intervention_compare"), f))
sc <- .read("intervention_scatter_mr.tsv"); cor <- .read("intervention_correlations.tsv")
# causal intermediates = the SAME cast as panel d (shared-language network), with
# each protein's primary causal disease (colocalized edge first, else strongest).
ne <- .read("shared_language_network_edges.tsv")[etype == "gen_fwd"]
ne[, coloc := tier == "colocalized"]
setorder(ne, from, -coloc, -weight)
causal <- ne[, .SD[1], by = from][, .(protein = from, dz_causal = to, coloc_confirmed = coloc)]
eff_col <- EFF[[INT]]
keep <- intersect(c("exposure_id","protein","beta_HEAP","olink_soma_r","GLP1_effect2","GLP1_effect1","HERITAGE_effect"), names(sc))

build_d <- function(i) {
  dd <- sc[exposure_id == EXPS$exposure_id[i], ..keep][is.finite(beta_HEAP) & is.finite(get(eff_col))]
  if (!nrow(dd)) return(NULL)
  dd <- merge(dd, causal, by = "protein", all.x = TRUE)
  dd[is.na(coloc_confirmed), coloc_confirmed := FALSE]
  dd[, `:=`(effect = get(eff_col), elab = EXPS$elab[i])]
  dd
}
d <- rbindlist(lapply(seq_len(nrow(EXPS)), build_d), fill = TRUE)
d[, causal := !is.na(dz_causal)]
d[, dz_causal := factor(dz_causal, levels = names(DZ_COL))]
d[, elab := factor(elab, levels = EXPS$elab)]
d <- d[order(causal)]                                      # causal on top

ann <- merge(EXPS, cor[intervention == eff_col, .(exposure_id, r = as.numeric(r), p = as.numeric(pval_BH))], by = "exposure_id")
ann[, elab := factor(elab, levels = EXPS$elab)]
ann[, lab := sprintf("r = %s\nadj. p = %s", ifelse(is.finite(r), formatC(r,2,format="f"), "NA"),
                     ifelse(is.finite(p), formatC(p,2,format="g"), "NA"))]
message(sprintf("panel c (causal): %d causal proteins across %d diseases (both exposures); %d colocalized",
                d[causal == TRUE, uniqueN(protein)], d[causal == TRUE, uniqueN(dz_causal)],
                d[coloc_confirmed == TRUE, uniqueN(protein)]))

p <- ggplot(d, aes(beta_HEAP, effect)) +
  geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey80") +
  geom_vline(xintercept = 0, linewidth = 0.25, colour = "grey80") +
  geom_point(aes(colour = dz_causal, alpha = causal), size = if (CELL) 1.7 else 2.4) +
  geom_point(data = d[coloc_confirmed == TRUE], shape = 21, fill = NA, colour = "grey10",
             size = if (CELL) 3.2 else 4.6, stroke = if (CELL) 0.6 else 0.85) +
  scale_alpha_manual(values = c(`TRUE` = 0.95, `FALSE` = 0.28), guide = "none") +
  scale_colour_manual(values = DZ_COL, breaks = names(DZ_COL), labels = DZ_SHORT,
                      name = if (CELL) "causal for" else "causal for disease", na.value = "grey80", drop = FALSE) +
  geom_text(data = ann, aes(x = -Inf, y = Inf, label = lab), inherit.aes = FALSE,
            hjust = -0.12, vjust = 1.15, size = if (CELL) 1.95 else 2.6, lineheight = 0.9, colour = "grey20") +
  facet_wrap(~ elab, ncol = 1, scales = "free") +
  labs(title = if (CELL) NULL else paste0("Causal cardiometabolic proteins vs ", ILAB[[INT]]),
       # the terse "HEAP exposure->protein beta" / "<trial> shift" did not say what
       # either axis measures or in which direction
       x = "Observational effect\n(exposure → protein, β per SD)",
       y = paste0("Treatment effect\n(protein change in ", ILAB[[INT]], ", β)"),
       caption = NULL) +
  theme_heap(base_size = BS) +
  theme(legend.position = if (CELL) "bottom" else "right",
        legend.key.height = grid::unit(if (CELL) 0.7 else 0.9, "lines"),
        legend.text = element_text(size = if (CELL) 6.0 else 7), legend.title = element_text(size = if (CELL) 6.6 else 8),
        axis.title = element_text(face = if (CELL) "plain" else "bold"),
        plot.title = element_text(size = if (CELL) TTL else rel(1.05), hjust = 0.5, face = "bold"),
        plot.title.position = "panel",
        plot.caption = element_text(size = if (CELL) 5 else 7, colour = "grey35", hjust = 0),
        strip.text = element_text(size = if (CELL) 7.8 else 9, face = "bold"),
        panel.spacing = grid::unit(if (CELL) 0.3 else 0.5, "lines"),
        plot.margin = margin(3, 4, 1, 3)) +
  guides(colour = guide_legend(nrow = if (CELL) 2 else 1,
                               override.aes = list(size = if (CELL) 2 else 3, alpha = 1)))

if (has_repel) {
  lab <- d[causal == TRUE][, .SD[order(-coloc_confirmed)][seq_len(min(if (CELL) 12 else 12, .N))], by = elab]
  if (nrow(lab)) p <- p + ggrepel::geom_text_repel(data = lab, aes(label = protein, fontface = ifelse(coloc_confirmed, "bold", "plain")),
            size = LBL, max.overlaps = 30, segment.colour = "grey70", min.segment.length = 0,
            box.padding = 0.18, colour = "#222222", bg.color = "white", bg.r = 0.12)
}

if (CELL) {
  ggsave(file.path(FIGDIR, "fig_m4_panel_d_cell.png"), p, width = B_W, height = B_H, dpi = 400, bg = "white", type = "cairo")
  ggsave(file.path(FIGDIR, "fig_m4_panel_d_cell.pdf"), p, width = B_W, height = B_H, bg = "white", device = cairo_pdf)
  message("panel d CELL done (tiered)")
} else {
  ggsave(file.path(FIGDIR, "fig_m4_panel_d.png"), p, width = 4.6, height = 6.2, dpi = 200, bg = "white", type = "cairo")
  message("panel d standalone done (tiered)")
}
