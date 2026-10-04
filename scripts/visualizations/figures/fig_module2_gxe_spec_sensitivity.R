#!/usr/bin/env Rscript
# ============================================================================
# fig_module2_gxe_spec_sensitivity.R  == Module-2 GxE covariate-spec sensitivity (supp) ==
# ----------------------------------------------------------------------------
# Robustness of the polygenic GxE interactions across covariate specifications.
# GxE F-tests are directionless (no per-pair beta), so this is p-value /
# replication based. Reads precomputed tables (analysis = scripts/analysis_
# summaries/module2_gxe_spec_sensitivity.R); NO analysis here.
#   a  # replicated GxE per spec, split cis / trans / joint-only
#   b  retention: % of base-replicated GxE kept per spec (+ Spearman of -log10 p)
#   c  base vs +clinical -log10 p_GxE_joint (base-replicated); FOLR3/CCL3 hubs
# ============================================================================
local({
  cand <- c(file.path(getwd(), "scripts", "visualizations", "common"),
            "/n/groups/patel/shakson_ukb/HEAP/scripts/visualizations/common")
  common <- cand[dir.exists(cand)][1]; if (is.na(common)) stop("cannot locate common/")
  for (f in c("figure_paths","load_heap_results","plot_theme","label_helpers","export_helpers"))
    source(file.path(common, paste0(f, ".R")))
})
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(ggrepel); library(patchwork) })

figure_id <- Sys.getenv("HEAP_FIGURE_ID", unset = "fig_module2_gxe_spec_sensitivity")
SS <- file.path(heap_project_output("module2"), "spec_sensitivity")
summ <- fread(file.path(SS, "gxe_spec_summary.tsv"))
att  <- fread(file.path(SS, "gxe_attenuation.tsv"))
GXE <- c(cis = "#1B6CA8", trans = "#D55E00", "joint-only" = "grey60")

specord <- summ[order(kind != "primary", -n_repl), lab]

## a -- # replicated GxE per spec, cis/trans/joint-only -----------------------
long <- melt(summ[, .(lab, n_cis, n_trans, n_jointonly)], id.vars = "lab",
             variable.name = "component", value.name = "n")
long[, component := factor(c(n_cis = "cis", n_trans = "trans", n_jointonly = "joint-only")[as.character(component)],
                           levels = c("cis","trans","joint-only"))]
long[, lab := factor(lab, levels = specord)]
pa <- ggplot(long, aes(n, lab, fill = component)) +
  geom_col(width = 0.72) +
  geom_text(data = summ[, .(lab = factor(lab, levels = specord), n_repl)], aes(n_repl, lab, label = n_repl),
            inherit.aes = FALSE, hjust = -0.25, size = 2.2, fontface = "bold", colour = "grey20") +
  scale_fill_manual(values = GXE, name = NULL) +
  coord_cartesian(xlim = c(0, max(summ$n_repl) * 1.12)) +
  labs(title = NULL, subtitle = NULL,
       x = "# replicated GxE (protein × exposure)", y = NULL) +
  theme_heap(base_size = 8) + theme(plot.subtitle = element_text(size = 6, colour = "grey35"),
        legend.position = c(0.98, 0.04), legend.justification = c(1, 0), legend.key.size = unit(0.3, "cm"), legend.text = element_text(size = 6))

## b -- retention of base-replicated GxE + concordance ------------------------
cb <- summ[exp != "M2_base_main"][order(pct_retained)]; cb[, lab := factor(lab, levels = lab)]
pb <- ggplot(cb, aes(pct_retained, lab)) +
  geom_col(width = 0.72, fill = "#4575B4") +
  geom_text(aes(label = sprintf("%.0f%%", pct_retained)), hjust = 1.1, size = 2.2, colour = "white", fontface = "bold") +
  geom_text(aes(x = 1, label = sprintf("rho=%.2f", spearman_repl)), hjust = 0, size = 2.0, colour = "grey92") +
  coord_cartesian(xlim = c(0, 100)) +
  labs(title = NULL, subtitle = NULL,
       x = "% of base-replicated GxE retained", y = NULL) +
  theme_heap(base_size = 8) + theme(plot.subtitle = element_text(size = 6, colour = "grey35"))

## c -- base vs +clinical -log10 p_GxE_joint ----------------------------------
cl <- att[spec == "M2_base_clinical_main"]
cl[, retained := ifelse(repl_s, "still replicated", "lost")]
HUB <- c("FOLR3","CCL3","CCL4","CXCL5","FOLR2","IL6","LEP","ADIPOQ","CRP","HGF","GDF15")
cl[, lab2 := ifelse(omicID %in% HUB | logp_base >= sort(logp_base, decreasing = TRUE)[min(.N, 12)], omicID, NA_character_)]
labset <- cl[!is.na(lab2)][order(-logp_base)][!duplicated(omicID)][1:14]
mx <- max(c(cl$logp_base, cl$logp_s), na.rm = TRUE)
pc <- ggplot(cl, aes(logp_base, logp_s)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey60", linewidth = 0.3) +
  geom_point(aes(colour = retained), size = 0.7, alpha = 0.5) +
  geom_point(data = labset, colour = "black", size = 1.0) +
  ggrepel::geom_text_repel(data = labset, aes(label = omicID), size = 2.0, max.overlaps = 20,
                           bg.color = "white", bg.r = 0.12, segment.size = 0.2, min.segment.length = 0) +
  scale_colour_manual(values = c("still replicated" = "grey55", "lost" = "#B2182B"), name = NULL) +
  coord_equal(xlim = c(0, mx), ylim = c(0, mx)) +
  labs(title = NULL, subtitle = NULL,
       x = expression(-log[10]~p[GxE]~"(base)"), y = expression(-log[10]~p[GxE]~"(+ clinical)")) +
  theme_heap(base_size = 8) + theme(plot.subtitle = element_text(size = 6, colour = "grey35"),
        legend.position = c(0.02, 0.98), legend.justification = c(0, 1), legend.key.size = unit(0.3, "cm"), legend.text = element_text(size = 6))

p <- (pa | pb) / pc + plot_layout(heights = c(1, 1.25)) +
  plot_annotation(tag_levels = "a") & theme(plot.tag = element_text(face = "bold", size = 11))

heap_emit_figure(p, figure_id, data = summ, category = "supplement", website = FALSE, width = 8.2, height = 7.2)
