#!/usr/bin/env Rscript
# ============================================================================
# fig_module2_spec_sensitivity.R   == Module-2 covariate-spec sensitivity (supp) ==
# ----------------------------------------------------------------------------
# Robustness of exposure->protein associations across covariate specifications.
# Reads precomputed tables (analysis = scripts/analysis_summaries/
# module2_spec_sensitivity.R); NO analysis here.
#   a  beta concordance with base (Spearman, base-replicated pairs) per spec
#   b  # replicated associations per spec
#   c  base vs +clinical beta scatter (base-replicated): BMI/clinical-mediated
#      attenuation, most-attenuated proteins labeled (LEP/FABP4 ...)
# ============================================================================
local({
  cand <- c(file.path(getwd(), "scripts", "visualizations", "common"),
            "/n/groups/patel/shakson_ukb/HEAP/scripts/visualizations/common")
  common <- cand[dir.exists(cand)][1]; if (is.na(common)) stop("cannot locate common/")
  for (f in c("figure_paths","load_heap_results","plot_theme","label_helpers","export_helpers"))
    source(file.path(common, paste0(f, ".R")))
})
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(ggrepel); library(patchwork) })

figure_id <- Sys.getenv("HEAP_FIGURE_ID", unset = "fig_module2_spec_sensitivity")
SS <- file.path(heap_project_output("module2"), "spec_sensitivity")
summ <- fread(file.path(SS, "spec_summary.tsv"))
att  <- fread(file.path(SS, "attenuation.tsv"))

KCOL <- c(primary = "#333333", covariate = "#2166AC", sample = "#1A9850",
          estimand = "#9E9E9E")   # grey: an alternative ESTIMAND, not a robustness check
summ[, lab := factor(lab, levels = lab[order(kind != "primary", -n_repl)])]

## a -- beta concordance with base (exclude base = trivial 1.0) --------------
ca <- summ[exp != "M2_base_main"][order(spearman_repl)]
ca[, lab := factor(lab, levels = lab)]
pa <- ggplot(ca, aes(spearman_repl, lab, fill = kind)) +
  geom_col(width = 0.72) +
  geom_text(aes(label = sprintf("%.3f", spearman_repl)), hjust = 1.1, size = 2.2, colour = "white", fontface = "bold") +
  geom_text(aes(x = 0.01, label = sprintf("%.0f%% kept", pct_retained)), hjust = 0, size = 2.0, colour = "grey92") +
  scale_fill_manual(values = KCOL, guide = "none") +
  coord_cartesian(xlim = c(0, 1)) +
  labs(title = NULL, subtitle = NULL, x = expression(Spearman~rho), y = NULL) +
  theme_heap(base_size = 8) + theme(plot.subtitle = element_text(size = 6, colour = "grey35"))

## b -- # replicated per spec ------------------------------------------------
cb <- copy(summ)[order(kind != "primary", -n_repl)]; cb[, lab := factor(lab, levels = lab)]
pb <- ggplot(cb, aes(n_repl, lab, fill = kind)) +
  geom_col(width = 0.72) +
  geom_text(aes(label = formatC(n_repl, big.mark = ",", format = "d")), hjust = 1.1, size = 2.2, colour = "white", fontface = "bold") +
  scale_fill_manual(values = KCOL, name = NULL,
                    labels = c(primary = "base", covariate = "+ covariate", sample = "sample",
                               estimand = "alt. estimand*")) +
  labs(title = NULL, subtitle = NULL, x = "# replicated (protein × exposure)", y = NULL) +
  theme_heap(base_size = 8) + theme(plot.subtitle = element_text(size = 6, colour = "grey35"),
        legend.position = c(0.99, 0.99), legend.justification = c(1, 1), legend.key.size = unit(0.3, "cm"), legend.text = element_text(size = 6))

## c -- base vs +clinical attenuation ---------------------------------------
cl <- att[spec == "M2_base_clinical_main"]
cl[, retained := ifelse(repl_s, "still replicated", "lost")]
ALLOW <- c("LEP","FABP4","RETN","LPL","IGFBP1","IGFBP2","GDF15","MMP12","ADIPOQ","CRP","HGF","WFDC2","CEACAM5","APOM","SELENOP","IL6","PLAUR","CXCL17")
cl[, drop := abs(beta_base) - abs(beta_s)]
labset <- cl[omicID %in% ALLOW][order(-drop)][!duplicated(omicID)][1:12]
rng <- range(c(cl$beta_base, cl$beta_s), na.rm = TRUE)
pc <- ggplot(cl, aes(beta_base, beta_s)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey60", linewidth = 0.3) +
  geom_point(aes(colour = retained), size = 0.5, alpha = 0.35) +
  geom_point(data = labset, colour = "black", size = 0.9) +
  ggrepel::geom_text_repel(data = labset, aes(label = omicID), size = 2.0, max.overlaps = 20,
                           bg.color = "white", bg.r = 0.12, segment.size = 0.2, min.segment.length = 0) +
  scale_colour_manual(values = c("still replicated" = "grey55", "lost" = "#B2182B"), name = NULL) +
  coord_equal(xlim = rng, ylim = rng) +
  labs(title = NULL, subtitle = NULL,
       x = expression(beta[base]), y = expression(beta["+ clinical"])) +
  theme_heap(base_size = 8) + theme(plot.subtitle = element_text(size = 6, colour = "grey35"),
        legend.position = c(0.02, 0.98), legend.justification = c(0, 1), legend.key.size = unit(0.3, "cm"), legend.text = element_text(size = 6))

p <- (pa | pb) / pc + plot_layout(heights = c(1, 1.25)) +
  plot_annotation(tag_levels = "a") & theme(plot.tag = element_text(face = "bold", size = 11))

heap_emit_figure(p, figure_id, data = summ, category = "supplement", website = FALSE, width = 8.2, height = 7.2)
