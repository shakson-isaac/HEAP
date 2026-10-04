#!/usr/bin/env Rscript
# ============================================================================
# fig_module3_spec_sensitivity.R   == Module-3 mediation spec sensitivity (supp) ==
# ----------------------------------------------------------------------------
# Robustness of the mediation (natural indirect effect, NIE) estimates across
# covariate specifications, model families and sample definitions. Reads
# precomputed tables (analysis = scripts/analysis_summaries/
# module3_spec_sensitivity.R); NO analysis here.
#   a  NIE concordance with base (Spearman, base-sig exposomic pairs) per spec
#   b  # FDR-significant exposomic mediations per spec
#   c  base vs +clinical NIE scatter (base-sig): BMI/clinical-mediated attenuation,
#      most-attenuated proteins labeled (LEP/FABP4 ...)
#   d  per-exposure-category |NIE| attenuation under +clinical (partitioned):
#      which exposure axes are BMI/clinical-mediated (cardiometabolic attenuates)
# ============================================================================
local({
  cand <- c(file.path(getwd(), "scripts", "visualizations", "common"),
            "/n/groups/patel/shakson_ukb/HEAP/scripts/visualizations/common")
  common <- cand[dir.exists(cand)][1]; if (is.na(common)) stop("cannot locate common/")
  for (f in c("figure_paths","load_heap_results","plot_theme","label_helpers","export_helpers"))
    source(file.path(common, paste0(f, ".R")))
})
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(ggrepel); library(patchwork) })

figure_id <- Sys.getenv("HEAP_FIGURE_ID", unset = "fig_module3_spec_sensitivity")
SS   <- file.path(heap_project_output("module3"), "spec_sensitivity")
summ <- fread(file.path(SS, "spec_summary.tsv"))
att  <- fread(file.path(SS, "attenuation.tsv"))
catt <- if (file.exists(file.path(SS, "category_attenuation.tsv"))) fread(file.path(SS, "category_attenuation.tsv")) else NULL

KCOL <- c(primary = "#333333", covariate = "#2166AC", sample = "#1A9850", model = "#B2182B")
ord  <- function(d) d[order(kind != "primary", -n_sig)]

## a -- NIE concordance with base (exclude base = trivial 1.0) ---------------
ca <- summ[exp != "M3_base_lasso_primary"][order(spearman_sig)]
ca[, lab := factor(lab, levels = lab)]
pa <- ggplot(ca, aes(spearman_sig, lab, fill = kind)) +
  geom_col(width = 0.72) +
  geom_text(aes(label = sprintf("%.3f", spearman_sig)), hjust = 1.1, size = 2.2, colour = "white", fontface = "bold") +
  geom_text(aes(x = 0.01, label = sprintf("%.0f%% kept", pct_retained)), hjust = 0, size = 2.0, colour = "grey92") +
  scale_fill_manual(values = KCOL, guide = "none") +
  coord_cartesian(xlim = c(0, 1)) +
  labs(title = "Mediation-effect concordance with base", subtitle = NULL,
       x = expression(Spearman~rho), y = NULL) +
  theme_heap(base_size = 8) + theme(plot.subtitle = element_text(size = 6, colour = "grey35"))

## b -- # FDR-significant mediations per spec (with disease/case context) -----
cb <- ord(copy(summ)); cb[, lab := factor(lab, levels = lab)]
cb[, ctx := sprintf("%d dz · %.0fk cases", n_diseases, total_cases/1000)]
pb <- ggplot(cb, aes(n_sig, lab, fill = kind)) +
  geom_col(width = 0.72) +
  geom_text(aes(label = formatC(n_sig, big.mark = ",", format = "d")), hjust = 1.1, size = 2.2, colour = "white", fontface = "bold") +
  geom_text(aes(x = max(n_sig) * 0.012, label = ctx), hjust = 0, size = 1.85, colour = "grey92") +
  scale_fill_manual(values = KCOL, name = NULL,
    labels = c(primary = "base", covariate = "+ covariate", sample = "sample", model = "model family")) +
  labs(title = "Significant exposomic mediations per spec",
       subtitle = NULL,
       x = "# FDR-significant (protein x disease)", y = NULL) +
  theme_heap(base_size = 8) + theme(plot.subtitle = element_text(size = 6, colour = "grey35"),
        legend.position = c(0.98, 0.98), legend.justification = c(1, 1), legend.key.size = unit(0.3, "cm"), legend.text = element_text(size = 6))

## c -- base vs +clinical NIE attenuation ------------------------------------
cl <- att[spec == "M3_base_clinical_lasso_primary"]
cl[, retained := ifelse(sig_s, "still significant", "lost")]
cl[, drop := abs(logHR_base) - abs(logHR_s)]
labset <- cl[order(-drop)][!duplicated(protID)][1:12]
rng <- range(c(cl$logHR_base, cl$logHR_s), na.rm = TRUE)
pc <- ggplot(cl, aes(logHR_base, logHR_s)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey60", linewidth = 0.3) +
  geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.2) + geom_vline(xintercept = 0, colour = "grey85", linewidth = 0.2) +
  geom_point(aes(colour = retained), size = 0.5, alpha = 0.35) +
  geom_point(data = labset, colour = "black", size = 0.9) +
  ggrepel::geom_text_repel(data = labset, aes(label = protID), size = 2.0, max.overlaps = 20,
                           bg.color = "white", bg.r = 0.12, segment.size = 0.2, min.segment.length = 0) +
  scale_colour_manual(values = c("still significant" = "grey55", "lost" = "#B2182B"), name = NULL) +
  coord_equal(xlim = rng, ylim = rng) +
  labs(title = "Attenuation under BMI/clinical adjustment", subtitle = NULL,
       x = "NIE log-HR  (base)", y = "NIE log-HR  (+clinical)") +
  theme_heap(base_size = 8) + theme(plot.subtitle = element_text(size = 6, colour = "grey35"),
        legend.position = c(0.02, 0.98), legend.justification = c(0, 1), legend.key.size = unit(0.3, "cm"), legend.text = element_text(size = 6))

## d -- per-category attenuation under +clinical -----------------------------
if (!is.null(catt) && nrow(catt)) {
  catt <- catt[n_base_sig >= 100][order(median_ratio)]
  catt[, catf := factor(heap_category_pretty(category), levels = heap_category_pretty(category[order(median_ratio)]))]
  pal <- HEAP_ECAT_COLORS; names(pal) <- heap_category_pretty(names(pal))
  pd <- ggplot(catt, aes(median_ratio, catf, fill = catf)) +
    geom_vline(xintercept = 1, linetype = 2, colour = "grey55", linewidth = 0.3) +
    geom_col(width = 0.72) +
    geom_text(aes(label = sprintf("%.2f", median_ratio)), hjust = -0.15, size = 2.0, colour = "grey25") +
    scale_fill_manual(values = pal, guide = "none") +
    coord_cartesian(xlim = c(0, max(catt$median_ratio) * 1.12)) +
    labs(title = "Which exposure axes attenuate under BMI/clinical adjustment", subtitle = NULL,
         x = "median |NIE| ratio  (+clinical / base)", y = NULL) +
    theme_heap(base_size = 8) + theme(plot.subtitle = element_text(size = 6, colour = "grey35"),
          axis.text.y = element_text(size = 6.5))
} else pd <- patchwork::plot_spacer()

p <- (pa | pb) / (pc | pd) + plot_layout(heights = c(1, 1.15)) +
  plot_annotation(tag_levels = "a") & theme(plot.tag = element_text(face = "bold", size = 11))

heap_emit_figure(p, figure_id, data = summ, category = "supplement", website = FALSE, width = 9.0, height = 7.6)
