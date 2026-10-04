#!/usr/bin/env Rscript

# ============================================================================
# fig_gxe_noise_floor.R  [figure_id: fig_gxe_noise_floor]
# ----------------------------------------------------------------------------
# CONSOLIDATED Module-1 GxE evidence. Replaces:
#   fig_traintest_stability_components  (train vs test scatter by component)
#   fig_interaction_sensitivity         (GxE across age/sex-interaction models)
# and RETIRES fig_gxe_proteins (a top-25 GxE leaderboard, which ranked the
# extreme right tail of a noise distribution as if the proteins were findings --
# it argued AGAINST the claim the sentence it supported was making).
#
# The claim is that the GxE component sits at a noise floor. Two independent
# failures of reproducibility make that case, and this figure shows both:
#   a  WITHIN-method: GxE does not survive train -> test. Genetic holds
#      (r = 0.99) and Exposome mostly holds (r = 0.96), but GxE degrades
#      (r = 0.76) with a systematic positive train-test gap -- i.e. the train
#      signal is fitted noise that does not carry out of fold.
#   b  The estimate is nonetheless STABLE across interaction models (rho 0.92-
#      0.98 vs base). This is the honest framing: GxE is a *consistently
#      estimated* noise floor, not a fragile artifact of unmodelled age/sex
#      interactions. It pre-empts the obvious reviewer question without
#      defending GxE's reality.
#
# The CROSS-method failure (GREML vs HEAP GxE rho = 0.06 -- the two methods
# emit GxE numbers that are uncorrelated with each other) is deliberately NOT
# duplicated here: it is panel c of fig_greml_concordance, which the same
# sentence already cites.
#
# Input : module1_predictive_r2_score_partition/M1_base_lasso{,_GxC,_ExC,_GxC_ExC}/
#           base/lasso/predictive_r2_coarse_*  (needs the r2_train column)
# Output: figures/supplement/module1/fig_gxe_noise_floor.{pdf,png} + data tsv
#
# Authored at 6.5in = the supplement's \textwidth (scale 1.0).
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
  library(data.table); library(ggplot2); library(patchwork)
})

figure_id <- Sys.getenv("HEAP_FIGURE_ID", unset = "fig_gxe_noise_floor")
BS  <- 8
THR <- 0.01

PC <- HEAP_PAL_COMPONENT
PAL <- c(Genetic = PC[["Genetic"]], Exposome = PC[["Exposome"]], GxE = PC[["GxE"]])

# ------------------------------------------- a: train vs test by component ---
pr <- load_module1_predictive_r2(covarType = "base", method = "lasso",
                                 level = "coarse", experiment = "M1_base_lasso")
if (!"r2_train" %in% names(pr))
  stop("predictive_r2_coarse has no r2_train column -- re-run Module 1.")

ud  <- pr[get("method") == "score_unique_drop" & block %in% c("G", "E", "GxE")]
ud  <- ud[is.finite(r2) & is.finite(r2_train)]
agg <- ud[, .(train_r2 = mean(r2_train), test_r2 = mean(r2)), by = .(omic, block)]
agg[, component := factor(block, levels = c("G", "E", "GxE"),
                          labels = c("Genetic", "Exposome", "GxE"))]

st <- agg[, .(r = cor(train_r2, test_r2, use = "complete.obs"),
              med_gap = median(train_r2 - test_r2, na.rm = TRUE),
              n = .N), by = component]
st[, lab := sprintf("r = %.2f", r)]   # median gap dropped from the panel (still printed below)

pa <- ggplot(agg, aes(train_r2, test_r2, colour = component)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              colour = "grey55", linewidth = .3) +
  geom_point(size = .5, alpha = .35) +
  geom_text(data = st, aes(x = -Inf, y = Inf, label = lab), inherit.aes = FALSE,
            hjust = -0.06, vjust = 1.15, size = 2.3, colour = "grey20", lineheight = 1.05) +
  facet_wrap(~ component, scales = "free", nrow = 1) +
  scale_colour_manual(values = PAL, guide = "none") +
  labs(x = expression("Train unique"~R^2~"(drop-one)"),
       y = expression("Test (out-of-fold)"~R^2)) +
  theme_heap(base_size = BS) +
  theme(panel.grid = element_blank(),
        strip.text = element_text(face = "bold", size = BS),
        plot.margin = margin(10, 4, 2, 4))

# ------------------------- b: GxE across age/sex-interaction models -----------
EXPS <- list(
  list(exp = "M1_base_lasso",        lab = "base"),
  list(exp = "M1_base_lasso_GxC",    lab = "+ G×C\n(gene × age/sex)"),
  list(exp = "M1_base_lasso_ExC",    lab = "+ E×C\n(exposure × age/sex)"),
  list(exp = "M1_base_lasso_GxC_ExC",lab = "+ both"))

getGxE <- function(e) {
  co <- load_module1_predictive_r2("base", "lasso", "coarse", experiment = e)
  co[get("method") == "score_unique_drop" & block == "GxE",
     .(GxE = pmax(0, mean(r2))), by = omic]
}
G <- lapply(EXPS, function(s) getGxE(s$exp))
names(G) <- sapply(EXPS, `[[`, "lab")
base_v <- G[[1]]

cnt <- rbindlist(lapply(names(G), function(n) {
  m <- merge(base_v, G[[n]], by = "omic", suffixes = c("_base", "_x"))
  data.table(lab = n,
             n   = sum(G[[n]]$GxE >= THR),
             rho = cor(m$GxE_base, m$GxE_x, method = "spearman"))
}))
cnt[, lab := factor(lab, levels = names(G))]

pb <- ggplot(cnt, aes(lab, n)) +
  geom_col(width = .62, fill = PC[["GxE"]], alpha = .85) +
  geom_text(aes(label = n), vjust = -0.4, size = 2.6, fontface = "bold", colour = "grey15") +
  # ASCII "rho": the Greek glyph (U+03C1) fails to encode in the pdf device
  geom_text(aes(label = sprintf("rho = %.2f", rho)), y = 16, size = 2.2, colour = "white") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.16))) +
  labs(x = NULL, y = expression("# proteins with GxE"~R^2 >= 0.01)) +
  theme_heap(base_size = BS) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(size = BS - 1),
        plot.margin = margin(10, 4, 2, 4))

# ------------------------------------------------------------- assemble ------
p <- (pa / pb) +
  plot_layout(heights = c(1.25, 0.75)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 9),
        plot.tag.position = c(0, 1))

out <- rbindlist(list(
  agg[, .(panel = "a", key = as.character(component), omic,
          value1 = train_r2, value2 = test_r2)],
  cnt[, .(panel = "b", key = as.character(lab), omic = NA_character_,
          value1 = as.numeric(n), value2 = rho)]), use.names = TRUE)

heap_emit_figure(p, figure_id, data = out, category = "supplement", subdir = "module1",
                 formats = c("pdf", "png"), width = 6.5, height = 5.4, website = TRUE)

message("fig_gxe_noise_floor: done.")
print(st[, .(component, r = round(r, 3), med_gap = round(med_gap, 4), n)])
print(cnt[, .(lab, n, rho = round(rho, 3))])
