#!/usr/bin/env Rscript

# ============================================================================
# fig_traintest_stability_categories.R  [figure_id: fig_traintest_stability_categories]
# ----------------------------------------------------------------------------
# Exposome generalisation, resolved to the 13 fine EXPOSURE CATEGORIES: per
# protein, TRAIN vs TEST (out-of-fold) unique drop-one R2 for each category's
# PXS, faceted by category. Points on y = x replicate; a category whose points
# sit below the line contributes train-only (non-replicating) variance. This is
# the per-category counterpart of fig_traintest_stability_components and the
# direct check that the exposure-type-routing fix removed the spurious
# train-only signal (e.g. the deprivation/IMD polynomial terms).
#
# Source: predictive_r2_exposure_categories (score_unique_drop per PXS_<cat>),
# which holds out-of-fold `r2` and in-sample `r2_train`; averaged over folds.
#
# Input : module1_predictive_r2_score_partition/<exp>/<covarType>/<method>/predictive_r2_exposure_categories_*
# Output: figures/supplement/fig_traintest_stability_categories.{pdf,png} + data tsv
#
# Run:
#   HEAP_PATHS_FILE=.../workflow/00_paths.R \
#     Rscript scripts/visualizations/figures/fig_traintest_stability_categories.R [covarType] [method] [experiment]
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
a <- a[!a %in% c("fig_traintest_stability_categories", "all_main", "all_supplement", "all", "website")]
covarType  <- if (length(a) >= 1) a[1] else "base"
method     <- if (length(a) >= 2) a[2] else "lasso"
experiment <- if (length(a) >= 3) a[3] else "M1_base_lasso"
figure_id  <- Sys.getenv("HEAP_FIGURE_ID", unset = "fig_traintest_stability_categories")

pr <- load_module1_predictive_r2(covarType = covarType, method = method,
                                 level = "exposure_categories", experiment = experiment)
if (!"r2_train" %in% names(pr))
  stop("predictive_r2_exposure_categories has no r2_train column — re-run Module 1.")
pr <- pr[is.finite(r2) & is.finite(r2_train)]
pr[, category := sub("^PXS_", "", block)]

# per protein x category: mean over folds
agg <- pr[, .(train_r2 = mean(r2_train), test_r2 = mean(r2)), by = .(omic, category)]
agg[, category := heap_category_factor(category)]

# per-category generalisation stats (correlation over proteins with any signal)
stats <- agg[, {
  nz <- abs(train_r2) > 1e-9 | abs(test_r2) > 1e-9
  .(rho = if (sum(nz) > 2) cor(train_r2[nz], test_r2[nz], use = "complete.obs") else NA_real_,
    med_gap = median(train_r2 - test_r2, na.rm = TRUE), n_signal = sum(nz))
}, by = category]
stats[, lab := sprintf("r=%.2f  n=%d", rho, n_signal)]
setorder(stats, category); message("per-category stability:"); print(stats[, .(category, rho = round(rho,3), n_signal)])

# equal x/y range per facet so y = x is the diagonal
rng <- agg[, .(lo = min(c(train_r2, test_r2)), hi = max(c(train_r2, test_r2))), by = category]
frame <- rbind(rng[, .(category, v = lo)], rng[, .(category, v = hi)])

p <- ggplot(agg, aes(train_r2, test_r2)) +
  geom_blank(data = frame, aes(v, v), inherit.aes = FALSE) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey45", linewidth = 0.3) +
  geom_point(aes(colour = category), size = 0.7, alpha = 0.5, stroke = 0, show.legend = FALSE) +
  geom_text(data = stats, aes(x = -Inf, y = Inf, label = lab), inherit.aes = FALSE,
            hjust = -0.06, vjust = 1.3, size = 2.5, colour = "grey20") +
  scale_colour_exposure(drop = FALSE) +
  facet_wrap(~ category, scales = "free", labeller = as_labeller(heap_category_pretty)) +
  labs(title = "Train vs test (out-of-fold) R2 by exposure category",
       subtitle = NULL,
       x = expression("Train"~R^2~"(PXS unique, drop-one)"),
       y = expression("Test (out-of-fold)"~R^2)) +
  theme_heap() +
  theme(panel.spacing = unit(0.6, "lines"),
        strip.text = element_text(size = 7.5),
        axis.text = element_text(size = 6))

heap_emit_figure(p, figure_id, data = agg, category = "supplement",
                 formats = c("pdf", "png"), width = 10, height = 8, website = TRUE)
message("fig_traintest_stability_categories: done (", experiment, "/", covarType, "/", method, ").")
