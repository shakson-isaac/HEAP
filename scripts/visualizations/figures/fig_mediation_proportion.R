#!/usr/bin/env Rscript

# ============================================================================
# fig_mediation_proportion.R  [figure_id: fig_mediation_proportion]
# ----------------------------------------------------------------------------
# CONSOLIDATED "how much of a disease effect actually flows through a protein?"
# Merges:
#   fig_mediation_pm_landscape    (PM vs total effect, coarse cis/trans/exposomic)
#   fig_mediation_pm_by_category  (PM violin per fine driver)
#
# The claim it serves, verbatim from the results: "Individual proteins explained
# only smaller fractions of a single exposure-disease association, suggestive of
# lifestyle influencing disease through shared proteomic responses rather than
# singular mediators."
#   a  the landscape: PM is bounded low for exposomic links across the whole
#      total-effect range -- no exposure-disease pair is carried by one protein
#   b  per driver: EVERY lifestyle/environmental driver sits at a median PM of
#      ~11-24%, while the protein's own cis genetic score runs at ~54% (the
#      positive control -- a cis-pQTL acts on disease essentially only through
#      its protein, so PM there SHOULD be high). That contrast is the point.
#
# Retired alongside this: fig_mediation_main (a volcano of mediated EFFECT SIZE,
# which is not what the sentence claims), fig_mediation_scale and
# fig_mediation_category_heatmap (both restate main Fig 3b; readers go to the
# per-disease supplementary table instead).
#
# Input : module3/<exp>/<covarType>/<family>/primary_total/       (exposomic PM)
#         module3/<exp>/<covarType>/<family>/partitioned_categories/ (genetic + per-category PM)
# Output: figures/supplement/module3/fig_mediation_proportion.{pdf,png} + data tsv
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
  library(data.table); library(ggplot2); library(ggrepel); library(patchwork)
})

a <- commandArgs(trailingOnly = TRUE); a <- a[!startsWith(a, "--")]
a <- a[!a %in% c("fig_mediation_proportion", "all_main", "all_supplement", "all", "website")]
covarType <- if (length(a) >= 1) a[1] else "base"
family    <- if (length(a) >= 2) a[2] else "lasso"
figure_id <- Sys.getenv("HEAP_FIGURE_ID", unset = "fig_mediation_proportion")

BS      <- 7.5
ALPHA   <- as.numeric(Sys.getenv("HEAP_MED_ALPHA", unset = "0.05"))
MIN_N   <- as.integer(Sys.getenv("HEAP_MED_MINN",  unset = "30"))
PM_LABELS <- "LEP:obesity,FABP4:chronic renal,IGSF9:non insulin,PRSS8:interstitial"

SEL <- c("protID", "DZ_ID", "predictor", "predictor_class", "effect_type",
         "effect_logHR", "effect_HR", "delta_p", "instrument_present",
         "n", "n_cases", "protein_HR", "protein_p")

# ---------------------------------------------------------------- data ------
# Exposomic PM comes from the primary model (the TOTAL exposome, PXS_total);
# genetic PM is split cis/trans in the partitioned model. FDR is taken over each
# model's own NIE family before the two are combined.
pri     <- load_module3_results(covarType = covarType, family = family,
                                mode = "primary_total")
pm_expo <- heap_proportion_mediated(pri)[predictor_class == "exposure_total"]

par <- load_module3_results(covarType = covarType, family = family,
                            mode = "partitioned_categories", select = SEL)
# Genetic rows with no cis/trans instrument are structural zeros (NIE=0), not
# evidence of "no mediation" -- drop them so they don't deflate the genetic PM.
par    <- par[!(predictor_class %in% c("genetic_cis", "genetic_trans") &
                instrument_present == FALSE)]
pm_par <- heap_proportion_mediated(par)

keep_sig <- function(d) d[is.finite(NIE_q) & NIE_q < ALPHA &
                          pm_consistent == TRUE & is.finite(pm_display)]

# ---- a: coarse landscape (PM vs total effect) -------------------------------
pmA <- rbind(pm_expo,
             pm_par[predictor_class %in% c("genetic_cis", "genetic_trans")],
             fill = TRUE)
pmA[, driver := heap_md_driver_group(predictor_class)]
pmA <- pmA[driver %in% c("Genetic (cis)", "Genetic (trans)", "Exposomic (PXS)")]
pmA[, driver := factor(driver, levels = c("Genetic (cis)", "Genetic (trans)",
                                          "Exposomic (PXS)"))]
pmA <- keep_sig(pmA)[is.finite(TE_HR)]
pmA[, disease := heap_pretty_disease(DZ_ID)]

share <- pmA[, .N, by = driver][, pct := round(100 * N / sum(N))][]
message(sprintf("PM landscape: %d links | %s", nrow(pmA),
                paste(sprintf("%s=%d(%d%%)", share$driver, share$N, share$pct), collapse = " ")))

want   <- strsplit(strsplit(PM_LABELS, ",")[[1]], ":")
lab_dt <- rbindlist(lapply(want, function(w)
  pmA[protID == trimws(w[1]) & grepl(trimws(w[2]), disease, ignore.case = TRUE)][
      order(-abs(NIE_logHR))][1]), fill = TRUE)
lab_dt <- lab_dt[!is.na(protID)]
lab_dt[, lbl := sprintf("%s → %s", protID, stringr::str_trunc(disease, 28))]

pal_drv <- heap_md_driver_group_colors()[c("Genetic (cis)", "Genetic (trans)",
                                           "Exposomic (PXS)")]

pa <- ggplot(pmA, aes(TE_HR, pm_display, colour = driver)) +
  geom_hline(yintercept = c(0.25, 0.5, 0.75), linetype = "dotted", colour = "grey88",
             linewidth = .25) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey60", linewidth = .25) +
  geom_point(aes(size = n_cases), alpha = .45, stroke = 0) +
  geom_text_repel(data = lab_dt, aes(label = lbl), size = 1.9, max.overlaps = 20,
                  min.segment.length = 0, segment.colour = "grey55", segment.size = .25,
                  box.padding = .4, bg.color = "white", bg.r = .12, seed = 5,
                  show.legend = FALSE) +
  scale_colour_manual(values = pal_drv, name = "Mediated driver") +
  scale_size_continuous(range = c(.4, 3.2), name = "Incident cases",
                        breaks = c(200, 500, 1000, 2000)) +
  scale_x_continuous(trans = "log10") +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(.01, .03)),
                     labels = scales::percent) +
  labs(x = "Total effect on disease (HR, log scale)",
       y = "Proportion mediated (NIE / total)") +
  theme_heap(base_size = BS) +
  theme(panel.grid = element_blank(),
        legend.position = "right", legend.box = "vertical",
        legend.key.size = unit(8, "pt"),
        legend.text  = element_text(size = BS - 1.5),
        legend.title = element_text(size = BS - 0.5),
        legend.spacing.y = unit(2, "pt"),
        plot.margin = margin(10, 4, 2, 2)) +
  guides(colour = guide_legend(override.aes = list(size = 2, alpha = .9), order = 1),
         size = guide_legend(order = 2))

# ---- b: PM per fine driver --------------------------------------------------
pmB <- keep_sig(pm_par)
pmB[, category := heap_md_category(predictor)]
pmB[, driver := fifelse(!is.na(category), category,
             fifelse(predictor_class == "genetic_cis",   "Genetic (cis)",
             fifelse(predictor_class == "genetic_trans", "Genetic (trans)", NA_character_)))]
pmB <- pmB[!is.na(driver)]
pmB <- pmB[driver %in% pmB[, .N, by = driver][N >= MIN_N, driver]]

ord <- pmB[, .(med = median(pm_display), n = .N), by = driver][order(med)]
pmB[, driver := factor(driver, levels = ord$driver)]

pal_fine <- c(HEAP_ECAT_COLORS,
              "Genetic (cis)"   = unname(pal_drv[["Genetic (cis)"]]),
              "Genetic (trans)" = unname(pal_drv[["Genetic (trans)"]]))

ann <- ord[, .(driver = factor(driver, levels = ord$driver),
               lab = sprintf("%d%%  (n=%s)", round(100 * med),
                             formatC(n, big.mark = ",", format = "d")))]

message(sprintf("PM by driver: %d drivers | %d links | exposomic medians %s",
                uniqueN(pmB$driver), nrow(pmB),
                paste(sprintf("%.0f%%", 100 * ord[!grepl("Genetic", driver), med]), collapse = "/")))

pb <- ggplot(pmB, aes(pm_display, driver, fill = as.character(driver))) +
  geom_violin(scale = "width", colour = NA, alpha = .55, width = .95) +
  geom_boxplot(width = .16, outlier.shape = NA, fill = "white",
               colour = "grey25", linewidth = .25) +
  geom_text(data = ann, aes(x = 1.03, y = driver, label = lab), inherit.aes = FALSE,
            hjust = 0, size = 1.9, colour = "grey25") +
  scale_fill_manual(values = pal_fine, guide = "none") +
  scale_y_discrete(labels = heap_category_pretty) +
  scale_x_continuous(limits = c(0, 1.30), breaks = c(0, .25, .5, .75, 1),
                     labels = scales::percent) +
  labs(x = "Proportion mediated (NIE / total effect)", y = NULL) +
  theme_heap(base_size = BS) +
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(size = BS - 1),
        plot.margin = margin(10, 4, 2, 2))

# ------------------------------------------------------------- assemble ------
p <- (pa / pb) + plot_layout(heights = c(1, 1.15)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 9),
        plot.tag.position = c(0, 1))

out <- rbindlist(list(
  pmA[, .(panel = "a", protID, DZ_ID, disease, driver = as.character(driver),
          TE_HR, pm_display, NIE_q, n_cases)],
  pmB[, .(panel = "b", protID, DZ_ID, disease = NA_character_,
          driver = as.character(driver), TE_HR = NA_real_, pm_display, NIE_q, n_cases)]),
  use.names = TRUE)

heap_emit_figure(p, figure_id, data = out, category = "supplement", subdir = "module3",
                 formats = c("pdf", "png"), width = 6.5, height = 7.2, website = TRUE)

message("fig_mediation_proportion: done (", covarType, "/", family, ").")
