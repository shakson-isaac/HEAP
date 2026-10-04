#!/usr/bin/env Rscript

# ============================================================================
# fig_pes_covariate_robustness.R  [figure_id: fig_pes_covariate_robustness]
# ----------------------------------------------------------------------------
# CONSOLIDATED Module-6 covariate-sensitivity supplement. Merges the three
# one-off builders that all argued the same thing:
#   build_covsens_supp1_adds.R      "what the PES adds beyond covariates"
#   build_covsens_supp2_survives.R  "the increment survives stricter covariates"
#   build_covsens_supp3_robust.R    "robust to covariates, not a prevalence artifact"
#
# supp1's dumbbell panels were just the top-12 slice of fig_pes_incremental_value
# (which already shows EVERY exposure), and its AUPR panel duplicated supp3's.
# supp2's spaghetti and supp3's concordance scatter are the same claim at two
# resolutions, so both are kept -- one as the per-exposure trajectory, one as its
# summary.
#
# The citing sentence has three clauses; this figure has one panel group each:
#   a-c  "sustained over multiple covariate specifications ... and evaluation
#         metrics including AUPR": the increment across the
#         base -> +draw -> +BMI -> +clinical ladder, on R2, AUC and AUPR, with a
#         detached healthy-at-baseline column after the dashed rule
#   d    the same claim summarized: base vs +clinical increment (Spearman ~0.96)
#   e    "not a prevalence artifact": held-out AUPR against its prevalence null
#   f    the DISEASE arm under the healthy-at-baseline restriction: PES-disease
#        log hazard ratio, all participants vs healthy at baseline
# The MAGNITUDE of the increment lives in fig_pes_incremental_value (every exposure).
#
# NB built from the sensitivity table DIRECTLY rather than by sourcing the atomic
# panel plotters (fig_pes_covsens_*.R): those bake geom-level text sizes sized for
# an 8in canvas, which at the supplement's 6.5in came out enormous and clipped.
#
# Input : module6_pes_longitudinal/covariate_sensitivity/covariate_sensitivity.tsv
#         module6_pes_longitudinal/multipes_disease/pes_disease_scale.tsv           (f)
#         module6_pes_longitudinal/multipes_disease/pes_disease_scale_base_exclprev.tsv (f)
# Output: figures/supplement/module6/fig_pes_covariate_robustness.{pdf,png} + data tsv
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

figure_id <- Sys.getenv("HEAP_FIGURE_ID", unset = "fig_pes_covariate_robustness")
BS  <- 7.5
LBL <- 1.7    # repel label size (mm) at true print scale

# ---------------------------------------------------------------- data ------
root <- dirname(heap_project_output("module6_pes_longitudinal", "base"))
d <- fread(file.path(root, "covariate_sensitivity", "covariate_sensitivity.tsv"))
d[, category := heap_category_factor(category)]

# ORDER IS AN ARGUMENT. base -> +draw -> +BMI -> +clinical keeps adjustment
# monotone AND puts the negative control immediately before the test rung:
# +draw adds fasting time and season, i.e. more covariates and a slightly smaller
# complete-case sample with no adiposity. It lands flat (median increment ratio
# to base 1.003) where +BMI dips (0.972), so the dip is attributable to adiposity
# rather than to adjustment per se -- a reading the figure could not support
# before, because the control rung was not drawn.
#
# The fifth column is NOT a rung: base_exclprev restricts to participants free of
# major disease at baseline AND uses the score retrained on that subset, so
# sample and score both move. Lines run through it because it reads better as one
# trajectory, which makes the dashed rule the only marker of that break -- the
# caption must therefore say so explicitly.
LAD <- c("base", "base_draw", "base_bmi", "base_clinical", "base_exclprev")
LL  <- c("base", "+ draw", "+ BMI", "+ clinical", "healthy at baseline")
RESTRICT_AT <- 4.5   # rule between the covariate ladder and the restriction

# ---- a-c: the increment across the covariate ladder -------------------------
# One line per exposure. Flat = the PES increment is not covariate-mediated;
# a steep drop = the covariates absorbed the signal (BMI is on the pathway).
ladder <- function(type, ycol, ylab, thr, show_cats = FALSE) {
  x <- d[spec %in% LAD & exposure_type == type & is.finite(get(ycol))]
  keep <- x[spec == "base" & get(ycol) > thr, exposure_id]
  x <- x[exposure_id %in% keep]
  x[, sx := factor(LL[match(spec, LAD)], levels = LL)]

  lab <- x[sx == "+ clinical"][order(-get(ycol))][seq_len(min(4, .N))]
  lab[, elab := heap_exposure_label(exposure_id)]

  ggplot(x, aes(sx, get(ycol), group = exposure_id, colour = category)) +
    geom_vline(xintercept = RESTRICT_AT, linetype = "22",
               colour = "grey70", linewidth = .3) +
    geom_line(alpha = .3, linewidth = .3) +
    geom_point(size = .7, alpha = .55) +
    geom_text_repel(data = lab, aes(label = elab), colour = "grey20", size = LBL,
                    direction = "y", hjust = 1, nudge_x = -.30, segment.colour = "grey70",
                    segment.size = .18, min.segment.length = 0, max.overlaps = 30,
                    seed = 6, show.legend = FALSE) +
    scale_colour_exposure(drop = FALSE,
                          guide = if (show_cats) "legend" else "none") +
    guides(colour = if (show_cats)
             guide_legend(override.aes = list(linewidth = 0, size = 1.4, alpha = 1))
           else "none") +
    scale_x_discrete(expand = expansion(mult = c(.14, .16))) +
    labs(x = NULL, y = ylab) +
    theme_heap(base_size = BS) +
    theme(panel.grid = element_blank(),
          # five rungs do not fit horizontally at 1/3 of \textwidth
          axis.text.x = element_text(size = BS - 2.4, angle = 38, hjust = 1),
          plot.margin = margin(12, 3, 2, 4))
}
pa <- ladder("continuous", "incremental",      "R² increment",   0.03)
pb <- ladder("binary",     "incremental",      "AUC increment",  0.05)
pc <- ladder("binary",     "incremental_aupr", "AUPR increment", 0.05)

# ---- d: base vs +clinical increment (the ladder, summarized) ----------------
w <- dcast(d[spec %in% c("base", "base_clinical")],
           exposure_id + exposure_type + category ~ spec, value.var = "incremental")
setnames(w, c("base", "base_clinical"), c("inc_base", "inc_clin"))
w <- w[is.finite(inc_base) & is.finite(inc_clin)]
rho <- cor(w$inc_base, w$inc_clin, method = "spearman")

EX <- c("alcohol_intake_frequency_f1558_0_0",
        "number_of_days_week_of_vigorous_physical_activity_10_plus_minutes_f904_0_0",
        "processed_meat_intake_f1349_0_0", "usual_walking_pace_f924_0_0",
        "coffee_intake_f1498_0_0")
EL <- c("Alcohol frequency", "Vigorous activity", "Processed meat", "Walking pace", "Coffee")
ex <- w[exposure_id %in% EX]; ex[, elab := EL[match(exposure_id, EX)]]
lim <- range(c(w$inc_base, w$inc_clin), na.rm = TRUE)

pd <- ggplot(w, aes(inc_base, inc_clin)) +
  geom_abline(slope = 1, intercept = 0, linetype = "22", colour = "grey60", linewidth = .3) +
  geom_hline(yintercept = 0, colour = "grey88", linewidth = .25) +
  geom_point(aes(colour = category, shape = exposure_type), size = 1.1, alpha = .75) +
  geom_point(data = ex, shape = 21, fill = "white", colour = "grey15", size = 1.5, stroke = .5) +
  geom_text_repel(data = ex, aes(label = elab), colour = "grey15", size = LBL,
                  box.padding = .45, min.segment.length = 0, segment.colour = "grey55",
                  segment.size = .18, max.overlaps = 30, seed = 4, show.legend = FALSE) +
  annotate("text", x = lim[1], y = lim[2], hjust = 0, vjust = 1,
           label = sprintf("Spearman = %.2f", rho), size = 2.1, colour = "grey30") +
  scale_colour_exposure(drop = TRUE) +
  guides(colour = guide_legend(override.aes = list(size = 1.4, alpha = 1),
                               order = 1, nrow = 3, byrow = TRUE),
         shape  = guide_legend(order = 2, nrow = 1)) +
  # NB: labels must be NAMED -- an unnamed vector binds positionally to the
  # alphabetically-sorted breaks (binary, continuous) and inverts the legend.
  scale_shape_manual(NULL, values = c(continuous = 16, binary = 17),
                     labels = c(binary = "binary (AUC)", continuous = "continuous (R²)")) +
  coord_equal(xlim = lim, ylim = lim) +
  labs(x = "Increment, base covariates",
       y = "Increment, + clinical") +
  theme_heap(base_size = BS) +
  theme(panel.grid = element_blank(),
        plot.margin = margin(10, 4, 2, 4))

# ---- e: AUPR against its prevalence null ------------------------------------
# For a rare binary exposure a high AUC can still be uninformative; the honest
# null for AUPR is the class prevalence itself (dashed curve).
e <- d[spec == "base" & exposure_type == "binary" &
       is.finite(aupr_full) & is.finite(prevalence)]
e[, signal := fifelse(aupr_lift > 0.02, "signal (AUPR > prevalence)",
                                        "null (AUPR ~ prevalence)")]
e[, elab := heap_exposure_label(exposure_id)]
lab_e <- e[order(-aupr_lift)][seq_len(min(4, .N))]
nullc <- data.table(x = exp(seq(log(min(e$prevalence)), log(max(e$prevalence)),
                                length.out = 200)))[, y := x]

pe <- ggplot(e, aes(prevalence, aupr_full)) +
  geom_line(data = nullc, aes(x, y), linetype = "22", colour = "grey55", linewidth = .35) +
  geom_segment(aes(xend = prevalence, y = aupr_cov, yend = aupr_full),
               colour = "grey82", linewidth = .3) +
  geom_point(aes(y = aupr_cov), colour = "grey65", size = .8, alpha = .6) +
  geom_point(aes(colour = signal), size = 1.2, alpha = .85) +
  geom_text_repel(data = lab_e, aes(label = elab), size = LBL, colour = "grey20",
                  max.overlaps = 40, min.segment.length = 0, segment.colour = "grey70",
                  segment.size = .18, box.padding = .4, seed = 7, show.legend = FALSE) +
  scale_colour_manual(NULL, values = c("signal (AUPR > prevalence)" = "#1B7837",
                                       "null (AUPR ~ prevalence)"   = "#9E9E9E")) +
  scale_x_log10() +
  labs(x = "Positive-class prevalence (log)", y = "Held-out AUPR, cov + PES") +
  theme_heap(base_size = BS) +
  theme(panel.grid = element_blank(), legend.position = "bottom",
        legend.key.size = unit(6, "pt"), legend.text = element_text(size = BS - 2),
        plot.margin = margin(10, 4, 2, 4))

# ---- f: the DISEASE arm under the same healthy-at-baseline restriction ------
# In HAZARD RATIO, not delta-C. Both were computed; the HR is a one-degree-of-
# freedom coefficient fit on the full cohort with a closed-form SE, whereas
# held-out delta-C comes from a single 70/30 split and is noisy per cell. The
# same comparison gives Spearman 0.95 in log HR and only 0.67 in delta-C, and
# that gap is estimator noise rather than biology (the restricted arm has fewer
# events per cell, which both attenuates delta-C and degrades its rank order).
# The delta-C columns ship in the deposit for anyone who wants them.
MDZ <- heap_project_output("module6_pes_longitudinal", "multipes_disease")
dzb <- fread(file.path(MDZ, "pes_disease_scale.tsv"),
             select = c("exposure_id", "disease", "hr_pes"))
dzx <- fread(file.path(MDZ, "pes_disease_scale_base_exclprev.tsv"),
             select = c("exposure_id", "disease", "hr_pes"))
DZ <- merge(dzb, dzx, by = c("exposure_id", "disease"), suffixes = c("_all", "_hb"))
DZ <- DZ[is.finite(hr_pes_all) & is.finite(hr_pes_hb) & hr_pes_all > 0 & hr_pes_hb > 0]
DZ[, `:=`(lp_all = log(hr_pes_all), lp_hb = log(hr_pes_hb))]
if (!nrow(DZ)) stop("panel f: no shared exposure-disease pairs across the two scans")
rho_dz  <- cor(DZ$lp_all, DZ$lp_hb, method = "spearman")
dzlim   <- quantile(c(DZ$lp_all, DZ$lp_hb), c(.002, .998), na.rm = TRUE)

pf <- ggplot(DZ, aes(lp_all, lp_hb)) +
  geom_abline(slope = 1, intercept = 0, linetype = "22", colour = "grey55", linewidth = .3) +
  geom_hline(yintercept = 0, colour = "grey88", linewidth = .25) +
  geom_vline(xintercept = 0, colour = "grey88", linewidth = .25) +
  geom_hex(bins = 46) +
  scale_fill_gradient("Pairs", low = "#DCE6F2", high = "#1F3A6E", trans = "log10",
                      breaks = c(1, 10, 100, 1000)) +
  annotate("text", x = dzlim[1], y = dzlim[2], hjust = 0, vjust = 1, size = 2.1,
           colour = "grey25", lineheight = .95,
           label = sprintf("Spearman = %.2f\n%s pairs", rho_dz,
                           format(nrow(DZ), big.mark = ","))) +
  coord_equal(xlim = dzlim, ylim = dzlim) +
  labs(x = "log HR, all participants", y = "log HR, healthy at baseline") +
  theme_heap(base_size = BS) +
  theme(panel.grid = element_blank(), plot.margin = margin(10, 4, 2, 4),
        legend.position = "bottom", legend.direction = "horizontal",
        legend.key.width = unit(.85, "cm"), legend.key.height = unit(.16, "cm"),
        legend.title = element_text(size = BS - 2),
        legend.text = element_text(size = BS - 2.5))

message(sprintf("PES covariate robustness: rho(base, +clinical) = %.3f | %d/%d binary exposures beat their prevalence null | disease rho(all, healthy) = %.3f on %d pairs",
                rho, e[aupr_lift > 0.02, .N], nrow(e), rho_dz, nrow(DZ)))

# ------------------------------------------------------------- assemble ------
# The exposure-category legend is shared by a-d; collect it (and d's shape key)
# into one right-hand column so the panels keep equal width.
# The category key runs as a horizontal strip along the bottom. As a right-hand
# column its 13 entries are taller than the panels and force a dead band between
# the rows, which at three panels per row leaves d-f too narrow to read.
# NB no plot.tag.position: it is in plot-relative units, so it tracks the margin
# and lands on top of the y-axis title once the panels narrow.
p <- ((pa | pb | pc) / (pd | pe | pf)) +
  plot_layout(heights = c(1, 1.25), guides = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(legend.position = "bottom", legend.box = "vertical",
        legend.key.size = unit(6, "pt"),
        legend.text  = element_text(size = BS - 2),
        legend.title = element_text(size = BS - 1),
        legend.margin = margin(1, 4, 0, 4),
        legend.spacing.y = unit(4, "pt"),
        plot.tag = element_text(face = "bold", size = 9),
        plot.tag.location = "margin")

out <- rbindlist(list(
  d[spec %in% LAD, .(panel = "a-c", exposure_id, exposure_type,
                     category = as.character(category), spec,
                     value = incremental, value_aupr = incremental_aupr)],
  w[, .(panel = "d", exposure_id, exposure_type, category = as.character(category),
        spec = "base_vs_clinical", value = inc_base, value_aupr = inc_clin)],
  e[, .(panel = "e", exposure_id, exposure_type, category = as.character(category),
        spec = "base", value = aupr_full, value_aupr = prevalence)],
  # panel f is per exposure x DISEASE, so `disease` is carried and the exposure
  # columns the other panels use do not apply
  DZ[, .(panel = "f", exposure_id, disease, spec = "base_vs_exclprev",
         value = lp_all, value_aupr = lp_hb)]), use.names = TRUE, fill = TRUE)

heap_emit_figure(p, figure_id, data = out, category = "supplement", subdir = "module6",
                 formats = c("pdf", "png"), width = 6.5, height = 5.6, website = TRUE)

message("fig_pes_covariate_robustness: done.")
