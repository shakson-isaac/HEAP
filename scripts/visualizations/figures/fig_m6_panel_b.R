#!/usr/bin/env Rscript
# Panel b (Q1): how well the proteome reads each exposure. Split sub-panels for
# the two metrics -- continuous (R²) and binary (AUC). Every exposure is a dot
# colored by its category (canonical HEAP_ECAT_COLORS); one representative
# exposure per category is labeled. The recurring cast (alcohol, smoking, oily
# fish, vigorous activity, PM2.5) is the labeled exemplar of its category.
local({
  cand <- c(file.path(getwd(),"scripts","visualizations","common"),
            "/n/groups/patel/shakson_ukb/HEAP/scripts/visualizations/common")
  cm <- cand[dir.exists(cand)][1]
  for (f in c("figure_paths.R","load_heap_results.R","plot_theme.R","label_helpers.R","export_helpers.R")) source(file.path(cm,f))
})
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(ggrepel); library(patchwork) })

# --- render profile: standalone (default) vs small composite CELL ----------
CELL  <- nzchar(Sys.getenv("HEAP_CELL"))
BS    <- if (CELL) 7   else 10     # base_size
PSZ_S <- if (CELL) 0.9 else 1.5    # cloud point
PSZ_L <- if (CELL) 2.3 else 2.9    # exemplar point (bolder than the cloud)
LBL   <- if (CELL) 1.8 else 2.5    # exemplar label (mm)
TTL   <- if (CELL) 9   else 13     # plot title (pt)
B_W   <- 5.93; B_H <- 2.45         # cell size (in); fills top-right row

CAST <- c("alcohol_intake_frequency_f1558_0_0","pack_years_of_smoking_f20161_0_0",
          "oily_fish_intake_f1329_0_0",
          "number_of_days_week_of_vigorous_physical_activity_10_plus_minutes_f904_0_0")
PRETTY <- c(Deprivation_Indices="Deprivation / income", Exercise_MET="Exercise (MET)",
            Exercise_Freq="Exercise (frequency)", Sun_Exposure="Sun exposure", Diet_Weekly="Diet",
            Internet_Usage="Internet use", Sexual_Factors="Sexual factors",
            Residential_Air_Pollution="Air pollution", Residential_Noise_Pollution="Noise pollution")
pcat <- function(x) ifelse(x %in% names(PRETTY), PRETTY[x], x)
shorten <- function(x) { x <- gsub("_", " ", x); ifelse(nchar(x) > 26, paste0(substr(x,1,24),"…"), x) }
# curated short labels for the one exemplar shown per category (avoids ugly auto-truncation)
EX_LABEL <- c(
  alcohol_intake_frequency_f1558_0_0 = "Alcohol frequency",
  average_total_household_income_before_tax_f738_0_0 = "Household income",
  oily_fish_intake_f1329_0_0 = "Oily fish",
  processed_meat_intake_f1349_0_0 = "Processed meat",
  number_of_days_week_of_vigorous_physical_activity_10_plus_minutes_f904_0_0 = "Vigorous activity",
  summed_days_activity_f22033_0_0 = "Active days/week",
  weekly_usage_of_mobile_phone_in_last_3_months_f1120_0_0 = "Mobile-phone use",
  pm2_5_mean = "PM2.5",
  particulate_matter_air_pollution_pm10._2007_f24019_0_0 = "PM10",
  no2_mean = "NO2",
  average_night_time_sound_level_of_noise_pollution_f24022_0_0 = "Night-time noise",
  age_first_had_sexual_intercourse_f2139_0_0 = "Age at first sex",
  nap_during_day_f1190_0_0 = "Daytime napping",
  pack_years_of_smoking_f20161_0_0 = "Pack-years smoking",
  time_spend_outdoors_in_summer_f1050_0_0 = "Time outdoors (summer)",
  alcohol_drinker_status_f20117_0_0_Never = "Never drinks",
  spread_type_f1428_0_0_Flora_Pro.Active.Benecol = "Flora/Benecol spread",
  types_of_physical_activity_in_last_4_weeks_f6164_0_0.multi_Strenuous_sports = "Strenuous sport",
  at_or_above_moderate_vigorous_walking_recommendation_f22036_0_0 = "Meets activity guideline",
  difference_in_mobile_phone_use_compared_to_two_years_previously_f1140_0_0_Yes._use_is_now_less_frequent = "Less phone use vs 2y",
  answered_sexual_history_questions_f2129_0_0 = "Answered sexual history",
  snoring_f1210_0_0_Yes = "Snores",
  current_tobacco_smoking_f1239_0_0_Yes._on_most_or_all_days = "Smokes most days",
  mineral_and_other_dietary_supplements_f6179_0_0.multi_Calcium = "Calcium supplement")
pretty_ex <- function(id) ifelse(id %in% names(EX_LABEL), EX_LABEL[id], shorten(heap_exposure_label(id)))
# exposures we never label as a category's exemplar (survey-completion / non-exposure proxies)
EX_EXCLUDE <- c("answered_sexual_history_questions_f2129_0_0")
# forced exemplar per category (overrides cast/highest-value pick) -- steers the headline item.
# The shared-category picks are kept identical to panel c so the cast thread holds across b->c.
PREFER <- c(Diet_Weekly = "processed_meat_intake_f1349_0_0",
            Sleep = "sleep_duration_f1160_0_0")

# Held-out generalization accuracy with a CI computed FROM the held-out data.
# point = held-out accuracy averaged over the 3 repeat visits; bar = 95% bootstrap
# CI (resampling held-out PEOPLE) -- the uncertainty of the held-out point itself,
# computed identically for every exposure (so even area-level exposures get a real
# interval) and large enough to see (~+/-0.025) because the held-out set is ~3.4k
# people, not the 50k training folds. Pre-computed by support/module6_holdout_ci.R.
ci <- fread(file.path(heap_project_output("module6_pes_longitudinal","base"),
                      "PESlong_base_HoldoutAccuracyCI.tsv"))
agg <- ci[, .(exposure_id, category, exposure_type,
              r2  = ifelse(exposure_type == "continuous", point, NA_real_),
              r2_lo = ifelse(exposure_type == "continuous", ci_lo, NA_real_),
              r2_hi = ifelse(exposure_type == "continuous", ci_hi, NA_real_),
              auc = ifelse(exposure_type == "binary", point, NA_real_),
              auc_lo = ifelse(exposure_type == "binary", ci_lo, NA_real_),
              auc_hi = ifelse(exposure_type == "binary", ci_hi, NA_real_))]
agg[, is_cast := exposure_id %in% CAST]

# SHARED category ordering so the two sub-panels keep the SAME relative order
# (binary is a subset of categories but must not reshuffle): rank by continuous
# held-out R2 median; any binary-only category is appended at the bottom.
.ordc <- agg[exposure_type == "continuous" & is.finite(r2), .(m = median(r2)), by = category][order(m)]$category
.ordb <- agg[exposure_type == "binary" & is.finite(auc), .(m = median(auc)), by = category][order(m)]$category
ORD_GLOBAL <- c(setdiff(.ordb, .ordc), .ordc)

build_sub <- function(typ, mcol, locol, hicol, xlab, ttl, xlim, ref) {
  d <- agg[exposure_type == typ]
  d[, `:=`(val = get(mcol), lo = get(locol), hi = get(hicol))]; d <- d[is.finite(val)]
  present <- ORD_GLOBAL[ORD_GLOBAL %in% unique(as.character(d$category))]   # shared order, present-only
  d[, category := factor(category, levels = present)]
  d[, ycat := as.integer(category)]
  set.seed(1); d[, yj := ycat + runif(.N, -0.17, 0.17)]      # precompute jitter so bar tracks its dot
  d[, prefid := ifelse(as.character(category) %in% names(PREFER), PREFER[as.character(category)], NA_character_)]
  d[, pref := !is.na(prefid) & exposure_id == prefid]
  ex <- d[!exposure_id %in% EX_EXCLUDE][order(category, -pref, -is_cast, -val)][, .SD[1], by = category]
  ex[, lab := pretty_ex(exposure_id)]
  ex[, yj := ycat]    # centre the highlighted exemplar on its row
  gp <- 0.045 * (xlim[2] - xlim[1])           # fixed gap: label sits this far RIGHT of the dot, on-row
  ex[, `:=`(lx = val + gp, ly = ycat)]
  ggplot(d, aes(val, yj, colour = category)) +
    geom_vline(xintercept = ref, colour = "grey75", linewidth = 0.4, linetype = "22") +   # no-skill
    geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, linewidth = 0.35, alpha = 0.30) +
    geom_point(size = PSZ_S, alpha = if (CELL) 0.55 else 0.85) +
    geom_errorbarh(data = ex, aes(xmin = lo, xmax = hi), height = 0, linewidth = 0.6, alpha = 0.9) +
    geom_point(data = ex, size = PSZ_L) +
    geom_point(data = ex, size = PSZ_L, shape = 1, colour = "grey20", stroke = if (CELL) 0.4 else 0.5) +  # dark outline ring -> the exemplar reads bolder than the cloud
    {if (CELL) list(
        geom_segment(data = ex, aes(x = val, xend = lx, y = ycat, yend = ly), colour = "grey55", linewidth = 0.25, inherit.aes = FALSE),
        geom_text(data = ex, aes(x = lx, y = ly, label = lab), hjust = 0, colour = "#222222", size = LBL, inherit.aes = FALSE))
     else ggrepel::geom_text_repel(data = ex, aes(label = lab), colour = "#222222", size = LBL,
        box.padding = 0.45, point.padding = 0.2, min.segment.length = 0, segment.colour = "grey65", max.overlaps = 30, seed = 1)} +
    scale_colour_exposure(drop = TRUE, guide = "none") +
    scale_y_continuous(breaks = seq_along(present), labels = pcat(present), expand = expansion(add = if (CELL) 0.55 else 0.6)) +
    scale_x_continuous(limits = xlim, expand = expansion(mult = c(0.02, 0.04))) +
    labs(title = ttl, x = xlab, y = NULL) +
    theme_heap(base_size = BS) +
    theme(panel.grid.major.y = if (CELL) element_blank() else element_line(colour = "grey93"),
          panel.grid.major.x = if (CELL) element_blank() else element_line(colour = "grey92"),
          axis.text.y = element_text(size = rel(0.85)),
          axis.title = element_text(face = if (CELL) "plain" else "bold"),
          plot.title = element_text(size = if (CELL) TTL-1 else rel(1.05), hjust = 0.5),
          plot.title.position = "panel", plot.margin = margin(4,3,2,3))
}
b1 <- build_sub("continuous", "r2", "r2_lo", "r2_hi", expression("proteome-only held-out "*R^2), "Continuous exposures", c(-0.08, 0.70), ref = 0)
b2 <- build_sub("binary", "auc", "auc_lo", "auc_hi", "proteome-only held-out AUC", "Binary exposures", c(0.43, 1.33), ref = 0.5)

p <- b1 | b2
FIGDIR <- file.path(Sys.getenv("HEAP_PROJECT_ROOT", "/n/groups/patel/IGLOO/UKB/HEAP"), "figures", "exploratory/module6")
if (CELL) {
  # NB: overall title removed -> lives in the separate figure legend file.
  ggsave(file.path(FIGDIR,"fig_m6_panel_b_cell.png"), p, width = B_W, height = B_H, dpi = 400, bg = "white")
  ggsave(file.path(FIGDIR,"fig_m6_panel_b_cell.pdf"), p, width = B_W, height = B_H, bg = "white")
  message("panel b CELL done")
} else {
  p <- p + plot_annotation(title = "The proteome captures some exposures far better than others",
          subtitle = "each dot = one exposure (proteome-only score); point = held-out accuracy (mean of 3 visits), bar = 95% bootstrap CI; dashed line = no skill; one labeled per category",
          theme = theme(plot.title = element_text(face = "bold", size = 13),
                        plot.subtitle = element_text(size = 9.5, colour = "grey35")))
  ggsave(file.path(FIGDIR,"fig_m6_panel_b.png"), p, width = 12, height = 5.6, dpi = 175, bg = "white")
  message("panel b (split) done")
}
