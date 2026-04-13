#!/usr/bin/env Rscript

# ============================================================
# Cox dumbbell plots (Type5 only)
# 1) Per-exposure: top 10 diseases (ranked by Δ M0->M1 and Δ M2->M3)
# 2) Disease-focused: for each disease code (e.g., E11/J43/F10),
#    show top N exposures ranked by Δ(M2->M3)
#
# FIXES INCLUDED:
# - Use Type5 only across the board
# - Disease-focused plots filter by stable disease_code (NOT label_short)
# - De-dup per exposure/disease deterministically (slice top by delta)
# - Factor levels made unique to avoid duplicated level errors
# - Robust x-limits so points/labels stay inside grid
# - limitsize=FALSE for tall plots
# - Filename sanitization applies only to basename
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(tibble)
})

# ----------------------------
# USER CONFIG
# ----------------------------
base_dir <- "/n/groups/patel/shakson_ukb/UK_Biobank/Data/Parallel/PES_test"
out_dir  <- "/n/groups/patel/shakson_ukb/UK_Biobank/Output/PES"

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
if (!dir.exists(out_dir)) stop("Output directory could not be created: ", out_dir)
testfile <- file.path(out_dir, paste0(".__write_test__", Sys.getpid()))
ok <- tryCatch({ writeLines("test", testfile); TRUE }, error = function(e) FALSE)
if (!ok) stop("No write permission in output directory: ", out_dir)
unlink(testfile)

pes_kind_main <- "PESprot"
pes_kind_supp <- "PESfull"

type_for_plots <- "Type5"
topN_per_exposure <- 10
min_events <- 50

# disease-focused plots:
disease_codes_focus <- c("E11", "J43", "F10")  # add/remove here
topN_exposures_for_disease <- 6
label_top_k_disease <- 6

# Label display:
use_disease_code_only <- TRUE

# Optional label maps (safe to leave as-is)
exposure_label_map <- tibble::tribble(
  ~exposure_id, ~exposure_label,
  "alcohol_intake_frequency_f1558_0_0", "Alcohol frequency",
  "alcohol_drinker_status_f20117_0_0_Current", "Alcohol (Current)",
  "beef_intake_f1369_0_0","Beef Intake",
  "bread_intake_f1438_0_0","Bread Intake",
  "fresh_fruit_intake_f1309_0_0", "Fresh fruit",
  "income_score_england_f26411_0_0", "Income Score (SES)",
  "met_minutes_per_week_for_vigorous_activity_f22039_0_0","Vigorous Activity (MET-min/wk)",
  "summed_met_minutes_per_week_for_all_activity_f22040_0_0", "Physical activity (MET-min/wk)",
  "pork_intake_f1389_0_0", "Pork Intake",
  "poultry_intake_f1359_0_0", "Poultry Intake",
  "processed_meat_intake_f1349_0_0", "Processed meat",
  "smoking_status_f20116_0_0_Current", "Smoking (Current)",
  "summed_days_activity_f22033_0_0","Summed Days Activity",
  "summed_minutes_activity_f22034_0_0", "Summed Minutes Activity",
  "time_spent_watching_television_tv_f1070_0_0", "TV time",
  "types_of_physical_activity_in_last_4_weeks_f6164_0_0.multi_Strenuous_sports", "Exercise: Strenuous Sports",
  "types_of_physical_activity_in_last_4_weeks_f6164_0_0.multi_Other_exercises_.eg._swimming._cycling._keep_fit._bowling.", "Exercise: Swimming/Cycling/etc.",
  "usual_walking_pace_f924_0_0", "Usual Walking Pace",
  "water_intake_f1528_0_0", "Water Intake"
)

disease_label_map <- tibble::tribble(
  ~disease_age_col, ~disease_label,
  "age_e11_first_reported_non_insulin_dependent_diabetes_mellitus_f130708_0_0", "T2D",
  "age_j43_first_reported_emphysema_f131490_0_0", "Emphysema",
  "age_j44_first_reported_other_chronic_obstructive_pulmonary_disease_f131492_0_0", "COPD",
  "age_f10_first_reported_mental_and_behavioural_disorders_due_to_use_of_alcohol_f130560_0_0", "Alcohol use disorder",
  "age_n18_first_reported_chronic_renal_failure_f132032_0_0", "Chronic renal failure",
  "age_i10_first_reported_essential_primary_hypertension_f131286_0_0", "Hypertension"
)

# ----------------------------
# Helpers
# ----------------------------
safe_fread <- function(path) tryCatch(fread(path), error = function(e) NULL)

safe_basename <- function(x) {
  x <- str_replace_all(x, "[^A-Za-z0-9_\\-\\.]", "_")
  x <- str_replace_all(x, "_+", "_")
  x <- str_replace(x, "^_+", "")
  substr(x, 1, 220)
}

extract_disease_code <- function(disease_age_col) {
  x <- tolower(as.character(disease_age_col))
  code <- str_match(x, "^age_([^_]+)_first_reported_")[,2]
  code <- ifelse(is.na(code), str_match(x, "^age_([^_]+)_")[,2], code)
  toupper(code)
}

read_type_files <- function(base_dir, type_for_plots, pattern) {
  ty_dir <- file.path(base_dir, type_for_plots)
  if (!dir.exists(ty_dir)) stop("Missing type directory: ", ty_dir)
  files <- list.files(ty_dir, full.names = TRUE)
  files <- files[str_detect(basename(files), pattern)]
  if (length(files) == 0) stop("No files matching pattern under: ", ty_dir)
  
  dt <- rbindlist(lapply(files, function(f) {
    x <- safe_fread(f)
    if (is.null(x)) return(NULL)
    x[, file := basename(f)]
    x[, path := f]
    x
  }), fill = TRUE)
  
  dt[, Type := type_for_plots]
  dt
}

filter_pes <- function(df, pes_kind_keep) {
  if (!"pes_kind" %in% names(df)) return(df)
  df %>% filter(is.na(pes_kind) | pes_kind == pes_kind_keep)
}

# ----------------------------
# Load Cox outputs (Type5 only)
# ----------------------------
cox_dt <- read_type_files(
  base_dir = base_dir,
  type_for_plots = type_for_plots,
  pattern = "^Cox4All_.*__PES(prot|full)\\.tsv$"
)

cox_dt[, pes_kind := fifelse(
  str_detect(file, "__PESprot\\.tsv$"), "PESprot",
  fifelse(str_detect(file, "__PESfull\\.tsv$"), "PESfull", NA_character_)
)]

cox_all <- as_tibble(cox_dt)

if (nrow(cox_all) == 0) stop("No Cox4All_*__PESprot/full.tsv found under: ", file.path(base_dir, type_for_plots))

cox_all <- cox_all %>%
  mutate(
    exposure_id     = as.character(exposure_id),
    disease_age_col = as.character(disease_age_col),
    disease_code    = extract_disease_code(disease_age_col)
  ) %>%
  left_join(exposure_label_map, by="exposure_id") %>%
  left_join(disease_label_map,  by="disease_age_col") %>%
  mutate(
    exposure_label = ifelse(is.na(exposure_label), exposure_id, exposure_label),
    disease_label  = ifelse(is.na(disease_label), disease_code, disease_label),
    disease_label_short = if (use_disease_code_only) disease_code else disease_label
  ) %>%
  mutate(
    events = suppressWarnings(as.numeric(events)),
    cindex_M0 = suppressWarnings(as.numeric(cindex_M0)),
    cindex_M1 = suppressWarnings(as.numeric(cindex_M1)),
    cindex_M2 = suppressWarnings(as.numeric(cindex_M2)),
    cindex_M3 = suppressWarnings(as.numeric(cindex_M3)),
    delta_cindex_M0_to_M1 = suppressWarnings(as.numeric(delta_cindex_M0_to_M1)),
    delta_cindex_M2_to_M3 = suppressWarnings(as.numeric(delta_cindex_M2_to_M3))
  )

cox_all_T5 <- cox_all %>% filter(Type == type_for_plots)

cox_main_T5 <- filter_pes(cox_all_T5, pes_kind_main)
cox_supp_T5 <- filter_pes(cox_all_T5, pes_kind_supp)

# ----------------------------
# Pick top N diseases per exposure (Type5)
# ----------------------------
pick_top_diseases_T5 <- function(df_T5, exposure_label_one,
                                 delta_col = c("delta_cindex_M0_to_M1", "delta_cindex_M2_to_M3"),
                                 topN = 10,
                                 min_events = 0) {
  delta_col <- match.arg(delta_col)
  
  dd <- df_T5 %>%
    filter(exposure_label == exposure_label_one) %>%
    filter(is.na(events) | events >= min_events) %>%
    mutate(delta = suppressWarnings(as.numeric(.data[[delta_col]]))) %>%
    filter(is.finite(delta)) %>%
    group_by(disease_label_short) %>%
    summarize(delta = max(delta, na.rm=TRUE), .groups="drop") %>%
    arrange(desc(delta)) %>%
    slice_head(n = topN)
  
  dd$disease_label_short
}

# ----------------------------
# Plot one exposure (two-panel dumbbell) for selected diseases
# ----------------------------
plot_one_exposure_T5 <- function(df_T5, exposure_label_one, diseases_keep, rank_tag, out_file,
                                 label_top_k = 3) {
  
  dd0 <- df_T5 %>%
    filter(exposure_label == exposure_label_one,
           disease_label_short %in% diseases_keep) %>%
    filter(is.na(events) | events >= min_events) %>%
    filter(is.finite(cindex_M0), is.finite(cindex_M1), is.finite(cindex_M2), is.finite(cindex_M3))
  
  if (nrow(dd0) == 0) return(NULL)
  
  # de-dup per disease deterministically (take row with largest delta23, then delta01)
  dd <- dd0 %>%
    group_by(disease_label_short) %>%
    arrange(desc(delta_cindex_M2_to_M3), desc(delta_cindex_M0_to_M1)) %>%
    slice(1) %>%
    ungroup() %>%
    mutate(
      ylab = ifelse(is.na(disease_label) | disease_label == disease_label_short,
                    disease_label_short,
                    paste0(disease_label, " (", disease_label_short, ")"))
    )
  
  # order y by input ranked disease codes
  ordered_ylabs <- dd$ylab[match(diseases_keep, dd$disease_label_short)]
  ordered_ylabs <- ordered_ylabs[!is.na(ordered_ylabs)]
  dd <- dd %>% mutate(ylab = factor(ylab, levels = rev(unique(ordered_ylabs))))
  
  long <- bind_rows(
    dd %>% transmute(panel = "Add PES beyond exposure model (M2 → M3)",
                     ylab,
                     base = cindex_M2, plus = cindex_M3, delta = delta_cindex_M2_to_M3),
    dd %>% transmute(panel = "Add PES to covariates (M0 → M1)",
                     ylab,
                     base = cindex_M0, plus = cindex_M1, delta = delta_cindex_M0_to_M1)
  ) %>% filter(is.finite(base), is.finite(plus), is.finite(delta))
  
  long <- long %>%
    mutate(panel = factor(
      panel,
      levels = c(
        "Add PES to covariates (M0 → M1)",              # TOP
        "Add PES beyond exposure model (M2 → M3)"       # BOTTOM
      )
    ))
  
  label_df <- long %>%
    group_by(panel) %>%
    arrange(desc(delta), .by_group = TRUE) %>%
    slice_head(n = label_top_k) %>%
    ungroup() %>%
    mutate(label = sprintf("Δ=%.3f", delta))
  
  vals <- c(long$base, long$plus)
  xmin <- min(vals, na.rm=TRUE)
  xmax <- max(vals, na.rm=TRUE)
  span <- xmax - xmin
  if (!is.finite(span) || span <= 0) span <- 0.02
  xlim <- c(xmin - 0.03*span, xmax + 0.20*span)
  
  p <- ggplot(long, aes(y = ylab)) +
    geom_segment(aes(x = base, xend = plus, yend = ylab), linewidth = 0.8, alpha = 0.7) +
    geom_point(aes(x = base), shape = 21, fill = "white", stroke = 1.1, size = 2.7) +
    geom_point(aes(x = plus), shape = 19, size = 2.7) +
    geom_text(data = label_df,
              aes(x = pmax(base, plus) + 0.01*span, label = label),
              size = 3.2, hjust = 0) +
    facet_wrap(~panel, ncol = 1, scales = "free_x") +
    coord_cartesian(xlim = xlim, clip = "on") +
    labs(
      title = paste0(exposure_label_one, " — top ", length(unique(dd$ylab)), " diseases (", rank_tag, ", Type5)"),
      subtitle = paste0("Open circle = baseline; filled = +PES. Filter: events ≥ ", min_events, "."),
      x = "c-index", y = NULL
    ) +
    theme_bw(base_size = 12) +
    theme(
      strip.background = element_rect(fill="grey95", color=NA),
      axis.text.y = element_text(size=10),
      plot.margin = margin(5.5, 14, 5.5, 5.5)
    )
  
  ggsave(out_file, p, width = 10.5, height = 6.5, dpi = 300, limitsize = FALSE)
  p
}

# ----------------------------
# Run per-exposure suite
# ----------------------------
run_per_exposure <- function(cox_df_T5, tag, out_dir) {
  if (nrow(cox_df_T5) == 0) return(invisible(NULL))
  
  plot_dir <- file.path(out_dir, paste0("CoxTop10_PerExposure_", tag, "_", type_for_plots))
  dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
  if (!dir.exists(plot_dir)) stop("Could not create plot_dir: ", plot_dir)
  
  exposures <- sort(unique(cox_df_T5$exposure_label))
  
  for (exp_lab in exposures) {
    
    keep01 <- pick_top_diseases_T5(
      cox_df_T5, exposure_label_one = exp_lab,
      delta_col = "delta_cindex_M0_to_M1",
      topN = topN_per_exposure,
      min_events = min_events
    )
    
    if (length(keep01) >= 2) {
      fname01 <- safe_basename(paste0("Top", topN_per_exposure, "_", tag, "_", type_for_plots, "_",
                                      exp_lab, "_rankM0toM1.png"))
      out01 <- file.path(plot_dir, fname01)
      plot_one_exposure_T5(cox_df_T5, exp_lab, keep01, "ranked by Δ(M0→M1)", out01)
    }
    
    keep23 <- pick_top_diseases_T5(
      cox_df_T5, exposure_label_one = exp_lab,
      delta_col = "delta_cindex_M2_to_M3",
      topN = topN_per_exposure,
      min_events = min_events
    )
    
    if (length(keep23) >= 2) {
      fname23 <- safe_basename(paste0("Top", topN_per_exposure, "_", tag, "_", type_for_plots, "_",
                                      exp_lab, "_rankM2toM3.png"))
      out23 <- file.path(plot_dir, fname23)
      plot_one_exposure_T5(cox_df_T5, exp_lab, keep23, "ranked by Δ(M2→M3)", out23)
    }
  }
  
  message("Wrote per-exposure Type5 plots to: ", plot_dir)
}

# ----------------------------
# Disease-focused plot: top N exposures for a disease code
# IMPORTANT: filters by disease_code (stable) so E11/J43/F10 are distinct.
# ----------------------------
plot_one_disease_across_exposures_T5 <- function(df_T5,
                                                 disease_code_focus = "E11",
                                                 out_file,
                                                 label_top_k = 6,
                                                 topN_exposures = 30) {
  
  dc <- toupper(disease_code_focus)
  
  dd0 <- df_T5 %>%
    filter(disease_code == dc) %>%
    filter(is.na(events) | events >= min_events) %>%
    filter(is.finite(cindex_M0), is.finite(cindex_M1), is.finite(cindex_M2), is.finite(cindex_M3)) %>%
    filter(is.finite(delta_cindex_M0_to_M1) | is.finite(delta_cindex_M2_to_M3))
  
  if (nrow(dd0) == 0) {
    message("No rows found for disease ", dc, " at Type5 after filtering.")
    return(NULL)
  }
  
  # de-dup per exposure deterministically (choose row with largest delta23 then delta01)
  dd <- dd0 %>%
    group_by(exposure_label) %>%
    arrange(desc(delta_cindex_M2_to_M3), desc(delta_cindex_M0_to_M1)) %>%
    slice(1) %>%
    ungroup()
  
  # keep top exposures
  dd <- dd %>%
    arrange(desc(delta_cindex_M2_to_M3), desc(delta_cindex_M0_to_M1)) %>%
    slice_head(n = topN_exposures)
  
  exp_levels <- unique(dd$exposure_label)
  dd <- dd %>% mutate(exposure_label = factor(exposure_label, levels = rev(exp_levels)))
  
  long <- bind_rows(
    dd %>% transmute(panel="Add PES beyond exposure model (M2 → M3)",
                     exposure_label,
                     base=cindex_M2, plus=cindex_M3, delta=delta_cindex_M2_to_M3),
    dd %>% transmute(panel="Add PES to covariates (M0 → M1)",
                     exposure_label,
                     base=cindex_M0, plus=cindex_M1, delta=delta_cindex_M0_to_M1)
  ) %>% filter(is.finite(base), is.finite(plus), is.finite(delta))
  
  long <- long %>%
    mutate(panel = factor(
      panel,
      levels = c(
        "Add PES to covariates (M0 → M1)",              # TOP
        "Add PES beyond exposure model (M2 → M3)"       # BOTTOM
      )
    ))
  
  label_df <- long %>%
    group_by(panel) %>%
    arrange(desc(delta), .by_group = TRUE) %>%
    slice_head(n = label_top_k) %>%
    ungroup() %>%
    mutate(label = sprintf("Δ=%.3f", delta))
  
  vals <- c(long$base, long$plus)
  xmin <- min(vals, na.rm=TRUE)
  xmax <- max(vals, na.rm=TRUE)
  span <- xmax - xmin
  if (!is.finite(span) || span <= 0) span <- 0.02
  xlim <- c(xmin - 0.03*span, xmax + 0.20*span)
  
  p <- ggplot(long, aes(y = exposure_label)) +
    geom_segment(aes(x=base, xend=plus, yend=exposure_label), linewidth=0.8, alpha=0.7) +
    geom_point(aes(x=base), shape=21, fill="white", stroke=1.1, size=2.7) +
    geom_point(aes(x=plus), shape=19, size=2.7) +
    geom_text(data=label_df,
              aes(x=pmax(base, plus) + 0.01*span, label=label),
              size=3.2, hjust=0) +
    facet_wrap(~panel, ncol=1, scales="free_x") +
    coord_cartesian(xlim=xlim, clip="on") +
    labs(
      title = paste0("Disease ", dc, " (Type5): discrimination gain from adding PES across exposures"),
      subtitle = paste0("Open circle = baseline; filled = +PES. Top ", nrow(dd),
                        " exposures by Δ(M2→M3). Filter: events ≥ ", min_events, "."),
      x="c-index", y=NULL
    ) +
    theme_bw(base_size = 12) +
    theme(
      strip.background = element_rect(fill="grey95", color=NA),
      axis.text.y = element_text(size=9),
      plot.margin = margin(5.5, 60, 5.5, 5.5)
    )
  
  height_in <- max(6, 0.28*nrow(dd) + 3)
  ggsave(out_file, p, width=10.5, height=height_in, dpi=300, limitsize=FALSE)
  p
}

# ----------------------------
# EXECUTE
# ----------------------------
# Disease-focused plots for each disease code (E11/J43/F10)
if (nrow(cox_main_T5) > 0) {
  for (dc in disease_codes_focus) {
    plot_one_disease_across_exposures_T5(
      df_T5 = cox_main_T5,
      disease_code_focus = dc,
      out_file = file.path(out_dir, paste0("Disease_", dc, "_AcrossExposures_", type_for_plots, "_PESprot_Top", topN_exposures_for_disease, ".png")),
      label_top_k = label_top_k_disease,
      topN_exposures = topN_exposures_for_disease
    )
  }
}

if (nrow(cox_supp_T5) > 0) {
  for (dc in disease_codes_focus) {
    plot_one_disease_across_exposures_T5(
      df_T5 = cox_supp_T5,
      disease_code_focus = dc,
      out_file = file.path(out_dir, paste0("Disease_", dc, "_AcrossExposures_", type_for_plots, "_PESfull_Top", topN_exposures_for_disease, ".png")),
      label_top_k = label_top_k_disease,
      topN_exposures = topN_exposures_for_disease
    )
  }
}

# Per-exposure top-10 disease plots
if (nrow(cox_main_T5) > 0) run_per_exposure(cox_main_T5, pes_kind_main, out_dir)
if (nrow(cox_supp_T5) > 0) run_per_exposure(cox_supp_T5, pes_kind_supp, out_dir)

# Save Type5 table for convenience
fwrite(as.data.table(cox_all_T5), file.path(out_dir, paste0("cox_all_", type_for_plots, ".tsv")), sep="\t")

message("\nDONE. Output in: ", out_dir, "\n")
