#!/usr/bin/env Rscript

# ============================================================
# Cox dumbbell / step plots (Type5 only)
# 1) Per-exposure: top 10 diseases (ranked by Δ M0->M1 and Δ M2->M3)  [unchanged]
# 2) Disease-focused (VERSION B): for each disease code (E11/J43/F10),
#    show 3-stage step plot per exposure: M0 (covariates) → M2 (+exposure) → M3 (+PES)
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
topN_exposures_for_disease <- 10
label_top_k_disease <- 10

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
  "types_of_physical_activity_in_last_4_weeks_f6164_0_0.multi_Strenuous_sports", "Strenuous Sports",
  "types_of_physical_activity_in_last_4_weeks_f6164_0_0.multi_Other_exercises_.eg._swimming._cycling._keep_fit._bowling.", "Swimming/Cycling/etc.",
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
# Disease-focused plot (VERSION B):
# For a disease code, show step plot per exposure: M0 → M2 → M3
# - exposures ranked by Δ(M2→M3) then Δ(M0→M1)
# - open circles for M0/M2, filled for M3
# - labels show Δ(M2→M3) for top-k exposures
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
    filter(is.finite(cindex_M0), is.finite(cindex_M2), is.finite(cindex_M3)) %>%
    filter(is.finite(delta_cindex_M2_to_M3) | is.finite(delta_cindex_M0_to_M1))
  
  if (nrow(dd0) == 0) {
    message("No rows found for disease ", dc, " at Type5 after filtering.")
    return(NULL)
  }
  
  dd <- dd0 %>%
    group_by(exposure_label) %>%
    arrange(desc(delta_cindex_M2_to_M3), desc(delta_cindex_M0_to_M1)) %>%
    slice(1) %>%
    ungroup() %>%
    arrange(desc(delta_cindex_M2_to_M3), desc(delta_cindex_M0_to_M1)) %>%
    slice_head(n = topN_exposures) %>%
    mutate(exposure_label = factor(exposure_label, levels = rev(unique(exposure_label))))
  
  # --- SHORT LABELS for compact legend ---
  stage_levels <- c("M1 (covariates)", "M2 (+exposure)", "M3 (+PES)")
  
  pts <- bind_rows(
    dd %>% transmute(exposure_label, model = stage_levels[1], cindex = cindex_M0),
    dd %>% transmute(exposure_label, model = stage_levels[2], cindex = cindex_M2),
    dd %>% transmute(exposure_label, model = stage_levels[3], cindex = cindex_M3)
  ) %>%
    filter(is.finite(cindex)) %>%
    mutate(model = factor(model, levels = stage_levels))
  
  seg <- bind_rows(
    dd %>% transmute(exposure_label, x = cindex_M0, xend = cindex_M2, inc = "M1 → M2"),
    dd %>% transmute(exposure_label, x = cindex_M2, xend = cindex_M3, inc = "M2 → M3")
  ) %>%
    mutate(inc = factor(inc, levels = c("M1 → M2", "M2 → M3")))
  
  label_df <- dd %>%
    arrange(desc(delta_cindex_M2_to_M3), desc(delta_cindex_M0_to_M1)) %>%
    slice_head(n = label_top_k) %>%
    transmute(
      exposure_label,
      x = cindex_M3,
      label = sprintf("Δ=%.3f", delta_cindex_M2_to_M3)
    )
  
  vals <- c(dd$cindex_M0, dd$cindex_M2, dd$cindex_M3)
  xmin <- min(vals, na.rm = TRUE)
  xmax <- max(vals, na.rm = TRUE)
  span <- xmax - xmin
  if (!is.finite(span) || span <= 0) span <- 0.02
  
  # compact x padding; labels close to M3
  xlim <- c(xmin - 0.02 * span, xmax + 0.10 * span)
  label_offset <- 0.015 * span
  
  # increment colors
  inc_cols <- c(
    "M1 → M2" = "#3B7DDD",
    "M2 → M3" = "#D55E00"
  )
  
  p <- ggplot() +
    geom_segment(
      data = seg,
      aes(x = x, xend = xend, y = exposure_label, yend = exposure_label, color = inc),
      linewidth = 0.9, alpha = 0.92, lineend = "round"
    ) +
    geom_point(
      data = pts,
      aes(x = cindex, y = exposure_label, shape = model),
      size = 2.9, stroke = 0.9, color = "black"
    ) +
    geom_text(
      data = label_df,
      aes(x = x + label_offset, y = exposure_label, label = label),
      size = 3.0, hjust = 0, color = "grey10"
    ) +
    coord_cartesian(xlim = xlim, clip = "on") +
    scale_color_manual(values = inc_cols) +
    scale_shape_manual(values = c(
      "M1 (covariates)" = 16,
      "M2 (+exposure)"  = 15,
      "M3 (+PES)"       = 17
    )) +
    labs(
      title = paste0("Disease ", dc),
      subtitle = paste0("Top ", nrow(dd), " exposures by Δ(M2 → M3)."),
      x = "c-index",
      y = NULL,
      color = "Increment",
      shape = "Model"
    ) +
    theme_classic(base_size = 11) +
    theme(
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.y  = element_text(size = 9.5, color = "grey20"),
      axis.text.x  = element_text(size = 10, color = "grey20"),
      axis.title.x = element_text(size = 11, margin = margin(t = 6)),
      
      plot.title    = element_text(size = 16, face = "bold", margin = margin(b = 4)),
      plot.subtitle = element_text(size = 11, color = "grey30", margin = margin(b = 6)),
      
      panel.grid.major.x = element_line(color = "grey92", linewidth = 0.45),
      panel.grid.major.y = element_line(color = "grey94", linewidth = 0.40),
      panel.grid.minor   = element_blank(),
      
      # --- KEY FIX: prevent legend clipping ---
      legend.position = "bottom",
      legend.box = "vertical",
      legend.box.just = "center",
      legend.title = element_text(size = 10.5),
      legend.text  = element_text(size = 10),
      
      legend.key.width = grid::unit(1.0, "lines"),
      legend.key.height = grid::unit(0.85, "lines"),
      legend.spacing.y = grid::unit(0.2, "lines"),
      legend.spacing.x = grid::unit(0.6, "lines"),
      
      # more bottom margin so legend isn't cut off
      plot.margin = margin(5, 5, 18, 5)
    ) +
    guides(
      color = guide_legend(order = 1, nrow = 1, byrow = TRUE, override.aes = list(linewidth = 2.5)),
      shape = guide_legend(order = 2, nrow = 1, byrow = TRUE, override.aes = list(size = 3.2))
    )
  
  # slightly taller save to accommodate legends cleanly
  height_in <- max(3.4, 0.22 * nrow(dd) + 2.2)
  ggsave(out_file, p, width = 6.5, height = height_in, dpi = 1000, limitsize = FALSE)
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
      out_file = file.path(out_dir, paste0(
        "Disease_", dc, "_AcrossExposures_STEP_M0_M2_M3_",
        type_for_plots, "_PESprot_Top", topN_exposures_for_disease, ".png"
      )),
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
      out_file = file.path(out_dir, paste0(
        "Disease_", dc, "_AcrossExposures_STEP_M0_M2_M3_",
        type_for_plots, "_PESfull_Top", topN_exposures_for_disease, ".png"
      )),
      label_top_k = label_top_k_disease,
      topN_exposures = topN_exposures_for_disease
    )
  }
}


# Save Type5 table for convenience
fwrite(as.data.table(cox_all_T5), file.path(out_dir, paste0("cox_all_", type_for_plots, ".tsv")), sep = "\t")

message("\nDONE. Output in: ", out_dir, "\n")