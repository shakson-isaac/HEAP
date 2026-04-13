#!/usr/bin/env Rscript

# ============================================================
# Cox Top-10 dumbbell plots (ONE exposure per file, Type5 only)
# + One-disease (E11) across exposures
#
# Key fixes:
#   1) Canonical row per exposure–disease chosen ONCE (no per-column max()).
#   2) Robust per-panel x-limits with padding so points/labels don't go off-grid.
#   3) Filename sanitization only on basename (keeps paths sane).
#   4) FIX: prevent duplicated factor levels in disease-across-exposures plot
#      (dedupe exposure_label and use unique() for levels)
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(forcats)
  library(ggplot2)
  library(tibble)
})

# ----------------------------
# USER CONFIG
# ----------------------------
base_dir <- "/n/groups/patel/shakson_ukb/UK_Biobank/Data/Parallel/PES_test"
types <- paste0("Type", 1:5)

out_path <- "/n/groups/patel/shakson_ukb/UK_Biobank/Output/PES/"
out_dir  <- file.path(out_path)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

if (!dir.exists(out_dir)) stop("Output directory could not be created: ", out_dir)
testfile <- file.path(out_dir, paste0(".__write_test__", Sys.getpid()))
ok <- tryCatch({ writeLines("test", testfile); TRUE }, error = function(e) FALSE)
if (!ok) stop("No write permission in output directory: ", out_dir)
unlink(testfile)

pes_kind_main <- "PESprot"
pes_kind_supp <- "PESfull"

type_for_plots    <- "Type5"
topN_per_exposure <- 10
min_events        <- 50

use_disease_code_only <- TRUE
disease_of_interest   <- "E11"

# Optional label maps
exposure_label_map <- tibble::tribble(
  ~exposure_id, ~exposure_label,
  "alcohol_intake_frequency_f1558_0_0", "Alcohol frequency",
  "smoking_status_f20116_0_0_Current", "Smoking (current)",
  "fresh_fruit_intake_f1309_0_0", "Fresh fruit",
  "processed_meat_intake_f1349_0_0", "Processed meat",
  "summed_met_minutes_per_week_for_all_activity_f22040_0_0", "Physical activity (MET-min/wk)",
  "types_of_physical_activity_in_last_4_weeks_f6164_0_0.multi_Strenuous_sports", "Strenuous sport (type)",
  "time_spent_watching_television_tv_f1070_0_0", "TV time",
  "types_of_physical_activity_in_last_4_weeks_f6164_0_0.multi_Other_exercises_.eg._swimming._cycling._keep_fit._bowling.", "Swimming/Cycling/etc."
)

disease_label_map <- tibble::tribble(
  ~disease_age_col, ~disease_label,
  "age_e11_first_reported_non_insulin_dependent_diabetes_mellitus_f130708_0_0", "T2D",
  "age_j43_first_reported_emphysema_f131490_0_0", "Emphysema",
  "age_j44_first_reported_other_chronic_obstructive_pulmonary_disease_f131492_0_0", "COPD",
  "age_n18_first_reported_chronic_renal_failure_f132032_0_0", "Chronic renal failure",
  "age_i10_first_reported_essential_primary_hypertension_f131286_0_0", "Hypertension"
)

# ----------------------------
# Helpers
# ----------------------------
safe_fread <- function(path) tryCatch(fread(path), error = function(e) NULL)

read_type_files <- function(base_dir, types, pattern, add_type=TRUE, add_from_filename=NULL) {
  out <- list()
  for (ty in types) {
    ty_dir <- file.path(base_dir, ty)
    if (!dir.exists(ty_dir)) next
    files <- list.files(ty_dir, full.names = TRUE)
    files <- files[str_detect(basename(files), pattern)]
    if (length(files) == 0) next
    
    dt <- rbindlist(lapply(files, function(f) {
      x <- safe_fread(f)
      if (is.null(x)) return(NULL)
      x[, file := basename(f)]
      x[, path := f]
      x
    }), fill = TRUE)
    
    if (add_type) dt[, Type := ty]
    if (!is.null(add_from_filename) && is.function(add_from_filename)) {
      dt <- as.data.table(add_from_filename(dt))
    }
    out[[ty]] <- dt
  }
  rbindlist(out, fill = TRUE)
}

filter_pes <- function(df, pes_kind) {
  if (!"pes_kind" %in% names(df)) return(df)
  df %>% filter(is.na(pes_kind) | pes_kind == !!pes_kind)
}

extract_disease_code <- function(disease_age_col) {
  x <- tolower(as.character(disease_age_col))
  code <- str_match(x, "^age_([^_]+)_first_reported_")[,2]
  code <- ifelse(is.na(code), str_match(x, "^age_([^_]+)_")[,2], code)
  toupper(code)
}

safe_basename <- function(x) {
  x <- str_replace_all(x, "[^A-Za-z0-9_\\-\\.]", "_")
  x <- str_replace_all(x, "_+", "_")
  x <- str_replace(x, "^_+", "")
  substr(x, 1, 220)
}

theme_set(theme_bw(base_size = 12))

choose_canonical_rows <- function(df) {
  df %>%
    mutate(
      cox_status = as.character(cox_status),
      n = suppressWarnings(as.numeric(n)),
      events = suppressWarnings(as.numeric(events))
    ) %>%
    group_by(pes_kind, Type, exposure_id, disease_age_col) %>%
    arrange(
      desc(cox_status %in% c("OK","ok","Success","SUCCESS")),
      desc(is.finite(events)), desc(events),
      desc(is.finite(n)), desc(n),
      .by_group = TRUE
    ) %>%
    slice(1) %>%
    ungroup()
}

panel_xlim <- function(v, right_pad_frac = 0.20, left_pad_frac = 0.03) {
  v <- v[is.finite(v)]
  if (length(v) == 0) return(c(0.5, 1.0))
  xmin <- min(v); xmax <- max(v)
  span <- xmax - xmin
  if (!is.finite(span) || span <= 0) span <- 0.01
  c(xmin - left_pad_frac*span, xmax + right_pad_frac*span)
}

# ----------------------------
# Load Cox outputs
# ----------------------------
cox_dt <- read_type_files(
  base_dir, types,
  "^Cox4All_.*__PES(prot|full)\\.tsv$",
  add_from_filename = function(dt) {
    dt <- as.data.table(dt)
    dt[, pes_kind := fifelse(
      str_detect(file, "__PESprot\\.tsv$"), "PESprot",
      fifelse(str_detect(file, "__PESfull\\.tsv$"), "PESfull", NA_character_)
    )]
    dt
  }
)

cox_all <- as_tibble(cox_dt) %>%
  mutate(Type = factor(Type, levels = types))

if (nrow(cox_all) == 0) stop("No Cox4All_*__PESprot/full.tsv found under: ", base_dir)

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

cox_T5 <- cox_all %>% filter(Type == type_for_plots)
cox_T5 <- choose_canonical_rows(cox_T5)

# ----------------------------
# Pick top N diseases per exposure at Type5
# ----------------------------
pick_top_diseases_T5 <- function(df_T5, exposure_label_one,
                                 delta_col = c("delta_cindex_M0_to_M1", "delta_cindex_M2_to_M3"),
                                 topN = 10,
                                 min_events = 0) {
  delta_col <- match.arg(delta_col)
  
  df_T5 %>%
    filter(exposure_label == exposure_label_one) %>%
    filter(is.na(events) | events >= min_events) %>%
    mutate(delta = suppressWarnings(as.numeric(.data[[delta_col]]))) %>%
    filter(is.finite(delta)) %>%
    arrange(desc(delta)) %>%
    slice_head(n = topN) %>%
    pull(disease_label_short)
}

# ----------------------------
# Plot one exposure (Type5 only)
# ----------------------------
plot_one_exposure_T5 <- function(df_T5, exposure_label_one, diseases_keep, rank_tag, out_file,
                                 label_top_k = 3) {
  
  dd <- df_T5 %>%
    filter(exposure_label == exposure_label_one,
           disease_label_short %in% diseases_keep) %>%
    transmute(
      disease_label_short,
      disease_label,
      cindex_M0, cindex_M1, cindex_M2, cindex_M3,
      delta01 = delta_cindex_M0_to_M1,
      delta23 = delta_cindex_M2_to_M3
    )
  
  if (nrow(dd) == 0) return(NULL)
  
  dd <- dd %>%
    mutate(
      ylab = ifelse(is.na(disease_label) | disease_label == disease_label_short,
                    disease_label_short,
                    paste0(disease_label, " (", disease_label_short, ")"))
    )
  
  y_levels <- dd %>%
    distinct(disease_label_short, ylab) %>%
    mutate(ord = match(disease_label_short, diseases_keep)) %>%
    arrange(ord) %>%
    pull(ylab)
  
  dd <- dd %>% mutate(ylab = factor(ylab, levels = rev(y_levels)))
  
  long <- bind_rows(
    dd %>% transmute(panel = "Add PES beyond exposure model (M2 → M3)",
                     ylab, base = cindex_M2, plus = cindex_M3, delta = delta23),
    dd %>% transmute(panel = "Add PES to covariates (M0 → M1)",
                     ylab, base = cindex_M0, plus = cindex_M1, delta = delta01)
  ) %>% filter(is.finite(base), is.finite(plus), is.finite(delta))
  
  label_df <- long %>%
    group_by(panel) %>%
    arrange(desc(delta), .by_group = TRUE) %>%
    slice_head(n = label_top_k) %>%
    ungroup() %>%
    mutate(label = sprintf("Δ=%.3f", delta))
  
  xlims <- long %>%
    group_by(panel) %>%
    summarize(xmin = panel_xlim(c(base, plus))[1],
              xmax = panel_xlim(c(base, plus))[2],
              .groups="drop")
  
  long <- long %>% left_join(xlims, by="panel")
  label_df <- label_df %>% left_join(xlims, by="panel")
  
  p <- ggplot(long, aes(y = ylab)) +
    geom_segment(aes(x = base, xend = plus, yend = ylab), linewidth = 0.8, alpha = 0.75) +
    geom_point(aes(x = base), shape = 21, fill = "white", stroke = 1.1, size = 2.7) +
    geom_point(aes(x = plus), shape = 19, size = 2.7) +
    geom_text(
      data = label_df,
      aes(x = plus + 0.01*(xmax-xmin), label = label),
      size = 3.2, hjust = 0
    ) +
    facet_wrap(~panel, ncol = 1, scales = "free_x") +
    labs(
      title = paste0(exposure_label_one, " — top ", length(diseases_keep), " diseases (", rank_tag, ", Type5)"),
      subtitle = paste0("Open circle = baseline model; filled = +PES. Filter: events ≥ ", min_events, "."),
      x = "c-index",
      y = NULL
    ) +
    coord_cartesian(clip = "off") +
    theme(
      strip.background = element_rect(fill="grey95", color=NA),
      axis.text.y = element_text(size=10),
      plot.margin = margin(5.5, 70, 5.5, 5.5)
    )
  
  ggsave(out_file, p, width = 10.8, height = 6.8, dpi = 300)
  p
}

# ----------------------------
# Per exposure runner
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
      fname01 <- safe_basename(paste0("Top", topN_per_exposure, "_", tag, "_", type_for_plots, "_", exp_lab, "_rankM0toM1.png"))
      plot_one_exposure_T5(cox_df_T5, exp_lab, keep01, "ranked by Δ(M0→M1)", file.path(plot_dir, fname01))
    }
    
    keep23 <- pick_top_diseases_T5(
      cox_df_T5, exposure_label_one = exp_lab,
      delta_col = "delta_cindex_M2_to_M3",
      topN = topN_per_exposure,
      min_events = min_events
    )
    if (length(keep23) >= 2) {
      fname23 <- safe_basename(paste0("Top", topN_per_exposure, "_", tag, "_", type_for_plots, "_", exp_lab, "_rankM2toM3.png"))
      plot_one_exposure_T5(cox_df_T5, exp_lab, keep23, "ranked by Δ(M2→M3)", file.path(plot_dir, fname23))
    }
  }
  
  message("Wrote per-exposure Type5 plots to: ", plot_dir)
}

# ----------------------------
# One disease across exposures (E11)  [FIXED duplicates]
# ----------------------------
plot_one_disease_across_exposures_T5 <- function(df_T5, disease_code = "E11",
                                                 out_file,
                                                 label_top_k = 6) {
  
  dd <- df_T5 %>%
    filter(disease_label_short == disease_code) %>%
    filter(is.na(events) | events >= min_events) %>%
    transmute(
      exposure_label,
      exposure_id,
      cindex_M0, cindex_M1, cindex_M2, cindex_M3,
      delta01 = delta_cindex_M0_to_M1,
      delta23 = delta_cindex_M2_to_M3,
      events
    ) %>%
    filter(is.finite(delta01) | is.finite(delta23))
  
  if (nrow(dd) == 0) {
    message("No rows found for disease ", disease_code, " at Type5 after filtering.")
    return(NULL)
  }
  
  # If multiple exposure_ids map to the same exposure_label, keep the "best" one
  dd <- dd %>%
    arrange(desc(delta23), desc(delta01), desc(events)) %>%
    group_by(exposure_label) %>%
    slice(1) %>%                     # keep single row per label
    ungroup()
  
  # sort exposures by delta23 primarily
  dd <- dd %>%
    arrange(desc(delta23), desc(delta01)) %>%
    mutate(exposure_label = factor(exposure_label, levels = rev(unique(exposure_label))))
  
  long <- bind_rows(
    dd %>% transmute(panel="Add PES beyond exposure model (M2 → M3)",
                     exposure_label, base=cindex_M2, plus=cindex_M3, delta=delta23),
    dd %>% transmute(panel="Add PES to covariates (M0 → M1)",
                     exposure_label, base=cindex_M0, plus=cindex_M1, delta=delta01)
  ) %>% filter(is.finite(base), is.finite(plus), is.finite(delta))
  
  label_df <- long %>%
    group_by(panel) %>%
    arrange(desc(delta), .by_group = TRUE) %>%
    slice_head(n = label_top_k) %>%
    ungroup() %>%
    mutate(label = sprintf("Δ=%.3f", delta))
  
  xlims <- long %>%
    group_by(panel) %>%
    summarize(xmin = panel_xlim(c(base, plus))[1],
              xmax = panel_xlim(c(base, plus))[2],
              .groups="drop")
  
  long <- long %>% left_join(xlims, by="panel")
  label_df <- label_df %>% left_join(xlims, by="panel")
  
  p <- ggplot(long, aes(y = exposure_label)) +
    geom_segment(aes(x=base, xend=plus, yend=exposure_label), linewidth=0.8, alpha=0.75) +
    geom_point(aes(x=base), shape=21, fill="white", stroke=1.1, size=2.7) +
    geom_point(aes(x=plus), shape=19, size=2.7) +
    geom_text(
      data=label_df,
      aes(x = plus + 0.01*(xmax-xmin), label=label),
      size=3.2, hjust=0
    ) +
    facet_wrap(~panel, ncol=1, scales="free_x") +
    labs(
      title = paste0("Disease ", disease_code, " (Type5): discrimination gain from adding PES across exposures"),
      subtitle = paste0("Open circle = baseline; filled = +PES. Exposures sorted by Δ(M2→M3). Filter: events ≥ ", min_events, "."),
      x="c-index",
      y=NULL
    ) +
    coord_cartesian(clip="off") +
    theme(
      strip.background = element_rect(fill="grey95", color=NA),
      axis.text.y = element_text(size=9),
      plot.margin = margin(5.5, 80, 5.5, 5.5)
    )
  
  ggsave(out_file, p, width=11.2, height=max(6, 0.35*nrow(dd) + 3), dpi=300)
  p
}

# ----------------------------
# Execute for MAIN and SUPP
# ----------------------------
cox_main_T5 <- filter_pes(cox_T5, pes_kind_main)
cox_supp_T5 <- filter_pes(cox_T5, pes_kind_supp)

if (nrow(cox_main_T5) > 0) {
  plot_one_disease_across_exposures_T5(
    df_T5 = cox_main_T5,
    disease_code = disease_of_interest,
    out_file = file.path(out_dir, paste0("Disease_", disease_of_interest, "_AcrossExposures_", type_for_plots, "_", pes_kind_main, ".png")),
    label_top_k = 6
  )
  run_per_exposure(cox_main_T5, pes_kind_main, out_dir)
}

if (nrow(cox_supp_T5) > 0) {
  plot_one_disease_across_exposures_T5(
    df_T5 = cox_supp_T5,
    disease_code = disease_of_interest,
    out_file = file.path(out_dir, paste0("Disease_", disease_of_interest, "_AcrossExposures_", type_for_plots, "_", pes_kind_supp, ".png")),
    label_top_k = 6
  )
  run_per_exposure(cox_supp_T5, pes_kind_supp, out_dir)
}

# Save Type5 canonical table for debugging
fwrite(as.data.table(cox_T5), file.path(out_dir, paste0("cox_canonical_", type_for_plots, ".tsv")), sep="\t")

message("\nDONE. Output in: ", out_dir, "\n")

