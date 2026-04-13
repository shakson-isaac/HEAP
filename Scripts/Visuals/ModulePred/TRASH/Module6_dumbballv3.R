#!/usr/bin/env Rscript

# ============================================================
# Cox Top-10 dumbbell plots (ONE exposure per file, Type5 only)
# FIXES:
#   (1) Δ labels now ALWAYS match the plotted dumbbell segments:
#       - We select ONE consistent row per group (disease or exposure),
#         then recompute Δ = (plus - base).
#       - No more mixing max(cindex_*) from one row with max(delta_*) from another.
#   (2) Δ text placement fixed for tight x-ranges:
#       - Labels placed at midpoint of segment (scale-free, stays aligned).
#   (3) Axis clipping/out-of-grid fixed:
#       - coord_cartesian(xlim=..., clip="on") with padded limits from data.
#   (4) Filenames sanitized ONLY on basename (paths preserved).
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

# Fail loudly if out_dir isn't usable (prevents silent writing to getwd())
if (!dir.exists(out_dir)) stop("Output directory could not be created: ", out_dir)
testfile <- file.path(out_dir, paste0(".__write_test__", Sys.getpid()))
ok <- tryCatch({ writeLines("test", testfile); TRUE }, error = function(e) FALSE)
if (!ok) stop("No write permission in output directory: ", out_dir)
unlink(testfile)

pes_kind_main <- "PESprot"
pes_kind_supp <- "PESfull"

# Only use this Type
type_for_plots    <- "Type5"
topN_per_exposure <- 10
min_events        <- 50

# Label display
use_disease_code_only <- TRUE

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

# Keep only Type5 for plotting & ranking
cox_all_T5 <- cox_all %>% filter(Type == type_for_plots)

# ----------------------------
# Pick top N diseases per exposure at Type5 (robust de-dup)
#   - We rank using a consistent choice per disease: max(delta)
#     (this is only for selection; plotting will choose 1 consistent row)
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
# Plot one exposure (Type5 only)
#   - Select ONE consistent row per disease (ties broken by events)
#   - Recompute delta = plus - base so label matches dumbbell
# ----------------------------
plot_one_exposure_T5 <- function(df_T5, exposure_label_one, diseases_keep, rank_tag, out_file,
                                 label_top_k = 2,
                                 trim_x = TRUE,
                                 trim_quantile = 0.98,
                                 rank_mode = c("M2toM3", "M0toM1")) {
  
  rank_mode <- match.arg(rank_mode)
  
  dd0 <- df_T5 %>%
    filter(exposure_label == exposure_label_one,
           disease_label_short %in% diseases_keep) %>%
    mutate(
      delta01 = suppressWarnings(as.numeric(delta_cindex_M0_to_M1)),
      delta23 = suppressWarnings(as.numeric(delta_cindex_M2_to_M3)),
      cindex_M0 = suppressWarnings(as.numeric(cindex_M0)),
      cindex_M1 = suppressWarnings(as.numeric(cindex_M1)),
      cindex_M2 = suppressWarnings(as.numeric(cindex_M2)),
      cindex_M3 = suppressWarnings(as.numeric(cindex_M3))
    )
  
  if (nrow(dd0) == 0) return(NULL)
  
  # Pick ONE consistent row per disease (important: avoid mixing max across columns)
  if (rank_mode == "M2toM3") {
    dd <- dd0 %>%
      group_by(disease_label_short) %>%
      arrange(desc(delta23), desc(delta01), desc(events)) %>%
      slice(1) %>%
      ungroup()
  } else {
    dd <- dd0 %>%
      group_by(disease_label_short) %>%
      arrange(desc(delta01), desc(delta23), desc(events)) %>%
      slice(1) %>%
      ungroup()
  }
  
  # nice y label: "T2D (E11)" if disease_label exists and differs from code
  dd <- dd %>%
    mutate(
      ylab = ifelse(is.na(disease_label) | disease_label == disease_label_short,
                    disease_label_short,
                    paste0(disease_label, " (", disease_label_short, ")")),
      # recompute deltas to guarantee label matches plotted segment
      delta01 = cindex_M1 - cindex_M0,
      delta23 = cindex_M3 - cindex_M2
    )
  
  # Keep y order consistent with diseases_keep ranking (desc)
  # Map disease_label_short -> ylab order
  y_by_code <- dd %>% select(disease_label_short, ylab) %>% distinct()
  y_levels <- y_by_code$ylab[match(diseases_keep, y_by_code$disease_label_short)]
  y_levels <- y_levels[!is.na(y_levels)]
  dd <- dd %>% mutate(ylab = factor(ylab, levels = rev(unique(y_levels))))
  
  # two panels
  long <- bind_rows(
    dd %>% transmute(panel = "Add PES beyond exposure model (M2 → M3)",
                     ylab,
                     base = cindex_M2, plus = cindex_M3, delta = delta23),
    dd %>% transmute(panel = "Add PES to covariates (M0 → M1)",
                     ylab,
                     base = cindex_M0, plus = cindex_M1, delta = delta01)
  ) %>%
    filter(is.finite(base), is.finite(plus), is.finite(delta))
  
  # label only the top K deltas per panel
  label_df <- long %>%
    group_by(panel) %>%
    arrange(desc(delta), .by_group = TRUE) %>%
    slice_head(n = label_top_k) %>%
    ungroup() %>%
    mutate(label = sprintf("Δ=%.3f", delta))
  
  # Robust x-limits: always include all points + modest padding
  vals <- c(long$base, long$plus)
  xmin <- min(vals, na.rm=TRUE)
  xmax <- max(vals, na.rm=TRUE)
  if (isTRUE(trim_x)) {
    xmax_suggest <- as.numeric(quantile(vals, probs = trim_quantile, na.rm=TRUE))
    xmax <- max(xmax, xmax_suggest)
  }
  span <- xmax - xmin
  if (!is.finite(span) || span <= 0) span <- 0.05
  xlim_min <- xmin - 0.05 * span
  xlim_max <- xmax + 0.10 * span
  
  p <- ggplot(long, aes(y = ylab)) +
    geom_segment(aes(x = base, xend = plus, yend = ylab), linewidth = 0.8, alpha = 0.7) +
    geom_point(aes(x = base), shape = 21, fill = "white", stroke = 1.1, size = 2.7) +
    geom_point(aes(x = plus), shape = 19, size = 2.7) +
    # midpoint labels => always aligned & scale-free
    geom_text(data = label_df,
              aes(x = (base + plus)/2, label = label),
              size = 3.2, vjust = -0.8) +
    facet_wrap(~panel, ncol = 1, scales = "free_x") +
    labs(
      title = paste0(exposure_label_one, " — top ", length(diseases_keep), " diseases (", rank_tag, ", Type5)"),
      subtitle = paste0("Open circle = baseline model; filled = +PES. Filter: events ≥ ", min_events, "."),
      x = "c-index",
      y = NULL
    ) +
    coord_cartesian(xlim = c(xlim_min, xlim_max), clip = "on") +
    theme(
      strip.background = element_rect(fill="grey95", color=NA),
      axis.text.y = element_text(size=10),
      plot.margin = margin(5.5, 12, 5.5, 5.5)
    )
  
  ggsave(out_file, p, width = 10.5, height = 6.5, dpi = 300)
  p
}

# ----------------------------
# Run per exposure, for a given PES kind
# ----------------------------
run_per_exposure <- function(cox_df_T5, tag, out_dir) {
  if (nrow(cox_df_T5) == 0) return(invisible(NULL))
  
  plot_dir <- file.path(out_dir, paste0("CoxTop10_PerExposure_", tag, "_", type_for_plots))
  dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
  if (!dir.exists(plot_dir)) stop("Could not create plot_dir: ", plot_dir)
  
  exposures <- sort(unique(cox_df_T5$exposure_label))
  
  for (exp_lab in exposures) {
    
    # rank by Δ M0->M1
    keep01 <- pick_top_diseases_T5(
      cox_df_T5, exposure_label_one = exp_lab,
      delta_col = "delta_cindex_M0_to_M1",
      topN = topN_per_exposure,
      min_events = min_events
    )
    
    if (length(keep01) >= 2) {
      fname01 <- safe_basename(paste0(
        "Top", topN_per_exposure, "_", tag, "_", type_for_plots, "_",
        exp_lab, "_rankM0toM1.png"
      ))
      out01 <- file.path(plot_dir, fname01)
      plot_one_exposure_T5(cox_df_T5, exp_lab, keep01, "ranked by Δ(M0→M1)", out01,
                           rank_mode = "M0toM1")
    }
    
    # rank by Δ M2->M3
    keep23 <- pick_top_diseases_T5(
      cox_df_T5, exposure_label_one = exp_lab,
      delta_col = "delta_cindex_M2_to_M3",
      topN = topN_per_exposure,
      min_events = min_events
    )
    
    if (length(keep23) >= 2) {
      fname23 <- safe_basename(paste0(
        "Top", topN_per_exposure, "_", tag, "_", type_for_plots, "_",
        exp_lab, "_rankM2toM3.png"
      ))
      out23 <- file.path(plot_dir, fname23)
      plot_one_exposure_T5(cox_df_T5, exp_lab, keep23, "ranked by Δ(M2→M3)", out23,
                           rank_mode = "M2toM3")
    }
  }
  
  message("Wrote per-exposure Type5 plots to: ", plot_dir)
}

# ----------------------------
# Plot one disease (Type5) across exposures
#   - Select ONE consistent row per exposure
#   - Recompute delta = plus - base
#   - Midpoint delta labels
#   - Padded x-lims, clip on
# ----------------------------
plot_one_disease_across_exposures_T5 <- function(df_T5, disease_code = "E11",
                                                 out_file,
                                                 label_top_k = 6) {
  
  dd <- df_T5 %>%
    filter(disease_label_short == disease_code) %>%
    filter(is.na(events) | events >= min_events) %>%
    mutate(
      delta01 = suppressWarnings(as.numeric(delta_cindex_M0_to_M1)),
      delta23 = suppressWarnings(as.numeric(delta_cindex_M2_to_M3)),
      cindex_M0 = suppressWarnings(as.numeric(cindex_M0)),
      cindex_M1 = suppressWarnings(as.numeric(cindex_M1)),
      cindex_M2 = suppressWarnings(as.numeric(cindex_M2)),
      cindex_M3 = suppressWarnings(as.numeric(cindex_M3))
    ) %>%
    group_by(exposure_label) %>%
    arrange(desc(delta23), desc(delta01), desc(events)) %>%
    slice(1) %>%
    ungroup() %>%
    mutate(
      delta01 = cindex_M1 - cindex_M0,
      delta23 = cindex_M3 - cindex_M2
    )
  
  if (nrow(dd) == 0) {
    message("No rows found for disease ", disease_code, " at Type5 after filtering.")
    return(NULL)
  }
  
  # order exposures by Δ(M2->M3)
  dd <- dd %>%
    arrange(desc(delta23), desc(delta01)) %>%
    mutate(exposure_label = factor(exposure_label, levels = rev(exposure_label)))
  
  long <- bind_rows(
    dd %>% transmute(panel="Add PES beyond exposure model (M2 → M3)",
                     exposure_label,
                     base=cindex_M2, plus=cindex_M3, delta=delta23),
    dd %>% transmute(panel="Add PES to covariates (M0 → M1)",
                     exposure_label,
                     base=cindex_M0, plus=cindex_M1, delta=delta01)
  ) %>% filter(is.finite(base), is.finite(plus), is.finite(delta))
  
  label_df <- long %>%
    group_by(panel) %>%
    arrange(desc(delta), .by_group = TRUE) %>%
    slice_head(n = label_top_k) %>%
    ungroup() %>%
    mutate(label = sprintf("Δ=%.3f", delta))
  
  # padded x-lims so labels/points stay inside
  vals <- c(long$base, long$plus)
  xmin <- min(vals, na.rm=TRUE)
  xmax <- max(vals, na.rm=TRUE)
  span <- xmax - xmin
  if (!is.finite(span) || span <= 0) span <- 0.01
  xlim <- c(xmin - 0.05*span, xmax + 0.12*span)
  
  p <- ggplot(long, aes(y = exposure_label)) +
    geom_segment(aes(x=base, xend=plus, yend=exposure_label), linewidth=0.8, alpha=0.7) +
    geom_point(aes(x=base), shape=21, fill="white", stroke=1.1, size=2.7) +
    geom_point(aes(x=plus), shape=19, size=2.7) +
    geom_text(data=label_df,
              aes(x=(base + plus)/2, label=label),
              size=3.2, vjust=-0.8) +
    facet_wrap(~panel, ncol=1, scales="free_x") +
    labs(
      title = paste0("Disease ", disease_code, " (Type5): discrimination gain from adding PES across exposures"),
      subtitle = paste0("Open circle = baseline; filled = +PES. Exposures sorted by Δ(M2→M3). Filter: events ≥ ", min_events, "."),
      x="c-index",
      y=NULL
    ) +
    coord_cartesian(xlim = xlim, clip="on") +
    theme(
      strip.background = element_rect(fill="grey95", color=NA),
      axis.text.y = element_text(size=9),
      plot.margin = margin(5.5, 12, 5.5, 5.5)
    )
  
  ggsave(out_file, p, width=10.5, height=max(6, 0.35*nrow(dd) + 3), dpi=300)
  p
}

# ----------------------------
# Execute for MAIN and SUPP
# ----------------------------
cox_main_T5 <- filter_pes(cox_all_T5, pes_kind_main)
cox_supp_T5 <- filter_pes(cox_all_T5, pes_kind_supp)

if (nrow(cox_main_T5) > 0) {
  plot_one_disease_across_exposures_T5(
    df_T5 = cox_main_T5,
    disease_code = "E11",
    out_file = file.path(out_dir, "Disease_E11_AcrossExposures_Type5_PESprot.png"),
    label_top_k = 6
  )
}

if (nrow(cox_supp_T5) > 0) {
  plot_one_disease_across_exposures_T5(
    df_T5 = cox_supp_T5,
    disease_code = "E11",
    out_file = file.path(out_dir, "Disease_E11_AcrossExposures_Type5_PESfull.png"),
    label_top_k = 6
  )
}

if (nrow(cox_main_T5) > 0) run_per_exposure(cox_main_T5, pes_kind_main, out_dir)
if (nrow(cox_supp_T5) > 0) run_per_exposure(cox_supp_T5, pes_kind_supp, out_dir)

# Save Type5 table for convenience
fwrite(as.data.table(cox_all_T5), file.path(out_dir, paste0("cox_all_", type_for_plots, ".tsv")), sep="\t")

message("\nDONE. Output in: ", out_dir, "\n")
