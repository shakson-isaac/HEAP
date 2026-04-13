#!/usr/bin/env Rscript

# ============================================================
# Cox Top-10 dumbbell plots (ONE exposure per file, Type5 only)
# FIXED: filename sanitization now applies ONLY to the basename,
#        so paths are not flattened into giant filenames.
#
# For each exposure:
#   - Select top N diseases (ranked at Type5) by:
#       (1) Δc-index M0->M1 (cov -> cov+PES)     [rankM0toM1]
#       (2) Δc-index M2->M3 (cov+E -> cov+E+PES) [rankM2toM3]
#   - Make TWO figures per exposure:
#       * ranked by Δ(M0->M1)
#       * ranked by Δ(M2->M3)
#   - Each figure has two panels:
#       Panel A: M0 vs M1
#       Panel B: M2 vs M3
#
# Inputs expected (under base_dir/Type*/):
#   Cox4All_<Type>_<exposure>__PESprot.tsv
#   Cox4All_<Type>_<exposure>__PESfull.tsv
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
type_for_plots   <- "Type5"
topN_per_exposure <- 10
min_events <- 50

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
# ----------------------------
plot_one_exposure_T5 <- function(df_T5, exposure_label_one, diseases_keep, rank_tag, out_file,
                                 label_top_k = 2,
                                 trim_x = TRUE,
                                 trim_quantile = 0.98) {
  
  # Need disease_label too for nicer y-axis labels
  # Expect df_T5 has: disease_label_short (code), disease_label (maybe name), disease_code
  dd <- df_T5 %>%
    filter(exposure_label == exposure_label_one,
           disease_label_short %in% diseases_keep) %>%
    group_by(disease_label_short) %>%
    summarize(
      disease_label = dplyr::first(na.omit(disease_label)),
      cindex_M0 = max(cindex_M0, na.rm=TRUE),
      cindex_M1 = max(cindex_M1, na.rm=TRUE),
      cindex_M2 = max(cindex_M2, na.rm=TRUE),
      cindex_M3 = max(cindex_M3, na.rm=TRUE),
      delta01   = max(delta_cindex_M0_to_M1, na.rm=TRUE),
      delta23   = max(delta_cindex_M2_to_M3, na.rm=TRUE),
      .groups="drop"
    )
  
  if (nrow(dd) == 0) return(NULL)
  
  # nice y label: "T2D (E11)" if disease_label exists and differs from code
  dd <- dd %>%
    mutate(
      ylab = ifelse(is.na(disease_label) | disease_label == disease_label_short,
                    disease_label_short,
                    paste0(disease_label, " (", disease_label_short, ")"))
    )
  
  # order y by the ranking set (largest gain at top)
  # diseases_keep is already in ranked order (desc), so keep that order
  dd <- dd %>% mutate(ylab = factor(ylab, levels = rev(unique(ylab[match(diseases_keep, disease_label_short)]))))
  
  # two panels: M2->M3 and M0->M1
  long <- bind_rows(
    dd %>% transmute(panel = "Add PES beyond exposure model (M2 → M3)",
                     ylab,
                     base = cindex_M2, plus = cindex_M3, delta = delta23),
    dd %>% transmute(panel = "Add PES to covariates (M0 → M1)",
                     ylab,
                     base = cindex_M0, plus = cindex_M1, delta = delta01)
  ) %>%
    filter(is.finite(base), is.finite(plus), is.finite(delta))
  
  # optionally label only the top K rows per panel
  label_df <- long %>%
    group_by(panel) %>%
    arrange(desc(delta), .by_group = TRUE) %>%
    slice_head(n = label_top_k) %>%
    ungroup() %>%
    mutate(label = sprintf("Δ=%.3f", delta))
  
  # choose x-limits (trim extreme) to avoid one outlier compressing the plot
  # --- Robust x-limits: always include all points + room for labels ---
  vals <- c(long$base, long$plus)
  
  # start from actual min/max (never exclude points)
  xmin <- min(vals, na.rm=TRUE)
  xmax <- max(vals, na.rm=TRUE)
  
  # optionally compress extreme *right* tail for readability, but never drop points:
  # Instead of trimming the axis, we only use trimming to define a "suggested" max,
  # then take max(actual, suggested) so nothing gets clipped.
  if (isTRUE(trim_x)) {
    xmax_suggest <- as.numeric(quantile(vals, probs = trim_quantile, na.rm=TRUE))
    xmax <- max(xmax, xmax_suggest)
  }
  
  # padding: small fraction of span
  span <- xmax - xmin
  if (!is.finite(span) || span <= 0) span <- 0.05
  
  pad_left  <- 0.03 * span
  pad_right <- 0.18 * span   # extra room for Δ labels to the right
  
  # if we label, ensure enough extra right padding for the label text
  # (if label_top_k==0, reduce padding)
  if (is.null(label_top_k) || label_top_k <= 0) pad_right <- 0.06 * span
  
  xlim_min <- xmin - pad_left
  xlim_max <- xmax + pad_right
  
  
  p <- ggplot(long, aes(y = ylab)) +
    # segment
    geom_segment(aes(x = base, xend = plus, yend = ylab), linewidth = 0.8, alpha = 0.7) +
    # baseline open
    geom_point(aes(x = base), shape = 21, fill = "white", stroke = 1.1, size = 2.7) +
    # +PES filled
    geom_point(aes(x = plus), shape = 19, size = 2.7) +
    # label only top K deltas
    geom_text(data = label_df,
              aes(x = pmax(base, plus) + 0.005, label = label),
              size = 3.2, hjust = 0) +
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
      plot_one_exposure_T5(cox_df_T5, exp_lab, keep01, "ranked by Δ(M0→M1)", out01)
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
      plot_one_exposure_T5(cox_df_T5, exp_lab, keep23, "ranked by Δ(M2→M3)", out23)
    }
  }
  
  message("Wrote per-exposure Type5 plots to: ", plot_dir)
}

plot_one_disease_across_exposures_T5 <- function(df_T5, disease_code = "E11",
                                                 out_file,
                                                 label_top_k = 5) {
  
  dd <- df_T5 %>%
    filter(disease_label_short == disease_code) %>%
    # de-dup any repeats
    group_by(exposure_label) %>%
    summarize(
      cindex_M0 = max(cindex_M0, na.rm=TRUE),
      cindex_M1 = max(cindex_M1, na.rm=TRUE),
      cindex_M2 = max(cindex_M2, na.rm=TRUE),
      cindex_M3 = max(cindex_M3, na.rm=TRUE),
      delta01   = max(delta_cindex_M0_to_M1, na.rm=TRUE),
      delta23   = max(delta_cindex_M2_to_M3, na.rm=TRUE),
      events    = max(events, na.rm=TRUE),
      .groups="drop"
    ) %>%
    filter(is.finite(delta01) | is.finite(delta23)) %>%
    filter(is.na(events) | events >= min_events)
  
  if (nrow(dd) == 0) {
    message("No rows found for disease ", disease_code, " at Type5 after filtering.")
    return(NULL)
  }
  
  # order exposures by Δ(M2->M3) primarily
  dd <- dd %>% arrange(desc(delta23), desc(delta01)) %>%
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
  
  p <- ggplot(long, aes(y = exposure_label)) +
    geom_segment(aes(x=base, xend=plus, yend=exposure_label), linewidth=0.8, alpha=0.7) +
    geom_point(aes(x=base), shape=21, fill="white", stroke=1.1, size=2.7) +
    geom_point(aes(x=plus), shape=19, size=2.7) +
    geom_text(data=label_df,
              aes(x=pmax(base, plus) + 0.005, label=label),
              size=3.2, hjust=0) +
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
      plot.margin = margin(5.5, 60, 5.5, 5.5)
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
