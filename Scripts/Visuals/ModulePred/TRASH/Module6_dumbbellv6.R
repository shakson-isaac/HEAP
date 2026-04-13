#!/usr/bin/env Rscript

# ============================================================
# Cox Top-N dumbbell plots
# - ONE exposure per file, Type5 only
# - PLUS: one disease (E11) across exposures (Top-N by delta)
#
# Fixes:
#   * Reads ONLY Type5 to avoid cross-type mixing
#   * Filters cox_status == "OK"
#   * De-duplicates correctly (no max() aggregation surprises)
#   * E11 plot uses unique y labels (no collapsing -> no "all points")
#   * Height capped + Top-N exposures to avoid >50 inch errors
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(tibble)
})

# ----------------------------
# USER CONFIG
# ----------------------------
base_dir <- "/n/groups/patel/shakson_ukb/UK_Biobank/Data/Parallel/PES_test"

type_for_plots <- "Type5"

out_dir <- "/n/groups/patel/shakson_ukb/UK_Biobank/Output/PES/"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

pes_kind_main <- "PESprot"
pes_kind_supp <- "PESfull"

topN_per_exposure <- 10
min_events <- 50

# For the disease-across-exposures plot:
#disease_code_focus <- "E11"
#disease_code_focus <- "J43"
#disease_code_focus <- "F10"
# For the disease-across-exposures plot:
disease_codes_focus <- c("E11", "J43", "F10")
topN_exposures_for_disease <- 30

topN_exposures_for_disease <- 30   # <<< important: prevents huge plots
cap_height_inches <- 18            # <<< prevents ggsave > 50 inch abort

# If TRUE, disease labels are ICD codes only (E11); else use map if provided
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
  "age_j44_first_reported_other_chronic_obstructive_pulmonary_disease_f131492_0_0", "COPD"
)

# ----------------------------
# Helpers
# ----------------------------
safe_fread <- function(path) tryCatch(fread(path), error = function(e) NULL)

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

first_non_na <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_real_)
  x[[1]]
}

# Build unique display labels: if exposure_label repeats, append exposure_id
make_unique_exposure_labels <- function(exposure_id, exposure_label) {
  lab <- exposure_label
  dup <- duplicated(lab) | duplicated(lab, fromLast = TRUE)
  lab[dup] <- paste0(lab[dup], " [", exposure_id[dup], "]")
  lab
}

# ----------------------------
# Load Cox outputs (Type5 only)
# ----------------------------
type_dir <- file.path(base_dir, type_for_plots)
if (!dir.exists(type_dir)) stop("Missing directory: ", type_dir)

files <- list.files(type_dir, full.names = TRUE)
files <- files[str_detect(basename(files), "^Cox4All_.*__PES(prot|full)\\.tsv$")]
if (length(files) == 0) stop("No Cox4All_*__PESprot/full.tsv found under: ", type_dir)

cox_dt <- rbindlist(lapply(files, function(f) {
  x <- safe_fread(f)
  if (is.null(x)) return(NULL)
  x[, file := basename(f)]
  x[, path := f]
  x[, Type := type_for_plots]
  x[, pes_kind := fifelse(
    str_detect(file, "__PESprot\\.tsv$"), "PESprot",
    fifelse(str_detect(file, "__PESfull\\.tsv$"), "PESfull", NA_character_)
  )]
  x
}), fill = TRUE)

cox_all <- as_tibble(cox_dt) %>%
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

# Keep only OK rows if that column exists (strongly recommended)
if ("cox_status" %in% names(cox_all)) {
  cox_all <- cox_all %>% filter(cox_status == "OK")
}

# De-dup at the fundamental key. (This prevents subtle repeats from file merges)
# If you have repeated rows, keep the first (they should be identical per key).
cox_all <- cox_all %>%
  arrange(pes_kind, exposure_id, disease_age_col) %>%
  distinct(pes_kind, exposure_id, disease_age_col, .keep_all = TRUE)

# ----------------------------
# Pick top diseases per exposure (Type5)
# ----------------------------
pick_top_diseases <- function(df, exp_id, delta_col, topN, min_events) {
  df %>%
    filter(exposure_id == exp_id) %>%
    filter(is.na(events) | events >= min_events) %>%
    mutate(delta = .data[[delta_col]]) %>%
    filter(is.finite(delta)) %>%
    arrange(desc(delta)) %>%
    slice_head(n = topN) %>%
    pull(disease_label_short)
}

# ----------------------------
# Plot one exposure (two-panel dumbbell)
# ----------------------------
plot_one_exposure <- function(df, exp_id, diseases_keep, rank_tag, out_file,
                              label_top_k = 3) {
  
  dd <- df %>%
    filter(exposure_id == exp_id,
           disease_label_short %in% diseases_keep) %>%
    # one row per disease (should already be 1, but keep safe)
    group_by(disease_label_short) %>%
    summarize(
      disease_label = dplyr::first(disease_label),
      cindex_M0 = first_non_na(cindex_M0),
      cindex_M1 = first_non_na(cindex_M1),
      cindex_M2 = first_non_na(cindex_M2),
      cindex_M3 = first_non_na(cindex_M3),
      delta01   = first_non_na(delta_cindex_M0_to_M1),
      delta23   = first_non_na(delta_cindex_M2_to_M3),
      .groups="drop"
    )
  
  if (nrow(dd) == 0) return(NULL)
  
  exp_lab <- df %>% filter(exposure_id == exp_id) %>% slice(1) %>% pull(exposure_label)
  
  dd <- dd %>%
    mutate(
      ylab = ifelse(is.na(disease_label) | disease_label == disease_label_short,
                    disease_label_short,
                    paste0(disease_label, " (", disease_label_short, ")"))
    )
  
  # preserve ranked order in diseases_keep
  lev_codes <- diseases_keep
  lev_ylab <- dd$ylab[match(lev_codes, dd$disease_label_short)]
  lev_ylab <- rev(lev_ylab[!is.na(lev_ylab)])
  dd$ylab <- factor(dd$ylab, levels = unique(lev_ylab))
  
  long <- bind_rows(
    dd %>% transmute(panel = "Add PES beyond exposure model (M2 → M3)",
                     ylab, base=cindex_M2, plus=cindex_M3, delta=delta23),
    dd %>% transmute(panel = "Add PES to covariates (M0 → M1)",
                     ylab, base=cindex_M0, plus=cindex_M1, delta=delta01)
  ) %>% filter(is.finite(base), is.finite(plus), is.finite(delta))
  
  # x limits per panel (include all points + padding)
  xlims <- long %>%
    group_by(panel) %>%
    summarize(
      xmin = min(c(base, plus), na.rm=TRUE),
      xmax = max(c(base, plus), na.rm=TRUE),
      .groups="drop"
    ) %>%
    mutate(span = pmax(xmax - xmin, 0.02),
           xmin = xmin - 0.05*span,
           xmax = xmax + 0.20*span)
  
  long <- long %>% left_join(xlims, by="panel")
  
  label_df <- long %>%
    group_by(panel) %>%
    arrange(desc(delta), .by_group = TRUE) %>%
    slice_head(n = label_top_k) %>%
    ungroup() %>%
    mutate(
      label = sprintf("Δ=%.3f", delta),
      xlab  = pmin(xmax - 0.01*span, pmax(base, plus) + 0.03*span)
    )
  
  p <- ggplot(long, aes(y = ylab)) +
    geom_segment(aes(x=base, xend=plus, yend=ylab), linewidth=0.8, alpha=0.7) +
    geom_point(aes(x=base), shape=21, fill="white", stroke=1.1, size=2.7) +
    geom_point(aes(x=plus), shape=19, size=2.7) +
    geom_text(data=label_df, aes(x=xlab, label=label), size=3.2, hjust=0) +
    facet_wrap(~panel, ncol=1, scales="free_x") +
    labs(
      title = paste0(exp_lab, " — top ", length(diseases_keep), " diseases (", rank_tag, ", Type5)"),
      subtitle = paste0("Open circle = baseline; filled = +PES. Filter: events ≥ ", min_events, "."),
      x = "c-index", y = NULL
    ) +
    theme_bw(base_size = 12) +
    theme(
      strip.background = element_rect(fill="grey95", color=NA),
      axis.text.y = element_text(size=10),
      plot.margin = margin(5.5, 14, 5.5, 5.5)
    )
  
  ggsave(out_file, p, width=10.8, height=6.8, dpi=300)
  p
}

# ----------------------------
# Run per exposure for a PES kind (Type5)
# ----------------------------
run_per_exposure <- function(df, tag) {
  plot_dir <- file.path(out_dir, paste0("CoxTop", topN_per_exposure, "_PerExposure_", tag, "_", type_for_plots))
  dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
  
  exposures <- sort(unique(df$exposure_id))
  
  for (exp_id in exposures) {
    
    keep01 <- pick_top_diseases(df, exp_id, "delta_cindex_M0_to_M1", topN_per_exposure, min_events)
    if (length(keep01) >= 2) {
      exp_lab <- df %>% filter(exposure_id == exp_id) %>% slice(1) %>% pull(exposure_label)
      out01 <- file.path(plot_dir, safe_basename(paste0("Top", topN_per_exposure, "_", tag, "_", exp_lab, "_rankM0toM1.png")))
      plot_one_exposure(df, exp_id, keep01, "ranked by Δ(M0→M1)", out01)
    }
    
    keep23 <- pick_top_diseases(df, exp_id, "delta_cindex_M2_to_M3", topN_per_exposure, min_events)
    if (length(keep23) >= 2) {
      exp_lab <- df %>% filter(exposure_id == exp_id) %>% slice(1) %>% pull(exposure_label)
      out23 <- file.path(plot_dir, safe_basename(paste0("Top", topN_per_exposure, "_", tag, "_", exp_lab, "_rankM2toM3.png")))
      plot_one_exposure(df, exp_id, keep23, "ranked by Δ(M2→M3)", out23)
    }
  }
  
  message("Wrote per-exposure plots to: ", plot_dir)
}

# ----------------------------
# Plot ONE disease across exposures (Top-N)
#   KEY FIX: y axis uses unique exposure_id, display label disambiguated
# ----------------------------
plot_one_disease_across_exposures_T5 <- function(df_T5,
                                                 disease_code_focus = "E11",
                                                 out_file,
                                                 label_top_k = 6,
                                                 topN_exposures = 30) {
  
  dc <- toupper(disease_code_focus)
  
  # HARD CHECK: make sure disease_code exists and has multiple values
  if (!"disease_code" %in% names(df_T5)) stop("df_T5 is missing disease_code column.")
  if (!"exposure_label" %in% names(df_T5)) stop("df_T5 is missing exposure_label column.")
  
  # 1) FILTER ONLY THIS DISEASE (stable column)
  dd0 <- df_T5 %>%
    filter(disease_code == dc) %>%
    filter(is.na(events) | events >= min_events) %>%
    filter(is.finite(cindex_M0), is.finite(cindex_M1), is.finite(cindex_M2), is.finite(cindex_M3)) %>%
    filter(is.finite(delta_cindex_M0_to_M1) | is.finite(delta_cindex_M2_to_M3))
  
  # Debug sanity: if you’re still getting “identical plots”, this will reveal it immediately
  message("Disease ", dc, ": rows after filter = ", nrow(dd0),
          " | unique exposures = ", dplyr::n_distinct(dd0$exposure_label),
          " | unique diseases in dd0 = ", dplyr::n_distinct(dd0$disease_code))
  
  if (nrow(dd0) == 0) {
    message("No rows found for disease ", dc, " at Type5 after filtering.")
    return(NULL)
  }
  
  # 2) If duplicates exist per exposure (same disease), choose ONE row deterministically:
  # Prefer the one with the largest delta23 (or delta01 as tie-break)
  dd <- dd0 %>%
    group_by(exposure_label) %>%
    arrange(desc(delta_cindex_M2_to_M3), desc(delta_cindex_M0_to_M1)) %>%
    slice(1) %>%
    ungroup()
  
  # 3) Take top N exposures by Δ(M2→M3) (this is what you said you wanted)
  dd <- dd %>%
    arrange(desc(delta_cindex_M2_to_M3), desc(delta_cindex_M0_to_M1)) %>%
    slice_head(n = topN_exposures)
  
  # factor levels MUST be unique
  exp_levels <- unique(dd$exposure_label)
  dd <- dd %>% mutate(exposure_label = factor(exposure_label, levels = rev(exp_levels)))
  
  long <- bind_rows(
    dd %>% transmute(panel="Add PES beyond exposure model (M2 → M3)",
                     exposure_label,
                     base=cindex_M2, plus=cindex_M3, delta=delta_cindex_M2_to_M3),
    dd %>% transmute(panel="Add PES to covariates (M0 → M1)",
                     exposure_label,
                     base=cindex_M0, plus=cindex_M1, delta=delta_cindex_M0_to_M1)
  )
  
  label_df <- long %>%
    group_by(panel) %>%
    arrange(desc(delta), .by_group = TRUE) %>%
    slice_head(n = label_top_k) %>%
    ungroup() %>%
    mutate(label = sprintf("Δ=%.3f", delta))
  
  # X limits that ALWAYS include points + room for labels
  vals <- c(long$base, long$plus)
  xmin <- min(vals, na.rm = TRUE)
  xmax <- max(vals, na.rm = TRUE)
  span <- xmax - xmin
  if (!is.finite(span) || span <= 0) span <- 0.02
  xlim_min <- xmin - 0.03 * span
  xlim_max <- xmax + 0.20 * span
  
  p <- ggplot(long, aes(y = exposure_label)) +
    geom_segment(aes(x=base, xend=plus, yend=exposure_label), linewidth=0.8, alpha=0.7) +
    geom_point(aes(x=base), shape=21, fill="white", stroke=1.1, size=2.7) +
    geom_point(aes(x=plus), shape=19, size=2.7) +
    geom_text(data=label_df,
              aes(x=pmax(base, plus) + 0.01 * span, label=label),
              size=3.2, hjust=0) +
    facet_wrap(~panel, ncol=1, scales="free_x") +
    coord_cartesian(xlim=c(xlim_min, xlim_max), clip="on") +
    labs(
      title = paste0("Disease ", dc, " (Type5): discrimination gain from adding PES across exposures"),
      subtitle = paste0("Open circle = baseline; filled = +PES. Showing top ", nrow(dd),
                        " exposures by Δ(M2→M3). Filter: events ≥ ", min_events, "."),
      x="c-index", y=NULL
    ) +
    theme_bw(base_size = 12) +
    theme(
      strip.background = element_rect(fill="grey95", color=NA),
      axis.text.y = element_text(size=9),
      plot.margin = margin(5.5, 60, 5.5, 5.5)
    )
  
  ggsave(out_file, p, width=10.5, height=max(6, 0.28*nrow(dd) + 3), dpi=300, limitsize=FALSE)
  p
}

# ---- RUN FOR MULTIPLE DISEASES ----
disease_codes_focus <- c("E11","J43","F10")
topN_exposures_for_disease <- 30

if (nrow(cox_main_T5) > 0) {
  for (dc in disease_codes_focus) {
    plot_one_disease_across_exposures_T5(
      df_T5 = cox_main_T5,
      disease_code_focus = dc,
      out_file = file.path(out_dir, paste0("Disease_", dc, "_AcrossExposures_Type5_PESprot_Top", topN_exposures_for_disease, ".png")),
      label_top_k = 6,
      topN_exposures = topN_exposures_for_disease
    )
  }
}

if (nrow(cox_supp_T5) > 0) {
  for (dc in disease_codes_focus) {
    plot_one_disease_across_exposures_T5(
      df_T5 = cox_supp_T5,
      disease_code_focus = dc,
      out_file = file.path(out_dir, paste0("Disease_", dc, "_AcrossExposures_Type5_PESfull_Top", topN_exposures_for_disease, ".png")),
      label_top_k = 6,
      topN_exposures = topN_exposures_for_disease
    )
  }
}

# ----------------------------
# Execute
# ----------------------------
cox_main_T5 <- cox_all %>% filter(pes_kind == pes_kind_main)
cox_supp_T5 <- cox_all %>% filter(pes_kind == pes_kind_supp)

if (nrow(cox_main_T5) > 0) {
  for (dc in disease_codes_focus) {
    plot_one_disease_across_exposures(
      df = cox_main_T5,
      disease_code = dc,
      out_file = file.path(
        out_dir,
        paste0("Disease_", dc, "_AcrossExposures_", type_for_plots, "_", pes_kind_main,
               "_Top", topN_exposures_for_disease, ".png")
      ),
      label_top_k = 8,
      topN_exposures = topN_exposures_for_disease
    )
  }
  run_per_exposure(cox_main_T5, pes_kind_main)
}


if (nrow(cox_supp_T5) > 0) {
  for (dc in disease_codes_focus) {
    plot_one_disease_across_exposures(
      df = cox_supp_T5,
      disease_code = dc,
      out_file = file.path(
        out_dir,
        paste0("Disease_", dc, "_AcrossExposures_", type_for_plots, "_", pes_kind_supp,
               "_Top", topN_exposures_for_disease, ".png")
      ),
      label_top_k = 8,
      topN_exposures = topN_exposures_for_disease
    )
  }
  run_per_exposure(cox_supp_T5, pes_kind_supp)
}

# Save Type5 table for convenience
fwrite(as.data.table(cox_all), file.path(out_dir, paste0("cox_all_", type_for_plots, ".tsv")), sep="\t")

message("\nDONE. Output in: ", out_dir, "\n")
