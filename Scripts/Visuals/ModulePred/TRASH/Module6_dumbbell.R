#!/usr/bin/env Rscript

# ============================================================
# Cox Top-10 "dumbbell" plots per exposure
# - Goal: show how much PES adds to disease discrimination
#   (baseline c-index vs +PES c-index), like your example figure.
#
# For each exposure:
#   (A) rank diseases by Δc-index (M0->M1) at Type5 (cov -> cov+PES)
#   (B) rank diseases by Δc-index (M2->M3) at Type5 (cov+E -> cov+E+PES)
#   Then plot those same diseases across Types (Type1..Type5) OR Type5-only.
#
# Inputs expected (per Type folder under base_dir):
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
  library(purrr)
})

# ----------------------------
# USER CONFIG
# ----------------------------
base_dir <- "/n/groups/patel/shakson_ukb/UK_Biobank/Data/Parallel/PES_test"
types <- paste0("Type", 1:5)

out_path <- "/n/groups/patel/shakson_ukb/UK_Biobank/Output/PES/"
out_dir  <- file.path(out_path)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Which PES to plot as "main" vs "supp"
pes_kind_main <- "PESprot"
pes_kind_supp <- "PESfull"

# Ranking settings
type_for_ranking <- "Type5"
topN_per_exposure <- 10

# Filtering (recommended to avoid tiny-event noise)
min_events <- 50

# Plot options
plot_types_mode <- c("ALL_TYPES", "TYPE5_ONLY")[1]   # change to "TYPE5_ONLY" if you want a compact main fig style
show_types <- if (plot_types_mode == "TYPE5_ONLY") type_for_ranking else types

# Label shortening (for long ICD10 derived strings)
use_disease_code_only <- TRUE  # TRUE => show ICD10 code extracted from age_... ; FALSE => show your curated labels if present

# ----------------------------
# Optional label maps (keep your existing)
# ----------------------------
exposure_label_map <- tibble::tribble(
  ~exposure_id, ~exposure_label,
  "alcohol_intake_frequency_f1558_0_0", "Alcohol frequency",
  "smoking_status_f20116_0_0_Current", "Smoking (current)",
  "fresh_fruit_intake_f1309_0_0", "Fresh fruit",
  "processed_meat_intake_f1349_0_0", "Processed meat",
  "summed_met_minutes_per_week_for_all_activity_f22040_0_0", "Physical activity (MET-min/wk)",
  "types_of_physical_activity_in_last_4_weeks_f6164_0_0.multi_Strenuous_sports", "Strenuous sport (type)",
  "time_spent_watching_television_tv_f1070_0_0", "TV time"
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

# Example: age_e11_first_reported_... -> E11
extract_disease_code <- function(disease_age_col) {
  x <- tolower(as.character(disease_age_col))
  code <- str_match(x, "^age_([^_]+)_first_reported_")[,2]
  code <- ifelse(is.na(code), str_match(x, "^age_([^_]+)_")[,2], code)
  toupper(code)
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
    exposure_id = as.character(exposure_id),
    disease_age_col = as.character(disease_age_col),
    disease_code = extract_disease_code(disease_age_col)
  ) %>%
  left_join(exposure_label_map, by="exposure_id") %>%
  left_join(disease_label_map, by="disease_age_col") %>%
  mutate(
    exposure_label = ifelse(is.na(exposure_label), exposure_id, exposure_label),
    disease_label  = ifelse(is.na(disease_label), disease_code, disease_label),
    disease_label_short = if (use_disease_code_only) disease_code else disease_label
  )

# Make sure numeric columns are numeric
cox_all <- cox_all %>%
  mutate(
    events = suppressWarnings(as.numeric(events)),
    cindex_M0 = suppressWarnings(as.numeric(cindex_M0)),
    cindex_M1 = suppressWarnings(as.numeric(cindex_M1)),
    cindex_M2 = suppressWarnings(as.numeric(cindex_M2)),
    cindex_M3 = suppressWarnings(as.numeric(cindex_M3)),
    delta_cindex_M0_to_M1 = suppressWarnings(as.numeric(delta_cindex_M0_to_M1)),
    delta_cindex_M2_to_M3 = suppressWarnings(as.numeric(delta_cindex_M2_to_M3))
  )

# ----------------------------
# Pick top diseases per exposure (rank once on Type5)
# ----------------------------
pick_top_diseases_per_exposure <- function(cox_df,
                                           delta_col = c("delta_cindex_M2_to_M3", "delta_cindex_M0_to_M1"),
                                           type_for_ranking = "Type5",
                                           topN = 10,
                                           min_events = 0) {
  delta_col <- match.arg(delta_col)
  
  dd <- cox_df %>%
    filter(Type == type_for_ranking) %>%
    mutate(
      delta = suppressWarnings(as.numeric(.data[[delta_col]])),
      events = suppressWarnings(as.numeric(events))
    ) %>%
    filter(is.finite(delta)) %>%
    filter(is.na(events) | events >= min_events) %>%
    group_by(exposure_label, disease_label_short) %>%
    summarize(delta = max(delta, na.rm=TRUE), .groups="drop")  # <-- de-dup
  
  dd %>%
    group_by(exposure_label) %>%
    arrange(desc(delta), .by_group = TRUE) %>%
    slice_head(n = topN) %>%
    ungroup() %>%
    select(exposure_label, disease_label_short) %>%
    distinct()
}

# ----------------------------
# Dumbbell plot (baseline vs +PES)
# ----------------------------
plot_dumbbell_topN <- function(cox_df, keep_tbl, out_file,
                               show_types = c("Type1","Type2","Type3","Type4","Type5"),
                               order_by = c("delta_cindex_M2_to_M3", "delta_cindex_M0_to_M1"),
                               title = NULL,
                               subtitle = NULL) {
  order_by <- match.arg(order_by)
  if (nrow(cox_df) == 0 || nrow(keep_tbl) == 0) return(NULL)
  
  dd <- cox_df %>%
    filter(Type %in% show_types) %>%
    inner_join(keep_tbl, by = c("exposure_label","disease_label_short")) %>%
    filter(is.finite(cindex_M0) & is.finite(cindex_M1) & is.finite(cindex_M2) & is.finite(cindex_M3))
  
  if (nrow(dd) == 0) return(NULL)
  
  # Order diseases within each exposure using Type5 ranking metric (stable y-axis)
  ord_tbl <- cox_df %>%
    filter(Type == type_for_ranking) %>%
    inner_join(keep_tbl, by=c("exposure_label","disease_label_short")) %>%
    mutate(ord = suppressWarnings(as.numeric(.data[[order_by]]))) %>%
    filter(is.finite(ord)) %>%
    group_by(exposure_label, disease_label_short) %>%
    summarize(ord = max(ord, na.rm=TRUE), .groups="drop") %>%     # <-- de-dup
    group_by(exposure_label) %>%
    arrange(ord, .by_group = TRUE) %>%
    summarize(disease_levels = list(unique(disease_label_short)), .groups="drop")  # <-- unique
  
  
  
  dd <- dd %>%
    left_join(ord_tbl, by="exposure_label") %>%
    rowwise() %>%
    mutate(disease_label_short = factor(disease_label_short, levels = unique(disease_levels))) %>%
    ungroup()
  
  
  long <- bind_rows(
    dd %>% transmute(exposure_label, Type, disease_label_short,
                     panel="Add PES to covariates",
                     x0=cindex_M0, x1=cindex_M1),
    dd %>% transmute(exposure_label, Type, disease_label_short,
                     panel="Add PES beyond exposure model",
                     x0=cindex_M2, x1=cindex_M3)
  )
  
  # One plot: facet by exposure (free y-space) and panel; columns = Type
  p <- ggplot(long, aes(y = disease_label_short)) +
    geom_segment(aes(x = x0, xend = x1, yend = disease_label_short),
                 linewidth = 0.7, alpha = 0.7) +
    geom_point(aes(x = x0), size = 2.3) +
    geom_point(aes(x = x1), size = 2.3) +
    facet_grid(exposure_label + panel ~ Type, scales="free_y", space="free_y") +
    labs(
      title = title %||% "Top disease gains in c-index from adding PES",
      subtitle = subtitle %||% paste0(
        "Top ", topN_per_exposure, " diseases per exposure selected by ",
        order_by, " at ", type_for_ranking, ". "
      ),
      x = "c-index", y = NULL
    ) +
    theme(
      strip.background = element_rect(fill="grey95", color=NA),
      axis.text.y = element_text(size=9)
    )
  
  # Dynamic height: exposures * diseases rows (reasonable scaling)
  n_exp <- length(unique(long$exposure_label))
  height <- min(4 + n_exp * 2.2, 40)
  
  ggsave(out_file, p, width = 16, height = height, dpi = 300)
  p
}

`%||%` <- function(x, y) if (!is.null(x) && length(x) > 0) x else y

# ----------------------------
# Runner: generate plots for MAIN and SUPP
# ----------------------------
run_top10_suite <- function(cox_df, tag, out_dir) {
  if (nrow(cox_df) == 0) return(invisible(NULL))
  
  plot_dir <- file.path(out_dir, paste0("CoxTop10_", tag))
  dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Rank sets (one stable set per exposure, based on Type5)
  keep_m0m1 <- pick_top_diseases_per_exposure(
    cox_df, delta_col = "delta_cindex_M0_to_M1",
    type_for_ranking = type_for_ranking,
    topN = topN_per_exposure,
    min_events = min_events
  )
  keep_m2m3 <- pick_top_diseases_per_exposure(
    cox_df, delta_col = "delta_cindex_M2_to_M3",
    type_for_ranking = type_for_ranking,
    topN = topN_per_exposure,
    min_events = min_events
  )
  
  # Plot: set ranked by ΔM0->M1 (cov -> cov+PES), show both panels (M0->M1 and M2->M3) anyway
  plot_dumbbell_topN(
    cox_df, keep_m0m1,
    out_file = file.path(plot_dir, paste0("Fig_Top10_byExposure_rankM0toM1_", plot_types_mode, ".png")),
    show_types = show_types,
    order_by = "delta_cindex_M0_to_M1",
    title = paste0("Top ", topN_per_exposure, " disease gains per exposure (ranked by Δc M0→M1)"),
    subtitle = paste0("Ranking at ", type_for_ranking, " using Δc-index (M0→M1). Points show baseline vs +PES. min events ≥ ", min_events, ".")
  )
  
  # Plot: set ranked by ΔM2->M3 (cov+E -> cov+E+PES)
  plot_dumbbell_topN(
    cox_df, keep_m2m3,
    out_file = file.path(plot_dir, paste0("Fig_Top10_byExposure_rankM2toM3_", plot_types_mode, ".png")),
    show_types = show_types,
    order_by = "delta_cindex_M2_to_M3",
    title = paste0("Top ", topN_per_exposure, " disease gains per exposure (ranked by Δc M2→M3)"),
    subtitle = paste0("Ranking at ", type_for_ranking, " using Δc-index (M2→M3). Points show baseline vs +PES. min events ≥ ", min_events, ".")
  )
  
  # Save the keep tables so you can cite them / reuse downstream
  fwrite(as.data.table(keep_m0m1), file.path(plot_dir, "Top10_keep_rankM0toM1.tsv"), sep="\t")
  fwrite(as.data.table(keep_m2m3), file.path(plot_dir, "Top10_keep_rankM2toM3.tsv"), sep="\t")
  
  message("Wrote Top10 plots + keep tables to: ", plot_dir)
}

# ----------------------------
# Run for MAIN and SUPP
# ----------------------------
cox_main <- filter_pes(cox_all, pes_kind_main)
cox_supp <- filter_pes(cox_all, pes_kind_supp)

if (nrow(cox_main) > 0) run_top10_suite(cox_main, tag = pes_kind_main, out_dir = out_dir)
if (nrow(cox_supp) > 0) run_top10_suite(cox_supp, tag = pes_kind_supp, out_dir = out_dir)

# Also save the combined Cox table for convenience
fwrite(as.data.table(cox_all), file.path(out_dir, "cox_all.tsv"), sep="\t")

message("\nDONE. Output in: ", out_dir, "\n")
