#!/usr/bin/env Rscript

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
types <- paste0("Type", 1:5)

out_dir  <- "/n/groups/patel/shakson_ukb/UK_Biobank/Output/PES/"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

type_for_plots    <- "Type5"
topN_per_exposure <- 10
min_events        <- 50

pes_kind_main <- "PESprot"
pes_kind_supp <- "PESfull"

use_disease_code_only <- TRUE

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
  "age_e11_first_reported_non_insulin_dependent_diabetes_mellitus_f130708_0_0", "T2D"
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

read_all_cox <- function(base_dir, types) {
  out <- list()
  for (ty in types) {
    ty_dir <- file.path(base_dir, ty)
    if (!dir.exists(ty_dir)) next
    files <- list.files(ty_dir, full.names = TRUE)
    files <- files[str_detect(basename(files), "^Cox4All_.*__PES(prot|full)\\.tsv$")]
    if (!length(files)) next
    
    dt <- rbindlist(lapply(files, function(f) {
      x <- safe_fread(f)
      if (is.null(x)) return(NULL)
      x[, file := basename(f)]
      x[, Type := ty]
      x[, pes_kind := fifelse(str_detect(file, "__PESprot\\.tsv$"), "PESprot",
                              fifelse(str_detect(file, "__PESfull\\.tsv$"), "PESfull", NA_character_))]
      x
    }), fill = TRUE)
    
    out[[ty]] <- dt
  }
  rbindlist(out, fill = TRUE)
}

# Core: enforce ONE row per key deterministically (NO max())
dedup_one_row <- function(df, key_cols) {
  df %>%
    arrange(desc(events), desc(n), file) %>%
    group_by(across(all_of(key_cols))) %>%
    slice(1) %>%
    ungroup()
}

# ----------------------------
# Load + clean
# ----------------------------
cox_dt <- read_all_cox(base_dir, types)
if (!nrow(cox_dt)) stop("No Cox4All files found under: ", base_dir)

cox <- as_tibble(cox_dt) %>%
  mutate(
    Type = as.character(Type),
    exposure_id = as.character(exposure_id),
    disease_age_col = as.character(disease_age_col),
    disease_code = extract_disease_code(disease_age_col)
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
    n      = suppressWarnings(as.numeric(n)),
    cindex_M0 = suppressWarnings(as.numeric(cindex_M0)),
    cindex_M1 = suppressWarnings(as.numeric(cindex_M1)),
    cindex_M2 = suppressWarnings(as.numeric(cindex_M2)),
    cindex_M3 = suppressWarnings(as.numeric(cindex_M3)),
    delta_cindex_M0_to_M1 = suppressWarnings(as.numeric(delta_cindex_M0_to_M1)),
    delta_cindex_M2_to_M3 = suppressWarnings(as.numeric(delta_cindex_M2_to_M3))
  )

# ----------------------------
# IMPORTANT FILTER: Type5 ONLY, OK ONLY, exact pes_kind ONLY
# ----------------------------
prep_df <- function(df, pes_kind_keep) {
  df %>%
    filter(Type == type_for_plots) %>%
    filter(pes_kind == pes_kind_keep) %>%          # <-- no NA rows
    filter(cox_status == "OK") %>%                 # <-- drop ERROR_OR_NULL
    filter(is.na(events) | events >= min_events) %>%
    # enforce uniqueness of (exposure_id, disease_age_col)
    dedup_one_row(c("exposure_id", "disease_age_col"))
}

# ----------------------------
# Plot functions
# ----------------------------
plot_one_disease_across_exposures <- function(df_T5, disease_code="E11", out_file, topN=NULL, label_top_k=6) {
  dd <- df_T5 %>%
    filter(disease_label_short == disease_code) %>%
    # enforce uniqueness (again, but now key is exposure_id + disease)
    dedup_one_row(c("exposure_id", "disease_age_col"))
  
  if (!nrow(dd)) return(NULL)
  
  dd <- dd %>%
    mutate(delta01 = delta_cindex_M0_to_M1,
           delta23 = delta_cindex_M2_to_M3) %>%
    arrange(desc(delta23), desc(delta01))
  
  if (!is.null(topN)) dd <- dd %>% slice_head(n = topN)
  
  # Use exposure_id as unique y key
  dd <- dd %>%
    mutate(exposure_id = factor(exposure_id, levels = rev(unique(exposure_id))))
  
  long <- bind_rows(
    dd %>% transmute(panel="Add PES beyond exposure model (M2 → M3)",
                     exposure_id, exposure_label,
                     base=cindex_M2, plus=cindex_M3, delta=delta23),
    dd %>% transmute(panel="Add PES to covariates (M0 → M1)",
                     exposure_id, exposure_label,
                     base=cindex_M0, plus=cindex_M1, delta=delta01)
  ) %>% filter(is.finite(base), is.finite(plus), is.finite(delta))
  
  vals <- c(long$base, long$plus)
  xmin <- min(vals); xmax <- max(vals)
  span <- xmax - xmin; if (!is.finite(span) || span <= 0) span <- 0.05
  xlim <- c(xmin - 0.03*span, xmax + 0.22*span)
  
  label_df <- long %>%
    group_by(panel) %>%
    arrange(desc(delta), .by_group=TRUE) %>%
    slice_head(n=label_top_k) %>%
    ungroup() %>%
    mutate(label=sprintf("Δ=%.3f", delta),
           xlab = plus + 0.02*span)
  
  p <- ggplot(long, aes(y=exposure_id)) +
    geom_segment(aes(x=base, xend=plus, yend=exposure_id), linewidth=0.8, alpha=0.7) +
    geom_point(aes(x=base), shape=21, fill="white", stroke=1.1, size=2.8) +
    geom_point(aes(x=plus), shape=19, size=2.8) +
    geom_text(data=label_df, aes(x=xlab, label=label), hjust=0, size=3.2) +
    scale_y_discrete(labels=setNames(dd$exposure_label, as.character(dd$exposure_id))) +
    facet_wrap(~panel, ncol=1, scales="fixed") +
    coord_cartesian(xlim=xlim, clip="on") +
    labs(
      title=paste0("Disease ", disease_code, " (", type_for_plots, "): discrimination gain from adding PES across exposures"),
      subtitle=paste0("Open circle = baseline; filled = +PES. Filter: events ≥ ", min_events, "."),
      x="c-index", y=NULL
    ) +
    theme_bw(base_size=12) +
    theme(strip.background=element_rect(fill="grey95", color=NA),
          axis.text.y=element_text(size=9),
          plot.margin=margin(5.5, 14, 5.5, 5.5))
  
  h <- max(6.5, 0.32*nrow(dd) + 3)
  h <- min(h, 18)  # cap at 18 inches (adjust as you like)
  
  
  ggsave(out_file, p, width=10.8, height=h, dpi=300)
  p
}

pick_top_diseases <- function(df_T5, exposure_id_one, delta_col, topN=10) {
  df_T5 %>%
    filter(exposure_id == exposure_id_one) %>%
    mutate(delta = .data[[delta_col]]) %>%
    arrange(desc(delta)) %>%
    slice_head(n=topN) %>%
    pull(disease_label_short)
}

plot_one_exposure_top10 <- function(df_T5, exposure_id_one, diseases_keep, rank_tag, out_file, label_top_k=3) {
  dd <- df_T5 %>%
    filter(exposure_id == exposure_id_one,
           disease_label_short %in% diseases_keep) %>%
    dedup_one_row(c("exposure_id", "disease_age_col"))
  
  if (!nrow(dd)) return(NULL)
  
  dd <- dd %>%
    mutate(ylab = ifelse(disease_label_short == disease_label, disease_label_short,
                         paste0(disease_label, " (", disease_label_short, ")")))
  
  # preserve rank order from diseases_keep
  y_order <- dd %>% mutate(idx=match(disease_label_short, diseases_keep)) %>%
    arrange(idx) %>% pull(ylab)
  dd$ylab <- factor(dd$ylab, levels=rev(unique(y_order)))
  
  dd <- dd %>%
    mutate(delta01 = delta_cindex_M0_to_M1,
           delta23 = delta_cindex_M2_to_M3)
  
  long <- bind_rows(
    dd %>% transmute(panel="Add PES beyond exposure model (M2 → M3)",
                     ylab, base=cindex_M2, plus=cindex_M3, delta=delta23),
    dd %>% transmute(panel="Add PES to covariates (M0 → M1)",
                     ylab, base=cindex_M0, plus=cindex_M1, delta=delta01)
  ) %>% filter(is.finite(base), is.finite(plus), is.finite(delta))
  
  vals <- c(long$base, long$plus)
  xmin <- min(vals); xmax <- max(vals)
  span <- xmax - xmin; if (!is.finite(span) || span <= 0) span <- 0.05
  xlim <- c(xmin - 0.03*span, xmax + 0.22*span)
  
  label_df <- long %>%
    group_by(panel) %>%
    arrange(desc(delta), .by_group=TRUE) %>%
    slice_head(n=label_top_k) %>%
    ungroup() %>%
    mutate(label=sprintf("Δ=%.3f", delta),
           xlab = plus + 0.02*span)
  
  exp_lab <- df_T5 %>% filter(exposure_id==exposure_id_one) %>% slice(1) %>% pull(exposure_label)
  
  p <- ggplot(long, aes(y=ylab)) +
    geom_segment(aes(x=base, xend=plus, yend=ylab), linewidth=0.8, alpha=0.7) +
    geom_point(aes(x=base), shape=21, fill="white", stroke=1.1, size=2.8) +
    geom_point(aes(x=plus), shape=19, size=2.8) +
    geom_text(data=label_df, aes(x=xlab, label=label), hjust=0, size=3.2) +
    facet_wrap(~panel, ncol=1, scales="fixed") +
    coord_cartesian(xlim=xlim, clip="on") +
    labs(
      title=paste0(exp_lab, " — top ", length(diseases_keep), " diseases (", rank_tag, ", ", type_for_plots, ")"),
      subtitle=paste0("Open circle = baseline; filled = +PES. Filter: events ≥ ", min_events, "."),
      x="c-index", y=NULL
    ) +
    theme_bw(base_size=12) +
    theme(strip.background=element_rect(fill="grey95", color=NA),
          axis.text.y=element_text(size=10),
          plot.margin=margin(5.5, 14, 5.5, 5.5))
  
  ggsave(out_file, p, width=10.8, height=6.8, dpi=300)
  p
}

run_per_exposure <- function(df_T5, tag) {
  plot_dir <- file.path(out_dir, paste0("CoxTop10_PerExposure_", tag, "_", type_for_plots))
  dir.create(plot_dir, showWarnings=FALSE, recursive=TRUE)
  
  for (eid in sort(unique(df_T5$exposure_id))) {
    keep01 <- pick_top_diseases(df_T5, eid, "delta_cindex_M0_to_M1", topN_per_exposure)
    keep23 <- pick_top_diseases(df_T5, eid, "delta_cindex_M2_to_M3", topN_per_exposure)
    
    exp_lab <- df_T5 %>% filter(exposure_id==eid) %>% slice(1) %>% pull(exposure_label)
    
    if (length(keep01) >= 2) {
      out01 <- file.path(plot_dir, safe_basename(paste0("Top", topN_per_exposure, "_", tag, "_", exp_lab, "_rankM0toM1.png")))
      plot_one_exposure_top10(df_T5, eid, keep01, "ranked by Δ(M0→M1)", out01)
    }
    if (length(keep23) >= 2) {
      out23 <- file.path(plot_dir, safe_basename(paste0("Top", topN_per_exposure, "_", tag, "_", exp_lab, "_rankM2toM3.png")))
      plot_one_exposure_top10(df_T5, eid, keep23, "ranked by Δ(M2→M3)", out23)
    }
  }
  message("Wrote per-exposure plots to: ", plot_dir)
}

# ----------------------------
# Execute
# ----------------------------
cox_main_T5 <- prep_df(cox, pes_kind_main)
cox_supp_T5 <- prep_df(cox, pes_kind_supp)

# E11 plots
if (nrow(cox_main_T5)) {
  plot_one_disease_across_exposures(
    cox_main_T5, "E11",
    file.path(out_dir, "Disease_E11_AcrossExposures_Type5_PESprot.png"),
    topN=NULL, label_top_k=6
  )
}
if (nrow(cox_supp_T5)) {
  plot_one_disease_across_exposures(
    cox_supp_T5, "E11",
    file.path(out_dir, "Disease_E11_AcrossExposures_Type5_PESfull.png"),
    topN=NULL, label_top_k=6
  )
}

# Per-exposure top10 disease plots
if (nrow(cox_main_T5)) run_per_exposure(cox_main_T5, pes_kind_main)
if (nrow(cox_supp_T5)) run_per_exposure(cox_supp_T5, pes_kind_supp)

# save the exact df used for plotting (crucial for sanity checks)
fwrite(as.data.table(cox_main_T5), file.path(out_dir, "cox_T5_used_PESprot.tsv"), sep="\t")
fwrite(as.data.table(cox_supp_T5), file.path(out_dir, "cox_T5_used_PESfull.tsv"), sep="\t")

message("\nDONE. Output in: ", out_dir, "\n")

