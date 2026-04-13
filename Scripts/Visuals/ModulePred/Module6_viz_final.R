#!/usr/bin/env Rscript

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

out_dir <- "/n/groups/patel/shakson_ukb/UK_Biobank/Output/PES/"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

pes_kind_main <- "PESprot"
pes_kind_supp <- "PESfull"

# Optional: exposure label map (keep if you want human-friendly legend labels)
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

# ----------------------------
# Helpers
# ----------------------------
safe_fread <- function(path) tryCatch(fread(path), error = function(e) NULL)

theme_set(theme_bw(base_size = 12))
theme_legend_compact <- function() {
  theme(
    legend.key.size = grid::unit(0.35, "cm"),
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 10),
    legend.spacing.y = grid::unit(0.05, "cm"),
    legend.box.spacing = grid::unit(0.1, "cm")
  )
}

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

choose_primary_metric <- function(overall_long) {
  overall_long %>%
    filter(metric_name %in% c("AUC","R2","R2_code"),
           model %in% c("cov_only","prot_plus_cov")) %>%
    group_by(exposure_id) %>%
    summarize(
      has_auc = any(metric_name=="AUC" & is.finite(value)),
      has_r2  = any(metric_name=="R2" & is.finite(value)),
      has_r2c = any(metric_name=="R2_code" & is.finite(value)),
      .groups="drop"
    ) %>%
    mutate(primary_metric = case_when(
      has_auc ~ "AUC",
      has_r2  ~ "R2",
      has_r2c ~ "R2_code",
      TRUE    ~ NA_character_
    ))
}

filter_pes <- function(df, pes_kind) {
  if (!"pes_kind" %in% names(df)) return(df)
  df %>% filter(is.na(pes_kind) | pes_kind == !!pes_kind)
}

p_to_star <- function(p) {
  p <- suppressWarnings(as.numeric(p))
  ifelse(is.na(p), "",
         ifelse(p < 1e-3, "***",
                ifelse(p < 1e-2, "**",
                       ifelse(p < 5e-2, "*", ""))))
}

extract_disease_code <- function(disease_age_col) {
  x <- tolower(as.character(disease_age_col))
  code <- str_match(x, "^age_([^_]+)_first_reported_")[,2]
  code <- ifelse(is.na(code), str_match(x, "^age_([^_]+)_")[,2], code)
  toupper(code)
}

# ----------------------------
# Load only what we need
# ----------------------------
overall_dt <- read_type_files(base_dir, types, "_OverallMetrics\\.tsv$")
fold_dt    <- read_type_files(base_dir, types, "_FoldMetrics\\.tsv$")

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

if (nrow(overall_dt) == 0) stop("No *_OverallMetrics.tsv found under: ", base_dir)
if (nrow(fold_dt)    == 0) stop("No *_FoldMetrics.tsv found under: ", base_dir)
if (nrow(cox_dt)     == 0) stop("No Cox4All_*__PESprot/full.tsv found under: ", base_dir)

overall_all <- as_tibble(overall_dt) %>% mutate(Type = factor(Type, levels=types))
fold_all    <- as_tibble(fold_dt)    %>% mutate(Type = factor(Type, levels=types))
cox_all     <- as_tibble(cox_dt)     %>% mutate(Type = factor(Type, levels=types))

# label exposures (if not found, fall back to exposure_id)
overall_all <- overall_all %>%
  left_join(exposure_label_map, by="exposure_id") %>%
  mutate(exposure_label = ifelse(is.na(exposure_label), exposure_id, exposure_label))

fold_all <- fold_all %>%
  left_join(exposure_label_map, by="exposure_id") %>%
  mutate(exposure_label = ifelse(is.na(exposure_label), exposure_id, exposure_label))

cox_all <- cox_all %>%
  left_join(exposure_label_map, by="exposure_id") %>%
  mutate(
    exposure_label = ifelse(is.na(exposure_label), exposure_id, exposure_label),
    disease_code = extract_disease_code(disease_age_col),
    pair_label = paste0(exposure_label, " \u2192 ", disease_code)
  )

# ----------------------------
# PREP for Fig1/Fig2 (primary metric)
# ----------------------------
overall_metric_cols <- c(
  "r2","rmse","auc","logloss","r2_code","mse_code",
  "delta_r2","delta_auc","delta_logloss","delta_r2_code","delta_mse_code"
)

overall_long <- overall_all %>%
  mutate(
    model = as.character(model),
    exposure_id = as.character(exposure_id),
    exposure_type = as.character(exposure_type)
  ) %>%
  pivot_longer(cols = any_of(overall_metric_cols),
               names_to="metric_key", values_to="value") %>%
  mutate(
    value = suppressWarnings(as.numeric(value)),
    metric_name = case_when(
      metric_key %in% c("r2","delta_r2") ~ "R2",
      metric_key %in% c("auc","delta_auc") ~ "AUC",
      metric_key %in% c("r2_code","delta_r2_code") ~ "R2_code",
      metric_key %in% c("rmse") ~ "RMSE",
      metric_key %in% c("logloss","delta_logloss") ~ "LogLoss",
      metric_key %in% c("mse_code","delta_mse_code") ~ "MSE_code",
      TRUE ~ NA_character_
    ),
    kind = ifelse(str_starts(metric_key, "delta_"), "delta", "metric")
  ) %>%
  filter(!is.na(metric_name))

primary_tbl <- choose_primary_metric(overall_long %>% filter(kind=="metric"))

fold_candidates <- c(
  "r2_full","auc_full","r2_code_full",
  "r2_prot","auc_prot","r2_code_prot",
  "r2_cov","auc_cov","r2_code_cov"
)

fold_long <- fold_all %>%
  mutate(
    exposure_id = as.character(exposure_id),
    exposure_type = as.character(exposure_type),
    fold = suppressWarnings(as.integer(fold))
  ) %>%
  pivot_longer(cols = any_of(fold_candidates), names_to="metric_key", values_to="value") %>%
  mutate(
    value = suppressWarnings(as.numeric(value)),
    metric_name = case_when(
      str_starts(metric_key, "auc_") ~ "AUC",
      str_starts(metric_key, "r2_code_") ~ "R2_code",
      str_starts(metric_key, "r2_") ~ "R2",
      TRUE ~ NA_character_
    ),
    model = case_when(
      str_ends(metric_key, "_full") ~ "prot_plus_cov",
      str_ends(metric_key, "_cov")  ~ "cov_only",
      str_ends(metric_key, "_prot") ~ "prot_only",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(metric_name), !is.na(model), is.finite(value))

# Fig1: fold-based Δ(primary metric) = full - cov
fold_primary_delta <- fold_long %>%
  filter(model %in% c("prot_plus_cov","cov_only")) %>%
  inner_join(primary_tbl %>% filter(!is.na(primary_metric)), by="exposure_id") %>%
  filter(metric_name == primary_metric) %>%
  select(Type, exposure_id, exposure_label, exposure_type, metric_name, fold, model, value) %>%
  pivot_wider(names_from = model, values_from = value) %>%
  mutate(delta_metric = prot_plus_cov - cov_only) %>%
  filter(is.finite(delta_metric)) %>%
  group_by(Type, exposure_id, exposure_label, exposure_type, metric_name) %>%
  summarize(
    mean = mean(delta_metric, na.rm=TRUE),
    sd   = sd(delta_metric, na.rm=TRUE),
    n    = sum(is.finite(delta_metric)),
    se   = sd / sqrt(pmax(n, 1)),
    .groups="drop"
  ) %>%
  filter(n >= 2) %>%
  mutate(exposure_label = factor(exposure_label))

# Fig2: fold stability (full model, primary metric)
fold_primary_full <- fold_long %>%
  filter(model=="prot_plus_cov") %>%
  inner_join(primary_tbl %>% filter(!is.na(primary_metric)), by="exposure_id") %>%
  filter(metric_name == primary_metric) %>%
  group_by(Type, exposure_id, exposure_label, metric_name) %>%
  summarize(
    mean = mean(value, na.rm=TRUE),
    sd   = sd(value, na.rm=TRUE),
    n    = sum(is.finite(value)),
    se   = sd / sqrt(pmax(n,1)),
    .groups="drop"
  ) %>%
  filter(n >= 2) %>%
  mutate(exposure_label = factor(exposure_label))

# Fig3: # selected proteins (requires n_proteins_selected_full in FoldMetrics)
if (!"n_proteins_selected_full" %in% names(fold_all)) {
  stop("Fig3 requested but FoldMetrics is missing column: n_proteins_selected_full")
}
sel_primary <- fold_all %>%
  mutate(
    exposure_id = as.character(exposure_id),
    n_proteins_selected_full = suppressWarnings(as.numeric(n_proteins_selected_full))
  ) %>%
  filter(is.finite(n_proteins_selected_full)) %>%
  group_by(Type, exposure_id, exposure_label) %>%
  summarize(
    mean = mean(n_proteins_selected_full),
    sd   = sd(n_proteins_selected_full),
    n    = n(),
    se   = sd / sqrt(pmax(n,1)),
    .groups="drop"
  ) %>%
  filter(n >= 2) %>%
  mutate(exposure_label = factor(exposure_label))

# ----------------------------
# Write Fig1 / Fig2 / Fig3
# ----------------------------
p1 <- ggplot(fold_primary_delta,
             aes(x=Type, y=mean, color=exposure_label, group=exposure_label)) +
  geom_hline(yintercept=0, linetype=2) +
  geom_line(alpha=0.7, linewidth=0.7) +
  geom_point(size=2) +
  geom_errorbar(aes(ymin=mean-1.96*se, ymax=mean+1.96*se),
                width=0.12, alpha=0.75) +
  facet_wrap(~metric_name, scales="free_y") +
  labs(
    title="Incremental predictive value of proteins beyond covariates",
    subtitle="Δ = (prot_plus_cov) − (cov_only); primary metric per exposure (AUC > R2 > R2_code). Mean±95% CI across folds.",
    x="Covariate specification (Type)",
    y="Δ primary metric",
    color="Exposure"
  ) +
  theme(legend.position="right") +
  theme_legend_compact()

ggsave(file.path(out_dir, "Fig1_DeltaPrimaryMetric_lines.png"), p1, width=13, height=6, dpi=300)

p2 <- ggplot(fold_primary_full,
             aes(x=Type, y=mean, color=exposure_label, group=exposure_label)) +
  geom_line(linewidth=0.7) +
  geom_point(size=2) +
  geom_errorbar(aes(ymin=mean-1.96*se, ymax=mean+1.96*se),
                width=0.12, alpha=0.75) +
  facet_wrap(~metric_name, scales="free_y") +
  labs(
    title="Fold-to-fold stability of exposure prediction (full model)",
    subtitle="Mean±95% CI across folds (primary metric per exposure)",
    x="Covariate specification (Type)",
    y="Primary metric",
    color="Exposure"
  ) +
  theme(legend.position="right") +
  theme_legend_compact()

ggsave(file.path(out_dir, "Fig2_FoldStability_primaryMetric_lines.png"), p2, width=13, height=6, dpi=300)

p3 <- ggplot(sel_primary,
             aes(x=Type, y=mean, color=exposure_label, group=exposure_label)) +
  geom_line(linewidth=0.7) +
  geom_point(size=2) +
  geom_errorbar(aes(ymin=mean-1.96*se, ymax=mean+1.96*se),
                width=0.12, alpha=0.75) +
  labs(
    title="Number of selected proteins across folds (full model)",
    subtitle="Mean±95% CI across folds",
    x="Covariate specification (Type)",
    y="# proteins selected",
    color="Exposure"
  ) +
  theme(legend.position="right") +
  theme_legend_compact()

ggsave(file.path(out_dir, "Fig3_SelectedProteins_lines.png"), p3, width=13, height=6, dpi=300)

# ----------------------------
# Fig4 + FigS4 only: Cox heatmap
# ----------------------------
make_cox_heatmap <- function(cox_df, out_file, star_source="WALD", topN=30) {
  if (nrow(cox_df) == 0) return(invisible(NULL))
  
  cox_hm <- bind_rows(
    cox_df %>%
      mutate(
        model = "M1: cov + PES",
        HR = suppressWarnings(as.numeric(HR_PES_M1_perSD)),
        pval = if (star_source == "LRT" && "p_LRT_M0_vs_M1" %in% names(.)) p_LRT_M0_vs_M1 else p_PES_M1
      ) %>%
      select(Type, pair_label, model, HR, pval),
    cox_df %>%
      mutate(
        model = "M3: cov + exposure + PES",
        HR = suppressWarnings(as.numeric(HR_PES_M3_perSD)),
        pval = if (star_source == "LRT" && "p_LRT_M2_vs_M3" %in% names(.)) p_LRT_M2_vs_M3 else p_PES_M3
      ) %>%
      select(Type, pair_label, model, HR, pval)
  ) %>%
    mutate(
      model = factor(model, levels=c("M1: cov + PES", "M3: cov + exposure + PES")),
      star = p_to_star(pval)
    ) %>%
    filter(is.finite(HR), HR > 0) %>%
    mutate(logHR = log(HR))
  
  keep_pairs <- cox_hm %>%
    filter(model == "M3: cov + exposure + PES") %>%
    group_by(pair_label) %>%
    summarize(score = max(abs(logHR), na.rm=TRUE), .groups="drop") %>%
    arrange(desc(score)) %>%
    slice_head(n=topN) %>%
    pull(pair_label)
  
  cox_hm_plot <- cox_hm %>%
    filter(pair_label %in% keep_pairs) %>%
    group_by(pair_label) %>%
    summarize(ord = median(abs(logHR), na.rm=TRUE), .groups="drop") %>%
    right_join(cox_hm, by="pair_label") %>%
    filter(pair_label %in% keep_pairs) %>%
    mutate(pair_label = fct_reorder(pair_label, ord)) %>%
    select(-ord)
  
  lim <- max(abs(cox_hm_plot$logHR), na.rm=TRUE)
  hr_breaks <- c(0.67, 0.8, 1, 1.25, 1.5, 2, 3)
  log_breaks <- log(hr_breaks)
  
  p <- ggplot(cox_hm_plot, aes(x=Type, y=pair_label, fill=logHR)) +
    geom_tile(color="white", linewidth=0.25) +
    geom_text(aes(label=star), size=4, color="black") +
    facet_wrap(~model, ncol=2) +
    scale_fill_gradient2(
      low  = "#2c7bb6",
      mid  = "white",
      high = "#d7191c",
      midpoint = 0,
      limits = c(-lim, lim),
      breaks = log_breaks,
      labels = hr_breaks
    ) +
    labs(
      title="PES association with incident disease",
      subtitle=paste0(
        "HR per SD of PES; stars from ",
        ifelse(star_source=="LRT","LRT","Wald"),
        " p-values (*<0.05, **<0.01, ***<0.001). Top ", topN, " pairs by |log(HR)| in M3."
      ),
      x="Covariate specification (Type)",
      y="Exposure \u2192 disease (ICD10)",
      fill="HR"
    ) +
    theme(axis.text.y = element_text(size=9),
          strip.background = element_rect(fill="grey90", color=NA))
  
  ggsave(out_file, p, width=14, height=9, dpi=300)
  invisible(p)
}

cox_main <- filter_pes(cox_all, pes_kind_main)
cox_supp <- filter_pes(cox_all, pes_kind_supp)

make_cox_heatmap(cox_main, file.path(out_dir, "Fig4_CoxHeatmap_MAIN_PESprot.png"),
                 star_source="WALD", topN=30)
make_cox_heatmap(cox_supp, file.path(out_dir, "FigS4_CoxHeatmap_SUPP_PESfull.png"),
                 star_source="WALD", topN=30)

message("\nDONE. Wrote ONLY Fig1, Fig2, Fig3, Fig4, FigS4 to: ", out_dir, "\n")


## TABLE Stuff:

# fold_primary_full: mean/se for full model primary metric
# fold_primary_delta: mean/se for delta primary metric

metric_tbl <- fold_primary_full %>%
  transmute(
    Type, exposure_id, exposure_label, metric_name,
    full_mean = mean,
    full_l95  = mean - 1.96*se,
    full_u95  = mean + 1.96*se
  ) %>%
  inner_join(
    fold_primary_delta %>%
      transmute(
        Type, exposure_id, exposure_label, metric_name,
        delta_mean = mean,
        delta_l95  = mean - 1.96*se,
        delta_u95  = mean + 1.96*se
      ),
    by = c("Type","exposure_id","exposure_label","metric_name")
  ) %>%
  mutate(
    full_fmt  = sprintf("%.3f [%.3f, %.3f]", full_mean,  full_l95,  full_u95),
    delta_fmt = sprintf("%.3f [%.3f, %.3f]", delta_mean, delta_l95, delta_u95)
  )

metric_tbl_type5 <- metric_tbl %>%   # metric_tbl from my earlier message
  filter(Type == "Type5") %>%
  mutate(
    full_mean  = as.numeric(full_mean),
    full_l95   = as.numeric(full_l95),
    full_u95   = as.numeric(full_u95),
    delta_mean = as.numeric(delta_mean),
    delta_l95  = as.numeric(delta_l95),
    delta_u95  = as.numeric(delta_u95)
  ) %>%
  arrange(metric_name, desc(delta_mean))

suppressPackageStartupMessages({
  library(grid)
  library(gridExtra)
  library(gtable)
})

# pick columns you want to show in the figure
tabA <- metric_tbl_type5 %>%
  transmute(
    Exposure = exposure_label,
    Metric   = metric_name,
    `Full (mean [95% CI])`  = sprintf("%.3f [%.3f, %.3f]", full_mean, full_l95, full_u95),
    `Δ Full−Cov (mean [95% CI])` = sprintf("%.3f [%.3f, %.3f]", delta_mean, delta_l95, delta_u95)
  )

# optional: keep top N rows per metric for the main panel (avoids a giant table in main)
topN <- 30
tabA_show <- tabA %>%
  group_by(Metric) %>%
  slice_head(n = topN) %>%
  ungroup()

tg <- tableGrob(tabA_show, rows = NULL,
                theme = ttheme_minimal(
                  base_size = 10,
                  core = list(fg_params = list(hjust = 0, x = 0.02)),
                  colhead = list(fg_params = list(fontface = "bold", hjust = 0, x = 0.02))
                ))

# Save as a figure (PowerPoint-friendly)
png(file.path(out_dir, "Fig1A_Table_Type5.png"), width = 2200, height = 2200, res = 300)
grid.newpage(); grid.draw(tg)
dev.off()

# Also save vector PDF (journal-friendly)
pdf(file.path(out_dir, "Fig1A_Table_Type5.pdf"), width = 7.2, height = 7)
grid.newpage(); grid.draw(tg)
dev.off()




#### Another one:
library(grid)
library(gridExtra)
library(dplyr)
library(data.table)

topN_per_metric <- 10

tabA <- metric_tbl_type5 %>%
  mutate(
    Metric = metric_name,
    `Full (mean [95% CI])` = sprintf("%.3f [%.3f, %.3f]", full_mean, full_l95, full_u95),
    `Δ (Full − Cov) (mean [95% CI])` = sprintf("%.3f [%.3f, %.3f]", delta_mean, delta_l95, delta_u95)
  ) %>%
  select(
    Exposure = exposure_label,
    Metric,
    `Full (mean [95% CI])`,
    `Δ (Full − Cov) (mean [95% CI])`
  ) %>%
  group_by(Metric) %>%
  # rank within metric by delta mean (extract the first number before space)
  mutate(delta_num = as.numeric(sub("^([0-9.]+).*", "\\1", `Δ (Full − Cov) (mean [95% CI])`))) %>%
  arrange(desc(delta_num), .by_group = TRUE) %>%
  slice_head(n = topN_per_metric) %>%
  ungroup() %>%
  select(-delta_num)

tg <- tableGrob(
  tabA, rows = NULL,
  theme = ttheme_minimal(
    base_size = 11,
    core = list(fg_params = list(hjust = 0, x = 0.02)),
    colhead = list(fg_params = list(fontface = "bold", hjust = 0, x = 0.02))
  )
)

png(file.path(out_dir, "Fig1A_Table_Type5.png"), width = 2400, height = 1600, res = 300)
grid.newpage(); grid.draw(tg)
dev.off()

pdf(file.path(out_dir, "Fig1A_Table_Type5.pdf"), width = 7.5, height = 5.0)
grid.newpage(); grid.draw(tg)
dev.off()
