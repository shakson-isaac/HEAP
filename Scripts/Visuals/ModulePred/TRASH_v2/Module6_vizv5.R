#!/usr/bin/env Rscript

# ============================================================
# Single end-to-end visualization script for PES Option1 outputs
# - Reads all Types (Type1..Type5)
# - Robust to missing metrics (R2 vs AUC vs r2_code etc.)
# - Produces coherent, publication-style plots with:
#     * exposure = color
#     * disease = shape (where applicable)
#     * model = linetype (M1 vs M3 where applicable)
# - Replaces boxplots with mean±CI line plots across folds
# - Includes Cox heatmap (readable alternative to huge forest)
# - Includes Δc-index plot with exposure colors + disease shapes
# - Includes Wald p-value plot (and LRT if present)
#
# Output:
#   <base_dir>/Viz2/*.png + summary TSVs
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(forcats)
  library(ggplot2)
})

# ----------------------------
# USER CONFIG
# ----------------------------
base_dir <- "/n/groups/patel/shakson_ukb/UK_Biobank/Data/Parallel/PES_Option1_fastv2"
types <- paste0("Type", 1:5)

out_path <- "/n/groups/patel/shakson_ukb/UK_Biobank/Output/PES/"
out_dir <- file.path(out_path)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Optional: short labels (highly recommended)
# If you don’t want label mapping, leave as is.
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

clip01 <- function(p, eps=1e-15) pmin(pmax(p, eps), 1-eps)

theme_set(theme_bw(base_size = 12))

# Reads all files of a given pattern under each Type folder
read_type_files <- function(base_dir, types, pattern, add_type=TRUE) {
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
      x
    }), fill = TRUE)
    if (add_type) dt[, Type := ty]
    out[[ty]] <- dt
  }
  rbindlist(out, fill = TRUE)
}

# Make sure required columns exist; create if missing
ensure_cols <- function(df, cols) {
  miss <- setdiff(cols, names(df))
  if (length(miss) > 0) {
    for (nm in miss) df[[nm]] <- NA
  }
  df
}

# Identify which primary metric is relevant per exposure_id based on what exists:
# Prefer AUC if present for binomial models; else R2; else R2_code.
choose_primary_metric <- function(overall_long) {
  # overall_long: Type, exposure_id, model, metric_name, value
  present <- overall_long %>%
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
  present
}

# ----------------------------
# Load outputs
# ----------------------------
overall_dt <- read_type_files(base_dir, types, "OverallMetrics\\.tsv$")
fold_dt    <- read_type_files(base_dir, types, "FoldMetrics\\.tsv$")
cox_dt     <- read_type_files(base_dir, types, "^Cox4_.*\\.tsv$")

if (nrow(overall_dt) == 0) stop("No OverallMetrics.tsv found under: ", base_dir)

overall_all <- as_tibble(overall_dt) %>% mutate(Type = factor(Type, levels=types))
fold_all    <- as_tibble(fold_dt)    %>% mutate(Type = factor(Type, levels=types))
cox_all     <- as_tibble(cox_dt)     %>% mutate(Type = factor(Type, levels=types))

# ----------------------------
# PREP: Overall metrics (long)
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
    kind = case_when(
      str_starts(metric_key, "delta_") ~ "delta",
      TRUE ~ "metric"
    )
  ) %>%
  filter(!is.na(metric_name))

# Choose primary metric per exposure
primary_tbl <- choose_primary_metric(
  overall_long %>% filter(kind=="metric")
)

# Build Δ metric from cov_only and prot_plus_cov (robust)
overall_primary_delta <- overall_long %>%
  filter(kind=="metric", model %in% c("cov_only","prot_plus_cov")) %>%
  inner_join(primary_tbl %>% filter(!is.na(primary_metric)),
             by="exposure_id") %>%
  filter(metric_name == primary_metric) %>%
  select(Type, exposure_id, exposure_type, metric_name, model, value) %>%
  group_by(Type, exposure_id, exposure_type, metric_name, model) %>%
  summarize(value=value[1], .groups="drop") %>%
  pivot_wider(names_from=model, values_from=value) %>%
  mutate(
    delta_metric = prot_plus_cov - cov_only
  ) %>%
  filter(is.finite(delta_metric))

# Add labels
overall_primary_delta <- overall_primary_delta %>%
  left_join(exposure_label_map, by="exposure_id") %>%
  mutate(exposure_label = ifelse(is.na(exposure_label), exposure_id, exposure_label),
         exposure_label = factor(exposure_label))

# ----------------------------
# PREP: Fold metrics for full model (primary metric)
# ----------------------------
# Build fold-long for possible metrics
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

# keep only full model and only primary metric per exposure
fold_primary_full <- fold_long %>%
  filter(model=="prot_plus_cov") %>%
  inner_join(primary_tbl %>% filter(!is.na(primary_metric)), by="exposure_id") %>%
  filter(metric_name == primary_metric) %>%
  group_by(Type, exposure_id, metric_name) %>%
  summarize(
    mean = mean(value, na.rm=TRUE),
    sd   = sd(value, na.rm=TRUE),
    n    = sum(is.finite(value)),
    se   = sd / sqrt(pmax(n,1)),
    .groups="drop"
  ) %>%
  filter(n >= 2) %>%
  left_join(exposure_label_map, by="exposure_id") %>%
  mutate(exposure_label = ifelse(is.na(exposure_label), exposure_id, exposure_label),
         exposure_label = factor(exposure_label))

# Selected proteins across folds (full model)
sel_primary <- NULL
if ("n_proteins_selected_full" %in% names(fold_all)) {
  sel_primary <- fold_all %>%
    mutate(
      exposure_id = as.character(exposure_id),
      n_proteins_selected_full = suppressWarnings(as.numeric(n_proteins_selected_full))
    ) %>%
    filter(is.finite(n_proteins_selected_full)) %>%
    group_by(Type, exposure_id) %>%
    summarize(
      mean = mean(n_proteins_selected_full),
      sd   = sd(n_proteins_selected_full),
      n    = n(),
      se   = sd / sqrt(pmax(n,1)),
      .groups="drop"
    ) %>%
    filter(n >= 2) %>%
    left_join(exposure_label_map, by="exposure_id") %>%
    mutate(exposure_label = ifelse(is.na(exposure_label), exposure_id, exposure_label),
           exposure_label = factor(exposure_label))
}

# ----------------------------
# PREP: Cox metrics for readable plots
# ----------------------------
cox_need <- c(
  "Type","disease_age_col","exposure_id","exposure_type","n","events",
  "HR_PES_M1_perSD","HR_PES_M1_L95","HR_PES_M1_U95","p_PES_M1",
  "HR_PES_M3_perSD","HR_PES_M3_L95","HR_PES_M3_U95","p_PES_M3",
  "p_LRT_M0_vs_M1","p_LRT_M2_vs_M3",
  "cindex_M0","cindex_M1","cindex_M2","cindex_M3",
  "delta_cindex_M0_to_M1","delta_cindex_M2_to_M3",
  "pes_used"
)

cox_sum <- as_tibble(cox_all)
if (nrow(cox_sum) > 0) {
  cox_sum <- ensure_cols(cox_sum, cox_need) %>%
    mutate(
      exposure_id = as.character(exposure_id),
      disease_age_col = as.character(disease_age_col)
    ) %>%
    left_join(exposure_label_map, by="exposure_id") %>%
    left_join(disease_label_map, by="disease_age_col") %>%
    mutate(
      exposure_label = ifelse(is.na(exposure_label), exposure_id, exposure_label),
      disease_label  = ifelse(is.na(disease_label), disease_age_col, disease_label),
      pair_label = paste0(exposure_label, " → ", disease_label)
    )
}

# ----------------------------
# PLOT 1: Δ metric (primary) vs Type, colored by exposure
# ----------------------------
p1 <- ggplot(overall_primary_delta,
             aes(x=Type, y=delta_metric, color=exposure_label, group=exposure_label)) +
  geom_hline(yintercept=0, linetype=2) +
  geom_line(alpha=0.7, linewidth=0.7) +
  geom_point(size=2) +
  facet_wrap(~metric_name, scales="free_y") +
  labs(
    title="Incremental predictive value of proteins beyond covariates",
    subtitle="Δ metric = (prot_plus_cov) − (cov_only); metric chosen per exposure (AUC > R2 > R2_code)",
    x="Covariate specification (Type)",
    y="Δ metric",
    color="Exposure"
  ) +
  theme(legend.position="right")

ggsave(file.path(out_dir, "Fig1_DeltaPrimaryMetric_lines.png"), p1, width=13, height=6, dpi=300)

# ----------------------------
# PLOT 2: Fold stability (full model, primary metric), mean±95% CI, colored by exposure
# ----------------------------
p2 <- ggplot(fold_primary_full,
             aes(x=Type, y=mean, color=exposure_label, group=exposure_label)) +
  geom_line(linewidth=0.7) +
  geom_point(size=2) +
  geom_errorbar(aes(ymin=mean-1.96*se, ymax=mean+1.96*se),
                width=0.12, alpha=0.75) +
  facet_wrap(~metric_name, scales="free_y") +
  labs(
    title="Fold-to-fold stability of exposure prediction (full model)",
    subtitle="Points/lines = mean across folds; error bars = 95% CI across folds",
    x="Covariate specification (Type)",
    y="Fold metric (mean ± 95% CI)",
    color="Exposure"
  ) +
  theme(legend.position="right")

ggsave(file.path(out_dir, "Fig2_FoldStability_primaryMetric_lines.png"), p2, width=13, height=6, dpi=300)

# ----------------------------
# PLOT 3: Selected proteins (full model), mean±95% CI, colored by exposure
# ----------------------------
if (!is.null(sel_primary) && nrow(sel_primary) > 0) {
  p3 <- ggplot(sel_primary,
               aes(x=Type, y=mean, color=exposure_label, group=exposure_label)) +
    geom_line(linewidth=0.7) +
    geom_point(size=2) +
    geom_errorbar(aes(ymin=mean-1.96*se, ymax=mean+1.96*se),
                  width=0.12, alpha=0.75) +
    labs(
      title="Number of selected proteins across folds (full model)",
      subtitle="Mean ± 95% CI across folds",
      x="Covariate specification (Type)",
      y="# proteins selected",
      color="Exposure"
    ) +
    theme(legend.position="right")
  
  ggsave(file.path(out_dir, "Fig3_SelectedProteins_lines.png"), p3, width=13, height=6, dpi=300)
}

# ----------------------------
# PLOT 4: Cox heatmap (readable) – log(HR) per Type
#         Two panels: M1 and M3
# ----------------------------
# if (nrow(cox_sum) > 0) {
#   cox_hm <- bind_rows(
#     cox_sum %>%
#       mutate(model="M1: cov + PES",
#              HR = suppressWarnings(as.numeric(HR_PES_M1_perSD))) %>%
#       select(Type, pair_label, model, HR),
#     cox_sum %>%
#       mutate(model="M3: cov + exposure + PES",
#              HR = suppressWarnings(as.numeric(HR_PES_M3_perSD))) %>%
#       select(Type, pair_label, model, HR)
#   ) %>%
#     mutate(
#       logHR = log(HR),
#       model = factor(model, levels=c("M1: cov + PES", "M3: cov + exposure + PES"))
#     ) %>%
#     filter(is.finite(logHR))
#   
#   # If too many pairs, show top N by abs(logHR) in M3 across Types (optional)
#   # Comment out if you want all.
#   topN <- 30
#   keep_pairs <- cox_hm %>%
#     filter(model=="M3: cov + exposure + PES") %>%
#     group_by(pair_label) %>%
#     summarize(score = max(abs(logHR), na.rm=TRUE), .groups="drop") %>%
#     arrange(desc(score)) %>%
#     slice_head(n=topN) %>%
#     pull(pair_label)
#   
#   cox_hm_plot <- cox_hm %>%
#     filter(pair_label %in% keep_pairs) %>%
#     mutate(pair_label = fct_reorder(pair_label, logHR, .fun = median))
#   
#   p4 <- ggplot(cox_hm_plot, aes(x=Type, y=pair_label, fill=logHR)) +
#     geom_tile() +
#     facet_wrap(~model, ncol=2) +
#     labs(
#       title="PES association with incident disease (readable summary)",
#       subtitle=paste0("Heatmap of log(HR per SD of PES); showing top ", topN, " exposure→disease pairs by |log(HR)| in M3"),
#       x="Covariate specification (Type)",
#       y="Exposure → disease",
#       fill="log(HR)"
#     ) +
#     theme(axis.text.y = element_text(size=8))
#   
#   ggsave(file.path(out_dir, "Fig4_CoxHeatmap_topPairs.png"), p4, width=14, height=9, dpi=300)
# }

# ----------------------------
# PLOT 4 (FIXED): Heatmap where white is EXACTLY HR=1
#   - Fill is log(HR), so 0 == HR 1 (guaranteed white at HR=1)
#   - Legend labels shown in HR units
# ----------------------------
if (nrow(cox_sum) > 0) {
  
  star_source <- "WALD"  # or "LRT" if you have p_LRT_* columns
  
  p_to_star <- function(p) {
    p <- suppressWarnings(as.numeric(p))
    ifelse(is.na(p), "",
           ifelse(p < 1e-3, "***",
                  ifelse(p < 1e-2, "**",
                         ifelse(p < 5e-2, "*", ""))))
  }
  
  cox_hm <- bind_rows(
    cox_sum %>%
      mutate(
        model = "M1: cov + PES",
        HR = suppressWarnings(as.numeric(HR_PES_M1_perSD)),
        pval = if (star_source == "LRT" && "p_LRT_M0_vs_M1" %in% names(.)) p_LRT_M0_vs_M1 else p_PES_M1
      ) %>%
      select(Type, pair_label, model, HR, pval),
    cox_sum %>%
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
    mutate(logHR = log(HR))  # <-- key: 0 means HR=1
  
  # Top N pairs by strongest deviation from HR=1 in M3
  topN <- 30
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
  
  # symmetric limits around 0 to keep palette balanced
  lim <- max(abs(cox_hm_plot$logHR), na.rm=TRUE)
  
  # HR breaks for legend (displayed in HR units)
  hr_breaks <- c(0.67, 0.8, 1, 1.25, 1.5, 2, 3)
  log_breaks <- log(hr_breaks)
  
  p4 <- ggplot(cox_hm_plot, aes(x=Type, y=pair_label, fill=logHR)) +
    geom_tile(color="white", linewidth=0.25) +
    geom_text(aes(label=star), size=4, color="black") +
    facet_wrap(~model, ncol=2) +
    scale_fill_gradient2(
      low  = "#2c7bb6",
      mid  = "white",
      high = "#d7191c",
      midpoint = 0,              # <-- EXACTLY log(HR)=0 => HR=1 is white
      limits = c(-lim, lim),     # symmetric so the midpoint is visually centered
      breaks = log_breaks,
      labels = hr_breaks
    ) +
    labs(
      title="PES association with incident disease",
      subtitle=paste0(
        "Heatmap of HR per SD of PES (stars: ",
        ifelse(star_source=="LRT","LRT p-values","Wald p-values"),
        "; * <0.05, ** <0.01, *** <0.001). Top ", topN, " pairs by |log(HR)| in M3."
      ),
      x="Covariate specification (Type)",
      y="Exposure → disease",
      fill="HR"
    ) +
    theme(
      axis.text.y = element_text(size=9),
      strip.background = element_rect(fill="grey90", color=NA)
    )
  
  ggsave(file.path(out_dir, "Fig4_CoxHeatmap_HR_withStars_whiteAt1.png"),
         p4, width=14, height=9, dpi=300)
}










# ----------------------------
# PLOT 5: Δc-index (M0→M1 and M2→M3) with exposure colors + disease shapes
# ----------------------------
if (nrow(cox_sum) > 0) {
  cidx_long <- cox_sum %>%
    select(Type, exposure_label, disease_label,
           delta_cindex_M0_to_M1, delta_cindex_M2_to_M3) %>%
    pivot_longer(cols=c(delta_cindex_M0_to_M1, delta_cindex_M2_to_M3),
                 names_to="delta_kind", values_to="delta_cindex") %>%
    mutate(
      delta_cindex = suppressWarnings(as.numeric(delta_cindex)),
      delta_kind = recode(delta_kind,
                          delta_cindex_M0_to_M1 = "Δc-index: M0→M1 (add PES)",
                          delta_cindex_M2_to_M3 = "Δc-index: M2→M3 (add PES to exposure model)")
    ) %>%
    filter(is.finite(delta_cindex))
  
  p5 <- ggplot(cidx_long, aes(x=Type, y=delta_cindex,
                              color=exposure_label, shape=disease_label,
                              group=interaction(exposure_label, disease_label))) +
    geom_hline(yintercept=0, linetype=2) +
    geom_point(position=position_jitter(width=0.12, height=0), alpha=0.85, size=2) +
    stat_summary(aes(group=exposure_label),
                 fun=mean, geom="line", linewidth=0.7, alpha=0.6) +
    facet_wrap(~delta_kind, ncol=1, scales="free_y") +
    labs(
      title="Incremental discrimination from PES",
      subtitle="Points = exposure–disease pairs; line = mean across diseases for each exposure",
      x="Covariate specification (Type)",
      y="Δc-index",
      color="Exposure",
      shape="Disease"
    ) +
    theme(legend.position="right")
  
  ggsave(file.path(out_dir, "Fig5_DeltaCindex_pointsLines.png"), p5, width=14, height=10, dpi=300)
}

# ----------------------------
# PLOT 6: Significance evidence
#   - Wald p-values for PES term (always present if p_PES_* exists)
#   - LRT p-values if p_LRT_* exists in your cox files
# ----------------------------
if (nrow(cox_sum) > 0) {
  
  # Wald (PES term)
  wald_long <- cox_sum %>%
    select(Type, exposure_label, disease_label, any_of(c("p_PES_M1","p_PES_M3"))) %>%
    pivot_longer(cols=any_of(c("p_PES_M1","p_PES_M3")),
                 names_to="test", values_to="pval") %>%
    mutate(
      pval = suppressWarnings(as.numeric(pval)),
      pval = ifelse(is.na(pval) | pval <= 0, NA_real_, pval),
      mlog10p = -log10(pval),
      test = recode(test,
                    p_PES_M1="Wald p: PES in M1 (cov + PES)",
                    p_PES_M3="Wald p: PES in M3 (cov + exposure + PES)")
    ) %>%
    filter(is.finite(mlog10p))
  
  # LRT (nested model comparisons) if present
  lrt_long <- cox_sum %>%
    select(Type, exposure_label, disease_label, any_of(c("p_LRT_M0_vs_M1","p_LRT_M2_vs_M3"))) %>%
    pivot_longer(cols=any_of(c("p_LRT_M0_vs_M1","p_LRT_M2_vs_M3")),
                 names_to="test", values_to="pval") %>%
    mutate(
      pval = suppressWarnings(as.numeric(pval)),
      pval = ifelse(is.na(pval) | pval <= 0, NA_real_, pval),
      mlog10p = -log10(pval),
      test = recode(test,
                    p_LRT_M0_vs_M1="LRT: M0 vs M1 (add PES)",
                    p_LRT_M2_vs_M3="LRT: M2 vs M3 (add PES to exposure model)")
    ) %>%
    filter(is.finite(mlog10p))
  
  sig_long <- bind_rows(
    wald_long %>% mutate(test_family="Wald"),
    lrt_long  %>% mutate(test_family="LRT")
  )
  
  if (nrow(sig_long) > 0) {
    p6 <- ggplot(sig_long,
                 aes(x=Type, y=mlog10p,
                     color=exposure_label, shape=disease_label)) +
      geom_hline(yintercept=-log10(0.05), linetype=2) +
      geom_point(position=position_jitter(width=0.12, height=0), alpha=0.85, size=2) +
      facet_grid(test_family ~ test, scales="free_y") +
      labs(
        title="Statistical evidence for PES contribution",
        subtitle="Points are exposure–disease pairs; dashed line is p=0.05",
        x="Covariate specification (Type)",
        y=expression(-log[10](p)),
        color="Exposure",
        shape="Disease"
      ) +
      theme(legend.position="right")
    
    ggsave(file.path(out_dir, "Fig6_Significance_Wald_and_LRT.png"), p6, width=16, height=9, dpi=300)
  }
}

# ----------------------------
# Save summary tables used for plotting
# ----------------------------
fwrite(as.data.table(overall_primary_delta), file.path(out_dir, "overall_primary_delta.tsv"), sep="\t")
fwrite(as.data.table(fold_primary_full), file.path(out_dir, "fold_primary_full.tsv"), sep="\t")
if (!is.null(sel_primary)) fwrite(as.data.table(sel_primary), file.path(out_dir, "selected_proteins_summary.tsv"), sep="\t")
if (nrow(cox_sum) > 0) fwrite(as.data.table(cox_sum), file.path(out_dir, "cox_all.tsv"), sep="\t")

message("\nDONE. Wrote plots + tables to: ", out_dir, "\n")
