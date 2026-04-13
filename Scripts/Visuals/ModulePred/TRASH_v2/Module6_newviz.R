#!/usr/bin/env Rscript

# ============================================================
# End-to-end visualization script for UPDATED PES outputs
# - Reads all Types (Type1..Type5)
# - Works with per-exposure output filenames:
#     PES_<Type>_<exposure>_OverallMetrics.tsv
#     PES_<Type>_<exposure>_FoldMetrics.tsv
#     Cox4All_<Type>_<exposure>__PESprot.tsv
#     Cox4All_<Type>_<exposure>__PESfull.tsv
#     PES_<Type>_<exposure>_WeightsSummary_ProtOnly.tsv
#     PES_<Type>_<exposure>_WeightsSummary_Full_ProtOnly.tsv
#
# Main vs Supplement:
#   - Main Cox figures use PESprot
#   - Supplement Cox figures use PESfull
#
# Adds:
#   - Protein signature plots: top 20 proteins per exposure with weights
#   - UMAP of proteomic signatures across exposures (Option C):
#       signature_p = selected_frac * mean_beta_among_selected
#     Uses cosine distance (uwot::umap metric="cosine")
#     One embedding across all Types; draws Type trajectories per exposure
#
# Output:
#   <out_dir>/*.png + summary TSVs + ProteinSignatures/ + UMAP/
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(forcats)
  library(ggplot2)
  library(ggrepel)
  library(tibble)
  library(purrr)
  library(uwot)
})

# Optional but recommended packages
has_uwot    <- requireNamespace("uwot", quietly = TRUE)
has_ggrepel <- requireNamespace("ggrepel", quietly = TRUE)

# ----------------------------
# USER CONFIG
# ----------------------------
base_dir <- "/n/groups/patel/shakson_ukb/UK_Biobank/Data/Parallel/PES_test"
types <- paste0("Type", 1:5)

out_path <- "/n/groups/patel/shakson_ukb/UK_Biobank/Output/PES/"
out_dir <- file.path(out_path)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Main vs Supp for Cox-based plots
pes_kind_main <- "PESprot"
pes_kind_supp <- "PESfull"

# Optional: short labels
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

# Make sure required columns exist; create if missing
ensure_cols <- function(df, cols) {
  miss <- setdiff(cols, names(df))
  if (length(miss) > 0) {
    for (nm in miss) df[[nm]] <- NA
  }
  df
}

# Identify which primary metric is relevant per exposure_id
choose_primary_metric <- function(overall_long) {
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

filter_pes <- function(df, pes_kind) {
  if (!"pes_kind" %in% names(df)) return(df)
  df %>% filter(is.na(pes_kind) | pes_kind == !!pes_kind)
}

`%||%` <- function(x, y) if (!is.null(x) && length(x) > 0) x else y

# Parse exposure_id from weights summary filenames
parse_exposure_from_weights_file <- function(file, type) {
  x <- str_replace(file, paste0("^PES_", type, "_"), "")
  x <- str_replace(x, "_WeightsSummary_ProtOnly\\.tsv$", "")
  x <- str_replace(x, "_WeightsSummary_Full_ProtOnly\\.tsv$", "")
  x
}

add_exposure_id_from_weights_filename <- function(df) {
  df %>%
    mutate(
      Type_chr = as.character(Type),
      exposure_id = map2_chr(file, Type_chr, parse_exposure_from_weights_file)
    ) %>%
    select(-Type_chr)
}

# Pretty stars for p-values
p_to_star <- function(p) {
  p <- suppressWarnings(as.numeric(p))
  ifelse(is.na(p), "",
         ifelse(p < 1e-3, "***",
                ifelse(p < 1e-2, "**",
                       ifelse(p < 5e-2, "*", ""))))
}

# ----------------------------
# Load outputs
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

w_prot_sum_dt <- read_type_files(base_dir, types, "_WeightsSummary_ProtOnly\\.tsv$")
w_full_sum_dt <- read_type_files(base_dir, types, "_WeightsSummary_Full_ProtOnly\\.tsv$")

if (nrow(overall_dt) == 0) stop("No *_OverallMetrics.tsv found under: ", base_dir)

overall_all <- as_tibble(overall_dt) %>% mutate(Type = factor(Type, levels=types))
fold_all    <- as_tibble(fold_dt)    %>% mutate(Type = factor(Type, levels=types))
cox_all     <- as_tibble(cox_dt)     %>% mutate(Type = factor(Type, levels=types))

w_prot_sum <- as_tibble(w_prot_sum_dt) %>% mutate(Type = factor(Type, levels=types))
w_full_sum <- as_tibble(w_full_sum_dt) %>% mutate(Type = factor(Type, levels=types))

# Add exposure_id to weights summaries (from filename)
if (nrow(w_prot_sum) > 0) w_prot_sum <- add_exposure_id_from_weights_filename(w_prot_sum)
if (nrow(w_full_sum) > 0) w_full_sum <- add_exposure_id_from_weights_filename(w_full_sum)

# Add exposure labels to everything that has exposure_id
overall_all <- overall_all %>%
  left_join(exposure_label_map, by="exposure_id") %>%
  mutate(exposure_label = ifelse(is.na(exposure_label), exposure_id, exposure_label))

fold_all <- fold_all %>%
  left_join(exposure_label_map, by="exposure_id") %>%
  mutate(exposure_label = ifelse(is.na(exposure_label), exposure_id, exposure_label))

if (nrow(cox_all) > 0) {
  cox_all <- cox_all %>%
    left_join(exposure_label_map, by="exposure_id") %>%
    left_join(disease_label_map, by="disease_age_col") %>%
    mutate(
      exposure_label = ifelse(is.na(exposure_label), exposure_id, exposure_label),
      disease_label  = ifelse(is.na(disease_label), disease_age_col, disease_label),
      pair_label = paste0(exposure_label, " \u2192 ", disease_label)
    )
}

if (nrow(w_prot_sum) > 0) {
  w_prot_sum <- w_prot_sum %>%
    left_join(exposure_label_map, by="exposure_id") %>%
    mutate(exposure_label = ifelse(is.na(exposure_label), exposure_id, exposure_label))
}
if (nrow(w_full_sum) > 0) {
  w_full_sum <- w_full_sum %>%
    left_join(exposure_label_map, by="exposure_id") %>%
    mutate(exposure_label = ifelse(is.na(exposure_label), exposure_id, exposure_label))
}

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

primary_tbl <- choose_primary_metric(overall_long %>% filter(kind=="metric"))

overall_primary_delta <- overall_long %>%
  filter(kind=="metric", model %in% c("cov_only","prot_plus_cov")) %>%
  inner_join(primary_tbl %>% filter(!is.na(primary_metric)),
             by="exposure_id") %>%
  filter(metric_name == primary_metric) %>%
  select(Type, exposure_id, exposure_type, exposure_label, metric_name, model, value) %>%
  group_by(Type, exposure_id, exposure_type, exposure_label, metric_name, model) %>%
  summarize(value=value[1], .groups="drop") %>%
  pivot_wider(names_from=model, values_from=value) %>%
  mutate(delta_metric = prot_plus_cov - cov_only) %>%
  filter(is.finite(delta_metric)) %>%
  mutate(exposure_label = factor(exposure_label))

# ----------------------------
# PREP: Fold metrics for full model (primary metric)
# ----------------------------
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

sel_primary <- NULL
if ("n_proteins_selected_full" %in% names(fold_all)) {
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
}

# ----------------------------
# PREP: Cox summary
# ----------------------------
cox_need <- c(
  "Type","pes_kind","cox_status",
  "disease_age_col","exposure_id","exposure_type","n","events",
  "HR_PES_M1_perSD","HR_PES_M1_L95","HR_PES_M1_U95","p_PES_M1",
  "HR_PES_M3_perSD","HR_PES_M3_L95","HR_PES_M3_U95","p_PES_M3",
  "p_LRT_M0_vs_M1","p_LRT_M2_vs_M3",
  "cindex_M0","cindex_M1","cindex_M2","cindex_M3",
  "delta_cindex_M0_to_M1","delta_cindex_M2_to_M3",
  "pes_used","pair_label","exposure_label","disease_label"
)

cox_sum <- as_tibble(cox_all)
if (nrow(cox_sum) > 0) {
  cox_sum <- ensure_cols(cox_sum, cox_need) %>%
    mutate(
      exposure_id = as.character(exposure_id),
      disease_age_col = as.character(disease_age_col),
      pes_kind = as.character(pes_kind),
      cox_status = as.character(cox_status),
      exposure_label = as.character(exposure_label),
      disease_label = as.character(disease_label),
      pair_label = as.character(pair_label)
    )
  # Keep only successful fits if you want:
  # cox_sum <- cox_sum %>% filter(cox_status == "OK")
}

# ----------------------------
# PLOT 1: Δ metric vs Type (exposure colored)
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
# PLOT 2: Fold stability (full model, primary metric), mean±95% CI
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
# PLOT 3: # selected proteins, mean±95% CI
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

# ============================================================
# Cox plots as functions (run twice: MAIN PESprot, SUPP PESfull)
# ============================================================

make_cox_heatmap <- function(cox_df, out_file, star_source="WALD", topN=30) {
  if (nrow(cox_df) == 0) return(NULL)
  
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
  
  p4 <- ggplot(cox_hm_plot, aes(x=Type, y=pair_label, fill=logHR)) +
    geom_tile(color="white", linewidth=0.25) +
    geom_text(aes(label=star), size=4, color="black") +
    facet_wrap(~model, ncol=2) +
    scale_fill_gradient2(
      low  = "#2c7bb6",
      mid  = "white",
      high = "#d7191c",
      midpoint = 0,              # log(HR)=0 => HR=1 is white
      limits = c(-lim, lim),
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
      y="Exposure \u2192 disease",
      fill="HR"
    ) +
    theme(
      axis.text.y = element_text(size=9),
      strip.background = element_rect(fill="grey90", color=NA)
    )
  
  ggsave(out_file, p4, width=14, height=9, dpi=300)
  p4
}

make_delta_cindex <- function(cox_df, out_file) {
  if (nrow(cox_df) == 0) return(NULL)
  
  cidx_long <- cox_df %>%
    select(Type, exposure_label, disease_label,
           delta_cindex_M0_to_M1, delta_cindex_M2_to_M3) %>%
    pivot_longer(cols=c(delta_cindex_M0_to_M1, delta_cindex_M2_to_M3),
                 names_to="delta_kind", values_to="delta_cindex") %>%
    mutate(
      delta_cindex = suppressWarnings(as.numeric(delta_cindex)),
      delta_kind = recode(delta_kind,
                          delta_cindex_M0_to_M1 = "\u0394c-index: M0\u2192M1 (add PES)",
                          delta_cindex_M2_to_M3 = "\u0394c-index: M2\u2192M3 (add PES to exposure model)")
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
      y="\u0394c-index",
      color="Exposure",
      shape="Disease"
    ) +
    theme(legend.position="right")
  
  ggsave(out_file, p5, width=14, height=10, dpi=300)
  p5
}

make_significance_plot <- function(cox_df, out_file) {
  if (nrow(cox_df) == 0) return(NULL)
  
  wald_long <- cox_df %>%
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
  
  lrt_long <- cox_df %>%
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
  
  if (nrow(sig_long) == 0) return(NULL)
  
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
  
  ggsave(out_file, p6, width=16, height=9, dpi=300)
  p6
}

# Run Cox plots twice: MAIN (PESprot) and SUPP (PESfull)
if (nrow(cox_sum) > 0) {
  cox_main <- filter_pes(cox_sum, pes_kind_main)
  cox_supp <- filter_pes(cox_sum, pes_kind_supp)
  
  if (nrow(cox_main) > 0) {
    make_cox_heatmap(cox_main, file.path(out_dir, "Fig4_CoxHeatmap_MAIN_PESprot.png"), star_source="WALD", topN=30)
    make_delta_cindex(cox_main, file.path(out_dir, "Fig5_DeltaCindex_MAIN_PESprot.png"))
    make_significance_plot(cox_main, file.path(out_dir, "Fig6_Significance_MAIN_PESprot.png"))
  }
  
  if (nrow(cox_supp) > 0) {
    make_cox_heatmap(cox_supp, file.path(out_dir, "FigS4_CoxHeatmap_SUPP_PESfull.png"), star_source="WALD", topN=30)
    make_delta_cindex(cox_supp, file.path(out_dir, "FigS5_DeltaCindex_SUPP_PESfull.png"))
    make_significance_plot(cox_supp, file.path(out_dir, "FigS6_Significance_SUPP_PESfull.png"))
  }
}

# ============================================================
# Protein signature plots: Top 20 proteins per exposure
# ============================================================
sig_dir <- file.path(out_dir, "ProteinSignatures")
dir.create(sig_dir, showWarnings = FALSE, recursive = TRUE)

plot_top_proteins <- function(w_df, exposure_id_one, out_file, topN=20, title_prefix="") {
  if (nrow(w_df) == 0) return(NULL)
  
  dd <- w_df %>%
    filter(exposure_id == exposure_id_one) %>%
    mutate(
      mean_beta_including_zeros = as.numeric(mean_beta_including_zeros),
      selected_frac = as.numeric(selected_frac)
    ) %>%
    filter(is.finite(mean_beta_including_zeros)) %>%
    mutate(abs_w = abs(mean_beta_including_zeros)) %>%
    arrange(desc(abs_w)) %>%
    slice_head(n = topN) %>%
    mutate(protein = fct_reorder(protein, mean_beta_including_zeros))
  
  if (nrow(dd) == 0) return(NULL)
  
  exp_lab <- unique(dd$exposure_label)[1] %||% exposure_id_one
  
  p <- ggplot(dd, aes(x = protein, y = mean_beta_including_zeros)) +
    geom_col() +
    coord_flip() +
    labs(
      title = paste0(title_prefix, exp_lab, ": top ", topN, " proteins by |mean weight|"),
      subtitle = "Weights are mean_beta_including_zeros (Lasso coefficients averaged across folds; zeros included).",
      x = NULL,
      y = "Mean weight (signed)"
    ) +
    theme_bw(base_size = 12)
  
  ggsave(out_file, p, width=9, height=7, dpi=300)
  p
}

# Make signature plots for BOTH weights summaries (Full and ProtOnly)
if (nrow(w_full_sum) > 0) {
  all_exposures_full <- sort(unique(w_full_sum$exposure_id))
  for (eid in all_exposures_full) {
    plot_top_proteins(
      w_df = w_full_sum,
      exposure_id_one = eid,
      out_file = file.path(sig_dir, paste0("Top20Proteins_FULL_", eid, ".png")),
      topN = 20,
      title_prefix = "Signature (Full model, protein terms): "
    )
  }
  
  top_tbl_full <- w_full_sum %>%
    group_by(Type, exposure_id, exposure_label) %>%
    mutate(abs_w = abs(as.numeric(mean_beta_including_zeros))) %>%
    arrange(desc(abs_w), .by_group = TRUE) %>%
    slice_head(n = 20) %>%
    ungroup()
  
  fwrite(as.data.table(top_tbl_full), file.path(sig_dir, "Top20Proteins_perExposure_FULL.tsv"), sep="\t")
}

if (nrow(w_prot_sum) > 0) {
  all_exposures_prot <- sort(unique(w_prot_sum$exposure_id))
  for (eid in all_exposures_prot) {
    plot_top_proteins(
      w_df = w_prot_sum,
      exposure_id_one = eid,
      out_file = file.path(sig_dir, paste0("Top20Proteins_PROT_", eid, ".png")),
      topN = 20,
      title_prefix = "Signature (Protein-only model): "
    )
  }
  
  top_tbl_prot <- w_prot_sum %>%
    group_by(Type, exposure_id, exposure_label) %>%
    mutate(abs_w = abs(as.numeric(mean_beta_including_zeros))) %>%
    arrange(desc(abs_w), .by_group = TRUE) %>%
    slice_head(n = 20) %>%
    ungroup()
  
  fwrite(as.data.table(top_tbl_prot), file.path(sig_dir, "Top20Proteins_perExposure_PROT.tsv"), sep="\t")
}

# ============================================================
# UMAP of proteomic signatures across exposures (Option C)
# signature_p = selected_frac * mean_beta_among_selected
# - Uses weights_summary_full by default (proteins from full model)
# - One embedding across all Types; shows Type trajectories per exposure
# ============================================================
umap_dir <- file.path(out_dir, "UMAP")
dir.create(umap_dir, showWarnings = FALSE, recursive = TRUE)

make_signature_matrix <- function(w_df, signature_mode=c("C_signed_stability","A_mean_beta_including_zeros","B_selected_frac")) {
  signature_mode <- match.arg(signature_mode)
  
  dd <- w_df %>%
    mutate(
      selected_frac = as.numeric(selected_frac),
      mean_beta_among_selected = as.numeric(mean_beta_among_selected),
      mean_beta_including_zeros = as.numeric(mean_beta_including_zeros)
    ) %>%
    filter(!is.na(exposure_id), !is.na(protein)) %>%
    filter(is.finite(selected_frac) | is.finite(mean_beta_among_selected) | is.finite(mean_beta_including_zeros))
  
  if (nrow(dd) == 0) return(NULL)
  
  if (signature_mode == "C_signed_stability") {
    dd <- dd %>% mutate(sig_value = selected_frac * mean_beta_among_selected)
  } else if (signature_mode == "A_mean_beta_including_zeros") {
    dd <- dd %>% mutate(sig_value = mean_beta_including_zeros)
  } else {
    dd <- dd %>% mutate(sig_value = selected_frac)
  }
  
  dd <- dd %>% mutate(sig_value = ifelse(is.finite(sig_value), sig_value, 0))
  
  # row id
  key_df <- dd %>%
    mutate(row_id = paste0(as.character(Type), "||", exposure_id)) %>%
    select(row_id, Type, exposure_id, exposure_label, protein, sig_value)
  
  # Add exposure_type if present; otherwise fill NA
  if ("exposure_type" %in% names(w_df)) {
    key_df <- key_df %>%
      left_join(w_df %>% select(Type, exposure_id, exposure_type) %>% distinct(),
                by = c("Type","exposure_id"))
  } else {
    key_df$exposure_type <- NA_character_
  }
  
  # wide matrix: rows = Type||exposure_id, cols = proteins
  mat <- key_df %>%
    select(row_id, protein, sig_value) %>%
    distinct() %>%
    pivot_wider(names_from = protein, values_from = sig_value, values_fill = 0)
  
  rowmeta <- key_df %>%
    select(row_id, Type, exposure_id, exposure_label, exposure_type) %>%
    distinct()
  
  X <- mat %>% select(-row_id) %>% as.matrix()
  rownames(X) <- mat$row_id
  
  list(X = X, meta = rowmeta)
}

l2_normalize_rows <- function(X) {
  nrm <- sqrt(rowSums(X^2))
  nrm[nrm == 0 | !is.finite(nrm)] <- 1
  X / nrm
}

plot_umap <- function(emb_df, out_file, label_top_n = 20, title="UMAP of proteomic signatures", subtitle=NULL) {
  p <- ggplot(emb_df, aes(x = UMAP1, y = UMAP2)) +
    geom_path(aes(group = exposure_id), alpha=0.25) +
    geom_point(aes(color = exposure_label, shape = Type), size=2, alpha=0.9) +
    labs(title = title, subtitle = subtitle, color = "Exposure", shape = "Type") +
    theme(legend.position = "right")
  
  # label a subset (largest signature norms)
  if (label_top_n > 0 && has_ggrepel) {
    lab <- emb_df %>%
      group_by(exposure_id) %>%
      summarize(exposure_label = first(exposure_label), max_norm = max(sig_norm, na.rm=TRUE), .groups="drop") %>%
      arrange(desc(max_norm)) %>%
      slice_head(n = label_top_n)
    
    lab_pts <- emb_df %>%
      group_by(exposure_id) %>%
      filter(sig_norm == max(sig_norm, na.rm=TRUE)) %>%
      ungroup() %>%
      inner_join(lab, by = c("exposure_id","exposure_label"))
    
    p <- p + ggrepel::geom_text_repel(
      data = lab_pts,
      aes(label = exposure_label),
      size = 3,
      max.overlaps = Inf,
      box.padding = 0.3,
      point.padding = 0.2
    )
  }
  
  ggsave(out_file, p, width=14, height=9, dpi=300)
  p
}

if (!has_uwot) {
  message("NOTE: Package 'uwot' not installed. Skipping UMAP. Install with: install.packages('uwot')")
} else if (nrow(w_full_sum) > 0) {
  
  # Use FULL weights summaries for signature space
  sig_mode <- "C_signed_stability"  # recommended
  
  sig_obj <- make_signature_matrix(w_full_sum, signature_mode = sig_mode)
  if (!is.null(sig_obj)) {
    
    X <- sig_obj$X
    meta <- sig_obj$meta
    
    # signature norm before normalization (for labeling)
    sig_norm <- sqrt(rowSums(X^2))
    sig_norm[!is.finite(sig_norm)] <- NA_real_
    
    # L2 normalize => cosine distance focuses on direction
    Xn <- l2_normalize_rows(X)
    
    set.seed(123)
    emb <- uwot::umap(
      Xn,
      n_neighbors = 25,
      min_dist = 0.15,
      metric = "cosine",
      n_components = 2,
      verbose = TRUE
    )
    
    emb_df <- meta %>%
      mutate(
        row_id = paste0(as.character(Type), "||", exposure_id)
      ) %>%
      mutate(
        UMAP1 = emb[match(row_id, rownames(Xn)), 1],
        UMAP2 = emb[match(row_id, rownames(Xn)), 2],
        sig_norm = sig_norm[match(row_id, rownames(Xn))]
      ) %>%
      filter(is.finite(UMAP1), is.finite(UMAP2)) %>%
      mutate(
        Type = factor(Type, levels = types),
        exposure_label = as.character(exposure_label)
      )
    
    # Save embedding table
    fwrite(as.data.table(emb_df), file.path(umap_dir, "UMAP_embeddings_FULL_signatures.tsv"), sep="\t")
    
    # Plot (trajectories + points)
    plot_umap(
      emb_df,
      file.path(umap_dir, "UMAP_FULL_signatures_Trajectories.png"),
      label_top_n = 20,
      title = "UMAP of proteomic signatures across exposures",
      subtitle = paste0(
        "Signature = selected_frac \u00D7 mean_beta_among_selected (Full model protein terms). ",
        "UMAP uses cosine distance on L2-normalized signatures; lines connect Types per exposure."
      )
    )
    
    # Optional: per-Type UMAP (sometimes easier to interpret)
    for (ty in types) {
      emb_ty <- emb_df %>% filter(Type == ty)
      if (nrow(emb_ty) < 10) next
      
      p_ty <- ggplot(emb_ty, aes(x = UMAP1, y = UMAP2)) +
        geom_point(aes(color = exposure_label), size=2, alpha=0.9) +
        labs(
          title = paste0("UMAP (subset): ", ty),
          subtitle = "Same global embedding, filtered to a single Type",
          color = "Exposure"
        ) +
        theme(legend.position = "right")
      
      ggsave(file.path(umap_dir, paste0("UMAP_FULL_signatures_", ty, ".png")),
             p_ty, width=14, height=9, dpi=300)
    }
  }
}

# ----------------------------
# Save summary tables used for plotting
# ----------------------------
fwrite(as.data.table(overall_primary_delta), file.path(out_dir, "overall_primary_delta.tsv"), sep="\t")
fwrite(as.data.table(fold_primary_full),    file.path(out_dir, "fold_primary_full.tsv"), sep="\t")
if (!is.null(sel_primary)) fwrite(as.data.table(sel_primary), file.path(out_dir, "selected_proteins_summary.tsv"), sep="\t")
if (nrow(cox_sum) > 0) fwrite(as.data.table(cox_sum), file.path(out_dir, "cox_all.tsv"), sep="\t")

message("\nDONE. Wrote plots + tables to: ", out_dir, "\n")
