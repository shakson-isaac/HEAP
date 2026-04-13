#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
  library(purrr)
  library(tidyr)
  library(ggplot2)
})

# ============================================================
# Config
# ============================================================
folder_name <- "/n/groups/patel/shakson_ukb/UK_Biobank/Data/Parallel/PES_Option1_fastv2/"
types <- c("Type1","Type2","Type3","Type4","Type5")

out_path <- "/n/groups/patel/shakson_ukb/UK_Biobank/Output/PES/"
out_dir <- file.path(out_path)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

theme_set(theme_bw(base_size = 12))
safe_fread <- function(path) tryCatch(fread(path), error = function(e) NULL)

# ============================================================
# Load all Types
# ============================================================
load_all_types <- function(types, folder_name) {
  overall_list <- list()
  fold_list <- list()
  cox_list <- list()
  
  for (ty in types) {
    ty_dir <- file.path(folder_name, ty)
    if (!dir.exists(ty_dir)) next
    files <- list.files(ty_dir, full.names = TRUE)
    
    overall_files <- files[str_detect(basename(files), "OverallMetrics\\.tsv$")]
    fold_files    <- files[str_detect(basename(files), "FoldMetrics\\.tsv$")]
    cox_files     <- files[str_detect(basename(files), "^Cox4_") & str_detect(basename(files), "\\.tsv$")]
    
    if (length(overall_files) > 0) {
      overall_dt <- rbindlist(lapply(overall_files, function(f) {
        dt <- safe_fread(f); if (is.null(dt)) return(NULL)
        dt[, `:=`(Type = ty, file = basename(f))]
        dt
      }), fill = TRUE)
      overall_list[[ty]] <- overall_dt
    }
    
    if (length(fold_files) > 0) {
      fold_dt <- rbindlist(lapply(fold_files, function(f) {
        dt <- safe_fread(f); if (is.null(dt)) return(NULL)
        dt[, `:=`(Type = ty, file = basename(f))]
        dt
      }), fill = TRUE)
      fold_list[[ty]] <- fold_dt
    }
    
    if (length(cox_files) > 0) {
      cox_dt <- rbindlist(lapply(cox_files, function(f) {
        dt <- safe_fread(f); if (is.null(dt)) return(NULL)
        dt[, `:=`(Type = ty, file = basename(f))]
        dt
      }), fill = TRUE)
      cox_list[[ty]] <- cox_dt
    }
  }
  
  list(
    overall = rbindlist(overall_list, fill = TRUE),
    fold    = rbindlist(fold_list, fill = TRUE),
    cox     = rbindlist(cox_list, fill = TRUE)
  )
}

dat <- load_all_types(types, folder_name)
overall_all <- dat$overall %>% as_tibble()
fold_all    <- dat$fold %>% as_tibble()
cox_all     <- dat$cox %>% as_tibble()

if (nrow(overall_all) == 0) stop("No OverallMetrics.tsv found under: ", folder_name)

# ============================================================
# Prep OverallMetrics (robust to missing columns)
# ============================================================
prep_overall <- function(overall_all, types) {
  metric_cols <- c("r2","rmse","auc","logloss","r2_code","mse_code",
                   "delta_r2","delta_auc","delta_logloss","delta_r2_code","delta_mse_code")
  
  ov_long <- overall_all %>%
    mutate(
      Type = factor(Type, levels = types),
      model = as.character(model),
      exposure_id = as.character(exposure_id),
      exposure_type = as.character(exposure_type)
    ) %>%
    pivot_longer(cols = any_of(metric_cols),
                 names_to = "metric_key", values_to = "value")
  
  key_map <- tibble::tribble(
    ~metric_key,        ~metric_name, ~kind,
    "r2",               "R2",         "metric",
    "auc",              "AUC",        "metric",
    "r2_code",          "R2_code",    "metric",
    "rmse",             "RMSE",       "aux",
    "logloss",          "LogLoss",    "aux",
    "mse_code",         "MSE_code",   "aux",
    "delta_r2",         "R2",         "delta",
    "delta_auc",        "AUC",        "delta",
    "delta_r2_code",    "R2_code",    "delta",
    "delta_logloss",    "LogLoss",    "delta",
    "delta_mse_code",   "MSE_code",   "delta"
  )
  
  ov_long <- ov_long %>%
    left_join(key_map, by = "metric_key") %>%
    filter(!is.na(metric_name)) %>%
    mutate(value = suppressWarnings(as.numeric(value)))
  
  primary_metrics <- c("R2","AUC","R2_code")
  
  cov_full <- ov_long %>%
    filter(kind == "metric",
           metric_name %in% primary_metrics,
           model %in% c("cov_only","prot_plus_cov")) %>%
    select(Type, exposure_id, exposure_type, metric_name, model, value) %>%
    group_by(Type, exposure_id, exposure_type, metric_name, model) %>%
    summarize(value = value[1], .groups = "drop") %>%
    pivot_wider(names_from = model, values_from = value) %>%
    rename(cov_metric = cov_only, full_metric = prot_plus_cov)
  
  delta_tbl <- ov_long %>%
    filter(kind == "delta", metric_name %in% primary_metrics) %>%
    select(Type, exposure_id, exposure_type, metric_name, value) %>%
    group_by(Type, exposure_id, exposure_type, metric_name) %>%
    summarize(delta_metric = value[1], .groups = "drop")
  
  out <- cov_full %>%
    left_join(delta_tbl, by = c("Type","exposure_id","exposure_type","metric_name")) %>%
    mutate(delta_metric = ifelse(is.na(delta_metric), full_metric - cov_metric, delta_metric)) %>%
    # Drop rows where we still can't compute anything
    filter(is.finite(delta_metric))
  
  out
}

overall_sum <- prep_overall(overall_all, types)

# ============================================================
# Prep FoldMetrics (robust)
# ============================================================
prep_fold <- function(fold_all, types) {
  if (nrow(fold_all) == 0) return(tibble())
  
  fold_all <- fold_all %>%
    mutate(
      Type = factor(Type, levels = types),
      exposure_id = as.character(exposure_id),
      exposure_type = as.character(exposure_type),
      family = as.character(family),
      fold = as.integer(fold)
    )
  
  candidates <- c("r2_prot","r2_cov","r2_full",
                  "auc_prot","auc_cov","auc_full",
                  "r2_code_prot","r2_code_cov","r2_code_full")
  
  fold_long <- fold_all %>%
    pivot_longer(cols = any_of(candidates),
                 names_to = "metric_key", values_to = "value") %>%
    mutate(
      metric_name = case_when(
        str_detect(metric_key, "^r2_") ~ "R2",
        str_detect(metric_key, "^auc_") ~ "AUC",
        str_detect(metric_key, "^r2_code_") ~ "R2_code",
        TRUE ~ NA_character_
      ),
      model = case_when(
        str_detect(metric_key, "_prot$") ~ "prot_only",
        str_detect(metric_key, "_cov$") ~ "cov_only",
        str_detect(metric_key, "_full$") ~ "prot_plus_cov",
        TRUE ~ NA_character_
      ),
      value = suppressWarnings(as.numeric(value))
    ) %>%
    filter(!is.na(metric_name), !is.na(model), is.finite(value))
  
  if (nrow(fold_long) == 0) return(tibble())
  
  fold_wide <- fold_long %>%
    select(Type, exposure_id, exposure_type, family, fold, metric_name, model, value,
           any_of(c("n_proteins_selected_full","n_proteins_selected_prot",
                    "n_selected_total_full","n_selected_total_prot"))) %>%
    distinct() %>%
    pivot_wider(names_from = model, values_from = value) %>%
    rename(metric_prot = prot_only, metric_cov = cov_only, metric_full = prot_plus_cov)
  
  # Ensure metric_full is present and finite
  fold_wide %>% filter(is.finite(metric_full))
}

fold_sum <- prep_fold(fold_all, types)

# ============================================================
# Prep Cox (robust)
# ============================================================
prep_cox <- function(cox_all, types) {
  if (nrow(cox_all) == 0) return(tibble())
  
  want <- c(
    "Type","exposure_id","exposure_type","disease_age_col","n","events",
    "HR_PES_M1_perSD","HR_PES_M1_L95","HR_PES_M1_U95","p_PES_M1",
    "HR_PES_M3_perSD","HR_PES_M3_L95","HR_PES_M3_U95","p_PES_M3",
    "p_LRT_M0_vs_M1","p_LRT_M2_vs_M3",
    "cindex_M0","cindex_M1","cindex_M2","cindex_M3",
    "delta_cindex_M0_to_M1","delta_cindex_M2_to_M3",
    "pes_used"
  )
  
  out <- cox_all %>%
    mutate(
      Type = factor(Type, levels = types),
      exposure_id = as.character(exposure_id),
      exposure_type = as.character(exposure_type),
      disease_age_col = as.character(disease_age_col)
    ) %>%
    select(any_of(want))
  
  # Add missing expected columns as NA
  missing_cols <- setdiff(want, names(out))
  if (length(missing_cols) > 0) {
    for (nm in missing_cols) out[[nm]] <- NA_real_
  }
  
  out
}

cox_sum <- prep_cox(cox_all, types)

# ============================================================
# Diagnostics (helps you see why plots might be empty)
# ============================================================
message("\n=== DIAGNOSTICS ===")
message("overall_sum rows: ", nrow(overall_sum))
if (nrow(overall_sum) > 0) print(overall_sum %>% count(metric_name, Type))
message("fold_sum rows: ", nrow(fold_sum))
if (nrow(fold_sum) > 0) print(fold_sum %>% count(metric_name, Type))
message("cox_sum rows: ", nrow(cox_sum))

# ============================================================
# PLOTS
# ============================================================

# ---------- Fig 1A (only if non-empty) ----------
if (nrow(overall_sum) > 0) {
  p1a <- ggplot(overall_sum, aes(x = Type, y = delta_metric, group = exposure_id)) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_line(alpha = 0.15) +
    geom_point(alpha = 0.35, size = 1) +
    facet_wrap(~metric_name, scales = "free_y") +
    labs(
      title = "Incremental predictive value of proteins beyond covariates",
      subtitle = "Δ metric = (prot_plus_cov) − (cov_only)",
      x = "Covariate specification", y = "Δ metric"
    )
  ggsave(file.path(out_dir, "Fig1A_delta_metric_vs_Type.png"), p1a, width = 11, height = 6, dpi = 300)
}

# ---------- Fig 2A / 2B ----------
# Only draw boxplots if there are at least 2 non-missing points per Type per facet
if (nrow(fold_sum) > 0) {
  
  fold_ok <- fold_sum %>%
    group_by(metric_name, Type) %>%
    mutate(n_group = sum(is.finite(metric_full))) %>%
    ungroup() %>%
    filter(n_group >= 2)
  
  if (nrow(fold_ok) > 0) {
    p2a <- ggplot(fold_ok, aes(x = Type, y = metric_full)) +
      geom_boxplot(outlier.alpha = 0.2) +
      facet_wrap(~metric_name, scales = "free_y") +
      labs(
        title = "Fold-to-fold stability of exposure prediction (prot_plus_cov)",
        x = "Covariate specification",
        y = "Fold metric (full model)"
      )
    ggsave(file.path(out_dir, "Fig2A_fold_metric_full_by_Type.png"), p2a, width = 11, height = 6, dpi = 300)
  } else {
    message("Skipping Fig2A: not enough non-missing fold metrics per group for boxplots.")
  }
  
  if ("n_proteins_selected_full" %in% names(fold_sum)) {
    sel_ok <- fold_sum %>%
      mutate(n_proteins_selected_full = suppressWarnings(as.numeric(n_proteins_selected_full))) %>%
      filter(is.finite(n_proteins_selected_full)) %>%
      group_by(Type) %>%
      mutate(n_group = n()) %>%
      ungroup() %>%
      filter(n_group >= 2)
    
    if (nrow(sel_ok) > 0) {
      p2b <- ggplot(sel_ok, aes(x = Type, y = n_proteins_selected_full)) +
        geom_boxplot(outlier.alpha = 0.2) +
        labs(
          title = "Number of selected proteins across folds (full model)",
          x = "Covariate specification",
          y = "# proteins selected"
        )
      ggsave(file.path(out_dir, "Fig2B_n_proteins_selected_full_by_Type.png"), p2b, width = 9, height = 5, dpi = 300)
    } else {
      message("Skipping Fig2B: not enough non-missing protein counts for boxplots.")
    }
  }
}

# ---------- Fig 3A forest ----------
if (nrow(cox_sum) > 0 &&
    all(c("HR_PES_M1_perSD","HR_PES_M1_L95","HR_PES_M1_U95",
          "HR_PES_M3_perSD","HR_PES_M3_L95","HR_PES_M3_U95") %in% names(cox_sum))) {
  
  cox_long <- cox_sum %>%
    mutate(row = paste0(exposure_id, " → ", disease_age_col)) %>%
    select(Type, row,
           HR_PES_M1_perSD, HR_PES_M1_L95, HR_PES_M1_U95,
           HR_PES_M3_perSD, HR_PES_M3_L95, HR_PES_M3_U95) %>%
    pivot_longer(
      cols = c(HR_PES_M1_perSD, HR_PES_M3_perSD),
      names_to = "model_key", values_to = "HR"
    ) %>%
    mutate(
      model = ifelse(model_key == "HR_PES_M1_perSD", "M1: cov + PES", "M3: cov + exposure + PES"),
      L95 = ifelse(model_key == "HR_PES_M1_perSD", HR_PES_M1_L95, HR_PES_M3_L95),
      U95 = ifelse(model_key == "HR_PES_M1_perSD", HR_PES_M1_U95, HR_PES_M3_U95),
      HR  = suppressWarnings(as.numeric(HR)),
      L95 = suppressWarnings(as.numeric(L95)),
      U95 = suppressWarnings(as.numeric(U95))
    ) %>%
    filter(is.finite(HR), is.finite(L95), is.finite(U95), HR > 0, L95 > 0, U95 > 0)
  
  if (nrow(cox_long) > 0) {
    cox_long <- cox_long %>%
      mutate(row = factor(row, levels = rev(unique(row))))
    
    p3a <- ggplot(cox_long, aes(x = HR, y = row)) +
      geom_vline(xintercept = 1, linetype = 2) +
      # ggplot2 >= 4.0: use geom_errorbar with orientation="y"
      geom_errorbar(aes(xmin = L95, xmax = U95), width = 0.2, alpha = 0.7, orientation = "y") +
      geom_point(size = 1.6) +
      scale_x_log10() +
      facet_grid(Type ~ model, scales = "free_y", space = "free_y") +
      labs(title = "PES association with incident disease (HR per SD)",
           x = "Hazard ratio (log scale)", y = "")
    ggsave(file.path(out_dir, "Fig3A_cox_forest_by_Type.png"), p3a, width = 14, height = 10, dpi = 300)
  } else {
    message("Skipping Fig3A: no finite HR/L95/U95 rows to plot.")
  }
}

# ---------- Fig 3B Δc-index ----------
if (nrow(cox_sum) > 0 && any(c("delta_cindex_M0_to_M1","delta_cindex_M2_to_M3") %in% names(cox_sum))) {
  
  cidx_long <- cox_sum %>%
    select(Type, any_of(c("delta_cindex_M0_to_M1","delta_cindex_M2_to_M3"))) %>%
    pivot_longer(cols = starts_with("delta_cindex"), names_to = "delta_kind", values_to = "delta_cindex") %>%
    mutate(delta_cindex = suppressWarnings(as.numeric(delta_cindex))) %>%
    filter(is.finite(delta_cindex))
  
  if (nrow(cidx_long) > 0) {
    cidx_long <- cidx_long %>%
      mutate(delta_kind = recode(delta_kind,
                                 delta_cindex_M0_to_M1 = "Δc-index: M0→M1",
                                 delta_cindex_M2_to_M3 = "Δc-index: M2→M3"))
    
    # guard: need at least 2 per group for boxplot
    cidx_ok <- cidx_long %>%
      group_by(delta_kind, Type) %>%
      mutate(n_group = n()) %>%
      ungroup() %>%
      filter(n_group >= 2)
    
    if (nrow(cidx_ok) > 0) {
      p3b <- ggplot(cidx_ok, aes(x = Type, y = delta_cindex)) +
        geom_hline(yintercept = 0, linetype = 2) +
        geom_point(alpha = 0.35, position = position_jitter(width = 0.12, height = 0)) +
        geom_boxplot(outlier.shape = NA, alpha = 0.4) +
        facet_wrap(~delta_kind, ncol = 1) +
        labs(title = "Incremental discrimination from PES", x = "Type", y = "Δc-index")
      ggsave(file.path(out_dir, "Fig3B_delta_cindex_by_Type.png"), p3b, width = 10, height = 8, dpi = 300)
    } else {
      message("Skipping Fig3B boxplots: not enough observations per (Type × delta_kind).")
    }
  }
}

# ============================================================
# Fig 3C alternative: Wald p-values for PES (what you have)
# ============================================================
if (nrow(cox_sum) > 0 && any(c("p_PES_M1","p_PES_M3") %in% names(cox_sum))) {
  
  wald_long <- cox_sum %>%
    select(Type, exposure_id, disease_age_col, any_of(c("p_PES_M1","p_PES_M3"))) %>%
    pivot_longer(cols = any_of(c("p_PES_M1","p_PES_M3")),
                 names_to = "test", values_to = "pval") %>%
    mutate(
      pval = suppressWarnings(as.numeric(pval)),
      pval = ifelse(is.na(pval) | pval <= 0, NA_real_, pval),
      mlog10p = -log10(pval),
      test = recode(test,
                    p_PES_M1 = "Wald p: PES in M1 (cov + PES)",
                    p_PES_M3 = "Wald p: PES in M3 (cov + exposure + PES)")
    ) %>%
    filter(is.finite(mlog10p))
  
  if (nrow(wald_long) > 0) {
    # guard: need >=2 per Type for boxplot
    wald_ok <- wald_long %>%
      group_by(test, Type) %>%
      mutate(n_group = n()) %>%
      ungroup() %>%
      filter(n_group >= 2)
    
    if (nrow(wald_ok) > 0) {
      p3c <- ggplot(wald_ok, aes(x = Type, y = mlog10p)) +
        geom_hline(yintercept = -log10(0.05), linetype = 2) +
        geom_point(alpha = 0.25, position = position_jitter(width = 0.12, height = 0)) +
        geom_boxplot(outlier.shape = NA, alpha = 0.4) +
        facet_wrap(~test, ncol = 1) +
        labs(title = "Evidence that PES adds information (Wald tests)",
             x = "Covariate specification",
             y = expression(-log[10](p)))
      ggsave(file.path(out_dir, "Fig3C_Wald_pvalues_by_Type.png"), p3c,
             width = 10, height = 8, dpi = 300)
    }
  }
}

# 
# # ---------- Optional Fig 3C LRT p-values ----------
# if (nrow(cox_sum) > 0 && any(c("p_LRT_M0_vs_M1","p_LRT_M2_vs_M3") %in% names(cox_sum))) {
#   lrt_long <- cox_sum %>%
#     select(Type, any_of(c("p_LRT_M0_vs_M1","p_LRT_M2_vs_M3"))) %>%
#     pivot_longer(cols = starts_with("p_LRT"), names_to = "lrt_kind", values_to = "pval") %>%
#     mutate(pval = suppressWarnings(as.numeric(pval)),
#            pval = ifelse(is.na(pval) | pval <= 0, NA_real_, pval),
#            mlog10p = -log10(pval)) %>%
#     filter(is.finite(mlog10p))
#   
#   if (nrow(lrt_long) > 0) {
#     lrt_long <- lrt_long %>%
#       mutate(lrt_kind = recode(lrt_kind,
#                                p_LRT_M0_vs_M1 = "LRT: M0 vs M1",
#                                p_LRT_M2_vs_M3 = "LRT: M2 vs M3"))
#     
#     lrt_ok <- lrt_long %>%
#       group_by(lrt_kind, Type) %>%
#       mutate(n_group = n()) %>%
#       ungroup() %>%
#       filter(n_group >= 2)
#     
#     if (nrow(lrt_ok) > 0) {
#       p3c <- ggplot(lrt_ok, aes(x = Type, y = mlog10p)) +
#         geom_hline(yintercept = -log10(0.05), linetype = 2) +
#         geom_point(alpha = 0.35, position = position_jitter(width = 0.12, height = 0)) +
#         geom_boxplot(outlier.shape = NA, alpha = 0.4) +
#         facet_wrap(~lrt_kind, ncol = 1) +
#         labs(title = "LRT evidence for added PES value", x = "Type", y = expression(-log[10](p)))
#       ggsave(file.path(out_dir, "Fig3C_LRT_pvalues_by_Type.png"), p3c, width = 10, height = 8, dpi = 300)
#     } else {
#       message("Skipping Fig3C boxplots: not enough observations per (Type × lrt_kind).")
#     }
#   }
# }

# ============================================================
# Save summaries
# ============================================================
fwrite(as.data.table(overall_sum), file.path(out_dir, "overall_summary_allTypes.tsv"), sep = "\t")
if (nrow(fold_sum) > 0) fwrite(as.data.table(fold_sum), file.path(out_dir, "fold_summary_allTypes.tsv"), sep = "\t")
if (nrow(cox_sum) > 0) fwrite(as.data.table(cox_sum), file.path(out_dir, "cox_summary_allTypes.tsv"), sep = "\t")

message("\nDONE. Output written to: ", out_dir)
