#!/usr/bin/env Rscript

# ============================================================
# Compare original GxE results (Type5) vs sensitivity analysis
# with ExCov and GxCov interactions (Type5_sens / Type5_sens_ext)
#
# INPUT:
#   - folders of batch .rds files produced by runProt_univar_assoc()
#   - each .rds contains list(train=..., test=...)
#
# OUTPUT:
#   1) merged comparison tables
#   2) summary metrics tables
#   3) visualizations:
#        - beta scatter (orig vs sens)
#        - -log10 p scatter
#        - volcano-style delta beta
#        - retention barplots
#        - QQ-like rank comparison
#   4) top retained / lost hits tables
#
# USAGE:
#   Rscript compare_type5_sensitivity.R Type5 Type5_sens
#
# OPTIONAL:
#   Rscript compare_type5_sensitivity.R Type5 Type5_sens_ext
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(ggplot2)
})

# ----------------------------
# Config
# ----------------------------
CFG <- list(
  out_root = "/n/groups/patel/shakson_ukb/UK_Biobank/Data/Parallel/Module2_sens/",
  out_subdir = "SensitivityCompare",
  eps = 1e-300
)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript compare_type5_sensitivity.R <orig_folder> <sens_folder>")
}

orig_folder <- "Type5" #args[1]   # e.g. "Type5"
sens_folder <- "Type5_sens_ext" #args[2]   # e.g. "Type5_sens"

orig_dir <- file.path(CFG$out_root, orig_folder)
sens_dir <- file.path(CFG$out_root, sens_folder)

out_dir <- file.path(CFG$out_root, CFG$out_subdir, paste0(orig_folder, "_vs_", sens_folder))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("Original folder: ", orig_dir)
message("Sensitivity folder: ", sens_dir)
message("Output folder: ", out_dir)

# ----------------------------
# Helpers
# ----------------------------
`%||%` <- function(x, y) if (!is.null(x)) x else y

safe_fread_rds <- function(path) {
  tryCatch(readRDS(path), error = function(e) NULL)
}

get_rds_files <- function(folder) {
  files <- list.files(folder, pattern = "\\.rds$", full.names = TRUE)
  files[order(files)]
}

safe_name <- function(x) {
  x <- gsub("[^A-Za-z0-9_\\-\\.]+", "_", x)
  x
}

pick_col <- function(df, candidates, required = FALSE) {
  hit <- intersect(candidates, names(df))[1]
  if (length(hit) == 0 || is.na(hit)) {
    if (required) stop("Missing required column. Tried: ", paste(candidates, collapse = ", "))
    return(NA_character_)
  }
  hit
}

coalesce_col <- function(df, candidates, out_name) {
  hit <- intersect(candidates, names(df))
  if (length(hit) == 0) {
    df[[out_name]] <- NA
  } else if (length(hit) == 1) {
    df[[out_name]] <- df[[hit]]
  } else {
    vals <- df[[hit[1]]]
    for (i in hit[-1]) vals <- dplyr::coalesce(vals, df[[i]])
    df[[out_name]] <- vals
  }
  df
}

standardize_statgxe <- function(df, source_label, split_label) {
  if (is.null(df) || nrow(df) == 0) return(NULL)
  df <- as.data.frame(df)
  
  # standardize key columns
  df <- coalesce_col(df, c("Estimate"), "beta")
  df <- coalesce_col(df, c("Std. Error", "Std_Error", "Std.Error"), "se")
  df <- coalesce_col(df, c("t value", "t_value", "t.value"), "tval")
  df <- coalesce_col(df, c("Pr(>|t|)", "p.value", "p_value", "P"), "pval")
  df <- coalesce_col(df, c("omicID"), "omicID_std")
  df <- coalesce_col(df, c("ID"), "termID")
  df <- coalesce_col(df, c("E_id"), "E_id_std")
  df <- coalesce_col(df, c("E_term"), "E_term_std")
  df <- coalesce_col(df, c("G_component"), "G_component_std")
  df <- coalesce_col(df, c("Category"), "Category_std")
  df <- coalesce_col(df, c("Eid"), "Eid_std")
  df <- coalesce_col(df, c("R2"), "R2_std")
  df <- coalesce_col(df, c("adj.R2", "adj_R2"), "adjR2_std")
  df <- coalesce_col(df, c("samplesize", "sample_size", "n"), "n_std")
  
  df$source <- source_label
  df$split <- split_label
  
  # make a stable exposure id
  df$exposure_id <- dplyr::coalesce(df$E_id_std, df$Eid_std, df$E_term_std)
  
  # stable omic id
  df$omicID <- df$omicID_std
  
  # stable G component
  df$G_component <- df$G_component_std
  
  # unique hit key
  df$hit_key <- paste(df$omicID, df$G_component, df$exposure_id, sep = "||")
  
  # numeric types
  num_cols <- c("beta", "se", "tval", "pval", "R2_std", "adjR2_std", "n_std")
  for (cc in intersect(num_cols, names(df))) df[[cc]] <- suppressWarnings(as.numeric(df[[cc]]))
  
  df
}

extract_statgxe_from_rds <- function(obj, source_label) {
  res <- list()
  
  for (split_label in c("train", "test")) {
    split_obj <- obj[[split_label]]
    if (is.null(split_obj) || length(split_obj) < 2) next
    
    # by construction:
    # [[1]] statE
    # [[2]] statGxE
    # [[3]] statR2
    # [[4]] statFblock
    statGxE <- split_obj[[2]]
    if (!is.null(statGxE) && nrow(statGxE) > 0) {
      res[[split_label]] <- standardize_statgxe(statGxE, source_label, split_label)
    }
  }
  
  bind_rows(res)
}

extract_fblock_from_rds <- function(obj, source_label) {
  out <- list()
  for (split_label in c("train", "test")) {
    split_obj <- obj[[split_label]]
    if (is.null(split_obj) || length(split_obj) < 4) next
    fblock <- split_obj[[4]]
    if (is.null(fblock) || nrow(fblock) == 0) next
    
    fblock <- as.data.frame(fblock)
    fblock$source <- source_label
    fblock$split <- split_label
    fblock <- coalesce_col(fblock, c("omicID"), "omicID_std")
    fblock <- coalesce_col(fblock, c("ID"), "exposure_id")
    fblock$omicID <- fblock$omicID_std
    fblock$key_FE <- paste(fblock$omicID, fblock$exposure_id, sep = "||")
    out[[split_label]] <- fblock
  }
  bind_rows(out)
}

load_folder_results <- function(folder, source_label) {
  files <- get_rds_files(folder)
  if (length(files) == 0) stop("No .rds files found in ", folder)
  
  message("Loading ", length(files), " files from ", folder)
  
  gxelist <- list()
  fblocklist <- list()
  
  for (i in seq_along(files)) {
    obj <- safe_fread_rds(files[i])
    if (is.null(obj)) next
    
    gxelist[[i]] <- extract_statgxe_from_rds(obj, source_label)
    fblocklist[[i]] <- extract_fblock_from_rds(obj, source_label)
  }
  
  gxe <- bind_rows(gxelist)
  fblock <- bind_rows(fblocklist)
  
  list(gxe = gxe, fblock = fblock)
}

bh_by_split <- function(df, pcol = "pval", outcol = "FDR") {
  df %>%
    group_by(split) %>%
    mutate(!!outcol := p.adjust(.data[[pcol]], method = "BH")) %>%
    ungroup()
}

summarize_counts <- function(comp_df, sig_col_orig, sig_col_sens, label) {
  tibble(
    summary_set = label,
    n_overlap = nrow(comp_df),
    n_sig_orig = sum(comp_df[[sig_col_orig]], na.rm = TRUE),
    n_sig_sens = sum(comp_df[[sig_col_sens]], na.rm = TRUE),
    retained_sig = sum(comp_df[[sig_col_orig]] & comp_df[[sig_col_sens]], na.rm = TRUE),
    lost_sig = sum(comp_df[[sig_col_orig]] & !comp_df[[sig_col_sens]], na.rm = TRUE),
    gained_sig = sum(!comp_df[[sig_col_orig]] & comp_df[[sig_col_sens]], na.rm = TRUE),
    concordant_direction = sum(sign(comp_df$beta_orig) == sign(comp_df$beta_sens), na.rm = TRUE),
    pct_retained_of_orig_sig = ifelse(sum(comp_df[[sig_col_orig]], na.rm = TRUE) == 0, NA_real_,
                                      100 * sum(comp_df[[sig_col_orig]] & comp_df[[sig_col_sens]], na.rm = TRUE) /
                                        sum(comp_df[[sig_col_orig]], na.rm = TRUE)),
    pct_concordant_direction = ifelse(nrow(comp_df) == 0, NA_real_,
                                      100 * mean(sign(comp_df$beta_orig) == sign(comp_df$beta_sens), na.rm = TRUE))
  )
}

theme_heap <- function() {
  theme_bw(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.2),
      strip.background = element_rect(fill = "grey95"),
      legend.position = "right"
    )
}

save_plot <- function(plot_obj, filename, width = 7, height = 6) {
  ggsave(
    filename = file.path(out_dir, filename),
    plot = plot_obj,
    width = width,
    height = height,
    dpi = 300
  )
}

# ----------------------------
# Load original + sensitivity
# ----------------------------
orig <- load_folder_results(orig_dir, "orig")
sens <- load_folder_results(sens_dir, "sens")

orig_gxe <- orig$gxe
sens_gxe <- sens$gxe

if (nrow(orig_gxe) == 0) stop("No GxE rows found in original results")
if (nrow(sens_gxe) == 0) stop("No GxE rows found in sensitivity results")

orig_gxe <- bh_by_split(orig_gxe, "pval", "FDR")
sens_gxe <- bh_by_split(sens_gxe, "pval", "FDR")

# ----------------------------
# Restrict to comparable rows
# ----------------------------
key_cols <- c("hit_key", "omicID", "G_component", "exposure_id", "split")

orig_keep <- orig_gxe %>%
  select(
    all_of(key_cols),
    beta_orig = beta,
    se_orig = se,
    tval_orig = tval,
    pval_orig = pval,
    FDR_orig = FDR,
    R2_orig = R2_std,
    adjR2_orig = adjR2_std,
    n_orig = n_std,
    termID_orig = termID,
    Category_orig = Category_std
  ) %>%
  distinct()

sens_keep <- sens_gxe %>%
  select(
    all_of(key_cols),
    beta_sens = beta,
    se_sens = se,
    tval_sens = tval,
    pval_sens = pval,
    FDR_sens = FDR,
    R2_sens = R2_std,
    adjR2_sens = adjR2_std,
    n_sens = n_std,
    termID_sens = termID,
    Category_sens = Category_std
  ) %>%
  distinct()

comp <- inner_join(orig_keep, sens_keep, by = key_cols)

if (nrow(comp) == 0) stop("No overlapping GxE hits between original and sensitivity results")

comp <- comp %>%
  mutate(
    delta_beta = beta_sens - beta_orig,
    abs_delta_beta = abs(delta_beta),
    beta_ratio = ifelse(abs(beta_orig) < 1e-12, NA_real_, beta_sens / beta_orig),
    same_direction = sign(beta_orig) == sign(beta_sens),
    neglog10p_orig = -log10(pmax(pval_orig, CFG$eps)),
    neglog10p_sens = -log10(pmax(pval_sens, CFG$eps)),
    sig_p_0.05_orig = pval_orig < 0.05,
    sig_p_0.05_sens = pval_sens < 0.05,
    sig_fdr_0.05_orig = FDR_orig < 0.05,
    sig_fdr_0.05_sens = FDR_sens < 0.05,
    retained_p = sig_p_0.05_orig & sig_p_0.05_sens,
    retained_fdr = sig_fdr_0.05_orig & sig_fdr_0.05_sens,
    rank_p_orig = rank(pval_orig, ties.method = "average"),
    rank_p_sens = rank(pval_sens, ties.method = "average")
  )

# ----------------------------
# Merge F-block exposure-level summaries
# ----------------------------
orig_f <- orig$fblock
sens_f <- sens$fblock

if (nrow(orig_f) > 0 && nrow(sens_f) > 0) {
  orig_f2 <- orig_f %>%
    select(
      key_FE,
      omicID,
      exposure_id,
      split,
      p_GxE_joint_orig = p_GxE_joint,
      p_GcisxE_orig = p_GcisxE,
      p_GtrxE_orig = p_GtrxE
    ) %>%
    distinct()
  
  sens_f2 <- sens_f %>%
    select(
      key_FE,
      omicID,
      exposure_id,
      split,
      p_GxE_joint_sens = p_GxE_joint,
      p_GcisxE_sens = p_GcisxE,
      p_GtrxE_sens = p_GtrxE
    ) %>%
    distinct()
  
  comp_FE <- inner_join(orig_f2, sens_f2, by = c("key_FE", "omicID", "exposure_id", "split"))
  fwrite(comp_FE, file.path(out_dir, "comparison_exposure_level_Fblock.tsv"), sep = "\t")
}

# ----------------------------
# Save main comparison table
# ----------------------------
fwrite(comp, file.path(out_dir, "comparison_GxE_hit_level.tsv"), sep = "\t")

# ----------------------------
# Summary tables
# ----------------------------
summary_split_p <- comp %>%
  group_by(split) %>%
  group_modify(~ summarize_counts(.x, "sig_p_0.05_orig", "sig_p_0.05_sens", "P<0.05")) %>%
  ungroup()

summary_split_fdr <- comp %>%
  group_by(split) %>%
  group_modify(~ summarize_counts(.x, "sig_fdr_0.05_orig", "sig_fdr_0.05_sens", "FDR<0.05")) %>%
  ungroup()

summary_by_component <- comp %>%
  group_by(split, G_component) %>%
  summarise(
    n_overlap = n(),
    cor_beta = suppressWarnings(cor(beta_orig, beta_sens, use = "pairwise.complete.obs")),
    cor_neglog10p = suppressWarnings(cor(neglog10p_orig, neglog10p_sens, use = "pairwise.complete.obs")),
    median_abs_delta_beta = median(abs_delta_beta, na.rm = TRUE),
    pct_same_direction = 100 * mean(same_direction, na.rm = TRUE),
    n_sig_orig_p = sum(sig_p_0.05_orig, na.rm = TRUE),
    n_sig_sens_p = sum(sig_p_0.05_sens, na.rm = TRUE),
    retained_sig_p = sum(retained_p, na.rm = TRUE),
    n_sig_orig_fdr = sum(sig_fdr_0.05_orig, na.rm = TRUE),
    n_sig_sens_fdr = sum(sig_fdr_0.05_sens, na.rm = TRUE),
    retained_sig_fdr = sum(retained_fdr, na.rm = TRUE),
    .groups = "drop"
  )

summary_global <- tibble(
  n_overlap = nrow(comp),
  cor_beta = suppressWarnings(cor(comp$beta_orig, comp$beta_sens, use = "pairwise.complete.obs")),
  cor_neglog10p = suppressWarnings(cor(comp$neglog10p_orig, comp$neglog10p_sens, use = "pairwise.complete.obs")),
  median_abs_delta_beta = median(comp$abs_delta_beta, na.rm = TRUE),
  pct_same_direction = 100 * mean(comp$same_direction, na.rm = TRUE)
)

fwrite(summary_split_p, file.path(out_dir, "summary_counts_p005.tsv"), sep = "\t")
fwrite(summary_split_fdr, file.path(out_dir, "summary_counts_fdr005.tsv"), sep = "\t")
fwrite(summary_by_component, file.path(out_dir, "summary_by_component.tsv"), sep = "\t")
fwrite(summary_global, file.path(out_dir, "summary_global.tsv"), sep = "\t")

# ----------------------------
# Top retained / lost / changed hits
# ----------------------------
top_retained_fdr <- comp %>%
  filter(sig_fdr_0.05_orig, sig_fdr_0.05_sens) %>%
  arrange(FDR_sens, FDR_orig, desc(abs(beta_sens))) %>%
  head(200)

top_lost_fdr <- comp %>%
  filter(sig_fdr_0.05_orig, !sig_fdr_0.05_sens) %>%
  arrange(FDR_orig, desc(abs(delta_beta))) %>%
  head(200)

top_changed_beta <- comp %>%
  arrange(desc(abs_delta_beta)) %>%
  head(200)

top_flip_direction <- comp %>%
  filter(!same_direction) %>%
  arrange(FDR_orig, FDR_sens, desc(abs_delta_beta)) %>%
  head(200)

fwrite(top_retained_fdr, file.path(out_dir, "top_retained_FDR_hits.tsv"), sep = "\t")
fwrite(top_lost_fdr, file.path(out_dir, "top_lost_FDR_hits.tsv"), sep = "\t")
fwrite(top_changed_beta, file.path(out_dir, "top_changed_beta_hits.tsv"), sep = "\t")
fwrite(top_flip_direction, file.path(out_dir, "top_direction_flip_hits.tsv"), sep = "\t")

# ----------------------------
# Plots: Beta scatter
# ----------------------------
plot_beta_scatter <- comp %>%
  ggplot(aes(x = beta_orig, y = beta_sens)) +
  geom_point(alpha = 0.45, size = 1) +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  facet_grid(split ~ G_component, scales = "free") +
  labs(
    title = paste0(orig_folder, " vs ", sens_folder, ": GxE beta comparison"),
    x = "Original beta",
    y = "Sensitivity beta"
  ) +
  theme_heap()

save_plot(plot_beta_scatter, "plot_beta_scatter.png", width = 9, height = 7)

# ----------------------------
# Plots: -log10 p scatter
# ----------------------------
plot_p_scatter <- comp %>%
  ggplot(aes(x = neglog10p_orig, y = neglog10p_sens)) +
  geom_point(alpha = 0.45, size = 1) +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  facet_grid(split ~ G_component, scales = "free") +
  labs(
    title = paste0(orig_folder, " vs ", sens_folder, ": GxE significance comparison"),
    x = expression(-log[10](p["original"])),
    y = expression(-log[10](p["sensitivity"]))
  ) +
  theme_heap()

save_plot(plot_p_scatter, "plot_neglog10p_scatter.png", width = 9, height = 7)

# ----------------------------
# Plots: delta beta
# ----------------------------
plot_delta_beta <- comp %>%
  mutate(status = case_when(
    sig_fdr_0.05_orig & sig_fdr_0.05_sens ~ "Retained FDR",
    sig_fdr_0.05_orig & !sig_fdr_0.05_sens ~ "Lost FDR",
    !sig_fdr_0.05_orig & sig_fdr_0.05_sens ~ "Gained FDR",
    TRUE ~ "Not FDR sig"
  )) %>%
  ggplot(aes(x = beta_orig, y = delta_beta)) +
  geom_point(aes(shape = status), alpha = 0.6, size = 1.2) +
  geom_hline(yintercept = 0, linetype = 2) +
  facet_grid(split ~ G_component, scales = "free") +
  labs(
    title = "Change in GxE beta after covariate-interaction adjustment",
    x = "Original beta",
    y = "Sensitivity beta - Original beta"
  ) +
  theme_heap()

save_plot(plot_delta_beta, "plot_delta_beta.png", width = 10, height = 7)

# ----------------------------
# Plots: retention barplot
# ----------------------------
retention_bar_df <- bind_rows(
  comp %>%
    group_by(split, G_component) %>%
    summarise(
      threshold = "P<0.05",
      category = "Original sig",
      n = sum(sig_p_0.05_orig, na.rm = TRUE),
      .groups = "drop"
    ),
  comp %>%
    group_by(split, G_component) %>%
    summarise(
      threshold = "P<0.05",
      category = "Retained sig",
      n = sum(retained_p, na.rm = TRUE),
      .groups = "drop"
    ),
  comp %>%
    group_by(split, G_component) %>%
    summarise(
      threshold = "FDR<0.05",
      category = "Original sig",
      n = sum(sig_fdr_0.05_orig, na.rm = TRUE),
      .groups = "drop"
    ),
  comp %>%
    group_by(split, G_component) %>%
    summarise(
      threshold = "FDR<0.05",
      category = "Retained sig",
      n = sum(retained_fdr, na.rm = TRUE),
      .groups = "drop"
    )
)

plot_retention_bar <- retention_bar_df %>%
  ggplot(aes(x = G_component, y = n, fill = category)) +
  geom_col(position = "dodge") +
  facet_grid(split ~ threshold) +
  labs(
    title = "Retention of significant GxE hits after covariate-interaction adjustment",
    x = "Genetic component",
    y = "Number of hits"
  ) +
  theme_heap()

save_plot(plot_retention_bar, "plot_retention_bar.png", width = 9, height = 7)

# ----------------------------
# Plots: retention percentage
# ----------------------------
retention_pct_df <- comp %>%
  group_by(split, G_component) %>%
  summarise(
    orig_sig_p = sum(sig_p_0.05_orig, na.rm = TRUE),
    retained_p = sum(retained_p, na.rm = TRUE),
    orig_sig_fdr = sum(sig_fdr_0.05_orig, na.rm = TRUE),
    retained_fdr = sum(retained_fdr, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    pct_retained_p = ifelse(orig_sig_p == 0, NA_real_, 100 * retained_p / orig_sig_p),
    pct_retained_fdr = ifelse(orig_sig_fdr == 0, NA_real_, 100 * retained_fdr / orig_sig_fdr)
  ) %>%
  pivot_longer(
    cols = c(pct_retained_p, pct_retained_fdr),
    names_to = "metric",
    values_to = "pct"
  ) %>%
  mutate(
    metric = recode(metric,
                    pct_retained_p = "P<0.05",
                    pct_retained_fdr = "FDR<0.05")
  )

plot_retention_pct <- retention_pct_df %>%
  ggplot(aes(x = G_component, y = pct, fill = metric)) +
  geom_col(position = "dodge") +
  facet_wrap(~ split) +
  labs(
    title = "Percent of original significant hits retained",
    x = "Genetic component",
    y = "Retained (%)"
  ) +
  ylim(0, 100) +
  theme_heap()

save_plot(plot_retention_pct, "plot_retention_percent.png", width = 8, height = 6)

# ----------------------------
# Plots: rank comparison
# ----------------------------
plot_rank <- comp %>%
  ggplot(aes(x = rank_p_orig, y = rank_p_sens)) +
  geom_point(alpha = 0.35, size = 1) +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  facet_grid(split ~ G_component, scales = "free") +
  labs(
    title = "Rank stability of GxE p-values",
    x = "Original p-value rank",
    y = "Sensitivity p-value rank"
  ) +
  theme_heap()

save_plot(plot_rank, "plot_rank_comparison.png", width = 9, height = 7)

# ----------------------------
# Optional: category-level summaries
# ----------------------------
if ("Category_orig" %in% names(comp)) {
  category_summary <- comp %>%
    mutate(Category = Category_orig %||% Category_sens) %>%
    group_by(split, G_component, Category) %>%
    summarise(
      n_hits = n(),
      cor_beta = suppressWarnings(cor(beta_orig, beta_sens, use = "pairwise.complete.obs")),
      pct_same_direction = 100 * mean(same_direction, na.rm = TRUE),
      retained_fdr = sum(retained_fdr, na.rm = TRUE),
      orig_fdr = sum(sig_fdr_0.05_orig, na.rm = TRUE),
      pct_retained_fdr = ifelse(orig_fdr == 0, NA_real_, 100 * retained_fdr / orig_fdr),
      .groups = "drop"
    ) %>%
    arrange(split, G_component, desc(orig_fdr), desc(n_hits))
  
  fwrite(category_summary, file.path(out_dir, "summary_by_category.tsv"), sep = "\t")
  
  category_plot_df <- category_summary %>%
    filter(!is.na(pct_retained_fdr), orig_fdr >= 5)
  
  if (nrow(category_plot_df) > 0) {
    plot_category_retention <- category_plot_df %>%
      ggplot(aes(x = reorder(Category, pct_retained_fdr), y = pct_retained_fdr)) +
      geom_col() +
      coord_flip() +
      facet_grid(split ~ G_component, scales = "free_y") +
      labs(
        title = "FDR-hit retention by exposure category",
        x = "Exposure category",
        y = "Retained FDR hits (%)"
      ) +
      theme_heap()
    
    save_plot(plot_category_retention, "plot_category_retention.png", width = 10, height = 9)
  }
}

# ----------------------------
# Write an interpretation-ready text summary
# ----------------------------
summary_text <- c(
  paste0("Original folder: ", orig_folder),
  paste0("Sensitivity folder: ", sens_folder),
  paste0("Number of overlapping GxE rows: ", nrow(comp)),
  paste0("Global beta correlation: ", round(summary_global$cor_beta, 4)),
  paste0("Global -log10(p) correlation: ", round(summary_global$cor_neglog10p, 4)),
  paste0("Median |delta beta|: ", round(summary_global$median_abs_delta_beta, 4)),
  paste0("Percent same direction: ", round(summary_global$pct_same_direction, 2), "%"),
  "",
  "By split / threshold:",
  capture.output(print(summary_split_p)),
  "",
  capture.output(print(summary_split_fdr)),
  "",
  "By split / G component:",
  capture.output(print(summary_by_component))
)

writeLines(summary_text, con = file.path(out_dir, "summary_readme.txt"))

message("Done. Outputs written to: ", out_dir)


test_comp <- comp %>%
              filter(split == "test")


