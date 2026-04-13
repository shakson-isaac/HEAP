#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(stringr)
})

# ============================================================
# Config
# ============================================================
IN_DIR  <- "/n/groups/patel/shakson_ukb/UK_Biobank/Output/MR_edges/summary"
OUTDIR  <- "/n/groups/patel/shakson_ukb/UK_Biobank/Output/MRplots"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

ADJ_METHOD <- "BH"
ALPHA_MR   <- 0.05    # q threshold for MR hit
ALPHA_Q    <- 0.05    # heterogeneity flag threshold (IVW Q_pval)
ALPHA_EG   <- 0.05    # Egger intercept pval flag threshold

# MR methods to consider as the "main" estimate
KEEP_METHODS <- c("Inverse variance weighted", "Wald ratio")

# How to treat missing sensitivity (common when nsnp small):
# "pass" => don't penalize if sensitivity not available (recommended)
# "fail" => require sensitivity available and non-flagged
MISSING_SENS_POLICY <- "pass"  # or "fail"


# ============================================================
# Helpers
# ============================================================

read_rds_safe <- function(path) {
  if (!file.exists(path)) stop("Missing: ", path)
  readRDS(path)
}

# Standardize the MR table inside each RDS:
# - keep IVW/Wald
# - compute q within edge_dir
prep_mr <- function(mr_dt) {
  if (is.null(mr_dt) || nrow(mr_dt) == 0) return(data.table())
  
  dt <- as.data.table(mr_dt)
  # Keep main methods only
  if ("method" %in% names(dt)) dt <- dt[method %in% KEEP_METHODS]
  
  if (nrow(dt) == 0) return(data.table())
  
  # Compute q within each edge_dir (cis/trans separated naturally by edge_dir)
  dt[, q := p.adjust(pval, method = ADJ_METHOD), by = .(edge_dir)]
  
  # Canonical keys
  # (src_id, tgt_id already exist from your loader; keep them)
  keep_cols <- intersect(
    c("edge_type","edge_dir","src_id","tgt_id",
      "id.exposure","id.outcome","exposure","outcome",
      "method","nsnp","b","se","pval","q"),
    names(dt)
  )
  dt <- dt[, ..keep_cols]
  
  dt
}

# Standardize heterogeneity:
# Use the IVW row’s Q_pval as a diagnostic for the instrument set
prep_het <- function(het_dt) {
  if (is.null(het_dt) || nrow(het_dt) == 0) return(data.table())
  dt <- as.data.table(het_dt)
  
  # Prefer IVW heterogeneity (most standard)
  if ("method" %in% names(dt)) {
    dt_ivw <- dt[method == "Inverse variance weighted"]
    if (nrow(dt_ivw) > 0) dt <- dt_ivw
  }
  
  # Keep key columns
  # Your example includes Q, Q_df, Q_pval
  keep_cols <- intersect(
    c("edge_type","edge_dir","src_id","tgt_id","Q","Q_df","Q_pval","method"),
    names(dt)
  )
  dt <- dt[, ..keep_cols]
  
  # De-duplicate (sometimes multiple rows remain)
  # keep first per edge
  if (all(c("edge_dir","src_id","tgt_id") %in% names(dt))) {
    setkeyv(dt, c("edge_dir","src_id","tgt_id"))
    dt <- unique(dt)
  }
  
  dt
}

# Standardize pleiotropy (Egger intercept):
prep_ple <- function(ple_dt) {
  if (is.null(ple_dt) || nrow(ple_dt) == 0) return(data.table())
  dt <- as.data.table(ple_dt)
  
  # Your example: egger_intercept, se, pval
  keep_cols <- intersect(
    c("edge_type","edge_dir","src_id","tgt_id","egger_intercept","se","pval"),
    names(dt)
  )
  dt <- dt[, ..keep_cols]
  setnames(dt, old = intersect("pval", names(dt)), new = "egger_pval", skip_absent = TRUE)
  
  if (all(c("edge_dir","src_id","tgt_id") %in% names(dt))) {
    setkeyv(dt, c("edge_dir","src_id","tgt_id"))
    dt <- unique(dt)
  }
  
  dt
}

# Convert edge_dir to nice label + ordering for plots
edge_label_from_dir <- function(edge_dir) {
  dplyr::case_when(
    edge_dir == "E_to_P"      ~ "E→P",
    edge_dir == "E_to_D"      ~ "E→D",
    edge_dir == "D_to_E"      ~ "D→E",
    edge_dir == "D_to_P"      ~ "D→P",
    edge_dir == "Pcis_to_D"   ~ "P→D (cis)",
    edge_dir == "Ptrans_to_D" ~ "P→D (trans)",
    edge_dir == "Pcis_to_E"   ~ "P→E (cis)",
    edge_dir == "Ptrans_to_E" ~ "P→E (trans)",
    TRUE ~ edge_dir
  )
}
pair_from_dir <- function(edge_dir) {
  dplyr::case_when(
    edge_dir %in% c("E_to_P","Pcis_to_E","Ptrans_to_E") ~ "E-P",
    edge_dir %in% c("E_to_D","D_to_E")                  ~ "E-D",
    edge_dir %in% c("Pcis_to_D","Ptrans_to_D","D_to_P") ~ "P-D",
    TRUE ~ "Other"
  )
}
direction_from_dir <- function(edge_dir) {
  dplyr::case_when(
    edge_dir %in% c("E_to_P","E_to_D","Pcis_to_D","Ptrans_to_D") ~ "forward",
    edge_dir %in% c("Pcis_to_E","Ptrans_to_E","D_to_E","D_to_P") ~ "reverse",
    TRUE ~ NA_character_
  )
}

# ============================================================
# Load all RDS
# ============================================================
RDS_FILES <- c(
  PD = file.path(IN_DIR, "PDres.rds"),
  EP = file.path(IN_DIR, "EPres.rds"),
  ED = file.path(IN_DIR, "EDres.rds"),
  DE = file.path(IN_DIR, "DEres.rds"),
  PE = file.path(IN_DIR, "PEres.rds"),
  DP = file.path(IN_DIR, "DPres.rds")
)

obj_list <- lapply(RDS_FILES, read_rds_safe)

# ============================================================
# Build long table of MR + sensitivity
# ============================================================
mr_all  <- rbindlist(lapply(obj_list, function(x) prep_mr(x$mr)),  fill = TRUE)
het_all <- rbindlist(lapply(obj_list, function(x) prep_het(x$het)), fill = TRUE)
ple_all <- rbindlist(lapply(obj_list, function(x) prep_ple(x$ple)), fill = TRUE)

if (nrow(mr_all) == 0) stop("No MR rows found after filtering methods.")

# Join sensitivity onto MR rows by (edge_dir, src_id, tgt_id)
# (These keys match what you created in the loader.)
setDT(mr_all); setDT(het_all); setDT(ple_all)
setkeyv(mr_all,  c("edge_dir","src_id","tgt_id"))
if (nrow(het_all) > 0) setkeyv(het_all, c("edge_dir","src_id","tgt_id"))
if (nrow(ple_all) > 0) setkeyv(ple_all, c("edge_dir","src_id","tgt_id"))

mr_sens <- copy(mr_all)
if (nrow(het_all) > 0) mr_sens <- het_all[mr_sens]  # left join
if (nrow(ple_all) > 0) mr_sens <- ple_all[mr_sens]  # left join

# Add plotting labels + flags
edge_levels_optA <- c(
  "E→P",
  "P→E (cis)", "P→E (trans)",
  "E→D",
  "D→E",
  "P→D (cis)", "P→D (trans)",
  "D→P"
)
edge_levels_display <- rev(edge_levels_optA)  # for y-axis order

mr_sens <- mr_sens %>%
  mutate(
    edge_label = edge_label_from_dir(edge_dir),
    pair = factor(pair_from_dir(edge_dir), levels = c("E-P","E-D","P-D")),
    direction = direction_from_dir(edge_dir),
    
    mr_hit = !is.na(q) & q < ALPHA_MR,
    
    # flags (only meaningful when present)
    het_flag = !is.na(Q_pval) & (Q_pval < ALPHA_Q),
    pleio_flag = !is.na(egger_pval) & (egger_pval < ALPHA_EG),
    
    # pass rules depending on policy
    het_pass = dplyr::case_when(
      MISSING_SENS_POLICY == "pass" & is.na(Q_pval) ~ TRUE,
      MISSING_SENS_POLICY == "fail" & is.na(Q_pval) ~ FALSE,
      TRUE ~ !het_flag
    ),
    pleio_pass = dplyr::case_when(
      MISSING_SENS_POLICY == "pass" & is.na(egger_pval) ~ TRUE,
      MISSING_SENS_POLICY == "fail" & is.na(egger_pval) ~ FALSE,
      TRUE ~ !pleio_flag
    ),
    
    sens_pass = het_pass & pleio_pass,
    mr_hit_sens = mr_hit & sens_pass,
    
    edge_label = factor(edge_label, levels = edge_levels_display)
  )

# ============================================================
# Summaries: hit counts before vs after sensitivity
# ============================================================
summ_counts <- mr_sens %>%
  group_by(pair, edge_label) %>%
  summarise(
    n_test = sum(!is.na(q)),                 # tested = has MR p/q
    n_hit  = sum(mr_hit, na.rm = TRUE),      # q<0.05
    n_hit_sens = sum(mr_hit_sens, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    n_removed = n_hit - n_hit_sens,
    pct_hit = 100 * n_hit / pmax(n_test, 1),
    pct_hit_sens = 100 * n_hit_sens / pmax(n_test, 1),
    
    # label strings
    lab_hit      = paste0(n_hit, "/", n_test),
    lab_hit_sens = paste0(n_hit_sens, "/", n_test),
    lab_removed  = paste0("−", n_removed, "/", n_hit)
  )


# Also: among hits, what fraction are flagged?
summ_flags <- mr_sens %>%
  filter(mr_hit) %>%
  group_by(pair, edge_label) %>%
  summarise(
    n_hit = n(),
    frac_het_flag = mean(het_flag, na.rm = TRUE),
    frac_pleio_flag = mean(pleio_flag, na.rm = TRUE),
    frac_any_flag = mean(het_flag | pleio_flag, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    pct_any_flag = 100 * frac_any_flag,
    pct_het_flag = 100 * frac_het_flag,
    pct_pleio_flag = 100 * frac_pleio_flag
  )

# Save tables for reference
fwrite(as.data.table(summ_counts), file.path(OUTDIR, "MR_hit_counts_before_after_sensitivity.tsv"), sep = "\t")
fwrite(as.data.table(summ_flags),  file.path(OUTDIR, "MR_hit_flag_rates_among_hits.tsv"), sep = "\t")

# ============================================================
# Plot 1: Hit-rate before vs after sensitivity (two bars per edge)
# ============================================================
plot_df <- summ_counts %>%
  select(pair, edge_label, n_test, n_hit, n_hit_sens, pct_hit, pct_hit_sens) %>%
  tidyr::pivot_longer(
    cols = c(pct_hit, pct_hit_sens),
    names_to = "stage",
    values_to = "pct"
  ) %>%
  mutate(
    stage = recode(stage,
                   pct_hit = "MR hits (q<0.05)",
                   pct_hit_sens = "After sensitivity filter"),
    stage = factor(stage, levels = c("MR hits (q<0.05)", "After sensitivity filter"))
  )

# Headroom for labels
ymax <- max(plot_df$pct, na.rm = TRUE)
ymax_pad <- ymax * 1.35

# Label: show removed counts on the "after" bar
removed_lab <- summ_counts %>%
  mutate(label = ifelse(n_removed > 0, paste0("−", n_removed), "0")) %>%
  select(pair, edge_label, label)

plot_df2 <- plot_df %>%
  left_join(removed_lab, by = c("pair","edge_label")) %>%
  mutate(
    lab2 = ifelse(stage == "After sensitivity filter", label, NA_character_)
  )

p1 <- ggplot(plot_df2, aes(x = edge_label, y = pct, fill = stage)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7) +
  geom_text(
    aes(label = lab2),
    position = position_dodge(width = 0.75),
    hjust = -0.15, size = 3, na.rm = TRUE
  ) +
  coord_flip(clip = "off") +
  facet_grid(pair ~ ., scales = "free_y", space = "free_y", drop = TRUE) +
  scale_y_continuous(limits = c(0, ymax_pad)) +
  theme_bw() +
  labs(
    x = NULL,
    y = "% significant",
    fill = NULL,
    title = "MR hit-rate before vs after sensitivity filtering",
    subtitle = paste0("Sensitivity filter: IVW Q_pval≥", ALPHA_Q,
                      " and Egger intercept p≥", ALPHA_EG,
                      " (missing=", MISSING_SENS_POLICY, ")")
  ) +
  theme(
    strip.background = element_rect(fill = NA),
    plot.margin = margin(t = 10, r = 55, b = 10, l = 10)
  )

ggsave(file.path(OUTDIR, "MR_hit_rate_before_after_sensitivity.png"), p1, width = 9, height = 6, dpi = 600)
ggsave(file.path(OUTDIR, "MR_hit_rate_before_after_sensitivity.svg"), p1, width = 9, height = 6)
print(p1)

# ============================================================
# Plot 2: Among MR hits, % flagged by heterogeneity/pleiotropy
# (nice supplemental panel)
# ============================================================
flag_long <- summ_flags %>%
  select(pair, edge_label, pct_het_flag, pct_pleio_flag) %>%
  tidyr::pivot_longer(
    cols = c(pct_het_flag, pct_pleio_flag),
    names_to = "flag",
    values_to = "pct"
  ) %>%
  mutate(
    flag = recode(flag,
                  pct_het_flag = "Heterogeneity (Q_pval<0.05)",
                  pct_pleio_flag = "Egger intercept (p<0.05)")
  )

ymax2 <- max(flag_long$pct, na.rm = TRUE)
p2 <- ggplot(flag_long, aes(x = edge_label, y = pct, fill = flag)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7) +
  coord_flip() +
  facet_grid(pair ~ ., scales = "free_y", space = "free_y", drop = TRUE) +
  scale_y_continuous(limits = c(0, max(10, ymax2 * 1.15))) +
  theme_bw() +
  labs(
    x = NULL,
    y = "% of MR hits flagged",
    fill = NULL,
    title = "Sensitivity flags among MR hits"
  ) +
  theme(strip.background = element_rect(fill = NA))

ggsave(file.path(OUTDIR, "MR_sensitivity_flag_rates_among_hits.png"), p2, width = 9, height = 6, dpi = 600)
ggsave(file.path(OUTDIR, "MR_sensitivity_flag_rates_among_hits.svg"), p2, width = 9, height = 6)
print(p2)

# ============================================================
# Plot 3: Distributions of Q_pval and Egger p-values (supp)
# ============================================================
# Note: only where available
p3a <- mr_sens %>%
  filter(!is.na(Q_pval)) %>%
  ggplot(aes(x = Q_pval)) +
  geom_histogram(bins = 50) +
  facet_grid(pair ~ edge_label, scales = "free_y", space = "free_y", drop = TRUE) +
  theme_bw() +
  labs(x = "IVW Cochran Q p-value", y = "Count", title = "Heterogeneity p-values (where available)") +
  theme(strip.text = element_text(size = 7))

ggsave(file.path(OUTDIR, "MR_heterogeneity_Qpval_hist.png"), p3a, width = 12, height = 6, dpi = 600)
print(p3a)

p3b <- mr_sens %>%
  filter(!is.na(egger_pval)) %>%
  ggplot(aes(x = egger_pval)) +
  geom_histogram(bins = 50) +
  facet_grid(pair ~ edge_label, scales = "free_y", space = "free_y", drop = TRUE) +
  theme_bw() +
  labs(x = "Egger intercept p-value", y = "Count", title = "Directional pleiotropy p-values (where available)") +
  theme(strip.text = element_text(size = 7))

ggsave(file.path(OUTDIR, "MR_pleiotropy_egger_pval_hist.png"), p3b, width = 12, height = 6, dpi = 600)
print(p3b)

cat("\nDone. Outputs written to: ", OUTDIR, "\n")
cat("Key tables:\n",
    " - MR_hit_counts_before_after_sensitivity.tsv\n",
    " - MR_hit_flag_rates_among_hits.tsv\n")
