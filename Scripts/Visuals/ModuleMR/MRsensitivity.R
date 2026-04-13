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
CFG <- list(
  rds_dir   = "/n/groups/patel/shakson_ukb/UK_Biobank/Output/MR_edges/summary",
  out_dir   = "/n/groups/patel/shakson_ukb/UK_Biobank/Output/MRplots",
  adj_method = "BH",
  alpha_q    = 0.05,  # MR hit threshold on adjusted p-value
  alpha_sens = 0.05   # sensitivity threshold for Q_pval and Egger pval
)

dir.create(CFG$out_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# Helpers
# ============================================================

read_res <- function(prefix, rds_dir) {
  path <- file.path(rds_dir, paste0(prefix, "res.rds"))
  if (!file.exists(path)) stop("Missing: ", path)
  x <- readRDS(path)
  # x is list(mr=..., het=..., ple=...)
  if (is.null(x$mr)) x$mr <- data.table()
  if (is.null(x$het)) x$het <- data.table()
  if (is.null(x$ple)) x$ple <- data.table()
  x$mr  <- as.data.table(x$mr)
  x$het <- as.data.table(x$het)
  x$ple <- as.data.table(x$ple)
  x
}

# ---- FIXED EDGE LABEL MAP (vectorized) ----
edge_label_map <- c(
  "E_to_P"      = "E→P",
  "Pcis_to_E"   = "P→E (cis)",
  "Ptrans_to_E" = "P→E (trans)",
  "E_to_D"      = "E→D",
  "D_to_E"      = "D→E",
  "Pcis_to_D"   = "P→D (cis)",
  "Ptrans_to_D" = "P→D (trans)",
  "D_to_P"      = "D→P"
)

pair_map <- c(
  "E_to_P"="E-P","Pcis_to_E"="E-P","Ptrans_to_E"="E-P",
  "E_to_D"="E-D","D_to_E"="E-D",
  "Pcis_to_D"="P-D","Ptrans_to_D"="P-D","D_to_P"="P-D"
)

# Option A ordering: forward first, reverse second; cis/trans nested under same direction
edge_levels <- c(
  # E-P: forward first
  "E→P",
  "P→E (cis)",
  "P→E (trans)",
  # E-D
  "E→D",
  "D→E",
  # P-D: forward P→D first (cis, trans), then reverse D→P
  "P→D (cis)",
  "P→D (trans)",
  "D→P"
)

# For ECDF collapse to 6 directions (cis/trans merged)
edge6_label_map <- c(
  "E_to_P"    = "E→P",
  "P_to_E"    = "P→E",
  "E_to_D"    = "E→D",
  "D_to_E"    = "D→E",
  "P_to_D"    = "P→D",
  "D_to_P"    = "D→P"
)

pair6_map <- c(
  "E_to_P"="E-P","P_to_E"="E-P",
  "E_to_D"="E-D","D_to_E"="E-D",
  "P_to_D"="P-D","D_to_P"="P-D"
)

edge6_levels <- c("E→P","P→E","E→D","D→E","P→D","D→P")

# pick the MR effect column if needed (not essential here)
pick_effect_col <- function(dt) {
  cand <- c("b","beta","beta.exposure","b.exposure","effect","estimate")
  hit <- cand[cand %in% names(dt)]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

safe_padj <- function(p) {
  if (all(is.na(p))) return(rep(NA_real_, length(p)))
  p.adjust(p, method = CFG$adj_method)
}

# ============================================================
# Load all edge-type bundles
# ============================================================
PD <- read_res("PD", CFG$rds_dir)
EP <- read_res("EP", CFG$rds_dir)
ED <- read_res("ED", CFG$rds_dir)
DE <- read_res("DE", CFG$rds_dir)
PE <- read_res("PE", CFG$rds_dir)
DP <- read_res("DP", CFG$rds_dir)

mr_all  <- rbindlist(list(PD$mr, EP$mr, ED$mr, DE$mr, PE$mr, DP$mr), fill = TRUE)
het_all <- rbindlist(list(PD$het,EP$het,ED$het,DE$het,PE$het,DP$het), fill = TRUE)
ple_all <- rbindlist(list(PD$ple,EP$ple,ED$ple,DE$ple,PE$ple,DP$ple), fill = TRUE)

# ============================================================
# Standardize + label MR table
# ============================================================

# Keep only IVW + Wald ratio for the "main MR" results used for hit calls
keep_methods <- c("Inverse variance weighted", "Wald ratio")

mr_all <- as.data.table(mr_all)

# edge label + pair label (FIXED vectorized lookup; no recycling)
mr_all[, edge := unname(edge_label_map[edge_dir])]
mr_all[, pair := unname(pair_map[edge_dir])]

# sanity check unmapped
bad <- mr_all[is.na(edge), unique(edge_dir)]
if (length(bad)) stop("Unmapped edge_dir values in MR table: ", paste(bad, collapse=", "))

# filter to main methods
mr_main <- mr_all[method %in% keep_methods]

# compute q-values *within each edge_dir* (you can change to within edge_type if you prefer)
mr_main[, pval_adj := safe_padj(pval), by = edge_dir]

# MR hit definition
mr_main[, mr_hit := !is.na(pval_adj) & (pval_adj < CFG$alpha_q)]

# ============================================================
# Sensitivity tables: IVW heterogeneity and Egger intercept pleiotropy
# ============================================================

# Heterogeneity: file includes both IVW and MR Egger rows; use IVW Cochran Q p-value
het_ivw <- as.data.table(het_all)
if (nrow(het_ivw)) {
  het_ivw <- het_ivw[method == "Inverse variance weighted"]
  # expected columns: edge_dir, src_id, tgt_id, Q_pval
  keep_cols <- intersect(c("edge_dir","src_id","tgt_id","Q_pval"), names(het_ivw))
  het_ivw <- het_ivw[, ..keep_cols]
  setnames(het_ivw, old = "Q_pval", new = "het_pval", skip_absent = TRUE)
}

# Pleiotropy: Egger intercept p-value lives in ple table column "pval"
ple_dt <- as.data.table(ple_all)
if (nrow(ple_dt)) {
  keep_cols <- intersect(c("edge_dir","src_id","tgt_id","pval"), names(ple_dt))
  ple_dt <- ple_dt[, ..keep_cols]
  setnames(ple_dt, old = "pval", new = "egger_pval", skip_absent = TRUE)
}

# Merge sensitivity onto main MR rows (by edge_dir, src_id, tgt_id)
setkeyv(mr_main,  c("edge_dir","src_id","tgt_id"))
if (nrow(het_ivw)) setkeyv(het_ivw, c("edge_dir","src_id","tgt_id"))
if (nrow(ple_dt))  setkeyv(ple_dt,  c("edge_dir","src_id","tgt_id"))

mr_sens <- copy(mr_main)
if (nrow(het_ivw)) mr_sens <- het_ivw[mr_sens]
if (nrow(ple_dt))  mr_sens <- ple_dt[mr_sens]

# Sensitivity pass rule:
# - missing diagnostics => pass (because Wald ratio / low nsnp)
# - otherwise: het_pval >= 0.05 AND egger_pval >= 0.05
mr_sens[, sens_pass :=
          (is.na(het_pval)   | het_pval   >= CFG$alpha_sens) &
          (is.na(egger_pval) | egger_pval >= CFG$alpha_sens)
]

mr_sens[, hit_after_sens := mr_hit & sens_pass]

# ============================================================
# Plot 1: Hit-rate before vs after sensitivity filtering
# (labels as num/denom; avoid overlap by nudging)
# ============================================================

hit_summary <- mr_sens[, .(
  denom = .N,
  n_hit = sum(mr_hit, na.rm=TRUE),
  n_after = sum(hit_after_sens, na.rm=TRUE)
), by = .(pair, edge)]

hit_long <- rbindlist(list(
  hit_summary[, .(pair, edge, stage="MR hits (q<0.05)", n = n_hit, denom)],
  hit_summary[, .(pair, edge, stage="After sensitivity filter", n = n_after, denom)]
), fill = TRUE)

hit_long[, pct := ifelse(denom > 0, 100 * n/denom, NA_real_)]
hit_long[, lab := paste0(n, "/", denom)]

# apply ordering
hit_long[, edge := factor(edge, levels = edge_levels)]
hit_long[, pair := factor(pair, levels = c("E-P","E-D","P-D"))]
hit_long[, stage := factor(stage, levels = c("MR hits (q<0.05)", "After sensitivity filter"))]

# label x position: slightly beyond bar end, with stage-specific nudges
hit_long[, lab_x := pmin(pct + 0.6, max(pct, na.rm=TRUE) + 1.5)]

p_hit_before_after <- ggplot(hit_long, aes(x = pct, y = edge, fill = stage)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(
    aes(x = pct, label = lab),
    position = position_dodge(width = 0.8),
    hjust = -0.15, size = 3
  ) +
  facet_grid(pair ~ ., scales = "free_y", space = "free_y") +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.2), add = c(0, 1))) +
  theme_bw() +
  labs(
    x = "% significant",
    y = NULL,
    title = "MR hit-rate before vs after sensitivity filtering",
    subtitle = "Sensitivity pass: IVW Q_pval≥0.05 and Egger intercept p≥0.05 (missing=pass).\nWald ratio has no Q/Egger."
  ) +
  theme(
    legend.title = element_blank(),
    strip.background = element_rect(fill = NA),
    plot.margin = margin(2, 2, 2, 2, "cm")
  ) +
  coord_cartesian(clip = "off")  # allow labels outside panel

ggsave(file.path(CFG$out_dir, "MR_hit_rate_before_after_sensitivity.png"),
       p_hit_before_after, width = 10, height = 6, dpi = 1000)
ggsave(file.path(CFG$out_dir, "MR_hit_rate_before_after_sensitivity.svg"),
       p_hit_before_after, width = 10, height = 6)

# ============================================================
# Plot 2: Sensitivity flags among MR hits (diagnostic availability + flagged)
# - Denominator: among MR hits where diagnostic exists (Wald ratio excluded automatically)
# - Show num/denom labels
# ============================================================

hits_only <- mr_sens[mr_hit == TRUE]

flag_summary <- rbindlist(list(
  # heterogeneity availability/flag among hits
  hits_only[!is.na(het_pval), .(
    num = sum(het_pval < CFG$alpha_sens, na.rm = TRUE),
    denom = .N
  ), by = .(pair, edge)][, .(pair, edge, diag="Heterogeneity (IVW Q_pval<0.05)", num, denom)],
  
  # egger availability/flag among hits
  hits_only[!is.na(egger_pval), .(
    num = sum(egger_pval < CFG$alpha_sens, na.rm = TRUE),
    denom = .N
  ), by = .(pair, edge)][, .(pair, edge, diag="Directional pleiotropy (Egger p<0.05)", num, denom)]
), fill = TRUE)

flag_summary[, pct := ifelse(denom > 0, 100 * num/denom, NA_real_)]
flag_summary[, lab := paste0(num, "/", denom)]

flag_summary[, edge := factor(edge, levels = edge_levels)]
flag_summary[, pair := factor(pair, levels = c("E-P","E-D","P-D"))]
flag_summary[, diag := factor(diag, levels = c("Directional pleiotropy (Egger p<0.05)",
                                               "Heterogeneity (IVW Q_pval<0.05)"))]

p_flags <- ggplot(flag_summary, aes(x = pct, y = edge, fill = diag)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(
    aes(x = pct, label = lab),
    position = position_dodge(width = 0.8),
    hjust = -0.15, size = 3
  ) +
  facet_grid(pair ~ ., scales = "free_y", space = "free_y") +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.2), add = c(0, 1))) +
  theme_bw() +
  labs(
    x = "% of MR hits flagged (among hits with diagnostic available)",
    y = NULL,
    title = "Sensitivity flags among MR hits",
    subtitle = "Percentages computed among MR hits where the diagnostic exists (Wald ratio excluded from denominators)."
  ) +
  theme(
    legend.title = element_blank(),
    strip.background = element_rect(fill = NA),
    plot.margin = margin(1, 1, 1, 1, "cm")
  ) +
  coord_cartesian(clip = "off")

ggsave(file.path(CFG$out_dir, "MR_sensitivity_flags_among_hits.png"),
       p_flags, width = 10, height = 6, dpi = 1000)
ggsave(file.path(CFG$out_dir, "MR_sensitivity_flags_among_hits.svg"),
       p_flags, width = 10, height = 6)

# ============================================================
# ECDF plots: collapse cis/trans into 6 edge directions
# ============================================================

edge8_levels <- c(
  # E-P (forward first)
  "E→P", "P→E (cis)", "P→E (trans)",
  # E-D
  "E→D", "D→E",
  # P-D (forward first)
  "P→D (cis)", "P→D (trans)",
  # reverse
  "D→P"
)


# Build ECDF dataset for heterogeneity p-values
ecdf_het <- hits_only[!is.na(het_pval), .(edge, pair, pval = het_pval)]
ecdf_het[, edge := factor(edge, levels = edge8_levels)]
ecdf_het[, pair := factor(pair, levels = c("E-P","E-D","P-D"))]


# Build ECDF dataset for Egger p-values
ecdf_ple <- hits_only[!is.na(egger_pval), .(edge, pair, pval = egger_pval)]
ecdf_ple[, edge := factor(edge, levels = edge8_levels)]
ecdf_ple[, pair := factor(pair, levels = c("E-P","E-D","P-D"))]


# Reference (uniform) line
ref_line <- data.table(pval = c(0,1), ecdf = c(0,1))

# ============================================================
# Plot 3: ECDF of heterogeneity p-values (6 facets only)
# Interpretation: curve above diagonal => enrichment of small p (more heterogeneity)
# ============================================================

p_ecdf_het <- ggplot(ecdf_het, aes(x = pval)) +
  stat_ecdf(geom = "step", linewidth = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  facet_wrap(~ edge, ncol = 4, drop = TRUE) +
  coord_cartesian(xlim = c(0,1), ylim = c(0,1)) +
  theme_bw() +
  labs(
    x = "IVW Cochran Q p-value",
    y = "ECDF",
    title = "Heterogeneity p-values (ECDF; among MR hits where available)",
    subtitle = "Dashed line: uniform reference (expected if no heterogeneity enrichment)."
  ) +
  theme(
    strip.background = element_rect(fill = "grey90"),
    panel.grid.minor = element_blank()
  )


ggsave(file.path(CFG$out_dir, "MR_ecdf_heterogeneity.png"),
       p_ecdf_het, width = 14, height = 6, dpi = 400)
ggsave(file.path(CFG$out_dir, "MR_ecdf_heterogeneity.svg"),
       p_ecdf_het, width = 14, height = 6)

# ============================================================
# Plot 4: ECDF of Egger intercept p-values (6 facets only)
# Interpretation: curve above diagonal => enrichment of small p (more pleiotropy)
# ============================================================

p_ecdf_ple <- ggplot(ecdf_ple, aes(x = pval)) +
  stat_ecdf(geom = "step", linewidth = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  facet_wrap(~ edge, ncol = 4, drop = TRUE) +
  coord_cartesian(xlim = c(0,1), ylim = c(0,1)) +
  theme_bw() +
  labs(
    x = "Egger intercept p-value",
    y = "ECDF",
    title = "Directional pleiotropy p-values (ECDF; among MR hits where available)",
    subtitle = "Dashed line: uniform reference (expected if no pleiotropy enrichment)."
  ) +
  theme(
    strip.background = element_rect(fill = "grey90"),
    panel.grid.minor = element_blank()
  )


ggsave(file.path(CFG$out_dir, "MR_ecdf_pleiotropy.png"),
       p_ecdf_ple, width = 14, height = 6, dpi = 400)
ggsave(file.path(CFG$out_dir, "MR_ecdf_pleiotropy.svg"),
       p_ecdf_ple, width = 14, height = 6)

# ============================================================
# Save a compact table for reporting (optional)
# ============================================================
fwrite(hit_summary, file.path(CFG$out_dir, "MR_hit_summary_before_after.tsv"), sep = "\t")
fwrite(flag_summary, file.path(CFG$out_dir, "MR_sensitivity_flag_summary.tsv"), sep = "\t")

message("Done. Plots saved to: ", CFG$out_dir)






