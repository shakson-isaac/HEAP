suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(ggplot2)
  library(ggpmisc)
  library(pbapply)
  library(dplyr)
  library(purrr)
  library(plotly)
  library(htmlwidgets)
  library(readxl)
  library(qs)
  library(ComplexHeatmap)
  library(circlize)
})

# LOAD MR/Assoc/GLP1 data together:
MRres <- fread(file = "/n/groups/patel/shakson_ukb/UK_Biobank/Output/MR_edges/summary/MRmotifs.csv")
#Load INT structure:
setwd("/n/groups/patel/shakson_ukb/UK_Biobank/")
HEAPint <- qread("./Output/HEAPres/HEAPintv2.qs")

#SomaScan vs Olink Comparison:
prot_rel <- fread("/n/groups/patel/IGLOO/UKB/OlinkSoma/OlinkSoma.csv", skip = 3)
prot_rel <- prot_rel %>% select(c("gene_name","olink_nonnorm_corr","olink_smpnorm_corr"))
colnames(prot_rel) <- c("EntrezGeneSymbol", "r_cross", "r_crossv2")
prot_rel <- na.omit(prot_rel)

head(MRres)
head(prot_rel)
head(HEAPint@sList$Model1, 5)


##### ScatterPlot Plotting: #####


# -----------------------------
# 0) Extract HEAP Model1 to a DT
# -----------------------------
heap <- as.data.table(HEAPint@sList$Model6)
setnames(heap,
         old = c("ID","EntrezGeneSymbol","Estimate","Std. Error"),
         new = c("Exposure","Protein","beta_HEAP","se_HEAP"),
         skip_absent = TRUE)

# Keep intervention columns if present (you showed these in head())
# HERITAGE_effect, HERITAGE_se, GLP1_effect1, GLP1_se1, GLP1_effect2, GLP1_se2
heap[, beta_GLP1_1 := GLP1_effect1]
heap[, se_GLP1_1   := GLP1_se1]
heap[, beta_GLP1_2 := GLP1_effect2]
heap[, se_GLP1_2   := GLP1_se2]
heap[, beta_HERIT  := HERITAGE_effect]
heap[, se_HERIT    := HERITAGE_se]

# -----------------------------
# 1) Prep Soma/Olink cross-platform reliability
# -----------------------------
prot_rel_dt <- as.data.table(prot_rel)
setnames(prot_rel_dt, "EntrezGeneSymbol", "Protein")

# -----------------------------
# 2) Prep MR: choose a single "best" PD effect per row (cis vs trans)
#    (pick the one with smaller adjusted p; fallback to non-NA)
# -----------------------------
MR <- as.data.table(MRres)

MR[, `:=`(
  padj_PD_best = fifelse(!is.na(padj_PDcis) & !is.na(padj_PDtrans),
                         pmin(padj_PDcis, padj_PDtrans),
                         fifelse(!is.na(padj_PDcis), padj_PDcis, padj_PDtrans)),
  beta_PD_best = fifelse(!is.na(padj_PDcis) & !is.na(padj_PDtrans),
                         fifelse(padj_PDcis <= padj_PDtrans, beta_PDcis, beta_PDtrans),
                         fifelse(!is.na(beta_PDcis), beta_PDcis, beta_PDtrans)),
  se_PD_best   = fifelse(!is.na(se_PDcis) & !is.na(se_PDtrans),
                         fifelse(padj_PDcis <= padj_PDtrans, se_PDcis, se_PDtrans),
                         fifelse(!is.na(se_PDcis), se_PDcis, se_PDtrans)),
  PD_best_type = fifelse(!is.na(padj_PDcis) & !is.na(padj_PDtrans),
                         fifelse(padj_PDcis <= padj_PDtrans, "cis", "trans"),
                         fifelse(!is.na(padj_PDcis), "cis", fifelse(!is.na(padj_PDtrans), "trans", NA_character_)))
)]

MR[, `:=`(
  sig_EP = !is.na(padj_EP) & padj_EP < 0.05,
  sig_PD = !is.na(padj_PD_best) & padj_PD_best < 0.05,
  sig_ED = !is.na(padj_ED) & padj_ED < 0.05,
  sig_DP = !is.na(padj_DP) & padj_DP < 0.05
)]

# -----------------------------
# 3) Join HEAP + MR + platform reliability
#    NOTE: MR is triplet-level; HEAP is exposure-protein level but has GLP1/HERIT per protein.
# -----------------------------
df <- merge(MR, heap, by = c("Exposure","Protein"), all.x = TRUE)
df <- merge(df, prot_rel_dt, by = "Protein", all.x = TRUE, allow.cartesian = TRUE)

# Convenience: a "reliability weight" you can use for size/alpha
df[, rel_w := pmax(0, pmin(1, abs(r_crossv2)))]  # clamp to [0,1]

# Optional: focus to rows where you actually have HEAP + GLP1 info
df_glp1 <- df[!is.na(beta_HEAP) & !is.na(beta_GLP1_1)]

# -----------------------------
# Helper: concordance flags (directional agreement)
# -----------------------------
df_glp1[, conc_HEAP_GLP1 := sign(beta_HEAP) == sign(beta_GLP1_1)]
df_glp1[, conc_MR_EP_HEAP := !is.na(beta_EP) & (sign(beta_EP) == sign(beta_HEAP))]

# ============================================================
# PLOT 1: HEAP (x) vs GLP1 (y) with MR annotations
# - Color: MR EP significance
# - Shape: MR PD significance (best of cis/trans)
# - Alpha: platform reliability (Soma/Olink)
# - Label: top outliers or proteins of interest
# ============================================================
# Choose a manageable subset for clarity (e.g., one exposure, or top-N by abs(beta_HEAP))
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
})

#-----------------------------
# Helper: format p-values nicely
#-----------------------------
fmt_p <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 1e-300) return("<1e-300")
  format.pval(p, digits = 2, eps = 1e-300)
}

#-----------------------------
# Main plotting function
#-----------------------------

plot_HEAP_GLP1_MR_onepanel <- function(
    EXPOSURE_TO_PLOT,
    arm = c("GLP1_1", "GLP1_2", "HERITAGE"),
    disease_for_arm = NULL,   # optional; if NULL we don't filter MR by Disease
    mr_alpha = 0.05,
    
    # text controls
    title = NULL,
    subtitle = NULL,
    xlab = "HEAP beta (Exposure \u2192 Protein)",
    ylab = NULL,
    show_title = TRUE,
    show_subtitle = TRUE,
    show_axis_titles = TRUE,
    
    # correlation annotation controls
    corr_model = "Model6",
    corr_loc = c("topleft","topright","bottomleft","bottomright"),
    corr_digits = 2,
    
    # label controls
    label_n = 12,
    label_only_mr = FALSE,
    
    # legend / compactness
    legend_position = c("right","bottom","none"),
    base_size = 12
) {
  
  # -----------------------------
  # Helpers
  # -----------------------------
  fmt_p <- function(p) {
    if (is.na(p)) return("NA")
    if (p < 1e-300) return("<1e-300")
    format.pval(p, digits = 2, eps = 1e-300)
  }
  canon_exposure <- function(x) sub("_[0-9]+$", "", x)
  cap_n <- function(x, n) x[seq_len(min(length(x), n))]
  
  # -----------------------------
  # Arg parsing
  # -----------------------------
  arm <- match.arg(arm)
  corr_loc <- match.arg(corr_loc)
  legend_position <- match.arg(legend_position)
  
  # Map arm -> column in HEAP tables
  y_col <- switch(
    arm,
    "GLP1_1"    = "GLP1_effect1",
    "GLP1_2"    = "GLP1_effect2",
    "HERITAGE"  = "HERITAGE_effect"
  )
  
  if (is.null(title)) title <- paste0("HEAP vs ", arm, " protein shifts (MR edge encoded)")
  if (is.null(ylab))  ylab  <- y_col
  
  exposure_key_plot <- canon_exposure(EXPOSURE_TO_PLOT)
  
  # -----------------------------
  # 1) HEAP subset (relaxed exposure match)
  # -----------------------------
  heap_dt <- data.table::as.data.table(HEAPint@sList[[corr_model]])
  if (!all(c("ID","EntrezGeneSymbol") %in% names(heap_dt))) {
    stop("HEAPint@sList[[corr_model]] must contain ID and EntrezGeneSymbol.")
  }
  
  data.table::setnames(
    heap_dt,
    old = c("ID","EntrezGeneSymbol","Estimate","Std. Error"),
    new = c("Exposure","Protein","beta_HEAP","se_HEAP"),
    skip_absent = TRUE
  )
  
  if (!y_col %in% names(heap_dt)) {
    stop("Column ", y_col, " not found in HEAPint@sList[[", corr_model, "]].")
  }
  
  heap_dt[, Exposure_key := canon_exposure(Exposure)]
  
  heap_pick <- if (any(heap_dt$Exposure == EXPOSURE_TO_PLOT, na.rm = TRUE)) {
    heap_dt[Exposure == EXPOSURE_TO_PLOT]
  } else {
    heap_dt[Exposure_key == exposure_key_plot]
  }
  
  # Collapse to 1 row per protein
  heap_sub <- heap_pick[
    , .(
      Exposure  = EXPOSURE_TO_PLOT,
      beta_HEAP = mean(beta_HEAP, na.rm = TRUE),
      beta_arm  = mean(get(y_col), na.rm = TRUE)
    ),
    by = .(Protein)
  ][!is.na(beta_arm)]
  
  if (nrow(heap_sub) == 0) {
    stop("No HEAP rows for exposure with available ", arm, " effects after relaxed matching.")
  }
  
  # -----------------------------
  # 2) MR subset (optional disease filter) + merge
  # -----------------------------
  MR <- data.table::as.data.table(MRres)
  MR[, Exposure_key := canon_exposure(Exposure)]
  
  mr_sub <- MR[Exposure_key == exposure_key_plot]
  if (!is.null(disease_for_arm) && "Disease" %in% names(mr_sub)) {
    mr_sub <- mr_sub[Disease == disease_for_arm]
  }
  
  # Merge by Protein (HEAP defines the plotted set)
  dt <- merge(heap_sub, mr_sub, by = "Protein", all.x = TRUE, allow.cartesian = TRUE)
  
  # -----------------------------
  # 3) Soma/Olink reliability -> size/alpha
  # -----------------------------
  prot_rel_dt <- data.table::as.data.table(prot_rel)
  if ("EntrezGeneSymbol" %in% names(prot_rel_dt)) {
    data.table::setnames(prot_rel_dt, "EntrezGeneSymbol", "Protein", skip_absent = TRUE)
  }
  if (!"r_crossv2" %in% names(prot_rel_dt) && "r_cross" %in% names(prot_rel_dt)) {
    prot_rel_dt[, r_crossv2 := r_cross]
  }
  
  dt <- merge(dt, prot_rel_dt[, .(Protein, r_crossv2)], by = "Protein", all.x = TRUE)
  
  dt[, olink_soma_r := r_crossv2]
  dt[, size_r  := pmax(0, olink_soma_r)]
  dt[is.na(size_r), size_r := 0]
  dt[, alpha_r := data.table::fifelse(!is.na(olink_soma_r) & olink_soma_r > 0, 0.9, 0.05)]
  
  # -----------------------------
  # 4) MR edge significance (safe if columns missing)
  # -----------------------------
  has <- function(x) x %in% names(dt)
  
  dt[, sig_PDcis   := if (has("padj_PDcis"))   !is.na(padj_PDcis)   & padj_PDcis   < mr_alpha else FALSE]
  dt[, sig_PDtrans := if (has("padj_PDtrans")) !is.na(padj_PDtrans) & padj_PDtrans < mr_alpha else FALSE]
  dt[, sig_DP      := if (has("padj_DP"))      !is.na(padj_DP)      & padj_DP      < mr_alpha else FALSE]
  
  dt[, mr_edge_sig := data.table::fifelse(
    sig_PDcis, "PDcis",
    data.table::fifelse(sig_PDtrans, "PDtrans",
                        data.table::fifelse(sig_DP, "DP", "None"))
  )]
  dt[, mr_edge_sig := factor(mr_edge_sig, levels = c("None","PDcis","PDtrans","DP"))]
  
  # -----------------------------
  # 5) Labels
  # -----------------------------
  pd_cis_prots   <- unique(dt[mr_edge_sig=="PDcis"][order(-abs(beta_arm))]$Protein)
  pd_trans_prots <- unique(dt[mr_edge_sig=="PDtrans"][order(-abs(beta_arm))]$Protein)
  dp_prots       <- unique(dt[mr_edge_sig=="DP"][order(-abs(beta_arm))]$Protein)
  
  n_cis   <- ceiling(label_n * 0.5)
  n_trans <- ceiling(label_n * 0.3)
  n_dp    <- ceiling(label_n * 0.2)
  
  lab_mr <- unique(c(
    cap_n(pd_cis_prots, n_cis),
    cap_n(pd_trans_prots, n_trans),
    cap_n(dp_prots, n_dp)
  ))
  
  n_left <- max(0, label_n - length(lab_mr))
  fill_prots <- unique(c(
    dt[order(-abs(beta_HEAP))]$Protein,
    dt[order(-abs(beta_arm))]$Protein
  ))
  
  lab_fill  <- setdiff(fill_prots, lab_mr)
  lab_prots <- unique(c(lab_mr, cap_n(lab_fill, n_left)))
  
  if (isTRUE(label_only_mr)) lab_prots <- cap_n(lab_mr, label_n)
  
  dt[, label_once := (Protein %in% lab_prots) & !duplicated(Protein)]
  
  # -----------------------------
  # 6) Correlation annotation (uses same arm column)
  # -----------------------------
  ctab <- data.table::as.data.table(HEAPint@cList[[corr_model]], keep.rownames = "Exposure")
  ptab <- data.table::as.data.table(HEAPint@pList[[corr_model]], keep.rownames = "Exposure")
  ctab[, Exposure_key := canon_exposure(Exposure)]
  ptab[, Exposure_key := canon_exposure(Exposure)]
  
  if (!y_col %in% names(ctab) || !y_col %in% names(ptab)) {
    r_val <- NA_real_
    p_val <- NA_real_
  } else if (any(ctab$Exposure == EXPOSURE_TO_PLOT, na.rm = TRUE)) {
    r_val <- ctab[Exposure == EXPOSURE_TO_PLOT, get(y_col)]
    p_val <- ptab[Exposure == EXPOSURE_TO_PLOT, get(y_col)]
  } else {
    r_val <- ctab[Exposure_key == exposure_key_plot, get(y_col)]
    p_val <- ptab[Exposure_key == exposure_key_plot, get(y_col)]
  }
  
  r_txt <- ifelse(is.na(r_val), "NA", formatC(r_val, digits = corr_digits, format = "f"))
  p_txt <- fmt_p(p_val)
  
  corr_label <- paste0(arm, " corr: r = ", r_txt, "\n", "p = ", p_txt)
  
  if (is.null(subtitle)) {
    subtitle <- if (is.null(disease_for_arm)) {
      paste0("Exposure: ", EXPOSURE_TO_PLOT)
    } else {
      paste0("Exposure: ", EXPOSURE_TO_PLOT, " | Disease: ", disease_for_arm)
    }
  }
  
  # -----------------------------
  # 7) Annotation location
  # -----------------------------
  xr <- range(dt$beta_HEAP, na.rm = TRUE)
  yr <- range(dt$beta_arm,  na.rm = TRUE)
  xpad <- diff(xr) * 0.02
  ypad <- diff(yr) * 0.04
  
  ann_x <- switch(corr_loc,
                  "topleft"     = xr[1] + xpad,
                  "bottomleft"  = xr[1] + xpad,
                  "topright"    = xr[2] - xpad,
                  "bottomright" = xr[2] - xpad
  )
  ann_y <- switch(corr_loc,
                  "topleft"     = yr[2] - ypad,
                  "topright"    = yr[2] - ypad,
                  "bottomleft"  = yr[1] + ypad,
                  "bottomright" = yr[1] + ypad
  )
  ann_hjust <- if (grepl("right", corr_loc)) 1 else 0
  ann_vjust <- if (grepl("bottom", corr_loc)) 0 else 1
  
  # -----------------------------
  # 8) Plot
  # -----------------------------
  p <- ggplot2::ggplot(dt, ggplot2::aes(x = beta_HEAP, y = beta_arm)) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.25, color = "grey75") +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.25, color = "grey75") +
    ggplot2::geom_point(ggplot2::aes(color = mr_edge_sig, size = size_r, alpha = alpha_r)) +
    
    ggrepel::geom_text_repel(
      data = dt[label_once == TRUE],
      ggplot2::aes(label = Protein),
      size = 3,
      box.padding = 0.2,
      point.padding = 0.12,
      min.segment.length = 0,
      max.overlaps = 60
    ) +
    
    ggplot2::annotate(
      "label",
      x = ann_x, y = ann_y,
      label = corr_label,
      hjust = ann_hjust, vjust = ann_vjust,
      size = 3.2,
      label.size = 0.25,
      fill = "white"
    ) +
    
    ggplot2::scale_alpha_identity(guide = "none") +
    ggplot2::scale_size_continuous(
      range = c(1.6, 5.2),
      breaks = c(0.25, 0.5, 0.75),
      limits = c(0, 1),
      name = "SomaScan–Olink\nr (pos only)"
    ) +
    ggplot2::scale_color_manual(
      values = c(
        "None"    = "grey75",
        "PDcis"   = "#1b9e77",
        "PDtrans" = "#7570b3",
        "DP"      = "#d95f02"
      ),
      name = paste0("MR significant edge\n(adj.p<", mr_alpha, ")")
    ) +
    ggplot2::labs(
      title = if (show_title) title else NULL,
      subtitle = if (show_subtitle) subtitle else NULL,
      x = if (show_axis_titles) xlab else NULL,
      y = if (show_axis_titles) ylab else NULL
    ) +
    ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      legend.position = if (legend_position == "none") "none" else legend_position,
      legend.title = ggplot2::element_text(face = "bold"),
      legend.key.height = grid::unit(0.8, "lines"),
      legend.spacing.y = grid::unit(0.2, "lines"),
      plot.title.position = "plot",
      plot.margin = ggplot2::margin(6, 6, 6, 6),
      plot.subtitle = ggplot2::element_text(size = base_size - 1)
    ) +
    ggplot2::coord_cartesian(clip = "off")
  
  return(p)
}




#-----------------------------
# Example usage
#-----------------------------
EXPOSURE_TO_PLOT <- "types_of_physical_activity_in_last_4_weeks_f6164_0_0.multi_Strenuous_sports"

p0 <- plot_HEAP_GLP1_MR_onepanel(
  EXPOSURE_TO_PLOT = "types_of_physical_activity_in_last_4_weeks_f6164_0_0.multi_Strenuous_sports",
  arm = "HERITAGE",
  disease_for_arm = "finngen_R12_E4_OBESITY",
  title = "Exercise → proteins vs HERITAGE (Endurance Exercise)",
  subtitle = "Strenuous sports (UKB) | MR edges vs Obesity",
  xlab = "UKB assoc (E→P) beta",
  ylab = "HERITAGE protein shift",
  legend_position = "right",
  corr_loc = "topleft"
)

p1 <- plot_HEAP_GLP1_MR_onepanel(
  EXPOSURE_TO_PLOT = EXPOSURE_TO_PLOT,
  arm = "GLP1_1",
  disease_for_arm = "finngen_R12_E4_OBESITY",
  title = "Exercise → proteins vs GLP1_1 (obesity cohort)",
  subtitle = "Strenuous sports (UKB) | MR edges vs obesity",
  xlab = "UKB assoc (E→P) beta",
  ylab = "GLP1_1 protein shift",
  legend_position = "right",
  corr_loc = "topleft"
)

p2 <- plot_HEAP_GLP1_MR_onepanel(
  EXPOSURE_TO_PLOT = EXPOSURE_TO_PLOT,
  arm = "GLP1_2",
  disease_for_arm = "finngen_R12_T2D",
  title = "Exercise → proteins vs GLP1_2 (diabetes-inclusive)",
  subtitle = "Strenuous sports (UKB) | MR edges vs T2D",
  xlab = "UKB assoc (E→P) beta",
  ylab = "GLP1_2 protein shift",
  legend_position = "right",
  corr_loc = "topleft"
)

print(p0)
print(p1)
print(p2)

ggsave("/n/groups/patel/shakson_ukb/UK_Biobank/Output/MRplots/MRHERITAGE_StrenSports.png",
       plot=p0, dpi = 1000,
       width = 7, height = 4, unit = "in")
ggsave("/n/groups/patel/shakson_ukb/UK_Biobank/Output/MRplots/MRGLP1_STEP1_StrenSports.png",
       plot=p1, dpi = 1000,
       width = 7, height = 4, unit = "in")
ggsave("/n/groups/patel/shakson_ukb/UK_Biobank/Output/MRplots/MRGLP1_STEP2_StrenSports.png",
       plot=p2, dpi = 1000,
       width = 7, height = 4, unit = "in")


p3 <- plot_HEAP_GLP1_MR_onepanel(
  EXPOSURE_TO_PLOT = "fresh_fruit_intake_f1309_0_0",
  arm = "GLP1_1",
  disease_for_arm = "finngen_R12_E4_OBESITY",
  title = "Diet → proteins vs GLP1_1 (obesity cohort)",
  subtitle = "Fruit Intake (UKB) | MR edges vs obesity",
  xlab = "UKB assoc (E→P) beta",
  ylab = "GLP1_1 protein shift",
  legend_position = "right",
  corr_loc = "topleft"
)

p4 <- plot_HEAP_GLP1_MR_onepanel(
  EXPOSURE_TO_PLOT = "fresh_fruit_intake_f1309_0_0",
  arm = "GLP1_2",
  disease_for_arm = "finngen_R12_T2D",
  title = "Diet → proteins vs GLP1_2 (diabetes-inclusive)",
  subtitle = "Fruit Intake (UKB) | MR edges vs T2D",
  xlab = "UKB assoc (E→P) beta",
  ylab = "GLP1_2 protein shift",
  legend_position = "right",
  corr_loc = "topleft"
)

print(p3)
print(p4)

ggsave("/n/groups/patel/shakson_ukb/UK_Biobank/Output/MRplots/MRGLP1_STEP1_FruitIntake.png",
       plot=p3, dpi = 1000,
       width = 7, height = 4, unit = "in")
ggsave("/n/groups/patel/shakson_ukb/UK_Biobank/Output/MRplots/MRGLP1_STEP2_FruitIntake.png",
       plot=p4, dpi = 1000,
       width = 7, height = 4, unit = "in")


p5 <- plot_HEAP_GLP1_MR_onepanel(
  EXPOSURE_TO_PLOT = "past_tobacco_smoking_f1249_0_0",
  arm = "GLP1_1",
  disease_for_arm = "finngen_R12_E4_OBESITY",
  title = "Quit Smoking → proteins vs GLP1_1 (obesity cohort)",
  subtitle = "Former Smoker | MR edges vs obesity",
  xlab = "UKB assoc (E→P) beta",
  ylab = "GLP1_1 protein shift",
  legend_position = "right",
  corr_loc = "topleft"
)

p6 <- plot_HEAP_GLP1_MR_onepanel(
  EXPOSURE_TO_PLOT = "past_tobacco_smoking_f1249_0_04",
  arm = "GLP1_2",
  disease_for_arm = "finngen_R12_T2D",
  title = "Quit Smoking → proteins vs GLP1_2 (diabetes-inclusive)",
  subtitle = "Former Smoker | MR edges vs T2D",
  xlab = "UKB assoc (E→P) beta",
  ylab = "GLP1_2 protein shift",
  legend_position = "right",
  corr_loc = "topleft"
)

print(p5)
print(p6)

ggsave("/n/groups/patel/shakson_ukb/UK_Biobank/Output/MRplots/MRGLP1_STEP1_QuitSmoking.png",
       plot=p5, dpi = 1000,
       width = 7, height = 4, unit = "in")
ggsave("/n/groups/patel/shakson_ukb/UK_Biobank/Output/MRplots/MRGLP1_STEP2_QuitSmoking.png",
       plot=p6, dpi = 1000,
       width = 7, height = 4, unit = "in")




## PRIOR ISSUE:
# Past Tobacco smoking wasn't merging well with MRres b/c the specific
# ordinal was turned into a continuous variable.
# FIXED THE ABOVE!!
# Smoking and Disease Context doesnt help as much SINCE the smoking protein seems to e CXCL17

#exposurelist <- sort(unique(df_glp1$Exposure))


sort(unique(sig6$Disease))

xxx <- MRres %>% filter(grepl("past_tobacco",Exposure))

p7 <- plot_HEAP_GLP1_MR_onepanel(
  EXPOSURE_TO_PLOT = "past_tobacco_smoking_f1249_0_04",
  GLP_ARM = "GLP1_1",
  disease_for_arm = "finngen_R12_SMOKING_DEPEND",
  title = "Quit Smoking → proteins vs GLP1_1 (obesity cohort)",
  subtitle = "Former Smoker | MR edges vs smoking depend",
  xlab = "UKB assoc (E→P) beta",
  ylab = "GLP1_1 protein shift",
  legend_position = "right",
  corr_source = "GLP1",
  corr_loc = "topleft"
)

p8 <- plot_HEAP_GLP1_MR_onepanel(
  EXPOSURE_TO_PLOT = "past_tobacco_smoking_f1249_0_04",
  GLP_ARM = "GLP1_2",
  disease_for_arm = "finngen_R12_SMOKING_DEPEND",
  title = "Quit Smoking → proteins vs GLP1_2 (diabetes-inclusive cohort)",
  subtitle = "Former Smoker | MR edges vs smoking depend",
  xlab = "UKB assoc (E→P) beta",
  ylab = "GLP1_2 protein shift",
  legend_position = "right",
  corr_source = "GLP1",
  corr_loc = "topleft"
)

#xxx <- MRres %>% filter(Disease == "finngen_R12_SMOKING_DEPEND")
xxx <- MRres %>% filter(Disease == "finngen_R12_J10_EMPHYSEMA")

print(p7)
print(p8)

ggsave("/n/groups/patel/shakson_ukb/UK_Biobank/Output/MRplots/MRGLP1_STEP1_QuitSmoking_DiffDisease.png",
       plot=p1, dpi = 1000,
       width = 6, height = 4, unit = "in")
ggsave("/n/groups/patel/shakson_ukb/UK_Biobank/Output/MRplots/MRGLP1_STEP2_QuitSmoking_DiffDisease.png",
       plot=p2, dpi = 1000,
       width = 6, height = 4, unit = "in")

sort(unique(MRres$Exposure))
sort(unique(MRres$Disease))

p9 <- plot_HEAP_GLP1_MR_onepanel(
  EXPOSURE_TO_PLOT = "alcohol_intake_versus_10_years_previously_f1628_0_03",
  GLP_ARM = "GLP1_1",
  disease_for_arm = "finngen_R12_J10_EMPHYSEMA",
  title = "Quit Smoking → proteins vs GLP1_1 (obesity cohort)",
  subtitle = "Former Smoker | MR edges vs obesity",
  xlab = "UKB assoc (E→P) beta",
  ylab = "GLP1_1 protein shift",
  legend_position = "right",
  corr_source = "GLP1",
  corr_loc = "topleft"
)

p10 <- plot_HEAP_GLP1_MR_onepanel(
  EXPOSURE_TO_PLOT = "alcohol_intake_versus_10_years_previously_f1628_0_03",
  GLP_ARM = "GLP1_2",
  disease_for_arm = "finngen_R12_J10_EMPHYSEMA",
  title = "Quit Smoking → proteins vs GLP1_2 (diabetes-inclusive cohort)",
  subtitle = "Former Smoker | MR edges vs T2D",
  xlab = "UKB assoc (E→P) beta",
  ylab = "GLP1_2 protein shift",
  legend_position = "right",
  corr_source = "GLP1",
  corr_loc = "topleft"
)

p9
p10


p11 <- plot_HEAP_GLP1_MR_onepanel(
  EXPOSURE_TO_PLOT = "types_of_physical_activity_in_last_4_weeks_f6164_0_0.multi_Strenuous_sports",
  arm = "HERITAGE",
  disease_for_arm = "finngen_R12_E4_OBESITY",
  title = "Exercise → proteins vs HERITAGE (Endurance Exercise)",
  subtitle = "Strenuous sports (UKB) | MR edges vs Obesity",
  xlab = "UKB assoc (E→P) beta",
  ylab = "HERITAGE protein shift",
  legend_position = "right",
  corr_loc = "topleft"
)

print(p11)

View(HEAPint@cList$Model6)
View(HEAPint@pList$Model6)


