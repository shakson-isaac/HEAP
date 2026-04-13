plot_HEAP_GLP1_MR_onepanel <- function(
    EXPOSURE_TO_PLOT,
    GLP_ARM = c("GLP1_1", "GLP1_2"),
    disease_for_arm = NULL,
    mr_alpha = 0.05,

    # text controls
    title = "HEAP vs GLP1 protein shifts (MR edge encoded)",
    subtitle = NULL,
    xlab = "HEAP beta (Exposure → Protein)",
    ylab = NULL,
    show_title = TRUE,
    show_subtitle = TRUE,
    show_axis_titles = TRUE,

    # correlation annotation controls
    corr_model = "Model6",
    corr_source = c("GLP1", "HERITAGE"),   # which precomputed correlation to annotate
    corr_loc = c("topleft","topright","bottomleft","bottomright"),
    corr_digits = 2,

    # label controls
    label_n = 12,
    label_only_mr = FALSE,

    # legend / compactness
    legend_position = c("right","bottom","none"),
    base_size = 12
) {

  GLP_ARM <- match.arg(GLP_ARM)
  corr_source <- match.arg(corr_source)
  corr_loc <- match.arg(corr_loc)
  legend_position <- match.arg(legend_position)

  if (is.null(disease_for_arm)) {
    disease_for_arm <- if (GLP_ARM == "GLP1_1") "finngen_R12_E4_OBESITY" else "finngen_R12_T2D"
  }

  glp_beta_col <- if (GLP_ARM == "GLP1_1") "GLP1_effect1" else "GLP1_effect2"
  if (is.null(ylab)) ylab <- paste0("GLP1 effect (", GLP_ARM, ")")

  #-----------------------------
  # 1) Build HEAP subset (already curated upstream for signif)
  #-----------------------------
  heap_dt <- as.data.table(HEAPint@sList[[corr_model]])
  if (!all(c("ID","EntrezGeneSymbol") %in% names(heap_dt))) {
    stop("HEAPint@sList[[modelType]] must contain ID and EntrezGeneSymbol.")
  }
  setnames(heap_dt,
           old = c("ID","EntrezGeneSymbol","Estimate","Std. Error"),
           new = c("Exposure","Protein","beta_HEAP","se_HEAP"),
           skip_absent = TRUE)

  heap_sub <- heap_dt[Exposure == EXPOSURE_TO_PLOT, .(
    Exposure, Protein,
    beta_HEAP,
    beta_GLP1 = get(glp_beta_col)
  )][!is.na(beta_GLP1)]

  if (nrow(heap_sub) == 0) stop("No HEAP rows for exposure with ", GLP_ARM, " available.")

  #-----------------------------
  # 2) Merge MR (restricted to ONE disease for interpretability)
  #-----------------------------
  MR <- as.data.table(MRres)
  mr_sub <- MR[Exposure == EXPOSURE_TO_PLOT & Disease == disease_for_arm]

  dt <- merge(heap_sub, mr_sub, by = c("Exposure","Protein"), all.x = TRUE)

  #-----------------------------
  # 3) Merge Soma–Olink correlation + size metric
  #-----------------------------
  prot_rel_dt <- as.data.table(prot_rel)
  # normalize column name
  if ("EntrezGeneSymbol" %in% names(prot_rel_dt)) setnames(prot_rel_dt, "EntrezGeneSymbol", "Protein", skip_absent = TRUE)
  if (!"r_crossv2" %in% names(prot_rel_dt)) {
    # fallback to r_cross if needed
    if ("r_cross" %in% names(prot_rel_dt)) prot_rel_dt[, r_crossv2 := r_cross]
  }

  dt <- merge(dt, prot_rel_dt[, .(Protein, r_crossv2)], by = "Protein", all.x = TRUE)
  dt[, olink_soma_absr := pmax(0, pmin(1, abs(r_crossv2)))]

  # r_crossv2 is SomaScan–Olink correlation (can be negative)
  dt[, olink_soma_r := r_crossv2]

  # size: use r directly, but only positive values contribute to size
  # (negative/0 -> size 0)
  dt[, size_r := pmax(0, olink_soma_r)]
  dt[is.na(size_r), size_r := 0]

  # alpha: make r<=0 (or NA) extremely faint
  dt[, alpha_r := fifelse(!is.na(olink_soma_r) & olink_soma_r > 0, 0.9, 0.01)]


  #-----------------------------
  # 4) MR edge significance + causal-priority assignment
  # PDcis > PDtrans > DP > None
  #-----------------------------
  dt[, sig_PDcis   := !is.na(padj_PDcis)   & padj_PDcis   < mr_alpha]
  dt[, sig_PDtrans := !is.na(padj_PDtrans) & padj_PDtrans < mr_alpha]
  dt[, sig_DP      := !is.na(padj_DP)      & padj_DP      < mr_alpha]

  dt[, mr_edge_sig := fifelse(
    sig_PDcis, "PDcis",
    fifelse(sig_PDtrans, "PDtrans",
            fifelse(sig_DP, "DP", "None"))
  )]
  dt[, mr_edge_sig := factor(mr_edge_sig, levels = c("None","PDcis","PDtrans","DP"))]

  # optional multi-hit info (useful for labeling)
  dt[, n_sig_edges := as.integer(sig_PDcis) + as.integer(sig_PDtrans) + as.integer(sig_DP)]
  print(dt)
  #-----------------------------
  # 5) Labels: compact + sane default
  #-----------------------------
  # if (label_only_mr) {
  #   lab_candidates <- dt[mr_edge_sig != "None"]
  # } else {
  #   lab_candidates <- dt
  # }
  #
  # lab_prots <- unique(c(
  #   lab_candidates[mr_edge_sig != "None"][1:min(label_n, .N), Protein],
  #   dt[order(-abs(beta_HEAP))][1:ceiling(label_n/2), Protein],
  #   dt[order(-abs(beta_GLP1))][1:ceiling(label_n/2), Protein]
  # ))
  # dt[, label_me := Protein %in% lab_prots]

  #-----------------------------
  # 5) Labels: prioritize by MR edge category (PDcis > PDtrans > DP > None)
  #-----------------------------

  # helper: safely cap n
  cap_n <- function(x, n) x[seq_len(min(length(x), n))]

  # (A) always label MR-significant proteins first, in priority order
  pd_cis_prots   <- dt[mr_edge_sig == "PDcis",   unique(Protein)]
  pd_trans_prots <- dt[mr_edge_sig == "PDtrans", unique(Protein)]
  dp_prots       <- dt[mr_edge_sig == "DP",      unique(Protein)]

  # Optionally: within each category, order by something sensible
  # e.g., strongest GLP1 shift among that category
  pd_cis_prots   <- dt[mr_edge_sig=="PDcis"][order(-abs(beta_GLP1))]$Protein |> unique()
  pd_trans_prots <- dt[mr_edge_sig=="PDtrans"][order(-abs(beta_GLP1))]$Protein |> unique()
  dp_prots       <- dt[mr_edge_sig=="DP"][order(-abs(beta_GLP1))]$Protein |> unique()

  # allocate label slots across categories (you can tune these)
  # Example: 50% PDcis, 30% PDtrans, 20% DP (then fill remainder)
  n_cis   <- ceiling(label_n * 0.5)
  n_trans <- ceiling(label_n * 0.3)
  n_dp    <- ceiling(label_n * 0.2)

  lab_mr <- unique(c(
    cap_n(pd_cis_prots, n_cis),
    cap_n(pd_trans_prots, n_trans),
    cap_n(dp_prots, n_dp)
  ))

  # (B) fill remaining slots with extremes (HEAP/GLP1), regardless of MR
  n_left <- max(0, label_n - length(lab_mr))

  fill_prots <- unique(c(
    dt[order(-abs(beta_HEAP))]$Protein,
    dt[order(-abs(beta_GLP1))]$Protein
  ))

  lab_fill <- setdiff(fill_prots, lab_mr)
  lab_prots <- unique(c(lab_mr, cap_n(lab_fill, n_left)))

  dt[, label_me := Protein %in% lab_prots]

  # Label One Time Only
  dt[, label_once := label_me & !duplicated(Protein)]


  #-----------------------------
  # 6) Pull precomputed correlation + p-value for this exposure
  # from HEAPint@cList / pList, Model6
  #-----------------------------
  ctab <- as.data.table(HEAPint@cList[[corr_model]], keep.rownames = "Exposure")
  ptab <- as.data.table(HEAPint@pList[[corr_model]], keep.rownames = "Exposure")

  # column depends on corr_source + which GLP arm
  # your cList/pList columns are named: GLP1_effect1, GLP1_effect2, HERITAGE_effect
  corr_col <- if (corr_source == "HERITAGE") "HERITAGE_effect" else glp_beta_col

  r_val <- ctab[Exposure == EXPOSURE_TO_PLOT, get(corr_col)]
  p_val <- ptab[Exposure == EXPOSURE_TO_PLOT, get(corr_col)]

  r_txt <- ifelse(is.na(r_val), "NA", formatC(r_val, digits = corr_digits, format = "f"))
  p_txt <- fmt_p(p_val)

  corr_label <- paste0(
    corr_source, " corr: r = ", r_txt, "\n",
    "p = ", p_txt
  )

  # subtitle default
  if (is.null(subtitle)) {
    subtitle <- paste0("Exposure: ", EXPOSURE_TO_PLOT, " | Disease: ", disease_for_arm)
  }

  #-----------------------------
  # 7) Compact plot theme + correlation placement
  #-----------------------------
  # choose annotation coordinates based on data range
  xr <- range(dt$beta_HEAP, na.rm = TRUE)
  yr <- range(dt$beta_GLP1, na.rm = TRUE)
  xpad <- diff(xr) * 0.02
  ypad <- diff(yr) * 0.04

  ann_x <- switch(corr_loc,
                  "topleft"     = xr[1] + xpad,
                  "bottomleft"  = xr[1] + xpad,
                  "topright"    = xr[2] - xpad,
                  "bottomright" = xr[2] - xpad)
  ann_y <- switch(corr_loc,
                  "topleft"     = yr[2] - ypad,
                  "topright"    = yr[2] - ypad,
                  "bottomleft"  = yr[1] + ypad,
                  "bottomright" = yr[1] + ypad)
  ann_hjust <- if (grepl("right", corr_loc)) 1 else 0
  ann_vjust <- if (grepl("bottom", corr_loc)) 0 else 1

  #-----------------------------
  # 8) Draw plot
  #-----------------------------
  p <- ggplot(dt, aes(x = beta_HEAP, y = beta_GLP1)) +
    geom_hline(yintercept = 0, linewidth = 0.25, color = "grey75") +
    geom_vline(xintercept = 0, linewidth = 0.25, color = "grey75") +

    geom_point(aes(color = mr_edge_sig, size = size_r, alpha = alpha_r)) +

    ggrepel::geom_text_repel(
      data = dt[label_once == TRUE],
      aes(label = Protein),
      size = 3,
      box.padding = 0.2,
      point.padding = 0.12,
      min.segment.length = 0,
      max.overlaps = 60
    ) +

    annotate("label",
             x = ann_x, y = ann_y,
             label = corr_label,
             hjust = ann_hjust, vjust = ann_vjust,
             size = 3.2,
             label.size = 0.25,
             fill = "white") +

    #scale_size_continuous(
    #  range = c(1.6, 5.2),
    #  breaks = c(0.25, 0.5, 0.75),
    #  name = "SomaScan–Olink\n|r|"
    # ) +

    scale_alpha_identity(guide = "none") +
    scale_size_continuous(
      range = c(1.6, 5.2),
      breaks = c(0.25, 0.5, 0.75),
      limits = c(0, 1),
      name = "SomaScan–Olink\nr (pos only)"
    ) +

    scale_color_manual(
      values = c(
        "None"    = "grey75",
        "PDcis"   = "#1b9e77",
        "PDtrans" = "#7570b3",
        "DP"      = "#d95f02"
      ),
      name = paste0("MR significant edge\n(adj.p<", mr_alpha, ")")
    ) +

    labs(
      title = if (show_title) title else NULL,
      subtitle = if (show_subtitle) subtitle else NULL,
      x = if (show_axis_titles) xlab else NULL,
      y = if (show_axis_titles) ylab else NULL
    ) +

    theme_classic(base_size = base_size) +
    theme(
      legend.position = if (legend_position == "none") "none" else legend_position,
      legend.title = element_text(face = "bold"),
      legend.key.height = unit(0.8, "lines"),
      legend.spacing.y = unit(0.2, "lines"),

      # compact spacing
      plot.title.position = "plot",
      plot.margin = margin(6, 6, 6, 6),

      # reduce excess whitespace from long subtitles
      plot.subtitle = element_text(size = base_size - 1)
    ) +
    coord_cartesian(clip = "off")

  return(p)
}


####### TRASH AGAIN ######
# -----------------------------
# Assumes you already created:
#   MR   (from MRres)  with columns Exposure, Protein, Disease, beta_PDcis, padj_PDcis, ...
#   heap (from HEAPint@sList$Model1) with Exposure, Protein, beta_HEAP, GLP1_effect1/2, GLP1_se1/2
#   prot_rel_dt with Protein, r_crossv2 (optional)
# -----------------------------
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
})

# ------------------------------------------------------------
# 0) Ensure heap / MR / prot_rel_dt exist in your environment
# ------------------------------------------------------------
# heap  <- as.data.table(HEAPint@sList$Model1) with cols Exposure, Protein, beta_HEAP, GLP1_effect1/2 ...
# MR    <- as.data.table(MRres) with cols Exposure, Protein, Disease, padj_PDcis, padj_PDtrans, padj_DP, etc.
# prot_rel_dt <- as.data.table(prot_rel) renamed with Protein + r_crossv2

# If you haven't already renamed/created these, do it once:
if (!"beta_HEAP" %in% names(heap)) {
  heap <- as.data.table(HEAPint@sList$Model1)
  setnames(heap,
           old = c("ID","EntrezGeneSymbol","Estimate","Std. Error"),
           new = c("Exposure","Protein","beta_HEAP","se_HEAP"),
           skip_absent = TRUE)
}

if (!"Protein" %in% names(prot_rel_dt)) {
  prot_rel_dt <- as.data.table(prot_rel)
  setnames(prot_rel_dt, "EntrezGeneSymbol", "Protein", skip_absent = TRUE)
}

MR <- as.data.table(MRres)

# ------------------------------------------------------------
# 1) Helper: build and plot one exposure + one GLP arm + one disease
# ------------------------------------------------------------
plot_onepanel_glp1_mr <- function(EXPOSURE_TO_PLOT,
                                  GLP_ARM = c("GLP1_1","GLP1_2"),
                                  disease_for_arm = NULL,
                                  mr_alpha = 0.05,
                                  corr_method = c("pearson","spearman"),
                                  label_n = 12,
                                  add_lm = TRUE) {
  
  GLP_ARM <- match.arg(GLP_ARM)
  corr_method <- match.arg(corr_method)
  
  # Disease mapping by your rule
  if (is.null(disease_for_arm)) {
    disease_for_arm <- if (GLP_ARM == "GLP1_1") "finngen_R12_E4_OBESITY" else "finngen_R12_T2D"
  }
  
  # Pick GLP column names
  glp_beta_col <- if (GLP_ARM == "GLP1_1") "GLP1_effect1" else "GLP1_effect2"
  glp_se_col   <- if (GLP_ARM == "GLP1_1") "GLP1_se1"     else "GLP1_se2"
  
  # ---- HEAP subset (already curated for significance upstream, per you) ----
  heap_sub <- heap[Exposure == EXPOSURE_TO_PLOT, .(
    Exposure, Protein,
    beta_HEAP,
    beta_GLP1 = get(glp_beta_col),
    se_GLP1   = get(glp_se_col)
  )][!is.na(beta_GLP1)]
  
  if (nrow(heap_sub) == 0) {
    stop("No HEAP rows for this exposure with ", GLP_ARM, " effects available.")
  }
  
  # ---- MR subset: one disease only ----
  mr_sub <- MR[Exposure == EXPOSURE_TO_PLOT & Disease == disease_for_arm]
  
  # Merge (proteins without MR will have NAs)
  dt <- merge(heap_sub, mr_sub, by = c("Exposure","Protein"), all.x = TRUE)
  
  # Merge cross-platform correlation
  dt <- merge(dt, prot_rel_dt[, .(Protein, r_crossv2)], by = "Protein", all.x = TRUE)
  
  # Size: |r| in [0,1]
  dt[, olink_soma_absr := pmax(0, pmin(1, abs(r_crossv2)))]
  
  # ---- Determine which MR edge is significant (color) ----
  dt[, sig_PDcis   := !is.na(padj_PDcis)   & padj_PDcis   < mr_alpha]
  dt[, sig_PDtrans := !is.na(padj_PDtrans) & padj_PDtrans < mr_alpha]
  dt[, sig_DP      := !is.na(padj_DP)      & padj_DP      < mr_alpha]
  
  # Priority rule if multiple edges are significant:
  # ---- Determine which MR edge is significant (color) ----
  dt[, sig_PDcis   := !is.na(padj_PDcis)   & padj_PDcis   < mr_alpha]
  dt[, sig_PDtrans := !is.na(padj_PDtrans) & padj_PDtrans < mr_alpha]
  dt[, sig_DP      := !is.na(padj_DP)      & padj_DP      < mr_alpha]
  
  # ---- Causal-priority rule (NOT min-p):
  # PDcis > PDtrans > DP > None
  dt[, mr_edge_sig := fifelse(
    sig_PDcis, "PDcis",
    fifelse(sig_PDtrans, "PDtrans",
            fifelse(sig_DP, "DP", "None"))
  )]
  
  dt[, mr_edge_sig := factor(mr_edge_sig, levels = c("None","PDcis","PDtrans","DP"))]
  
  # ---- Labels: MR-significant + extremes ----
  # label MR sig plus top extremes in HEAP/GLP
  lab_prots <- unique(c(
    dt[mr_edge_sig != "None", Protein][1:min(label_n, sum(dt$mr_edge_sig != "None"))],
    dt[order(-abs(beta_HEAP))][1:ceiling(label_n/2), Protein],
    dt[order(-abs(beta_GLP1))][1:ceiling(label_n/2), Protein]
  ))
  dt[, label_me := Protein %in% lab_prots]
  
  # ---- Correlation between HEAP and GLP1 ----
  # (uses proteins present in this plot)
  ok <- dt[!is.na(beta_HEAP) & !is.na(beta_GLP1)]
  cor_txt <- "r=NA, p=NA"
  if (nrow(ok) >= 3) {
    ct <- suppressWarnings(cor.test(ok$beta_HEAP, ok$beta_GLP1, method = corr_method))
    r  <- unname(ct$estimate)
    p  <- ct$p.value
    cor_txt <- paste0(corr_method, " r=", formatC(r, digits = 2, format = "f"),
                      ", p=", format.pval(p, digits = 2, eps = 1e-300))
  }
  
  # ---- Plot ----
  p <- ggplot(dt, aes(x = beta_HEAP, y = beta_GLP1)) +
    geom_hline(yintercept = 0, linewidth = 0.25, color = "grey70") +
    geom_vline(xintercept = 0, linewidth = 0.25, color = "grey70") +
    
    # optional trend line (subtle)
    {if (add_lm) geom_smooth(method = "lm", se = FALSE, linewidth = 0.4, color = "grey40") else NULL} +
    
    # points: color by MR edge type, size by |r|
    geom_point(aes(color = mr_edge_sig, size = olink_soma_absr),
               alpha = 0.9) +
    
    ggrepel::geom_text_repel(
      data = dt[label_me == TRUE],
      aes(label = Protein),
      size = 3,
      box.padding = 0.25,
      point.padding = 0.15,
      min.segment.length = 0,
      max.overlaps = 60
    ) +
    
    scale_size_continuous(
      range = c(1.7, 5.0),
      breaks = c(0, 0.25, 0.5, 0.75, 1),
      name = "SomaScan–Olink\n|r|"
    ) +
    
    # MR edge colors (explicitly set for stability)
    scale_color_manual(
      values = c(
        "None"    = "grey75",
        "PDcis"   = "#1b9e77",
        "PDtrans" = "#7570b3",
        "DP"      = "#d95f02"
      ),
      name = paste0("MR significant edge\n(adj.p<", mr_alpha, ")")
    ) +
    
    labs(
      title = "HEAP exposure–protein vs GLP1 protein shift (single-panel MR encoding)",
      subtitle = paste0(
        "Exposure: ", EXPOSURE_TO_PLOT,
        " | GLP arm: ", GLP_ARM,
        " | Disease: ", disease_for_arm,
        " | HEAP↔GLP1: ", cor_txt
      ),
      x = "HEAP association beta (Exposure → Protein)",
      y = paste0("GLP1 effect (", GLP_ARM, "; Protein shift)")
    ) +
    
    theme_classic(base_size = 12) +
    theme(
      legend.title = element_text(face = "bold"),
      legend.key.height = unit(0.8, "lines"),
      legend.spacing.y = unit(0.25, "lines"),
      plot.margin = margin(7, 25, 7, 7)
    ) +
    coord_cartesian(clip = "off")
  
  return(p)
}

# ------------------------------------------------------------
# 2) Example calls
# ------------------------------------------------------------
EXPOSURE_TO_PLOT <- "types_of_physical_activity_in_last_4_weeks_f6164_0_0.multi_Strenuous_sports"

p_obesity <- plot_onepanel_glp1_mr(
  EXPOSURE_TO_PLOT = EXPOSURE_TO_PLOT,
  GLP_ARM = "GLP1_1",
  disease_for_arm = "finngen_R12_E4_OBESITY",
  corr_method = "pearson",
  label_n = 14
)

p_t2d <- plot_onepanel_glp1_mr(
  EXPOSURE_TO_PLOT = EXPOSURE_TO_PLOT,
  GLP_ARM = "GLP1_2",
  disease_for_arm = "finngen_R12_T2D",
  corr_method = "pearson",
  label_n = 14
)

print(p_obesity)
print(p_t2d)



head(HEAPint@cList$Model6)
head(HEAPint@pList$Model6)











####### TRASH ######
# Helper to make one plot for one exposure + one GLP arm
plot_exposure_glp1_mr <- function(EXPOSURE_TO_PLOT,
                                  GLP_ARM = c("GLP1_1","GLP1_2"),
                                  disease_for_arm = NULL,
                                  label_n = 12) {
  
  GLP_ARM <- match.arg(GLP_ARM)
  
  # map disease default by arm
  if (is.null(disease_for_arm)) {
    disease_for_arm <- if (GLP_ARM == "GLP1_1") "finngen_R12_E4_OBESITY" else "finngen_R12_T2D"
  }
  
  # pick correct GLP1 effect column names from HEAP table
  glp_beta_col <- if (GLP_ARM == "GLP1_1") "GLP1_effect1" else "GLP1_effect2"
  glp_se_col   <- if (GLP_ARM == "GLP1_1") "GLP1_se1"     else "GLP1_se2"
  
  # ---- 1) HEAP exposure-protein rows for this exposure (one row per protein) ----
  heap_sub <- heap[Exposure == EXPOSURE_TO_PLOT, .(
    Exposure, Protein,
    beta_HEAP,
    se_HEAP,
    beta_GLP1 = get(glp_beta_col),
    se_GLP1   = get(glp_se_col)
  )]
  
  heap_sub <- heap_sub[!is.na(beta_GLP1)]  # must have the relevant GLP arm effect
  
  if (nrow(heap_sub) == 0) {
    message("No HEAP rows found for this exposure with ", GLP_ARM, " effects.")
    return(NULL)
  }
  
  # ---- 2) MR triplets restricted to one disease for interpretability ----
  mr_sub <- MR[Exposure == EXPOSURE_TO_PLOT & Disease == disease_for_arm]
  
  # Merge: now each (Protein) should appear at most once per disease (if not, we’ll collapse)
  dt <- merge(heap_sub, mr_sub, by = c("Exposure","Protein"), all.x = TRUE)
  
  # Optional reliability alpha
  if (exists("prot_rel_dt")) {
    dt <- merge(dt, prot_rel_dt, by = "Protein", all.x = TRUE)
    dt[, rel_w := pmax(0, pmin(1, abs(r_crossv2)))]
  } else {
    dt[, rel_w := 1]
  }
  
  # If there are duplicate Protein rows (rare but can happen), keep the "best" MR evidence row
  # (min padj across PDcis/PDtrans/DP, ignoring NAs)
  dt[, best_p := pmin(padj_PDcis, padj_PDtrans, padj_DP, na.rm = TRUE)]
  dt[is.infinite(best_p), best_p := NA_real_]
  setorder(dt, Protein, best_p)
  dt <- dt[, .SD[1], by = .(Exposure, Protein)]
  
  # ---- 3) Long format for separate MR edges ----
  long <- rbindlist(list(
    dt[, .(Exposure, Protein, Disease = disease_for_arm, beta_HEAP, se_HEAP, beta_GLP1, se_GLP1, rel_w,
           edge = "PDcis", beta = beta_PDcis, se = se_PDcis, padj = padj_PDcis)],
    dt[, .(Exposure, Protein, Disease = disease_for_arm, beta_HEAP, se_HEAP, beta_GLP1, se_GLP1, rel_w,
           edge = "PDtrans", beta = beta_PDtrans, se = se_PDtrans, padj = padj_PDtrans)],
    dt[, .(Exposure, Protein, Disease = disease_for_arm, beta_HEAP, se_HEAP, beta_GLP1, se_GLP1, rel_w,
           edge = "DP", beta = beta_DP, se = se_DP, padj = padj_DP)]
  ), use.names = TRUE, fill = TRUE)
  
  long[, sig := !is.na(padj) & padj < 0.05]
  long[, edge := factor(edge, levels = c("PDcis","PDtrans","DP"))]
  
  # ---- 4) Labels: choose a few proteins (avoid repeated labels across facets) ----
  # Label proteins with strongest HEAP/GLP1 or any significant MR in any panel
  lab_prots <- unique(c(
    long[order(-abs(beta_HEAP))][1:label_n, Protein],
    long[order(-abs(beta_GLP1))][1:label_n, Protein],
    long[sig == TRUE][1:label_n, Protein]
  ))
  long[, label_me := Protein %in% lab_prots]
  
  # ---- 5) Plot ----
  p <- ggplot(long, aes(x = beta_HEAP, y = beta_GLP1)) +
    geom_hline(yintercept = 0, linewidth = 0.3) +
    geom_vline(xintercept = 0, linewidth = 0.3) +
    geom_point(aes(color = sig, alpha = rel_w), size = 2.6, stroke = 0.25) +
    ggrepel::geom_text_repel(
      data = long[label_me == TRUE],
      aes(label = Protein),
      size = 3, max.overlaps = 50,
      box.padding = 0.25, point.padding = 0.15,
      min.segment.length = 0
    ) +
    facet_wrap(~edge, nrow = 1) +
    scale_alpha_continuous(range = c(0.25, 1), guide = "none") +
    scale_color_manual(values = c(`FALSE` = "grey60", `TRUE` = "#1b9e77"),
                       name = "MR edge\n(adj.p<0.05)") +
    labs(
      title = "HEAP exposure–protein vs GLP1 protein shift (MR shown separately)",
      subtitle = paste0(
        "Exposure: ", EXPOSURE_TO_PLOT,
        " | GLP arm: ", GLP_ARM,
        " | Disease for MR: ", disease_for_arm
      ),
      x = "HEAP association beta (Exposure → Protein)",
      y = paste0("GLP1 effect (", GLP_ARM, "; Protein shift)")
    ) +
    theme_bw(base_size = 12) +
    theme(panel.grid = element_blank(),
          strip.background = element_rect(fill = "grey95", color = NA))
  
  return(p)
}

# -----------------------------
# Example usage
# -----------------------------
EXPOSURE_TO_PLOT <- "types_of_physical_activity_in_last_4_weeks_f6164_0_0.multi_Strenuous_sports"

p_obesity <- plot_exposure_glp1_mr(
  EXPOSURE_TO_PLOT = EXPOSURE_TO_PLOT,
  GLP_ARM = "GLP1_1",
  disease_for_arm = "finngen_R12_E4_OBESITY"
)

p_t2d <- plot_exposure_glp1_mr(
  EXPOSURE_TO_PLOT = EXPOSURE_TO_PLOT,
  GLP_ARM = "GLP1_2",
  disease_for_arm = "finngen_R12_T2D"
)

print(p_obesity)
print(p_t2d)


exposurelist <- sort(unique(df_glp1$Exposure))















# ============================================================
# PLOT 2: GLP1 (y) vs MR P→D (x) for one disease (triangulation plot)
# - Color: whether HEAP beta sign matches GLP1 sign (optional)
# - Shape: PD_best_type (cis/trans)
# ============================================================
DISEASE_TO_PLOT <- "finngen_R12_T2D"  # <-- change to your T2D code
p2_dt <- df[ Disease == DISEASE_TO_PLOT &
               !is.na(beta_GLP1_1) & !is.na(beta_PD_best) ]

if (nrow(p2_dt) == 0) {
  message("No rows for that disease with GLP1 + MR PD; pick a Disease that exists in MRres and proteins with GLP1 effects.")
} else {
  # “protective” definition depends on coding; keep as directional concordance plot
  p2_dt[, conc_GLP1_MRPD := sign(beta_GLP1_1) == sign(beta_PD_best)]
  p2_dt[, label_me := rank(-abs(beta_PD_best), ties.method="first") <= 10 |
          rank(-abs(beta_GLP1_1), ties.method="first") <= 10]
  
  p2 <- ggplot(p2_dt, aes(x = beta_PD_best, y = beta_GLP1_1)) +
    geom_hline(yintercept = 0, linewidth = 0.3) +
    geom_vline(xintercept = 0, linewidth = 0.3) +
    geom_point(aes(color = conc_GLP1_MRPD, shape = PD_best_type, alpha = rel_w),
               size = 2.6, stroke = 0.3) +
    ggrepel::geom_text_repel(
      data = p2_dt[label_me == TRUE],
      aes(label = Protein),
      max.overlaps = 50, size = 3, box.padding = 0.3, point.padding = 0.2
    ) +
    scale_alpha_continuous(range = c(0.25, 1), guide = "none") +
    scale_color_manual(values = c(`FALSE` = "grey55", `TRUE` = "#7570b3"),
                       name = "Concordant\n(GLP1 vs MR P→D)") +
    labs(
      title = "GLP1 protein shift vs MR protein→disease effect",
      subtitle = paste("Disease:", DISEASE_TO_PLOT),
      x = "MR beta (Protein → Disease) (best cis/trans)",
      y = "GLP1 effect (STEP; Protein shift)"
    ) +
    theme_bw(base_size = 12) +
    theme(panel.grid = element_blank())
  
  print(p2)
}
