#!/usr/bin/env Rscript

# ============================================================
# HEAP vs GLP1/HERITAGE scatter with MR overlay
# Shapes encode MR replication category:
#   - UKB only
#   - DECODE only
#   - Both (replicated)
#
# MR edge colors encode which MR edge is significant (PDcis/PDtrans/DP/None)
# Point size encodes SomaScan–Olink reliability (positive r only; 0 otherwise)
#
# NOTE:
# - This script assumes UKB + DECODE MRmotifs.csv have compatible columns.
# - Replication ("Both") is defined at the triplet level (Exposure+Protein+Disease),
#   optionally restricted to any_sig==TRUE within each dataset.
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(ggplot2)
  library(ggrepel)
  library(qs)
})

# -----------------------------
# File paths
# -----------------------------
MR_UKB_FP <- "/n/groups/patel/shakson_ukb/UK_Biobank/Output/MR_edges/summary/MRmotifs.csv"
MR_DEC_FP <- "/n/groups/patel/shakson_ukb/UK_Biobank/Output/MR_edges/summary/DECODE/MRmotifs.csv"

HEAPINT_FP <- "/n/groups/patel/shakson_ukb/UK_Biobank/Output/HEAPres/HEAPintv2.qs"
PROTREL_FP <- "/n/groups/patel/IGLOO/UKB/OlinkSoma/OlinkSoma.csv"

OUTDIR <- "/n/groups/patel/shakson_ukb/UK_Biobank/Output/MRplots/"
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

# Toggle: should "replicated" mean "present in both" OR "significant in both"?
# Most people mean "significant in both", so default TRUE.
KEEP_ONLY_ANY_SIG_PER_DATASET <- TRUE

# -----------------------------
# Utilities
# -----------------------------
fmt_p <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 1e-300) return("<1e-300")
  format.pval(p, digits = 2, eps = 1e-300)
}

canon_exposure <- function(x) sub("_[0-9]+$", "", x)

# pmin that ignores NAs
pmin_na <- function(...) {
  xs <- list(...)
  if (length(xs) == 1) return(xs[[1]])
  m <- do.call(cbind, xs)
  apply(m, 1, function(v) {
    v <- v[!is.na(v)]
    if (length(v) == 0) NA_real_ else min(v)
  })
}

# pick value (beta/se) from the dataset with smaller p-value; if only one exists use it
pick_by_best_p <- function(p1, v1, p2, v2) {
  out <- rep(NA_real_, length(p1))
  has1 <- !is.na(p1) & !is.na(v1)
  has2 <- !is.na(p2) & !is.na(v2)
  
  out[has1 & !has2] <- v1[has1 & !has2]
  out[has2 & !has1] <- v2[has2 & !has1]
  
  both <- has1 & has2
  out[both] <- ifelse(p1[both] <= p2[both], v1[both], v2[both])
  out
}

# Robust column getter
get_or_na <- function(dt, col) {
  if (col %in% names(dt)) dt[[col]] else rep(NA_real_, nrow(dt))
}

# -----------------------------
# Load HEAPint
# -----------------------------
HEAPint <- qread(HEAPINT_FP)

# -----------------------------
# Load SomaScan vs Olink reliability
# -----------------------------
prot_rel <- fread(PROTREL_FP, skip = 3) %>%
  select(c("gene_name","olink_nonnorm_corr","olink_smpnorm_corr")) %>%
  setNames(c("EntrezGeneSymbol", "r_cross", "r_crossv2")) %>%
  na.omit()

prot_rel_dt <- as.data.table(prot_rel)
setnames(prot_rel_dt, "EntrezGeneSymbol", "Protein", skip_absent = TRUE)
if (!"r_crossv2" %in% names(prot_rel_dt) && "r_cross" %in% names(prot_rel_dt)) {
  prot_rel_dt[, r_crossv2 := r_cross]
}

# -----------------------------
# Load MR (UKB + DECODE) and build replication category per triplet
# -----------------------------
read_mr <- function(fp, tag) {
  x <- fread(fp)
  x[, mr_dataset := tag]
  # stable triplet key (use triplet col if exists + non-empty)
  if ("triplet" %in% names(x)) {
    x[, trip_key := ifelse(!is.na(triplet) & triplet != "",
                           triplet,
                           paste(Exposure, Protein, Disease, sep = "||"))]
  } else {
    x[, trip_key := paste(Exposure, Protein, Disease, sep = "||")]
  }
  x
}

MR_ukb <- read_mr(MR_UKB_FP, "UKB")
MR_dec <- read_mr(MR_DEC_FP, "DECODE")

if (KEEP_ONLY_ANY_SIG_PER_DATASET) {
  if ("any_sig" %in% names(MR_ukb)) MR_ukb <- MR_ukb[any_sig == TRUE]
  if ("any_sig" %in% names(MR_dec)) MR_dec <- MR_dec[any_sig == TRUE]
}

# Determine replication status by trip_key membership
keys_ukb <- unique(MR_ukb$trip_key)
keys_dec <- unique(MR_dec$trip_key)

all_keys <- unique(c(keys_ukb, keys_dec))
rep_dt <- data.table(
  trip_key = all_keys,
  in_ukb = all_keys %in% keys_ukb,
  in_dec = all_keys %in% keys_dec
)
rep_dt[, mr_support := fifelse(in_ukb & in_dec, "Both",
                               fifelse(in_ukb, "UKB only",
                                       fifelse(in_dec, "DECODE only", "None")))]
rep_dt[, mr_support := factor(mr_support, levels = c("UKB only","DECODE only","Both","None"))]

# Merge UKB+DECODE MR rows into one row per triplet (trip_key),
# keeping *combined* evidence for edge p-values and a "best" beta/se
# (best = dataset with smaller p for that edge).
combine_two <- function(u, d) {
  # Full outer join on trip_key; suffix columns
  uu <- copy(u)
  dd <- copy(d)
  setkey(uu, trip_key); setkey(dd, trip_key)
  m <- merge(uu, dd, by = "trip_key", all = TRUE, suffixes = c("_ukb","_dec"))
  
  # Carry key columns
  # Prefer UKB values when present; else DECODE
  m[, Exposure := fifelse(!is.na(Exposure_ukb), Exposure_ukb, Exposure_dec)]
  m[, Protein  := fifelse(!is.na(Protein_ukb),  Protein_ukb,  Protein_dec)]
  m[, Disease  := fifelse(!is.na(Disease_ukb),  Disease_ukb,  Disease_dec)]
  
  # Attach replication category
  m <- merge(m, rep_dt[, .(trip_key, mr_support)], by = "trip_key", all.x = TRUE)
  
  # Combine edge p-values: take min non-NA across datasets
  # (these columns exist in your UKB MRmotifs.csv; if DECODE lacks some, it still works)
  m[, padj_PDcis   := pmin_na(get_or_na(m, "padj_PDcis_ukb"),   get_or_na(m, "padj_PDcis_dec"))]
  m[, padj_PDtrans := pmin_na(get_or_na(m, "padj_PDtrans_ukb"), get_or_na(m, "padj_PDtrans_dec"))]
  m[, padj_DP      := pmin_na(get_or_na(m, "padj_DP_ukb"),      get_or_na(m, "padj_DP_dec"))]
  
  # Combine betas/se by picking dataset with better (smaller) p-value for that edge
  m[, beta_PDcis := pick_by_best_p(get_or_na(m, "padj_PDcis_ukb"),   get_or_na(m, "beta_PDcis_ukb"),
                                   get_or_na(m, "padj_PDcis_dec"),   get_or_na(m, "beta_PDcis_dec"))]
  m[, se_PDcis   := pick_by_best_p(get_or_na(m, "padj_PDcis_ukb"),   get_or_na(m, "se_PDcis_ukb"),
                                   get_or_na(m, "padj_PDcis_dec"),   get_or_na(m, "se_PDcis_dec"))]
  
  m[, beta_PDtrans := pick_by_best_p(get_or_na(m, "padj_PDtrans_ukb"), get_or_na(m, "beta_PDtrans_ukb"),
                                     get_or_na(m, "padj_PDtrans_dec"), get_or_na(m, "beta_PDtrans_dec"))]
  m[, se_PDtrans   := pick_by_best_p(get_or_na(m, "padj_PDtrans_ukb"), get_or_na(m, "se_PDtrans_ukb"),
                                     get_or_na(m, "padj_PDtrans_dec"), get_or_na(m, "se_PDtrans_dec"))]
  
  m[, beta_DP := pick_by_best_p(get_or_na(m, "padj_DP_ukb"), get_or_na(m, "beta_DP_ukb"),
                                get_or_na(m, "padj_DP_dec"), get_or_na(m, "beta_DP_dec"))]
  m[, se_DP   := pick_by_best_p(get_or_na(m, "padj_DP_ukb"), get_or_na(m, "se_DP_ukb"),
                                get_or_na(m, "padj_DP_dec"), get_or_na(m, "se_DP_dec"))]
  
  # If you also want EP/ED columns, add them similarly.
  # Here we keep beta_EP/padj_EP if present, preferring min pval (optional).
  if ("padj_EP_ukb" %in% names(m) || "padj_EP_dec" %in% names(m)) {
    m[, padj_EP := pmin_na(get_or_na(m, "padj_EP_ukb"), get_or_na(m, "padj_EP_dec"))]
    m[, beta_EP := pick_by_best_p(get_or_na(m, "padj_EP_ukb"), get_or_na(m, "beta_EP_ukb"),
                                  get_or_na(m, "padj_EP_dec"), get_or_na(m, "beta_EP_dec"))]
  }
  
  # keep a lightweight MR table
  keep <- c("trip_key","Exposure","Protein","Disease","mr_support",
            "padj_PDcis","beta_PDcis","se_PDcis",
            "padj_PDtrans","beta_PDtrans","se_PDtrans",
            "padj_DP","beta_DP","se_DP",
            "padj_EP","beta_EP")
  keep <- keep[keep %in% names(m)]
  m[, ..keep]
}

MR_combined <- combine_two(MR_ukb, MR_dec)

# -----------------------------
# Plot function: HEAP vs GLP1/HERIT with MR overlay
# -----------------------------
plot_HEAP_GLP1_MR_onepanel <- function(
    EXPOSURE_TO_PLOT,
    arm = c("GLP1_1", "GLP1_2", "HERITAGE"),
    disease_for_arm = NULL,
    mr_alpha = 0.05,
    
    corr_model = "Model6",
    corr_loc = c("topleft","topright","bottomleft","bottomright"),
    corr_digits = 2,
    
    label_n = 20,
    label_only_mr = FALSE,
    
    legend_position = c("right","bottom","none"),
    base_size = 12,
    
    title = NULL,
    subtitle = NULL,
    xlab = "UKB assoc (E→P) beta",
    ylab = NULL
) {
  arm <- match.arg(arm)
  corr_loc <- match.arg(corr_loc)
  legend_position <- match.arg(legend_position)
  
  y_col <- switch(
    arm,
    "GLP1_1"    = "GLP1_effect1",
    "GLP1_2"    = "GLP1_effect2",
    "HERITAGE"  = "HERITAGE_effect"
  )
  if (is.null(ylab)) ylab <- y_col
  if (is.null(title)) title <- paste0("HEAP vs ", arm, " protein shifts (MR support encoded)")
  
  exposure_key_plot <- canon_exposure(EXPOSURE_TO_PLOT)
  
  # ---- HEAP table subset ----
  heap_dt <- as.data.table(HEAPint@sList[[corr_model]])
  setnames(heap_dt,
           old = c("ID","EntrezGeneSymbol","Estimate","Std. Error"),
           new = c("Exposure","Protein","beta_HEAP","se_HEAP"),
           skip_absent = TRUE)
  
  stopifnot(all(c("Exposure","Protein","beta_HEAP") %in% names(heap_dt)))
  if (!y_col %in% names(heap_dt)) stop("Column ", y_col, " not found in HEAPint@sList[[", corr_model, "]].")
  
  heap_dt[, Exposure_key := canon_exposure(Exposure)]
  
  heap_pick <- if (any(heap_dt$Exposure == EXPOSURE_TO_PLOT, na.rm = TRUE)) {
    heap_dt[Exposure == EXPOSURE_TO_PLOT]
  } else {
    heap_dt[Exposure_key == exposure_key_plot]
  }
  
  heap_sub <- heap_pick[
    , .(
      Exposure  = EXPOSURE_TO_PLOT,
      beta_HEAP = mean(beta_HEAP, na.rm = TRUE),
      beta_arm  = mean(get(y_col), na.rm = TRUE)
    ),
    by = .(Protein)
  ][is.finite(beta_HEAP) & is.finite(beta_arm)]
  
  if (nrow(heap_sub) == 0) stop("No HEAP rows found for exposure with available ", arm, " effects.")
  
  # ---- MR subset (by exposure key, optional disease) ----
  MR <- copy(MR_combined)
  MR[, Exposure_key := canon_exposure(Exposure)]
  mr_sub <- MR[Exposure_key == exposure_key_plot]
  if (!is.null(disease_for_arm)) mr_sub <- mr_sub[Disease == disease_for_arm]
  
  # Merge: keep plotted proteins, attach MR info if present
  dt <- merge(heap_sub, mr_sub, by = c("Protein"), all.x = TRUE, allow.cartesian = TRUE)
  
  # ---- Reliability: size encodes positive r only ----
  dt <- merge(dt, prot_rel_dt[, .(Protein, r_crossv2)], by = "Protein", all.x = TRUE)
  dt[, olink_soma_r := r_crossv2]
  dt[, size_r := fifelse(is.na(olink_soma_r) | olink_soma_r <= 0, 0.05,
                         fifelse(olink_soma_r < 0.1, 0.1,
                                 pmin(olink_soma_r, 1)))]
  #dt[, size_r := fifelse(!is.na(olink_soma_r) & olink_soma_r > 0,
  #                       pmin(olink_soma_r, 1),
  #                       0.05)]
  #dt[, size_r := pmax(0, olink_soma_r)]
  #dt[is.na(size_r), size_r := 0]
  #dt[, alpha_r := fifelse(!is.na(olink_soma_r) & olink_soma_r > 0, 0.9, 0.08)]
  
  # ---- MR edge significance category (priority: PDcis > PDtrans > DP) ----
  # Safe when missing: treat as not significant
  dt[, sig_PDcis   := !is.na(padj_PDcis)   & padj_PDcis   < mr_alpha]
  dt[, sig_PDtrans := !is.na(padj_PDtrans) & padj_PDtrans < mr_alpha]
  dt[, sig_DP      := !is.na(padj_DP)      & padj_DP      < mr_alpha]
  
  dt[, mr_edge_sig := fifelse(sig_PDcis, "PDcis",
                              fifelse(sig_PDtrans, "PDtrans",
                                      fifelse(sig_DP, "DP", "None")))]
  dt[, mr_edge_sig := factor(mr_edge_sig, levels = c("None","PDcis","PDtrans","DP"))]
  
  # ---- MR replication shape (UKB only / DECODE only / Both / None) ----
  if (!"mr_support" %in% names(dt)) dt[, mr_support := "None"]
  dt[, mr_support := fct_explicit_na(as.factor(mr_support), na_level = "None")]
  dt[, mr_support := factor(mr_support, levels = c("None","UKB only","DECODE only","Both"))]
  
  # ---- Labels (priority: edge PDcis >> PDtrans >> DP >> None;
  #              within: Both >> single >> None) ----
  cap_n <- function(x, n) x[seq_len(min(length(x), n))]
  
  # Define support tier: Both (2) > single (1) > None (0)
  dt[, support_tier := fifelse(mr_support == "Both", 2L,
                               fifelse(mr_support %in% c("UKB only","DECODE only"), 1L, 0L))]
  
  # Define edge tier: PDcis (3) > PDtrans (2) > DP (1) > None (0)
  dt[, edge_tier := fifelse(mr_edge_sig == "PDcis", 3L,
                            fifelse(mr_edge_sig == "PDtrans", 2L,
                                    fifelse(mr_edge_sig == "DP", 1L, 0L)))]
  
  # Ranking score: first by edge_tier, then support_tier, then effect magnitude
  # (use abs(beta_arm) as primary; tie-break by abs(beta_HEAP))
  dt[, label_rank := 1e6*edge_tier + 1e3*support_tier + 10*abs(beta_arm) + abs(beta_HEAP)]
  
  # Candidate pools:
  # - If label_only_mr=TRUE: only label points that have a non-None edge
  # - Else: fill remaining labels with biggest absolute effects
  if (isTRUE(label_only_mr)) {
    cand <- dt[mr_edge_sig != "None"][order(-label_rank)]
  } else {
    cand <- dt[order(-label_rank)]
  }
  
  lab_prots <- unique(cand$Protein)
  lab_prots <- cap_n(lab_prots, label_n)
  
  dt[, label_once := (Protein %in% lab_prots) & !duplicated(Protein)]
  
  # ---- Correlation annotation (from HEAPint cList/pList) ----
  ctab <- as.data.table(HEAPint@cList[[corr_model]], keep.rownames = "Exposure")
  ptab <- as.data.table(HEAPint@pList[[corr_model]], keep.rownames = "Exposure")
  
  ctab[, Exposure_key := canon_exposure(Exposure)]
  ptab[, Exposure_key := canon_exposure(Exposure)]
  
  if (!y_col %in% names(ctab) || !y_col %in% names(ptab)) {
    r_val <- NA_real_; p_val <- NA_real_
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
      paste0("Exposure: ", EXPOSURE_TO_PLOT, " | MR support shape: UKB vs DECODE replication")
    } else {
      paste0("Exposure: ", EXPOSURE_TO_PLOT, " | Disease: ", disease_for_arm, " | MR support shape: replication")
    }
  }
  
  # ---- Annotation placement ----
  xr <- range(dt$beta_HEAP, na.rm = TRUE)
  yr <- range(dt$beta_arm,  na.rm = TRUE)
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
  
  # ---- Plot ----
  p <- ggplot(dt, aes(x = beta_HEAP, y = beta_arm)) +
    geom_hline(yintercept = 0, linewidth = 0.25, color = "grey75") +
    geom_vline(xintercept = 0, linewidth = 0.25, color = "grey75") +
    
    geom_point(aes(color = mr_edge_sig, shape = mr_support, size = size_r)) + #, alpha = alpha_r)) +
    
    ggrepel::geom_text_repel(
      data = dt[label_once == TRUE],
      aes(label = Protein),
      size = 3,
      box.padding = 0.2,
      point.padding = 0.12,
      min.segment.length = 0,
      max.overlaps = 60
    ) +
    
    annotate(
      "label",
      x = ann_x, y = ann_y,
      label = corr_label,
      hjust = ann_hjust, vjust = ann_vjust,
      size = 3.2,
      label.size = 0.25,
      fill = "white"
    ) +
    
    scale_alpha_identity(guide = "none") +
    
    scale_size_continuous(
      range = c(1.6, 6.0),
      limits = c(0.05, 1),
      breaks = c(0.05, 0.25, 0.5, 0.75, 1),
      labels = c("≤0 / NA", "0.25", "0.5", "0.75", "1.0"),
      name = "SomaScan–Olink\nr (pos only)"
    ) +
    
    #scale_size_continuous(
    #  range = c(1.6, 5.2),
    #  breaks = c(0.25, 0.5, 0.75),
    #  limits = c(0, 1),
    #  name = "SomaScan–Olink\nr (pos only)"
    #) +
    
    scale_color_manual(
      values = c(
        "None"    = "grey75",
        "PDcis"   = "#1b9e77",
        "PDtrans" = "#7570b3",
        "DP"      = "#d95f02"
      ),
      name = paste0("MR significant edge\n(adj.p<", mr_alpha, ")")
    ) +
    
    scale_shape_manual(
      values = c(
        "None"      = 16,
        "UKB only"  = 17,
        "DECODE only" = 15,
        "Both"      = 18
      ),
      name = "MR support\n(UKB vs DECODE)"
    ) +
    
    labs(
      title = title,
      subtitle = subtitle,
      x = xlab,
      y = ylab
    ) +
    
    theme_classic(base_size = base_size) +
    theme(
      legend.position = if (legend_position == "none") "none" else legend_position,
      legend.title = element_text(face = "bold"),
      legend.key.height = grid::unit(0.8, "lines"),
      legend.spacing.y = grid::unit(0.2, "lines"),
      plot.title.position = "plot",
      plot.margin = margin(6, 6, 6, 6),
      plot.subtitle = element_text(size = base_size - 1)
    ) +
    coord_cartesian(clip = "off")
  
  return(p)
}

# ============================================================
# Example usage (your same examples)
# ============================================================

# 1) Exercise vs HERITAGE (obesity MR disease context)
EXPOSURE_TO_PLOT <- "types_of_physical_activity_in_last_4_weeks_f6164_0_0.multi_Strenuous_sports"

p0 <- plot_HEAP_GLP1_MR_onepanel(
  EXPOSURE_TO_PLOT = EXPOSURE_TO_PLOT,
  arm = "HERITAGE",
  disease_for_arm = "finngen_R12_E4_OBESITY",
  title = "Exercise → proteins vs HERITAGE (Endurance Exercise)",
  subtitle = "Strenuous sports (UKB) | Shape = MR replication (UKB/DECODE/Both)",
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
  subtitle = "Strenuous sports (UKB) | Shape = MR replication (UKB/DECODE/Both)",
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
  subtitle = "Strenuous sports (UKB) | Shape = MR replication (UKB/DECODE/Both)",
  xlab = "UKB assoc (E→P) beta",
  ylab = "GLP1_2 protein shift",
  legend_position = "right",
  corr_loc = "topleft"
)

print(p0); print(p1); print(p2)

ggsave(file.path(OUTDIR, "MRHERITAGE_StrenSports_replShape.png"),
       plot = p0, dpi = 1000, width = 7, height = 4, units = "in")
ggsave(file.path(OUTDIR, "MRGLP1_STEP1_StrenSports_replShape.png"),
       plot = p1, dpi = 1000, width = 7, height = 4, units = "in")
ggsave(file.path(OUTDIR, "MRGLP1_STEP2_StrenSports_replShape.png"),
       plot = p2, dpi = 1000, width = 7, height = 4, units = "in")

# 2) Diet (fruit intake)
p3 <- plot_HEAP_GLP1_MR_onepanel(
  EXPOSURE_TO_PLOT = "fresh_fruit_intake_f1309_0_0",
  arm = "GLP1_1",
  disease_for_arm = "finngen_R12_E4_OBESITY",
  title = "Diet → proteins vs GLP1_1 (obesity cohort)",
  subtitle = "Fruit intake (UKB) | Shape = MR replication (UKB/DECODE/Both)",
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
  subtitle = "Fruit intake (UKB) | Shape = MR replication (UKB/DECODE/Both)",
  xlab = "UKB assoc (E→P) beta",
  ylab = "GLP1_2 protein shift",
  legend_position = "right",
  corr_loc = "topleft"
)

print(p3); print(p4)

ggsave(file.path(OUTDIR, "MRGLP1_STEP1_FruitIntake_replShape.png"),
       plot = p3, dpi = 1000, width = 7, height = 4, units = "in")
ggsave(file.path(OUTDIR, "MRGLP1_STEP2_FruitIntake_replShape.png"),
       plot = p4, dpi = 1000, width = 7, height = 4, units = "in")

# 3) Smoking example (note your relaxed matching already handled by canon_exposure)
p5 <- plot_HEAP_GLP1_MR_onepanel(
  EXPOSURE_TO_PLOT = "past_tobacco_smoking_f1249_0_0",
  arm = "GLP1_1",
  disease_for_arm = "finngen_R12_E4_OBESITY",
  title = "Quit Smoking → proteins vs GLP1_1 (obesity cohort)",
  subtitle = "Former smoker | Shape = MR replication (UKB/DECODE/Both)",
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
  subtitle = "Former smoker | Shape = MR replication (UKB/DECODE/Both)",
  xlab = "UKB assoc (E→P) beta",
  ylab = "GLP1_2 protein shift",
  legend_position = "right",
  corr_loc = "topleft"
)

print(p5); print(p6)

ggsave(file.path(OUTDIR, "MRGLP1_STEP1_QuitSmoking_replShape.png"),
       plot = p5, dpi = 1000, width = 7, height = 4, units = "in")
ggsave(file.path(OUTDIR, "MRGLP1_STEP2_QuitSmoking_replShape.png"),
       plot = p6, dpi = 1000, width = 7, height = 4, units = "in")

cat("\nDone. Wrote plots to:\n  ", OUTDIR, "\n\n", sep = "")