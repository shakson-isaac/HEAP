#!/usr/bin/env Rscript

# ============================================================
# HEAP vs GLP1/HERITAGE scatter with MR overlay
#
# Shapes encode MR replication category (EDGE-TYPE MATCH REQUIRED):
#   - UKB only      (this triplet has a significant edge in UKB, but not DECODE)
#   - DECODE only   (significant edge in DECODE, but not UKB)
#   - Both          (SAME MR EDGE TYPE is significant in BOTH UKB and DECODE)
#   - None          (no significant MR edge in either)
#
# MR edge colors encode which MR edge is significant (PDcis/PDtrans/DP/None)
# Point size encodes SomaScan–Olink reliability:
#   - r <= 0 or NA  -> size_r = 0.05 (minimum visible)
#   - 0 < r < 0.1   -> size_r = 0.1  ("any positive gets credit")
#   - r >= 0.1      -> size_r = min(r, 1) (continuous)
#
# NOTE:
# - "Both" now REQUIRES replication of the SAME edge type (PDcis, PDtrans, or DP).
# - If a triplet is significant in BOTH datasets but for DIFFERENT edge types,
#   it will be treated as SINGLE-SOURCE support (UKB only / DECODE only) based on
#   which dataset provides the smaller p-value for the edge type used in coloring.
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

# Filter each MR file to any_sig==TRUE BEFORE replication logic (recommended)
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
# Load MR (UKB + DECODE) and build EDGE-TYPE replicated support per triplet
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

# Add edge-specific significance flags (assumes these columns exist; if not, they become FALSE)
for (xname in c("MR_ukb","MR_dec")) {
  x <- get(xname)
  x[, sig_PDcis   := (!is.na(get_or_na(x, "padj_PDcis"))   & get_or_na(x, "padj_PDcis")   < 0.05)]
  x[, sig_PDtrans := (!is.na(get_or_na(x, "padj_PDtrans")) & get_or_na(x, "padj_PDtrans") < 0.05)]
  x[, sig_DP      := (!is.na(get_or_na(x, "padj_DP"))      & get_or_na(x, "padj_DP")      < 0.05)]
  assign(xname, x)
}

# Summarize per triplet: does each dataset have each significant edge type?
ukb_keys <- MR_ukb[, .(
  in_ukb = TRUE,
  ukb_PDcis   = any(sig_PDcis),
  ukb_PDtrans = any(sig_PDtrans),
  ukb_DP      = any(sig_DP)
), by = trip_key]

dec_keys <- MR_dec[, .(
  in_dec = TRUE,
  dec_PDcis   = any(sig_PDcis),
  dec_PDtrans = any(sig_PDtrans),
  dec_DP      = any(sig_DP)
), by = trip_key]

rep_dt <- merge(ukb_keys, dec_keys, by = "trip_key", all = TRUE)
for (cc in setdiff(names(rep_dt), "trip_key")) rep_dt[is.na(get(cc)), (cc) := FALSE]

# SAME EDGE TYPE replicated?
rep_dt[, rep_PDcis   := ukb_PDcis   & dec_PDcis]
rep_dt[, rep_PDtrans := ukb_PDtrans & dec_PDtrans]
rep_dt[, rep_DP      := ukb_DP      & dec_DP]
rep_dt[, replicated_same_edge := rep_PDcis | rep_PDtrans | rep_DP]

# We'll compute "mr_support" later AFTER we decide which edge type we are coloring by
# (because if edges differ across datasets, we want to label support based on the chosen edge)

# -----------------------------
# Merge UKB+DECODE MR rows into one row per triplet (trip_key),
# keeping combined p-values (min) and "best" beta/se by smaller p
# -----------------------------
combine_two <- function(u, d) {
  uu <- copy(u)
  dd <- copy(d)
  setkey(uu, trip_key); setkey(dd, trip_key)
  m <- merge(uu, dd, by = "trip_key", all = TRUE, suffixes = c("_ukb","_dec"))
  
  m[, Exposure := fifelse(!is.na(Exposure_ukb), Exposure_ukb, Exposure_dec)]
  m[, Protein  := fifelse(!is.na(Protein_ukb),  Protein_ukb,  Protein_dec)]
  m[, Disease  := fifelse(!is.na(Disease_ukb),  Disease_ukb,  Disease_dec)]
  
  # attach replication summaries (edge-type flags)
  m <- merge(m,
             rep_dt[, .(trip_key, in_ukb, in_dec,
                        rep_PDcis, rep_PDtrans, rep_DP,
                        ukb_PDcis, ukb_PDtrans, ukb_DP,
                        dec_PDcis, dec_PDtrans, dec_DP)],
             by = "trip_key", all.x = TRUE)
  
  # Combine edge p-values: min across datasets (used for edge-significance coloring)
  m[, padj_PDcis   := pmin_na(get_or_na(m, "padj_PDcis_ukb"),   get_or_na(m, "padj_PDcis_dec"))]
  m[, padj_PDtrans := pmin_na(get_or_na(m, "padj_PDtrans_ukb"), get_or_na(m, "padj_PDtrans_dec"))]
  m[, padj_DP      := pmin_na(get_or_na(m, "padj_DP_ukb"),      get_or_na(m, "padj_DP_dec"))]
  
  # Combine betas/se by picking dataset with better (smaller) p-value for that edge
  m[, beta_PDcis := pick_by_best_p(get_or_na(m, "padj_PDcis_ukb"), get_or_na(m, "beta_PDcis_ukb"),
                                   get_or_na(m, "padj_PDcis_dec"), get_or_na(m, "beta_PDcis_dec"))]
  m[, se_PDcis   := pick_by_best_p(get_or_na(m, "padj_PDcis_ukb"), get_or_na(m, "se_PDcis_ukb"),
                                   get_or_na(m, "padj_PDcis_dec"), get_or_na(m, "se_PDcis_dec"))]
  
  m[, beta_PDtrans := pick_by_best_p(get_or_na(m, "padj_PDtrans_ukb"), get_or_na(m, "beta_PDtrans_ukb"),
                                     get_or_na(m, "padj_PDtrans_dec"), get_or_na(m, "beta_PDtrans_dec"))]
  m[, se_PDtrans   := pick_by_best_p(get_or_na(m, "padj_PDtrans_ukb"), get_or_na(m, "se_PDtrans_ukb"),
                                     get_or_na(m, "padj_PDtrans_dec"), get_or_na(m, "se_PDtrans_dec"))]
  
  m[, beta_DP := pick_by_best_p(get_or_na(m, "padj_DP_ukb"), get_or_na(m, "beta_DP_ukb"),
                                get_or_na(m, "padj_DP_dec"), get_or_na(m, "beta_DP_dec"))]
  m[, se_DP   := pick_by_best_p(get_or_na(m, "padj_DP_ukb"), get_or_na(m, "se_DP_ukb"),
                                get_or_na(m, "padj_DP_dec"), get_or_na(m, "se_DP_dec"))]
  
  # Optional EP
  if ("padj_EP_ukb" %in% names(m) || "padj_EP_dec" %in% names(m)) {
    m[, padj_EP := pmin_na(get_or_na(m, "padj_EP_ukb"), get_or_na(m, "padj_EP_dec"))]
    m[, beta_EP := pick_by_best_p(get_or_na(m, "padj_EP_ukb"), get_or_na(m, "beta_EP_ukb"),
                                  get_or_na(m, "padj_EP_dec"), get_or_na(m, "beta_EP_dec"))]
  }
  
  keep <- c("trip_key","Exposure","Protein","Disease",
            "in_ukb","in_dec",
            "rep_PDcis","rep_PDtrans","rep_DP",
            "ukb_PDcis","ukb_PDtrans","ukb_DP",
            "dec_PDcis","dec_PDtrans","dec_DP",
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
    base_size = 10,
    
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
  
  # ---- MR edge significance category (priority: PDcis > PDtrans > DP) ----
  dt[, sig_PDcis   := !is.na(padj_PDcis)   & padj_PDcis   < mr_alpha]
  dt[, sig_PDtrans := !is.na(padj_PDtrans) & padj_PDtrans < mr_alpha]
  dt[, sig_DP      := !is.na(padj_DP)      & padj_DP      < mr_alpha]
  
  dt[, mr_edge_sig := fifelse(sig_PDcis, "PDcis",
                              fifelse(sig_PDtrans, "PDtrans",
                                      fifelse(sig_DP, "DP", "None")))]
  dt[, mr_edge_sig := factor(mr_edge_sig, levels = c("None","PDcis","PDtrans","DP"))]
  
  # ---- MR support (SHAPE): require SAME edge-type replication for "Both" ----
  # If mr_edge_sig is PDcis, "Both" only if rep_PDcis==TRUE, etc.
  # If edges differ across datasets, label as single-source based on which dataset has the smaller p for THIS edge type.
  dt[, mr_support := "None"]
  
  # helper: choose which dataset supports the chosen edge type better (smaller p)
  choose_side <- function(p_ukb, p_dec) {
    ifelse(!is.na(p_ukb) & (is.na(p_dec) | p_ukb <= p_dec), "UKB only",
           ifelse(!is.na(p_dec), "DECODE only", "None"))
  }
  
  # Build per-row dataset-specific p-values for the chosen edge type
  # (May be absent if you filtered columns upstream; safe if missing -> NA)
  # We do this by re-reading from MR_ukb/MR_dec is expensive; instead rely on availability of *_ukb/*_dec in combined
  # Here: we approximate "dataset support" using the summarized flags (ukb_PDcis/dec_PDcis etc.)
  # If you want p-based tie-break precisely, keep padj_*_ukb/_dec columns in MR_combined.
  dt[, mr_support := fifelse(mr_edge_sig == "PDcis" & rep_PDcis, "Both",
                             fifelse(mr_edge_sig == "PDtrans" & rep_PDtrans, "Both",
                                     fifelse(mr_edge_sig == "DP" & rep_DP, "Both",
                                             fifelse(mr_edge_sig == "PDcis" & ukb_PDcis & !dec_PDcis, "UKB only",
                                                     fifelse(mr_edge_sig == "PDcis" & dec_PDcis & !ukb_PDcis, "DECODE only",
                                                             fifelse(mr_edge_sig == "PDtrans" & ukb_PDtrans & !dec_PDtrans, "UKB only",
                                                                     fifelse(mr_edge_sig == "PDtrans" & dec_PDtrans & !ukb_PDtrans, "DECODE only",
                                                                             fifelse(mr_edge_sig == "DP" & ukb_DP & !dec_DP, "UKB only",
                                                                                     fifelse(mr_edge_sig == "DP" & dec_DP & !ukb_DP, "DECODE only",
                                                                                             "None")))))))))]
  
  # If there is MR evidence in both datasets but edges differ, assign "UKB only"/"DECODE only" based on which has the edge for the chosen color.
  # Remaining cases: keep None.
  dt[, mr_support := fct_explicit_na(as.factor(mr_support), na_level = "None")]
  dt[, mr_support := factor(mr_support, levels = c("None","UKB only","DECODE only","Both"))]
  
  # ---- Labels (priority: edge PDcis >> PDtrans >> DP >> None;
  #              within: Both >> single >> None) ----
  cap_n <- function(x, n) x[seq_len(min(length(x), n))]
  
  dt[, support_tier := fifelse(mr_support == "Both", 2L,
                               fifelse(mr_support %in% c("UKB only","DECODE only"), 1L, 0L))]
  
  dt[, edge_tier := fifelse(mr_edge_sig == "PDcis", 3L,
                            fifelse(mr_edge_sig == "PDtrans", 2L,
                                    fifelse(mr_edge_sig == "DP", 1L, 0L)))]
  
  dt[, label_rank := 1e6*edge_tier + 1e3*support_tier + 10*abs(beta_arm) + abs(beta_HEAP)]
  
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
      paste0("Exposure: ", EXPOSURE_TO_PLOT, " | Shape = edge-replication (same edge type in UKB+DECODE)")
    } else {
      paste0("Exposure: ", EXPOSURE_TO_PLOT, " | Disease: ", disease_for_arm,
             " | Shape = edge-replication (same edge type in UKB+DECODE)")
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
    
    geom_point(aes(color = mr_edge_sig, shape = mr_support, size = size_r)) +
    
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
    
    scale_size_continuous(
      range  = c(1.6, 6.0),
      limits = c(0.05, 1),
      breaks = c(0.05, 0.25, 0.5, 0.75, 1),
      labels = c("≤0 / NA", "0.25", "0.5", "0.75", "1.0"),
      name   = "SomaScan–Olink\nr (pos only)"
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
    
    scale_shape_manual(
      values = c(
        "None"       = 16,
        "UKB only"   = 17,
        "DECODE only"= 15,
        "Both"       = 18
      ),
      name = "MR support\n(edge-replicated)"
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
      plot.subtitle = element_text(size = base_size - 2)
    ) +
    coord_cartesian(clip = "off")
  
  return(p)
}


# ============================================================
# Example usage (same as before)
# ============================================================

suppressPackageStartupMessages({
  library(cowplot)
})

EXPOSURE_TO_PLOT <- "types_of_physical_activity_in_last_4_weeks_f6164_0_0.multi_Strenuous_sports"

p_master <- plot_HEAP_GLP1_MR_onepanel(
  EXPOSURE_TO_PLOT = EXPOSURE_TO_PLOT,
  arm = "HERITAGE",
  disease_for_arm = "finngen_R12_E4_OBESITY",
  legend_position = "right"
)

leg <- cowplot::get_legend(p_master + theme(legend.position = "right"))

ggsave(file.path(OUTDIR, "LEGEND_MR_GLP1.png"),
       plot = cowplot::ggdraw(leg),
       width = 3.0, height = 4.5, units = "in", dpi = 600, bg = "white")


p0 <- plot_HEAP_GLP1_MR_onepanel(
  EXPOSURE_TO_PLOT = EXPOSURE_TO_PLOT,
  arm = "HERITAGE",
  disease_for_arm = "finngen_R12_E4_OBESITY",
  title = "Exercise → proteins vs HERITAGE (Endurance Exercise)",
  subtitle = "Strenuous sports (UKB) | Shape = MR replication (same edge type UKB+DECODE)",
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
  subtitle = "Strenuous sports (UKB) | Shape = MR replication (same edge type UKB+DECODE)",
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
  subtitle = "Strenuous sports (UKB) | Shape = MR replication (same edge type UKB+DECODE)",
  xlab = "UKB assoc (E→P) beta",
  ylab = "GLP1_2 protein shift",
  legend_position = "right",
  corr_loc = "topleft"
)

print(p0); print(p1); print(p2)

save_no_legend <- function(p, filename, w=5, h=4, dpi=1000) {
  ggsave(file.path(OUTDIR, filename),
         plot = p + theme(legend.position = "none")
           ,
         width = w, height = h, units = "in", dpi = dpi, bg = "white")
}

save_no_legend(p0, "MRHERITAGE_StrenSports_noLegend.png")
save_no_legend(p1, "MRGLP1_STEP1_StrenSports_noLegend.png")
save_no_legend(p2, "MRGLP1_STEP2_StrenSports_noLegend.png")


ggsave(file.path(OUTDIR, "MRHERITAGE_StrenSports_edgeRepShape.png"),
       plot = p0, dpi = 1000, width = 7, height = 4, units = "in")
ggsave(file.path(OUTDIR, "MRGLP1_STEP1_StrenSports_edgeRepShape.png"),
       plot = p1, dpi = 1000, width = 7, height = 4, units = "in")
ggsave(file.path(OUTDIR, "MRGLP1_STEP2_StrenSports_edgeRepShape.png"),
       plot = p2, dpi = 1000, width = 7, height = 4, units = "in")

# Diet
p3 <- plot_HEAP_GLP1_MR_onepanel(
  EXPOSURE_TO_PLOT = "fresh_fruit_intake_f1309_0_0",
  arm = "GLP1_1",
  disease_for_arm = "finngen_R12_E4_OBESITY",
  title = "Diet → proteins vs GLP1_1 (obesity cohort)",
  subtitle = "Fruit intake (UKB) | Shape = MR replication (same edge type UKB+DECODE)",
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
  subtitle = "Fruit intake (UKB) | Shape = MR replication (same edge type UKB+DECODE)",
  xlab = "UKB assoc (E→P) beta",
  ylab = "GLP1_2 protein shift",
  legend_position = "right",
  corr_loc = "topleft"
)

print(p3); print(p4)

save_no_legend(p3, "MRGLP1_STEP1_FruitIntake_noLegend.png")
save_no_legend(p4, "MRGLP1_STEP2_FruitIntake_noLegend.png")

ggsave(file.path(OUTDIR, "MRGLP1_STEP1_FruitIntake_edgeRepShape.png"),
       plot = p3, dpi = 1000, width = 7, height = 4, units = "in")
ggsave(file.path(OUTDIR, "MRGLP1_STEP2_FruitIntake_edgeRepShape.png"),
       plot = p4, dpi = 1000, width = 7, height = 4, units = "in")

# Smoking
p5 <- plot_HEAP_GLP1_MR_onepanel(
  EXPOSURE_TO_PLOT = "past_tobacco_smoking_f1249_0_0",
  arm = "GLP1_1",
  disease_for_arm = "finngen_R12_E4_OBESITY",
  title = "Quit Smoking → proteins vs GLP1_1 (obesity cohort)",
  subtitle = "Former smoker | Shape = MR replication (same edge type UKB+DECODE)",
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
  subtitle = "Former smoker | Shape = MR replication (same edge type UKB+DECODE)",
  xlab = "UKB assoc (E→P) beta",
  ylab = "GLP1_2 protein shift",
  legend_position = "right",
  corr_loc = "topleft"
)

print(p5); print(p6)

save_no_legend(p5, "MRGLP1_STEP1_QuitSmoking_noLegend.png")
save_no_legend(p6, "MRGLP1_STEP2_QuitSmoking_noLegend.png")


ggsave(file.path(OUTDIR, "MRGLP1_STEP1_QuitSmoking_edgeRepShape.png"),
       plot = p5, dpi = 1000, width = 7, height = 4, units = "in")
ggsave(file.path(OUTDIR, "MRGLP1_STEP2_QuitSmoking_edgeRepShape.png"),
       plot = p6, dpi = 1000, width = 7, height = 4, units = "in")

cat("\nDone. Wrote plots to:\n  ", OUTDIR, "\n\n", sep = "")