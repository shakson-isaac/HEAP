
#!/usr/bin/env Rscript

# ============================================================
# Plot shared motif paths A–E between UKB and DECODE
# - For each motif letter (A,B,C,D,E):
#   - find triplets present in BOTH cohorts within that motif
#   - (optional) require any_sig in both cohorts
#   - (optional) require at least 1 shared significant edge (q<alpha) in both cohorts
#   - save 2 diagrams per triplet: UKB + DECODE
# - Outputs folders: OUTDIR/A/UKB, OUTDIR/A/DECODE, ..., OUTDIR/E/DECODE
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(ggforce)
  library(grid)
  library(tibble)
})

# ============================================================
# Helpers + plotting function (your exact function)
# ============================================================
`%||%` <- function(a,b) if (!is.null(a)) a else b

pstars <- function(p) {
  ifelse(is.na(p), "",
         ifelse(p < 0.001, "***",
                ifelse(p < 0.01, "**",
                       ifelse(p < 0.05, "*", ""))))
}

plot_triplet_mr_diagram_nature <- function(
    DT,
    triplet_id,
    motif = NULL,
    title = NULL,
    subtitle = NULL,
    alpha = 0.05,
    show = c("all","sig_only"),
    pd_mode = c("cis","trans","both"),
    pe_mode = c("cis","trans","both"),
    digits = 3,
    label_mode = c("sig_only","all","none"),
    label_box = TRUE,
    label_pad = 0.20,
    node_radius = 0.95,
    node_size = 32,
    node_text_size = 5,
    arrow_mm = 3.0,
    edge_lwd = 1.5,
    xlim = c(-0.6, 10.6),
    ylim = c(-0.9, 3.9),
    legend = TRUE
) {
  
  show <- match.arg(show)
  pd_mode <- match.arg(pd_mode)
  pe_mode <- match.arg(pe_mode)
  label_mode <- match.arg(label_mode)
  
  row <- DT %>% as.data.frame() %>% filter(triplet == triplet_id)
  if (!is.null(motif)) row <- row %>% filter(motif_label == motif)
  if (nrow(row) != 1) stop("Expected exactly 1 row after filtering; got n=", nrow(row))
  
  nodes <- tibble(
    node = c("Exposure","Protein","Disease"),
    x    = c(0, 5, 10),
    y    = c(0, 3, 0)
  )
  
  edges_def <- tibble(
    edge_type = c("EP","PEcis","PEtrans",
                  "PDcis","PDtrans","DP",
                  "ED","DE"),
    from      = c("Exposure","Protein","Protein",
                  "Protein","Protein","Disease",
                  "Exposure","Disease"),
    to        = c("Protein","Exposure","Exposure",
                  "Disease","Disease","Protein",
                  "Disease","Exposure"),
    curve_mag = c(0.00, 0.34, 0.46,
                  0.00, 0.00, 0.34,
                  0.00, 0.24),
    side      = c(0,  -1,  -1,
                  0,   0,  -1,
                  0,  -1),
    ly_nudge  = c(+0.22, -0.28, -0.36,
                  +0.22, +0.22, +0.26,
                  -0.22, +0.22),
    
    beta_col  = c("beta_EP","beta_PEcis","beta_PEtrans",
                  "beta_PDcis","beta_PDtrans","beta_DP",
                  "beta_ED","beta_DE"),
    se_col    = c("se_EP","se_PEcis","se_PEtrans",
                  "se_PDcis","se_PDtrans","se_DP",
                  "se_ED","se_DE"),
    p_col     = c("padj_EP","padj_PEcis","padj_PEtrans",
                  "padj_PDcis","padj_PDtrans","padj_DP",
                  "padj_ED","padj_DE")
  )
  
  if (pd_mode != "both") {
    edges_def <- edges_def %>%
      filter(!(edge_type %in% c("PDcis","PDtrans")) | edge_type == paste0("PD", pd_mode))
  }
  if (pe_mode != "both") {
    edges_def <- edges_def %>%
      filter(!(edge_type %in% c("PEcis","PEtrans")) | edge_type == paste0("PE", pe_mode))
  }
  
  edges <- edges_def %>%
    rowwise() %>%
    mutate(beta = row[[beta_col]], se = row[[se_col]], padj = row[[p_col]]) %>%
    ungroup() %>%
    mutate(
      sig = !is.na(padj) & padj < alpha,
      lo  = ifelse(is.na(beta) | is.na(se), NA_real_, beta - 1.96 * se),
      hi  = ifelse(is.na(beta) | is.na(se), NA_real_, beta + 1.96 * se),
      label = case_when(
        is.na(beta) ~ "",
        is.na(se)   ~ paste0("β=", formatC(beta, format="f", digits=digits), pstars(padj)),
        TRUE        ~ paste0(
          "β=", formatC(beta, format="f", digits=digits),
          " [", formatC(lo, format="f", digits=digits), ", ",
          formatC(hi, format="f", digits=digits), "]",
          pstars(padj)
        )
      )
    )
  
  if (show == "sig_only") edges <- edges %>% filter(sig)
  
  edges <- edges %>%
    left_join(nodes %>% rename(from=node, x_from=x, y_from=y), by="from") %>%
    left_join(nodes %>% rename(to=node,   x_to=x,   y_to=y), by="to") %>%
    mutate(
      dx = x_to - x_from,
      dy = y_to - y_from,
      L  = sqrt(dx^2 + dy^2),
      x_from2 = x_from + node_radius * dx / L,
      y_from2 = y_from + node_radius * dy / L,
      x_to2   = x_to   - node_radius * dx / L,
      y_to2   = y_to   - node_radius * dy / L,
      x_mid = (x_from2 + x_to2)/2,
      y_mid = (y_from2 + y_to2)/2
    )
  
  edges_lab <- edges
  if (label_mode == "none") edges_lab <- edges_lab %>% filter(FALSE)
  if (label_mode == "sig_only") edges_lab <- edges_lab %>% filter(sig)
  
  edges_straight <- edges %>% filter(curve_mag == 0)
  edges_curved   <- edges %>% filter(curve_mag > 0)
  
  if (nrow(edges_curved) > 0) {
    edges_curved <- edges_curved %>%
      mutate(
        dx2 = x_to2 - x_from2,
        dy2 = y_to2 - y_from2,
        L2  = sqrt(dx2^2 + dy2^2),
        ux  = -dy2 / L2,
        uy  =  dx2 / L2,
        bend = curve_mag * side * L2,
        x_ctrl = x_mid + bend * ux,
        y_ctrl = y_mid + bend * uy
      )
    
    bez <- bind_rows(
      edges_curved %>% transmute(edge_type, sig, t=1, x=x_from2, y=y_from2),
      edges_curved %>% transmute(edge_type, sig, t=2, x=x_ctrl,  y=y_ctrl),
      edges_curved %>% transmute(edge_type, sig, t=3, x=x_to2,   y=y_to2)
    )
  } else {
    bez <- tibble(edge_type=character(), sig=logical(), t=integer(), x=double(), y=double())
  }
  
  g <- ggplot() +
    geom_segment(
      data = edges_straight,
      aes(x=x_from2, y=y_from2, xend=x_to2, yend=y_to2,
          color=edge_type, linetype=sig, alpha=sig),
      arrow = arrow(type="closed", length = unit(arrow_mm, "mm")),
      linewidth = edge_lwd
    ) +
    ggforce::geom_bezier(
      data = bez,
      aes(x=x, y=y, group=edge_type,
          color=edge_type, linetype=sig, alpha=sig),
      arrow = arrow(type="closed", length = unit(arrow_mm, "mm")),
      linewidth = edge_lwd
    ) +
    scale_linetype_manual(
      values = c(`TRUE`="solid", `FALSE`="dashed"),
      breaks = "TRUE",
      labels = "TRUE",
      name = paste0("Significant (q<", alpha, ")")
    ) +
    scale_alpha_manual(
      values = c(`TRUE`=1, `FALSE`=0.18),
      breaks = "TRUE",
      labels = "TRUE",
      name = paste0("Significant (q<", alpha, ")")
    )
  
  if (label_box && nrow(edges_lab) > 0) {
    g <- g +
      geom_label(
        data = edges_lab,
        aes(x=x_mid, y=y_mid + ly_nudge, label=label),
        size = 5,
        label.size = 0.25,
        label.padding = unit(label_pad, "lines"),
        fill = "white",
        alpha = 0.98
      )
  } else if (!label_box && nrow(edges_lab) > 0) {
    g <- g +
      geom_text(
        data = edges_lab,
        aes(x=x_mid, y=y_mid + ly_nudge, label=label),
        size = 5
      )
  }
  
  g <- g +
    geom_point(data=nodes, aes(x=x, y=y),
               size=node_size, shape=21, stroke=1.6, fill="white") +
    geom_text(data=nodes, aes(x=x, y=y, label=node), size=node_text_size) +
    coord_equal(xlim=xlim, ylim=ylim, clip="off") +
    theme_void() +
    theme(
      plot.title = element_text(hjust=0.5, size=16, face="bold"),
      plot.subtitle = element_text(hjust=0.5, size=13),
      legend.position = if (legend) "bottom" else "none",
      legend.box = "vertical",
      legend.title = element_text(size=12),
      legend.text  = element_text(size=11),
      plot.margin = margin(6, 6, 6, 6)
    ) +
    guides(
      color = guide_legend(title="Edge type", nrow=2, override.aes=list(alpha=1, linetype="solid")),
      linetype = guide_legend(
        title = paste0("Significant (q<", alpha, ")"),
        override.aes = list(color="black", alpha=1)
      ),
      alpha = "none"
    ) +
    labs(title = title %||% triplet_id, subtitle = subtitle)
  
  g
}

# ============================================================
# I/O
# ============================================================
MRfiles <- "/n/groups/patel/shakson_ukb/UK_Biobank/Output/MR_edges"
TRIPLET_UKB_FP <- file.path(MRfiles, "summary", "MRmotifs.csv")
TRIPLET_DEC_FP <- file.path(MRfiles, "summary", "DECODE", "MRmotifs.csv")

output_dir <- "/n/groups/patel/shakson_ukb/UK_Biobank/Output/MRplots/"
OUTDIR <- file.path(output_dir, "COMPARE_UKB_vs_DECODE", "SharedPaths_Motif_AtoE")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# Settings
# ============================================================
alpha   <- 0.05
pd_mode <- "both"   # "cis" / "trans" / "both"
pe_mode <- "both"   # "cis" / "trans" / "both"

# TRUE = require at least one edge significant in BOTH cohorts
REQUIRE_SHARED_SIG_EDGE <- TRUE

# TRUE = only consider triplets with any_sig in BOTH
REQUIRE_ANY_SIG_BOTH <- TRUE

# limit for debugging; set Inf for all
MAX_TRIPLETS_PER_MOTIF <- Inf

# plot aesthetics
SHOW_EDGES  <- "all"       # "all" or "sig_only"
LABEL_MODE  <- "sig_only"  # "sig_only" / "all" / "none"

# ============================================================
# Read
# ============================================================
ukb <- fread(TRIPLET_UKB_FP) %>% as_tibble()
dec <- fread(TRIPLET_DEC_FP) %>% as_tibble()

# ============================================================
# Build stable keys + motif letter
# ============================================================

make_trip_key <- function(df) {
  df %>%
    mutate(
      trip_key = ifelse(!is.na(triplet) & triplet != "",
                        as.character(triplet),
                        paste(Exposure, Protein, Disease, sep="||")),
      # plot function filters on DT$triplet
      triplet = trip_key
    )
}

assign_motif_letter <- function(df) {
  df <- df %>% mutate(motif_letter = NA_character_)
  
  if ("motif_label" %in% names(df)) {
    df <- df %>%
      mutate(motif_letter = str_extract(as.character(motif_label), "\\b[A-E]\\b"))
  }
  
  has_cols <- all(c("motif_A_mediator","motif_B_biomarker","motif_C_exposure_marker",
                    "motif_D_P_to_E","motif_E_disease_liability") %in% names(df))
  
  if (has_cols) {
    df <- df %>%
      mutate(
        motif_letter = ifelse(is.na(motif_letter) & (motif_A_mediator %in% c(TRUE,1)), "A", motif_letter),
        motif_letter = ifelse(is.na(motif_letter) & (motif_B_biomarker %in% c(TRUE,1)), "B", motif_letter),
        motif_letter = ifelse(is.na(motif_letter) & (motif_C_exposure_marker %in% c(TRUE,1)), "C", motif_letter),
        motif_letter = ifelse(is.na(motif_letter) & (motif_D_P_to_E %in% c(TRUE,1)), "D", motif_letter),
        motif_letter = ifelse(is.na(motif_letter) & (motif_E_disease_liability %in% c(TRUE,1)), "E", motif_letter)
      )
  }
  
  df
}

ukb <- ukb %>% make_trip_key() %>% assign_motif_letter()
dec <- dec %>% make_trip_key() %>% assign_motif_letter()

# keep only A-E
ukb <- ukb %>% filter(motif_letter %in% c("A","B","C","D","E"))
dec <- dec %>% filter(motif_letter %in% c("A","B","C","D","E"))

# de-dup by triplet within cohort (safest)
ukb <- ukb %>% distinct(triplet, .keep_all = TRUE)
dec <- dec %>% distinct(triplet, .keep_all = TRUE)

# ============================================================
# Determine which edge p-columns are "in play"
# ============================================================
edge_types <- c("EP","ED","DE","DP")
if (pd_mode == "both") edge_types <- c(edge_types, "PDcis","PDtrans")
if (pd_mode == "cis")  edge_types <- c(edge_types, "PDcis")
if (pd_mode == "trans")edge_types <- c(edge_types, "PDtrans")
if (pe_mode == "both") edge_types <- c(edge_types, "PEcis","PEtrans")
if (pe_mode == "cis")  edge_types <- c(edge_types, "PEcis")
if (pe_mode == "trans")edge_types <- c(edge_types, "PEtrans")

# ============================================================
# Output helpers
# ============================================================
safe_slug <- function(x) {
  x %>%
    str_replace_all("[^A-Za-z0-9]+", "_") %>%
    str_replace_all("_+", "_") %>%
    str_replace_all("^_|_$", "") %>%
    substr(1, 200)
}

# ============================================================
# Main loop: motif A-E
# ============================================================
all_summary <- list()

for (m in c("A","B","C","D","E")) {
  
  out_m <- file.path(OUTDIR, m)
  dir.create(file.path(out_m, "UKB"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(out_m, "DECODE"), recursive = TRUE, showWarnings = FALSE)
  
  ukb_m <- ukb %>% filter(motif_letter == m)
  dec_m <- dec %>% filter(motif_letter == m)
  
  shared_triplets <- intersect(unique(ukb_m$triplet), unique(dec_m$triplet))
  
  if (length(shared_triplets) == 0) {
    message(sprintf("[Motif %s] No shared triplets between UKB and DECODE.", m))
    next
  }
  
  # Build metadata table for shared triplets
  both <- ukb_m %>%
    select(triplet, Exposure, Protein, Disease,
           any_sig_ukb = any_sig,
           starts_with("padj_")) %>%
    inner_join(
      dec_m %>%
        select(triplet,
               any_sig_dec = any_sig,
               starts_with("padj_")),
      by = "triplet",
      suffix = c("_ukb", "_dec")
    )
  
  if (REQUIRE_ANY_SIG_BOTH && ("any_sig_ukb" %in% names(both)) && ("any_sig_dec" %in% names(both))) {
    both <- both %>% filter(any_sig_ukb, any_sig_dec)
  }
  
  # compute shared significant edges for each triplet
  shared_sig_tbl <- both %>%
    rowwise() %>%
    mutate(
      shared_edges = {
        keep <- c()
        for (e in edge_types) {
          pu_nm <- paste0("padj_", e, "_ukb")
          pd_nm <- paste0("padj_", e, "_dec")
          p_u <- if (pu_nm %in% names(cur_data())) get(pu_nm) else NA_real_
          p_d <- if (pd_nm %in% names(cur_data())) get(pd_nm) else NA_real_
          if (!is.na(p_u) && !is.na(p_d) && (p_u < alpha) && (p_d < alpha)) keep <- c(keep, e)
        }
        paste(keep, collapse = ",")
      },
      n_shared_edges = ifelse(shared_edges == "", 0L, length(strsplit(shared_edges, ",", fixed = TRUE)[[1]]))
    ) %>%
    ungroup()
  
  if (REQUIRE_SHARED_SIG_EDGE) {
    shared_sig_tbl <- shared_sig_tbl %>% filter(n_shared_edges > 0)
  }
  
  if (nrow(shared_sig_tbl) == 0) {
    message(sprintf("[Motif %s] Shared triplets exist, but none pass shared-edge filters.", m))
    next
  }
  
  shared_sig_tbl <- shared_sig_tbl %>% arrange(desc(n_shared_edges))
  
  # cap count if requested
  n_to_take <- min(
    nrow(shared_sig_tbl),
    if (is.finite(MAX_TRIPLETS_PER_MOTIF)) as.integer(MAX_TRIPLETS_PER_MOTIF) else nrow(shared_sig_tbl)
  )
  triplets_to_plot <- head(shared_sig_tbl$triplet, n_to_take)
  
  # save summary CSV for this motif
  motif_csv <- file.path(out_m, paste0("SharedTriplets_Motif_", m, "_summary.csv"))
  write.csv(shared_sig_tbl, motif_csv, row.names = FALSE)
  
  # store for combined summary
  all_summary[[m]] <- shared_sig_tbl %>% mutate(motif_letter = m)
  
  message(sprintf("[Motif %s] Plotting %d triplets.", m, length(triplets_to_plot)))
  
  # Plot each triplet: UKB + DECODE
  for (tid in triplets_to_plot) {
    
    row_meta <- shared_sig_tbl %>% filter(triplet == tid) %>% slice(1)
    
    short_title <- paste0(row_meta$Exposure, " → ", row_meta$Protein, " → ", row_meta$Disease)
    sub <- paste0(
      "Motif ", m, " | Shared sig edges (q<", alpha, "): ",
      ifelse(row_meta$shared_edges == "", "None", row_meta$shared_edges)
    )
    
    p_ukb <- plot_triplet_mr_diagram_nature(
      DT = ukb_m,               # motif-restricted DT so you only plot A-E motifs you asked for
      triplet_id = tid,
      title = paste0("UKB • ", short_title),
      subtitle = sub,
      alpha = alpha,
      show = SHOW_EDGES,
      pd_mode = pd_mode,
      pe_mode = pe_mode,
      label_mode = LABEL_MODE,
      legend = TRUE
    )
    
    p_dec <- plot_triplet_mr_diagram_nature(
      DT = dec_m,
      triplet_id = tid,
      title = paste0("DECODE • ", short_title),
      subtitle = sub,
      alpha = alpha,
      show = SHOW_EDGES,
      pd_mode = pd_mode,
      pe_mode = pe_mode,
      label_mode = LABEL_MODE,
      legend = TRUE
    )
    
    # filename uses E/P/D for readability
    fn_base <- safe_slug(paste0(row_meta$Exposure, "__", row_meta$Protein, "__", row_meta$Disease))
    
    ggsave(file.path(out_m, "UKB", paste0(fn_base, ".png")),
           p_ukb, width = 8.5, height = 5.5, units = "in", dpi = 450)
    
    ggsave(file.path(out_m, "DECODE", paste0(fn_base, ".png")),
           p_dec, width = 8.5, height = 5.5, units = "in", dpi = 450)
  }
}

# ============================================================
# Combined summary across motifs
# ============================================================
if (length(all_summary) > 0) {
  all_df <- bind_rows(all_summary)
  write.csv(all_df, file.path(OUTDIR, "SharedTriplets_AllMotifs_AtoE.csv"), row.names = FALSE)
  message("Wrote combined summary: ", file.path(OUTDIR, "SharedTriplets_AllMotifs_AtoE.csv"))
}

message("Done. Output root: ", OUTDIR)
