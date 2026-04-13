#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(scales)
  library(pbapply)
  library(ggforce)
  library(grid)
  library(tibble)
})

# ============================================================
# PATHS
# ============================================================

MRfiles <- "/n/groups/patel/shakson_ukb/UK_Biobank/Output/MR_edges"
PLOTBASE <- "/n/groups/patel/shakson_ukb/UK_Biobank/Output/MRplots"

TRIPLET_UKB_FP <- file.path(MRfiles, "summary", "MRmotifs.csv")
TRIPLET_DEC_FP <- file.path(MRfiles, "summary", "DECODE", "MRmotifs.csv")

# sensitivity output from your script
SENS_FP <- file.path(PLOTBASE, "SENSITIVITY_UKB_vs_DECODE", "MR_sensitivity_table_ALL_edges_UKB_and_DECODE.tsv")

# outputs
OUTDIR_COMPARE <- file.path(PLOTBASE, "COMPARE_UKB_vs_DECODE")
OUTDIR_SENS_TRIPLET <- file.path(PLOTBASE, "COMPARE_UKB_vs_DECODE", "Sensitivity_TripletLevel")
OUTDIR_SHARED_PATHS <- file.path(PLOTBASE, "COMPARE_UKB_vs_DECODE", "SharedPaths_Motif_AtoE")

dir.create(OUTDIR_COMPARE, recursive = TRUE, showWarnings = FALSE)
dir.create(OUTDIR_SENS_TRIPLET, recursive = TRUE, showWarnings = FALSE)
dir.create(OUTDIR_SHARED_PATHS, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# GLOBAL SETTINGS
# ============================================================

KEEP_ONLY_ANY_SIG <- FALSE            # for motif overlap plots
KEEP_ONLY_ANY_SIG_PER_DATASET <- TRUE # for shared paths plots (recommended)
ALPHA <- 0.05                         # MRmotifs padj threshold (used for "sig edge" calls in triplets)

# For shared-path diagram batching
REQUIRE_SHARED_SIG_EDGE <- TRUE
REQUIRE_ANY_SIG_BOTH <- TRUE
MAX_TRIPLETS_PER_MOTIF <- Inf
SHOW_EDGES  <- "all"        # "all" or "sig_only"
LABEL_MODE  <- "sig_only"   # "sig_only" / "all" / "none"
PD_MODE <- "both"           # "cis" / "trans" / "both"
PE_MODE <- "both"           # "cis" / "trans" / "both"

# ============================================================
# UTILS
# ============================================================

make_trip_key <- function(df) {
  df %>%
    mutate(
      trip_key = ifelse(!is.na(triplet) & triplet != "",
                        as.character(triplet),
                        paste(Exposure, Protein, Disease, sep="||"))
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

# edge_dir mapping used by sensitivity script
edge_dirs <- c("E_to_P","Pcis_to_E","Ptrans_to_E",
               "Pcis_to_D","Ptrans_to_D","D_to_P",
               "E_to_D","D_to_E")

# For a triplet (E,P,D), define which (edge_dir, src, tgt) tuples correspond
triplet_edge_keys <- function(E, P, D) {
  data.table(
    edge_dir = edge_dirs,
    src_id = c(E, P, P, P, P, D, E, D),
    tgt_id = c(P, E, E, D, D, P, D, E)
  )
}

# ============================================================
# LOAD INPUTS
# ============================================================

message("Reading MRmotifs UKB + DECODE ...")
ukb <- fread(TRIPLET_UKB_FP) %>% as_tibble()
dec <- fread(TRIPLET_DEC_FP) %>% as_tibble()

ukb <- ukb %>% mutate(dataset="UKB") %>% make_trip_key() %>% assign_motif_letter()
dec <- dec %>% mutate(dataset="DECODE") %>% make_trip_key() %>% assign_motif_letter()

# Optionally filter for overlap plots
if (KEEP_ONLY_ANY_SIG) {
  if ("any_sig" %in% names(ukb)) ukb <- ukb %>% filter(any_sig)
  if ("any_sig" %in% names(dec)) dec <- dec %>% filter(any_sig)
}

# Dedup
ukb <- ukb %>% distinct(trip_key, .keep_all = TRUE)
dec <- dec %>% distinct(trip_key, .keep_all = TRUE)

message("Reading sensitivity table ...")
sens <- fread(SENS_FP)
sens <- as.data.table(sens)

# Keep just what we need; the rest is huge
keep_cols <- intersect(
  c("dataset","edge_dir","src_id","tgt_id","method","nsnp","pval","pval_adj","mr_hit",
    "het_pval","egger_pval","het_flag","pleio_flag","sens_pass","hit_after_sens"),
  names(sens)
)
sens <- sens[, ..keep_cols]

# Collapse to one row per edge (dataset, edge_dir, src_id, tgt_id):
# - any hit_after_sens across IVW/Wald rows
# - any mr_hit
# - and we keep min pval_adj (so you can rank if needed)
sens_edge <- sens[, .(
  any_mr_hit = any(mr_hit %in% TRUE, na.rm=TRUE),
  any_hit_after = any(hit_after_sens %in% TRUE, na.rm=TRUE),
  any_het_flag = any(het_flag %in% TRUE, na.rm=TRUE),
  any_pleio_flag = any(pleio_flag %in% TRUE, na.rm=TRUE),
  min_padj = suppressWarnings(min(pval_adj, na.rm=TRUE))
), by = .(dataset, edge_dir, src_id, tgt_id)]
sens_edge[is.infinite(min_padj), min_padj := NA_real_]

setkeyv(sens_edge, c("dataset","edge_dir","src_id","tgt_id"))

# ============================================================
# PART 1) MOTIF OVERLAP PLOTS (your current script, packaged)
# ============================================================

message("Building motif overlap / transition plots ...")

# Inner join on triplets that exist in both
both <- ukb %>%
  select(trip_key,
         Exposure, Protein, Disease,
         motif_label_ukb = motif_label,
         any_sig_ukb = any_sig,
         n_motifs_ukb = n_motifs,
         motif_A_ukb = motif_A_mediator,
         motif_B_ukb = motif_B_biomarker,
         motif_C_ukb = motif_C_exposure_marker,
         motif_D_ukb = motif_D_P_to_E,
         motif_E_ukb = motif_E_disease_liability) %>%
  inner_join(
    dec %>%
      select(trip_key,
             motif_label_dec = motif_label,
             any_sig_dec = any_sig,
             n_motifs_dec = n_motifs,
             motif_A_dec = motif_A_mediator,
             motif_B_dec = motif_B_biomarker,
             motif_C_dec = motif_C_exposure_marker,
             motif_D_dec = motif_D_P_to_E,
             motif_E_dec = motif_E_disease_liability),
    by = "trip_key"
  ) %>%
  mutate(same_label = (motif_label_ukb == motif_label_dec))

write.csv(both, file.path(OUTDIR_COMPARE, "Triplet_motif_calls_UKB_vs_DECODE_from_MRmotifs.csv"),
          row.names = FALSE)

labels_all <- sort(unique(c(both$motif_label_ukb, both$motif_label_dec)))
labels_all <- labels_all[!is.na(labels_all)]

motif_label_counts <- tibble(motif_label = labels_all) %>%
  mutate(
    n_UKB = sapply(labels_all, \(m) sum(both$motif_label_ukb == m, na.rm=TRUE)),
    n_DECODE = sapply(labels_all, \(m) sum(both$motif_label_dec == m, na.rm=TRUE)),
    n_shared_same = sapply(labels_all, \(m) sum(both$motif_label_ukb == m & both$motif_label_dec == m, na.rm=TRUE)),
    n_unique_UKB = sapply(labels_all, \(m) sum(both$motif_label_ukb == m & both$motif_label_dec != m, na.rm=TRUE)),
    n_unique_DECODE = sapply(labels_all, \(m) sum(both$motif_label_dec == m & both$motif_label_ukb != m, na.rm=TRUE))
  ) %>%
  arrange(desc(n_shared_same))

write.csv(motif_label_counts, file.path(OUTDIR_COMPARE, "MotifLabel_overlap_unique_counts.csv"),
          row.names = FALSE)

motif_label_xtab <- both %>%
  count(motif_label_ukb, motif_label_dec, name = "n") %>%
  complete(motif_label_ukb = labels_all, motif_label_dec = labels_all, fill = list(n = 0))

write.csv(motif_label_xtab, file.path(OUTDIR_COMPARE, "MotifLabel_crosstab_UKB_to_DECODE.csv"),
          row.names = FALSE)

positions <- c("Null (no MR hits)", "Other (has hits)",
               "E","D","C","B","A","Multiple motifs")

MotifUnion <- motif_label_counts %>%
  mutate(
    n_union = n_shared_same + n_unique_UKB + n_unique_DECODE,
    pct = n_union / sum(n_union),
    label = paste0(percent(pct, accuracy = 0.01), " (n=", comma(n_union), ")")
  )

p_union <- ggplot(MotifUnion, aes(x = reorder(motif_label, n_union), y = n_union)) +
  geom_col() +
  geom_text(aes(label = label), hjust = -0.1, size = 3) +
  coord_flip(clip = "off") +
  scale_y_log10(labels = comma, expand = expansion(mult = c(0, 0.28))) +
  scale_x_discrete(limits = positions) +
  theme_classic() +
  labs(x = NULL, y = "Triad count (log10)") +
  theme(plot.margin = margin(t=10, r=40, b=10, l=10))

ggsave(file.path(OUTDIR_COMPARE, "MotifLabel_Union_DECODEUKB.png"),
       p_union, width = 6, height = 4, units = "in", dpi = 1000)

comp_long <- motif_label_counts %>%
  select(motif_label, n_shared_same, n_unique_UKB, n_unique_DECODE) %>%
  pivot_longer(-motif_label, names_to = "component", values_to = "n") %>%
  mutate(component = recode(component,
                            n_shared_same   = "Shared",
                            n_unique_UKB    = "UKB pQTLs",
                            n_unique_DECODE = "deCODE pQTLs")) %>%
  mutate(pct = n / sum(n),
         label = paste0(percent(pct, 0.01), " (n=", comma(n), ")")) %>%
  ungroup()

p_components <- ggplot(
  comp_long,
  aes(x = reorder(motif_label, n, FUN = sum), y = n + 1, fill = component)
) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(aes(label = label),
            position = position_dodge(width = 0.8),
            hjust = -0.1, size = 3) +
  coord_flip(clip = "off") +
  scale_y_log10(labels = comma, expand = expansion(mult = c(0, 0.30))) +
  scale_x_discrete(limits = positions) +
  theme_classic() +
  labs(x = NULL, y = "Triad count (log10)", title = "Shared vs Unique Triad Hits (UKB vs deCODE)") +
  theme(legend.position = "bottom",
        plot.margin = margin(t = 10, r = 45, b = 10, l = 10))

ggsave(file.path(OUTDIR_COMPARE, "MotifLabel_Separate_DECODEUKB.png"),
       p_components, width = 7, height = 5, units = "in", dpi = 1000)

p2 <- ggplot(motif_label_xtab, aes(x = motif_label_dec, y = motif_label_ukb, fill = n)) +
  geom_tile() +
  theme_bw() +
  labs(x = "deCODE motif label", y = "UKB motif label",
       title = "Motif-label transitions (UKB → deCODE)", fill = "# triplets") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(OUTDIR_COMPARE, "MotifLabel_transition_heatmap.png"),
       p2, width = 8.5, height = 7, units = "in", dpi = 400)

motif_bool_summary <- tibble(
  motif_type = c("A_mediator","B_biomarker","C_exposure_marker","D_P_to_E","E_disease_liability"),
  n_UKB = c(sum(both$motif_A_ukb), sum(both$motif_B_ukb), sum(both$motif_C_ukb), sum(both$motif_D_ukb), sum(both$motif_E_ukb)),
  n_DECODE = c(sum(both$motif_A_dec), sum(both$motif_B_dec), sum(both$motif_C_dec), sum(both$motif_D_dec), sum(both$motif_E_dec)),
  n_shared = c(sum(both$motif_A_ukb & both$motif_A_dec),
               sum(both$motif_B_ukb & both$motif_B_dec),
               sum(both$motif_C_ukb & both$motif_C_dec),
               sum(both$motif_D_ukb & both$motif_D_dec),
               sum(both$motif_E_ukb & both$motif_E_dec)),
  n_unique_UKB = c(sum(both$motif_A_ukb & !both$motif_A_dec),
                   sum(both$motif_B_ukb & !both$motif_B_dec),
                   sum(both$motif_C_ukb & !both$motif_C_dec),
                   sum(both$motif_D_ukb & !both$motif_D_dec),
                   sum(both$motif_E_ukb & !both$motif_E_dec)),
  n_unique_DECODE = c(sum(both$motif_A_dec & !both$motif_A_ukb),
                      sum(both$motif_B_dec & !both$motif_B_ukb),
                      sum(both$motif_C_dec & !both$motif_C_ukb),
                      sum(both$motif_D_dec & !both$motif_D_ukb),
                      sum(both$motif_E_dec & !both$motif_E_ukb))
) %>%
  mutate(jaccard = n_shared / (n_UKB + n_DECODE - n_shared))

write.csv(motif_bool_summary, file.path(OUTDIR_COMPARE, "MotifTypeBoolean_overlap_unique_counts.csv"),
          row.names = FALSE)

# ============================================================
# PART 2) NEW: SENSITIVITY AT TRIPLET LEVEL + PLOTS BY MOTIF
# ============================================================

message("Building triplet-level sensitivity summaries...")

# Work on the intersection triplets (like your overlap script), to keep apples-to-apples
shared_trip_keys <- intersect(ukb$trip_key, dec$trip_key)

ukb_shared <- ukb %>% filter(trip_key %in% shared_trip_keys)
dec_shared <- dec %>% filter(trip_key %in% shared_trip_keys)

trip_master <- bind_rows(ukb_shared, dec_shared) %>%
  select(dataset, trip_key, Exposure, Protein, Disease, motif_label, motif_letter, any_sig) %>%
  distinct(dataset, trip_key, .keep_all = TRUE)

trip_master_dt <- as.data.table(trip_master)

# Attach triplet->edge sensitivity info:
# For each triplet, we look up the 8 directed edges (E_to_P, Pcis_to_E, ...)
# and summarize:
#   - any edge is MR hit (any_mr_hit)
#   - any edge survives sensitivity (any_hit_after)
#   - any heterogeneity flagged (any_het_flag) among those edges
#   - any pleiotropy flagged (any_pleio_flag) among those edges
trip_sens_list <- pbapply::pblapply(seq_len(nrow(trip_master_dt)), function(i) {
  row <- trip_master_dt[i]
  E <- row$Exposure; P <- row$Protein; D <- row$Disease
  ds <- row$dataset
  
  keys <- triplet_edge_keys(E, P, D)
  keys[, dataset := ds]
  setkeyv(keys, c("dataset","edge_dir","src_id","tgt_id"))
  
  # join sensitivity per edge
  joined <- sens_edge[keys]
  
  # summarize per triplet
  out <- data.table(
    dataset = ds,
    trip_key = row$trip_key,
    any_edge_mr_hit = any(joined$any_mr_hit %in% TRUE, na.rm=TRUE),
    any_edge_hit_after = any(joined$any_hit_after %in% TRUE, na.rm=TRUE),
    any_edge_het_flag = any(joined$any_het_flag %in% TRUE, na.rm=TRUE),
    any_edge_pleio_flag = any(joined$any_pleio_flag %in% TRUE, na.rm=TRUE)
  )
  
  # retention among hits: if no edge hit, leave NA
  if (out$any_edge_mr_hit) {
    out[, triplet_retained := any_edge_hit_after]
  } else {
    out[, triplet_retained := NA]
  }
  
  out
})

trip_sens <- rbindlist(trip_sens_list, fill=TRUE)

trip_sens <- merge(trip_master_dt, trip_sens, by=c("dataset","trip_key"), all.x=TRUE)

# Write the big triplet-level sensitivity table
fwrite(trip_sens, file.path(OUTDIR_SENS_TRIPLET, "TripletLevel_sensitivity_UKB_vs_DECODE.tsv"), sep="\t")

# ---- Plot: retention by motif_letter (A-E) ----
tmp <- trip_sens[!is.na(motif_letter)]

ret_by_letter <- tmp[, .(
  denom = sum(any_edge_mr_hit %in% TRUE, na.rm=TRUE),
  num = sum(any_edge_mr_hit %in% TRUE & triplet_retained %in% TRUE, na.rm=TRUE)
), by=.(dataset, motif_letter)]

ret_by_letter[, pct := ifelse(denom > 0, 100*num/denom, NA_real_)]
ret_by_letter[, lab := paste0(num, "/", denom)]

p_ret_letter <- ggplot(ret_by_letter, aes(x=motif_letter, y=pct, fill=dataset)) +
  geom_col(position = position_dodge(width=0.8), width=0.7) +
  geom_text(aes(label=lab),
            position = position_dodge(width=0.8),
            vjust = -0.25, size=3) +
  theme_bw() +
  labs(
    x = "Motif letter",
    y = "% triplets retained after sensitivity (among triplets with any MR-hit edge)",
    title = "Triplet-level MR retention after sensitivity filtering (UKB vs DECODE)"
  ) +
  theme(legend.title = element_blank()) +
  coord_cartesian(clip="off")

ggsave(file.path(OUTDIR_SENS_TRIPLET, "Triplet_retention_by_motifLetter.png"),
       p_ret_letter, width=8, height=4.5, dpi=600)

# ---- Plot: violation rate by motif_letter among triplets with any MR-hit edge ----
viol_long <- rbindlist(list(
  tmp[, .(dataset, motif_letter, flag = any_edge_het_flag, which="Heterogeneity")],
  tmp[, .(dataset, motif_letter, flag = any_edge_pleio_flag, which="Pleiotropy")]
), fill=TRUE)

viol_sum <- viol_long[, .(
  denom = .N,
  num = sum(flag %in% TRUE, na.rm=TRUE)
), by=.(dataset, motif_letter, which)]

viol_sum[, pct := ifelse(denom>0, 100*num/denom, NA_real_)]
viol_sum[, lab := paste0(num, "/", denom)]

p_viol_letter <- ggplot(viol_sum, aes(x=motif_letter, y=pct, fill=dataset)) +
  geom_col(position = position_dodge(width=0.8), width=0.7) +
  geom_text(aes(label=lab),
            position = position_dodge(width=0.8),
            vjust = -0.25, size=3) +
  facet_wrap(~which, ncol=1) +
  theme_bw() +
  labs(
    x = "Motif letter",
    y = "% triplets flagged (triplet has ANY flagged edge)",
    title = "Triplet-level sensitivity flags by motif letter (UKB vs DECODE)"
  ) +
  theme(legend.title = element_blank()) +
  coord_cartesian(clip="off")

ggsave(file.path(OUTDIR_SENS_TRIPLET, "Triplet_flags_by_motifLetter.png"),
       p_viol_letter, width=8, height=7, dpi=600)

# Optional: do the same by motif_label (can be many; still useful)
if ("motif_label" %in% names(trip_sens)) {
  tmp2 <- trip_sens[!is.na(motif_label)]
  ret_by_label <- tmp2[, .(
    denom = sum(any_edge_mr_hit %in% TRUE, na.rm=TRUE),
    num = sum(any_edge_mr_hit %in% TRUE & triplet_retained %in% TRUE, na.rm=TRUE)
  ), by=.(dataset, motif_label)]
  ret_by_label[, pct := ifelse(denom > 0, 100*num/denom, NA_real_)]
  fwrite(ret_by_label, file.path(OUTDIR_SENS_TRIPLET, "Triplet_retention_by_motifLabel.tsv"), sep="\t")
}

# ============================================================
# PART 3) SHARED PATH DIAGRAMS (batch; OPTIONAL HEAVY)
# ============================================================

message("Preparing shared-path diagram batching...")

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
    scale_linetype_manual(values = c(`TRUE`="solid", `FALSE`="dashed")) +
    scale_alpha_manual(values = c(`TRUE`=1, `FALSE`=0.18)) +
    geom_point(data=nodes, aes(x=x, y=y), size=node_size, shape=21, stroke=1.6, fill="white") +
    geom_text(data=nodes, aes(x=x, y=y, label=node), size=node_text_size) +
    coord_equal(xlim=xlim, ylim=ylim, clip="off") +
    theme_void() +
    theme(
      plot.title = element_text(hjust=0.5, size=16, face="bold"),
      plot.subtitle = element_text(hjust=0.5, size=13),
      legend.position = if (legend) "bottom" else "none",
      plot.margin = margin(6, 6, 6, 6)
    ) +
    labs(title = title %||% triplet_id, subtitle = subtitle)
  
  if (nrow(edges_lab) > 0) {
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
  }
  
  g
}

# For diagram script compatibility: it filters on DT$triplet
ukb_diag <- fread(TRIPLET_UKB_FP) %>% as_tibble()
dec_diag <- fread(TRIPLET_DEC_FP) %>% as_tibble()

if (KEEP_ONLY_ANY_SIG_PER_DATASET) {
  if ("any_sig" %in% names(ukb_diag)) ukb_diag <- ukb_diag %>% filter(any_sig)
  if ("any_sig" %in% names(dec_diag)) dec_diag <- dec_diag %>% filter(any_sig)
}

ukb_diag <- ukb_diag %>% make_trip_key() %>% mutate(triplet = trip_key) %>% assign_motif_letter()
dec_diag <- dec_diag %>% make_trip_key() %>% mutate(triplet = trip_key) %>% assign_motif_letter()

ukb_diag <- ukb_diag %>% filter(motif_letter %in% c("A","B","C","D","E")) %>% distinct(triplet, .keep_all = TRUE)
dec_diag <- dec_diag %>% filter(motif_letter %in% c("A","B","C","D","E")) %>% distinct(triplet, .keep_all = TRUE)

# Shared triplets by motif letter
motifs <- c("A","B","C","D","E")

for (m in motifs) {
  message("Shared-path diagrams motif ", m, " ...")
  
  ukb_m <- ukb_diag %>% filter(motif_letter == m)
  dec_m <- dec_diag %>% filter(motif_letter == m)
  
  shared_triplets <- intersect(ukb_m$triplet, dec_m$triplet)
  if (length(shared_triplets) == 0) next
  
  # Optionally require any_sig in BOTH
  if (REQUIRE_ANY_SIG_BOTH) {
    ukb_m2 <- ukb_m %>% filter(triplet %in% shared_triplets, any_sig %in% TRUE)
    dec_m2 <- dec_m %>% filter(triplet %in% shared_triplets, any_sig %in% TRUE)
    shared_triplets <- intersect(ukb_m2$triplet, dec_m2$triplet)
  }
  
  if (length(shared_triplets) == 0) next
  
  # Optionally require at least one shared sig edge in BOTH (padj_* < alpha)
  if (REQUIRE_SHARED_SIG_EDGE) {
    has_any_sig_edge <- function(df_row) {
      pcols <- names(df_row)[grepl("^padj_", names(df_row))]
      if (length(pcols) == 0) return(FALSE)
      any(suppressWarnings(as.numeric(df_row[, pcols, drop=TRUE])) < ALPHA, na.rm=TRUE)
    }
    
    ok <- sapply(shared_triplets, function(tid) {
      urow <- ukb_m %>% filter(triplet == tid)
      drow <- dec_m %>% filter(triplet == tid)
      if (nrow(urow) != 1 || nrow(drow) != 1) return(FALSE)
      has_any_sig_edge(urow) && has_any_sig_edge(drow)
    })
    shared_triplets <- shared_triplets[ok]
  }
  
  if (is.finite(MAX_TRIPLETS_PER_MOTIF)) {
    shared_triplets <- head(shared_triplets, MAX_TRIPLETS_PER_MOTIF)
  }
  
  if (length(shared_triplets) == 0) next
  
  # Make subfolders
  dir.create(file.path(OUTDIR_SHARED_PATHS, m, "UKB"), recursive=TRUE, showWarnings=FALSE)
  dir.create(file.path(OUTDIR_SHARED_PATHS, m, "DECODE"), recursive=TRUE, showWarnings=FALSE)
  
  for (tid in shared_triplets) {
    # sanitize filename
    fn_safe <- gsub("[/\\\\:;\\*\\?\\\"\\<\\>\\|]", "_", tid)
    
    p_ukb <- plot_triplet_mr_diagram_nature(
      DT = ukb_m,
      triplet_id = tid,
      alpha = ALPHA,
      show = SHOW_EDGES,
      pd_mode = PD_MODE,
      pe_mode = PE_MODE,
      label_mode = LABEL_MODE,
      title = paste0("Motif ", m, " (UKB)"),
      subtitle = tid
    )
    p_dec <- plot_triplet_mr_diagram_nature(
      DT = dec_m,
      triplet_id = tid,
      alpha = ALPHA,
      show = SHOW_EDGES,
      pd_mode = PD_MODE,
      pe_mode = PE_MODE,
      label_mode = LABEL_MODE,
      title = paste0("Motif ", m, " (deCODE)"),
      subtitle = tid
    )
    
    ggsave(file.path(OUTDIR_SHARED_PATHS, m, "UKB", paste0(fn_safe, ".png")),
           p_ukb, width=8.5, height=5.5, units="in", dpi=600)
    ggsave(file.path(OUTDIR_SHARED_PATHS, m, "DECODE", paste0(fn_safe, ".png")),
           p_dec, width=8.5, height=5.5, units="in", dpi=600)
  }
}

message("\nDONE.\nWrote:\n  ", OUTDIR_COMPARE,
        "\n  ", OUTDIR_SENS_TRIPLET,
        "\n  ", OUTDIR_SHARED_PATHS, "\n", sep="")