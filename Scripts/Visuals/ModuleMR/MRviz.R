suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(ggplot2)
  library(svglite)
  library(scales)
})

#Read in files from HEAP and MR:
MRfiles <- "/n/groups/patel/shakson_ukb/UK_Biobank/Output/MR_edges"

HEAP <- fread(file=file.path(MRfiles, "global_edges","HEAPres.tsv"))
PD <- fread(file=file.path(MRfiles, "summary","PDres.csv"))
EP <- fread(file=file.path(MRfiles, "summary","EPres.csv"))
ED <- fread(file=file.path(MRfiles, "summary","EDres.csv"))
DE <- fread(file=file.path(MRfiles, "summary","DEres.csv"))
PE <- fread(file=file.path(MRfiles, "summary","PEres.csv"))
DP <- fread(file=file.path(MRfiles, "summary","DPres.csv"))
  
#Filter for INW or Wald Ratio for each type of test
unique(PD$method)

EP <- EP %>% filter(method == "Inverse variance weighted" | 
                      method == "Wald ratio")
PDcis <- PD %>% 
            filter(method == "Inverse variance weighted" | 
             method == "Wald ratio") %>%
            filter(edge_dir == "Pcis_to_D")
PDtrans <- PD %>% 
            filter(method == "Inverse variance weighted" | 
                     method == "Wald ratio") %>%
            filter(edge_dir == "Ptrans_to_D")

ED <- ED %>% filter(method == "Inverse variance weighted" | 
                      method == "Wald ratio")
DE <- DE %>% filter(method == "Inverse variance weighted" | 
                      method == "Wald ratio")

PEcis <- PE %>% 
            filter(method == "Inverse variance weighted" | 
                     method == "Wald ratio") %>%
            filter(edge_dir == "Pcis_to_E")
PEtrans <- PE %>% 
            filter(method == "Inverse variance weighted" | 
                     method == "Wald ratio") %>%
            filter(edge_dir == "Ptrans_to_E")

DP <- DP %>% filter(method == "Inverse variance weighted" | 
                      method == "Wald ratio")


# Adjust p-values by test type:
# pick one: "BH" (FDR) or "bonferroni"
# pick one: "BH" (FDR) or "bonferroni"
ADJ_METHOD <- "BH"

EP      <- EP      %>% mutate(pval_adj = p.adjust(pval, method = ADJ_METHOD))
PDcis   <- PDcis   %>% mutate(pval_adj = p.adjust(pval, method = ADJ_METHOD))
PDtrans <- PDtrans %>% mutate(pval_adj = p.adjust(pval, method = ADJ_METHOD))
ED      <- ED      %>% mutate(pval_adj = p.adjust(pval, method = ADJ_METHOD))
DE      <- DE      %>% mutate(pval_adj = p.adjust(pval, method = ADJ_METHOD))
PEcis   <- PEcis   %>% mutate(pval_adj = p.adjust(pval, method = ADJ_METHOD))
PEtrans <- PEtrans %>% mutate(pval_adj = p.adjust(pval, method = ADJ_METHOD))
DP      <- DP      %>% mutate(pval_adj = p.adjust(pval, method = ADJ_METHOD))


#Plotting and Inference/Visualization

EP_sig <- EP %>% filter(pval_adj < 0.05)
PDcis_sig   <- PDcis   %>% filter(pval_adj < 0.05)
PDtrans_sig <- PDtrans %>% filter(pval_adj < 0.05)
ED_sig      <- ED      %>% filter(pval_adj < 0.05)
DE_sig      <- DE      %>% filter(pval_adj < 0.05)
PEcis_sig   <- PEcis   %>% filter(pval_adj < 0.05)
PEtrans_sig <- PEtrans %>% filter(pval_adj < 0.05)
DP_sig      <- DP      %>% filter(pval_adj < 0.05)


#EP)sig: just count stuff:



#'*Relating HEAP to MR hits*
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(stringr)
})

alpha <- 0.05

# ---- helper: find the effect-size column in each MR table ----
pick_effect_col <- function(df) {
  cand <- c("b", "beta", "beta.exposure", "b.exposure", "effect", "estimate")
  hit <- cand[cand %in% names(df)]
  if (length(hit) == 0) stop("No effect column found. Add your effect column name to pick_effect_col().")
  hit[1]
}

# ---- helper: find SE column in each MR table ----
pick_se_col <- function(df) {
  cand <- c("se", "se.exposure", "se.outcome", "sebeta", "bse", "stderr", "std_error", "StdError", "Std.Err", "SE")
  hit <- cand[cand %in% names(df)]
  if (length(hit) == 0) stop("No SE column found. Add your SE column name to pick_se_col().")
  hit[1]
}

# compute 95% CI from (beta, se)
add_ci <- function(beta, se, z = 1.96) {
  tibble::tibble(
    lo = beta - z * se,
    hi = beta + z * se
  )
}


# ---- helper: map (beta, padj) -> {+, −, 0, NA} ----
# edge_state <- function(beta, padj, alpha = 0.05) {
#   dplyr::case_when(
#     is.na(beta) | is.na(padj) ~ NA_character_,
#     padj < alpha & beta > 0   ~ "+",
#     padj < alpha & beta < 0   ~ "−",
#     TRUE                      ~ "0"
#   )
# }

edge_state <- function(beta, padj, alpha = 0.05) {
  dplyr::case_when(
    # if missing beta or padj, treat as "no evidence"
    is.na(beta) | is.na(padj) ~ "0",
    padj < alpha & beta > 0   ~ "+",
    padj < alpha & beta < 0   ~ "−",
    TRUE                      ~ "0"
  )
}


# ---- standardize each MR edge table into keys + beta/padj ----
std_EP <- function(df) {
  bcol <- pick_effect_col(df)
  scol <- pick_se_col(df)
  df %>%
    transmute(Exposure = id.exposure,
              Protein  = id.outcome,
              beta = .data[[bcol]],
              se   = .data[[scol]],
              padj = pval_adj)
}

std_PD <- function(df) {
  bcol <- pick_effect_col(df)
  scol <- pick_se_col(df)
  df %>%
    transmute(Protein = id.exposure,
              Disease = id.outcome,
              beta = .data[[bcol]],
              se   = .data[[scol]],
              padj = pval_adj)
}

std_ED <- function(df) {
  bcol <- pick_effect_col(df)
  scol <- pick_se_col(df)
  df %>%
    transmute(Exposure = id.exposure,
              Disease  = id.outcome,
              beta = .data[[bcol]],
              se   = .data[[scol]],
              padj = pval_adj)
}

std_PE <- function(df) {
  bcol <- pick_effect_col(df)
  scol <- pick_se_col(df)
  df %>%
    transmute(Protein  = id.exposure,
              Exposure = id.outcome,
              beta = .data[[bcol]],
              se   = .data[[scol]],
              padj = pval_adj)
}

std_DP <- function(df) {
  bcol <- pick_effect_col(df)
  scol <- pick_se_col(df)
  df %>%
    transmute(Disease = id.exposure,
              Protein = id.outcome,
              beta = .data[[bcol]],
              se   = .data[[scol]],
              padj = pval_adj)
}

std_DE <- function(df) {
  bcol <- pick_effect_col(df)
  scol <- pick_se_col(df)
  df %>%
    transmute(Disease  = id.exposure,
              Exposure = id.outcome,
              beta = .data[[bcol]],
              se   = .data[[scol]],
              padj = pval_adj)
}

# ---- master triplets from HEAP ----
trip <- HEAP %>%
  distinct(Exposure, Protein, Disease)

# ---- join in each link (left_join keeps HEAP as master; missing edges => NA) ----
trip2 <- trip %>%
  left_join(std_EP(EP) %>% rename(beta_EP = beta, se_EP = se, padj_EP = padj),
            by = c("Exposure","Protein")) %>%
  left_join(std_PD(PDcis) %>% rename(beta_PDcis = beta, se_PDcis = se, padj_PDcis = padj),
            by = c("Protein","Disease")) %>%
  left_join(std_PD(PDtrans) %>% rename(beta_PDtrans = beta, se_PDtrans = se, padj_PDtrans = padj),
            by = c("Protein","Disease")) %>%
  left_join(std_ED(ED) %>% rename(beta_ED = beta, se_ED = se, padj_ED = padj),
            by = c("Exposure","Disease")) %>%
  left_join(std_PE(PEcis) %>% rename(beta_PEcis = beta, se_PEcis = se, padj_PEcis = padj),
            by = c("Protein","Exposure")) %>%
  left_join(std_PE(PEtrans) %>% rename(beta_PEtrans = beta, se_PEtrans = se, padj_PEtrans = padj),
            by = c("Protein","Exposure")) %>%
  left_join(std_DP(DP) %>% rename(beta_DP = beta, se_DP = se, padj_DP = padj),
            by = c("Disease","Protein")) %>%
  left_join(std_DE(DE) %>% rename(beta_DE = beta, se_DE = se, padj_DE = padj),
            by = c("Disease","Exposure"))


# ---- call +/−/0/NA for each link ----
sig <- trip2 %>%
  mutate(
    state_EP      = edge_state(beta_EP,      padj_EP,      alpha),
    state_PDcis   = edge_state(beta_PDcis,   padj_PDcis,   alpha),
    state_PDtrans = edge_state(beta_PDtrans, padj_PDtrans, alpha),
    state_ED      = edge_state(beta_ED,      padj_ED,      alpha),
    state_PEcis   = edge_state(beta_PEcis,   padj_PEcis,   alpha),
    state_PEtrans = edge_state(beta_PEtrans, padj_PEtrans, alpha),
    state_DP      = edge_state(beta_DP,      padj_DP,      alpha),
    state_DE      = edge_state(beta_DE,      padj_DE,      alpha)
  )

link_order <- c("state_EP","state_PDcis","state_PDtrans","state_ED",
                "state_PEcis","state_PEtrans","state_DP","state_DE")

# build a compact signature string per triplet in a fixed column order
sig_ord <- sig %>%
  mutate(triplet = str_c(Exposure, Protein, Disease, sep = " | ")) %>%
  unite("signature", all_of(link_order), sep = "", remove = FALSE, na.rm = FALSE) %>%
  arrange(signature)   # simple ordering; see below for true clustering

row_levels <- sig_ord$triplet

#Think logically about the ordering:
table(sig_ord$signature)
unique(sig_ord$signature)
2^8


#'*MOTIF Matching!*

library(dplyr)
library(stringr)

is_present <- function(x) !is.na(x) & x != "0"

sig6 <- sig %>%
  mutate(
    # collapse cis/trans to "any evidence"
    state_PDany = case_when(
      state_PDcis == "+" | state_PDtrans == "+" ~ "+",
      state_PDcis == "−" | state_PDtrans == "−" ~ "−",
      TRUE ~ "0"
    ),
    state_PEany = case_when(
      state_PEcis == "+" | state_PEtrans == "+" ~ "+",
      state_PEcis == "−" | state_PEtrans == "−" ~ "−",
      TRUE ~ "0"
    ),
    
    # presence flags (✓ vs ✗ in your motif table)
    pres_EP  = is_present(state_EP),
    pres_PD  = is_present(state_PDany),
    pres_ED  = is_present(state_ED),
    pres_PE  = is_present(state_PEany),
    pres_DP  = is_present(state_DP),
    pres_DE  = is_present(state_DE),
    
    triplet = str_c(Exposure, Protein, Disease, sep=" | ")
  )

#Motif A — “Exposure influences Disease through Protein mediator”
#✓ on 1 (EP), ✓ on 2 (PD), ✓ on 3 (ED), ✗ on 4 (PE), ✗ on 5 (DP), ✗ on 6 (DE)
sig6 <- sig6 %>%
  mutate(
    motif_A_mediator =
      pres_EP & pres_PD & pres_ED & !pres_PE & !pres_DP & !pres_DE
  )

# Motif B — “Exposure influences a protein that is a disease biomarker”
# ✓ EP, ✗ PD, ✓ ED, ✗ PE, ✓ DP, ✗ DE
sig6 <- sig6 %>%
  mutate(
    motif_B_biomarker =
      pres_EP & !pres_PD & pres_ED & !pres_PE & pres_DP & !pres_DE
  )

# Motif C — “Exposure responsive marker (protein not causal w/ disease)”
# ✓ EP, ✗ PD, ✓ ED, ✗ PE, ✗ DP, ? DE (ignore DE)
sig6 <- sig6 %>%
  mutate(
    motif_C_exposure_marker =
      pres_EP & !pres_PD & pres_ED & !pres_PE & !pres_DP
    # DE can be anything
  )

# Motif D — “Protein influences exposure (biology drives behavior)”
# ✗ EP, ✓ PE, other edges are ? (ignore)
sig6 <- sig6 %>%
  mutate(
    motif_D_P_to_E =
      !pres_EP & pres_PE
  )

# Motif E — “Disease liability impacts both protein and exposure”
# ✓ DP and ✓ DE, others are ? (ignore)
sig6 <- sig6 %>%
  mutate(
    motif_E_disease_liability =
      pres_DP & pres_DE
  )

# if you still have your 8-char signature column:
# signature8 <- paste0(state_EP, state_PDcis, state_PDtrans, state_ED,
#                      state_PEcis, state_PEtrans, state_DP, state_DE)

sig6 <- sig6 %>%
  mutate(signature8 = paste0(state_EP, state_PDcis, state_PDtrans, state_ED,
                             state_PEcis, state_PEtrans, state_DP, state_DE))

count_sigs <- function(df) {
  df %>% count(signature8, sort = TRUE)
}

MotifList <- list(
  A = count_sigs(filter(sig6, motif_A_mediator)),
  B = count_sigs(filter(sig6, motif_B_biomarker)),
  C = count_sigs(filter(sig6, motif_C_exposure_marker)),
  D = count_sigs(filter(sig6, motif_D_P_to_E)),
  E = count_sigs(filter(sig6, motif_E_disease_liability))
)

#'*Plotting Time*
#Look at Triplets in each Motif:
MotifView <- data.frame(
  A = sum(MotifList[["A"]]$n),
  B = sum(MotifList[["B"]]$n),
  C = sum(MotifList[["C"]]$n),
  D = sum(MotifList[["D"]]$n),
  E = sum(MotifList[["E"]]$n)
)
head(MotifView)

# total triplets in your universe
N_total <- nrow(sig6)

MotifView <- tibble::tibble(
  motif = c("A","B","C","D","E"),
  n = c(
    sum(MotifList[["A"]]$n),
    sum(MotifList[["B"]]$n),
    sum(MotifList[["C"]]$n),
    sum(MotifList[["D"]]$n),
    sum(MotifList[["E"]]$n)
  )
) %>%
  mutate(frac = n / N_total)

# If motifs are mutually exclusive, this is fine:
MotifView <- MotifView %>%
  bind_rows(tibble::tibble(motif = "Other", n = N_total - sum(.$n), frac = (N_total - sum(.$n))/N_total))

MotifView


### alternative:
library(dplyr)

link_cols <- c("state_EP","state_PDcis","state_PDtrans","state_ED",
               "state_PEcis","state_PEtrans","state_DP","state_DE")

sig6_labeled <- sig6 %>%
  mutate(
    # any significant evidence anywhere in the 8 links?
    any_sig = rowSums(across(all_of(link_cols), ~ .x %in% c("+","−"))) > 0,
    
    # how many motifs does this triplet match?
    n_motifs = (motif_A_mediator + motif_B_biomarker + motif_C_exposure_marker +
                  motif_D_P_to_E + motif_E_disease_liability),
    
    # exclusive label with a sensible priority order
    motif_label = case_when(
      !any_sig ~ "Null (no MR hits)",
      
      # if you want to keep overlaps explicit:
      n_motifs > 1 ~ "Multiple motifs",
      
      motif_A_mediator          ~ "A",
      motif_B_biomarker         ~ "B",
      motif_C_exposure_marker   ~ "C",
      motif_D_P_to_E            ~ "D",
      motif_E_disease_liability ~ "E",
      
      TRUE ~ "Other (has hits)"
    )
  )

MotifView <- sig6_labeled %>%
  count(motif_label, name = "n") %>%
  mutate(frac = n / sum(n)) %>%
  arrange(desc(n))

MotifView

library(ggplot2)

ggplot(MotifView, aes(x = reorder(motif_label, n), y = n)) +
  geom_col() +
  coord_flip() +
  scale_y_log10() +
  theme_bw() +
  labs(x = NULL, y = "Triplet count (log10)")


positions = c("Null (no MR hits)", "Other (has hits)",
              "E","D","C","B","A","Multiple motifs")
  #c("Multiple motifs","A","B","C","D","E","Other (has hits)","Null (no MR hits)")
gg1 <- ggplot(MotifView, aes(x = reorder(motif_label, n), y = n)) +
  geom_col() +
  coord_flip() +
  scale_y_log10() +
  scale_x_discrete(limits = positions) +
  theme_bw() +
  labs(x = NULL, y = "Triad count (log10)") + theme(
    plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "cm")
  )
gg1

MotifView_lab <- MotifView %>%
  mutate(pct = n / sum(n),
         label_pct = paste0(scales::percent(pct, 0.01), " (n=", scales::comma(n), ")"))#percent(pct, accuracy = 0.01))

gg2 <- ggplot(MotifView_lab, aes(x = reorder(motif_label, n), y = n)) +
  geom_col() +
  geom_text(aes(label = label_pct),
            hjust = -0.1, size = 3) +
  coord_flip() +
  scale_y_log10() +
  scale_x_discrete(limits = positions) +
  theme_bw() +
  labs(x = NULL, y = "Triad count (log10)") +
  theme(plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "cm")) +
  expand_limits(y = max(MotifView_lab$n)^1.4)  # extra room for labels

gg2

ggsave("/n/groups/patel/shakson_ukb/UK_Biobank/Output/MRplots/motifcount.svg", 
       plot = gg2)
ggsave("/n/groups/patel/shakson_ukb/UK_Biobank/Output/MRplots/motifcount.png", 
       plot = gg2,
       dpi = 1000,
       width = 6.5, height = 4, units = "in")


fwrite(sig6_labeled, 
      file = "/n/groups/patel/shakson_ukb/UK_Biobank/Output/MR_edges/summary/MRmotifs.csv")




# Count features within each MR set & Determine prioritization material:

#'*Updated Version*
############################################################
# Panel B (compact): Hit-rate + Effect-size regime
# Fixes applied:
#   (1) Hit-rate: add Wilson 95% CI + n_sig/n_total labels + tidy axis
#   (2) Ridgeline: use |z| = |beta/se| among significant hits (comparable)
#       and clip to 99.5% to avoid extreme outliers destroying the axis
#
# Assumes you already created:
#   EP, PDcis, PDtrans, ED, DE, PEcis, PEtrans, DP
# and you already filtered method to IVW/Wald and created pval_adj in each.
############################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

HAS_GGRIDGES <- requireNamespace("ggridges", quietly = TRUE)

alpha <- 0.05

# ------------------------------------------------------------
# 0) Helpers
# ------------------------------------------------------------

pick_effect_col <- function(df) {
  cand <- c("b", "beta", "beta.exposure", "b.exposure", "effect", "estimate")
  hit <- cand[cand %in% names(df)]
  if (length(hit) == 0) stop("No effect column found. Add your effect column name to pick_effect_col().")
  hit[1]
}

pick_se_col <- function(df) {
  cand <- c("se", "se.exposure", "se.outcome", "se.b", "stderr", "SE")
  hit <- cand[cand %in% names(df)]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}


# Wilson 95% CI for binomial proportion
wilson_ci <- function(x, n, conf = 0.95) {
  if (is.na(n) || n == 0) return(c(NA_real_, NA_real_))
  z <- qnorm((1 + conf) / 2)
  p <- x / n
  den <- 1 + (z^2 / n)
  center <- (p + (z^2 / (2 * n))) / den
  half <- (z * sqrt((p * (1 - p) / n) + (z^2 / (4 * n^2)))) / den
  c(max(0, center - half), min(1, center + half))
}

# Standardize a table into edge-long format with beta/se/z
std_edge <- function(df, edge,
                     pair = c("E-P","E-D","P-D"),
                     direction = c("forward","reverse"),
                     instrument = NA_character_) {
  pair <- match.arg(pair)
  direction <- match.arg(direction)
  
  bcol  <- pick_effect_col(df)
  secol <- pick_se_col(df)
  
  df %>%
    transmute(
      edge = edge,
      pair = pair,
      direction = direction,
      instrument = instrument,
      method = if ("method" %in% names(df)) method else NA_character_,
      beta = .data[[bcol]],
      se   = if (!is.na(secol)) .data[[secol]] else NA_real_,
      q    = if ("pval_adj" %in% names(df)) pval_adj else NA_real_
    ) %>%
    mutate(
      abs_beta = abs(beta),
      z = ifelse(!is.na(se) & is.finite(se) & se > 0, beta / se, NA_real_),
      abs_z = abs(z),
      sig = !is.na(q) & (q < alpha)
    )
}

# ------------------------------------------------------------
# 1) Build unified "edges" table
# ------------------------------------------------------------

edges <- bind_rows(
  std_edge(EP,      edge="E→P",         pair="E-P", direction="forward", instrument="IV_E"),
  std_edge(PEcis,   edge="P→E (cis)",   pair="E-P", direction="reverse", instrument="IV_Pcis"),
  std_edge(PEtrans, edge="P→E (trans)", pair="E-P", direction="reverse", instrument="IV_Ptrans"),
  
  std_edge(ED,      edge="E→D",         pair="E-D", direction="forward", instrument="IV_E"),
  std_edge(DE,      edge="D→E",         pair="E-D", direction="reverse", instrument="IV_D"),
  
  std_edge(PDcis,   edge="P→D (cis)",   pair="P-D", direction="forward", instrument="IV_Pcis"),
  std_edge(PDtrans, edge="P→D (trans)", pair="P-D", direction="forward", instrument="IV_Ptrans"),
  std_edge(DP,      edge="D→P",         pair="P-D", direction="reverse", instrument="IV_D")
) %>%
  mutate(
    edge = factor(edge, levels = c(
      "E→P", "P→E (cis)", "P→E (trans)",
      "E→D", "D→E",
      "P→D (cis)", "P→D (trans)", "D→P"
    )),
    pair = factor(pair, levels = c("E-P","E-D","P-D")),
    direction = factor(direction, levels = c("forward","reverse"))
  )


# If you already have hit_rate2, skip this block
hit_rate2 <- edges %>%
  group_by(edge, pair, direction) %>%
  summarise(
    n_total = sum(!is.na(q)),
    n_sig   = sum(sig, na.rm = TRUE),
    prop    = ifelse(n_total > 0, n_sig / n_total, NA_real_),
    .groups = "drop"
  ) %>%
  mutate(
    pct = 100 * prop,
    lab = paste0(n_sig, "/", n_total)
  )

# -----------------------------
# Option A ordering: forward first, reverse second
# - E-P: E→P (forward), then P→E (cis), P→E (trans)
# - E-D: E→D (forward), then D→E (reverse)
# - P-D: P→D (cis), P→D (trans) (forward), then D→P (reverse)
# -----------------------------
edge_levels_optA <- c(
  "E→P",
  "P→E (cis)", "P→E (trans)",
  "E→D",
  "D→E",
  "P→D (cis)", "P→D (trans)",
  "D→P"
)

# For coord_flip: first level appears at bottom, so reverse levels for display
edge_levels_display <- rev(edge_levels_optA)

hit_rate2 <- hit_rate2 %>%
  mutate(
    edge = factor(as.character(edge), levels = edge_levels_display),
    pair = factor(as.character(pair), levels = c("E-P","E-D","P-D"))
  )

ymax <- max(hit_rate2$pct, na.rm = TRUE)
ymax_pad <- ymax * 1.30

p_hit_final <- ggplot(hit_rate2, aes(x = edge, y = pct)) +
  geom_col(width = 0.8) +
  geom_text(aes(label = lab), hjust = -0.1, size = 3) +
  facet_grid(pair ~ ., scales = "free_y", space = "free_y", drop = TRUE) +  # drop unused levels
  scale_y_continuous(
    limits = c(0, ymax_pad),
    breaks = seq(0, ceiling(ymax_pad/5)*5, by = 5)
  ) +
  coord_flip(clip = "off") +
  theme_bw() +
  labs(
    x = NULL,
    y = paste0("% significant (q < ", alpha, ")"),
    title = "MR hit-rate by edge type"
  ) +
  theme(
    strip.background = element_rect(fill = NA),
    plot.margin = margin(t = 10, r = 45, b = 10, l = 10)
  )

p_hit_final


ggsave("/n/groups/patel/shakson_ukb/UK_Biobank/Output/MRplots/MRperchits.svg", 
       plot = p_hit_final)
ggsave("/n/groups/patel/shakson_ukb/UK_Biobank/Output/MRplots/MRperchits.png", 
       plot = p_hit_final,
       dpi = 1000, width = 5, height = 4, units = "in")


# -----------------------------
# B) Forest summary of |z| among significant hits
# Use median + IQR 
# -----------------------------

edge_levels_optA <- c(
  "E→P",
  "P→E (cis)", "P→E (trans)",
  "E→D",
  "D→E",
  "P→D (cis)", "P→D (trans)",
  "D→P"
)

# For y-axis: first level appears at bottom, so reverse for display
edge_levels_display <- rev(edge_levels_optA)

sig_edges <- edges %>%
  filter(sig, !is.na(abs_z), is.finite(abs_z))

forest_z <- sig_edges %>%
  group_by(edge, pair) %>%
  summarise(
    n = n(),
    med = median(abs_z, na.rm = TRUE),
    q25 = quantile(abs_z, 0.25, na.rm = TRUE),
    q75 = quantile(abs_z, 0.75, na.rm = TRUE),
    q10 = quantile(abs_z, 0.10, na.rm = TRUE),
    q90 = quantile(abs_z, 0.90, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    lab = paste0("n=", n),
    pair = factor(as.character(pair), levels = c("E-P","E-D","P-D")),
    edge = factor(as.character(edge), levels = edge_levels_display)
  )

p_forest <- ggplot(forest_z, aes(y = edge, x = med)) +
  geom_segment(aes(x = q10, xend = q90, yend = edge), linewidth = 0.5) +
  geom_segment(aes(x = q25, xend = q75, yend = edge), linewidth = 2) +
  geom_point(size = 2) +
  facet_grid(pair ~ ., scales = "free_y", space = "free_y", drop = TRUE) +
  theme_bw() +
  labs(
    x = "|z| (median; thick=IQR, thin=10–90%)",
    y = NULL,
    title = "Strength of MR evidence"
  ) +
  theme(
    strip.background = element_rect(fill = NA),
    plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "cm")
  )

p_forest

ggsave("/n/groups/patel/shakson_ukb/UK_Biobank/Output/MRplots/MReffsize.svg", 
       plot = p_forest)
ggsave("/n/groups/patel/shakson_ukb/UK_Biobank/Output/MRplots/MReffsize.png", 
       plot = p_forest,
       dpi = 1000,
       width = 4, height = 4, units = "in")


#### Pcis and Ptrans effect ####

PDeffs <- sig6_labeled %>% select(c("Protein","Disease",
                          "beta_PDcis","padj_PDcis",
                          "beta_PDtrans","padj_PDtrans")) %>%
                          unique()
cor.test(PDeffs$beta_PDcis, PDeffs$beta_PDtrans)


PDeffs_sig <- sig6_labeled %>% 
                    select(c("Protein","Disease",
                           "beta_PDcis","padj_PDcis",
                           "beta_PDtrans","padj_PDtrans")) %>%
                    unique() %>%
                    filter(padj_PDcis < 0.05 & padj_PDtrans < 0.05)

cor.test(PDeffs_sig$beta_PDcis, PDeffs_sig$beta_PDtrans)

#### MOTIF DAGs ####

#'*Look through Motifs*
#Pick Relevant DAGs for Motifs
MotifA <- sig6_labeled %>% filter(motif_label == "A")
MotifB <- sig6_labeled %>% filter(motif_label == "B")
MotifC <- sig6_labeled %>% filter(motif_label == "C")
MotifD <- sig6_labeled %>% filter(motif_label == "D")
MotifE <- sig6_labeled %>% filter(motif_label == "E")


head(MotifA)


# A is missing E->P, P->D but no E->D answer 



#### PLOT a SPECIFIC TRIAD ####

# ============================================================
# Triplet MR diagram — "Nature-ish" styling
# Fixes:
#  1) user-defined title (short)
#  2) edge labels in white boxes + drawn on top
#  3) DP curve bows outward (same side as PE/DE)
#  4) bigger nodes + less wasted whitespace
#  5) cleaner aesthetics: subtle non-sig edges, compact legend
#  - Bigger nodes AND bigger node text (so labels fit)
#  - Legend "Significant" is now LINESTYLE ONLY (no arrows) to avoid confusion
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggforce)
  library(grid)
  library(tibble)
})

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
    
    # labeling
    label_mode = c("sig_only","all","none"),
    label_box = TRUE,
    label_pad = 0.20,
    
    # geometry (BIGGER DEFAULTS)
    node_radius = 0.95,        # arrow shortening
    node_size = 32,            # bigger circles
    node_text_size = 5,        # bigger node label text
    arrow_mm = 3.0,
    edge_lwd = 1.5,
    
    # layout
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
    
    # forward straight; reverse curved outward
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
      
      # shorten so arrowheads don't enter nodes
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
    # IMPORTANT: make the "Significant" legend use lines, not arrows
    guides(
      color = guide_legend(title="Edge type", nrow=2, override.aes=list(alpha=1, linetype="solid")),
      linetype = guide_legend(
        title = paste0("Significant (q<", alpha, ")"),
        override.aes = list(color="black", alpha=1)  # no arrows here
      ),
      alpha = "none"
    ) +
    labs(title = title %||% triplet_id, subtitle = subtitle)
  
  g
}

#---- Examples ----
p <- plot_triplet_mr_diagram_nature(
  DT = sig6_labeled,
  triplet_id = "time_spent_watching_television_tv_f1070_0_0 | FURIN | finngen_R12_I9_HYPTENSESS",
  title = "TV time \u2192 FURIN \u2192 Hypertension",
  alpha = 0.05,
  show = "all",
  pd_mode = "cis",
  pe_mode = "cis",
  label_mode = "sig_only"
)
print(p)


p <-  plot_triplet_mr_diagram_nature(
  DT = sig6_labeled,
  triplet_id = "time_spent_watching_television_tv_f1070_0_0 | LEP | finngen_R12_E4_OBESITYCAL",
  title = "TV time \u2192 LEP  & Obesity \u2192 LEP ",
  alpha = 0.05,
  show = "all",
  pd_mode = "cis",
  pe_mode = "cis",
  digits = 3,
  label_mode = "sig_only"  # recommended
)
print(p)


p <- plot_triplet_mr_diagram_nature(
  DT = sig6_labeled,
  triplet_id = "current_tobacco_smoking_f1239_0_0_Yes._on_most_or_all_days | CXCL17 | finngen_R12_J10_EMPHYSEMA",
  title = "Every Day Smoker \u2192 CXCL17  & Emphysema \u2192 CXCL17 ",
  motif = "B",
  alpha = 0.05,
  show = "all",
  pd_mode = "cis",
  pe_mode = "cis",
  digits = 3,
  label_mode = "sig_only"  # recommended
)
print(p)

p <- plot_triplet_mr_diagram_nature(
  DT = sig6_labeled,
  triplet_id = "usual_walking_pace_f924_0_0 | IL1RN | finngen_R12_I9_HEARTFAIL",
  title = "Walking Pace \u2192 IL1RN  & Walking Pace \u2192 HeartFailure",
  motif = "C",
  alpha = 0.05,
  show = "all",
  pd_mode = "cis",
  pe_mode = "cis",
  digits = 3,
  label_mode = "sig_only"  # recommended
)
print(p)


p <- plot_triplet_mr_diagram_nature(
  DT = sig6_labeled,
  triplet_id = "cereal_intake_f1458_0_0 | APOC1 | finngen_R12_I9_IHD",
  title = "APOC1 \u2192 cereal intake  & Ischemic Heart Disease \u2192 APOC1",
  motif = "D",
  alpha = 0.05,
  show = "all",
  pd_mode = "cis",
  pe_mode = "cis",
  digits = 3,
  label_mode = "sig_only"  # recommended
)
print(p)


p <- plot_triplet_mr_diagram_nature(
  DT = sig6_labeled,
  triplet_id = "cereal_intake_f1458_0_0 | APOC1 | finngen_R12_T2D",
  title = "T2D \u2192 APOC1 & APOC1 \u2192 cereal intake",
  motif = "D",
  alpha = 0.05,
  show = "all",
  pd_mode = "cis",
  pe_mode = "cis",
  digits = 3,
  label_mode = "sig_only"  # recommended
)
print(p)


#PRSS8 is some sodium transporter
# Hypertension --> influence dont take as much salt. 
# Increase PRSS8 
# high levels of prostasin promote sodium retention
p <- plot_triplet_mr_diagram_nature(
  DT = sig6_labeled,
  triplet_id = "salt_added_to_food_f1478_0_0 | PRSS8 | finngen_R12_I9_HYPTENSESS",
  title = "Hypertension \u2192 Salt Intake &  Hypertension \u2192 PRSS8",
  motif = "E",
  alpha = 0.05,
  show = "all",
  pd_mode = "cis",
  pe_mode = "cis",
  digits = 3,
  label_mode = "sig_only"  # recommended
)
print(p)


p <- plot_triplet_mr_diagram_nature(
  DT = sig6_labeled,
  triplet_id = "major_dietary_changes_in_the_last_5_years_f1538_0_0_Yes._because_of_illness | GDF15 | finngen_R12_T2D",
  title = "T2D \u2192 Prior Dietary Change (within 5 years) &  T2D \u2192 GDF15",
  motif = "E",
  alpha = 0.05,
  show = "all",
  pd_mode = "cis",
  pe_mode = "cis",
  digits = 3,
  label_mode = "sig_only"  # recommended
)
print(p)


p <- plot_triplet_mr_diagram_nature(
  DT = sig6_labeled,
  triplet_id = "major_dietary_changes_in_the_last_5_years_f1538_0_0_Yes._because_of_illness | GUSB | finngen_R12_T2D_WIDE",
  title = "T2D \u2192 Prior Dietary Change (within 5 years) &  T2D \u2192 GUSB",
  motif = "E",
  alpha = 0.05,
  show = "all",
  pd_mode = "cis",
  pe_mode = "cis",
  digits = 3,
  label_mode = "sig_only"  # recommended
)
print(p)



















#### OLD VERSION OF TRYING TO PLOT TRIAD !!!! ####


# ============================================================
# Triplet MR diagram — forward edges straight, reverse edges curved
#  Forward:  EP, PDcis/PDtrans, ED
#  Reverse:  PEcis/PEtrans, DP, DE
#  - No overlap / bidirectional confusion
#  - Single bezier/segment layer per edge type (no double painting)
#  - Non-sig edges dashed + faded; sig edges solid + opaque
#  - Labels: β [95% CI] + stars
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggforce)
  library(grid)
  library(tibble)
})

`%||%` <- function(a,b) if (!is.null(a)) a else b

pstars <- function(p) {
  ifelse(is.na(p), "",
         ifelse(p < 0.001, "***",
                ifelse(p < 0.01, "**",
                       ifelse(p < 0.05, "*", ""))))
}

plot_triplet_mr_diagram_straight_forward <- function(DT,
                                                     triplet_id,
                                                     motif = NULL,
                                                     alpha = 0.05,
                                                     show = c("all","sig_only"),
                                                     pd_mode = c("cis","trans","both"),
                                                     pe_mode = c("cis","trans","both"),
                                                     digits = 3,
                                                     title = NULL,
                                                     label_mode = c("sig_only","all","none"),
                                                     node_radius = 0.55,
                                                     arrow_mm = 2.6,
                                                     edge_lwd = 1.25) {
  
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
  
  # Edges: forward straight (curve_mag=0), reverse curved (curve_mag>0)
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
    
    # forward edges straight; reverse edges curved
    curve_mag = c(0.00, 0.34, 0.46,
                  0.00, 0.00, 0.28,
                  0.00, 0.24),
    
    # only matters when curve_mag>0
    # choose a consistent bend "side" for reverse edges
    side      = c(0, -1, -1,
                  0,  0, +1,
                  0, +1),
    
    # label nudges (forward baseline label nudges down slightly)
    ly_nudge  = c(+0.20, -0.25, -0.35,
                  +0.22, +0.22, +0.28,
                  -0.22, +0.25),
    
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
  
  # drop cis/trans edges if requested
  if (pd_mode != "both") {
    edges_def <- edges_def %>%
      filter(!(edge_type %in% c("PDcis","PDtrans")) | edge_type == paste0("PD", pd_mode))
  }
  if (pe_mode != "both") {
    edges_def <- edges_def %>%
      filter(!(edge_type %in% c("PEcis","PEtrans")) | edge_type == paste0("PE", pe_mode))
  }
  
  # edge table for this row
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
  
  # join coords + shorten to avoid arrowheads inside nodes
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
  
  # labels selection
  edges_lab <- edges
  if (label_mode == "none") edges_lab <- edges_lab %>% filter(FALSE)
  if (label_mode == "sig_only") edges_lab <- edges_lab %>% filter(sig)
  
  # split forward (straight) and reverse (curved)
  edges_straight <- edges %>% filter(curve_mag == 0)
  edges_curved   <- edges %>% filter(curve_mag > 0)
  
  # build bezier points for curved edges only
  if (nrow(edges_curved) > 0) {
    edges_curved <- edges_curved %>%
      mutate(
        # perpendicular unit vector for curved control point
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
  
  ggplot() +
    # ---- straight edges (forward) ----
  geom_segment(
    data = edges_straight,
    aes(x=x_from2, y=y_from2, xend=x_to2, yend=y_to2, color=edge_type, linetype=sig, alpha=sig),
    arrow = arrow(type="closed", length = unit(arrow_mm, "mm")),
    linewidth = edge_lwd
  ) +
    # ---- curved edges (reverse) ----
  ggforce::geom_bezier(
    data = bez,
    aes(x=x, y=y, group=edge_type, color=edge_type, linetype=sig, alpha=sig),
    arrow = arrow(type="closed", length = unit(arrow_mm, "mm")),
    linewidth = edge_lwd
  ) +
    scale_linetype_manual(
      values = c(`TRUE`="solid", `FALSE`="dashed"),
      name = paste0("Significant (q<", alpha, ")")
    ) +
    scale_alpha_manual(
      values = c(`TRUE`=1, `FALSE`=0.22),
      name = paste0("Significant (q<", alpha, ")")
    ) +
    geom_text(
      data = edges_lab,
      aes(x=x_mid, y=y_mid + ly_nudge, label=label),
      size = 4
    ) +
    geom_point(data=nodes, aes(x=x, y=y), size=18, shape=21, stroke=1.2, fill="white") +
    geom_text(data=nodes, aes(x=x, y=y, label=node), size=5) +
    coord_equal(xlim=c(-1, 11), ylim=c(-1.2, 4.9), clip="off") +
    theme_void() +
    theme(
      plot.title = element_text(hjust=0.5, size=14, face="bold"),
      legend.position = "bottom",
      legend.box = "vertical",
      legend.title = element_text(size=11),
      legend.text  = element_text(size=10)
    ) +
    guides(
      color = guide_legend(title="Edge type", nrow=2, override.aes=list(alpha=1, linetype="solid")),
      linetype = guide_legend(override.aes=list(color="black")),
      alpha = "none"
    ) +
    labs(title = title %||% triplet_id)
}


suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggforce)
  library(grid)
  library(tibble)
})

pstars <- function(p) {
  ifelse(is.na(p), "",
         ifelse(p < 0.001, "***",
                ifelse(p < 0.01, "**",
                       ifelse(p < 0.05, "*", ""))))
}

`%||%` <- function(a,b) if (!is.null(a)) a else b

plot_triplet_mr_diagram2 <- function(DT,
                                     triplet_id,
                                     motif = NULL,
                                     alpha = 0.05,
                                     show = c("all","sig_only"),
                                     pd_mode = c("cis","trans","both"),
                                     pe_mode = c("cis","trans","both"),
                                     digits = 3,
                                     title = NULL,
                                     label_mode = c("sig_only","all","none"),
                                     node_radius = 0.55,   # controls arrow/node separation in data units
                                     arrow_mm = 2.6) {
  
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
    curve     = c(+0.28, -0.28, -0.40,
                  -0.22, -0.34, +0.22,
                  0.00, +0.18),
    ly_nudge  = c(+0.25, -0.25, -0.35,
                  +0.28, +0.15, +0.28,
                  -0.30, +0.30),
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
  
  if (pd_mode != "both") edges_def <- edges_def %>% filter(!(edge_type %in% c("PDcis","PDtrans")) | edge_type == paste0("PD", pd_mode))
  if (pe_mode != "both") edges_def <- edges_def %>% filter(!(edge_type %in% c("PEcis","PEtrans")) | edge_type == paste0("PE", pe_mode))
  
  edges <- edges_def %>%
    rowwise() %>%
    mutate(beta = row[[beta_col]], se = row[[se_col]], padj = row[[p_col]]) %>%
    ungroup() %>%
    mutate(
      sig = !is.na(padj) & padj < alpha,
      lo = ifelse(is.na(beta) | is.na(se), NA_real_, beta - 1.96*se),
      hi = ifelse(is.na(beta) | is.na(se), NA_real_, beta + 1.96*se),
      label = case_when(
        is.na(beta) ~ "",
        is.na(se)   ~ paste0("β=", formatC(beta, format="f", digits=digits), pstars(padj)),
        TRUE        ~ paste0("β=", formatC(beta, format="f", digits=digits),
                             " [", formatC(lo, format="f", digits=digits), ", ",
                             formatC(hi, format="f", digits=digits), "]",
                             pstars(padj))
      )
    )
  
  if (show == "sig_only") edges <- edges %>% filter(sig)
  
  # join coords
  edges <- edges %>%
    left_join(nodes %>% rename(from=node, x_from=x, y_from=y), by="from") %>%
    left_join(nodes %>% rename(to=node,   x_to=x,   y_to=y), by="to") %>%
    mutate(
      dx = x_to - x_from,
      dy = y_to - y_from,
      L  = sqrt(dx^2 + dy^2),
      # shorten line so arrowheads don't enter the node circles
      x_from2 = x_from + node_radius * dx / L,
      y_from2 = y_from + node_radius * dy / L,
      x_to2   = x_to   - node_radius * dx / L,
      y_to2   = y_to   - node_radius * dy / L,
      x_mid = (x_from2 + x_to2)/2,
      y_mid = (y_from2 + y_to2)/2,
      # perpendicular unit vector for curve
      ux = -dy / L,
      uy =  dx / L,
      bend = curve * L,
      x_ctrl = x_mid + bend * ux,
      y_ctrl = y_mid + bend * uy
    )
  
  bez <- bind_rows(
    edges %>% transmute(edge_type, sig, t=1, x=x_from2, y=y_from2),
    edges %>% transmute(edge_type, sig, t=2, x=x_ctrl,  y=y_ctrl),
    edges %>% transmute(edge_type, sig, t=3, x=x_to2,   y=y_to2)
  )
  
  # label selection
  edges_lab <- edges
  if (label_mode == "none") edges_lab <- edges_lab %>% filter(FALSE)
  if (label_mode == "sig_only") edges_lab <- edges_lab %>% filter(sig)
  
  ggplot() +
    ggforce::geom_bezier(
      data = bez,
      aes(x=x, y=y, group=edge_type, color=edge_type, linetype=sig, alpha=sig),
      arrow = arrow(type="closed", length = unit(arrow_mm, "mm")),
      linewidth = 1.2
    ) +
    scale_linetype_manual(values = c(`TRUE`="solid", `FALSE`="dashed"),
                          name = paste0("Significant (q<", alpha, ")")) +
    scale_alpha_manual(values = c(`TRUE`=1, `FALSE`=0.25),
                       name = paste0("Significant (q<", alpha, ")")) +
    geom_text(
      data = edges_lab,
      aes(x=x_mid, y=y_mid + ly_nudge, label=label),
      size = 4
    ) +
    geom_point(data=nodes, aes(x=x, y=y), size=18, shape=21, stroke=1.2, fill="white") +
    geom_text(data=nodes, aes(x=x, y=y, label=node), size=5) +
    coord_equal(xlim=c(-1, 11), ylim=c(-1.2, 4.9), clip="off") +
    theme_void() +
    theme(
      plot.title = element_text(hjust=0.5, size=14, face="bold"),
      legend.position = "bottom",
      legend.box = "vertical",
      legend.title = element_text(size=11),
      legend.text  = element_text(size=10)
    ) +
    guides(
      color = guide_legend(title="Edge type", nrow=2, override.aes=list(alpha=1)),
      linetype = guide_legend(override.aes=list(color="black")),
      alpha = "none"
    ) +
    labs(title = title %||% triplet_id)
}

p <- plot_triplet_mr_diagram2(
  DT = sig6_labeled,
  triplet_id = "time_spent_watching_television_tv_f1070_0_0 | LEP | finngen_R12_E4_OBESITYCAL",
  motif = "B",
  alpha = 0.05,
  show = "all",
  pd_mode = "cis",
  pe_mode = "cis",
  digits = 3,
  label_mode = "sig_only"   # <- huge readability win
)
print(p)


suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggforce)
  library(grid)
})

pstars <- function(p) {
  ifelse(is.na(p), "",
         ifelse(p < 0.001, "***",
                ifelse(p < 0.01, "**",
                       ifelse(p < 0.05, "*", ""))))
}

plot_triplet_mr_diagram <- function(DT,
                                    triplet_id,
                                    motif = NULL,
                                    alpha = 0.05,
                                    show = c("all","sig_only"),
                                    pd_mode = c("cis","trans","both"),
                                    pe_mode = c("cis","trans","both"),
                                    digits = 3,
                                    title = NULL) {
  
  show <- match.arg(show)
  pd_mode <- match.arg(pd_mode)
  pe_mode <- match.arg(pe_mode)
  
  row <- DT %>% as.data.frame() %>% filter(triplet == triplet_id)
  if (!is.null(motif)) row <- row %>% filter(motif_label == motif)
  if (nrow(row) != 1) stop("Expected exactly 1 row after filtering; got n=", nrow(row))
  
  nodes <- tibble::tibble(
    node = c("Exposure","Protein","Disease"),
    x    = c(0, 5, 10),
    y    = c(0, 3, 0)
  )
  
  # 8-edge version (EP, ED, DP, DE, and cis/trans for PD & PE)
  edges_def <- tibble::tibble(
    edge_type = c("EP","PEcis","PEtrans",
                  "PDcis","PDtrans","DP",
                  "ED","DE"),
    from      = c("Exposure","Protein","Protein",
                  "Protein","Protein","Disease",
                  "Exposure","Disease"),
    to        = c("Protein","Exposure","Exposure",
                  "Disease","Disease","Protein",
                  "Disease","Exposure"),
    
    # Curves are the key: opposite signs for reverse edges, and different magnitudes for cis/trans
    curve     = c( +0.28, -0.28, -0.40,     # E<->P lanes (cis vs trans further apart)
                   -0.22, -0.34, +0.22,      # P<->D lanes (cis vs trans + reverse)
                   0.00, +0.18),            # E->D straight, D->E small arc (so they don't overlap)
    
    # nudge labels away from nodes / baseline; tune if needed
    ly_nudge  = c( +0.25, -0.25, -0.35,
                   +0.28, +0.15, +0.28,
                   -0.28, +0.28),
    
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
  
  
  # optionally drop PD/PE cis/trans
  if (pd_mode != "both") edges_def <- edges_def %>% filter(!(edge_type %in% c("PDcis","PDtrans")) | edge_type == paste0("PD", pd_mode))
  if (pe_mode != "both") edges_def <- edges_def %>% filter(!(edge_type %in% c("PEcis","PEtrans")) | edge_type == paste0("PE", pe_mode))
  
  # build edge table for this row
  edges <- edges_def %>%
    rowwise() %>%
    mutate(
      beta = row[[beta_col]],
      se   = row[[se_col]],
      padj = row[[p_col]]
    ) %>%
    ungroup() %>%
    mutate(
      sig = !is.na(padj) & padj < alpha,
      lo = ifelse(is.na(beta) | is.na(se), NA_real_, beta - 1.96*se),
      hi = ifelse(is.na(beta) | is.na(se), NA_real_, beta + 1.96*se),
      label = dplyr::case_when(
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
  
  # join coords
  edges <- edges %>%
    left_join(nodes %>% rename(from=node, x_from=x, y_from=y), by="from") %>%
    left_join(nodes %>% rename(to=node,   x_to=x,   y_to=y), by="to") %>%
    mutate(
      x_mid = (x_from + x_to)/2,
      y_mid = (y_from + y_to)/2
    ) %>%
    mutate(
      dx = x_to - x_from,
      dy = y_to - y_from,
      L  = sqrt(dx^2 + dy^2),
      ux = ifelse(L == 0, 0, -dy / L),
      uy = ifelse(L == 0, 0,  dx / L),
      bend = curve * L,
      x_ctrl = x_mid + bend * ux,
      y_ctrl = y_mid + bend * uy
    )
  
  bez <- bind_rows(
    edges %>% transmute(edge_type, t=1, x=x_from, y=y_from),
    edges %>% transmute(edge_type, t=2, x=x_ctrl, y=y_ctrl),
    edges %>% transmute(edge_type, t=3, x=x_to,   y=y_to)
  )
  
  ggplot() +
    ggforce::geom_bezier(
      data = bez,
      aes(x=x, y=y, group=edge_type, color=edge_type, linetype=edge_type),
      arrow = arrow(type="closed", length = unit(2.8, "mm")),
      linewidth = 1.2
    ) +
    # linetype mapped to significance by overriding in a second scale:
    # easiest is to set it directly:
    ggforce::geom_bezier(
      data = bez %>% left_join(edges %>% select(edge_type, sig), by="edge_type"),
      aes(x=x, y=y, group=edge_type, color=edge_type, linetype=sig),
      arrow = arrow(type="closed", length = unit(2.8, "mm")),
      linewidth = 1.2,
      show.legend = TRUE
    ) +
    scale_linetype_manual(values = c(`TRUE`="solid", `FALSE`="dashed"), name = paste0("Significant (q<", alpha, ")")) +
    geom_text(
      data = edges,
      aes(x=x_mid, y=y_mid + ly_nudge, label=label),
      size = 4
    ) +
    geom_point(data=nodes, aes(x=x, y=y), size=18, shape=21, stroke=1.2) +
    geom_text(data=nodes, aes(x=x, y=y, label=node), size=5) +
    coord_equal(xlim=c(-1, 11), ylim=c(-1, 4.8), clip="off") +
    theme_void() +
    theme(
      plot.title = element_text(hjust=0.5, size=14, face="bold"),
      legend.position = "bottom",
      legend.box = "vertical"
    ) +
    guides(color = guide_legend(title="Edge type", nrow=2),
           linetype = guide_legend(title=paste0("Significant (q<", alpha, ")"))) +
    labs(title = title %||% triplet_id)
}

`%||%` <- function(a,b) if (!is.null(a)) a else b

p <- plot_triplet_mr_diagram(
  DT = sig6_labeled,
  triplet_id = "time_spent_watching_television_tv_f1070_0_0 | LEP | finngen_R12_E4_OBESITYCAL",
  motif = "B",
  alpha = 0.05,
  show = "all",     # "sig_only" if you want to drop non-sig edges
  pd_mode = "cis",  # or "trans" or "both"
  pe_mode = "cis",  # or "trans" or "both"
  digits = 3
)
print(p)



#sig6_labeled$triplet[1]
#sig6_labeled %>% filter(triplet == "time_spent_watching_television_tv_f1070_0_0 | LEP | finngen_R12_E4_OBESITYCAL")





#### TRASH #####


hit_rate <- edges %>%
  group_by(edge, pair, direction) %>%
  summarise(
    n_total = sum(!is.na(q)),
    n_sig   = sum(sig, na.rm = TRUE),
    prop    = ifelse(n_total > 0, n_sig / n_total, NA_real_),
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(
    ci = list(wilson_ci(n_sig, n_total, conf = 0.95)),
    lo = ci[[1]][1],
    hi = ci[[1]][2]
  ) %>%
  ungroup() %>%
  mutate(
    pct = 100 * prop,
    lo_pct = 100 * lo,
    hi_pct = 100 * hi,
    lab = paste0(n_sig, "/", n_total)
  )

ymax <- max(hit_rate$pct, na.rm = TRUE) * 1.15

p_hit <- ggplot(hit_rate, aes(x = edge, y = pct)) +
  geom_col(width = 0.8) +
  geom_errorbar(aes(ymin = lo_pct, ymax = hi_pct), width = 0.2) +
  geom_text(aes(label = lab), hjust = -0.1, size = 3) +
  coord_flip() +
  facet_grid(pair ~ ., scales = "free_y", space = "free_y") +
  scale_y_continuous(
    breaks = seq(0, ceiling(ymax/5)*5, by = 5),
    limits = c(0, ymax)
  ) +
  theme_bw() +
  labs(
    x = NULL,
    y = paste0("% significant (q < ", alpha, ")"),
    title = "MR hit-rate by edge type (Wilson 95% CI)"
  ) +
  theme(
    strip.background = element_rect(fill = NA),
    plot.margin = margin(10, 10, 10, 10)
  )

p_hit


sig_edges <- edges %>%
  filter(sig, !is.na(abs_z), is.finite(abs_z))

if (nrow(sig_edges) == 0) {
  warning("No significant hits (q < alpha). Ridgeline plot will be empty.")
}

# Clip extreme |z| so the distribution is readable
CLIP_Q <- 0.995
xmax <- as.numeric(quantile(sig_edges$abs_z, CLIP_Q, na.rm = TRUE))
if (!is.finite(xmax) || is.na(xmax)) xmax <- max(sig_edges$abs_z, na.rm = TRUE)

sig_edges <- sig_edges %>%
  mutate(abs_z_clip = pmin(abs_z, xmax))

xlab <- paste0("|z| among significant hits (clipped at ", CLIP_Q*100, "% = ", signif(xmax, 3), ")")

if (HAS_GGRIDGES) {
  p_ridge <- ggplot(sig_edges, aes(x = abs_z_clip, y = edge)) +
    ggridges::geom_density_ridges(scale = 1.2, rel_min_height = 0.01, alpha = 0.85) +
    facet_grid(pair ~ ., scales = "free_y", space = "free_y") +
    theme_bw() +
    labs(
      x = xlab,
      y = NULL,
      title = "Effect-size regime among significant MR hits"
    ) +
    theme(
      strip.background = element_rect(fill = NA),
      plot.margin = margin(10, 10, 10, 10)
    )
  p_ridge
} else {
  message("Package 'ggridges' not installed. Using violin fallback.")
  p_ridge <- ggplot(sig_edges, aes(x = abs_z_clip, y = edge)) +
    geom_violin(trim = TRUE) +
    facet_grid(pair ~ ., scales = "free_y", space = "free_y") +
    theme_bw() +
    labs(
      x = xlab,
      y = NULL,
      title = "Effect-size regime among significant MR hits (violin fallback)"
    )
  p_ridge
}

###'*REDO Plot:*
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})








