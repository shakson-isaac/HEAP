# ============================================================
# Motif-type overlap / uniqueness using triplet master MRmotifs.csv
# (Uses Exposure/Protein/Disease + motif_* columns already computed)
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

MRfiles <- "/n/groups/patel/shakson_ukb/UK_Biobank/Output/MR_edges"
output_dir <- "/n/groups/patel/shakson_ukb/UK_Biobank/Output/MRplots/"
OUTDIR <- file.path(output_dir , "COMPARE_UKB_vs_DECODE")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

# Triplet master paths (UPDATED)
TRIPLET_UKB_FP <- file.path(MRfiles, "summary", "MRmotifs.csv")
TRIPLET_DEC_FP <- file.path(MRfiles, "summary", "DECODE", "MRmotifs.csv")

# Toggle this depending on whether you want to include "None"/non-sig triplets
KEEP_ONLY_ANY_SIG <- FALSE   # TRUE = only triplets with any_sig == TRUE

# ----------------------------
# Read + validate
# ----------------------------
read_trip <- function(fp, tag) {
  df <- fread(fp) %>% as_tibble()
  needed <- c("Exposure","Protein","Disease","triplet","motif_label",
              "motif_A_mediator","motif_B_biomarker","motif_C_exposure_marker",
              "motif_D_P_to_E","motif_E_disease_liability",
              "any_sig","n_motifs")
  miss <- setdiff(needed, names(df))
  if (length(miss) > 0) stop(tag, " triplet master missing columns: ", paste(miss, collapse=", "))
  df %>% mutate(dataset = tag)
}

ukb <- read_trip(TRIPLET_UKB_FP, "UKB_pQTL")
dec <- read_trip(TRIPLET_DEC_FP, "DECODE_pQTL")

if (KEEP_ONLY_ANY_SIG) {
  ukb <- ukb %>% filter(any_sig)
  dec <- dec %>% filter(any_sig)
}

# Make a stable join key (triplet column is already there, but we'll be safe)
ukb <- ukb %>% mutate(trip_key = ifelse(!is.na(triplet) & triplet != "", triplet,
                                        paste(Exposure, Protein, Disease, sep="||")))
dec <- dec %>% mutate(trip_key = ifelse(!is.na(triplet) & triplet != "", triplet,
                                        paste(Exposure, Protein, Disease, sep="||")))

# Keep one row per triplet key
ukb <- ukb %>% distinct(trip_key, .keep_all = TRUE)
dec <- dec %>% distinct(trip_key, .keep_all = TRUE)

# Inner join on the triplets that exist in BOTH (so overlap/unique is apples-to-apples)
both <- ukb %>%
  select(trip_key, Exposure, Protein, Disease,
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
  mutate(
    same_label = (motif_label_ukb == motif_label_dec)
  )

write.csv(both, file.path(OUTDIR, "Triplet_motif_calls_UKB_vs_DECODE_from_MRmotifs.csv"),
          row.names = FALSE)

# ============================================================
# 1) motif_label overlap / unique counts
# ============================================================
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

write.csv(motif_label_counts, file.path(OUTDIR, "MotifLabel_overlap_unique_counts.csv"),
          row.names = FALSE)

# Cross-tab transitions (UKB -> deCODe)
motif_label_xtab <- both %>%
  count(motif_label_ukb, motif_label_dec, name = "n") %>%
  complete(motif_label_ukb = labels_all, motif_label_dec = labels_all, fill = list(n = 0))

write.csv(motif_label_xtab, file.path(OUTDIR, "MotifLabel_crosstab_UKB_to_DECODE.csv"),
          row.names = FALSE)

motif_label_xtab_mat <- motif_label_xtab %>%
  pivot_wider(names_from = motif_label_dec, values_from = n) %>%
  arrange(factor(motif_label_ukb, levels = labels_all))

write.csv(motif_label_xtab_mat, file.path(OUTDIR, "MotifLabel_crosstab_matrix.csv"),
          row.names = FALSE)

# Plot: overlap/unique stacked by motif_label

##PLOT 1: Motif
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)

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
  coord_flip() +
  scale_y_log10(labels = comma) +
  scale_x_discrete(limits = positions) +
  theme_classic() +
  labs(x = NULL, y = "Triad count (log10)") +
  expand_limits(y = 2010440)#max(MotifUnion$n_union) * 20) +
  theme(plot.margin = margin(t=1, r=1, b=1, l=1, unit="cm"))

p_union <- ggplot(MotifUnion, aes(x = reorder(motif_label, n_union), y = n_union)) +
    geom_col() +
    geom_text(aes(label = label), hjust = -0.1, size = 3) +
    coord_flip(clip = "off") +   # key
    scale_y_log10(labels = scales::comma,
                  expand = expansion(mult = c(0, 0.28))) +  # extra space for text on the right
    scale_x_discrete(limits = positions) +
    theme_classic() +
    labs(x = NULL, y = "Triad count (log10)") +
    theme(
      plot.margin = margin(t=10, r=40, b=10, l=10),  # more right margin (pt)
      axis.text.x = element_text(size = 10),
      axis.text.y = element_text(size = 12)
    )

p_union

ggsave(file.path(OUTDIR, "MotifLabel_Union_DECODEUKB.png"),
       p_union, width = 6, height = 4, units = "in", dpi = 1000)

#PLOT 2:
comp_long <- motif_label_counts %>%
  select(motif_label, n_shared_same, n_unique_UKB, n_unique_DECODE) %>%
  pivot_longer(-motif_label, names_to = "component", values_to = "n") %>%
  mutate(component = recode(component,
                            n_shared_same   = "Shared",
                            n_unique_UKB    = "UKB pQTLs",
                            n_unique_DECODE = "deCODE pQTLs")) %>%
  #group_by(motif_label) %>%
  mutate(pct = n / sum(n), 
         label = paste0(percent(pct, 0.01), " (n=", comma(n), ")")) %>%
  ungroup()


p_components <- ggplot(
  comp_long,
  aes(x = reorder(motif_label, n, FUN = sum), y = n + 1, fill = component)
) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(
    aes(label = label),
    position = position_dodge(width = 0.8),
    hjust = -0.1,
    size = 3
  ) +
  coord_flip(clip = "off") +
  scale_y_log10(labels = comma,
                expand = expansion(mult = c(0, 0.30))) +
  scale_x_discrete(limits = positions) +
  theme_classic() +
  labs(x = NULL, y = "Triad count (log10)", title = "Shared vs Unique Triad Hits (UKB vs deCODE)") +
  theme(
    legend.position = "bottom",
    plot.margin = margin(t = 10, r = 45, b = 10, l = 10),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 12)
  )

p_components

ggsave(file.path(OUTDIR, "MotifLabel_Separate_DECODEUKB.png"),
       p_components, width = 7, height = 5, units = "in", dpi = 1000)

# Plot: transition heatmap
p2 <- ggplot(motif_label_xtab, aes(x = motif_label_dec, y = motif_label_ukb, fill = n)) +
  geom_tile() +
  theme_bw() +
  labs(x = "deCODe motif label", y = "UKB motif label",
       title = "Motif-label transitions (UKB → deCODe)", fill = "# triplets") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(OUTDIR, "MotifLabel_transition_heatmap.png"),
       p2, width = 8.5, height = 7, units = "in", dpi = 400)

# ============================================================
# 2) Overlap/unique per motif TYPE (A–E booleans), independent of label
# ============================================================
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

write.csv(motif_bool_summary, file.path(OUTDIR, "MotifTypeBoolean_overlap_unique_counts.csv"),
          row.names = FALSE)

cat("\nWrote motif overlap outputs to:\n  ", OUTDIR, "\n\n", sep="")
print(motif_label_counts)
print(motif_bool_summary)

