#!/usr/bin/env Rscript

# ============================================================
# Compare MR hit sets: UKB pQTLs vs deCODe pQTLs
# Outputs overlap/unique counts for:
#   - EP
#   - PDcis, PDtrans
#   - PEcis, PEtrans
#   - DP
#   - (and ED, DE as a sanity check — should be similar)
#
# Assumptions (match your current scripts):
#   - Files exist at:
#       MRfiles/summary/*res.csv           (UKB pQTLs)
#       MRfiles/summary/DECODE/*res.csv    (deCODe pQTLs)
#   - Each table has: method, pval, id.exposure, id.outcome
#   - PD/PE tables also have: edge_dir with Pcis_to_D, Ptrans_to_D, etc.
#
# What counts as a "hit" here:
#   - method in {IVW, Wald ratio}
#   - pval adjusted within each edge type via ADJ_METHOD
#   - significant if pval_adj < ALPHA
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(ggplot2)
})

# ----------------------------
# Config
# ----------------------------
MRfiles    <- "/n/groups/patel/shakson_ukb/UK_Biobank/Output/MR_edges"
ADJ_METHOD <- "BH"     # "BH" or "bonferroni"
ALPHA      <- 0.05
output_dir <- "/n/groups/patel/shakson_ukb/UK_Biobank/Output/MRplots/"

OUTDIR <- file.path(output_dir , "COMPARE_UKB_vs_DECODE")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

# ----------------------------
# Helpers
# ----------------------------
read_mr <- function(path) {
  # fread is fast; if you ever have TSVs swap to fread(..., sep="\t")
  fread(path) %>% as_tibble()
}

keep_method <- function(df) {
  if (!"method" %in% names(df)) stop("Missing 'method' column.")
  df %>%
    filter(method %in% c("Inverse variance weighted", "Wald ratio"))
}

add_q <- function(df, adj_method = "BH") {
  if (!"pval" %in% names(df)) stop("Missing 'pval' column.")
  df %>% mutate(pval_adj = p.adjust(pval, method = adj_method))
}

require_cols <- function(df, cols, label) {
  miss <- setdiff(cols, names(df))
  if (length(miss) > 0) stop(label, " is missing columns: ", paste(miss, collapse = ", "))
  df
}

mk_key <- function(df, type) {
  # standard key: exposure|outcome plus edge_dir when relevant
  # this ensures PDcis vs PDtrans don’t accidentally merge.
  df <- df %>%
    mutate(
      key_pair = paste0(id.exposure, "||", id.outcome),
      key = if ("edge_dir" %in% names(df)) paste0(key_pair, "||", edge_dir) else key_pair,
      edge_type = type
    )
  df
}

get_hit_keys <- function(df, type) {
  # returns distinct significant keys for that edge type
  df %>%
    distinct(key) %>%
    mutate(edge_type = type)
}

set_compare <- function(keysA, keysB) {
  # keysA/keysB are character vectors
  A <- unique(keysA)
  B <- unique(keysB)
  
  inter <- intersect(A, B)
  onlyA <- setdiff(A, B)
  onlyB <- setdiff(B, A)
  
  list(
    n_A = length(A),
    n_B = length(B),
    n_intersect = length(inter),
    n_onlyA = length(onlyA),
    n_onlyB = length(onlyB),
    intersect = inter,
    onlyA = onlyA,
    onlyB = onlyB
  )
}

# ----------------------------
# Load UKB (baseline) + DECODE (alt)
# ----------------------------
paths_ukb <- list(
  EP = file.path(MRfiles, "summary", "EPres.csv"),
  PD = file.path(MRfiles, "summary", "PDres.csv"),
  ED = file.path(MRfiles, "summary", "EDres.csv"),
  DE = file.path(MRfiles, "summary", "DEres.csv"),
  PE = file.path(MRfiles, "summary", "PEres.csv"),
  DP = file.path(MRfiles, "summary", "DPres.csv")
)

paths_dec <- list(
  EP = file.path(MRfiles, "summary", "DECODE", "EPres.csv"),
  PD = file.path(MRfiles, "summary", "DECODE", "PDres.csv"),
  ED = file.path(MRfiles, "summary", "DECODE", "EDres.csv"),
  DE = file.path(MRfiles, "summary", "DECODE", "DEres.csv"),
  PE = file.path(MRfiles, "summary", "DECODE", "PEres.csv"),
  DP = file.path(MRfiles, "summary", "DECODE", "DPres.csv")
)

# ----------------------------
# Preprocess function per dataset
# ----------------------------
prep_dataset <- function(paths, tag) {
  EP <- read_mr(paths$EP) %>%
    require_cols(c("id.exposure", "id.outcome", "method", "pval"), paste0(tag, " EP")) %>%
    keep_method() %>% add_q(ADJ_METHOD) %>% mk_key("EP")
  
  PD <- read_mr(paths$PD) %>%
    require_cols(c("id.exposure", "id.outcome", "method", "pval", "edge_dir"), paste0(tag, " PD")) %>%
    keep_method() %>% add_q(ADJ_METHOD) %>% mk_key("PD")
  
  ED <- read_mr(paths$ED) %>%
    require_cols(c("id.exposure", "id.outcome", "method", "pval"), paste0(tag, " ED")) %>%
    keep_method() %>% add_q(ADJ_METHOD) %>% mk_key("ED")
  
  DE <- read_mr(paths$DE) %>%
    require_cols(c("id.exposure", "id.outcome", "method", "pval"), paste0(tag, " DE")) %>%
    keep_method() %>% add_q(ADJ_METHOD) %>% mk_key("DE")
  
  PE <- read_mr(paths$PE) %>%
    require_cols(c("id.exposure", "id.outcome", "method", "pval", "edge_dir"), paste0(tag, " PE")) %>%
    keep_method() %>% add_q(ADJ_METHOD) %>% mk_key("PE")
  
  DP <- read_mr(paths$DP) %>%
    require_cols(c("id.exposure", "id.outcome", "method", "pval"), paste0(tag, " DP")) %>%
    keep_method() %>% add_q(ADJ_METHOD) %>% mk_key("DP")
  
  # Split cis/trans for PD and PE (matching your definitions)
  PDcis   <- PD %>% filter(edge_dir == "Pcis_to_D")  %>% mutate(edge_type = "PDcis")
  PDtrans <- PD %>% filter(edge_dir == "Ptrans_to_D")%>% mutate(edge_type = "PDtrans")
  PEcis   <- PE %>% filter(edge_dir == "Pcis_to_E")  %>% mutate(edge_type = "PEcis")
  PEtrans <- PE %>% filter(edge_dir == "Ptrans_to_E")%>% mutate(edge_type = "PEtrans")
  
  # For these splits, keep keys that include edge_dir so they stay disambiguated
  # but we also label edge_type so we can summarize cleanly.
  out <- list(
    EP = EP, ED = ED, DE = DE, DP = DP,
    PDcis = PDcis, PDtrans = PDtrans,
    PEcis = PEcis, PEtrans = PEtrans
  )
  
  # Add dataset tag for optional downstream merges
  out <- lapply(out, function(x) x %>% mutate(dataset = tag))
  out
}

ukb <- prep_dataset(paths_ukb, "UKB_pQTL")
dec <- prep_dataset(paths_dec, "DECODE_pQTL")

edge_types <- c("PDcis","PDtrans","PEcis","PEtrans","DP","EP","ED","DE")

# ----------------------------
# Build hit key-sets per edge type
# ----------------------------
hit_keys <- function(lst, edge_type) {
  lst[[edge_type]] %>%
    filter(pval_adj < ALPHA) %>%
    get_hit_keys(edge_type)
}

ukb_hits <- lapply(edge_types, function(et) hit_keys(ukb, et)) %>% setNames(edge_types)
dec_hits <- lapply(edge_types, function(et) hit_keys(dec, et)) %>% setNames(edge_types)

# ----------------------------
# Compare + summarize counts
# ----------------------------
summ_list <- lapply(edge_types, function(et) {
  A <- ukb_hits[[et]]$key
  B <- dec_hits[[et]]$key
  cmp <- set_compare(A, B)
  tibble(
    edge_type = et,
    n_hits_UKB = cmp$n_A,
    n_hits_DECODE = cmp$n_B,
    n_shared = cmp$n_intersect,
    n_unique_UKB = cmp$n_onlyA,
    n_unique_DECODE = cmp$n_onlyB,
    jaccard = ifelse((cmp$n_A + cmp$n_B - cmp$n_intersect) > 0,
                     cmp$n_intersect / (cmp$n_A + cmp$n_B - cmp$n_intersect),
                     NA_real_)
  )
})

summary_tbl <- bind_rows(summ_list) %>%
  mutate(edge_type = factor(edge_type, levels = edge_types))

# Write summary
write.csv(summary_tbl, file.path(OUTDIR, "MR_hit_overlap_summary.csv"), row.names = FALSE)

# ----------------------------
# Write the actual key lists (so you can inspect differences)
# ----------------------------
for (et in edge_types) {
  A <- ukb_hits[[et]]$key
  B <- dec_hits[[et]]$key
  cmp <- set_compare(A, B)
  
  fwrite(data.table(key = cmp$intersect), file.path(OUTDIR, paste0("shared_", et, ".tsv")), sep = "\t")
  fwrite(data.table(key = cmp$onlyA),     file.path(OUTDIR, paste0("unique_UKB_", et, ".tsv")), sep = "\t")
  fwrite(data.table(key = cmp$onlyB),     file.path(OUTDIR, paste0("unique_DECODE_", et, ".tsv")), sep = "\t")
}

# ----------------------------
# Optional: decode keys back into columns for convenience
# ----------------------------
# key format:
#   id.exposure||id.outcome          (EP/ED/DE/DP)
#   id.exposure||id.outcome||edge_dir (PDcis/PDtrans/PEcis/PEtrans)
decode_key <- function(dt) {
  parts <- str_split_fixed(dt$key, "\\|\\|", 3)
  out <- dt %>%
    mutate(
      id.exposure = parts[,1],
      id.outcome  = parts[,2],
      edge_dir    = ifelse(nchar(parts[,3]) == 0, NA_character_, parts[,3])
    )
  out
}

# Example: write a decoded table of shared keys for each edge type
for (et in edge_types) {
  fp <- file.path(OUTDIR, paste0("shared_", et, ".tsv"))
  dt <- fread(fp) %>% as_tibble()
  if (nrow(dt) == 0) next
  dt2 <- decode_key(dt) %>% mutate(edge_type = et)
  write.csv(dt2, file.path(OUTDIR, paste0("shared_", et, "_decoded.csv")), row.names = FALSE)
}

# ----------------------------
# Plot: shared vs unique bars by edge type
# ----------------------------
plot_df <- summary_tbl %>%
  select(edge_type, n_shared, n_unique_UKB, n_unique_DECODE) %>%
  pivot_longer(-edge_type, names_to = "category", values_to = "n") %>%
  mutate(category = recode(category,
                           n_shared = "Shared",
                           n_unique_UKB = "Unique to UKB pQTL",
                           n_unique_DECODE = "Unique to deCODe pQTL"))

p1 <- ggplot(plot_df, aes(x = edge_type, y = n, fill = category)) +
  geom_col(position = "stack") +
  coord_flip() +
  theme_bw() +
  labs(
    x = NULL,
    y = paste0("# significant hits (q<", ALPHA, " ; ", ADJ_METHOD, ")"),
    title = "MR hit overlap: UKB pQTL vs deCODe pQTL"
  )

ggsave(file.path(OUTDIR, "MR_hit_overlap_stacked.png"),
       p1, width = 7, height = 4.5, units = "in", dpi = 400)
ggsave(file.path(OUTDIR, "MR_hit_overlap_stacked.svg"),
       p1, width = 7, height = 4.5, units = "in")

# Plot: Jaccard similarity per edge type
p2 <- ggplot(summary_tbl, aes(x = edge_type, y = jaccard)) +
  geom_point(size = 2) +
  coord_flip() +
  theme_bw() +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    x = NULL,
    y = "Jaccard (shared / union)",
    title = "Similarity of MR hit sets by edge type"
  )

ggsave(file.path(OUTDIR, "MR_hit_overlap_jaccard.png"),
       p2, width = 6.5, height = 4.0, units = "in", dpi = 400)
ggsave(file.path(OUTDIR, "MR_hit_overlap_jaccard.svg"),
       p2, width = 6.5, height = 4.0, units = "in")

# ----------------------------
# Print to console
# ----------------------------
cat("\nWrote outputs to:\n  ", OUTDIR, "\n\n", sep = "")
print(summary_tbl)

cat("\nNotes:\n")
cat(" - PD/PE are split by edge_dir; keys include edge_dir so cis/trans stay separated.\n")
cat(" - ED/DE should be similar across pQTL sources; if not, it may indicate different filtering/IDs.\n")
cat(" - If you want to compare *any* PD (cis OR trans), you can union PDcis+PDtrans before comparing.\n")
