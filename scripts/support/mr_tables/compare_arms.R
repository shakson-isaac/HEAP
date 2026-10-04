#!/usr/bin/env Rscript
# ============================================================================
# support/mr_tables/compare_arms.R
# ----------------------------------------------------------------------------
# Cross-arm (UKB pQTL vs deCODE pQTL) comparison of MR triad motif calls.
# Reads the per-arm wide motif tables written by build_mr_tables.R
#   mr_edges/summary/MRmotifs.tsv          (UKB)
#   mr_edges/summary/DECODE/MRmotifs.tsv   (deCODE)
# and emits the shared-vs-unique tables the figures consume (this is the cross-arm
# STATISTIC, kept out of the plotter):
#   mr_edges/summary/arm_comparison_motif_counts.tsv    long (motif_label, component, n)
#   mr_edges/summary/arm_comparison_label_wide.tsv      wide per-label counts
#   mr_edges/summary/arm_comparison_motif_transition.tsv UKB->deCODE label crosstab
#   mr_edges/summary/arm_comparison_motif_type_jaccard.tsv per-motif-type (A-E) overlap
#
# Comparison is on the INNER set of triads present in BOTH arms (apples-to-apples),
# ported from ModuleMR/COMBO/comparePQTL_viz2.R.
#
# Usage:
#   module load gcc/14.2.0 R/4.4.2
#   HEAP_PATHS_FILE=.../workflow/00_paths.R \
#     Rscript scripts/support/mr_tables/compare_arms.R
# ============================================================================
local({
  cand <- c(Sys.getenv("HEAP_PATHS_FILE", unset = ""),
            file.path(getwd(), "workflow", "00_paths.R"),
            "/n/groups/patel/shakson_ukb/HEAP/workflow/00_paths.R")
  hit <- cand[nzchar(cand) & file.exists(cand)][1]
  if (is.na(hit)) stop("Could not locate workflow/00_paths.R (set HEAP_PATHS_FILE).")
  source(hit)
})
suppressPackageStartupMessages({ library(data.table) })

summ_dir  <- heap_project_output("mr_edges", "summary")
ukb_fp    <- file.path(summ_dir, "MRmotifs.tsv")
dec_fp    <- file.path(summ_dir, "DECODE", "MRmotifs.tsv")
for (p in c(ukb_fp, dec_fp))
  if (!file.exists(p)) stop("Missing MRmotifs: ", p,
    "\nRun build_mr_tables.R for both arms first.", call. = FALSE)

MOTIF_BOOL <- c(A = "motif_A_mediator", B = "motif_B_biomarker",
                C = "motif_C_exposure_marker", D = "motif_D_P_to_E",
                E = "motif_E_disease_liability")
keep <- c("triplet", "Exposure", "Protein", "Disease", "motif_label",
          "any_sig", "n_motifs", unname(MOTIF_BOOL))

read_arm <- function(fp) {
  x <- fread(fp, select = keep)
  x <- unique(x, by = "triplet")
  x
}
ukb <- read_arm(ukb_fp)
dec <- read_arm(dec_fp)
message(sprintf("UKB triads=%d  deCODE triads=%d", nrow(ukb), nrow(dec)))

# inner set present in both arms
both <- merge(
  ukb[, c("triplet", "motif_label", "any_sig", unname(MOTIF_BOOL)), with = FALSE],
  dec[, c("triplet", "motif_label", "any_sig", unname(MOTIF_BOOL)), with = FALSE],
  by = "triplet", suffixes = c("_ukb", "_dec"))
message(sprintf("triads in BOTH arms: %d", nrow(both)))

labels_all <- sort(unique(c(both$motif_label_ukb, both$motif_label_dec)))
labels_all <- labels_all[!is.na(labels_all)]

# ---- per-label shared / unique counts --------------------------------------
lab_wide <- rbindlist(lapply(labels_all, function(m) {
  data.table(
    motif_label     = m,
    n_UKB           = sum(both$motif_label_ukb == m, na.rm = TRUE),
    n_DECODE        = sum(both$motif_label_dec == m, na.rm = TRUE),
    n_shared_same   = sum(both$motif_label_ukb == m & both$motif_label_dec == m, na.rm = TRUE),
    n_unique_UKB    = sum(both$motif_label_ukb == m & both$motif_label_dec != m, na.rm = TRUE),
    n_unique_DECODE = sum(both$motif_label_dec == m & both$motif_label_ukb != m, na.rm = TRUE))
}))
fwrite(lab_wide, file.path(summ_dir, "arm_comparison_label_wide.tsv"), sep = "\t")

# long form for the grouped bar (Shared / UKB pQTLs / deCODE pQTLs)
counts_long <- melt(
  lab_wide[, .(motif_label,
               Shared          = n_shared_same,
               `UKB pQTLs`     = n_unique_UKB,
               `deCODE pQTLs`  = n_unique_DECODE)],
  id.vars = "motif_label", variable.name = "component", value.name = "n")
counts_long[, pct := n / sum(n)]
fwrite(counts_long, file.path(summ_dir, "arm_comparison_motif_counts.tsv"), sep = "\t")

# ---- UKB -> deCODE label transition crosstab -------------------------------
xtab <- both[, .N, by = .(motif_label_ukb, motif_label_dec)]
fwrite(xtab, file.path(summ_dir, "arm_comparison_motif_transition.tsv"), sep = "\t")

# ---- per-motif-TYPE (A-E boolean) overlap + Jaccard ------------------------
type_tab <- rbindlist(lapply(names(MOTIF_BOOL), function(k) {
  cu <- as.logical(both[[paste0(MOTIF_BOOL[[k]], "_ukb")]])
  cd <- as.logical(both[[paste0(MOTIF_BOOL[[k]], "_dec")]])
  cu[is.na(cu)] <- FALSE; cd[is.na(cd)] <- FALSE
  ns <- sum(cu & cd); nu <- sum(cu & !cd); nd <- sum(cd & !cu)
  data.table(motif_type = k, n_UKB = sum(cu), n_DECODE = sum(cd),
             n_shared = ns, n_unique_UKB = nu, n_unique_DECODE = nd,
             jaccard = if ((sum(cu) + sum(cd) - ns) > 0) ns / (sum(cu) + sum(cd) - ns) else NA_real_)
}))
fwrite(type_tab, file.path(summ_dir, "arm_comparison_motif_type_jaccard.tsv"), sep = "\t")

# ============================================================================
# CROSS-ARM TIERING: replication -> Tier 1+; per-lane attrition funnel
# (mr_tier is computed per-arm by build_mr_tables.R; Tier 1+ needs both arms)
# ============================================================================
read_long <- function(co) {
  f <- file.path(if (co == "UKB") summ_dir else file.path(summ_dir, "DECODE"),
                 "mr_sensitivity_long.tsv")
  if (!file.exists(f)) stop("Missing ", f, " — run build_mr_tables.R ", co, " first.")
  fread(f)
}
ukb_l <- read_long("UKB"); dec_l <- read_long("DECODE")
key <- c("edge_dir", "src_id", "tgt_id")

# replicated = same edge reaches Tier 1 in BOTH arms; protein-free edges
# (E_to_D/D_to_E) are byte-identical across arms -> not eligible for replication.
mt <- merge(ukb_l[, c(key, "mr_tier"), with = FALSE],
            dec_l[, c(key, "mr_tier"), with = FALSE],
            by = key, suffixes = c("_ukb", "_dec"), all = TRUE)
mt[, replicated := !(edge_dir %in% c("E_to_D", "D_to_E")) &
                   mr_tier_ukb == "Tier1" & mr_tier_dec == "Tier1"]
mt[is.na(replicated), replicated := FALSE]

finalize <- function(long, tag) {
  x <- merge(long, mt[, c(key, "replicated"), with = FALSE], by = key, all.x = TRUE)
  x[is.na(replicated), replicated := FALSE]
  x[, mr_tier_final := fifelse(mr_tier == "Tier1" & replicated, "Tier1plus", mr_tier)]
  x[, dataset := tag]; x
}
tiered <- rbind(finalize(ukb_l, "UKB"), finalize(dec_l, "DECODE"), fill = TRUE)

# --- fold in targeted colocalization (run_coloc_shortlist.R), if present ------
# cis-pQTL Tier 1/1+ that are NOT coloc-confirmed demote to Tier 2 (flowchart).
tiered[, coloc_pph4 := NA_real_]
coloc_fp <- file.path(heap_project_output("support", "coloc"), "coloc_shortlist_results.tsv")
if (file.exists(coloc_fp)) {
  cz <- fread(coloc_fp)[, .(dataset, edge_dir, src_id = protein, tgt_id = outcome,
                            cz_pph4 = `PP.H4`, cz_status = coloc_status)]
  tiered <- merge(tiered, cz, by = c("dataset", "edge_dir", "src_id", "tgt_id"), all.x = TRUE)
  tiered[!is.na(cz_status), coloc_status := cz_status]
  tiered[!is.na(cz_pph4),  coloc_pph4 := cz_pph4]
  # demote ONLY cis Tier1 that were coloc-ASSESSED and not confirmed (distinct/
  # ambiguous/failed). cis Tier1 not yet coloc'd stay Tier1 with status "pending"
  # (coloc not-attempted != refuted) -> they form the next coloc shortlist.
  tiered[lane == "pQTL" & edge_class == "cis" &
         mr_tier_final %in% c("Tier1", "Tier1plus") &
         !is.na(cz_status) & cz_status != "confirmed", mr_tier_final := "Tier2"]
  tiered[, c("cz_pph4", "cz_status") := NULL]
  message("coloc folded in: ", cz[, sum(cz_status == "confirmed")], " confirmed of ", nrow(cz))
}
tiered[lane == "pQTL" & edge_class == "cis" &
       mr_tier_final %in% c("Tier1", "Tier1plus") & is.na(coloc_pph4),
       coloc_status := "pending"]

keep_t <- intersect(c("dataset", key, "lane", "edge_class", "nsnp", "mr_tier",
                      "mr_tier_final", "replicated", "tier_reason", "het_status",
                      "coloc_status", "coloc_pph4"), names(tiered))
fwrite(tiered[, ..keep_t], file.path(summ_dir, "mr_tiered_edges.tsv"), sep = "\t")
fwrite(tiered[, .N, by = .(dataset, lane, mr_tier_final)],
       file.path(summ_dir, "tier_counts.tsv"), sep = "\t")

TIER_RANK <- c(Null = 0L, Suggestive = 1L, Reverse = 1L, Tier2 = 2L, Tier1 = 3L, Tier1plus = 4L)
tiered[, trank := TIER_RANK[mr_tier_final]]
funnel <- tiered[, .(tested = .N,
                     significant = sum(mr_tier_final != "Null"),
                     qualified   = sum(trank >= 2L),   # Tier 2+
                     tier1       = sum(trank >= 3L),    # Tier 1 / 1+
                     replicated  = sum(trank >= 4L)),   # Tier 1+
                 by = .(dataset, lane)]
fwrite(funnel, file.path(summ_dir, "tier_funnel.tsv"), sep = "\t")

cat("\n================ TIER FUNNEL (per lane) ================\n")
print(funnel[order(dataset, lane)])
cat("\nFinal tier counts:\n")
print(dcast(tiered, lane + dataset ~ mr_tier_final, value.var = "trank", fun.aggregate = length))

cat("\n================ ARM COMPARISON ================\n")
print(lab_wide)
cat("\nMotif-type overlap (A-E):\n"); print(type_tab)
cat("\nWrote arm_comparison_* to: ", summ_dir, "\nDONE.\n", sep = "")
