#!/usr/bin/env Rscript
# ============================================================================
# summarize_module1_component_r2.R
# ----------------------------------------------------------------------------
# Per-protein variance-component R2 with confidence intervals -- the per-entity
# table behind Supplementary Table S1 (see docs/SUPPLEMENT_DATA_STANDARD.md §2.1:
# supplementary tables report one row per protein, not a distributional summary).
#
# Two granularities, both one row per protein:
#   component_r2_per_protein.tsv       C / G / E / GxE  + whole-model total
#   component_r2_per_protein_fine.tsv  C / Gcis / Gtrans / E-by-category
#                                      / GxE-by-category   (29 components)
#
# Each component contributes TWO columns: a numeric `R2 <comp>` (sortable) and a
# text `CI <comp>` of the form "(lo, hi)".
#
# R2 is the mean over the K out-of-fold estimates of the leave-one-out drop
# decomposition; the interval is a t-based 95% CI on that mean (K-1 df). Unique
# R2 can be NEGATIVE when a component adds nothing -- the reduced model beats the
# full model out-of-fold. Values are NOT clipped at zero: clipping would bias
# every component upward and make the interval meaningless.
#
# Run:
#   Rscript scripts/analysis_summaries/summarize_module1_component_r2.R \
#          [covarType] [method] [experiment]
#   Defaults: base / lasso / M1_base_lasso   (the base specification -> Tier 1)
#   Alternative specs (-> Supplementary Data) e.g.:
#          base_bmi lasso M1_base_bmi_lasso
# ============================================================================
suppressPackageStartupMessages(library(data.table))
options(scipen = 999)   # write 0.0009, not 9e-04, in the shipped TSVs

local({
  cand <- c(Sys.getenv("HEAP_PATHS_FILE", unset = ""),
            "/n/groups/patel/shakson_ukb/HEAP/workflow/00_paths.R")
  hit <- cand[nzchar(cand) & file.exists(cand)][1]; if (!is.na(hit)) source(hit)
})
source(file.path(heap_root, "scripts", "visualizations", "common", "load_heap_results.R"))
source(file.path(heap_root, "scripts", "visualizations", "common", "label_helpers.R"))  # heap_broad_category

a          <- commandArgs(trailingOnly = TRUE)
covarType  <- if (length(a) >= 1) a[1] else "base"
method     <- if (length(a) >= 2) a[2] else "lasso"
experiment <- if (length(a) >= 3) a[3] else "M1_base_lasso"

outdir <- file.path(heap_root, "docs", "manuscript_stats", "module1_variance")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
suffix <- if (identical(experiment, "M1_base_lasso")) "" else paste0("_", experiment)

message("Module 1 component R2 (", experiment, "/", covarType, "/", method, ") ...")
ld <- function(lev) load_module1_predictive_r2(covarType, method, level = lev,
                                               experiment = experiment)

## mean over folds + t-based 95% CI on the mean
agg <- function(d, by) {
  d <- d[is.finite(r2)]
  d[, .(K = .N, r2 = mean(r2), sd = sd(r2)), by = by
    ][, `:=`(lo = r2 - qt(0.975, K - 1) * sd / sqrt(K),
             hi = r2 + qt(0.975, K - 1) * sd / sqrt(K))][]
}

## append a numeric "R2 <lab>" + text "CI <lab>" pair onto `dst`, keyed by protein
add_pair <- function(dst, src, lab, dgt = 4) {
  i <- match(dst$Protein, src$omic)
  set(dst, j = paste("R2", lab), value = round(src$r2[i], dgt))
  set(dst, j = paste("CI", lab),
      value = ifelse(is.na(i), NA_character_,
                     sprintf("(%.*f, %.*f)", dgt, src$lo[i], dgt, src$hi[i])))
  dst
}

co   <- ld("coarse")
uniq <- agg(co[method == "score_unique_drop"], c("omic", "block"))
tot  <- agg(co[method == "score_model_total" & block == "C+G+E+GxE"], "omic")
prot <- sort(unique(co$omic))

# ---- sheet A: coarse -------------------------------------------------------
A <- data.table(Protein = prot)
A <- add_pair(A, tot, "whole model")
lab_map <- c(Covars = "covariates", G = "genetic", E = "exposomic", GxE = "GxE")
for (b in names(lab_map)) A <- add_pair(A, uniq[block == b], lab_map[[b]])

# ---- sheet B: fine ---------------------------------------------------------
gs <- agg(ld("genetic_subblocks")[method == "score_unique_drop"], c("omic", "block"))
ec <- agg(ld("exposure_categories")[method == "score_unique_drop"], c("omic", "category"))
gc <- agg(ld("gxe_categories")[method == "score_unique_drop"], c("omic", "category"))

B <- data.table(Protein = prot)
B <- add_pair(B, uniq[block == "Covars"], "covariates")
for (b in c("Gcis", "Gtrans")) B <- add_pair(B, gs[block == b], b)
for (k in sort(unique(ec$category))) B <- add_pair(B, ec[category == k], paste0("E: ", k))
for (k in sort(unique(gc$category))) B <- add_pair(B, gc[category == k], paste0("GxE: ", k))

# ---- dominant GxE category (carried over from the retired gxe_proteins) -----
# argmax over the 13 per-category GxE R2. Computed on UNROUNDED values: at the
# 4-dp display rounding 214/2686 proteins tie spuriously. Genuine ties (a protein
# whose top categories are exactly equal, typically both 0) break alphabetically,
# and n_tied records how many categories shared the maximum.
domg <- gc[order(omic, -r2, category),
           .(dominant_GxE_category = category[1],
             n_tied = sum(r2 == max(r2, na.rm = TRUE), na.rm = TRUE)), by = omic]
i <- match(A$Protein, domg$omic)
A[, `dominant GxE category` := domg$dominant_GxE_category[i]]
A[, `dominant GxE ties`     := domg$n_tied[i]]
message(sprintf("  dominant GxE category: %d proteins, %d with a tie at the maximum",
                nrow(domg), sum(domg$n_tied > 1L, na.rm = TRUE)))

# ---- sheet C: distributional summary, DERIVED from the per-protein tables ---
# One place for every "across the proteome" number: replaces the old
# component_summary.tsv + exposure_categories.tsv, and is by construction
# consistent with sheets A and B because it is computed from them.
# n_proteins_pos and n_proteins_gt_0.005 are carried over from the retired
# category_r2 table; they are counted on UNROUNDED r2, so they are not
# recoverable from the 4-dp per-protein sheets.
dist <- function(d, section, label_col) {
  d[, .(section = section, n_proteins = sum(is.finite(r2)),
        mean = round(mean(r2, na.rm = TRUE), 4),
        median = round(median(r2, na.rm = TRUE), 4),
        q25 = round(quantile(r2, .25, na.rm = TRUE), 4),
        q75 = round(quantile(r2, .75, na.rm = TRUE), 4),
        max = round(max(r2, na.rm = TRUE), 4),
        n_proteins_pos      = sum(r2 > 0, na.rm = TRUE),
        n_proteins_gt_0.005 = sum(r2 > 0.005, na.rm = TRUE),
        pct_gt_0.01 = round(100 * mean(r2 > 0.01, na.rm = TRUE), 1),
        pct_gt_0.05 = round(100 * mean(r2 > 0.05, na.rm = TRUE), 1),
        pct_negative = round(100 * mean(r2 < 0, na.rm = TRUE), 1)),
    by = c(label_col)]
}
tot2 <- copy(tot)[, block := "Whole model"]
C <- rbindlist(list(
  dist(rbind(uniq[, .(omic, block, r2)], tot2[, .(omic, block, r2)]), "Component", "block"),
  dist(gs, "Genetic sub-block",  "block"),
  dist(ec, "Exposure category",  "category"),
  dist(gc, "GxE category",       "category")), fill = TRUE)
C[, component := fifelse(is.na(block), category, block)][, c("block", "category") := NULL]
# broad grouping (Lifestyle: … / Environmental) for the exposure categories,
# from the canonical config via heap_broad_category() -- NOT copied out of the
# category_r2 table this replaces.
C[, broad := NA_character_]
isec <- C$section %in% c("Exposure category", "GxE category")
C[isec, broad := tryCatch(heap_broad_category(component),
                          error = function(e) NA_character_)]
setcolorder(C, c("section", "component", "broad"))
C[, component := c(Covars = "Covariates", G = "Genetic", E = "Exposomic (PXS)",
                   GxE = "GxE")[component] |> (\(x) fifelse(is.na(x), C$component, x))()]

fa <- file.path(outdir, sprintf("component_r2_per_protein%s.tsv", suffix))
fb <- file.path(outdir, sprintf("component_r2_per_protein_fine%s.tsv", suffix))
fc <- file.path(outdir, sprintf("component_r2_summary%s.tsv", suffix))
fwrite(A, fa, sep = "\t"); fwrite(B, fb, sep = "\t"); fwrite(C, fc, sep = "\t")
message(sprintf("  %s  %d x %d", basename(fc), nrow(C), ncol(C)))

message(sprintf("  %s  %d x %d", basename(fa), nrow(A), ncol(A)))
message(sprintf("  %s  %d x %d", basename(fb), nrow(B), ncol(B)))
message(sprintf("  K = %d folds; %d negative unique-R2 cells retained (not clipped)",
                max(uniq$K),
                sum(uniq$r2 < 0) + sum(gs$r2 < 0) + sum(ec$r2 < 0) + sum(gc$r2 < 0)))
