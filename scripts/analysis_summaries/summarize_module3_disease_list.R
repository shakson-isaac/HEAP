#!/usr/bin/env Rscript
# ============================================================================
# summarize_module3_disease_list.R
# ----------------------------------------------------------------------------
# Regenerates the Module 3 incident-disease inventory -- the "181 diseases" set
# behind Supplementary Table S20 (`disease_list`) and behind
# multipes_disease/disease_set_181.tsv, which module6_pes_disease_scale.R and
# module6_pes_disease_dc_holdout.R both read.
#
# WHY THIS EXISTS. Neither file had a generator: both were flat files that
# nothing in the repo wrote (SUPP_TABLE_REVIEW.md item E6). That is the same
# shape as ladder_candidates.tsv (L9) and the retired mediation builder (L3) --
# readable, plausible, and impossible to check against a rerun.
#
# WHAT THE SET IS. Not a standalone cohort calculation. Module3.R applies
# MIN_CASES_PER_DZ = 100 per (protein, disease) AFTER complete-casing on
# protein, survival time, predictors and covariates (Module3.R:589), so the
# disease inventory is an EMERGENT property of the run: a disease is in the set
# iff at least one protein cleared 100 cases for it. This reads that set back
# off the canonical base-specification link table rather than recomputing a
# cohort, which is why it agrees with the mediation output by construction.
#
# WHAT IS VERIFIED vs WHAT IS NOT (checked against the shipped file):
#   * the 181-disease SET       -- reproduces EXACTLY, and under every
#                                  aggregation tried, so it is not tie-break
#                                  dependent
#   * n_cases = max over proteins -- reproduces EXACTLY for all 181
#   ! n_participants            -- the rule below matches 140/181 exactly. The
#                                  other 41 come out HIGHER than the shipped
#                                  file by 1..191 participants (<=1.5%), always
#                                  in that direction.
#
# WHY, run to ground 2026-08-08 (PRV-1). Do not spend another pass guessing rules:
#   * NOT a different statistic. Seven aggregations tested (max/min/median over all
#     proteins, min among argmax-cases, n at min cases, first protein alphabetically).
#     The rule below is the best of them; every alternative scores worse than 140.
#   * NOT a tie-break. Five orderings tested (cache order first/last, protID sort
#     first/last, modal n among ties). All worse.
#   * NOT a predictor_class subset. All five classes give byte-identical results --
#     n and n_cases do not vary by predictor_class.
#   * The shipped value IS a genuine per-protein n from the tied-at-max-cases rows,
#     for 181/181 diseases. So it is a real complete-case sample, not an artifact.
#   => the shipped column came from a NARROWER PROTEIN ROSTER than this cache holds.
#     443 proteins carry the excess across the 41 diseases; no roster explains all 41
#     (the widest covers 21), and the excess proteins have systematically larger n
#     (median 12,258 vs 10,550) -- a panel/batch missingness signature.
#   The original roster is not recoverable from .med_links_base_cache.rds. The
#   discrepancy is <=1.5% in a column nothing cites, so the shipped file stands.
#
#   Note what the column MEANS: not a per-disease cohort, but the complete-case n of
#   the protein attaining that disease's maximum case count. S20's legend implies the
#   former; a sharpening is proposed at HEAP_manuscript/docs/PENDING_CITATIONS.txt #8.
#
# So this script does NOT overwrite anything. It writes *_regen.tsv beside the
# shipped files and prints the diff. Promoting a regenerated file over a cited
# one is a rule-5 decision (SUPPLEMENT_README.md) and belongs to the author.
#
# Input: the base-spec link cache written by summarize_module3_mediation.R.
#        Run that first if it is missing -- it derives from the canonical run
#        directories, and this script deliberately does not duplicate that
#        loading logic.
#
#   Run:  Rscript scripts/analysis_summaries/summarize_module3_disease_list.R
# ============================================================================
local({ cand <- c(Sys.getenv("HEAP_PATHS_FILE",""),"/n/groups/patel/shakson_ukb/HEAP/workflow/00_paths.R")
  h <- cand[nzchar(cand)&file.exists(cand)][1]; if(!is.na(h)) source(h) })
suppressPackageStartupMessages(library(data.table))
msg <- function(...) cat(format(Sys.time(),"[%H:%M:%S] "),...,"\n",sep="")

MIN_CASES <- as.integer(Sys.getenv("HEAP_MIN_CASES_PER_DZ", unset = "100"))
STAT  <- file.path(heap_root, "docs", "manuscript_stats", "module3")
M6DIR <- heap_project_output("module6_pes_longitudinal", "multipes_disease")
CACHE <- heap_project_intermediate(".med_links_base_cache.rds")

if (!file.exists(CACHE))
  stop("missing ", CACHE, "\n  run scripts/analysis_summaries/summarize_module3_mediation.R first")

msg("loading base-spec link table")
L <- readRDS(CACHE); setDT(L)
msg(format(nrow(L), big.mark=","), " links | ", uniqueN(L$protID), " proteins | ",
    uniqueN(L$DZ_ID), " diseases")

# Module3.R already applied the >=100 filter per (protein, disease); re-assert it
# so the rule is visible here and the script stays correct if the cache ever
# carries a lower-threshold run.
L <- L[n_cases >= MIN_CASES]

# n_cases: the maximum across proteins -- the count attained by the protein with
#          the most complete data. Verified exact against the shipped file.
# n_participants: the largest complete-case sample among the rows attaining that
#          maximum. Matches 140/181; see the header note.
D <- L[, {
  mx <- max(n_cases)
  .(n_participants = max(n[n_cases == mx]), n_cases = mx)
}, by = DZ_ID]
setorder(D, DZ_ID)
msg(nrow(D), " diseases at >=", MIN_CASES, " incident cases")

dir.create(STAT, recursive = TRUE, showWarnings = FALSE)
out_list <- file.path(STAT,  "disease_list_regen.tsv")
out_set  <- file.path(M6DIR, "disease_set_181_regen.tsv")
fwrite(D, out_list, sep = "\t")
fwrite(D[, .(DZ_ID)], out_set, sep = "\t")
msg("wrote ", out_list)
msg("wrote ", out_set)

# ---- diff against the shipped files; never overwrite them -------------------
ref_f <- file.path(STAT, "disease_list.tsv")
if (file.exists(ref_f)) {
  R <- fread(ref_f)
  cat("\n=== diff vs shipped disease_list.tsv ===\n")
  cat("  diseases: shipped ", nrow(R), " | regenerated ", nrow(D), "\n", sep="")
  extra <- setdiff(D$DZ_ID, R$DZ_ID); missing <- setdiff(R$DZ_ID, D$DZ_ID)
  cat("  set identical: ", length(extra) == 0 && length(missing) == 0, "\n", sep="")
  if (length(extra))   cat("    only in regenerated: ", paste(substr(extra,1,44),   collapse=", "), "\n")
  if (length(missing)) cat("    only in shipped:     ", paste(substr(missing,1,44), collapse=", "), "\n")
  M <- merge(R, D, by = "DZ_ID", suffixes = c("_shipped","_regen"))
  dc <- M[n_cases_shipped != n_cases_regen]
  dn <- M[n_participants_shipped != n_participants_regen]
  cat("  n_cases differing:        ", nrow(dc), " / ", nrow(M), "\n", sep="")
  cat("  n_participants differing: ", nrow(dn), " / ", nrow(M), "\n", sep="")
  if (nrow(dn)) {
    dn[, delta := n_participants_regen - n_participants_shipped]
    cat("    delta range: ", min(dn$delta), " .. ", max(dn$delta),
        "  (regen - shipped; expected all >= 0, <= ~191)\n", sep="")
    print(head(dn[order(-abs(delta)), .(DZ_ID = substr(DZ_ID,5,40),
              shipped = n_participants_shipped, regen = n_participants_regen, delta)], 5))
  }
  if (nrow(dc)) stop("n_cases drifted -- that is a cited number; investigate before using this output")
  cat("\n  NOT overwriting the shipped files. Promoting *_regen.tsv is a rule-5 decision.\n")
}
