#!/usr/bin/env Rscript
# ============================================================================
# summarize_module2_supp_tables.R
# ----------------------------------------------------------------------------
# Supplementary Tables S8 (exposure-protein) and S11 (polygenic GxE), built from
# the SAME code path the manuscript figures use -- load_module2_results() over
# the experiment-nested output module2/M2_base_main/base/univar_assoc_*.rds.
#
# Replaces the flat root CSVs (module2/Replicated{E,GxE}assoc_base.csv), which
# no current script writes and which are 11 rows behind the nested computation:
# they were produced before the experiment-nesting restructure with a slightly
# different Bonferroni denominator, flipping borderline pairs.
#
# UNIT: one row per (model TERM x protein). Ordinal exposures contribute one row
# per polynomial-contrast level, so every row carries an interpretable beta.
#
# IMPORTANT -- two different replication tests exist, and this table is the
# term-level one:
#   * per-TERM t-test  (Pr(>|t|) on the term)      -> the rows of this table
#   * per-PAIR block F-test (p_E_block, omnibus)   -> the manuscript's headline
#     count of 22,240 replicated exposure-protein associations (\nReplAssoc)
# The block p-values are carried as columns so a reader can see both, BUT the
# 22,240 figure is NOT recoverable from this table: 1,747 pairs are block-
# significant with no individual term passing the stricter term threshold, and a
# further 36 have no statE rows at all. The legend must say so; do not claim the
# headline count can be derived here.
#
#   Run:  Rscript scripts/analysis_summaries/summarize_module2_supp_tables.R
# ============================================================================
suppressPackageStartupMessages(library(data.table))

local({
  cand <- c(Sys.getenv("HEAP_PATHS_FILE", unset = ""),
            "/n/groups/patel/shakson_ukb/HEAP/workflow/00_paths.R")
  hit <- cand[nzchar(cand) & file.exists(cand)][1]; if (!is.na(hit)) source(hit)
})
source(file.path(heap_root, "scripts", "visualizations", "common", "load_heap_results.R"))

a          <- commandArgs(trailingOnly = TRUE)
covarType  <- if (length(a) >= 1) a[1] else "base"
experiment <- if (length(a) >= 2) a[2] else "M2_base_main"
OUT <- file.path(heap_root, "docs", "manuscript_stats", "module2_associations")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

message("Loading Module 2 (", experiment, "/", covarType, ") train + test ...")
tr <- load_module2_results(covarType, "train", FALSE, experiment)
te <- load_module2_results(covarType, "test",  FALSE, experiment)

## ---- per-pair block F-test (context columns; defines the headline count) ----
eb <- merge(tr$statFblock[, .(Eid = ID, omicID, pb_tr = p_E_block, pgx_tr = p_GxE_joint)],
            te$statFblock[, .(Eid = ID, omicID, pb_te = p_E_block, pgx_te = p_GxE_joint)],
            by = c("Eid", "omicID"))
ebE <- eb[is.finite(pb_tr) & is.finite(pb_te)]
thr_block <- 0.05 / nrow(ebE)
message(sprintf("  block F-test: %d pairs tested, thr %.3g, %d replicated",
                nrow(ebE), thr_block, ebE[pb_tr < thr_block & pb_te < thr_block, .N]))

## ---- S8: term-level replicated exposure-protein associations ---------------
m   <- load_module2_replicated(covarType, experiment = experiment)
thrT <- attr(m, "pval_thresh")
s8  <- merge(m[replicated == TRUE], ebE[, .(Eid, omicID, pb_tr, pb_te)],
             by = c("Eid", "omicID"), all.x = TRUE)
setorder(s8, p_train)
S8 <- s8[, .(Exposure = Eid, Term = ID, Protein = omicID, Category,
             `N train` = samplesize_train, `N test` = samplesize_test,
             `beta train` = signif(beta_train, 4), `SE train` = signif(se_train, 4),
             `t train` = signif(t_train, 4), `p train` = p_train,
             `beta test` = signif(beta_test, 4), `SE test` = signif(se_test, 4),
             `t test` = signif(t_test, 4), `p test` = p_test,
             `p block train` = pb_tr, `p block test` = pb_te)]
fwrite(S8, file.path(OUT, "assoc_replicated_terms.tsv"), sep = "\t")
message(sprintf("  assoc_replicated_terms.tsv  %d rows x %d cols (term thr %.3g)",
                nrow(S8), ncol(S8), thrT))

## ---- S11: term-level replicated polygenic GxE ------------------------------
gx_pull <- function(r) r$statGxE[, .(ID, omicID, Eid, Category, E_term, G_component,
                                     beta = Estimate, se = `Std. Error`,
                                     tval = `t value`, p = `Pr(>|t|)`, samplesize)]
gtr <- gx_pull(tr); gte <- gx_pull(te)
key <- c("ID", "omicID", "G_component")
g <- merge(gtr, gte[, c(key, "beta", "se", "tval", "p", "samplesize"), with = FALSE],
           by = key, suffixes = c("_train", "_test"))
thrG <- 0.05 / nrow(g)
g <- g[is.finite(p_train) & is.finite(p_test) & p_train < thrG & p_test < thrG]
g <- merge(g, ebE[, .(Eid, omicID, pgx_tr, pgx_te)], by = c("Eid", "omicID"), all.x = TRUE)
setorder(g, p_train)
S11 <- g[, .(Exposure = Eid, Term = E_term, Protein = omicID, Category,
             `Genetic component` = G_component,
             `N train` = samplesize_train, `N test` = samplesize_test,
             `beta train` = signif(beta_train, 4), `SE train` = signif(se_train, 4),
             `p train` = p_train,
             `beta test` = signif(beta_test, 4), `SE test` = signif(se_test, 4),
             `p test` = p_test,
             `p GxE block train` = pgx_tr, `p GxE block test` = pgx_te)]
fwrite(S11, file.path(OUT, "gxe_replicated_terms.tsv"), sep = "\t")
message(sprintf("  gxe_replicated_terms.tsv    %d rows x %d cols (GxE term thr %.3g)",
                nrow(S11), ncol(S11), thrG))
