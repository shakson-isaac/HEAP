#!/usr/bin/env Rscript
# ============================================================================
# export_assoc_spec_deposit.R   (DEPOSIT -- Supplementary Data, associations)
# ----------------------------------------------------------------------------
# Writes the per-covariate-specification exposure->protein (E) and polygenic
# gene-by-environment (GxE) association tables that ship in
# HEAP_Supplementary_Data.zip under associations_E/ and associations_GxE/.
#
# WHY THIS EXISTS. Module 1 already deposits its full specification grid
# (varcomp_specs*.tsv, 11 specs). Module 2 ran the same sweep -- all seven
# experiments are on disk under output/module2/ -- but only BMI and clinical
# were ever aggregated, and only into the paired attenuation.tsv the robustness
# figure needs. A reader could not pull "the E associations under +blood draw".
#
# CANONICAL SOURCE RUNS -- this header is the contract:
#   E    output/module2/<experiment>/<covariate_set>/  via load_module2_replicated()
#   GxE  the same runs, statGxE, via load_module2_results()
# Both go through the shared loaders in visualizations/common/load_heap_results.R
# ON PURPOSE. Those loaders carry the unique() guard against Module2.R's
# first-protein-per-batch double-count bug (load_heap_results.R:186-189). Do NOT
# hand-roll a reader over the run_config shards -- it will silently double-count.
#
# CONSTRUCTION, matched to the already-shipped base tables so the deposit and
# the supplement cannot disagree:
#   E   = load_module2_replicated(), the same call fig_module2_spec_sensitivity
#         and \Tref{assoc_full} (S8) use. Replication is Bonferroni 0.05/n_terms
#         in BOTH train and test, with n_terms recomputed within each spec.
#   GxE = term-level statGxE merged train/test on (ID, omicID, G_component) with
#         Bonferroni 0.05/n_terms -- the exact recipe of \Tref{gxe_assoc} (S13)
#         at summarize_module2_supp_tables.R:76-95.
#
# ROW SET. Every file carries the SAME key set: the union, across all specs, of
# every pair replicating in at least one specification. That makes the files
# directly joinable and, crucially, keeps a base hit visible in the spec where
# it FAILS -- the attenuation signal (LEP/FABP4 under BMI) would be invisible if
# each file held only its own replicated rows. `replicated` is per-spec, so
# survival, loss and novelty are all readable off the flag.
#
# SPEC COVERAGE. E ships 7, GxE ships 6. M2_base_ses is deliberately ABSENT from
# GxE: base_ses calls E_to_Covariate(), deleting the Deprivation_Indices category
# and moving its 9 variables into the covariate matrix, so 16 of the 108
# base-replicated GxE pairs become untestable by construction rather than "lost
# to adjustment"; it is also mis-specified for interactions absent GxC/ExC
# controls. See module2_gxe_spec_sensitivity.R:28-40. It IS shipped for E, where
# that reasoning does not apply, labelled as an alternative estimand.
#
# Run: module load gcc/14.2.0 R/4.4.2
#      export HEAP_PATHS_FILE=/n/groups/patel/shakson_ukb/HEAP/workflow/00_paths.R
#      Rscript scripts/analysis_summaries/export_assoc_spec_deposit.R
# ============================================================================
suppressPackageStartupMessages({ library(data.table) })
source(Sys.getenv("HEAP_PATHS_FILE"))
local({ c <- "/n/groups/patel/shakson_ukb/HEAP/scripts/visualizations/common"
  for (f in c("figure_paths", "load_heap_results", "label_helpers")) source(file.path(c, paste0(f, ".R"))) })

DEP <- file.path(heap_project_output("supp_deposit"))
OUT_E <- file.path(DEP, "assoc_E_specs"); OUT_G <- file.path(DEP, "assoc_GxE_specs")
for (d in c(OUT_E, OUT_G)) { unlink(d, recursive = TRUE); dir.create(d, recursive = TRUE, showWarnings = FALSE) }
msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")

# spec table -- `file` is the reader-facing name, not the internal run id
# NOTE: prevalent disease is released as an EXCLUSION only. Conditioning on it
# would adjust for a variable downstream of exposure, opening the collider path the
# reviews raised; the exclusion answers the same question without that risk. The
# SES/deprivation estimand is likewise withheld: moving deprivation out of the
# exposome into the covariate set deletes a whole exposure domain, so it is a
# re-specification of the estimand rather than a robustness check.
SPECS <- data.table(
  exp  = c("M2_base_main", "M2_base_bmi", "M2_base_clinical_main", "M2_base_draw",
           "M2_base_exclprev"),
  cov  = c("base", "base_bmi", "base_clinical", "base_draw", "base"),
  file = c("base", "base_plus_bmi", "base_plus_clinical", "base_plus_blood_draw",
           "exclude_prevalent_disease"),
  kind = c("primary", "covariate", "covariate", "covariate", "sample"),
  gxe  = c(TRUE, TRUE, TRUE, TRUE, TRUE))

# ---- E: two passes. Pass 1 collects replicated keys only (cheap to hold), so
# every output file can share one key set; pass 2 subsets and writes.
key_of <- function(d) paste0(d$ID, "\r", d$omicID)

msg("E pass 1/2 -- collecting replicated keys across ", nrow(SPECS), " specs")
repl_keys <- list()
for (i in seq_len(nrow(SPECS))) {
  m <- load_module2_replicated(SPECS$cov[i], experiment = SPECS$exp[i])
  repl_keys[[SPECS$exp[i]]] <- key_of(m[replicated == TRUE])
  msg("  ", SPECS$file[i], ": ", nrow(m), " terms, ", length(repl_keys[[SPECS$exp[i]]]), " replicated")
  rm(m); gc(FALSE)
}
UNION <- unique(unlist(repl_keys, use.names = FALSE))
msg("E union across specs: ", length(UNION), " term x protein pairs")

msg("E pass 2/2 -- writing")
e_rows <- integer(0)
for (i in seq_len(nrow(SPECS))) {
  m <- load_module2_replicated(SPECS$cov[i], experiment = SPECS$exp[i])
  m <- m[key_of(m) %in% UNION]
  setorder(m, p_train)
  out <- m[, .(Exposure = Eid, Term = ID, Protein = omicID, Category,
               spec = SPECS$file[i], spec_kind = SPECS$kind[i], run_id = SPECS$exp[i],
               covariate_set = SPECS$cov[i],
               N_train = samplesize_train, N_test = samplesize_test,
               beta_train = signif(beta_train, 4), SE_train = signif(se_train, 4), p_train,
               beta_test  = signif(beta_test, 4),  SE_test  = signif(se_test, 4),  p_test,
               replicated)]
  fwrite(out, file.path(OUT_E, paste0(SPECS$file[i], ".tsv")), sep = "\t")
  e_rows[SPECS$file[i]] <- nrow(out)
  msg("  wrote associations_E/", SPECS$file[i], ".tsv  ", nrow(out), " rows, ",
      sum(out$replicated), " replicated here")
  rm(m, out); gc(FALSE)
}

# ---- GxE: term level, the S13 recipe --------------------------------------
gx_terms <- function(cov, exp) {
  pull <- function(split) load_module2_results(cov, split, experiment = exp)$statGxE[
    , .(ID, omicID, Eid, Category, E_term, G_component,
        beta = Estimate, se = `Std. Error`, p = `Pr(>|t|)`, samplesize)]
  tr <- pull("train"); te <- pull("test")
  k <- c("ID", "omicID", "G_component")
  g <- merge(tr, te[, c(k, "beta", "se", "p", "samplesize"), with = FALSE], by = k,
             suffixes = c("_train", "_test"))
  thr <- 0.05 / nrow(g)
  g[, replicated := is.finite(p_train) & is.finite(p_test) & p_train < thr & p_test < thr]
  g[]
}
GS <- SPECS[gxe == TRUE]
gkey <- function(d) paste0(d$ID, "\r", d$omicID, "\r", d$G_component)

msg("GxE pass 1/2 -- ", nrow(GS), " specs (M2_base_ses excluded by design)")
gk <- list()
for (i in seq_len(nrow(GS))) {
  g <- gx_terms(GS$cov[i], GS$exp[i])
  gk[[GS$exp[i]]] <- gkey(g[replicated == TRUE])
  msg("  ", GS$file[i], ": ", nrow(g), " terms, ", length(gk[[GS$exp[i]]]), " replicated")
  rm(g); gc(FALSE)
}
GUNION <- unique(unlist(gk, use.names = FALSE))
msg("GxE union across specs: ", length(GUNION), " term x protein x component")

msg("GxE pass 2/2 -- writing")
g_rows <- integer(0)
for (i in seq_len(nrow(GS))) {
  g <- gx_terms(GS$cov[i], GS$exp[i])
  g <- g[gkey(g) %in% GUNION]
  setorder(g, p_train)
  out <- g[, .(Exposure = Eid, Term = E_term, Protein = omicID, Category,
               Genetic_component = G_component,
               spec = GS$file[i], spec_kind = GS$kind[i], run_id = GS$exp[i],
               covariate_set = GS$cov[i],
               N_train = samplesize_train, N_test = samplesize_test,
               beta_train = signif(beta_train, 4), SE_train = signif(se_train, 4), p_train,
               beta_test  = signif(beta_test, 4),  SE_test  = signif(se_test, 4),  p_test,
               replicated)]
  fwrite(out, file.path(OUT_G, paste0(GS$file[i], ".tsv")), sep = "\t")
  g_rows[GS$file[i]] <- nrow(out)
  msg("  wrote associations_GxE/", GS$file[i], ".tsv  ", nrow(out), " rows, ",
      sum(out$replicated), " replicated here")
  rm(g, out); gc(FALSE)
}

# ---- per-folder manifests --------------------------------------------------
wm <- function(dir, S, n) fwrite(data.table(
  file = paste0(S$file, ".tsv"), specification = S$file, spec_kind = S$kind,
  covariate_set = S$cov, run_id = S$exp, rows = as.integer(n[S$file])),
  file.path(dir, "manifest.tsv"), sep = "\t")
wm(OUT_E, SPECS, e_rows); wm(OUT_G, GS, g_rows)
msg("done. ", OUT_E, " and ", OUT_G)
