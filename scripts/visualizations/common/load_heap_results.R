#!/usr/bin/env Rscript

# ============================================================================
# load_heap_results.R — canonical loaders for HEAP module outputs
# ----------------------------------------------------------------------------
# One loader per module. Each loader:
#   * reads the *canonical per-idx module outputs* (the source of truth),
#   * aggregates them into a tidy object for plotting,
#   * fails with an informative message naming the module/experiment to run
#     first when the upstream output is missing.
#
# These replace the legacy two-stage pattern where analy*.R scripts hand-rolled
# HEAPres/*.qs aggregate objects that *_main.R plot scripts then read. Plotting
# scripts should call these loaders instead of re-implementing aggregation.
#
# Canonical output layout (see scripts/module*/ and 00_paths.R):
#   module1/<covarType>/{R2groups,ShapMain,lassofit}_<idx>.txt   (+ *_<idx>.rds)
#   module1_predictive_r2_final/<covarType>/<method>/predictive_r2_*_<idx>.txt
#   module2/<covarType>/univar_assoc_<idx>.rds  (train/test: statE,statGxE,statR2,statFblock)
#   module3/<covarType>/<family>/<mode>/MDres_<idx>.txt
#   mr_edges/global_edges/edges_*.tsv, HEAPres.tsv, MR_priority_table.tsv
#   module6_pes_test/<covarType>/{PES_*,Cox4All_*}.tsv
#   population_architecture/summary/*.tsv
# ============================================================================

local({
  if (exists("heap_resolve_output", mode = "function")) return(invisible())
  here <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) NA)
  cand <- c(
    if (!is.na(here)) file.path(here, "figure_paths.R") else character(0),
    file.path(getwd(), "scripts", "visualizations", "common", "figure_paths.R"),
    "/n/groups/patel/shakson_ukb/HEAP/scripts/visualizations/common/figure_paths.R"
  )
  hit <- cand[file.exists(cand)][1]
  if (is.na(hit)) stop("load_heap_results.R: cannot find figure_paths.R")
  source(hit)
})

suppressPackageStartupMessages({
  library(data.table)
})

# Internal: list per-idx files matching a glob in a resolved module dir.
.heap_list_idx_files <- function(subdir, pattern, recursive = FALSE) {
  d <- heap_resolve_output(subdir, must_exist = TRUE)
  f <- list.files(d, pattern = pattern, full.names = TRUE, recursive = recursive)
  if (length(f) == 0L)
    stop("No files matching '", pattern, "' under ", d,
         "\nThe module directory exists but holds no results — re-run the module.",
         call. = FALSE)
  f
}

# ---------------------------------------------------------------------------
# Module 1 — variance decomposition (Shapley R2) + predictive R2
# ---------------------------------------------------------------------------

#' Load Module 1 LEGACY Shapley results (R2groups / ShapMain).
#'
#' DEPRECATED for figures: the CURRENT Module 1 (Module1_suggested.R) does NOT
#' produce R2groups/ShapMain — those come only from the trashed old module
#' (scripts/TRASH/Module1_oldv1.R), so `output/module1/R2groups_*` is stale
#' legacy output. The canonical Module 1 variance decomposition is the
#' predictive-R2 tables; use `load_module1_predictive_r2()` instead (the coarse
#' table's `method == "sequential_incremental"` rows are the block increments).
#' Retained only for inspecting legacy runs.
#'
#' @param covarType covariate set name, e.g. "Type5" (default) or "Type3"
#' @param what one or more of: "r2groups" (group-level Shapley R2),
#'   "shap" (per-feature Shapley R2), "lassofit" (lasso fit summaries)
#' @return named list of data.tables, each with an added `omic` protein column
load_module1_results <- function(covarType = "Type5",
                                  what = c("r2groups", "shap")) {
  what <- match.arg(what, c("r2groups", "shap", "lassofit"), several.ok = TRUE)
  sub <- file.path("module1", covarType)
  map <- c(r2groups = "^R2groups_.*\\.txt$",
           shap     = "^ShapMain_.*\\.txt$",
           lassofit = "^lassofit_.*\\.txt$")
  out <- lapply(what, function(w) {
    files <- .heap_list_idx_files(sub, map[[w]])
    rbindlist(lapply(files, fread), fill = TRUE)
  })
  names(out) <- what
  out
}

#' Load Module 1 *predictive* R2 partition tables (cross-validated OOF R2).
#'
#' Canonical (post-covariate-restructure) layout is experiment-nested:
#'   module1_predictive_r2_score_partition/<experiment>/<covarType>/<method>/
#'     predictive_r2_<level>_<idx>.txt
#' e.g. the primary base run is .../M1_base_lasso/base/lasso/. Pass `experiment`
#' to read that layout. When `experiment` is NULL the loader reads the legacy
#' flat layout <root>/<covarType>/<method>/ (kept for old Type* runs).
#'
#' METHOD COLUMN (current Module1_suggested.R, score_partition mode):
#'   coarse:               method in {score_model_total, score_unique_drop}
#'     - score_model_total  -> nested cumulative TEST R2 (block C, C+G, C+G+E, C+G+E+GxE)
#'     - score_unique_drop  -> unique block R2 = r2_full - r2_reduced
#'                             (block Covars, G, E, GxE)
#'   genetic_subblocks:    method == score_unique_drop (block Gcis, Gtrans; +present)
#'   exposure_categories:  method == score_unique_drop (per fine exposure category)
#'   gxe_categories:       method == score_unique_drop (per fine GxE category)
#'
#' @param covarType covariate-set name, e.g. "base" (default; legacy "Type5")
#' @param method    regularization family folder: "lasso" (default), "enet", "ridge"
#' @param level     which partition table: "coarse" (default), "genetic_subblocks",
#'   "exposure_categories", "gxe_categories"
#' @param experiment experiment name (e.g. "M1_base_lasso"); NULL -> legacy flat layout
#' @param final     TRUE -> module1_predictive_r2_final;
#'                  FALSE (default) -> module1_predictive_r2_score_partition
load_module1_predictive_r2 <- function(covarType = "base", method = "lasso",
                                       level = "coarse",
                                       experiment = "M1_base_lasso",
                                       final = FALSE) {
  level <- match.arg(level, c("coarse", "genetic_subblocks",
                              "exposure_categories", "gxe_categories"))
  root <- if (final) "module1_predictive_r2_final" else
                     "module1_predictive_r2_score_partition"
  sub <- if (is.null(experiment) || !nzchar(experiment))
           file.path(root, covarType, method)
         else file.path(root, experiment, covarType, method)
  files <- .heap_list_idx_files(sub, paste0("^predictive_r2_", level, "_.*\\.txt$"))
  rbindlist(lapply(files, fread), fill = TRUE)
}

#' Load Module 1 per-fold lasso fit summaries (whole-model train/test R2, lambda,
#' sparsity). One row per (protein, fold). Same experiment-nested layout as
#' load_module1_predictive_r2().
#'
#' Columns: omic, fold, family, model_label, model_class, train_r2, test_r2,
#'   alpha, lambda, cvm, n_nonzero, n_design_cols, n_train, n_test, inner_kfold.
load_module1_fit_summary <- function(covarType = "base", method = "lasso",
                                     experiment = "M1_base_lasso", final = FALSE) {
  root <- if (final) "module1_predictive_r2_final" else
                     "module1_predictive_r2_score_partition"
  sub <- if (is.null(experiment) || !nzchar(experiment))
           file.path(root, covarType, method)
         else file.path(root, experiment, covarType, method)
  files <- .heap_list_idx_files(sub, "^fit_summary_.*\\.txt$")
  rbindlist(lapply(files, fread), fill = TRUE)
}

# ---------------------------------------------------------------------------
# Module 2 — E / G / GxE univariate associations
# ---------------------------------------------------------------------------

#' Load Module 2 association results, aggregated across all proteins.
#'
#' Each univar_assoc_<idx>.rds holds list(train=list(statE,statGxE,statR2,
#' statFblock), test=list(...)). Elements are positional; this loader names them.
#'
#' Canonical (post-covariate-restructure) layout is experiment-nested:
#'   module2/<experiment>/<covariate_set>/univar_assoc_<idx>.rds
#' e.g. the primary base run is module2/M2_base_main/base/. Pass `experiment`
#' to read that layout. When `experiment` is NULL the loader reads the legacy
#' flat layout module2/<covarType>/ (kept for old Type* runs).
#'
#' @param covarType covariate-set name, e.g. "base" (default) or "base_clinical"
#'   (legacy: "Type3"/"Type5")
#' @param split "train" (default) or "test"
#' @param sens  if TRUE read module2_sens instead of module2
#' @param experiment experiment name (e.g. "M2_base_main"); NULL -> legacy flat
#'   layout module2/<covarType>/
#' @return named list with elements statE, statGxE, statR2, statFblock
#'   (data.tables, row-bound across proteins)
load_module2_results <- function(covarType = "base", split = c("train", "test"),
                                 sens = FALSE, experiment = "M2_base_main") {
  split <- match.arg(split)
  base <- if (sens) "module2_sens" else "module2"
  sub <- if (is.null(experiment) || !nzchar(experiment))
           file.path(base, covarType)
         else file.path(base, experiment, covarType)
  files <- .heap_list_idx_files(sub, "^univar_assoc_.*\\.rds$")
  comp_names <- c("statE", "statGxE", "statR2", "statFblock")
  acc <- setNames(vector("list", length(comp_names)), comp_names)
  for (f in files) {
    obj <- readRDS(f)
    part <- obj[[split]]
    if (is.null(part)) next
    for (i in seq_along(comp_names)) {
      if (length(part) >= i && !is.null(part[[i]]))
        acc[[i]][[length(acc[[i]]) + 1L]] <- as.data.table(part[[i]])
    }
  }
  # unique(): each stat row (term x protein) is unique by construction; this
  # guards against any exact-duplicate rows in older outputs (e.g. the fixed
  # first-protein double-count bug in Module2.R's batch accumulator).
  lapply(acc, function(x) if (length(x)) unique(rbindlist(x, fill = TRUE)) else data.table())
}

#' Module 2 statE merged across the 80/20 TRAIN/TEST split, with a replication
#' flag — the basis for the manuscript Figure-3 association panels.
#'
#' Merges the per-TERM E main-effect tables (statE) from the train and test
#' splits on (ID, omicID), suffixing stats `_train`/`_test`, and flags each
#' exposure-term x protein pair as `replicated` when it passes Bonferroni in
#' BOTH splits (the legacy "significant in train AND test" criterion; the
#' replication filter also naturally drops the ordered-factor high-degree
#' polynomial-contrast terms — e.g. deprivation indices — whose raw coefficients
#' explode but do not reproduce across splits). Bonferroni threshold is
#' 0.05 / (number of merged term x protein pairs), matching Module2's legacy
#' miami plot and the manuscript's p<7e-8.
#'
#' Each row is one (exposure TERM, protein). `ID` is the model term name
#' (treatment-coded factor levels keep an interpretable per-level coefficient,
#' e.g. alcohol_intake_frequency_f1558_0_06); `Eid` is the base exposure;
#' `Category` is the fine exposure category. `spec` labels the source
#' experiment so several specs can be row-bound and pooled (Fig 3D/3E pool the
#' effect sizes "across covariate specifications").
#'
#' @param covarType covariate-set name (default "base")
#' @param experiment experiment name (default "M2_base_main"); NULL -> legacy flat
#' @param sens read module2_sens instead of module2
#' @return data.table with columns ID, omicID, Eid, Category, spec, AssocID,
#'   beta_train, se_train, t_train, p_train, beta_test, se_test, t_test, p_test,
#'   samplesize_train/test, replicated; attr "pval_thresh" holds the Bonferroni
#'   threshold used.
load_module2_replicated <- function(covarType = "base",
                                    experiment = "M2_base_main", sens = FALSE) {
  pull <- function(split) {
    s <- load_module2_results(covarType, split, sens, experiment)$statE
    if (!nrow(s)) stop("Empty statE (", split, ") for ", experiment, "/", covarType)
    s[, .(ID, omicID, Eid, Category,
          beta = Estimate, se = `Std. Error`, tval = `t value`,
          p = `Pr(>|t|)`, samplesize)]
  }
  tr <- pull("train"); te <- pull("test")
  trn <- tr[, .(ID, omicID, Eid, Category,
                beta_train = beta, se_train = se, t_train = tval,
                p_train = p, samplesize_train = samplesize)]
  ten <- te[, .(ID, omicID,
                beta_test = beta, se_test = se, t_test = tval,
                p_test = p, samplesize_test = samplesize)]
  m <- merge(trn, ten, by = c("ID", "omicID"))
  m[, spec := experiment]
  m[, AssocID := paste0(omicID, ":", ID)]
  thr <- 0.05 / nrow(m)
  m[, replicated := is.finite(p_train) & is.finite(p_test) &
                    p_train < thr & p_test < thr]
  setattr(m, "pval_thresh", thr)
  m[]
}

# ---------------------------------------------------------------------------
# Module 3 — mediation (GEM / disease-specific)
# ---------------------------------------------------------------------------

#' Load Module 3 mediation result tables (MDres_<idx>.txt), row-bound.
#'
#' @param covarType e.g. "Type3"
#' @param family    family_type folder, e.g. "lasso"
#' @param mode      mediation mode folder: "primary_total" (default) or
#'                  "partitioned_categories"
#' Default experiment-name for a Module 3 mediation mode.
#'
#' Canonical post-restructure runs are experiment-nested as
#'   module3/<experiment>/<covarType>/<family>/<mode>/MDres_*.txt
#' where the experiment encodes the base run + mode, e.g.
#'   primary_total                  -> M3_base_lasso_primary
#'   partitioned_categories         -> M3_base_lasso_partitioned
#'   partitioned_grouped_categories -> M3_base_lasso_grouped
#' Override by passing `experiment` explicitly to load_module3_results().
heap_module3_default_experiment <- function(mode = "primary_total",
                                            covarType = "base",
                                            family = "lasso") {
  short <- c(primary_total = "primary",
             partitioned_categories = "partitioned",
             partitioned_grouped_categories = "grouped")[mode]
  if (is.na(short)) return(NA_character_)
  sprintf("M3_%s_%s_%s", covarType, family, short)
}

#' Load Module 3 generalized-mediation (GEM) results, row-bound across proteins.
#'
#' Canonical (post-covariate-restructure) layout is experiment-nested:
#'   module3/<experiment>/<covarType>/<family>/<mode>/MDres_<idx>.txt
#' e.g. the primary base run is module3/M3_base_lasso_primary/base/lasso/
#' primary_total/. When `experiment` is NULL the loader derives the default
#' experiment from `mode` (see heap_module3_default_experiment); if that path is
#' absent it falls back to the legacy flat layout module3/<covarType>/<family>/
#' <mode>/ (kept for old Type* pilots).
#'
#' Each MDres_<idx>.txt row is one (protID, DZ_ID, predictor, effect_type) record.
#' effect_type is "NDE" (direct) or "NIE" (mediated/indirect); predictor_class is
#' genetic_total/exposure_total (primary) or genetic_cis/genetic_trans/
#' exposure_category (partitioned). Key columns: effect_logHR, effect_HR,
#' delta_se/l95/u95/p (delta-method CI + p on the log-HR scale),
#' delta_l95_HR/delta_u95_HR (HR-scale CI), protein_HR (protein->disease Cox HR),
#' n, n_cases, instrument_present (FALSE = structural-0 genetic component).
#'
#' @param covarType covariate-set name, e.g. "base" (default; legacy "Type3"/"Type5")
#' @param family    regularization family folder: "lasso" (default), "enet", "ridge"
#' @param mode      "primary_total" (default), "partitioned_categories",
#'   "partitioned_grouped_categories"
#' @param experiment experiment name; NULL -> derive from mode then fall back to flat
#' @param select    optional character vector of columns to read (passed to fread);
#'   use to keep memory down on the large partitioned tables (~12M rows)
#' @return data.table row-bound across all per-idx MDres files
load_module3_results <- function(covarType = "base", family = "lasso",
                                 mode = "primary_total",
                                 experiment = NULL, select = NULL) {
  flat <- file.path("module3", covarType, family, mode)
  sub  <- flat
  if (is.null(experiment)) {
    exp_default <- heap_module3_default_experiment(mode, covarType, family)
    if (!is.na(exp_default)) {
      nested <- file.path("module3", exp_default, covarType, family, mode)
      if (heap_output_exists(nested)) sub <- nested
    }
  } else if (nzchar(experiment)) {
    sub <- file.path("module3", experiment, covarType, family, mode)
  }
  files <- .heap_list_idx_files(sub, "^MDres_.*\\.txt$")
  dt <- rbindlist(lapply(files, function(f) fread(f, select = select)), fill = TRUE)
  # Defensive: a partially-written MDres file (the Module 3 run may still be in
  # progress) can make a numeric column parse as character in one file, which
  # rbindlist then coerces across the WHOLE column -> is.finite() blanks every
  # row. Coerce the known numeric columns back to numeric (bad tokens -> NA, i.e.
  # the malformed rows drop out) so one in-flight file can't sink the figure.
  num_cols <- c("effect_logHR", "effect_HR", "delta_se", "delta_l95", "delta_u95",
                "delta_p", "delta_l95_HR", "delta_u95_HR", "contrast_sd_total",
                "protein_HR", "protein_HR_l95", "protein_HR_u95", "protein_p",
                "mediator_adjR2", "cox_cindex", "cox_cindex_se", "n", "n_cases")
  for (cc in intersect(num_cols, names(dt)))
    if (!is.numeric(dt[[cc]]))
      suppressWarnings(dt[, (cc) := as.numeric(get(cc))])
  dt[]
}

#' Add proportion-mediated columns to a Module 3 table.
#'
#' Reshapes the long (NDE, NIE) rows to one row per (protID, DZ_ID, predictor)
#' and computes the mediation decomposition on the log-HR (additive) scale:
#'   TE_logHR  = NDE_logHR + NIE_logHR        (total effect of the driver)
#'   PM        = NIE_logHR / TE_logHR         (proportion mediated through protein)
#' PM is the standard VanderWeele proportion-mediated. It is only interpretable
#' when NDE and NIE share sign (0<=PM<=1); we flag the rest as `pm_consistent`
#' (FALSE = inconsistent mediation / suppression, where PM can be <0 or >1) and
#' clamp a reported `pm_display` to [0,1] for plotting while keeping raw `pm`.
#' A tiny |TE| guard avoids division blow-ups. NIE significance carries the
#' delta-method p (`NIE_p`) AND a multiple-testing-adjusted q (`NIE_q`): figures
#' should threshold on `NIE_q`, not the raw p. The q is computed over EVERY NIE
#' test in the passed table (one per protein x disease x predictor), so pass the
#' FULL mode table here, then filter — that makes the adjustment family the whole
#' mediation screen for that analysis. Structural-0 genetic rows have NA p and are
#' (correctly) excluded from the family.
#'
#' @param md a data.table from load_module3_results() (must hold both NDE & NIE)
#' @param te_floor minimum |TE_logHR| for a defined PM (default 1e-4)
#' @param padjust p.adjust method for NIE_q over the NIE family (default "BH" = FDR)
#' @return data.table, one row per (protID, DZ_ID, predictor, predictor_class)
heap_proportion_mediated <- function(md, te_floor = 1e-4, padjust = "BH") {
  stopifnot(all(c("protID","DZ_ID","predictor","effect_type","effect_logHR") %in% names(md)))
  keep <- c("protID","DZ_ID","predictor","predictor_class","effect_type",
            "effect_logHR","effect_HR","delta_p","instrument_present",
            "n","n_cases","protein_HR","protein_p")
  keep <- intersect(keep, names(md))
  long <- md[effect_type %in% c("NDE","NIE"), ..keep]
  idv  <- intersect(c("protID","DZ_ID","predictor","predictor_class",
                      "instrument_present","n","n_cases","protein_HR","protein_p"),
                    names(long))
  w <- dcast(long, as.formula(paste(paste(idv, collapse = " + "),
                                    "~ effect_type")),
             value.var = c("effect_logHR","effect_HR","delta_p"))
  setnames(w,
    old = c("effect_logHR_NDE","effect_logHR_NIE","effect_HR_NDE","effect_HR_NIE",
            "delta_p_NDE","delta_p_NIE"),
    new = c("NDE_logHR","NIE_logHR","NDE_HR","NIE_HR","NDE_p","NIE_p"),
    skip_absent = TRUE)
  w[, TE_logHR := NDE_logHR + NIE_logHR]
  w[, TE_HR    := exp(TE_logHR)]
  w[, pm := fifelse(abs(TE_logHR) >= te_floor, NIE_logHR / TE_logHR, NA_real_)]
  # PM is well-behaved only when direct and indirect effects agree in sign.
  w[, pm_consistent := is.finite(pm) & sign(NDE_logHR) == sign(NIE_logHR) &
                        pm >= 0 & pm <= 1]
  w[, pm_display := pmin(pmax(pm, 0), 1)]
  # Multiple-testing-adjusted NIE significance over the whole NIE family.
  w[, NIE_q := NA_real_]
  if ("NIE_p" %in% names(w))
    w[is.finite(NIE_p), NIE_q := stats::p.adjust(NIE_p, padjust)]
  w[]
}

#' Add a multiple-testing-adjusted q column to a long Module 3 table.
#'
#' For figures that threshold delta_p directly (not via heap_proportion_mediated):
#' computes p.adjust over every FINITE p in `dt[p_col]` and writes `q_col`. Pass
#' the rows that form the test family (e.g. ALL NIE rows for a mode) BEFORE
#' filtering to the significant/plotted subset, so the family is the full screen.
#'
#' @param dt data.table (modified by reference; also returned)
#' @param p_col source p-value column (default "delta_p")
#' @param q_col destination column (default "delta_q")
#' @param method p.adjust method (default "BH" = Benjamini-Hochberg FDR)
heap_md_fdr <- function(dt, p_col = "delta_p", q_col = "delta_q", method = "BH") {
  stopifnot(p_col %in% names(dt))
  dt[, (q_col) := NA_real_]
  dt[is.finite(get(p_col)), (q_col) := stats::p.adjust(get(p_col), method)]
  dt[]
}

# ---------------------------------------------------------------------------
# Module 5 — Mendelian randomization edge tables
# ---------------------------------------------------------------------------

#' Load Module 5 MR global edge tables.
#'
#' @param which one or more edge-table stems (without .tsv), e.g.
#'   "edges_PD", "edges_EP", "MR_priority_table", "HEAPres". Default loads all
#'   *.tsv in mr_edges/global_edges.
#' @return named list of data.tables keyed by file stem
load_module5_results <- function(which = NULL) {
  d <- heap_resolve_output(file.path("mr_edges", "global_edges"), must_exist = TRUE)
  files <- list.files(d, pattern = "\\.tsv$", full.names = TRUE)
  if (length(files) == 0L)
    stop("No edge .tsv files under ", d, " — run Module 5 (MR).", call. = FALSE)
  stems <- tools::file_path_sans_ext(basename(files))
  if (!is.null(which)) {
    keep <- stems %in% which
    if (!any(keep))
      stop("None of requested tables found. Available: ",
           paste(stems, collapse = ", "), call. = FALSE)
    files <- files[keep]; stems <- stems[keep]
  }
  setNames(lapply(files, fread), stems)
}

# ---------------------------------------------------------------------------
# Module 5 — MR summary/sensitivity/motif tables (support/mr_tables aggregation)
#
# These read the canonical per-arm aggregates written by
# scripts/support/mr_tables/build_mr_tables.R (NOT the global_edges triad lists):
#   mr_edges/summary/<which>.tsv            (UKB primary arm)
#   mr_edges/summary/DECODE/<which>.tsv     (deCODE replication arm)
# where <which> is one of: MRmotifs, mr_sensitivity_long, sensitivity_by_edgedir,
# sensitivity_overall (+ arm_comparison_* written by compare_arms.R).
# ---------------------------------------------------------------------------

#' Directory holding the MR summary tables for one instrument arm.
heap_mr_summary_dir <- function(cohort = "UKB") {
  base <- heap_resolve_output(file.path("mr_edges", "summary"), must_exist = FALSE)
  if (toupper(cohort) == "DECODE") file.path(base, "DECODE") else base
}

#' Load one MR summary table for an arm.
#' @param which file stem, e.g. "MRmotifs", "mr_sensitivity_long",
#'   "sensitivity_by_edgedir", "sensitivity_overall", "arm_comparison_motif_counts"
#' @param cohort "UKB" | "DECODE"
load_mr_table <- function(which, cohort = "UKB") {
  f <- file.path(heap_mr_summary_dir(cohort), paste0(which, ".tsv"))
  if (!file.exists(f))
    stop("Missing MR table: ", f, "\nRun: Rscript ",
         "scripts/support/mr_tables/build_mr_tables.R ", cohort, call. = FALSE)
  fread(f)
}

#' Load an MR summary table for BOTH arms, row-bound (tagged by `dataset`).
#' Returns whatever arms are present (NULL-skips a missing arm).
load_mr_table_both <- function(which) {
  out <- rbindlist(lapply(c("UKB", "DECODE"), function(co) {
    f <- file.path(heap_mr_summary_dir(co), paste0(which, ".tsv"))
    if (!file.exists(f)) return(NULL)
    dt <- fread(f)
    if (!"dataset" %in% names(dt)) dt[, dataset := co]
    dt[]
  }), fill = TRUE)
  if (!nrow(out))
    stop("No MR '", which, "' tables found for either arm. Run ",
         "scripts/support/mr_tables/build_mr_tables.R {UKB,DECODE}.", call. = FALSE)
  out
}

#' Load the systematic colocalization results (one row per cis-pQTL x outcome
#' locus, both instrument arms) written by
#' scripts/support/coloc/run_coloc_systematic.R.
#'
#' Canonical table: support/coloc/coloc_results.tsv with columns
#'   arm, protID, target, edge_dir, lead_snp, chr, pos, nsnps, PP.H3, PP.H4, status.
#' `target` is a FinnGen disease code for Pcis_to_D edges and a UKB exposure
#' field for Pcis_to_E edges. Adds a `colocalized` flag at the canonical PP.H4
#' threshold (the same gate build_mr_tables.R applies to demote LD-confounded
#' Tier-1 cis edges).
#'
#' @param pp_h4_thresh canonical colocalization cutoff (default 0.8).
#' @return data.table of coloc loci with a logical `colocalized` column.
load_coloc_results <- function(pp_h4_thresh = 0.8) {
  f <- file.path(heap_resolve_output(file.path("support", "coloc"),
                                     must_exist = FALSE), "coloc_results.tsv")
  if (!file.exists(f))
    stop("Missing coloc results: ", f, "\nRun: Rscript ",
         "scripts/support/coloc/run_coloc_systematic.R", call. = FALSE)
  dt <- fread(f)
  if (!"PP.H4" %in% names(dt))
    stop("coloc_results.tsv lacks a PP.H4 column: ", f, call. = FALSE)
  dt[, PP.H4 := as.numeric(`PP.H4`)]
  dt[, colocalized := is.finite(`PP.H4`) & `PP.H4` >= pp_h4_thresh]
  dt[]
}

# ---------------------------------------------------------------------------
# Module 6 — PES prediction + longitudinal Cox validation
# ---------------------------------------------------------------------------

#' Load Module 6 PES prediction / Cox validation tables for a covariate set.
#'
#' @param covarType e.g. "Type5"
#' @param pattern   filename glob to collect (default Cox summary tables).
#'   Use "^PES_.*overall.*\\.tsv$" etc. for other artifacts.
#' @return data.table row-bound across proteins/exposures (file stem in `source_file`)
load_module6_results <- function(covarType = "Type5",
                                 pattern = "^Cox4All_.*\\.tsv$") {
  sub <- file.path("module6_pes_test", covarType)
  files <- .heap_list_idx_files(sub, pattern, recursive = TRUE)
  rbindlist(lapply(files, function(f) {
    dt <- fread(f); dt[, source_file := basename(f)][]
  }), fill = TRUE)
}

# ---------------------------------------------------------------------------
# Module 6 — longitudinal PES (prod sub-workflow): training + repeat-visit
# holdout + within-person change metric tables.
#
# Module6_prod_longitudinal.R writes, per exposure/level, a family of metric
# TSVs under module6_pes_longitudinal/<covarType>/ with the stem
#   PESlong_<covarType>_<exposure_id>_<Suffix>.tsv
# These are ALREADY-aggregated metric tables (cross-validated OOF performance,
# repeat-visit generalization, within-person change) — the figures consume them
# directly. The heavy per-person artefacts (FinalModelArtifact / HoldoutScores /
# TrainOOF .rds) are intentionally NOT read here.
#
# `model` levels across the tables: prot_only (proteome-only PES), cov_only
# (covariate baseline), prot_plus_cov (proteome + covariates).
# ---------------------------------------------------------------------------

# Internal: parse a Module 6 exposure_id into (base_variable, level) and join
# the canonical fine category. Categorical exposures arrive one row per level
# with the level appended after the UKB field suffix, e.g.
#   alcohol_drinker_status_f20117_0_0_Current -> base=..._f20117_0_0 level=Current
# Continuous / non-field ids (no2_mean, pm2_5_mean, *.multi_*) keep base == id.
.heap_m6_exposure_meta <- function(exposure_id) {
  eid  <- as.character(exposure_id)
  base <- sub("(_f[0-9]+_[0-9]+_[0-9]+)_.+$", "\\1", eid)   # strip trailing _Level
  level <- ifelse(base == eid, NA_character_,
                  sub("^.*_f[0-9]+_[0-9]+_[0-9]+_", "", eid))
  meta <- data.table(exposure_id = eid, base_variable = base, level = level)
  # category: prefer the curated exposure_labels.tsv, fall back to the analysis
  # manifest. Both config tables key categorical exposures by their FULL one-hot
  # id (e.g. alcohol_drinker_status_f20117_0_0_Current), so join on exposure_id
  # directly; fall back to the stripped base for any continuous id that differs.
  cat_lut <- NULL
  lf <- tryCatch(heap_config("exposure_sets", "exposure_labels.tsv"),
                 error = function(e) NA_character_)
  if (!is.na(lf) && file.exists(lf)) {
    lm <- fread(lf)
    if (all(c("variable", "category") %in% names(lm)))
      cat_lut <- setNames(lm$category, lm$variable)
  }
  if (is.null(cat_lut)) {
    af <- tryCatch(heap_config("exposure_sets", "analysis_exposures.tsv"),
                   error = function(e) NA_character_)
    if (!is.na(af) && file.exists(af)) {
      am <- fread(af)
      if (all(c("variable", "category") %in% names(am)))
        cat_lut <- setNames(am$category, am$variable)
    }
  }
  if (is.null(cat_lut)) {
    meta[, category := NA_character_]
  } else {
    cat_by_id   <- unname(cat_lut[eid])
    cat_by_base <- unname(cat_lut[base])
    meta[, category := fifelse(!is.na(cat_by_id), cat_by_id, cat_by_base)]
  }
  meta[]
}

#' Load Module 6 longitudinal PES metric tables (prod sub-workflow).
#'
#' Row-binds the per-exposure metric TSVs under
#' module6_pes_longitudinal/<covarType>/ and joins canonical exposure metadata
#' (base_variable, level, category). The `exposure_id` is taken from the file
#' stem (reliable across every table, including SelectedProteins which has no
#' exposure column).
#'
#' @param covarType covariate-set subdir (default "base")
#' @param which one or more of:
#'   "overall"    -> *_TrainOverallMetrics.tsv   (CV OOF r2/correlation/rmse per model)
#'   "fold"       -> *_TrainFoldMetrics.tsv       (per-fold r2/corr for prot/cov/full)
#'   "holdout"    -> *_HoldoutByVisitMetrics.tsv  (repeat-visit generalization)
#'   "within"     -> *_WithinPersonChangeMetrics.tsv (within-person change tracking)
#'   "transition" -> *_WithinPersonTransitionSummary.tsv (categorical level moves)
#'   "selected"   -> *_SelectedProteins.tsv       (panel size + selected protein list)
#' @return a single data.table when `which` is length 1, else a named list of
#'   data.tables. Each has added columns exposure_id, base_variable, level,
#'   category (and source_file).
load_module6_pes_longitudinal <- function(covarType = "base",
                                          which = "overall") {
  suffix_map <- c(
    overall    = "TrainOverallMetrics",
    fold       = "TrainFoldMetrics",
    holdout    = "HoldoutByVisitMetrics",
    within     = "WithinPersonChangeMetrics",
    transition = "WithinPersonTransitionSummary",
    selected   = "SelectedProteins")
  which <- match.arg(which, names(suffix_map), several.ok = TRUE)
  sub   <- file.path("module6_pes_longitudinal", covarType)
  stem_re <- paste0("^PESlong_", covarType, "_(.*)_%s\\.tsv$")

  read_one <- function(w) {
    suf <- suffix_map[[w]]
    pat <- sprintf("^PESlong_%s_.*_%s\\.tsv$", covarType, suf)
    files <- .heap_list_idx_files(sub, pat, recursive = FALSE)
    dt <- rbindlist(lapply(files, function(f) {
      d <- fread(f)
      eid <- sub(sprintf(stem_re, suf), "\\1", basename(f))
      d[, exposure_id := eid]
      d[, source_file := basename(f)]
      d[]
    }), fill = TRUE, use.names = TRUE)
    # join exposure metadata once per distinct id (cheap, avoids per-row regex)
    meta <- .heap_m6_exposure_meta(unique(dt$exposure_id))
    dt <- merge(dt, meta, by = "exposure_id", all.x = TRUE, sort = FALSE)
    dt[]
  }

  out <- lapply(which, read_one)
  names(out) <- which
  if (length(out) == 1L) out[[1]] else out
}

# ---------------------------------------------------------------------------
# Module 6 — frozen-risk Cox validation (frozenrisk sub-workflow).
#
# Module6_longitudinal_cox_frozenrisk.R freezes each trained PES and fits Cox
# disease-risk models at landmark visits, writing per-exposure tables into the
# SAME directory as the prod-longitudinal output (module6_pes_longitudinal/
# <covarType>/) with stems:
#   PESlong_<covarType>_<exposure>_FrozenCox{TrainFits,Eval,EvalLong,PersonDelta}.tsv
#   PESlong_<covarType>_<exposure>_RiskCalculator{Metrics,Deciles}.tsv
# (an optional --output-tag inserts a tag after "FrozenCox" / before the
# RiskCalculator suffix; default is empty.) Each row is keyed by
# (disease_age_col, model/landmark) — the disease outcome is `disease_age_col`.
#
# These figures are DEFERRED until the frozen-risk array runs; the loader and
# its plotters are in place so they activate automatically once the tables land.
# ---------------------------------------------------------------------------

#' Load Module 6 frozen-risk Cox validation tables.
#'
#' @param covarType covariate-set subdir (default "base")
#' @param which one or more of:
#'   "eval"         -> *_FrozenCox<tag>Eval.tsv      (WIDE: cindex_M0..M3 + delta_cindex_*)
#'   "eval_long"    -> *_FrozenCox<tag>EvalLong.tsv  (long per model/landmark)
#'   "fits"         -> *_FrozenCox<tag>TrainFits.tsv (HR_PES_perSD, p_PES per model)
#'   "person_delta" -> *_FrozenCox<tag>PersonDelta.tsv
#'   "risk_metrics" -> *_RiskCalculator<tag>Metrics.tsv (brier, calibration_slope, ...)
#'   "risk_deciles" -> *_RiskCalculator<tag>Deciles.tsv (observed vs predicted by decile)
#' @param tag optional output-tag inserted in the file stem (default "" = none)
#' @return a single data.table when `which` is length 1, else a named list. Each
#'   has added exposure_id, base_variable, level, category (and source_file).
load_module6_frozenrisk <- function(covarType = "base", which = "eval", tag = "") {
  suffix_map <- c(
    fits         = paste0("FrozenCox", tag, "TrainFits"),
    eval         = paste0("FrozenCox", tag, "Eval"),
    eval_long    = paste0("FrozenCox", tag, "EvalLong"),
    person_delta = paste0("FrozenCox", tag, "PersonDelta"),
    # RiskCalculator files put the tag at the END (RiskCalculatorMetrics_<tag>), unlike
    # the FrozenCox files which put it in the middle — match the writer in
    # Module6_longitudinal_cox_frozenrisk.R (risk_tag suffix). Empty tag is unaffected.
    risk_metrics = paste0("RiskCalculatorMetrics", if (nzchar(tag)) paste0("_", tag) else ""),
    risk_deciles = paste0("RiskCalculatorDeciles", if (nzchar(tag)) paste0("_", tag) else ""))
  which <- match.arg(which, names(suffix_map), several.ok = TRUE)
  sub   <- file.path("module6_pes_longitudinal", covarType)
  stem_re <- paste0("^PESlong_", covarType, "_(.*)_%s\\.tsv$")

  read_one <- function(w) {
    suf <- suffix_map[[w]]
    pat <- sprintf("^PESlong_%s_.*_%s\\.tsv$", covarType, suf)
    files <- .heap_list_idx_files(sub, pat, recursive = FALSE)
    dt <- rbindlist(lapply(files, function(f) {
      d <- fread(f)
      d[, exposure_id := sub(sprintf(stem_re, suf), "\\1", basename(f))]
      d[, source_file := basename(f)]
      d[]
    }), fill = TRUE, use.names = TRUE)
    # Drop censoring/administrative columns that the Cox disease selection (grep ^age_)
    # historically captured as fake "diseases" (age_of_lastfollowup, age_of_removal_0_0);
    # they are degenerate, not real outcomes. Harmless once the Cox script also excludes them.
    if ("disease_age_col" %in% names(dt))
      dt <- dt[!disease_age_col %chin% c("age_of_lastfollowup", "age_of_removal_0_0")]
    meta <- .heap_m6_exposure_meta(unique(dt$exposure_id))
    merge(dt, meta, by = "exposure_id", all.x = TRUE, sort = FALSE)[]
  }

  out <- lapply(which, read_one)
  names(out) <- which
  if (length(out) == 1L) out[[1]] else out
}

# ---------------------------------------------------------------------------
# Population architecture — REML/HE variance components
# ---------------------------------------------------------------------------

#' Load population-architecture summary tables (REML/HE variance components).
#'
#' @param pattern filename glob under population_architecture/summary
#' @return named list of data.tables keyed by file stem
load_population_architecture_results <- function(pattern = "\\.tsv$") {
  # Summaries may be written directly under the output root or a summary/ subdir.
  root <- heap_resolve_output("population_architecture", must_exist = TRUE)
  files <- list.files(root, pattern = pattern, full.names = TRUE, recursive = TRUE)
  if (length(files) == 0L)
    stop("No summary tables under ", root,
         " — run scripts/population_architecture (run + summarize).", call. = FALSE)
  setNames(lapply(files, fread), tools::file_path_sans_ext(basename(files)))
}

# ---------------------------------------------------------------------------
# Exposure GWAS (REGENIE step 2) — two-sample MR summary statistics
#
# The gwas_regenie stage writes one all-chromosome REGENIE summary per exposure
# at gwas/regenie_step2/<exposure>/<exposure>.regenie (legacy doubled name
# regenie_step2_<exp>_<exp>.regenie still read as a fallback; see io_map.yml).
# Columns: CHROM GENPOS ID ALLELE0 ALLELE1 A1FREQ INFO N TEST BETA SE CHISQ
#          LOG10P EXTRA  (LOG10P = -log10 p; ALLELE1 = effect/tested allele).
# These loaders + small QC helpers back the GWAS Manhattan / QQ / summary
# figures and reuse the same heap_resolve_output() path discipline as the rest.
# ---------------------------------------------------------------------------

# Internal: resolve the .regenie summary file for one exposure (canonical name
# first, then the legacy doubled name). Returns "" when neither exists.
.exposure_gwas_file <- function(exposure, root = NULL) {
  if (is.null(root))
    root <- heap_resolve_output(file.path("gwas", "regenie_step2"), must_exist = TRUE)
  cand <- c(
    file.path(root, exposure, paste0(exposure, ".regenie")),
    file.path(root, exposure, paste0("regenie_step2_", exposure, "_", exposure, ".regenie"))
  )
  hit <- cand[file.exists(cand)]
  if (!length(hit)) "" else hit[1]
}

#' List exposures that have a (completed) REGENIE step-2 summary statistic file.
#'
#' @param completed_only TRUE (default) keeps only exposures whose .regenie file
#'   exists; FALSE returns every exposure sub-directory under gwas/regenie_step2.
#' @return sorted character vector of exposure ids (the directory names)
list_exposure_gwas <- function(completed_only = TRUE) {
  root <- heap_resolve_output(file.path("gwas", "regenie_step2"), must_exist = TRUE)
  exps <- list.dirs(root, full.names = FALSE, recursive = FALSE)
  if (!length(exps)) return(character(0))
  if (!completed_only) return(sort(exps))
  has <- vapply(exps, function(e) nzchar(.exposure_gwas_file(e, root)), logical(1))
  sort(exps[has])
}

#' Load one exposure's REGENIE step-2 GWAS summary statistics.
#'
#' Reads only the columns needed for plotting / QC by default (the file is ~7.8M
#' variants, ~750 MB; a full read of all columns is rarely needed). Chromosome
#' and position columns are renamed to `chr`/`pos`; everything else keeps its
#' REGENIE name. Optionally pre-filters to LOG10P >= `min_log10p` so callers that
#' only need the signal (e.g. lead-variant tables) avoid materializing the cloud.
#'
#' @param exposure exposure id (a directory name under gwas/regenie_step2)
#' @param select REGENIE columns to read
#' @param min_log10p optional minimum -log10(p) filter applied during load
#' @return data.table with at least chr,pos,ID,LOG10P (+ requested columns)
load_exposure_gwas <- function(exposure,
                               select = c("CHROM", "GENPOS", "ID", "ALLELE0", "ALLELE1",
                                          "A1FREQ", "N", "BETA", "SE", "CHISQ", "LOG10P"),
                               min_log10p = NULL) {
  root <- heap_resolve_output(file.path("gwas", "regenie_step2"), must_exist = TRUE)
  f <- .exposure_gwas_file(exposure, root)
  if (!nzchar(f))
    stop("No REGENIE summary stats for exposure '", exposure, "' under ", root,
         "\nRun the gwas_regenie stage for this exposure first ",
         "(see config/io_map.yml: gwas_regenie).", call. = FALSE)
  DT <- fread(f, select = select, showProgress = FALSE)
  if ("CHROM"  %in% names(DT)) setnames(DT, "CHROM", "chr")
  if ("GENPOS" %in% names(DT)) setnames(DT, "GENPOS", "pos")
  if (!"LOG10P" %in% names(DT))
    stop("REGENIE file lacks a LOG10P column (unexpected format): ", f, call. = FALSE)
  if (!is.null(min_log10p)) DT <- DT[LOG10P >= min_log10p]
  DT[]
}

#' Genomic-inflation factor (lambda_GC) from a GWAS summary data.table.
#'
#' Uses the REGENIE CHISQ column when present, else reconstructs the 1-df chi-square
#' from LOG10P. lambda_GC = median(chi^2) / qchisq(0.5, 1).
#'
#' @param dt data.table from load_exposure_gwas (needs CHISQ or LOG10P)
#' @return single numeric lambda_GC
heap_gwas_lambda <- function(dt) {
  chisq <- if ("CHISQ" %in% names(dt)) dt[["CHISQ"]]
           else qchisq(-dt[["LOG10P"]] * log(10), df = 1, lower.tail = FALSE, log.p = TRUE)
  unname(median(chisq, na.rm = TRUE) / qchisq(0.5, df = 1))
}

#' Distance-based lead-variant clumping of genome-wide-significant variants.
#'
#' Greedy LD-free clumping: take the most significant variant, exclude every
#' variant within +/- `window` bp on the same chromosome, repeat. No LD panel
#' needed — adequate for marking independent loci on a Manhattan plot.
#'
#' @param dt data.table with chr,pos,LOG10P (from load_exposure_gwas)
#' @param log10p_thresh significance threshold on -log10(p) (default 5e-8)
#' @param window half-window in bp for the exclusion radius (default 500 kb)
#' @return data.table of lead variants (one per locus), ordered by significance
heap_gwas_lead_variants <- function(dt, log10p_thresh = -log10(5e-8), window = 5e5) {
  sig <- dt[LOG10P >= log10p_thresh][order(-LOG10P)]
  if (!nrow(sig)) return(sig)
  pos <- sig$pos; chr <- sig$chr
  used <- logical(nrow(sig)); leads <- integer(0)
  repeat {
    k <- which(!used)[1]
    if (is.na(k)) break
    leads <- c(leads, k)
    used[chr == chr[k] & abs(pos - pos[k]) <= window] <- TRUE
  }
  sig[leads]
}

#' Thin a GWAS summary table for plotting.
#'
#' Keeps every variant at or above `keep_full` (so no real signal is dropped) and
#' grid-deduplicates the dense sub-threshold cloud: at most one point per
#' (chromosome, `pos_bin`-bp column, `y_bin`-unit -log10p row). This turns ~7.8M
#' points into a few hundred thousand that render the same envelope at a fraction
#' of the file size.
#'
#' @param dt data.table with chr,pos,LOG10P
#' @param keep_full -log10(p) above which every variant is retained (default 4)
#' @param pos_bin position bin width in bp for the cloud grid (default 200 kb)
#' @param y_bin -log10(p) bin height for the cloud grid (default 0.1)
#' @return thinned data.table (a row subset of `dt`)
heap_gwas_thin <- function(dt, keep_full = 4, pos_bin = 2e5, y_bin = 0.1) {
  full <- dt[LOG10P >= keep_full]
  lo   <- dt[LOG10P <  keep_full]
  if (nrow(lo)) {
    # Group by a (chr, pos-bin, y-bin) grid and keep one representative per cell.
    # Use .I row indices rather than .SD[1L]: data.table drops columns referenced
    # in `by` expressions (here pos, LOG10P) from .SD, which would lose them.
    lo <- copy(lo)
    lo[, `:=`(.gx = floor(pos / pos_bin), .gy = round(LOG10P / y_bin))]
    lo <- lo[lo[, .I[1L], by = .(chr, .gx, .gy)]$V1]
    lo[, c(".gx", ".gy") := NULL]
  }
  rbind(full, lo[, names(full), with = FALSE])
}

# ===========================================================================
# LDSC (LD Score Regression) loaders — SNP heritability + genetic correlation
# of the exposure GWAS. Read the collected summary tables produced by
# scripts/ldsc/collect_ldsc_h2.R and scripts/ldsc/collect_ldsc_rg.R (the
# parsing of raw ldsc.py logs lives there, NOT in the plotters). See
# config/io_map.yml: ldsc and the memory project_ldsc_workflow.
# ===========================================================================

#' Locate an LDSC summary file, preferring the canonical output tree then the
#' staged figures/data copy (so a plotter still works if only the latter exists).
.ldsc_summary_file <- function(filename) {
  cand <- c(
    tryCatch(file.path(heap_resolve_output(file.path("gwas", "ldsc"),
                                           must_exist = FALSE), filename),
             error = function(e) NA_character_),
    tryCatch(heap_figure_data_path(filename), error = function(e) NA_character_)
  )
  hit <- cand[!is.na(cand) & file.exists(cand)]
  if (!length(hit)) "" else hit[1]
}

#' Load the per-exposure LDSC heritability summary (ldsc_h2_summary.tsv).
#'
#' One row per exposure with SNP h2 (+ s.e.), the LDSC intercept (model-based
#' genomic-inflation estimate), lambda_GC, mean chi^2, the inflation ratio
#' (intercept-1)/(mean_chi2-1), regression-SNP count, and exposure category/label.
#' Re-coerces numerics (collect script writes characters for NA-laden cols),
#' (re)derives the h2 z-score, and attaches the canonical category factor + short
#' label so callers don't re-implement annotation.
#'
#' @param min_h2_z optional: drop exposures whose h2 is < this many s.e. above 0
#'   (NULL keeps all). h2_z >= 4 is the usual "well-powered enough for rg" cutoff.
#' @return data.table ordered by descending h2
load_ldsc_h2 <- function(min_h2_z = NULL) {
  f <- .ldsc_summary_file("ldsc_h2_summary.tsv")
  if (!nzchar(f))
    stop("LDSC h2 summary not found (ldsc_h2_summary.tsv).\n",
         "Run the LDSC h2 array (slurm/ldsc/submit_ldsc.sh) then ",
         "scripts/ldsc/collect_ldsc_h2.R.", call. = FALSE)
  dt <- fread(f)
  num <- c("h2", "h2_se", "intercept", "intercept_se", "lambda_gc",
           "mean_chi2", "ratio", "ratio_se", "n_snp", "h2_z")
  for (cl in intersect(num, names(dt)))
    dt[[cl]] <- suppressWarnings(as.numeric(dt[[cl]]))
  if (!"failed" %in% names(dt)) dt[, failed := FALSE]
  dt[, failed := as.logical(failed)]
  dt[is.na(failed), failed := FALSE]
  dt[, h2_z := h2 / h2_se]
  if (!is.null(min_h2_z)) dt <- dt[is.finite(h2_z) & h2_z >= min_h2_z]
  if ("label" %in% names(dt)) dt[is.na(label) | label == "", label := exposure]
  else dt[, label := exposure]
  setorder(dt, -h2)
  dt[]
}

#' Load the per-exposure GWAS lead-locus summary (gwas_locus_summary.tsv).
#'
#' One row per exposure with the count of independent genome-wide-significant
#' lead loci (n_lead; 500 kb distance-clumped at p < 5e-8), all GW-significant
#' variants (n_gwsig), suggestive variants (n_suggestive), lambda_GC, mean chi^2,
#' and the top hit (rsid / locus / p). Built by the GWAS lead-variant annotation
#' step that writes output/gwas/gwas_locus_summary.tsv. Prefers the canonical
#' output tree, falling back to the staged figures/data copy.
#'
#' @return data.table ordered by descending n_lead, with category factor + label.
load_gwas_locus_summary <- function() {
  cand <- c(
    tryCatch(file.path(heap_resolve_output("gwas", must_exist = FALSE),
                       "gwas_locus_summary.tsv"), error = function(e) NA_character_),
    tryCatch(heap_figure_data_path("gwas_locus_summary.tsv"), error = function(e) NA_character_))
  f <- cand[!is.na(cand) & file.exists(cand)][1]
  if (is.na(f) || !nzchar(f))
    stop("GWAS locus summary not found (gwas_locus_summary.tsv).\n",
         "Run the GWAS lead-variant annotation step first.", call. = FALSE)
  dt <- fread(f)
  num <- c("n_variants", "lambda_gc", "mean_chi2", "n_suggestive", "n_gwsig",
           "n_lead", "min_f", "median_f", "max_log10p", "top_p")
  for (cl in intersect(num, names(dt))) dt[[cl]] <- suppressWarnings(as.numeric(dt[[cl]]))
  if (!"label" %in% names(dt)) dt[, label := exposure]
  dt[is.na(label) | label == "", label := exposure]
  setorder(dt, -n_lead)
  dt[]
}

#' Load the pairwise LDSC genetic-correlation summary (ldsc_rg_summary.tsv).
#'
#' Long format, one row per unordered exposure pair: p1, p2, rg, se, z, p, plus
#' category/label for each member (attached from the h2 summary when available).
#' Optionally returns a symmetric square rg matrix instead of the edge table.
#'
#' @param max_se optional: drop noisy pairs with rg s.e. above this (NULL keeps all)
#' @param as_matrix TRUE returns a symmetric exposure x exposure matrix (diag = 1)
#' @return data.table (edge list) or numeric matrix
load_ldsc_rg <- function(max_se = NULL, as_matrix = FALSE) {
  f <- .ldsc_summary_file("ldsc_rg_summary.tsv")
  if (!nzchar(f))
    stop("LDSC rg summary not found (ldsc_rg_summary.tsv).\n",
         "Run the LDSC h2 array, then slurm/ldsc/ldsc_rg.sh, then ",
         "scripts/ldsc/collect_ldsc_rg.R.", call. = FALSE)
  dt <- fread(f)
  for (cl in intersect(c("rg", "se", "z", "p"), names(dt)))
    dt[[cl]] <- suppressWarnings(as.numeric(dt[[cl]]))
  if (!is.null(max_se)) dt <- dt[is.finite(se) & se <= max_se]
  # attach category/label per member from the h2 summary (best effort)
  meta <- tryCatch(load_ldsc_h2(), error = function(e) NULL)
  if (!is.null(meta)) {
    lab <- setNames(meta$label, meta$exposure)
    cat <- setNames(as.character(meta$category), meta$exposure)
    dt[, `:=`(label1 = lab[p1], label2 = lab[p2],
              category1 = cat[p1], category2 = cat[p2])]
  }
  if (!as_matrix) return(dt[])
  exps <- sort(unique(c(dt$p1, dt$p2)))
  m <- matrix(NA_real_, length(exps), length(exps), dimnames = list(exps, exps))
  diag(m) <- 1
  for (i in seq_len(nrow(dt))) {
    a <- dt$p1[i]; b <- dt$p2[i]
    m[a, b] <- dt$rg[i]; m[b, a] <- dt$rg[i]
  }
  m
}

#' Per-exposure PES archetype (single source of truth; matches main Fig 6 panel d
#' and fig_pes_disease_specificity). reads = the proteome holds incremental info
#' beyond covariates (incr > 0.06); tracks = the score follows within-person change
#' (prot_lo > 0); maxgain = best held-out Cox C-index gain of the PES over the base
#' covariates. Returns data.table(exposure_id, archetype) with levels
#' modifiable / fixed / classifier / confounded / other.
heap_pes_archetype <- function(covarType = "base") {
  od  <- heap_project_output("module6_pes_longitudinal", "base")
  mpd <- heap_project_output("module6_pes_longitudinal", "multipes_disease")
  hdr <- as.data.table(load_module6_pes_longitudinal(covarType, "holdout"))
  hdr[, r2n := suppressWarnings(as.numeric(r2))]
  hdr[, aucn := suppressWarnings(as.numeric(auc))]
  inc <- hdr[, {
    typ <- if (all(is.na(r2n))) "binary" else "continuous"
    mv  <- function(mod) if (typ == "continuous") mean(r2n[model == mod], na.rm = TRUE)
                         else mean(aucn[model == mod], na.rm = TRUE)
    flo <- if (typ == "continuous") 0 else 0.5
    .(incr = mv("prot_plus_cov") - max(mv("cov_only"), flo))
  }, by = exposure_id]
  wc <- fread(file.path(od, "PESlong_base_WithinDeltaCorCI.tsv"))[, .(exposure_id, tracks = prot_lo > 0)]
  sc <- fread(file.path(mpd, "quadrant_scan.tsv"))[is.finite(dC_pes) & events >= 500]
  sc <- sc[!disease %in% c("Obesity", "Alcoholic liver disease")]
  ds <- sc[, .(maxgain = max(dC_pes)), by = exposure_id]
  ex <- merge(inc, wc, by = "exposure_id", all.x = TRUE); ex[is.na(tracks), tracks := FALSE]
  ex <- merge(ex, ds, by = "exposure_id", all.x = TRUE); ex[is.na(maxgain), maxgain := 0]
  ex[, reads := is.finite(incr) & incr > 0.06]
  ex[, archetype := fcase(reads &  tracks & maxgain > 0.03, "modifiable",
                          reads & !tracks & maxgain > 0.03, "fixed",
                          reads & maxgain <= 0.03,          "classifier",
                          !reads & maxgain > 0.03,          "confounded",
                          default = "other")]
  ex[, .(exposure_id, archetype)]
}
