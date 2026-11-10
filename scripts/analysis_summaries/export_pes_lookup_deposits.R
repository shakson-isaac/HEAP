#!/usr/bin/env Rscript
# ============================================================================
# export_pes_lookup_deposits.R
# ----------------------------------------------------------------------------
# Builds the three Supplementary Data lookup tables for Module 6, one per PES
# capability, so a reader can look up ANY exposure (or exposure-disease pair)
# rather than only the exemplars the main figure draws.
#
#   pes_read_specs.tsv     what the proteome adds BEYOND covariates for each
#                          exposure, under all five covariate specifications
#                          (panel b / the covariate-robustness supplement)
#   pes_tracking_specs.tsv within-person Delta-correlation between change in the
#                          PES and change in the reported exposure, with the
#                          covariate-only benchmark, under all five
#                          specifications (panel c)
#   pes_disease_lookup.tsv per exposure-disease pair: hazard ratio per SD of PES
#                          with CIs, and the held-out C-index ladder (panel d)
#
# NOTHING IS RECOMPUTED HERE. Each output is a stack or join of tables the
# analysis already wrote; this script only assembles and documents them.
#
# WHICH COLUMNS CARRY INTERVALS -- state this in the legends, it is not uniform:
#   read     intervals on all three nested models AND on the increment, for both
#            the primary metric (AUC/R2) and, for binary exposures, AUPR. All are
#            bootstrap over held-out PEOPLE (B=1000); the increment is a PAIRED
#            difference so its interval accounts for model correlation.
#   track    intervals on all three nested models (bootstrap over people, B=1000).
#   disease  intervals on the HAZARD RATIOS (bootstrap). The held-out C-index
#            columns are single 70/30-split point estimates with NO intervals;
#            \Tref{pes_cindex} carries bootstrapped Delta-C for a curated slate.
#
# WHICH METRIC IS WHICH: the `metric` column names what the heldout_* columns
# hold -- AUC for binary exposures, R2 for continuous. The aupr_* block is
# populated for binary exposures ONLY and is NA for continuous ones by
# construction. AUPR must be read against `prevalence` (its no-skill floor),
# not against 0.5.
#
#   Usage:  Rscript export_pes_lookup_deposits.R
#   Writes: <output>/supp_deposit/pes_{read,tracking}_specs.tsv,
#           <output>/supp_deposit/pes_disease_lookup.tsv
# ============================================================================
local({ cand <- c(Sys.getenv("HEAP_PATHS_FILE",""), "/n/groups/patel/shakson_ukb/HEAP/workflow/00_paths.R")
  h <- cand[nzchar(cand) & file.exists(cand)][1]; if (!is.na(h)) source(h) })
suppressPackageStartupMessages({ library(data.table) })
msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")

O   <- heap_project_output("module6_pes_longitudinal")
OUT <- heap_project_output("supp_deposit"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
SPECS <- c("base","base_bmi","base_clinical","base_draw","base_exclprev")

# ---- 1. READ ---------------------------------------------------------------
f <- file.path(O, "covariate_sensitivity", "covariate_sensitivity.tsv")
r <- fread(f)
# Join the paired-bootstrap intervals produced by module6_holdout_ci.R. The
# increment CI is a PAIRED difference -- covariate-only and covariate-plus-proteome
# are recomputed inside each resample -- so it accounts for the correlation between
# the two models rather than differencing two independent intervals.
ci <- rbindlist(lapply(SPECS, function(sp) {
  fp <- file.path(O, sp, sprintf("PESlong_%s_HoldoutAccuracyCI.tsv", sp))
  if (!file.exists(fp)) return(NULL)
  d <- fread(fp); d[, covariate_spec := sp][]
}), fill = TRUE)
if (nrow(ci)) {
  keep <- intersect(c("exposure_id","covariate_spec","n_person","n_obs",
                      "point","ci_lo","ci_hi","cov_point","cov_lo","cov_hi",
                      "full_point","full_lo","full_hi",
                      "increment","increment_lo","increment_hi"), names(ci))
  # NB covariate_sensitivity is computed from the SAME held-out visit-level scores
  # as module6_holdout_ci.R -- verified byte-identical on all 814 rows for the
  # covariate-only, covariate-plus-proteome and increment columns. They are one
  # estimator, not two, so the duplicates are dropped below and only the version
  # carrying bootstrap intervals is kept.
  setnames(ci,
    c("point","ci_lo","ci_hi","cov_point","cov_lo","cov_hi","full_point","full_lo","full_hi",
      "increment","increment_lo","increment_hi",
      "aupr_prot","aupr_prot_lo","aupr_prot_hi",
      "aupr_cov","aupr_cov_lo","aupr_cov_hi",
      "aupr_full","aupr_full_lo","aupr_full_hi",
      "aupr_increment","aupr_increment_lo","aupr_increment_hi",
      "prevalence_holdout","prevalence_holdout_lo","prevalence_holdout_hi"),
    c("heldout_proteome_only","heldout_proteome_only_lo","heldout_proteome_only_hi",
      "heldout_covariates_only","heldout_covariates_only_lo","heldout_covariates_only_hi",
      "heldout_covariates_plus_pes","heldout_covariates_plus_pes_lo","heldout_covariates_plus_pes_hi",
      "heldout_increment","heldout_increment_lo","heldout_increment_hi",
      # boot_ prefix keeps these clear of the point-only AUPR columns already in
      # covariate_sensitivity.tsv until the duplicate audit below resolves them
      "boot_aupr_proteome_only","boot_aupr_proteome_only_lo","boot_aupr_proteome_only_hi",
      "boot_aupr_covariates_only","boot_aupr_covariates_only_lo","boot_aupr_covariates_only_hi",
      "boot_aupr_covariates_plus_pes","boot_aupr_covariates_plus_pes_lo","boot_aupr_covariates_plus_pes_hi",
      "boot_aupr_increment","boot_aupr_increment_lo","boot_aupr_increment_hi",
      "boot_prevalence","boot_prevalence_lo","boot_prevalence_hi"), skip_absent = TRUE)
  keep <- c("exposure_id","covariate_spec","heldout_proteome_only","heldout_proteome_only_lo",
            "heldout_proteome_only_hi","heldout_covariates_only","heldout_covariates_only_lo",
            "heldout_covariates_only_hi","heldout_covariates_plus_pes",
            "heldout_covariates_plus_pes_lo","heldout_covariates_plus_pes_hi",
            "heldout_increment","heldout_increment_lo","heldout_increment_hi",
            # anchor on the AUPR/prevalence blocks specifically -- a bare "^boot_"
            # also catches boot_sd (the bootstrap SD of the proteome metric), which
            # the rename step below would strip to a meaningless column called "sd"
            # and which is redundant with the lo/hi interval already carried.
            grep("^boot_(aupr|prevalence)", names(ci), value = TRUE))
  keep <- intersect(keep, names(ci))
  r <- merge(r, ci[, keep, with = FALSE], by.x = c("exposure_id","spec"),
             by.y = c("exposure_id","covariate_spec"), all.x = TRUE)
}
setnames(r,
  c("spec","cov_only","cov_plus_pes","incremental","incremental_aupr","aupr_cov","aupr_full","aupr_lift"),
  c("covariate_spec","skill_covariates_only","skill_covariates_plus_pes","increment",
    "increment_aupr","aupr_covariates_only","aupr_covariates_plus_pes","aupr_lift"),
  skip_absent = TRUE)
# identical to their heldout_ twins (checked: max abs diff 0) -- drop, keep the
# columns that carry intervals
for (dup in c("skill_covariates_only","skill_covariates_plus_pes","increment"))
  if (dup %in% names(r) && paste0("heldout_", sub("^skill_","",dup)) %in% names(r)) r[[dup]] <- NULL
if ("increment" %in% names(r) && "heldout_increment" %in% names(r)) r[, increment := NULL]

# --- duplicate audit for the AUPR block -------------------------------------
# covariate_sensitivity.tsv already carries point-estimate AUPR for two of the
# three models plus the prevalence. The bootstrap run recomputes them from the
# SAME held-out scores, so they are expected to coincide -- verify rather than
# assume, then keep only the copy carrying intervals. Shipping both without
# checking is exactly how the AUC block ended up silently duplicated.
for (p in list(c("aupr_covariates_only",     "boot_aupr_covariates_only"),
               c("aupr_covariates_plus_pes", "boot_aupr_covariates_plus_pes"),
               c("increment_aupr",           "boot_aupr_increment"),
               c("prevalence",               "boot_prevalence"))) {
  if (!all(p %in% names(r))) next
  d    <- suppressWarnings(max(abs(r[[p[1]]] - r[[p[2]]]), na.rm = TRUE))
  same <- is.finite(d) && d < 1e-9
  msg(sprintf("  dup-check %-26s max|diff| = %.3e  %s", p[1], d,
              if (same) "IDENTICAL -> dropping point-only copy" else "DIFFERS -> keeping both"))
  if (same) r[[p[1]]] <- NULL
}
for (b in grep("^boot_", names(r), value = TRUE)) {
  tgt <- sub("^boot_", "", b)
  if (tgt %in% names(r)) msg("  ! kept BOTH ", tgt, " and ", b) else setnames(r, b, tgt)
}

# Name the metric the heldout_* columns actually hold, and rename aupr_lift so it
# cannot be misread as the PES increment -- upstream defines it as
# (covariates+PES AUPR) minus prevalence, i.e. lift over a RANDOM classifier.
r[, metric := fifelse(exposure_type == "binary", "AUC", "R2")]
setnames(r, "aupr_lift", "aupr_covariates_plus_pes_over_prevalence", skip_absent = TRUE)

setcolorder(r, intersect(c("exposure_id","category","exposure_type","metric","covariate_spec",
                           "heldout_covariates_only","heldout_covariates_only_lo","heldout_covariates_only_hi",
                           "heldout_proteome_only","heldout_proteome_only_lo","heldout_proteome_only_hi",
                           "heldout_covariates_plus_pes","heldout_covariates_plus_pes_lo","heldout_covariates_plus_pes_hi",
                           "heldout_increment","heldout_increment_lo","heldout_increment_hi",
                           "prevalence","prevalence_lo","prevalence_hi",
                           "aupr_covariates_only","aupr_covariates_only_lo","aupr_covariates_only_hi",
                           "aupr_proteome_only","aupr_proteome_only_lo","aupr_proteome_only_hi",
                           "aupr_covariates_plus_pes","aupr_covariates_plus_pes_lo","aupr_covariates_plus_pes_hi",
                           "aupr_increment","aupr_increment_lo","aupr_increment_hi",
                           "aupr_covariates_plus_pes_over_prevalence"), names(r)))
setorder(r, exposure_id, covariate_spec)
fwrite(r, file.path(OUT, "pes_read_specs.tsv"), sep = "\t")
msg("pes_read_specs.tsv        ", nrow(r), " rows | ", uniqueN(r$covariate_spec), " specs | ",
    uniqueN(r$exposure_id), " exposures")

# ---- 2. TRACK --------------------------------------------------------------
tl <- rbindlist(lapply(SPECS, function(s) {
  f <- file.path(O, s, sprintf("PESlong_%s_WithinDeltaCorCI.tsv", s))
  if (!file.exists(f)) { msg("  ! missing tracking file for ", s); return(NULL) }
  d <- fread(f); d[, covariate_spec := s][]
}), fill = TRUE)
setnames(tl,
  c("dcor_prot","prot_lo","prot_hi","dcor_cov","cov_lo","cov_hi","dcor_full","full_lo","full_hi"),
  c("dcor_pes","dcor_pes_lo","dcor_pes_hi","dcor_covariates","dcor_covariates_lo","dcor_covariates_hi",
    "dcor_covariates_plus_pes","dcor_covariates_plus_pes_lo","dcor_covariates_plus_pes_hi"),
  skip_absent = TRUE)
setcolorder(tl, intersect(c("exposure_id","category","exposure_type","covariate_spec",
                            "n_person","n_pairs","n_change",
                            "dcor_covariates","dcor_covariates_lo","dcor_covariates_hi",
                            "dcor_pes","dcor_pes_lo","dcor_pes_hi",
                            "dcor_covariates_plus_pes","dcor_covariates_plus_pes_lo",
                            "dcor_covariates_plus_pes_hi"), names(tl)))
setorder(tl, exposure_id, covariate_spec)
fwrite(tl, file.path(OUT, "pes_tracking_specs.tsv"), sep = "\t")
msg("pes_tracking_specs.tsv    ", nrow(tl), " rows | ", uniqueN(tl$covariate_spec), " specs | ",
    uniqueN(tl$exposure_id), " exposures")

# ---- 3. PREDICT ------------------------------------------------------------
M  <- file.path(O, "multipes_disease")
sc <- fread(file.path(M, "pes_disease_scale.tsv"))
# Prefer the interval-carrying version. It is the SAME 70/30 split re-scored with
# a test-set bootstrap, and collect_dc_holdout_ci.R refuses to write unless the
# point estimates reproduce the screen exactly (verified max|diff| = 0), so this
# adds columns without moving a single published number.
HO_CI <- file.path(M, "pes_disease_dc_holdout_ci.tsv")
ho <- if (file.exists(HO_CI)) {
  msg("  using pes_disease_dc_holdout_ci.tsv (point estimates + bootstrap intervals)")
  fread(HO_CI)
} else {
  msg("  ! no _ci file; falling back to point estimates only")
  fread(file.path(M, "pes_disease_dc_holdout.tsv"))
}
# `n`/`events` are identical in both; keep one copy
d  <- merge(sc, ho[, setdiff(names(ho), c("n","events")), with = FALSE],
            by = c("exposure_id","disease"), all = TRUE)
# The in-sample ("apparent") Delta-C columns are DROPPED. They are biased toward
# zero because the ~46-covariate baseline overfits in sample, and they disagree in
# SIGN with the held-out columns (e.g. -0.0044 apparent vs +0.0079 held-out on the
# first row). Shipping both invites a reader to quote the wrong one.
d[, c("dC_overcov_app","dC_beyondE_app") := NULL]
# readable disease label alongside the raw UK Biobank field name
d[, disease_label := gsub("_", " ", sub("_f[0-9]+_0_0$", "", sub("^age_[a-z][0-9]+_first_reported_", "", disease)))]
d[, disease_label := paste0(toupper(substr(disease_label,1,1)), substr(disease_label,2,nchar(disease_label)))]
setnames(d, "disease", "disease_field_id")
setnames(d,
  c("C0_ho","CP_ho","CE_ho","CEP_ho",
    "dC_overcov_ho","dC_beyondE_ho","dC_e_ho"),
  c("C_covariates_heldout","C_covariates_pes_heldout","C_covariates_selfreport_heldout",
    "C_covariates_selfreport_pes_heldout",
    "dC_over_covariates_heldout","dC_beyond_selfreport_heldout","dC_selfreport_heldout"),
  skip_absent = TRUE)
# intervals on the same split (test-set bootstrap, B=400)
setnames(d,
  c("nboot_ho","C0_ho_lo","C0_ho_hi","CP_ho_lo","CP_ho_hi","CE_ho_lo","CE_ho_hi",
    "CEP_ho_lo","CEP_ho_hi","dC_overcov_ho_lo","dC_overcov_ho_hi","dC_overcov_ho_pgt0",
    "dC_beyondE_ho_lo","dC_beyondE_ho_hi","dC_beyondE_ho_pgt0"),
  c("nboot_heldout","C_covariates_heldout_lo","C_covariates_heldout_hi",
    "C_covariates_pes_heldout_lo","C_covariates_pes_heldout_hi",
    "C_covariates_selfreport_heldout_lo","C_covariates_selfreport_heldout_hi",
    "C_covariates_selfreport_pes_heldout_lo","C_covariates_selfreport_pes_heldout_hi",
    "dC_over_covariates_heldout_lo","dC_over_covariates_heldout_hi",
    "dC_over_covariates_heldout_pgt0",
    "dC_beyond_selfreport_heldout_lo","dC_beyond_selfreport_heldout_hi",
    "dC_beyond_selfreport_heldout_pgt0"),
  skip_absent = TRUE)
setnames(d, c("hr_pesE","hr_pesE_lo","hr_pesE_hi","wald_pesE","lrt_pesE","q_pesE"),
         c("hr_pes_beyond_selfreport","hr_pes_beyond_selfreport_lo","hr_pes_beyond_selfreport_hi",
           "wald_pes_beyond_selfreport","lrt_pes_beyond_selfreport","q_pes_beyond_selfreport"),
         skip_absent = TRUE)
# --- bootstrap C-index ladder for the pairs that carry one ------------------
# Two estimators coexist here and are NOT interchangeable:
#   * the C_*_heldout columns above are ONE 70/30 split (seed 42), the only
#     estimator that covers all 29,865 pairs;
#   * the boot_* columns below refit all four Cox models inside each of B=400
#     resamples and evaluate out-of-bag, which is what Fig 6d plots and what
#     \Tref{pes_cindex} ships -- but only for the ladder slate.
# The bootstrap runs ~1,600 Cox fits per pair against the split's 4, so it is
# not feasible grid-wide. Because training uses only ~63% unique rows, the
# bootstrap is mildly pessimistic on the ABSOLUTE C-index (mean -0.003 across
# the shared pairs, and larger for the richer models) while the INCREMENTS
# agree to ~0.0005 -- both models pay the same penalty, so it cancels.
LADF <- file.path(M, "pes_disease_ladder.tsv")
if (file.exists(LADF)) {
  L <- fread(LADF)
  bcols <- c(boot_C_covariates="C_cov", boot_C_covariates_pes="C_covP",
             boot_C_covariates_selfreport="C_covE", boot_C_covariates_selfreport_pes="C_covEP",
             boot_dC_over_covariates="dC_overCov", boot_dC_over_covariates_lo="overCov_lo",
             boot_dC_over_covariates_hi="overCov_hi",
             boot_dC_beyond_selfreport="dC_beyondE", boot_dC_beyond_selfreport_lo="beyondE_lo",
             boot_dC_beyond_selfreport_hi="beyondE_hi", boot_n="nboot")
  bcols <- bcols[bcols %in% names(L)]
  Lb <- L[, c("exposure_id","dz_id", unname(bcols)), with = FALSE]
  setnames(Lb, unname(bcols), names(bcols))
  d <- merge(d, Lb, by.x = c("exposure_id","disease_field_id"), by.y = c("exposure_id","dz_id"),
             all.x = TRUE)
  msg("  bootstrap ladder: ", nrow(L), " rows; ", sum(!is.na(d$boot_dC_over_covariates)),
      " joined into the grid (", nrow(L) - sum(!is.na(d$boot_dC_over_covariates)),
      " are for diseases outside the 181-disease set)")
}
# Flag the exemplars actually DRAWN in Fig 6d, so a reader can find the figure's
# points without guessing. Keep in sync with EX in fig_m6_panel_d.R.
FIGP <- data.table(
  exposure_id = c("current_tobacco_smoking_f1239_0_0_Yes._on_most_or_all_days",
                  "pack_years_of_smoking_f20161_0_0",
                  "alcohol_intake_frequency_f1558_0_0",
                  "number_of_days_week_of_vigorous_physical_activity_10_plus_minutes_f904_0_0",
                  "coffee_intake_f1498_0_0", "usual_walking_pace_f924_0_0"),
  dzpat = c("j44_first_reported_other_chronic","j44_first_reported_other_chronic",
            "f10_first_reported_mental","e11_first_reported_non_insulin",
            "e11_first_reported_non_insulin","e11_first_reported_non_insulin"))
d[, in_main_figure := FALSE]
for (i in seq_len(nrow(FIGP)))
  d[exposure_id == FIGP$exposure_id[i] & grepl(FIGP$dzpat[i], disease_field_id), in_main_figure := TRUE]
msg("  in_main_figure flagged on ", sum(d$in_main_figure), " rows (expect ", nrow(FIGP), ")")
if (sum(d$in_main_figure) != nrow(FIGP))
  msg("  ! WARNING: Fig 6d exemplar list is out of sync with fig_m6_panel_d.R")

setcolorder(d, intersect(c("exposure_id","disease_label","disease_field_id","in_main_figure","n","events",
                           "hr_pes","hr_pes_lo","hr_pes_hi","wald_pes","lrt_pes","q_pes",
                           "hr_pes_beyond_selfreport","hr_pes_beyond_selfreport_lo",
                           "hr_pes_beyond_selfreport_hi","wald_pes_beyond_selfreport",
                           "lrt_pes_beyond_selfreport","q_pes_beyond_selfreport",
                           "C_covariates_heldout","C_covariates_pes_heldout",
                           "C_covariates_selfreport_heldout","C_covariates_selfreport_pes_heldout",
                           "dC_over_covariates_heldout","dC_beyond_selfreport_heldout",
                           "dC_selfreport_heldout",
                           "boot_n","boot_C_covariates","boot_C_covariates_pes",
                           "boot_C_covariates_selfreport","boot_C_covariates_selfreport_pes",
                           "boot_dC_over_covariates","boot_dC_over_covariates_lo",
                           "boot_dC_over_covariates_hi","boot_dC_beyond_selfreport",
                           "boot_dC_beyond_selfreport_lo","boot_dC_beyond_selfreport_hi"), names(d)))
# --- drop sex-limited outcomes ----------------------------------------------
# The Cox model adjusts for sex, age*sex and age^2*sex. On an outcome only one sex
# can have, those terms perfectly predict who is at risk: the fit is separated, the
# sex coefficient diverges, and coxph stops wherever the convergence tolerance
# lands. The estimates are then numerically unstable rather than merely imprecise --
# chronotype -> ovarian disorders moved HR 0.860 to 1.103 between two runs of
# identical code on identical data. 11 such outcomes x 165 exposures = 1,815 rows,
# 435 of which carried q < 0.05 and would read as findings. Excluded until they are
# refit sex-stratified without the sex terms.
SEXLIMITED <- paste("penis","prostat","testis","testic","scrotum","male_genital",
                    "ovar","uter","cervi","vagin","vulv","breast","menstr",
                    "menopaus","pregnan","female_genital","fallopian", sep = "|")
n_before <- nrow(d)
sexsp <- grepl(SEXLIMITED, d$disease_field_id, ignore.case = TRUE)
msg("  dropping ", sum(sexsp), " sex-limited rows across ",
    uniqueN(d$disease_field_id[sexsp]), " outcomes (separated Cox fits)")
d <- d[!sexsp]

setorder(d, exposure_id, disease_field_id)
fwrite(d, file.path(OUT, "pes_disease_lookup.tsv"), sep = "\t")
msg("pes_disease_lookup.tsv    ", nrow(d), " rows | ", uniqueN(d$exposure_id), " exposures x ",
    uniqueN(d$disease_field_id), " diseases")
# --- alternative-specification HR scans -------------------------------------
# base_prevalent is deliberately NOT released: adjusting for prevalent disease
# conditions on a descendant of exposure, the collider path the reviews raised,
# and it was withdrawn from every other deposit for that reason.
# The first three vary ONLY the Cox adjustment: same PES, same participants.
# base_exclprev is different in kind -- it restricts to participants free of major
# disease at baseline AND uses the PES retrained on that subset, so sample and
# score both move. It is the reverse-causation check, not a covariate spec, and
# its exposure list differs from base in both directions (162 shared of 165/164),
# so cross-spec comparison is valid only on the intersection.
# NAMED BY SPECIFICATION, not by an expanded phrase. pes_scores/reading/ and
# pes_scores/tracking/ are split on the covariate_spec column, so their files are
# base_bmi.tsv, base_draw.tsv, base_exclprev.tsv ... A reader opening pes_scores/
# should not find base_bmi.tsv beside base_plus_bmi.tsv and have to work out that
# they are the same specification. The expansion belongs in the manifest, which
# has room for a sentence.
SPEC_SCANS <- c(base_bmi      = "pes_disease_scale_base_bmi.tsv",
                base_draw     = "pes_disease_scale_base_draw.tsv",
                base_clinical = "pes_disease_scale_base_clinical.tsv",
                base_exclprev = "pes_disease_scale_base_exclprev.tsv")
sdir <- file.path(OUT, "pes_disease_specs"); dir.create(sdir, recursive = TRUE, showWarnings = FALSE)
# The zip ships whatever is in this folder. A file left behind by an earlier naming
# scheme would be shipped alongside its own replacement, so clear before writing.
stale <- setdiff(basename(Sys.glob(file.path(sdir, "*.tsv"))),
                 c(paste0(names(SPEC_SCANS), ".tsv"), "manifest.tsv"))
if (length(stale)) {
  msg("  removing ", length(stale), " stale spec file(s): ", paste(stale, collapse = ", "))
  invisible(file.remove(file.path(sdir, stale)))
}
# Where a held-out C-index screen was run for the same specification, carry its
# columns alongside the hazard ratios so one file answers both questions for a
# pair -- as pes_disease_lookup.tsv already does for base. The apparent dC the
# scan writes is NOT a substitute: across the base grid it correlates with the
# held-out value at Pearson -0.02, so it is dropped wherever held-out exists.
DC_SCANS <- c(base_exclprev = "pes_disease_dc_holdout_base_exclprev.tsv")
for (nm in names(SPEC_SCANS)) {
  f <- file.path(M, SPEC_SCANS[[nm]])
  if (!file.exists(f)) { msg("  ! missing spec scan: ", nm); next }
  x <- fread(f)
  x <- x[!grepl(SEXLIMITED, disease, ignore.case = TRUE)]
  dcf <- if (nm %in% names(DC_SCANS)) file.path(M, DC_SCANS[[nm]]) else NA_character_
  if (!is.na(dcf) && file.exists(dcf)) {
    dc <- fread(dcf, select = c("exposure_id","disease","C0_ho","CP_ho","CE_ho","CEP_ho",
                                "dC_overcov_ho","dC_beyondE_ho","dC_e_ho"))
    n0 <- nrow(x)
    x <- merge(x, dc, by = c("exposure_id","disease"), all.x = TRUE)
    if (nrow(x) != n0) stop("dC merge changed row count for ", nm, ": ", n0, " -> ", nrow(x))
    x[, c("dC_overcov_app","dC_beyondE_app") := NULL]
    msg("    + held-out C-index for ", nm, ": ", sum(is.finite(x$dC_overcov_ho)), " pairs")
  }
  fwrite(x, file.path(sdir, paste0(nm, ".tsv")), sep = "\t")
  msg("  spec scan ", nm, ": ", nrow(x), " pairs")
}
# every other split folder in the archive ships a manifest; this one should too
sm <- rbindlist(lapply(names(SPEC_SCANS), function(nm) {
  fp <- file.path(sdir, paste0(nm, ".tsv")); if (!file.exists(fp)) return(NULL)
  x <- fread(fp, select = c("exposure_id","disease","q_pes"))
  data.table(file = paste0(nm, ".tsv"), specification = nm,
             covariate_set = c(base_bmi = "base_bmi", base_draw = "base_draw",
                               base_clinical = "base_clinical",
                               base_exclprev = "base")[[nm]],
             description = c(
               base_bmi      = "base adjustment plus body mass index",
               base_draw     = "base adjustment plus blood-draw conditions (fasting time, season)",
               base_clinical = "base adjustment plus body mass index, draw conditions and five medication classes",
               base_exclprev = "base adjustment, restricted to participants free of major disease at baseline")[[nm]],
             what_varies = c(base_bmi = "Cox adjustment only",
                             base_draw = "Cox adjustment only",
                             base_clinical = "Cox adjustment only",
                             base_exclprev =
                               "analysis sample and score: restricted to participants free of major disease at baseline, PES retrained on that subset")[[nm]],
             pairs = nrow(x), exposures = uniqueN(x$exposure_id),
             diseases = uniqueN(x$disease), significant_q05 = sum(x$q_pes < 0.05, na.rm = TRUE))
}), fill = TRUE)
if (nrow(sm)) { fwrite(sm, file.path(sdir, "manifest.tsv"), sep = "\t")
                msg("  spec scan manifest: ", nrow(sm), " rows") }

msg("  HR intervals present: ", sum(is.finite(d$hr_pes_lo)), " rows | held-out dC present: ",
    sum(is.finite(d$dC_over_covariates_heldout)), " rows")
