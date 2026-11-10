#!/usr/bin/env Rscript
# ============================================================================
# export_med_spec_deposit.R   (DEPOSIT -- Supplementary Data, mediation)
# ----------------------------------------------------------------------------
# Per-covariate-specification mediation tables for HEAP_Supplementary_Data.zip,
# under mediation/exposome/<spec>.tsv and mediation/genetic/<spec>.tsv.
#
# WHY. summarize_module3_mediation.R deposits the BASE specification only. The
# other seven specs ran -- output/module3/M3_*_primary all exist -- but were only
# ever aggregated into the paired attenuation.tsv that the robustness figure
# needs, and that file covers just +BMI and +clinical. Mediation is where the
# headline caveat lives (NIE roughly halves under BMI adjustment, which this
# design cannot separate from confounding), so a reader should be able to pull
# any specification directly.
#
# CANONICAL SOURCE RUNS -- this header is the contract:
#   output/module3/<experiment>/<covariate_set>/<family>/primary_total/MDres_*.txt
# The spec grid is the one in module3_spec_sensitivity.R:28-35, unchanged.
# Only primary_total exists for the alternative specs, so these files carry the
# aggregate exposomic driver (PXS_total) and the combined genetic driver
# (PGS_total). The 13 per-category drivers are base-only and stay in
# med_exposure_categories.tsv.
#
# CONSTRUCTION is copied line-for-line from summarize_module3_mediation.R --
# same KEEP columns, same NIE/NDE pivot, same total = NIE + NDE, the same
# proportion-mediated rule INCLUDING its clamp of PM>1 to 1 (the retired script's
# header said it excluded them; its code clamped, and the code is what produced
# the shipped numbers), the same BH-q within predictor_class, and the same
# estimable rule that separates a structural zero from an unmeasured cell.
# Deviating anywhere here would make the base file disagree with the shipped
# med_exposure_total.tsv, which --finalize checks and refuses to accept.
#
# TWO PASSES, because reading one spec is ~800 files and several GB:
#   --spec <experiment>   read one run, cache links to intermediate/
#   --finalize            write every cached spec out, verify base, manifest
#
# Run (normally via slurm/module3/HEAPmed_spec_deposit.sh):
#   Rscript export_med_spec_deposit.R --spec M3_base_bmi_lasso_primary
#   Rscript export_med_spec_deposit.R --finalize
# ============================================================================
suppressPackageStartupMessages(library(data.table))
setDTthreads(0)
local({
  cand <- c(Sys.getenv("HEAP_PATHS_FILE", unset = ""),
            "/n/groups/patel/shakson_ukb/HEAP/workflow/00_paths.R")
  hit <- cand[nzchar(cand) & file.exists(cand)][1]; if (!is.na(hit)) source(hit)
})

M3    <- heap_project_output("module3")
DEP   <- heap_project_output("supp_deposit")
CACHE <- heap_project_intermediate("med_spec_cache")
dir.create(CACHE, recursive = TRUE, showWarnings = FALSE)
msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")

# spec grid -- identical to module3_spec_sensitivity.R:28-35
# NOTE: prevalent disease is released as an EXCLUSION only. Conditioning on it
# would adjust for a variable downstream of exposure, opening the collider path the
# reviews raised; the exclusion answers the same question without that risk. The
# SES/deprivation estimand is likewise withheld: moving deprivation out of the
# exposome into the covariate set deletes a whole exposure domain, so it is a
# re-specification of the estimand rather than a robustness check.
SPECS <- data.table(
  exp  = c("M3_base_lasso_primary", "M3_base_bmi_lasso_primary", "M3_base_clinical_lasso_primary",
           "M3_base_exclprev_lasso_primary",
           "M3_base_draw_lasso_primary", "M3_base_ridge_primary", "M3_base_enet_primary"),
  cov  = c("base", "base_bmi", "base_clinical", "base", "base_draw", "base", "base"),
  fam  = c("lasso", "lasso", "lasso", "lasso", "lasso", "ridge", "enet"),
  file = c("base", "base_plus_bmi", "base_plus_clinical",
           "exclude_prevalent_disease", "base_plus_blood_draw", "estimator_ridge",
           "estimator_elastic_net"),
  kind = c("primary", "covariate", "covariate", "sample", "covariate", "model", "model"))

KEEP <- c("protID","DZ_ID","predictor","predictor_class","effect_type",
          "effect_logHR","effect_HR","delta_se","delta_l95_HR","delta_u95_HR",
          "delta_p","protein_HR","protein_p","n","n_cases","instrument_present")

parse_dz <- function(x) {
  s <- sub("^age_", "", x); s <- sub("_f[0-9]+_[0-9]+_[0-9]+$", "", s)
  code <- toupper(sub("^([a-z][0-9]+)_.*$", "\\1", s))
  nm   <- sub("^[a-z][0-9]+_", "", s); nm <- sub("^first_reported_", "", nm)
  list(icd = code, disease = gsub("_", " ", nm))
}

read_run <- function(dir) {
  f <- list.files(dir, "^MDres_.*\\.txt$", full.names = TRUE)
  if (!length(f)) stop("no MDres files under ", dir)
  msg("  reading ", length(f), " files from ", dir)
  acc <- vector("list", length(f))
  for (i in seq_along(f)) {
    d <- fread(f[i], select = KEEP, showProgress = FALSE)
    nie <- d[effect_type == "NIE",
             .(protID, DZ_ID, predictor, predictor_class, n, n_cases,
               protein_HR, protein_p, driver_present = instrument_present,
               nie_logHR = effect_logHR, nie_HR = effect_HR, nie_se = delta_se,
               nie_l95 = delta_l95_HR, nie_u95 = delta_u95_HR, nie_p = delta_p)]
    nde <- d[effect_type == "NDE",
             .(protID, DZ_ID, predictor, nde_logHR = effect_logHR, nde_HR = effect_HR,
               nde_l95 = delta_l95_HR, nde_u95 = delta_u95_HR, nde_p = delta_p)]
    acc[[i]] <- merge(nie, nde, by = c("protID", "DZ_ID", "predictor"))
    if (i %% 200 == 0) msg("    ", i, "/", length(f))
  }
  rbindlist(acc, fill = TRUE)
}

# the derived quantities. The rule itself lives in med_link_rules.R (NUM-6) so
# this script and module3_spec_sensitivity.R cannot drift apart again.
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])),
                 "med_link_rules.R"))
finish_links <- function(L) {
  numcols <- c("nie_logHR","nie_HR","nie_l95","nie_u95","nie_p",
               "nde_logHR","nde_HR","protein_HR","protein_p","n","n_cases")
  for (cc in intersect(numcols, names(L))) if (!is.numeric(L[[cc]]))
    set(L, j = cc, value = suppressWarnings(as.numeric(L[[cc]])))
  L <- med_derive_links(L)
  dz <- parse_dz(L$DZ_ID); L[, `:=`(icd = dz$icd, disease = dz$disease)]
  L[]
}

COLS <- function(lab) c(lab, "protID", "icd", "disease", "DZ_ID", "n_cases",
  "nie_HR", "nie_l95", "nie_u95", "nde_HR", "nde_l95", "nde_u95",
  "total_HR", "prop_mediated", "protein_HR", "protein_p",
  "nie_p", "nde_p", "q", "sig", "estimable", "driver_present", "pm_clamped")
rnd <- function(d) { for (c in c("nie_HR","nie_l95","nie_u95","nde_HR","total_HR",
                                 "prop_mediated","protein_HR"))
                       if (c %in% names(d)) set(d, j = c, value = signif(d[[c]], 4)); d }

# ---------------------------------------------------------------------------
args <- commandArgs(TRUE)
if ("--spec" %in% args) {
  exp <- args[which(args == "--spec") + 1L]
  i <- which(SPECS$exp == exp); if (!length(i)) stop("unknown spec: ", exp)
  dir <- file.path(M3, SPECS$exp[i], SPECS$cov[i], SPECS$fam[i], "primary_total")
  msg("PASS 1  ", exp, "  (", SPECS$file[i], ")")
  L <- finish_links(read_run(dir))
  out <- file.path(CACHE, paste0(SPECS$exp[i], ".rds"))
  saveRDS(L, out)
  msg("cached ", nrow(L), " links -> ", out, "  | significant: ", sum(L$sig))
  quit(status = 0)
}

if (!"--finalize" %in% args) stop("give --spec <experiment> or --finalize")

# ---- PASS 2 ---------------------------------------------------------------
OUT_E <- file.path(DEP, "med_specs_exposome"); OUT_G <- file.path(DEP, "med_specs_genetic")
for (d in c(OUT_E, OUT_G)) { unlink(d, recursive = TRUE); dir.create(d, recursive = TRUE, showWarnings = FALSE) }
have <- SPECS[file.exists(file.path(CACHE, paste0(exp, ".rds")))]
if (!nrow(have)) stop("no cached specs under ", CACHE, " -- run pass 1 first")
if (!"M3_base_lasso_primary" %in% have$exp) stop("base spec missing; cannot verify")
msg("PASS 2  finalizing ", nrow(have), " of ", nrow(SPECS), " specs")

man_e <- man_g <- list()
for (i in seq_len(nrow(have))) {
  L <- readRDS(file.path(CACHE, paste0(have$exp[i], ".rds")))
  ex <- L[predictor_class == "exposure_total"]; ex[, exposure_category := "PXS_total"]
  ex[, `:=`(spec = have$file[i], spec_kind = have$kind[i], run_id = have$exp[i])]
  fwrite(rnd(ex)[, c(COLS("exposure_category"), "spec", "spec_kind", "run_id"), with = FALSE],
         file.path(OUT_E, paste0(have$file[i], ".tsv")), sep = "\t")
  gn <- L[predictor_class == "genetic_total"]; gn[, component := "PGS_total"]
  gn[, `:=`(spec = have$file[i], spec_kind = have$kind[i], run_id = have$exp[i])]
  fwrite(rnd(gn)[, c(COLS("component"), "spec", "spec_kind", "run_id"), with = FALSE],
         file.path(OUT_G, paste0(have$file[i], ".tsv")), sep = "\t")
  man_e[[i]] <- data.table(file = paste0(have$file[i], ".tsv"), specification = have$file[i],
                           spec_kind = have$kind[i], run_id = have$exp[i],
                           links = nrow(ex), significant = sum(ex$sig))
  man_g[[i]] <- data.table(file = paste0(have$file[i], ".tsv"), specification = have$file[i],
                           spec_kind = have$kind[i], run_id = have$exp[i],
                           links = nrow(gn), significant = sum(gn$sig))
  msg("  ", have$file[i], ": exposome ", nrow(ex), " (", sum(ex$sig), " sig) | genetic ",
      nrow(gn), " (", sum(gn$sig), " sig)")
  rm(L, ex, gn); gc(FALSE)
}
fwrite(rbindlist(man_e), file.path(OUT_E, "manifest.tsv"), sep = "\t")
fwrite(rbindlist(man_g), file.path(OUT_G, "manifest.tsv"), sep = "\t")

# ---- rule 5: the base file must reproduce what already ships ---------------
# med_exposure_total.tsv is CITED. If this pipeline's base disagrees with it, the
# per-spec files are not comparable to the published numbers and the difference
# is an author decision, not something to overwrite silently.
ship <- file.path(DEP, "med_exposure_total.tsv")
if (file.exists(ship)) {
  a <- fread(ship, showProgress = FALSE)
  b <- fread(file.path(OUT_E, "base.tsv"), showProgress = FALSE)
  b <- b[, intersect(names(a), names(b)), with = FALSE]
  setcolorder(b, names(a)); setkeyv(a, c("protID","DZ_ID")); setkeyv(b, c("protID","DZ_ID"))
  same_rows <- nrow(a) == nrow(b)
  diffs <- if (!same_rows) NA_integer_ else sum(!mapply(function(x, y)
    isTRUE(all.equal(x, y, tolerance = 1e-8)), a, b))
  if (!same_rows || (!is.na(diffs) && diffs > 0)) {
    msg("!! BASE DOES NOT REPRODUCE med_exposure_total.tsv")
    msg("   shipped rows ", nrow(a), " vs regenerated ", nrow(b),
        if (!is.na(diffs)) paste0(" | differing columns: ", diffs) else "")
    stop("base mismatch -- do NOT ship these files; this is a rule-5 decision")
  }
  msg("base reproduces the shipped med_exposure_total.tsv exactly (", nrow(a), " links)")
} else msg("NOTE: ", ship, " absent -- base not verified")
msg("done. ", OUT_E, " and ", OUT_G)
