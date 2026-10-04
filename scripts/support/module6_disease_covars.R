#!/usr/bin/env Rscript
# ============================================================================
# module6_disease_covars.R  (shared: Cox covariate frame for the disease models)
# ----------------------------------------------------------------------------
# The Module 6 disease scripts hardcoded their Cox adjustment as
#   age0 + sex + age2 + age_sex + age2_sex + centre + 20 PCs
# which is exactly the `base` set in config/covariates/covariate_sets.yml. That
# hardcoding is why passing an alternative covarType to the quadrant scripts
# silently changed nothing: the score directory moved, the ADJUSTMENT did not.
#
# This helper drives the frame and the formula off the YAML instead, so
# "does the PES still predict disease after adjusting for BMI / the clinical
# set?" becomes answerable rather than assumed.
#
# CONTRACT: with covar_set = "base" the returned formula is term-for-term the
# old hardcoded string, so existing outputs must reproduce exactly. Callers
# should assert that before trusting any other set.
#
# CAUTION -- the added covariates are not missing-at-random. BMI and fasting
# time carry real missingness, so a richer set SHRINKS the analysis sample. A
# drop in C-index between sets therefore confounds "adjustment absorbed the
# signal" with "different people". build_disease_covars() reports the complete
# case count for exactly this reason; compare sets on the intersection, or
# report n alongside every estimate.
#
#   source(".../module6_disease_covars.R")
#   z <- build_disease_covars(cvl, "base_bmi")   # -> list(dt, formula, covars, n)
# ============================================================================
suppressPackageStartupMessages({ library(data.table) })

# Columns the disease scripts synthesise/rename rather than read verbatim.
# Everything else is taken from covars_long under its YAML name.
.M6_RENAME <- c(age_when_attended_assessment_centre_f21003_0_0 = "age0",
                sex_f31_0_0                                    = "sex",
                uk_biobank_assessment_centre_f54_0_0           = "centre")

build_disease_covars <- function(cvl, covar_set = "base", instance = 0L, verbose = TRUE) {
  if (!exists("load_covariate_set")) stop("source workflow/config_helpers.R first")
  want <- load_covariate_set(covar_set)

  have <- names(cvl)
  miss <- setdiff(want, have)
  if (length(miss)) stop("covariate_set '", covar_set, "' names ", length(miss),
                         " column(s) absent from covars_long: ", paste(miss, collapse=", "))

  # Evaluate the instance filter OUTSIDE the data.table frame. Writing
  #   cvl[get("instance") == instance, ...]
  # silently resolves BOTH sides to the column, so the test is always TRUE and
  # every instance is kept -- which duplicated the 61,321 participants who have
  # imaging visits and inflated every downstream Cox fit. Caught only by the
  # base-reproduction gate.
  if (!"instance" %in% names(cvl)) stop("covars_long has no `instance` column")
  keep <- cvl[["instance"]] == instance
  d <- cvl[keep, c("eid", want), with = FALSE]
  if (anyDuplicated(d[["eid"]]))
    stop("build_disease_covars: ", sum(duplicated(d[["eid"]])),
         " duplicated eid(s) after filtering to instance ", instance,
         " -- one row per participant is required for the Cox models")
  setnames(d, names(.M6_RENAME), unname(.M6_RENAME), skip_absent = TRUE)
  cols <- setdiff(names(d), "eid")

  # Character/low-cardinality columns become factors; the rest numeric. centre
  # and assessment_season are categorical; medication flags and BMI are numeric.
  for (cc in cols) {
    v <- d[[cc]]
    if (is.character(v) || is.factor(v) || (is.logical(v))) {
      d[, (cc) := factor(get(cc))]
    } else {
      d[, (cc) := as.numeric(get(cc))]
      # a numeric column with <=2 distinct non-NA values is an indicator; leaving
      # it numeric is fine for coxph and keeps the term count predictable
    }
  }
  # drop factor levels with no data, which otherwise make coxph rank-deficient
  for (cc in cols) if (is.factor(d[[cc]])) d[, (cc) := droplevels(get(cc))]

  n_all <- nrow(d)
  d_cc  <- d[stats::complete.cases(d)]
  fml   <- paste(cols, collapse = "+")

  if (verbose)
    message(sprintf("[covars] %-14s %2d terms | %d rows, %d complete (%.1f%% lost to missingness)",
                    covar_set, length(cols), n_all, nrow(d_cc), 100*(1 - nrow(d_cc)/max(n_all,1))))

  list(dt = d, dt_complete = d_cc, formula = fml, covars = cols,
       n = n_all, n_complete = nrow(d_cc), covar_set = covar_set)
}

# Assert that `base` still reproduces the string the scripts used to hardcode.
assert_base_formula_unchanged <- function(fml) {
  PCS <- paste0("genetic_principal_components_f22009_0_", 1:20)
  old <- paste("age0+sex+age2+age_sex+age2_sex+centre", paste(PCS, collapse="+"), sep="+")
  new_terms <- sort(strsplit(fml, "+", fixed=TRUE)[[1]])
  old_terms <- sort(strsplit(old, "+", fixed=TRUE)[[1]])
  if (!identical(new_terms, old_terms))
    stop("base formula drifted from the hardcoded original:\n  missing: ",
         paste(setdiff(old_terms, new_terms), collapse=","),
         "\n  extra:   ", paste(setdiff(new_terms, old_terms), collapse=","))
  invisible(TRUE)
}
