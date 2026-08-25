#!/usr/bin/env Rscript

# HEAP Module 3 -- Generalized Mediation Analysis
#
# Reads:
#   HEAP.rds            -> heap_filter_exposures() -> covariates, disease data
#   Module 1 mediation_scores/<covarType>/<family>/mediation_scores/
#             mediation_scores_<idx>.txt  -> per-protein OOF score table
#
# Mediation modes:
#   primary_total              : G_total + PXS_total
#   partitioned_categories     : Gcis + Gtrans + PXS_<category> ...
#   partitioned_grouped_categories : Gcis + Gtrans + PXSgrp_<group> ...
#
# Usage:
#   Rscript Module3.R <idx> <split_num> <covarType> [<family>] [<mediation_mode>]
#   family        default = lasso
#   mediation_mode default = primary_total

local({
  candidates <- c(
    Sys.getenv("HEAP_PATHS_FILE", unset = ""),
    file.path(getwd(), "workflow", "00_paths.R"),
    file.path(getwd(), "..", "workflow", "00_paths.R"),
    file.path(getwd(), "..", "..", "workflow", "00_paths.R"),
    file.path(getwd(), "..", "..", "..", "workflow", "00_paths.R")
  )
  candidates <- candidates[nzchar(candidates)]
  hit <- candidates[file.exists(candidates)][1]
  if (!is.na(hit)) source(hit)
})

# Config helpers: load_sample_filter / apply_sample_filter (sample-filter axis).
# Sourced from the same workflow dir as 00_paths.R.
local({
  candidates <- c(
    if (exists("HEAP_PATHS") && !is.null(HEAP_PATHS$heap_root))
      file.path(HEAP_PATHS$heap_root, "workflow", "config_helpers.R") else "",
    file.path(getwd(), "workflow", "config_helpers.R"),
    file.path(getwd(), "..", "workflow", "config_helpers.R"),
    file.path(getwd(), "..", "..", "workflow", "config_helpers.R"),
    file.path(getwd(), "..", "..", "..", "workflow", "config_helpers.R")
  )
  candidates <- candidates[nzchar(candidates)]
  hit <- candidates[file.exists(candidates)][1]
  if (!is.na(hit)) source(hit)
})

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(pbapply)
  library(survival)
  library(tibble)
})

# -------------------------
# Global settings
# -------------------------
DELTA_SD_TOTAL      <- 1
PRINT_NA_COEF_MSG   <- FALSE
MIN_CASES_PER_DZ    <- 100L
MIN_EID_OVERLAP     <- 500L

prot_clean <- function(x) gsub("-", "_", x)

# ============================================================
# Time-to-event utilities (unchanged from original Module3)
# ============================================================
survival_time <- function(Time2Event_df, event_age, recode_status, recode_survtime) {
  T2E_df <- as.data.frame(Time2Event_df)
  T2E_df <- T2E_df[!is.na(T2E_df$recode_age_of_assessment_0_0), ]
  T2E_df <- T2E_df[
    is.na(T2E_df[[event_age]]) |
      (T2E_df$recode_age_of_assessment_0_0 < T2E_df[[event_age]]),
  ]

  censor_age <- pmin(
    T2E_df$recode_age_of_death_0_0,
    T2E_df$age_of_removal_0_0,
    T2E_df$age_of_lastfollowup,
    na.rm = TRUE
  )

  T2E_df[[recode_status]] <- as.integer(
    !is.na(T2E_df[[event_age]]) & (T2E_df[[event_age]] <= censor_age)
  )
  T2E_df[[recode_survtime]] <- ifelse(
    T2E_df[[recode_status]] == 1,
    T2E_df[[event_age]],
    censor_age
  )
  T2E_df
}

mediation_analysis_DZ <- function(DZ_df, DZ_ID, Orig_df) {
  DZ_df <- DZ_df %>%
    dplyr::select(all_of(c(
      "eid", DZ_ID,
      "recode_age_of_assessment_0_0",
      "recode_age_of_death_0_0",
      "age_of_removal_0_0",
      "age_of_lastfollowup"
    )))

  icd10_df <- survival_time(
    DZ_df,
    event_age       = DZ_ID,
    recode_status   = "DZ_status",
    recode_survtime = "DZ_survtime"
  )
  icd10_df$survtime_standard <-
    icd10_df$DZ_survtime - icd10_df$recode_age_of_assessment_0_0

  suppressMessages(dplyr::inner_join(Orig_df, icd10_df, by = "eid"))
}

zscore_numeric <- function(df, skip = c("eid")) {
  num_cols <- names(df)[vapply(df, is.numeric, logical(1))]
  num_cols <- setdiff(num_cols, skip)
  for (cc in num_cols) {
    x <- df[[cc]]
    s <- sd(x, na.rm = TRUE)
    if (is.na(s) || s == 0) {
      df[[cc]] <- 0
    } else {
      df[[cc]] <- as.numeric((x - mean(x, na.rm = TRUE)) / s)
    }
  }
  df
}

# ============================================================
# QC checks
# ============================================================

qc_required_cols <- function(score_df, mediation_mode) {
  # Raw PGS scores (G_raw, Gcis_raw, Gtrans_raw) are used for genetics because
  # lasso can penalize G_score to zero when exposed to competition from many
  # exposure features. Raw scores from OmicsPred are external GWAS-derived
  # predispositions that are not subject to this penalization.
  base <- c("eid", "protID", "protein_value", "G_raw", "Gcis_raw", "Gtrans_raw", "PXS_total")
  missing_base <- setdiff(base, names(score_df))
  if (length(missing_base) > 0)
    stop("Score table missing required columns: ", paste(missing_base, collapse = ", "),
         "\nExpected score table from Module1_suggested.R.")

  if (mediation_mode == "primary_total") {
    need <- c("G_raw", "PXS_total")
  } else if (mediation_mode == "partitioned_categories") {
    cat_cols <- grep("^PXS_", names(score_df), value = TRUE)
    cat_cols <- setdiff(cat_cols, "PXS_total")
    if (length(cat_cols) == 0)
      stop("partitioned_categories mode requires at least one PXS_<category> column beyond PXS_total.")
    need <- c("Gcis_raw", "Gtrans_raw")
  } else if (mediation_mode == "partitioned_grouped_categories") {
    grp_cols <- grep("^PXSgrp_", names(score_df), value = TRUE)
    if (length(grp_cols) == 0)
      stop("partitioned_grouped_categories mode requires at least one PXSgrp_<group> column.")
    need <- c("Gcis_raw", "Gtrans_raw")
  } else {
    stop("Unknown mediation_mode: ", mediation_mode)
  }

  missing_mode <- setdiff(need, names(score_df))
  if (length(missing_mode) > 0)
    stop("Mode '", mediation_mode, "' requires columns: ",
         paste(missing_mode, collapse = ", "), " (not found in score table).")
}

# ============================================================
# Predictor definition by mode
# ============================================================

define_predictors <- function(score_df, mediation_mode) {
  if (mediation_mode == "primary_total") {
    preds <- c("G_raw", "PXS_total")
    classes <- c(G_raw = "genetic_total", PXS_total = "exposure_total")
    return(list(predictors = preds, predictor_classes = classes))
  }

  if (mediation_mode == "partitioned_categories") {
    pxs_cat_cols <- grep("^PXS_", names(score_df), value = TRUE)
    pxs_cat_cols <- setdiff(pxs_cat_cols, "PXS_total")
    preds   <- c("Gcis_raw", "Gtrans_raw", pxs_cat_cols)
    classes <- c(
      Gcis_raw   = "genetic_cis",
      Gtrans_raw = "genetic_trans",
      setNames(rep("exposure_category", length(pxs_cat_cols)), pxs_cat_cols)
    )
    return(list(predictors = preds, predictor_classes = classes))
  }

  if (mediation_mode == "partitioned_grouped_categories") {
    pxs_grp_cols <- grep("^PXSgrp_", names(score_df), value = TRUE)
    preds   <- c("Gcis_raw", "Gtrans_raw", pxs_grp_cols)
    classes <- c(
      Gcis_raw   = "genetic_cis",
      Gtrans_raw = "genetic_trans",
      setNames(rep("exposure_group", length(pxs_grp_cols)), pxs_grp_cols)
    )
    return(list(predictors = preds, predictor_classes = classes))
  }

  stop("Unknown mediation_mode: ", mediation_mode)
}

# ============================================================
# Generalized model fitting
# ============================================================

fit_mediation_models_general <- function(df, predictors, covars) {
  names(df) <- gsub("-", "_", names(df))
  covars    <- intersect(covars, names(df))

  rhs_med <- paste(c(predictors, covars), collapse = " + ")
  m_form  <- as.formula(paste("protein_value ~", rhs_med))
  m_df    <- df[, c("protein_value", predictors, covars), drop = FALSE]
  m_df    <- m_df[stats::complete.cases(m_df), , drop = FALSE]
  m_fit   <- lm(m_form, data = m_df)

  rhs_out <- paste(c("protein_value", predictors, covars), collapse = " + ")
  o_form  <- as.formula(paste0("Surv(survtime_standard, DZ_status) ~ ", rhs_out))
  o_fit   <- coxph(o_form, data = df)

  list(
    mediator   = m_fit,
    outcome    = o_fit,
    predictors = predictors,
    covars     = covars,
    prot_col   = "protein_value"
  )
}

# ============================================================
# Generalized g-computation (average-over-individuals)
#
# For each focal predictor Xk (continuous):
#   low_k  = mean(Xk) - 0.5 * DELTA_SD_TOTAL * sd(Xk)
#   high_k = mean(Xk) + 0.5 * DELTA_SD_TOTAL * sd(Xk)
#
# NDE_k : change Xk low->high, mediator fixed at low-Xk counterfactual,
#          all other predictors at their observed values.
# NIE_k : hold Xk at high, change mediator from low- to high-Xk cf.
# ============================================================

.hi <- function(x, d = DELTA_SD_TOTAL) mean(x, na.rm = TRUE) + (d / 2) * sd(x, na.rm = TRUE)
.lo <- function(x, d = DELTA_SD_TOTAL) mean(x, na.rm = TRUE) - (d / 2) * sd(x, na.rm = TRUE)

make_newdata_replace <- function(fit, col, val) {
  mf   <- model.frame(fit)
  resp <- attr(terms(fit), "response")
  if (!is.null(resp) && resp > 0) mf <- mf[, -resp, drop = FALSE]
  mf[[col]] <- val
  as.data.frame(mf)
}

g_compute_nde_nie <- function(m_fit, o_fit, focal_col, delta_sd = DELTA_SD_TOTAL) {
  mf_m <- model.frame(m_fit)
  resp <- attr(terms(m_fit), "response")
  if (!is.null(resp) && resp > 0) mf_m <- mf_m[, -resp, drop = FALSE]

  x_obs <- mf_m[[focal_col]]
  lo_v  <- .lo(x_obs, delta_sd)
  hi_v  <- .hi(x_obs, delta_sd)

  # Counterfactual mediator predictions (only Xk changes, others observed)
  nd_m_lo <- make_newdata_replace(m_fit, focal_col, lo_v)
  nd_m_hi <- make_newdata_replace(m_fit, focal_col, hi_v)

  m_hat_lo <- as.numeric(predict(m_fit, newdata = nd_m_lo))
  m_hat_hi <- as.numeric(predict(m_fit, newdata = nd_m_hi))

  # Outcome predictions: NDE
  # lp(Xk=hi, M=m_hat_lo) - lp(Xk=lo, M=m_hat_lo)
  nd_o_hi_mlo <- make_newdata_replace(o_fit, focal_col, hi_v)
  nd_o_lo_mlo <- make_newdata_replace(o_fit, focal_col, lo_v)
  nd_o_hi_mlo[["protein_value"]] <- m_hat_lo
  nd_o_lo_mlo[["protein_value"]] <- m_hat_lo

  lp_hi_mlo <- as.numeric(predict(o_fit, newdata = nd_o_hi_mlo, type = "lp"))
  lp_lo_mlo <- as.numeric(predict(o_fit, newdata = nd_o_lo_mlo, type = "lp"))
  NDE_logHR <- mean(lp_hi_mlo - lp_lo_mlo, na.rm = TRUE)

  # Outcome predictions: NIE
  # lp(Xk=hi, M=m_hat_hi) - lp(Xk=hi, M=m_hat_lo)
  nd_o_hi_mhi <- make_newdata_replace(o_fit, focal_col, hi_v)
  nd_o_hi_mhi[["protein_value"]] <- m_hat_hi

  lp_hi_mhi <- as.numeric(predict(o_fit, newdata = nd_o_hi_mhi, type = "lp"))
  NIE_logHR <- mean(lp_hi_mhi - lp_hi_mlo, na.rm = TRUE)

  c(NDE = NDE_logHR, NIE = NIE_logHR)
}

med_effects_general <- function(fits, delta_sd = DELTA_SD_TOTAL) {
  out_rows <- list()

  for (focal in fits$predictors) {
    # Skip if column not in outcome model frame
    mf_o <- tryCatch(model.frame(fits$outcome), error = function(e) NULL)
    if (is.null(mf_o) || !focal %in% names(mf_o)) next

    eff <- tryCatch(
      g_compute_nde_nie(fits$mediator, fits$outcome, focal, delta_sd),
      error = function(e) {
        message("g_compute_nde_nie failed for focal=", focal, ": ", conditionMessage(e))
        c(NDE = NA_real_, NIE = NA_real_)
      }
    )

    out_rows[[length(out_rows) + 1L]] <- tibble::tibble(
      predictor     = focal,
      effect_type   = c("NDE",      "NIE"),
      effect_logHR  = c(eff["NDE"], eff["NIE"]),
      effect_HR     = exp(c(eff["NDE"], eff["NIE"])),
      contrast_sd_total = delta_sd
    )
  }

  if (length(out_rows) == 0) return(tibble::tibble())
  bind_rows(out_rows)
}

# ============================================================
# Generalized delta method
#
# For continuous predictor Xk:
#   delta_k = hi_k - lo_k  (scalar, same for all individuals since we
#             set Xk to a scalar value)
#
#   NDE_k = delta_k * beta_o[Xk]
#   NIE_k = beta_o[M] * delta_k * beta_m[Xk]
#
#   Var(NDE_k) = delta_k^2 * Var(beta_o[Xk])
#   Var(NIE_k) = delta_k^2 * (beta_m[Xk]^2 * Var(beta_o[M])
#                             + beta_o[M]^2 * Var(beta_m[Xk]))
#   (independence of m_fit and o_fit coefs assumed)
# ============================================================

delta_nde_nie_continuous <- function(m_fit, o_fit, focal_col,
                                     delta_sd = DELTA_SD_TOTAL) {
  b_m  <- coef(m_fit)
  b_o  <- coef(o_fit)
  V_m  <- vcov(m_fit)
  V_o  <- vcov(o_fit)

  if (!focal_col %in% names(b_m) || !focal_col %in% names(b_o)) {
    return(list(
      NDE_se = NA_real_, NDE_l95 = NA_real_, NDE_u95 = NA_real_, NDE_p = NA_real_,
      NIE_se = NA_real_, NIE_l95 = NA_real_, NIE_u95 = NA_real_, NIE_p = NA_real_
    ))
  }
  if (!"protein_value" %in% names(b_o)) {
    return(list(
      NDE_se = NA_real_, NDE_l95 = NA_real_, NDE_u95 = NA_real_, NDE_p = NA_real_,
      NIE_se = NA_real_, NIE_l95 = NA_real_, NIE_u95 = NA_real_, NIE_p = NA_real_
    ))
  }

  mf_m <- model.frame(m_fit)
  resp <- attr(terms(m_fit), "response")
  if (!is.null(resp) && resp > 0) mf_m <- mf_m[, -resp, drop = FALSE]

  x_obs   <- mf_m[[focal_col]]
  delta_k <- .hi(x_obs, delta_sd) - .lo(x_obs, delta_sd)

  beta_o_k <- b_o[[focal_col]]
  beta_o_M <- b_o[["protein_value"]]
  beta_m_k <- b_m[[focal_col]]

  # NDE point + SE
  NDE_est  <- delta_k * beta_o_k
  if (focal_col %in% rownames(V_o)) {
    NDE_var <- delta_k^2 * V_o[focal_col, focal_col]
  } else {
    NDE_var <- NA_real_
  }
  NDE_se <- sqrt(pmax(NDE_var, 0))

  # NIE point + SE (product rule, independent models)
  NIE_est  <- beta_o_M * delta_k * beta_m_k
  var_o_M  <- if ("protein_value" %in% rownames(V_o)) V_o["protein_value", "protein_value"] else NA_real_
  var_m_k  <- if (focal_col %in% rownames(V_m)) V_m[focal_col, focal_col] else NA_real_
  NIE_var  <- delta_k^2 * (beta_m_k^2 * var_o_M + beta_o_M^2 * var_m_k)
  NIE_se   <- sqrt(pmax(NIE_var, 0))

  make_ci <- function(est, se) {
    if (!is.finite(se)) return(list(l95 = NA_real_, u95 = NA_real_, p = NA_real_))
    list(
      l95 = est - 1.96 * se,
      u95 = est + 1.96 * se,
      p   = 2 * pnorm(-abs(est / se))
    )
  }

  nde_ci <- make_ci(NDE_est, NDE_se)
  nie_ci <- make_ci(NIE_est, NIE_se)

  list(
    NDE_se  = NDE_se,  NDE_l95 = nde_ci$l95, NDE_u95 = nde_ci$u95, NDE_p = nde_ci$p,
    NIE_se  = NIE_se,  NIE_l95 = nie_ci$l95, NIE_u95 = nie_ci$u95, NIE_p = nie_ci$p
  )
}

delta_mediation_general <- function(fits, point_est_df, delta_sd = DELTA_SD_TOTAL) {
  rows <- list()

  for (focal in fits$predictors) {
    d <- tryCatch(
      delta_nde_nie_continuous(fits$mediator, fits$outcome, focal, delta_sd),
      error = function(e) {
        message("delta method failed for focal=", focal, ": ", conditionMessage(e))
        list(
          NDE_se = NA_real_, NDE_l95 = NA_real_, NDE_u95 = NA_real_, NDE_p = NA_real_,
          NIE_se = NA_real_, NIE_l95 = NA_real_, NIE_u95 = NA_real_, NIE_p = NA_real_
        )
      }
    )

    for (et in c("NDE", "NIE")) {
      idx_row <- which(point_est_df$predictor == focal & point_est_df$effect_type == et)
      if (length(idx_row) == 0) next

      rows[[length(rows) + 1L]] <- tibble::tibble(
        predictor     = focal,
        effect_type   = et,
        delta_se      = d[[paste0(et, "_se")]],
        delta_l95     = d[[paste0(et, "_l95")]],
        delta_u95     = d[[paste0(et, "_u95")]],
        delta_p       = d[[paste0(et, "_p")]]
      )
    }
  }

  if (length(rows) == 0) return(tibble::tibble())
  bind_rows(rows)
}

# ============================================================
# Single disease runner
# ============================================================

run_one_disease_general <- function(p, DZ_ID, df, fits, pred_def,
                                    mediation_mode, delta_sd = DELTA_SD_TOTAL) {
  eff_df <- tryCatch(
    med_effects_general(fits, delta_sd),
    error = function(e) {
      message("med_effects_general failed: prot=", p, " DZ=", DZ_ID, " ", conditionMessage(e))
      tibble::tibble()
    }
  )

  if (nrow(eff_df) == 0) return(NULL)

  delta_df <- delta_mediation_general(fits, eff_df, delta_sd)

  if (nrow(delta_df) > 0) {
    eff_df <- suppressMessages(dplyr::left_join(eff_df, delta_df, by = c("predictor", "effect_type")))
  } else {
    eff_df$delta_se  <- NA_real_
    eff_df$delta_l95 <- NA_real_
    eff_df$delta_u95 <- NA_real_
    eff_df$delta_p   <- NA_real_
  }

  eff_df$delta_l95_HR <- exp(eff_df$delta_l95)
  eff_df$delta_u95_HR <- exp(eff_df$delta_u95)

  # Protein HR from Cox model
  sm  <- summary(fits$outcome)
  rn  <- rownames(sm$coefficients)
  prot_col <- fits$prot_col

  protein_HR <- protein_HR_l95 <- protein_HR_u95 <- protein_p <- NA_real_
  if (prot_col %in% rn) {
    beta_p <- sm$coefficients[prot_col, "coef"]
    protein_HR <- exp(beta_p)
    ci <- tryCatch(
      exp(confint(fits$outcome))[prot_col, ],
      error = function(e) c(NA_real_, NA_real_)
    )
    protein_HR_l95 <- ci[1]
    protein_HR_u95 <- ci[2]
    protein_p      <- sm$coefficients[prot_col, "Pr(>|z|)"]
  }

  mediator_adjR2 <- summary(fits$mediator)$adj.r.squared
  cox_cindex     <- fits$outcome$concordance["concordance"]
  cox_cindex_se  <- fits$outcome$concordance["std"]
  n              <- nrow(df)
  n_cases        <- sum(df$DZ_status, na.rm = TRUE)

  eff_df %>%
    mutate(
      protID          = p,
      DZ_ID           = DZ_ID,
      mediation_mode  = mediation_mode,
      predictor_class = pred_def$predictor_classes[predictor],
      protein_HR      = protein_HR,
      protein_HR_l95  = protein_HR_l95,
      protein_HR_u95  = protein_HR_u95,
      protein_p       = protein_p,
      mediator_adjR2  = mediator_adjR2,
      cox_cindex      = cox_cindex,
      cox_cindex_se   = cox_cindex_se,
      n               = n,
      n_cases         = n_cases
    ) %>%
    dplyr::select(
      protID, DZ_ID, mediation_mode,
      predictor, predictor_class, effect_type,
      effect_logHR, effect_HR,
      delta_se, delta_l95, delta_u95, delta_p,
      delta_l95_HR, delta_u95_HR,
      contrast_sd_total,
      protein_HR, protein_HR_l95, protein_HR_u95, protein_p,
      mediator_adjR2, cox_cindex, cox_cindex_se,
      n, n_cases
    )
}

# ============================================================
# Per-protein runner
# ============================================================

run_protein <- function(p, score_df, DZ_df, DZ_ids, covars_df, covars_list,
                        mediation_mode, delta_sd = DELTA_SD_TOTAL) {
  score_p <- score_df[score_df$protID == prot_clean(p), ]

  if (nrow(score_p) == 0) {
    message("No score rows for protein: ", p)
    return(NULL)
  }

  # Merge: scores <- covariates <- (disease handled per DZ below)
  base_df <- suppressMessages(
    score_p %>%
      dplyr::inner_join(covars_df, by = "eid")
  )

  if (nrow(base_df) < MIN_EID_OVERLAP) {
    message("eid overlap too small for protein: ", p, " (n=", nrow(base_df), ")")
    return(NULL)
  }

  base_df <- zscore_numeric(base_df, skip = c("eid"))

  pred_def <- define_predictors(score_p, mediation_mode)

  # Per-protein genetic-instrument availability, from Module1's gcis_present /
  # gtrans_present flags (default TRUE for older score tables without them). A
  # protein with no cis/trans variants has an all-zero Gcis_raw/Gtrans_raw: not a
  # valid instrument. Left in the model it is a constant (zscore_numeric maps it to
  # 0), perfectly collinear -> aliased NA coefficient and a meaningless 0 NDE/NIE.
  # Drop the absent genetic predictor(s) and carry the flags onto the output so the
  # summary reports "no instrument" rather than a spurious 0 estimate.
  .is_absent <- function(x) length(x) > 0 && !is.na(x[1]) &&
    (isFALSE(x[1]) || identical(tolower(as.character(x[1])), "false"))
  gcis_present   <- !.is_absent(score_p$gcis_present)
  gtrans_present <- !.is_absent(score_p$gtrans_present)
  pc <- pred_def$predictor_classes
  drop_genetic <- character(0)
  if (!gcis_present)                   drop_genetic <- c(drop_genetic, names(pc)[pc == "genetic_cis"])
  if (!gtrans_present)                 drop_genetic <- c(drop_genetic, names(pc)[pc == "genetic_trans"])
  if (!gcis_present && !gtrans_present) drop_genetic <- c(drop_genetic, names(pc)[pc == "genetic_total"])
  drop_genetic <- unique(drop_genetic)
  if (length(drop_genetic) > 0)
    message("Protein ", p, ": no ",
            paste(c(if (!gcis_present) "cis", if (!gtrans_present) "trans"), collapse = "/"),
            " genetic instrument; dropping predictor(s): ", paste(drop_genetic, collapse = ", "))

  predictors_needed <- setdiff(pred_def$predictors, drop_genetic)
  missing_preds <- setdiff(predictors_needed, names(base_df))
  if (length(missing_preds) > 0) {
    message("Missing predictor columns for ", p, ": ", paste(missing_preds, collapse = ", "))
    return(NULL)
  }

  results <- lapply(DZ_ids, function(DZ_ID) {
    tryCatch({
      df_dz <- mediation_analysis_DZ(DZ_df, DZ_ID, base_df)
      df_dz <- df_dz[stats::complete.cases(
        df_dz[, c("protein_value", "survtime_standard", "DZ_status",
                  predictors_needed, covars_list), drop = FALSE]
      ), ]

      if (nrow(df_dz) == 0) return(NULL)
      if (sum(df_dz$DZ_status, na.rm = TRUE) < MIN_CASES_PER_DZ) return(NULL)

      fits <- fit_mediation_models_general(df_dz, predictors_needed, covars_list)

      if (PRINT_NA_COEF_MSG) {
        na_m <- sum(is.na(coef(fits$mediator)))
        na_o <- sum(is.na(coef(fits$outcome)))
        if (na_m > 0 || na_o > 0)
          message("NA coefs: prot=", p, " DZ=", DZ_ID,
                  " lm_na=", na_m, " cox_na=", na_o)
      }

      run_one_disease_general(p, DZ_ID, df_dz, fits, pred_def,
                              mediation_mode, delta_sd)
    }, error = function(e) {
      message("Error prot=", p, " DZ=", DZ_ID, ": ", conditionMessage(e))
      NULL
    })
  })

  res <- bind_rows(Filter(Negate(is.null), results))
  if (nrow(res) > 0) {
    # Per-protein genetic-instrument availability on every row.
    res$gcis_present       <- gcis_present
    res$gtrans_present     <- gtrans_present
    # TRUE = effect was genuinely estimated from a valid instrument.
    res$instrument_present <- TRUE

    # An absent cis/trans component has no valid instrument, so its direct (NDE)
    # and indirect (NIE) effects are a structural 0 (logHR = 0, HR = 1). Emit
    # explicit labelled-0 rows per disease, flagged instrument_present = FALSE, so
    # the component is visible in the summary rather than silently missing —
    # mirroring Module1's "r2 = 0 + present = FALSE" convention. SE/CI/p stay NA
    # (nothing was fit).
    if (length(drop_genetic) > 0) {
      add <- list()
      for (dz in unique(res$DZ_ID)) {
        tmpl <- res[res$DZ_ID == dz, , drop = FALSE][1, , drop = FALSE]
        for (gp in drop_genetic) {
          for (et in c("NDE", "NIE")) {
            r <- tmpl
            r$predictor       <- gp
            r$predictor_class <- unname(pc[gp])
            r$effect_type     <- et
            r$effect_logHR    <- 0
            r$effect_HR       <- 1
            r$delta_se <- NA_real_; r$delta_l95 <- NA_real_; r$delta_u95 <- NA_real_
            r$delta_p  <- NA_real_; r$delta_l95_HR <- NA_real_; r$delta_u95_HR <- NA_real_
            r$instrument_present <- FALSE
            add[[length(add) + 1L]] <- r
          }
        }
      }
      if (length(add) > 0) res <- dplyr::bind_rows(res, add)
    }
  }
  res
}

# ============================================================
# CLI + loading
# ============================================================

# ============================================================
# CLI: manifest-driven or legacy positional arguments
#
# Manifest mode (preferred):
#   Rscript Module3.R --manifest <path> --array-index <N>
#
# Legacy positional mode (backward compatibility):
#   Rscript Module3.R <idx> <split_num> <covarType> [<family>] [<mediation_mode>]
# ============================================================

args <- commandArgs(trailingOnly = TRUE)

.parse_flag_m3 <- function(args, flag, default = NULL) {
  i <- which(args == flag)
  if (length(i) == 0 || i[1] >= length(args)) return(default)
  args[i[1] + 1L]
}

.is_manifest_mode_m3 <- length(args) >= 2 && args[1] == "--manifest"

allowed_families <- c("lasso", "ridge", "enet")
allowed_modes    <- c("primary_total", "partitioned_categories",
                      "partitioned_grouped_categories")

if (.is_manifest_mode_m3) {
  .manifest_path <- .parse_flag_m3(args, "--manifest")
  .array_idx_str <- .parse_flag_m3(args, "--array-index")

  if (is.null(.manifest_path))
    stop("--manifest <path> is required")
  if (is.null(.array_idx_str))
    stop("--array-index <N> is required in manifest mode")
  if (!file.exists(.manifest_path))
    stop("Manifest not found: ", .manifest_path)

  .mf      <- utils::read.delim(.manifest_path, stringsAsFactors = FALSE, check.names = FALSE)
  .arr_idx <- as.integer(.array_idx_str)
  .row     <- .mf[.mf$array_index == .arr_idx, , drop = FALSE]

  if (nrow(.row) == 0)
    stop("No manifest row for array_index=", .arr_idx, " in ", .manifest_path)
  if (nrow(.row) > 1)
    stop("Duplicate array_index=", .arr_idx, " in ", .manifest_path)
  .row <- as.list(.row[1L, ])

  idx            <- as.integer(.row$chunk_id)
  split_num      <- as.integer(.row$n_chunks)
  covarType      <- as.character(.row$covariate_set)
  family_type    <- as.character(.row$family)
  mediation_mode <- as.character(.row$mediation_mode)
  experiment_name <- as.character(.row$experiment_name %||% "unknown")
  .sample_filter  <- as.character(.row$sample_filter %||% "none")

  # Optional manifest overrides
  if (!is.null(.row$min_cases_per_disease) && !is.na(.row$min_cases_per_disease))
    MIN_CASES_PER_DZ  <- as.integer(.row$min_cases_per_disease)
  if (!is.null(.row$min_eid_overlap) && !is.na(.row$min_eid_overlap))
    MIN_EID_OVERLAP   <- as.integer(.row$min_eid_overlap)
  if (!is.null(.row$delta_sd_total) && !is.na(.row$delta_sd_total))
    DELTA_SD_TOTAL    <- as.numeric(.row$delta_sd_total)

  # Upstream score directory from manifest (IGLOO-rooted)
  .score_dir_from_manifest <- if (!is.null(.row$upstream_score_dir) &&
                                   !is.na(.row$upstream_score_dir) &&
                                   nzchar(.row$upstream_score_dir)) {
    as.character(.row$upstream_score_dir)
  } else NULL

  # Output path from manifest (IGLOO-rooted)
  .out_root_from_manifest <- if (!is.null(.row$output_path) &&
                                   !is.na(.row$output_path) &&
                                   nzchar(.row$output_path)) {
    as.character(.row$output_path)
  } else NULL

  message(sprintf("[manifest] experiment=%s  array_index=%d  chunk=%d/%d  covar=%s  family=%s  mode=%s",
                  experiment_name, .arr_idx, idx, split_num, covarType, family_type, mediation_mode))

} else {
  # ---- Legacy positional mode ----
  if (length(args) < 3) {
    stop(
      "Usage (manifest):   Rscript Module3.R --manifest <path> --array-index <N>\n",
      "Usage (positional): Rscript Module3.R <idx> <split_num> <covarType> [<family>] [<mediation_mode>]\n",
      "  family        : lasso | ridge | enet  (default: lasso)\n",
      "  mediation_mode: primary_total | partitioned_categories | partitioned_grouped_categories\n",
      "                  (default: primary_total)"
    )
  }
  idx            <- as.integer(args[1])
  split_num      <- as.integer(args[2])
  covarType      <- as.character(args[3])
  family_type    <- if (length(args) >= 4) as.character(args[4]) else "lasso"
  mediation_mode <- if (length(args) >= 5) as.character(args[5]) else "primary_total"
  experiment_name <- paste0("positional_", covarType, "_", family_type, "_", mediation_mode)
  .sample_filter  <- "none"
  .score_dir_from_manifest <- NULL
  .out_root_from_manifest  <- NULL
}

if (!family_type %in% allowed_families)
  stop("family must be one of: ", paste(allowed_families, collapse = ", "))
if (!mediation_mode %in% allowed_modes)
  stop("mediation_mode must be one of: ", paste(allowed_modes, collapse = ", "))

# ============================================================
# Load HEAP
# ============================================================

if (!exists("heap_filter_exposures", mode = "function"))
  stop("heap_filter_exposures() not found. Ensure 00_paths.R is sourced correctly.")

heap <- readRDS(heap_loader_rds)
heap <- heap_filter_exposures(heap)

DZ_df     <- heap$disease$DZ_df
DZ_ids    <- heap$disease$DZ_ids
covars_df <- heap$covars_baseline

if (is.null(DZ_df) || is.null(DZ_ids) || length(DZ_ids) == 0)
  stop("Disease object is missing or empty. Was HEAP built with HEAP_SKIP_DISEASE=TRUE?")

message("Loaded HEAP: ", length(DZ_ids), " diseases, ",
        nrow(covars_df), " subjects with covariates.")

# ============================================================
# Resolve covariate set
# Prefer centralized config (covariate_sets.yml); fall back to inline.
# ============================================================

.resolve_covartype_m3 <- function(covarType, heap) {
  .cfg_path <- heap_config("covariates", "covariate_sets.yml")
  if (file.exists(.cfg_path) && requireNamespace("yaml", quietly = TRUE)) {
    tryCatch({
      .sets <- yaml::read_yaml(.cfg_path)$covariate_sets
      if (covarType %in% names(.sets)) {
        .entry <- .sets[[covarType]]
        .covars <- .entry$covariates
        if (is.null(.covars) || (length(.covars) == 1L && is.na(.covars[[1L]]))) {
          return(heap$covars_list)
        }
        return(as.character(.covars))
      }
    }, error = function(e) NULL)
  }
  # Inline fallback (used only if covariate_sets.yml is unreadable). Mirrors the
  # descriptive sets in covariate_sets.yml EXACTLY (base = PRIMARY; others are
  # supplementary layers on base). No BMI/fasting in base.
  .pcs  <- paste0("genetic_principal_components_f22009_0_", 1:20)
  .base <- c("age_when_attended_assessment_centre_f21003_0_0", "sex_f31_0_0",
             "age2", "age_sex", "age2_sex",
             "uk_biobank_assessment_centre_f54_0_0", .pcs)
  CovarSpec <- list(
    base          = .base,
    base_bmi      = c(.base, "body_mass_index_bmi_f23104_0_0"),
    base_draw     = c(.base, "fasting_time_f74_0_0", "assessment_season"),
    base_clinical = c(.base, "body_mass_index_bmi_f23104_0_0", "fasting_time_f74_0_0",
                      "assessment_season",
                      "combined_Blood_pressure_medication",
                      "combined_Hormone_replacement_therapy",
                      "combined_Oral_contraceptive_pill_or_minipill",
                      "combined_Insulin",
                      "combined_Cholesterol_lowering_medication"),
    base_ses       = .base,
    base_prevalent = c(.base, "prevalent_major_disease")
  )
  if (!covarType %in% names(CovarSpec))
    stop("covarType '", covarType, "' not found in covariate_sets.yml or inline CovarSpec.\n",
         "Available: ", paste(names(CovarSpec), collapse = ", "))
  CovarSpec[[covarType]]
}

covars_vec  <- .resolve_covartype_m3(covarType, heap)
covars_list <- intersect(covars_vec, names(covars_df))

# Sample filter (the WHO-IS-IN axis), applied to the FULL covariate frame BEFORE
# the covariate-column narrowing below (otherwise the filter column may already be
# gone). The model frame is built by inner joins on covars_df (see run_protein),
# so restricting covars_df rows here restricts the analysis sample for the run.
.sf_spec <- load_sample_filter(.sample_filter)
if (!is.null(.sf_spec)) {
  .sf <- apply_sample_filter(covars_df, .sf_spec)
  covars_df <- .sf$df
  message(sprintf("[sample_filter] %s: %d -> %d (dropped %d)",
                  .sf_spec$name, .sf$n_before, .sf$n_after, .sf$n_dropped))
}

covars_df   <- covars_df[, c("eid", covars_list), drop = FALSE]

# ============================================================
# Resolve upstream Module 1 score directory
# ============================================================

score_base_dir <- if (!is.null(.score_dir_from_manifest)) {
  # IGLOO-rooted path from manifest
  .score_dir_from_manifest
} else {
  # Legacy: derive from local heap_output pattern (backward compat)
  file.path(
    heap_output("module1_predictive_r2_score_partition"),
    covarType, family_type, "mediation_scores"
  )
}

if (!dir.exists(score_base_dir))
  stop("Module 1 mediation score directory not found: ", score_base_dir,
       "\nRun Module1_suggested.R first with covarType=", covarType,
       " family=", family_type,
       "\nIf using manifest mode, check upstream_score_dir in the manifest.")

score_files <- list.files(score_base_dir, pattern = "^mediation_scores_.*\\.txt$",
                          full.names = TRUE)
if (length(score_files) == 0)
  stop("No mediation score files found in: ", score_base_dir)

message("Loading ", length(score_files), " score file(s) from ", score_base_dir)
score_df <- rbindlist(lapply(score_files, fread), fill = TRUE) %>% as.data.frame()
score_df$protID <- prot_clean(score_df$protID)

# QC: validate required columns exist for chosen mediation mode
qc_required_cols(score_df, mediation_mode)

# Check eid overlap
overlap <- length(intersect(score_df$eid, covars_df$eid))
if (overlap < MIN_EID_OVERLAP)
  stop("eid overlap between score table and covariates is too small: n=", overlap)

all_prots <- unique(score_df$protID)
message("Score table: ", nrow(score_df), " rows, ", length(all_prots), " proteins, ",
        "mediation_mode=", mediation_mode)

# ============================================================
# Split proteins into batches
# ============================================================

if (!is.finite(split_num) || split_num < 1)
  stop("split_num must be a positive integer")

split_vectors <- if (split_num >= length(all_prots)) {
  unname(as.list(all_prots))
} else {
  groups <- cut(seq_along(all_prots), breaks = split_num, labels = FALSE)
  split(all_prots, groups)
}

if (idx < 1 || idx > length(split_vectors))
  stop("idx=", idx, " out of range [1, ", length(split_vectors), "]")

prot_batch <- split_vectors[[idx]]
message("Processing batch idx=", idx, ": ", length(prot_batch), " proteins")

# ============================================================
# Resolve output directory
# ============================================================

out_dir <- if (!is.null(.out_root_from_manifest)) {
  # IGLOO-rooted; experiment_name already in path from manifest
  file.path(.out_root_from_manifest, covarType, family_type, mediation_mode)
} else {
  # Non-manifest default: IGLOO canonical (was repo-local heap_output).
  file.path(heap_project_output("module3"), covarType, family_type, mediation_mode)
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Guard: output must not be under scratch
if (grepl(HEAP_PATHS$scratch_root, normalizePath(out_dir, mustWork = FALSE), fixed = TRUE))
  stop("REPRODUCIBILITY VIOLATION: Module3 out_dir points to scratch: ", out_dir)

# ============================================================
# Write resolved run config artifact
# ============================================================

tryCatch({
  .rc <- list(
    experiment_name  = experiment_name,
    module           = "module3",
    array_index      = if (.is_manifest_mode_m3) .arr_idx else NA_integer_,
    chunk_id         = idx,
    n_chunks         = split_num,
    covariate_set    = covarType,
    family           = family_type,
    mediation_mode   = mediation_mode,
    delta_sd_total   = DELTA_SD_TOTAL,
    min_cases_per_dz = MIN_CASES_PER_DZ,
    min_eid_overlap  = MIN_EID_OVERLAP,
    score_base_dir   = score_base_dir,
    out_dir          = out_dir,
    heap_rds         = heap_loader_rds,
    manifest_path    = if (.is_manifest_mode_m3) .manifest_path else NA_character_,
    run_timestamp    = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    run_host         = Sys.info()[["nodename"]]
  )
  if (requireNamespace("yaml", quietly = TRUE)) {
    yaml::write_yaml(.rc, file.path(out_dir, paste0("run_config_", idx, ".yml")))
  } else {
    saveRDS(.rc, file.path(out_dir, paste0("run_config_", idx, ".rds")))
  }
}, error = function(e) {
  warning("Could not write run_config: ", conditionMessage(e))
})

# ============================================================
# Run mediation
# ============================================================

MDres_list <- pblapply(prot_batch, function(p) {
  tryCatch(
    run_protein(
      p              = p,
      score_df       = score_df,
      DZ_df          = DZ_df,
      DZ_ids         = DZ_ids,
      covars_df      = covars_df,
      covars_list    = covars_list,
      mediation_mode = mediation_mode
    ),
    error = function(e) {
      message("run_protein failed for ", p, ": ", conditionMessage(e))
      NULL
    }
  )
})

MDfinal <- bind_rows(Filter(Negate(is.null), MDres_list))

if (nrow(MDfinal) == 0) {
  message("No results produced for idx=", idx, " (all proteins skipped or no cases).")
} else {
  required_out <- c("protID", "DZ_ID", "mediation_mode",
                    "predictor", "predictor_class", "effect_type",
                    "effect_logHR", "effect_HR")
  missing_out <- setdiff(required_out, names(MDfinal))
  if (length(missing_out) > 0)
    warning("Output table missing columns: ", paste(missing_out, collapse = ", "))
}

out_file <- file.path(out_dir, paste0("MDres_", idx, ".txt"))
fwrite(as.data.table(MDfinal), file = out_file, sep = "\t")
message("Saved: ", out_file, "  (", nrow(MDfinal), " rows)")
