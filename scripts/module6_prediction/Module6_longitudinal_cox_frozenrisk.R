#!/usr/bin/env Rscript

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

# Frozen baseline Cox risk scoring for Module 6 longitudinal PES outputs.
# Fits disease-risk Cox models in baseline training data using OOF PES scores,
# then applies the same frozen coefficients to held-out i0/i2/i3 visits.

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(survival)
  library(tibble)
  library(progress)
})

setDTthreads(1)

cfg <- list(
  out_dir = heap_project_output("module6_pes_longitudinal"),
  paths = list(
    # Reads HEAP.rds directly; derives the longitudinal PXS via as_pxs_longitudinal()
    # and the disease time-to-event table from its $disease$DZ_df slot (unified loader;
    # no separate legacy UKB_MDstore dependency).
    heap_rds = heap_loader_rds
  ),
  # Mirrors config/covariates/covariate_sets.yml (descriptive base-primary scheme).
  CovarSpec = list(
    base = c(
      "age_when_attended_assessment_centre_f21003_0_0", "sex_f31_0_0",
      "age2", "age_sex", "age2_sex",
      "uk_biobank_assessment_centre_f54_0_0",
      paste0("genetic_principal_components_f22009_0_", 1:20)
    ),
    base_bmi = c(
      "age_when_attended_assessment_centre_f21003_0_0", "sex_f31_0_0",
      "age2", "age_sex", "age2_sex",
      "uk_biobank_assessment_centre_f54_0_0",
      paste0("genetic_principal_components_f22009_0_", 1:20),
      "body_mass_index_bmi_f23104_0_0"
    ),
    base_draw = c(
      "age_when_attended_assessment_centre_f21003_0_0", "sex_f31_0_0",
      "age2", "age_sex", "age2_sex",
      "uk_biobank_assessment_centre_f54_0_0",
      paste0("genetic_principal_components_f22009_0_", 1:20),
      "fasting_time_f74_0_0", "assessment_season"
    ),
    base_clinical = c(
      "age_when_attended_assessment_centre_f21003_0_0", "sex_f31_0_0",
      "age2", "age_sex", "age2_sex",
      "uk_biobank_assessment_centre_f54_0_0",
      paste0("genetic_principal_components_f22009_0_", 1:20),
      "body_mass_index_bmi_f23104_0_0", "fasting_time_f74_0_0", "assessment_season",
      "combined_Blood_pressure_medication",
      "combined_Hormone_replacement_therapy",
      "combined_Oral_contraceptive_pill_or_minipill",
      "combined_Insulin",
      "combined_Cholesterol_lowering_medication"
    ),
    base_prevalent = c(
      "age_when_attended_assessment_centre_f21003_0_0", "sex_f31_0_0",
      "age2", "age_sex", "age2_sex",
      "uk_biobank_assessment_centre_f54_0_0",
      paste0("genetic_principal_components_f22009_0_", 1:20),
      "prevalent_major_disease"
    ),
    # exclude_prevalent SENSITIVITY: `base` covariates on the healthy-at-baseline
    # subset (sample filter drops prevalent_major_disease==1). Mirrors
    # Module6_prod_longitudinal.R; the frozen-risk Cox must run on the SAME subset.
    base_exclprev = c(
      "age_when_attended_assessment_centre_f21003_0_0", "sex_f31_0_0",
      "age2", "age_sex", "age2_sex",
      "uk_biobank_assessment_centre_f54_0_0",
      paste0("genetic_principal_components_f22009_0_", 1:20)
    )
  )
)

# covar_types whose cohort excludes prevalent major disease at baseline.
EXCLUDE_PREVALENT_COVAR_TYPES <- c("base_exclprev")

message_ts <- function(...) {
  message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), ...)
}

safe_name <- function(x) gsub("[^A-Za-z0-9_.-]+", "_", x)
bt <- function(x) paste0("`", gsub("`", "\\\\`", x), "`")

scalar_num <- function(x, default = NA_real_) {
  if (is.null(x) || length(x) == 0) return(default)
  out <- suppressWarnings(as.numeric(x[[1]]))
  ifelse(is.na(out), default, out)
}

append_table <- function(x, path) {
  if (is.null(x) || nrow(x) == 0) return(invisible(FALSE))
  append <- file.exists(path) && file.info(path)$size > 0
  fwrite(as.data.table(x), path, sep = "\t", append = append, col.names = !append)
  invisible(TRUE)
}

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 2) {
    stop(paste0(
      "Usage:\n",
      "  Rscript Module6_longitudinal_cox_frozenrisk.R <covarType> <exposure_id> [options]\n\n",
      "Options:\n",
      "  --result-dir <dir>       Default: <out_dir>/<covarType>\n",
      "  --instances 0,2,3        Held-out visit instances to score\n",
      "  --score-types prot,full  PES score types to run\n",
      "  --min-n 200              Minimum training sample size\n",
      "  --min-events 20          Minimum training incident events\n",
      "  --eval-min-events 5      Minimum evaluation events for c-index reporting\n",
      "  --max-diseases N         Optional pilot cap\n",
      "  --disease-regex REGEX    Optional disease-age-column filter\n",
      "  --output-tag TAG         Optional tag inserted after FrozenCox in output names\n",
      "  --pes-scale final|oof    Holdout PES scale. Default: final, use oof for strict OOF-Cox scale matching\n",
      "  --horizons 1,3,5,10      Fixed-year absolute-risk horizons\n",
      "  --risk-horizons 1,3,5,10 Legacy alias for --horizons\n",
      "  --save-risk-scores       Save deployable per-person absolute risk scores\n",
      "  --save-person-scores     Save per-person visit risk scores\n",
      "  --save-models            Save fitted Cox model RDS; memory-heavy for full runs\n"
    ))
  }

  out <- list(
    covar_type = trimws(args[[1]]),  # guard against stray trailing whitespace in COVAR_TYPE
    exposure_id = args[[2]],
    result_dir = NA_character_,
    instances = c(0L, 2L, 3L),
    score_types = c("prot", "full"),
    min_n = 200L,
    min_events = 20L,
    eval_min_events = 5L,
    max_diseases = NA_integer_,
    disease_regex = NA_character_,
    output_tag = "",
    pes_scale = "final",
    risk_horizons = c(1, 3, 5, 10),
    save_risk_scores = FALSE,
    save_person_scores = FALSE,
    save_models = FALSE
  )

  i <- 3
  while (i <= length(args)) {
    key <- args[[i]]
    val <- if (i + 1 <= length(args)) args[[i + 1]] else NA_character_
    if (key == "--result-dir") {
      out$result_dir <- val
      i <- i + 2
    } else if (key == "--instances") {
      out$instances <- as.integer(strsplit(val, ",", fixed = TRUE)[[1]])
      i <- i + 2
    } else if (key == "--score-types") {
      out$score_types <- strsplit(val, ",", fixed = TRUE)[[1]]
      i <- i + 2
    } else if (key == "--min-n") {
      out$min_n <- as.integer(val)
      i <- i + 2
    } else if (key == "--min-events") {
      out$min_events <- as.integer(val)
      i <- i + 2
    } else if (key == "--eval-min-events") {
      out$eval_min_events <- as.integer(val)
      i <- i + 2
    } else if (key == "--max-diseases") {
      out$max_diseases <- as.integer(val)
      i <- i + 2
    } else if (key == "--disease-regex") {
      out$disease_regex <- val
      i <- i + 2
    } else if (key == "--output-tag") {
      out$output_tag <- safe_name(val)
      i <- i + 2
    } else if (key == "--pes-scale") {
      out$pes_scale <- val
      i <- i + 2
    } else if (key %in% c("--horizons", "--risk-horizons")) {
      out$risk_horizons <- suppressWarnings(as.numeric(strsplit(val, ",", fixed = TRUE)[[1]]))
      i <- i + 2
    } else if (key == "--save-risk-scores") {
      out$save_risk_scores <- TRUE
      i <- i + 1
    } else if (key == "--save-person-scores") {
      out$save_person_scores <- TRUE
      i <- i + 1
    } else if (key == "--save-models") {
      out$save_models <- TRUE
      i <- i + 1
    } else if (key == "--no-save-models") {
      out$save_models <- FALSE
      i <- i + 1
    } else {
      stop("Unknown option: ", key)
    }
  }

  out$score_types <- intersect(out$score_types, c("prot", "full"))
  if (length(out$score_types) == 0) stop("--score-types must include prot and/or full")
  if (!out$pes_scale %in% c("final", "oof")) stop("--pes-scale must be final or oof")
  out$risk_horizons <- out$risk_horizons[is.finite(out$risk_horizons) & out$risk_horizons > 0]
  if (length(out$risk_horizons) == 0) stop("--horizons/--risk-horizons must include positive numeric years")
  out
}

as_long_pxs <- function(x) {
  if (is.list(x) && !isS4(x)) return(x)
  if (!isS4(x)) stop("PXS object must be S4 or list")
  list(
    covars_df = x@covars_df,
    covars_list = x@covars_list
  )
}

read_score_table <- function(result_dir, covar_type, exposure_id, suffix) {
  prefix <- file.path(result_dir, paste0("PESlong_", covar_type, "_", safe_name(exposure_id)))
  rds <- paste0(prefix, "_", suffix, ".rds")
  tsv <- paste0(prefix, "_", suffix, ".tsv")
  if (file.exists(tsv)) return(as.data.frame(fread(tsv)))
  if (file.exists(rds)) return(as.data.frame(readRDS(rds)))
  stop("Missing ", suffix, " artifact for exposure: ", exposure_id, "\nTried:\n  ", rds, "\n  ", tsv)
}

# Disease time-to-event table, sourced from the unified HEAP.rds $disease$DZ_df slot
# (built by HEAP_loader.R::load_disease()). Disease outcomes are covariate-scheme
# independent, so this no longer depends on covar_type or the legacy UKB_MDstore.
load_t2e <- function(heap) {
  dz <- heap$disease$DZ_df
  if (is.null(dz)) {
    stop("HEAP.rds has no $disease$DZ_df (was it built with HEAP_SKIP_DISEASE=TRUE?)")
  }
  out <- as.data.frame(dz)
  needed <- c(
    "eid", "recode_age_of_assessment_0_0", "recode_age_of_death_0_0",
    "age_of_removal_0_0", "age_of_lastfollowup"
  )
  miss <- setdiff(needed, names(out))
  if (length(miss) > 0) stop("DZ_df missing columns: ", paste(miss, collapse = ", "))
  out
}

make_censor_age <- function(df) {
  censor_age <- pmin(
    suppressWarnings(as.numeric(df$recode_age_of_death_0_0)),
    suppressWarnings(as.numeric(df$age_of_removal_0_0)),
    suppressWarnings(as.numeric(df$age_of_lastfollowup)),
    na.rm = TRUE
  )
  censor_age[!is.finite(censor_age)] <- NA_real_
  censor_age
}

survival_from_landmark <- function(df, event_age_col, landmark_age_col = "landmark_age") {
  out <- as.data.frame(df)
  out$event_age_col <- event_age_col
  out$landmark_age <- suppressWarnings(as.numeric(out[[landmark_age_col]]))
  out$event_age <- suppressWarnings(as.numeric(out[[event_age_col]]))
  out$censor_age <- make_censor_age(out)
  out <- out[is.finite(out$landmark_age) & is.finite(out$censor_age), , drop = FALSE]
  out <- out[out$censor_age > out$landmark_age, , drop = FALSE]
  out <- out[is.na(out$event_age) | out$event_age > out$landmark_age, , drop = FALSE]
  out$DZ_status <- as.integer(!is.na(out$event_age) & out$event_age <= out$censor_age)
  out$DZ_survtime <- ifelse(out$DZ_status == 1, out$event_age, out$censor_age)
  out$survtime_standard <- out$DZ_survtime - out$landmark_age
  out <- out[is.finite(out$survtime_standard) & out$survtime_standard > 0, , drop = FALSE]
  out
}

format_exposure_for_cox <- function(y_raw, exposure_type, train_levels = NULL) {
  if (exposure_type == "continuous") return(suppressWarnings(as.numeric(y_raw)))
  out <- as.factor(y_raw)
  if (!is.null(train_levels)) out <- factor(as.character(out), levels = train_levels)
  out
}

safe_cindex_lp <- function(time, status, lp) {
  ok <- is.finite(time) & !is.na(status) & is.finite(lp)
  time <- time[ok]
  status <- status[ok]
  lp <- lp[ok]
  if (length(time) < 2 || sum(status) == 0 || length(unique(lp)) < 2) return(NA_real_)
  out <- tryCatch(
    survival::concordance(Surv(time, status) ~ lp, reverse = TRUE)$concordance,
    error = function(e) NA_real_
  )
  scalar_num(out)
}

horizon_col <- function(prefix, h) {
  paste0(prefix, "_", gsub("\\.", "p", as.character(h)), "y")
}

make_baseline_survival <- function(fit, horizons) {
  bh <- tryCatch(basehaz(fit, centered = FALSE), error = function(e) NULL)
  if (is.null(bh) || nrow(bh) == 0) {
    return(list(table = tibble(time = numeric(), hazard = numeric(), survival = numeric()),
                s0 = setNames(rep(NA_real_, length(horizons)), as.character(horizons))))
  }
  bh <- as_tibble(bh) %>%
    transmute(time = as.numeric(time), hazard = as.numeric(hazard), survival = exp(-hazard)) %>%
    filter(is.finite(time), is.finite(hazard)) %>%
    arrange(time)
  s0 <- vapply(horizons, function(h) {
    if (nrow(bh) == 0) return(NA_real_)
    stats::approx(x = bh$time, y = bh$survival, xout = h, method = "constant",
                  f = 0, rule = 2, ties = "ordered")$y
  }, FUN.VALUE = numeric(1))
  names(s0) <- as.character(horizons)
  list(table = bh, s0 = s0)
}

event_status_at_horizon <- function(time, status, horizon) {
  out <- rep(NA_integer_, length(time))
  event_by_h <- !is.na(status) & status == 1L & is.finite(time) & time <= horizon
  event_free_at_h <- !is.na(status) & is.finite(time) & time >= horizon
  out[event_by_h] <- 1L
  out[event_free_at_h] <- 0L
  out
}

risk_group_from_percentile <- function(p) {
  dplyr::case_when(
    is.na(p) ~ NA_character_,
    p < 0.50 ~ "low",
    p < 0.80 ~ "moderate",
    p < 0.95 ~ "high",
    TRUE ~ "very_high"
  )
}

calibration_metrics <- function(event, risk, horizon, cindex = NA_real_) {
  n_total <- length(event)
  ok <- !is.na(event) & is.finite(risk)
  event_cc <- event[ok]
  risk_cc <- pmin(pmax(risk[ok], .Machine$double.eps), 1 - .Machine$double.eps)
  n_informative <- length(event_cc)
  events_horizon <- sum(event_cc == 1L, na.rm = TRUE)
  nonevents_observed <- sum(event_cc == 0L, na.rm = TRUE)
  prop_informative <- if (n_total > 0) n_informative / n_total else NA_real_

  eval_status <- dplyr::case_when(
    n_total == 0 ~ "SKIP_no_rows",
    n_informative == 0 ~ "SKIP_no_informative_followup",
    length(unique(event_cc)) < 2 ~ "SKIP_one_event_class",
    TRUE ~ "OK"
  )

  lp_risk <- qlogis(risk_cc)
  cal_fit <- if (eval_status == "OK" && stats::sd(lp_risk) > 0) {
    tryCatch(glm(event_cc ~ lp_risk, family = binomial()), error = function(e) NULL)
  } else NULL
  observed_cc <- if (n_informative > 0) mean(event_cc) else NA_real_
  predicted_cc <- if (n_informative > 0) mean(risk_cc) else NA_real_

  tibble(
    horizon = horizon,
    n_total = n_total,
    n_informative = n_informative,
    events_horizon = events_horizon,
    nonevents_observed_to_horizon = nonevents_observed,
    prop_informative = prop_informative,
    mean_predicted_risk = predicted_cc,
    observed_risk_complete_case = observed_cc,
    cindex = cindex,
    brier = if (n_informative > 0) mean((event_cc - risk_cc)^2) else NA_real_,
    calibration_slope = if (!is.null(cal_fit)) unname(coef(cal_fit)[["lp_risk"]]) else NA_real_,
    calibration_intercept = if (!is.null(cal_fit)) unname(coef(cal_fit)[["(Intercept)"]]) else NA_real_,
    eval_status = eval_status,
    horizon_years = horizon,
    n_calibration = n_informative,
    observed_risk = observed_cc,
    calibration_in_large = observed_cc - predicted_cc
  )
}

decile_metrics <- function(scores, horizon, risk_col, event_col) {
  risk <- scores[[risk_col]]
  event <- scores[[event_col]]
  ok <- is.finite(risk) & !is.na(event)
  if (sum(ok) == 0) return(tibble())
  dd <- scores[ok, , drop = FALSE]
  dd$risk_bin <- dplyr::ntile(dd[[risk_col]], 10)
  dd %>%
    group_by(risk_bin) %>%
    summarise(
      horizon = horizon,
      n_bin = n(),
      events_bin = sum(.data[[event_col]] == 1L, na.rm = TRUE),
      observed_risk = mean(.data[[event_col]], na.rm = TRUE),
      mean_predicted_risk = mean(.data[[risk_col]], na.rm = TRUE),
      min_predicted_risk = min(.data[[risk_col]], na.rm = TRUE),
      max_predicted_risk = max(.data[[risk_col]], na.rm = TRUE),
      bin_type = "decile",
      horizon_years = horizon,
      risk_decile = first(risk_bin),
      n = n_bin,
      .groups = "drop"
    )
}

extract_term <- function(fit, term = "pes_z") {
  sm <- summary(fit)
  if (!(term %in% rownames(sm$coefficients))) {
    return(list(HR = NA_real_, L95 = NA_real_, U95 = NA_real_, p = NA_real_))
  }
  beta <- sm$coefficients[term, "coef"]
  se <- sm$coefficients[term, "se(coef)"]
  list(
    HR = exp(beta),
    L95 = exp(beta - 1.96 * se),
    U95 = exp(beta + 1.96 * se),
    p = sm$coefficients[term, "Pr(>|z|)"]
  )
}

prepare_base <- function(score_df, covars_df, t2e_df, covars_used, pes_col,
                         landmark_label, exposure_levels = NULL) {
  covars_present <- intersect(covars_used, names(covars_df))
  score_df <- as.data.frame(score_df)
  score_df$eid <- as.integer(score_df$eid)
  score_df$pes_z <- suppressWarnings(as.numeric(score_df[[pes_col]]))
  exposure_type <- unique(score_df$exposure_type)[1]
  score_df$y_model <- format_exposure_for_cox(score_df$y_raw, exposure_type, exposure_levels)

  age_col <- "age_when_attended_assessment_centre_f21003_0_0"
  if (!age_col %in% names(covars_df)) stop("covars_df missing visit age column: ", age_col)

  base <- score_df %>%
    select(eid, instance, exposure_id, exposure_type, y_raw, y_model, pes_z) %>%
    inner_join(covars_df[, c("eid", age_col, covars_present), drop = FALSE], by = "eid") %>%
    inner_join(t2e_df, by = "eid") %>%
    mutate(
      landmark = landmark_label,
      landmark_age = suppressWarnings(as.numeric(.data[[age_col]]))
    )

  list(base = base, covars_present = covars_present)
}

fit_one_disease <- function(train_base, disease_age_col, exposure_id, exposure_type,
                            pes_used, pes_scale, covars_used, min_n, min_events,
                            risk_horizons) {
  train_surv <- survival_from_landmark(train_base, disease_age_col)
  cov_part <- if (length(covars_used) > 0) paste(bt(covars_used), collapse = " + ") else "1"
  formulas <- list(
    M0_covars = as.formula(paste0("Surv(survtime_standard, DZ_status) ~ ", cov_part)),
    M1_covars_PES = as.formula(paste0("Surv(survtime_standard, DZ_status) ~ pes_z + ", cov_part)),
    M2_covars_exposure = as.formula(paste0("Surv(survtime_standard, DZ_status) ~ y_model + ", cov_part)),
    M3_covars_exposure_PES = as.formula(paste0("Surv(survtime_standard, DZ_status) ~ y_model + pes_z + ", cov_part))
  )
  required <- list(
    M0_covars = c("survtime_standard", "DZ_status", covars_used),
    M1_covars_PES = c("survtime_standard", "DZ_status", "pes_z", covars_used),
    M2_covars_exposure = c("survtime_standard", "DZ_status", "y_model", covars_used),
    M3_covars_exposure_PES = c("survtime_standard", "DZ_status", "y_model", "pes_z", covars_used)
  )

  fits <- list()
  risk_specs <- list()
  fit_rows <- list()
  for (model_name in names(formulas)) {
    req <- intersect(required[[model_name]], names(train_surv))
    dt <- train_surv[complete.cases(train_surv[, req, drop = FALSE]), , drop = FALSE]
    n_i <- nrow(dt)
    ev_i <- sum(dt$DZ_status)
    status <- "OK"
    fit <- NULL
    if (n_i < min_n || ev_i < min_events) {
      status <- ifelse(n_i < min_n, "SKIP_low_n", "SKIP_low_events")
    } else if (model_name %in% c("M2_covars_exposure", "M3_covars_exposure_PES") &&
               exposure_type != "continuous" && length(unique(dt$y_model)) < 2) {
      status <- "SKIP_one_exposure_level"
    } else {
      fit <- tryCatch(coxph(formulas[[model_name]], data = dt, x = FALSE, y = FALSE, model = FALSE),
                      error = function(e) structure(list(error = conditionMessage(e)), class = "fit_error"))
      if (inherits(fit, "fit_error")) {
        status <- paste0("ERROR: ", fit$error)
        fit <- NULL
      }
    }
    fits[[model_name]] <- fit
    risk_specs[[model_name]] <- NULL
    if (!is.null(fit)) {
      bs <- make_baseline_survival(fit, risk_horizons)
      risk_specs[[model_name]] <- list(
        disease_age_col = disease_age_col,
        exposure_id = exposure_id,
        exposure_type = exposure_type,
        pes_used = pes_used,
        pes_scale = pes_scale,
        model = model_name,
        formula = deparse(formulas[[model_name]]),
        covars_used = covars_used,
        coefficients = coef(fit),
        xlevels = fit$xlevels,
        contrasts = fit$contrasts,
        baseline_survival = bs$table,
        s0_horizons = bs$s0,
        horizons = risk_horizons,
        train_n = n_i,
        train_events = ev_i,
        fit_status = status
      )
    }
    pes_term <- if (!is.null(fit) && model_name %in% c("M1_covars_PES", "M3_covars_exposure_PES")) extract_term(fit) else list(HR = NA_real_, L95 = NA_real_, U95 = NA_real_, p = NA_real_)
    fit_rows[[model_name]] <- tibble(
      disease_age_col = disease_age_col,
      exposure_id = exposure_id,
      exposure_type = exposure_type,
      pes_used = pes_used,
      pes_scale = pes_scale,
      model = model_name,
      train_n = n_i,
      train_events = ev_i,
      train_cindex_apparent = if (!is.null(fit)) safe_cindex_lp(dt$survtime_standard, dt$DZ_status, predict(fit, newdata = dt, type = "lp", reference = "zero")) else NA_real_,
      HR_PES_perSD = pes_term$HR,
      HR_PES_L95 = pes_term$L95,
      HR_PES_U95 = pes_term$U95,
      p_PES = pes_term$p,
      fit_status = status
    )
  }

  list(fits = fits, risk_specs = risk_specs, train_surv = train_surv, fit_summary = bind_rows(fit_rows))
}

score_one_landmark <- function(fit_obj, eval_base, disease_age_col, landmark_label,
                               eval_min_events, risk_horizons, save_person_scores = FALSE,
                               save_risk_scores = FALSE) {
  eval_surv <- survival_from_landmark(eval_base, disease_age_col)
  rows <- list()
  person_rows <- list()
  risk_rows <- list()
  risk_metric_rows <- list()
  decile_rows <- list()
  for (model_name in names(fit_obj$fits)) {
    fit <- fit_obj$fits[[model_name]]
    if (is.null(fit)) {
      rows[[model_name]] <- tibble(
        landmark = landmark_label,
        model = model_name,
        n = nrow(eval_surv),
        events = sum(eval_surv$DZ_status),
        cindex = NA_real_,
        eval_status = "SKIP_no_training_fit"
      )
      next
    }
    lp <- tryCatch(
      suppressWarnings(as.numeric(predict(fit, newdata = eval_surv, type = "lp", reference = "zero"))),
      error = function(e) structure(list(error = conditionMessage(e)), class = "pred_error")
    )
    if (inherits(lp, "pred_error")) {
      rows[[model_name]] <- tibble(
        landmark = landmark_label,
        model = model_name,
        n = nrow(eval_surv),
        events = sum(eval_surv$DZ_status),
        cindex = NA_real_,
        eval_status = paste0("ERROR: ", lp$error)
      )
      next
    }
    ok <- is.finite(lp)
    dt <- eval_surv[ok, , drop = FALSE]
    lp <- lp[ok]
    n_i <- nrow(dt)
    ev_i <- sum(dt$DZ_status)
    status <- ifelse(n_i == 0, "SKIP_no_rows",
                     ifelse(ev_i < eval_min_events, "SKIP_low_events", "OK"))
    cindex_i <- if (status == "OK") safe_cindex_lp(dt$survtime_standard, dt$DZ_status, lp) else NA_real_

    risk_tbl <- NULL
    spec <- fit_obj$risk_specs[[model_name]]
    if (!is.null(spec) && n_i > 0) {
      risk_tbl <- tibble(
        eid = dt$eid,
        instance = dt$instance,
        landmark = landmark_label,
        disease_age_col = disease_age_col,
        exposure_id = spec$exposure_id,
        exposure_type = spec$exposure_type,
        model = model_name,
        pes_used = spec$pes_used,
        pes_scale = spec$pes_scale,
        pes_z = dt$pes_z,
        y_raw = dt$y_raw,
        landmark_age = dt$landmark_age,
        survtime_standard = dt$survtime_standard,
        DZ_status = dt$DZ_status,
        event_age = dt$event_age,
        censor_age = dt$censor_age,
        LP = lp,
        lp = lp
      )
      for (h in risk_horizons) {
        risk_col <- horizon_col("risk", h)
        event_col <- horizon_col("event_status", h)
        pct_col <- horizon_col("risk_percentile", h)
        grp_col <- horizon_col("risk_group", h)
        s0 <- spec$s0_horizons[[as.character(h)]]
        risk_tbl[[risk_col]] <- if (is.finite(s0)) 1 - (s0 ^ exp(lp)) else NA_real_
        risk_tbl[[risk_col]] <- pmin(pmax(risk_tbl[[risk_col]], 0), 1)
        risk_tbl[[event_col]] <- event_status_at_horizon(dt$survtime_standard, dt$DZ_status, h)
        ok_risk <- is.finite(risk_tbl[[risk_col]])
        risk_tbl[[pct_col]] <- NA_real_
        if (sum(ok_risk) > 0) {
          risk_tbl[[pct_col]][ok_risk] <- percent_rank(risk_tbl[[risk_col]][ok_risk])
        }
        risk_tbl[[grp_col]] <- risk_group_from_percentile(risk_tbl[[pct_col]])
      }
      if (save_risk_scores) risk_rows[[model_name]] <- risk_tbl

      for (h in risk_horizons) {
        risk_col <- horizon_col("risk", h)
        event_col <- horizon_col("event_status", h)
        cal <- calibration_metrics(risk_tbl[[event_col]], risk_tbl[[risk_col]], h, cindex_i) %>%
          mutate(landmark = landmark_label, disease_age_col = disease_age_col,
                 exposure_id = spec$exposure_id, exposure_type = spec$exposure_type,
                 pes_used = spec$pes_used, pes_scale = spec$pes_scale,
                 model = model_name, .before = 1)
        risk_metric_rows[[paste(model_name, h, sep = "_")]] <- cal
        dec <- decile_metrics(risk_tbl, h, risk_col, event_col)
        if (nrow(dec) > 0) {
          decile_rows[[paste(model_name, h, sep = "_")]] <- dec %>%
            mutate(landmark = landmark_label, disease_age_col = disease_age_col,
                   exposure_id = spec$exposure_id, exposure_type = spec$exposure_type,
                   pes_used = spec$pes_used, pes_scale = spec$pes_scale,
                   model = model_name, .before = 1)
        }
      }
    }

    rows[[model_name]] <- tibble(
      landmark = landmark_label,
      model = model_name,
      n = n_i,
      events = ev_i,
      cindex = cindex_i,
      eval_status = status
    )
    if (save_person_scores) {
      person_rows[[model_name]] <- tibble(
        eid = dt$eid,
        instance = dt$instance,
        landmark = landmark_label,
        disease_age_col = disease_age_col,
        model = model_name,
        lp = lp,
        y_raw = dt$y_raw,
        DZ_status = dt$DZ_status,
        survtime_standard = dt$survtime_standard,
        landmark_age = dt$landmark_age,
        event_age = dt$event_age,
        censor_age = dt$censor_age
      )
    }
  }
  list(
    eval = bind_rows(rows),
    person = bind_rows(person_rows),
    risk_scores = bind_rows(risk_rows),
    risk_metrics = bind_rows(risk_metric_rows),
    risk_deciles = bind_rows(decile_rows)
  )
}

add_deltas <- function(eval_dt) {
  if (nrow(eval_dt) == 0) return(eval_dt)
  wide <- as.data.table(eval_dt)
  wide <- dcast(
    wide,
    disease_age_col + exposure_id + exposure_type + pes_used + pes_scale + landmark ~ model,
    value.var = c("cindex", "n", "events", "eval_status")
  )
  for (cc in c("cindex_M0_covars", "cindex_M1_covars_PES", "cindex_M2_covars_exposure", "cindex_M3_covars_exposure_PES")) {
    if (!cc %in% names(wide)) wide[, (cc) := NA_real_]
  }
  wide[, delta_cindex_M0_to_M1_PES := cindex_M1_covars_PES - cindex_M0_covars]
  wide[, delta_cindex_M0_to_M2_exposure := cindex_M2_covars_exposure - cindex_M0_covars]
  wide[, delta_cindex_M2_to_M3_add_PES := cindex_M3_covars_exposure_PES - cindex_M2_covars_exposure]
  as_tibble(wide)
}

summarise_person_deltas <- function(person_dt) {
  if (nrow(person_dt) == 0) return(tibble())
  dt <- as_tibble(person_dt) %>%
    filter(grepl("^holdout_i", landmark)) %>%
    select(eid, instance, disease_age_col, exposure_id, exposure_type, pes_used, model,
           lp, DZ_status, survtime_standard)
  pairs <- list(c(0L, 2L), c(0L, 3L), c(2L, 3L))
  out <- list()
  for (pp in pairs) {
    a <- pp[[1]]
    b <- pp[[2]]
    left <- dt %>% filter(instance == a) %>% rename(lp_from = lp)
    right <- dt %>% filter(instance == b) %>% rename(lp_to = lp, DZ_status_to = DZ_status, survtime_to = survtime_standard)
    joined <- inner_join(
      left %>% select(eid, disease_age_col, exposure_id, exposure_type, pes_used, model, lp_from),
      right %>% select(eid, disease_age_col, exposure_id, exposure_type, pes_used, model, lp_to, DZ_status_to, survtime_to),
      by = c("eid", "disease_age_col", "exposure_id", "exposure_type", "pes_used", "model")
    ) %>%
      mutate(pair = paste0("i", a, "_to_i", b), delta_lp = lp_to - lp_from)
    if (nrow(joined) == 0) next
    sm <- joined %>%
      group_by(disease_age_col, exposure_id, exposure_type, pes_used, model, pair) %>%
      summarise(
        n_pair = n(),
        events_after_later = sum(DZ_status_to, na.rm = TRUE),
        mean_delta_lp = mean(delta_lp, na.rm = TRUE),
        sd_delta_lp = sd(delta_lp, na.rm = TRUE),
        median_delta_lp = median(delta_lp, na.rm = TRUE),
        q10_delta_lp = quantile(delta_lp, 0.10, na.rm = TRUE),
        q90_delta_lp = quantile(delta_lp, 0.90, na.rm = TRUE),
        mean_delta_lp_event = mean(delta_lp[DZ_status_to == 1], na.rm = TRUE),
        mean_delta_lp_nonevent = mean(delta_lp[DZ_status_to == 0], na.rm = TRUE),
        cindex_delta_lp = safe_cindex_lp(survtime_to, DZ_status_to, delta_lp),
        .groups = "drop"
      )
    out[[length(out) + 1L]] <- sm
  }
  bind_rows(out)
}

main <- function() {
  args <- parse_args()
  covar_type <- args$covar_type
  exposure_id <- args$exposure_id
  if (is.na(args$result_dir)) args$result_dir <- file.path(cfg$out_dir, covar_type)
  if (!dir.exists(args$result_dir)) stop("Missing result_dir: ", args$result_dir)

  message_ts("Loading HEAP.rds + deriving longitudinal PXS: ", cfg$paths$heap_rds)
  heap <- readRDS(cfg$paths$heap_rds)
  message_ts("Loading disease time-to-event data from HEAP.rds $disease$DZ_df")
  t2e_df <- load_t2e(heap)
  pxs <- as_long_pxs(as_pxs_longitudinal(heap))
  rm(heap); invisible(gc())
  covars_used <- cfg$CovarSpec[[covar_type]]
  if (is.null(covars_used) || length(covars_used) == 0) stop("No covariates found for covarType: ", covar_type)
  covars_long <- as.data.frame(pxs$covars_df)
  if (!all(c("eid", "instance") %in% names(covars_long))) stop("covars_df must include eid and instance")

  # SAMPLE filter (healthy-at-baseline exclusion). Restricting covars_long to
  # prevalent_major_disease==0 propagates through cov_i0/cov_i -> prepare_base's
  # inner_join, so the frozen-risk Cox runs on the same healthy subset as the PES
  # training (Module6_prod_longitudinal.R). prevalent_major_disease is a sample
  # selector here, never a Cox covariate.
  if (covar_type %in% EXCLUDE_PREVALENT_COVAR_TYPES) {
    if (!"prevalent_major_disease" %in% names(covars_long))
      stop("exclude_prevalent requires prevalent_major_disease in covars_df ",
           "(rerun HEAP_loader with the prevalent-disease feature)")
    keep_eids <- unique(covars_long$eid[!is.na(covars_long$prevalent_major_disease) &
                                          covars_long$prevalent_major_disease == 0])
    n_before <- length(unique(covars_long$eid))
    covars_long <- covars_long[covars_long$eid %in% keep_eids, , drop = FALSE]
    message_ts("[exclude_prevalent] kept ", length(unique(covars_long$eid)), "/", n_before,
               " participants (prevalent_major_disease==0)")
  }

  disease_age_cols <- grep("^age_", names(t2e_df), value = TRUE)
  # Drop the censoring/administrative age_* columns — they start with "age_" but are
  # NOT disease outcomes; treating age_of_lastfollowup (everyone has it) as a disease
  # produces a degenerate, perfectly-separated Cox fit.
  disease_age_cols <- setdiff(disease_age_cols, c("age_of_removal_0_0", "age_of_lastfollowup"))
  if (!is.na(args$disease_regex)) disease_age_cols <- disease_age_cols[grepl(args$disease_regex, disease_age_cols)]
  if (is.finite(args$max_diseases)) disease_age_cols <- head(disease_age_cols, args$max_diseases)
  if (length(disease_age_cols) == 0) stop("No disease columns selected.")
  message_ts("Selected disease columns: ", length(disease_age_cols))

  train_oof <- read_score_table(args$result_dir, covar_type, exposure_id, "TrainOOF")
  train_oof$instance <- 0L
  holdout <- read_score_table(args$result_dir, covar_type, exposure_id, "HoldoutScores")
  holdout$instance <- as.integer(holdout$instance)
  exposure_type <- unique(train_oof$exposure_type)[1]
  exposure_levels <- if (exposure_type == "continuous") NULL else levels(format_exposure_for_cox(train_oof$y_raw, exposure_type))

  out_prefix <- file.path(args$result_dir, paste0("PESlong_", covar_type, "_", safe_name(exposure_id)))
  frozen_stem <- paste0("_FrozenCox", args$output_tag)
  fit_tsv <- paste0(out_prefix, frozen_stem, "TrainFits.tsv")
  eval_long_tsv <- paste0(out_prefix, frozen_stem, "EvalLong.tsv")
  eval_tsv <- paste0(out_prefix, frozen_stem, "Eval.tsv")
  delta_tsv <- paste0(out_prefix, frozen_stem, "PersonDelta.tsv")
  model_rds <- paste0(out_prefix, frozen_stem, "Models.rds")
  person_tsv <- paste0(out_prefix, frozen_stem, "PersonScores.tsv")
  risk_tag <- if (nzchar(args$output_tag)) paste0("_", args$output_tag) else ""
  risk_model_rds <- paste0(out_prefix, "_RiskCalculatorModels", risk_tag, ".rds")
  risk_scores_tsv <- paste0(out_prefix, "_RiskCalculatorScores", risk_tag, ".tsv")
  risk_metrics_tsv <- paste0(out_prefix, "_RiskCalculatorMetrics", risk_tag, ".tsv")
  risk_deciles_tsv <- paste0(out_prefix, "_RiskCalculatorDeciles", risk_tag, ".tsv")

  out_files <- c(fit_tsv, eval_long_tsv, eval_tsv, delta_tsv, person_tsv, model_rds,
                 risk_model_rds, risk_scores_tsv, risk_metrics_tsv, risk_deciles_tsv)
  invisible(file.remove(out_files[file.exists(out_files)]))

  saved_models <- list()
  risk_model_specs <- list()

  pb <- progress_bar$new(
    format = "  Frozen Cox :current/:total [:bar] :percent eta: :eta",
    total = length(disease_age_cols) * length(args$score_types),
    clear = FALSE
  )

  for (score_type in args$score_types) {
    train_pes_col <- if (score_type == "prot") "pes_prot_z" else "pes_full_z"
    holdout_pes_col <- if (args$pes_scale == "oof") {
      if (score_type == "prot") "pes_prot_z_oofscale" else "pes_full_z_oofscale"
    } else {
      if (score_type == "prot") "pes_prot_z_finalscale" else "pes_full_z_finalscale"
    }
    fallback_holdout_pes_col <- if (score_type == "prot") "pes_prot_z" else "pes_full_z"
    if (!holdout_pes_col %in% names(holdout)) {
      message_ts("Holdout column ", holdout_pes_col, " absent; falling back to ", fallback_holdout_pes_col)
      holdout_pes_col <- fallback_holdout_pes_col
    }
    cov_i0 <- covars_long[covars_long$instance == 0L, , drop = FALSE]
    train_prep <- prepare_base(train_oof, cov_i0, t2e_df, covars_used, train_pes_col, "train_i0", exposure_levels)

    hold_preps <- list()
    for (inst in args$instances) {
      hold_i <- holdout[holdout$instance == inst, , drop = FALSE]
      if (nrow(hold_i) == 0) next
      cov_i <- covars_long[covars_long$instance == inst, , drop = FALSE]
      hold_preps[[paste0("holdout_i", inst)]] <- prepare_base(
        hold_i, cov_i, t2e_df, covars_used, holdout_pes_col, paste0("holdout_i", inst), exposure_levels
      )$base
    }

    for (disease_age_col in disease_age_cols) {
      pb$tick()
      fit_obj <- fit_one_disease(
        train_base = train_prep$base,
        disease_age_col = disease_age_col,
        exposure_id = exposure_id,
        exposure_type = exposure_type,
        pes_used = train_pes_col,
        pes_scale = args$pes_scale,
        covars_used = train_prep$covars_present,
        min_n = args$min_n,
        min_events = args$min_events,
        risk_horizons = args$risk_horizons
      )
      append_table(fit_obj$fit_summary, fit_tsv)
      if (args$save_models) saved_models[[paste(disease_age_col, train_pes_col, sep = "__")]] <- fit_obj$fits
      risk_model_specs[[paste(disease_age_col, train_pes_col, sep = "__")]] <- fit_obj$risk_specs

      train_score <- score_one_landmark(
        fit_obj = fit_obj,
        eval_base = train_prep$base,
        disease_age_col = disease_age_col,
        landmark_label = "train_i0_apparent",
        eval_min_events = args$eval_min_events,
        risk_horizons = args$risk_horizons,
        save_person_scores = FALSE
      )$eval %>%
        mutate(disease_age_col = disease_age_col, exposure_id = exposure_id,
               exposure_type = exposure_type, pes_used = train_pes_col, pes_scale = args$pes_scale)
      append_table(train_score, eval_long_tsv)

      disease_person <- list()
      for (landmark_label in names(hold_preps)) {
        scored <- score_one_landmark(
          fit_obj = fit_obj,
          eval_base = hold_preps[[landmark_label]],
          disease_age_col = disease_age_col,
          landmark_label = landmark_label,
          eval_min_events = args$eval_min_events,
          risk_horizons = args$risk_horizons,
          save_person_scores = args$save_person_scores,
          save_risk_scores = args$save_risk_scores
        )
        eval_block <- scored$eval %>%
          mutate(disease_age_col = disease_age_col, exposure_id = exposure_id,
                 exposure_type = exposure_type, pes_used = train_pes_col, pes_scale = args$pes_scale)
        append_table(eval_block, eval_long_tsv)
        append_table(scored$risk_metrics, risk_metrics_tsv)
        append_table(scored$risk_deciles, risk_deciles_tsv)
        if (args$save_risk_scores && nrow(scored$risk_scores) > 0) {
          append_table(scored$risk_scores, risk_scores_tsv)
        }
        if (args$save_person_scores && nrow(scored$person) > 0) {
          person_block <- scored$person %>%
            mutate(exposure_id = exposure_id, exposure_type = exposure_type, pes_used = train_pes_col,
                   pes_scale = args$pes_scale)
          append_table(person_block, person_tsv)
          disease_person[[length(disease_person) + 1L]] <- person_block
        }
      }

      if (args$save_person_scores && length(disease_person) > 0) {
        person_delta <- summarise_person_deltas(bind_rows(disease_person))
        append_table(person_delta, delta_tsv)
      }
    }
  }

  eval_long <- fread(eval_long_tsv)
  eval_wide <- add_deltas(eval_long)
  fwrite(as.data.table(eval_wide), eval_tsv, sep = "\t")

  fit_summary_for_risk <- fread(fit_tsv)
  saveRDS(
    list(
      exposure_id = exposure_id,
      covar_type = covar_type,
      score_types = args$score_types,
      pes_scale = args$pes_scale,
      risk_horizons = args$risk_horizons,
      instances = args$instances,
      covars_used = covars_used,
      disease_age_cols = disease_age_cols,
      model_specs = risk_model_specs,
      fit_summary = fit_summary_for_risk,
      note = "Deployable Cox risk artifact: use model_specs baseline_survival and s0_horizons with Risk = 1 - S0(t)^exp(LP)."
    ),
    risk_model_rds
  )

  if (args$save_models) {
    saveRDS(
      list(
        exposure_id = exposure_id,
        covar_type = covar_type,
        score_types = args$score_types,
        instances = args$instances,
        covars_used = covars_used,
        models = saved_models,
        fit_summary = fit_summary_for_risk
      ),
      model_rds
    )
  }

  message_ts("Saved frozen Cox train fit summary: ", fit_tsv)
  message_ts("Saved frozen Cox eval long table: ", eval_long_tsv)
  message_ts("Saved frozen Cox eval summary: ", eval_tsv)
  if (args$save_person_scores) {
    message_ts("Saved frozen Cox person scores: ", person_tsv)
    message_ts("Saved frozen Cox person delta summary: ", delta_tsv)
  }
  message_ts("Saved risk calculator model artifact: ", risk_model_rds)
  if (file.exists(risk_metrics_tsv)) message_ts("Saved risk calculator metrics: ", risk_metrics_tsv)
  if (file.exists(risk_deciles_tsv)) message_ts("Saved risk calculator deciles: ", risk_deciles_tsv)
  if (args$save_risk_scores) message_ts("Saved risk calculator scores: ", risk_scores_tsv)
  if (args$save_models) message_ts("Saved frozen Cox model artifact: ", model_rds)
  message_ts("DONE")
}

if (identical(environment(), globalenv())) {
  main()
}
