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

# Longitudinal out-of-sample validation for Module 6 protein-exposure scores.
# Trains the PES models and validates exposure prediction on repeat-visit proteomics;
# the disease/Cox stage is handled separately (Module6_longitudinal_cox_frozenrisk.R).

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(purrr)
  library(tidyr)
  library(glmnet)
  library(caret)
  library(Matrix)
  library(pROC)
  library(progress)
  library(tibble)
})

cfg <- list(
  seed = 123,
  kfold = 5,
  miss_rate_prot = 0.20,
  out_dir = heap_project_output("module6_pes_longitudinal"),
  paths = list(
    # Module 6 reads HEAP.rds directly and derives the longitudinal PXS in-process
    # via as_pxs_longitudinal() (00_paths.R) -- no separate longitudinal loader RDS.
    heap_rds = heap_loader_rds,
    exposure_manifest = heap_project_output("module6_pes_test", "exposure_specs.tsv")
  ),
  # covariate sets mirror config/covariates/covariate_sets.yml (v2, base-primary)
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
      "body_mass_index_bmi_f23104_0_0", "fasting_time_f74_0_0",
      "assessment_season",
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
    # exclude_prevalent SENSITIVITY (reviewer exclusion): `base` covariate
    # adjustment on the healthy-at-baseline subset. Covariates are IDENTICAL to
    # `base`; the difference is the SAMPLE -- assemble_for_exposure_long() drops
    # participants with prevalent_major_disease==1 (sample_filter below). Outputs
    # land under module6_pes_longitudinal/base_exclprev/ (no collision with base).
    base_exclprev = c(
      "age_when_attended_assessment_centre_f21003_0_0", "sex_f31_0_0",
      "age2", "age_sex", "age2_sex",
      "uk_biobank_assessment_centre_f54_0_0",
      paste0("genetic_principal_components_f22009_0_", 1:20)
    )
  )
)

# covar_types that run on the healthy-at-baseline subset (exclude prevalent
# major disease). Covariates come from CovarSpec; the sample filter is applied
# in assemble_for_exposure_long().
EXCLUDE_PREVALENT_COVAR_TYPES <- c("base_exclprev")

set.seed(cfg$seed)
dir.create(cfg$out_dir, showWarnings = FALSE, recursive = TRUE)

`%||%` <- function(x, y) if (!is.null(x)) x else y

message_ts <- function(...) {
  message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), ...)
}

safe_name <- function(x) {
  gsub("[^A-Za-z0-9_.-]+", "_", x)
}

canonicalize_exposure_id <- function(x) {
  gsub("(_f[0-9]+)_[23]_", "\\1_0_", x)
}

as_long_pxs <- function(x) {
  if (is.list(x) && !isS4(x)) return(x)
  if (!isS4(x)) stop("PXS object must be S4 or list")
  list(
    Elist = x@Elist,
    Elist_names = x@Elist_names,
    Eid_cat = x@Eid_cat,
    ordinalIDs = x@ordinalIDs,
    UKBprot_df = x@UKBprot_df,
    protIDs = x@protIDs,
    covars_df = x@covars_df,
    covars_list = x@covars_list
  )
}

missing_cols <- function(df, miss_rate = 0.2) {
  na_rate <- colMeans(is.na(df))
  names(na_rate[na_rate > miss_rate])
}

zscore_with <- function(x, center, scale) {
  if (!is.finite(scale) || scale == 0) return(rep(0, length(x)))
  (x - center) / scale
}

fit_zscore_params <- function(x) {
  list(center = mean(x, na.rm = TRUE), scale = stats::sd(x, na.rm = TRUE))
}

clip01 <- function(p, eps = 1e-15) pmin(pmax(p, eps), 1 - eps)

fit_outcome_encoder <- function(y_raw, exposure_type) {
  if (exposure_type == "continuous") {
    return(list(type = exposure_type, family = "gaussian", kind = "continuous"))
  }

  if (exposure_type == "binary") {
    if (is.numeric(y_raw) || is.integer(y_raw) || is.logical(y_raw)) {
      vals <- sort(unique(stats::na.omit(as.numeric(y_raw))))
      if (length(vals) != 2) stop("Binary exposure does not have exactly two observed values in training.")
      return(list(type = exposure_type, family = "binomial", kind = "binary_numeric", values = vals))
    }
    levs <- sort(unique(stats::na.omit(as.character(y_raw))))
    if (length(levs) != 2) stop("Binary exposure does not have exactly two observed levels in training.")
    return(list(type = exposure_type, family = "binomial", kind = "binary_factor", levels = levs))
  }

  if (exposure_type %in% c("ordinal", "categorical")) {
    levs <- sort(unique(stats::na.omit(as.character(y_raw))))
    if (length(levs) < 2) stop("Categorical/ordinal exposure has fewer than two observed levels in training.")
    if (length(levs) == 2) {
      return(list(type = exposure_type, family = "binomial", kind = paste0("binary_from_", exposure_type), levels = levs))
    }
    codes <- seq_along(levs) - 1L
    names(codes) <- levs
    return(list(type = exposure_type, family = "multinomial", kind = exposure_type, levels = levs, codes = codes))
  }

  stop("Unknown exposure_type: ", exposure_type)
}

apply_outcome_encoder <- function(y_raw, enc) {
  if (enc$family == "gaussian") {
    if (is.factor(y_raw)) {
      lev_num <- suppressWarnings(as.numeric(levels(y_raw)))
      if (all(is.finite(lev_num))) {
        return(lev_num[as.integer(y_raw)])
      }
      return(as.numeric(y_raw) - 1)
    }
    return(suppressWarnings(as.numeric(y_raw)))
  }

  if (enc$family == "binomial") {
    if (enc$kind == "binary_numeric") {
      y <- suppressWarnings(as.numeric(y_raw))
      out <- rep(NA_integer_, length(y))
      out[y == enc$values[1]] <- 0L
      out[y == enc$values[2]] <- 1L
      return(out)
    }
    yy <- as.character(y_raw)
    out <- rep(NA_integer_, length(yy))
    out[yy == enc$levels[1]] <- 0L
    out[yy == enc$levels[2]] <- 1L
    return(out)
  }

  yy <- as.character(y_raw)
  factor(ifelse(yy %in% enc$levels, yy, NA_character_), levels = enc$levels)
}

outcome_code <- function(y_raw, enc) {
  if (enc$family %in% c("gaussian", "binomial")) return(as.numeric(apply_outcome_encoder(y_raw, enc)))
  yy <- as.character(y_raw)
  as.numeric(enc$codes[yy])
}

fit_median_imputer <- function(X) {
  med <- apply(X, 2, function(v) stats::median(v, na.rm = TRUE))
  med[!is.finite(med)] <- 0
  list(median = med)
}

apply_median_imputer <- function(X, imp) {
  X2 <- as.matrix(X)
  idx <- which(is.na(X2), arr.ind = TRUE)
  if (nrow(idx) > 0) X2[idx] <- imp$median[idx[, 2]]
  X2
}

fit_covariate_recipe <- function(df, covars_used) {
  recipes <- lapply(covars_used, function(cc) {
    x <- df[[cc]]
    if (is.numeric(x) || is.integer(x) || is.logical(x)) {
      med <- stats::median(as.numeric(x), na.rm = TRUE)
      if (!is.finite(med)) med <- 0
      list(name = cc, type = "numeric", median = med)
    } else {
      levs <- sort(unique(stats::na.omit(as.character(x))))
      list(name = cc, type = "factor", levels = levs, missing_level = ".__MISSING__")
    }
  })
  names(recipes) <- covars_used
  recipes
}

apply_covariate_recipe <- function(df, recipe) {
  out <- data.frame(row.names = seq_len(nrow(df)))
  for (cc in names(recipe)) {
    rr <- recipe[[cc]]
    if (rr$type == "numeric") {
      x <- suppressWarnings(as.numeric(df[[cc]]))
      x[!is.finite(x)] <- rr$median
      out[[cc]] <- x
    } else {
      x <- as.character(df[[cc]])
      x[is.na(x) | !(x %in% rr$levels)] <- rr$missing_level
      out[[cc]] <- factor(x, levels = c(rr$levels, rr$missing_level))
    }
  }
  out
}

fit_covariate_matrix <- function(df, covars_used) {
  if (length(covars_used) == 0) {
    return(list(recipe = list(), terms = NULL, columns = character()))
  }
  recipe <- fit_covariate_recipe(df, covars_used)
  cov_df <- apply_covariate_recipe(df, recipe)
  terms_obj <- terms(reformulate(covars_used), data = cov_df)
  X <- model.matrix(terms_obj, data = cov_df, na.action = stats::na.pass)
  if (ncol(X) > 0 && colnames(X)[1] == "(Intercept)") X <- X[, -1, drop = FALSE]
  list(recipe = recipe, terms = terms_obj, columns = colnames(X))
}

apply_covariate_matrix <- function(df, covar_plan) {
  if (length(covar_plan$columns) == 0) return(matrix(nrow = nrow(df), ncol = 0))
  cov_df <- apply_covariate_recipe(df, covar_plan$recipe)
  X <- model.matrix(covar_plan$terms, data = cov_df, na.action = stats::na.pass)
  if (ncol(X) > 0 && colnames(X)[1] == "(Intercept)") X <- X[, -1, drop = FALSE]
  missing <- setdiff(covar_plan$columns, colnames(X))
  if (length(missing) > 0) {
    X <- cbind(X, matrix(0, nrow = nrow(X), ncol = length(missing), dimnames = list(NULL, missing)))
  }
  X <- X[, covar_plan$columns, drop = FALSE]
  X[is.na(X)] <- 0
  X
}

predict_scalar <- function(fit, newx, family, enc, lambda = "lambda.min") {
  if (is.null(fit)) stop("predict_scalar received NULL fit")
  if (family == "multinomial") {
    pred_obj <- predict(fit, newx = newx, s = lambda, type = "response")
    p <- if (is.array(pred_obj)) pred_obj[, , 1, drop = TRUE] else as.matrix(pred_obj)
    if (is.null(colnames(p))) colnames(p) <- enc$levels
    p <- p[, enc$levels, drop = FALSE]
    score <- as.numeric(p %*% matrix(enc$codes[colnames(p)], ncol = 1))
    cls <- colnames(p)[max.col(p, ties.method = "first")]
    return(list(score = score, class = cls, prob = p))
  }
  list(score = as.numeric(predict(fit, newx = newx, s = lambda, type = "response")), class = NULL, prob = NULL)
}

null_prediction <- function(y_train, n, enc) {
  if (enc$family == "gaussian") return(rep(mean(as.numeric(y_train), na.rm = TRUE), n))
  if (enc$family == "binomial") return(rep(mean(as.numeric(y_train), na.rm = TRUE), n))
  tab <- prop.table(table(y_train))
  levs <- names(tab)
  rep(sum(as.numeric(tab) * enc$codes[levs]), n)
}

metric_row <- function(y_raw, enc, pred) {
  if (enc$family == "gaussian") {
    y <- apply_outcome_encoder(y_raw, enc)
    ok <- is.finite(y) & is.finite(pred)
    y <- y[ok]
    p <- pred[ok]
    ss_res <- sum((y - p)^2)
    ss_tot <- sum((y - mean(y))^2)
    return(tibble(
      n = length(y),
      r2 = ifelse(ss_tot == 0, NA_real_, 1 - ss_res / ss_tot),
      correlation = ifelse(length(y) > 2 && stats::sd(y) > 0 && stats::sd(p) > 0, stats::cor(y, p), NA_real_),
      rmse = sqrt(mean((y - p)^2))
    ))
  }

  if (enc$family == "binomial") {
    y <- apply_outcome_encoder(y_raw, enc)
    ok <- !is.na(y) & is.finite(pred)
    y <- y[ok]
    p <- clip01(pred[ok])
    auc <- if (length(unique(y)) == 2) as.numeric(pROC::auc(pROC::roc(y, p, quiet = TRUE))) else NA_real_
    return(tibble(
      n = length(y),
      auc = auc,
      logloss = -mean(y * log(p) + (1 - y) * log(1 - p))
    ))
  }

  y_code <- outcome_code(y_raw, enc)
  ok <- is.finite(y_code) & is.finite(pred)
  y_code <- y_code[ok]
  p <- pred[ok]
  tibble(
    n = length(y_code),
    correlation_code = ifelse(length(y_code) > 2 && stats::sd(y_code) > 0 && stats::sd(p) > 0, stats::cor(y_code, p), NA_real_),
    rmse_code = sqrt(mean((y_code - p)^2)),
    n_levels = length(enc$levels),
    metric_note = "Multiclass/ordinal proxy uses expected class-code score."
  )
}

extract_nonzero_predictors <- function(fit, family, lambda = "lambda.min") {
  if (is.null(fit)) return(character())
  if (family != "multinomial") {
    b <- as.matrix(coef(fit, s = lambda))
    return(setdiff(rownames(b)[as.numeric(b[, 1]) != 0], "(Intercept)"))
  }
  coefs <- coef(fit, s = lambda)
  setdiff(unique(unlist(lapply(coefs, function(cm) rownames(cm)[as.numeric(cm) != 0]))), "(Intercept)")
}

fit_cv_glmnet <- function(X, y, family) {
  cv.glmnet(x = X, y = y, family = family, alpha = 1, standardize = TRUE)
}

fit_fold_artifact <- function(train_df, exposure_id, enc, prot_cols, covars_used) {
  y <- apply_outcome_encoder(train_df[[exposure_id]], enc)
  ok <- !is.na(y)
  train_df <- train_df[ok, , drop = FALSE]
  y <- y[ok]

  Xprot_raw <- as.matrix(train_df[, prot_cols, drop = FALSE])
  prot_imp <- fit_median_imputer(Xprot_raw)
  Xprot <- apply_median_imputer(Xprot_raw, prot_imp)

  cov_plan <- fit_covariate_matrix(train_df, covars_used)
  Xcov <- apply_covariate_matrix(train_df, cov_plan)
  Xfull <- cbind(Xprot, Xcov)

  fit_prot <- fit_cv_glmnet(Xprot, y, enc$family)
  fit_cov <- if (ncol(Xcov) > 0) fit_cv_glmnet(Xcov, y, enc$family) else NULL
  fit_full <- fit_cv_glmnet(Xfull, y, enc$family)

  list(
    exposure_id = exposure_id,
    encoder = enc,
    prot_cols = prot_cols,
    covars_used = covars_used,
    protein_imputer = prot_imp,
    covariate_plan = cov_plan,
    fit_prot = fit_prot,
    fit_cov = fit_cov,
    fit_full = fit_full,
    null_y = y
  )
}

score_artifact <- function(df, artifact) {
  Xprot_raw <- as.matrix(df[, artifact$prot_cols, drop = FALSE])
  Xprot <- apply_median_imputer(Xprot_raw, artifact$protein_imputer)
  Xcov <- apply_covariate_matrix(df, artifact$covariate_plan)
  Xfull <- cbind(Xprot, Xcov)

  pred_prot <- predict_scalar(artifact$fit_prot, Xprot, artifact$encoder$family, artifact$encoder)$score
  pred_cov <- if (!is.null(artifact$fit_cov)) {
    predict_scalar(artifact$fit_cov, Xcov, artifact$encoder$family, artifact$encoder)$score
  } else {
    null_prediction(artifact$null_y, nrow(df), artifact$encoder)
  }
  pred_full <- predict_scalar(artifact$fit_full, Xfull, artifact$encoder$family, artifact$encoder)$score

  tibble(
    eid = df$eid,
    instance = df$instance,
    y_raw = df[[artifact$exposure_id]],
    pred_prot = pred_prot,
    pred_cov = pred_cov,
    pred_full = pred_full
  )
}

oof_predictors_glmnet_long <- function(df, exposure_id, exposure_type, prot_cols, covars_used,
                                       kfold = 5, seed = 1) {
  set.seed(seed)
  enc <- fit_outcome_encoder(df[[exposure_id]], exposure_type)
  y_all <- apply_outcome_encoder(df[[exposure_id]], enc)
  ok <- !is.na(y_all)
  df <- df[ok, , drop = FALSE]
  y_all <- y_all[ok]
  exposure_values <- df[[exposure_id]]

  kfold <- min(kfold, nrow(df))
  folds <- if (enc$family == "gaussian") {
    caret::createFolds(y_all, k = kfold, list = TRUE)
  } else {
    caret::createFolds(as.factor(y_all), k = kfold, list = TRUE)
  }

  oof_prot <- rep(NA_real_, nrow(df))
  oof_cov <- rep(NA_real_, nrow(df))
  oof_full <- rep(NA_real_, nrow(df))
  oof_fold <- rep(NA_integer_, nrow(df))
  fold_metrics <- vector("list", length(folds))

  pb <- progress_bar$new(
    format = "  Fold :current/:total [:bar] :percent | elapsed: :elapsed | eta: :eta",
    total = length(folds), clear = FALSE, width = 60
  )

  for (ff in seq_along(folds)) {
    pb$tick()
    idx_te <- folds[[ff]]
    idx_tr <- setdiff(seq_len(nrow(df)), idx_te)
    fold_art <- fit_fold_artifact(df[idx_tr, , drop = FALSE], exposure_id, enc, prot_cols, covars_used)
    fold_score <- score_artifact(df[idx_te, , drop = FALSE], fold_art)
    oof_prot[idx_te] <- fold_score$pred_prot
    oof_cov[idx_te] <- fold_score$pred_cov
    oof_full[idx_te] <- fold_score$pred_full
    oof_fold[idx_te] <- ff

    met_prot <- metric_row(exposure_values[idx_te], enc, oof_prot[idx_te])
    met_cov <- metric_row(exposure_values[idx_te], enc, oof_cov[idx_te])
    met_full <- metric_row(exposure_values[idx_te], enc, oof_full[idx_te])
    fold_metrics[[ff]] <- bind_cols(
      tibble(fold = ff, family = enc$family, n_train = length(idx_tr), n_test = length(idx_te)),
      met_prot %>% rename_with(~ paste0(.x, "_prot")),
      met_cov %>% rename_with(~ paste0(.x, "_cov")),
      met_full %>% rename_with(~ paste0(.x, "_full"))
    )
  }

  z_prot <- fit_zscore_params(oof_prot)
  z_full <- fit_zscore_params(oof_full)

  oof_tbl <- tibble(
    eid = df$eid,
    instance = df$instance,
    exposure_id = rep(exposure_id, nrow(df)),
    exposure_type = rep(exposure_type, nrow(df)),
    y_raw = exposure_values,
    pred_prot = oof_prot,
    pred_cov = oof_cov,
    pred_full = oof_full,
    pes_prot_z = zscore_with(oof_prot, z_prot$center, z_prot$scale),
    pes_full_z = zscore_with(oof_full, z_full$center, z_full$scale),
    fold = oof_fold
  )

  overall_tbl <- bind_rows(
    metric_row(exposure_values, enc, oof_prot) %>% mutate(model = "prot_only"),
    metric_row(exposure_values, enc, oof_cov) %>% mutate(model = "cov_only"),
    metric_row(exposure_values, enc, oof_full) %>% mutate(model = "prot_plus_cov")
  ) %>%
    mutate(exposure_id = exposure_id, exposure_type = exposure_type)

  list(
    oof = oof_tbl,
    fold_metrics = bind_rows(fold_metrics) %>% mutate(exposure_id = exposure_id, exposure_type = exposure_type),
    overall_metrics = overall_tbl,
    encoder = enc,
    oof_z_params = list(prot = z_prot, full = z_full)
  )
}

fit_final_artifact <- function(df, exposure_id, exposure_type, prot_cols, covars_used, oof_fit) {
  enc <- fit_outcome_encoder(df[[exposure_id]], exposure_type)
  y <- apply_outcome_encoder(df[[exposure_id]], enc)
  ok <- !is.na(y)
  df <- df[ok, , drop = FALSE]

  art <- fit_fold_artifact(df, exposure_id, enc, prot_cols, covars_used)
  train_scores <- score_artifact(df, art)
  z_prot <- fit_zscore_params(train_scores$pred_prot)
  z_full <- fit_zscore_params(train_scores$pred_full)

  art$training_n <- nrow(df)
  art$exposure_type <- exposure_type
  art$lambda <- list(
    prot = art$fit_prot$lambda.min,
    cov = if (!is.null(art$fit_cov)) art$fit_cov$lambda.min else NA_real_,
    full = art$fit_full$lambda.min
  )
  art$selected_predictors <- list(
    prot = extract_nonzero_predictors(art$fit_prot, enc$family),
    cov = extract_nonzero_predictors(art$fit_cov, enc$family),
    full = extract_nonzero_predictors(art$fit_full, enc$family)
  )
  art$selected_proteins <- list(
    prot = intersect(art$selected_predictors$prot, prot_cols),
    full = intersect(art$selected_predictors$full, prot_cols)
  )
  art$z_params <- list(prot = z_prot, full = z_full)
  art$oof_z_params <- oof_fit$oof_z_params
  art$training_score_summary <- train_scores %>%
    summarise(
      n = n(),
      mean_pred_prot = mean(pred_prot, na.rm = TRUE),
      sd_pred_prot = sd(pred_prot, na.rm = TRUE),
      mean_pred_full = mean(pred_full, na.rm = TRUE),
      sd_pred_full = sd(pred_full, na.rm = TRUE)
    )
  art
}

assemble_for_exposure_long <- function(pxs, exposure_id, covars_used, miss_rate_prot = 0.2,
                                       sample_filter = "none") {
  pxs <- as_long_pxs(pxs)
  join_key <- c("eid", "instance")

  E_df <- pxs$Elist %>% purrr::reduce(full_join, by = join_key)
  if (!exposure_id %in% names(E_df)) stop("Exposure not found in longitudinal Elist: ", exposure_id)

  prot_df <- as.data.frame(pxs$UKBprot_df)
  cov_df <- as.data.frame(pxs$covars_df)
  split_df <- as.data.frame(pxs$split_df %||% data.frame(eid = unique(prot_df$eid)))

  if (!all(join_key %in% names(prot_df))) stop("UKBprot_df must include eid and instance")
  if (!all(join_key %in% names(cov_df))) stop("covars_df must include eid and instance")

  covars_present <- intersect(covars_used, names(cov_df))
  missing_covars <- setdiff(covars_used, names(cov_df))

  df <- prot_df %>%
    inner_join(E_df[, c(join_key, exposure_id), drop = FALSE], by = join_key) %>%
    left_join(cov_df[, c(join_key, covars_present), drop = FALSE], by = join_key) %>%
    left_join(split_df, by = "eid")

  # SAMPLE filter (healthy-at-baseline exclusion). Drops participants with
  # prevalent_major_disease==1 from BOTH train and hold-out. The flag is NOT a
  # model covariate -- it is used only to subset the sample. prevalent_major_disease
  # is broadcast across visits by the loader, so a single keep-set by eid applies.
  if (identical(sample_filter, "exclude_prevalent")) {
    if (!"prevalent_major_disease" %in% names(cov_df))
      stop("sample_filter='exclude_prevalent' requires prevalent_major_disease in covars_df ",
           "(rerun HEAP_loader with the prevalent-disease feature)")
    keep_eids <- unique(cov_df$eid[!is.na(cov_df$prevalent_major_disease) &
                                     cov_df$prevalent_major_disease == 0])
    n_before <- length(unique(df$eid))
    df <- df[df$eid %in% keep_eids, , drop = FALSE]
    message("[exclude_prevalent] kept ", length(unique(df$eid)), "/", n_before,
            " participants (prevalent_major_disease==0)")
  } else if (!identical(sample_filter, "none")) {
    stop("Unknown sample_filter: ", sample_filter)
  }

  prot_cols <- setdiff(names(prot_df), join_key)

  train_idx <- df$analysis_set == "train_baseline_only" & df$instance == 0L
  train_df0 <- df[train_idx & !is.na(df[[exposure_id]]), , drop = FALSE]
  if (nrow(train_df0) < 50) stop("Too few training rows after baseline-only split: ", nrow(train_df0))

  drop_prot <- missing_cols(train_df0[, prot_cols, drop = FALSE], miss_rate = miss_rate_prot)
  keep_prot <- setdiff(prot_cols, drop_prot)

  train_df <- train_df0[, c(join_key, exposure_id, covars_present, keep_prot), drop = FALSE]
  holdout_df <- df[df$analysis_set == "holdout_repeat_proteomics" & df$instance %in% c(0L, 2L, 3L), , drop = FALSE]
  holdout_df <- holdout_df[, c(join_key, exposure_id, covars_present, keep_prot), drop = FALSE]

  list(
    train_df = train_df,
    holdout_df = holdout_df,
    prot_cols = keep_prot,
    covars_used = covars_present,
    missing_covars = missing_covars,
    dropped_proteins = drop_prot,
    split_counts = as.data.frame(table(df$analysis_set, df$instance, useNA = "ifany"))
  )
}

cross_sectional_metrics <- function(score_tbl, enc, exposure_id) {
  models <- c(pred_prot = "prot_only", pred_cov = "cov_only", pred_full = "prot_plus_cov")
  bind_rows(lapply(names(models), function(pred_col) {
    bind_rows(lapply(split(score_tbl, score_tbl$instance), function(dd) {
      metric_row(dd$y_raw, enc, dd[[pred_col]]) %>%
        mutate(instance = unique(dd$instance), model = models[[pred_col]])
    }))
  })) %>%
    mutate(exposure_id = exposure_id) %>%
    relocate(exposure_id, instance, model)
}

change_metrics_continuous <- function(pair_df, pred_col, enc) {
  dy <- outcome_code(pair_df$y_follow, enc) - outcome_code(pair_df$y_0, enc)
  dp <- pair_df[[paste0(pred_col, "_follow")]] - pair_df[[paste0(pred_col, "_0")]]
  ok <- is.finite(dy) & is.finite(dp)
  dy <- dy[ok]
  dp <- dp[ok]
  if (length(dy) == 0) {
    return(tibble(n = 0, delta_cor = NA_real_, delta_rmse = NA_real_, beta_delta_pes = NA_real_, se_delta_pes = NA_real_, p_delta_pes = NA_real_))
  }
  fit <- if (length(dy) >= 3 && stats::sd(dp) > 0) summary(lm(dy ~ dp)) else NULL
  tibble(
    n = length(dy),
    delta_cor = ifelse(length(dy) > 2 && stats::sd(dy) > 0 && stats::sd(dp) > 0, stats::cor(dy, dp), NA_real_),
    delta_rmse = sqrt(mean((dy - dp)^2)),
    beta_delta_pes = if (!is.null(fit)) fit$coefficients["dp", "Estimate"] else NA_real_,
    se_delta_pes = if (!is.null(fit)) fit$coefficients["dp", "Std. Error"] else NA_real_,
    p_delta_pes = if (!is.null(fit)) fit$coefficients["dp", "Pr(>|t|)"] else NA_real_
  )
}

change_transition_summary <- function(pair_df, pred_col, enc) {
  y0 <- outcome_code(pair_df$y_0, enc)
  yf <- outcome_code(pair_df$y_follow, enc)
  p0 <- pair_df[[paste0(pred_col, "_0")]]
  pf <- pair_df[[paste0(pred_col, "_follow")]]
  transition <- ifelse(is.na(y0) | is.na(yf), NA_character_, paste0(y0, "->", yf))
  tibble(
    transition = transition,
    pred_0 = p0,
    pred_follow = pf,
    delta_pred = pf - p0
  ) %>%
    filter(!is.na(transition), is.finite(delta_pred)) %>%
    group_by(transition) %>%
    summarise(
      n = n(),
      mean_pred_0 = mean(pred_0, na.rm = TRUE),
      mean_pred_follow = mean(pred_follow, na.rm = TRUE),
      mean_delta_pred = mean(delta_pred, na.rm = TRUE),
      sd_delta_pred = sd(delta_pred, na.rm = TRUE),
      median_delta_pred = median(delta_pred, na.rm = TRUE),
      .groups = "drop"
    )
}

within_person_change_eval <- function(score_tbl, enc, exposure_id) {
  models <- c(pred_prot = "prot_only", pred_cov = "cov_only", pred_full = "prot_plus_cov")
  continuous_family <- enc$family == "gaussian"

  out <- list(metrics = list(), transitions = list())
  for (follow_inst in c(2L, 3L)) {
    base <- score_tbl %>%
      filter(instance == 0L) %>%
      select(eid, y_0 = y_raw, pred_prot_0 = pred_prot, pred_cov_0 = pred_cov, pred_full_0 = pred_full)
    follow <- score_tbl %>%
      filter(instance == follow_inst) %>%
      select(eid, y_follow = y_raw, pred_prot_follow = pred_prot, pred_cov_follow = pred_cov, pred_full_follow = pred_full)
    pair_df <- inner_join(base, follow, by = "eid")
    if (nrow(pair_df) == 0) next

    for (pred_col in names(models)) {
      if (continuous_family) {
        out$metrics[[paste(follow_inst, pred_col, sep = "_")]] <- change_metrics_continuous(pair_df, pred_col, enc) %>%
          mutate(exposure_id = exposure_id, followup_instance = follow_inst, model = models[[pred_col]], metric_note = NA_character_)
      } else {
        out$transitions[[paste(follow_inst, pred_col, sep = "_")]] <- change_transition_summary(pair_df, pred_col, enc) %>%
          mutate(exposure_id = exposure_id, followup_instance = follow_inst, model = models[[pred_col]])
      }
    }
  }

  if (continuous_family && length(out$metrics) == 0) {
    out$metrics <- list(tibble(
      exposure_id = exposure_id,
      followup_instance = integer(),
      model = character(),
      n = integer(),
      delta_cor = numeric(),
      delta_rmse = numeric(),
      beta_delta_pes = numeric(),
      se_delta_pes = numeric(),
      p_delta_pes = numeric(),
      metric_note = character()
    ))
  }
  if (!continuous_family) {
    out$metrics <- list(tibble(
      exposure_id = exposure_id,
      followup_instance = integer(),
      model = character(),
      metric_note = character()
    ))
  }
  if (length(out$transitions) == 0) {
    out$transitions <- list(tibble(
      exposure_id = exposure_id,
      followup_instance = integer(),
      model = character(),
      transition = character(),
      n = integer(),
      mean_pred_0 = numeric(),
      mean_pred_follow = numeric(),
      mean_delta_pred = numeric(),
      sd_delta_pred = numeric(),
      median_delta_pred = numeric()
    ))
  }

  list(
    change_metrics = bind_rows(out$metrics) %>% relocate(exposure_id, followup_instance, model),
    transition_summary = bind_rows(out$transitions) %>% relocate(exposure_id, followup_instance, model, transition)
  )
}

main <- function() {
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop(paste0(
    "Usage:\n",
    "  Rscript Module6_prod_longitudinal.R <covarType> <exposure_id>\n\n",
    "Example:\n",
    "  Rscript Module6_prod_longitudinal.R base_clinical smoking_status_f20116_0_0_Current\n"
  ))
}

covarType <- trimws(as.character(args[1]))  # guard against stray trailing whitespace in COVAR_TYPE
exposure_id_input <- as.character(args[2])
exposure_id <- canonicalize_exposure_id(exposure_id_input)
if (!identical(exposure_id, exposure_id_input)) {
  message_ts("Canonicalized exposure id ", exposure_id_input, " -> ", exposure_id)
}

if (!covarType %in% names(cfg$CovarSpec)) stop("Unknown covariate_set: ", covarType)

if (!file.exists(cfg$paths$exposure_manifest)) {
  stop("Missing exposure manifest TSV: ", cfg$paths$exposure_manifest)
}
exposure_specs <- fread(cfg$paths$exposure_manifest)
if (!all(c("exposure_id", "exposure_type") %in% names(exposure_specs))) {
  stop("Exposure manifest missing required columns exposure_id, exposure_type")
}
if (!exposure_id %in% exposure_specs$exposure_id) {
  stop("exposure_id not found in exposure manifest after canonicalization: ", exposure_id)
}
exposure_type <- as.character(exposure_specs$exposure_type[match(exposure_id, exposure_specs$exposure_id)])
if (!exposure_type %in% c("continuous", "binary", "ordinal", "categorical")) {
  stop("Invalid exposure_type for ", exposure_id, ": ", exposure_type)
}

message_ts("Loading HEAP.rds + deriving longitudinal PXS: ", cfg$paths$heap_rds)
pxs0 <- as_long_pxs(as_pxs_longitudinal(readRDS(cfg$paths$heap_rds)))

covars_used <- cfg$CovarSpec[[covarType]]
if (is.null(covars_used) || length(covars_used) == 0) stop("No covariates found for covarType: ", covarType)

# Sample filter: exclude-prevalent covar_types run on the healthy-at-baseline subset.
sample_filter <- if (covarType %in% EXCLUDE_PREVALENT_COVAR_TYPES) "exclude_prevalent" else "none"
message_ts("Sample filter: ", sample_filter)

covar_out_dir <- file.path(cfg$out_dir, covarType)
dir.create(covar_out_dir, showWarnings = FALSE, recursive = TRUE)
out_prefix <- file.path(covar_out_dir, paste0("PESlong_", covarType, "_", safe_name(exposure_id)))

message_ts("Assembling train and longitudinal hold-out data for ", exposure_id, " (", exposure_type, ")")
assembled <- assemble_for_exposure_long(
  pxs = pxs0,
  exposure_id = exposure_id,
  covars_used = covars_used,
  miss_rate_prot = cfg$miss_rate_prot,
  sample_filter = sample_filter
)

message_ts("Training rows: ", nrow(assembled$train_df), " | hold-out rows: ", nrow(assembled$holdout_df),
           " | proteins used: ", length(assembled$prot_cols))
if (length(assembled$missing_covars) > 0) {
  message_ts("Covariates absent from loader and omitted: ", paste(assembled$missing_covars, collapse = ", "))
}

oof_fit <- oof_predictors_glmnet_long(
  df = assembled$train_df,
  exposure_id = exposure_id,
  exposure_type = exposure_type,
  prot_cols = assembled$prot_cols,
  covars_used = assembled$covars_used,
  kfold = cfg$kfold,
  seed = cfg$seed
)

message_ts("Fitting final models on full training set")
artifact <- fit_final_artifact(
  df = assembled$train_df,
  exposure_id = exposure_id,
  exposure_type = exposure_type,
  prot_cols = assembled$prot_cols,
  covars_used = assembled$covars_used,
  oof_fit = oof_fit
)

holdout_scores <- score_artifact(assembled$holdout_df, artifact) %>%
  mutate(
    exposure_id = exposure_id,
    exposure_type = exposure_type,
    pes_prot_z_finalscale = zscore_with(pred_prot, artifact$z_params$prot$center, artifact$z_params$prot$scale),
    pes_full_z_finalscale = zscore_with(pred_full, artifact$z_params$full$center, artifact$z_params$full$scale),
    pes_prot_z_oofscale = zscore_with(pred_prot, artifact$oof_z_params$prot$center, artifact$oof_z_params$prot$scale),
    pes_full_z_oofscale = zscore_with(pred_full, artifact$oof_z_params$full$center, artifact$oof_z_params$full$scale),
    # Backward-compatible names retain the historical final-model training scale.
    pes_prot_z = pes_prot_z_finalscale,
    pes_full_z = pes_full_z_finalscale
  ) %>%
  relocate(eid, instance, exposure_id, exposure_type, y_raw)

cross_metrics <- cross_sectional_metrics(holdout_scores, artifact$encoder, exposure_id)
change_eval <- within_person_change_eval(holdout_scores, artifact$encoder, exposure_id)

artifact$run_config <- cfg
artifact$data_summary <- list(
  train_n = nrow(assembled$train_df),
  holdout_n = nrow(assembled$holdout_df),
  prot_n = length(assembled$prot_cols),
  covars_used = assembled$covars_used,
  missing_covars = assembled$missing_covars,
  dropped_proteins = assembled$dropped_proteins
)
artifact$pes_scale_contract <- list(
  train_oof_columns = c(prot = "pes_prot_z", full = "pes_full_z"),
  holdout_default_columns = c(prot = "pes_prot_z", full = "pes_full_z"),
  holdout_finalscale_columns = c(prot = "pes_prot_z_finalscale", full = "pes_full_z_finalscale"),
  holdout_oofscale_columns = c(prot = "pes_prot_z_oofscale", full = "pes_full_z_oofscale"),
  note = paste(
    "TrainOOF pes_*_z columns are scaled with OOF prediction moments.",
    "Holdout pes_*_z columns are retained as final-model scale for backward compatibility.",
    "For Cox models trained on OOF PES, use holdout *_oofscale columns for strict scale matching."
  )
)

saveRDS(oof_fit$oof, paste0(out_prefix, "_TrainOOF.rds"))
fwrite(as.data.table(oof_fit$oof), paste0(out_prefix, "_TrainOOF.tsv"), sep = "\t")
fwrite(as.data.table(oof_fit$fold_metrics), paste0(out_prefix, "_TrainFoldMetrics.tsv"), sep = "\t")
fwrite(as.data.table(oof_fit$overall_metrics), paste0(out_prefix, "_TrainOverallMetrics.tsv"), sep = "\t")
saveRDS(artifact, paste0(out_prefix, "_FinalModelArtifact.rds"))
fwrite(as.data.table(holdout_scores), paste0(out_prefix, "_HoldoutScores.tsv"), sep = "\t")
saveRDS(holdout_scores, paste0(out_prefix, "_HoldoutScores.rds"))
fwrite(as.data.table(cross_metrics), paste0(out_prefix, "_HoldoutByVisitMetrics.tsv"), sep = "\t")
fwrite(as.data.table(change_eval$change_metrics), paste0(out_prefix, "_WithinPersonChangeMetrics.tsv"), sep = "\t")
fwrite(as.data.table(change_eval$transition_summary), paste0(out_prefix, "_WithinPersonTransitionSummary.tsv"), sep = "\t")

selected_tbl <- tibble(
  model = c("prot_only", "prot_plus_cov"),
  n_selected_proteins = c(length(artifact$selected_proteins$prot), length(artifact$selected_proteins$full)),
  selected_proteins = c(
    paste(artifact$selected_proteins$prot, collapse = ";"),
    paste(artifact$selected_proteins$full, collapse = ";")
  )
)
fwrite(as.data.table(selected_tbl), paste0(out_prefix, "_SelectedProteins.tsv"), sep = "\t")

message_ts("Saved outputs with prefix: ", out_prefix)
message_ts("DONE")
}

if (identical(environment(), globalenv()) && Sys.getenv("M6_LONG_SOURCE_ONLY") != "1") {
  main()
}
