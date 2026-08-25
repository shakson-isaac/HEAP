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

# Build compact/deployable Module 6 PES scores from existing longitudinal final
# model artifacts. This does not refit PES models; it freezes the baseline-trained
# LASSO weights, ranks proteins by final-model absolute coefficient, and scores
# held-out repeat visits using top-K protein panels.

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(glmnet)
  library(pROC)
  library(tibble)
})

cfg <- list(
  heap_rds = heap_loader_rds,  # read HEAP.rds directly; derive longitudinal PXS via as_pxs_longitudinal()
  result_dir = heap_project_output("module6_pes_longitudinal"),
  exposure_manifest = heap_project_output("module6_pes_test", "exposure_specs.tsv"),
  portability_file = igloo_path("UKB", "OlinkSoma", "OlinkSoma.csv"),
  out_tag = "CompactPES",
  ks = c(10, 25, 50, 100, 200, 500),
  selection_modes = c("lasso", "portable", "portable_weighted"),
  model_types = c("prot"),
  instances = c(0L, 2L, 3L),
  portability_corr_col = "olink_smpnorm_corr",
  portability_threshold = 0.70,
  portability_thresholds = c(0.50, 0.60, 0.70, 0.80, 0.90),
  portability_gamma = 1,
  include_train = TRUE
)

`%||%` <- function(x, y) if (!is.null(x)) x else y

message_ts <- function(...) {
  message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), ...)
}

safe_name <- function(x) gsub("[^A-Za-z0-9_.-]+", "_", x)

canonicalize_exposure_id <- function(x) {
  gsub("(_f[0-9]+)_[23]_", "\\1_0_", x)
}

parse_csv <- function(x, numeric = FALSE) {
  if (is.null(x) || !nzchar(x)) return(character())
  out <- trimws(unlist(strsplit(x, ",")))
  out <- out[nzchar(out)]
  if (numeric) return(as.numeric(out))
  out
}

parse_args <- function(args) {
  if (length(args) < 2) {
    stop(paste0(
      "Usage:\n",
      "  Rscript Module6_compact_pes_deployability.R <covarType> <exposure_id> [options]\n\n",
      "Common options:\n",
      "  --ks 10,50,100,200,500,all\n",
      "  --selection-modes lasso,portable,portable_weighted\n",
      "  --model-types prot\n",
      "  --portability-file /n/groups/patel/IGLOO/UKB/OlinkSoma/OlinkSoma.csv\n",
      "  --portability-threshold 0.7\n",
      "  --portability-thresholds 0.5,0.6,0.7,0.8,0.9\n",
      "  --output-tag CompactPES\n",
      "  --save-scores            (off by default; writes the large ~0.8GB/exposure per-person Scores.tsv)\n"
    ))
  }
  out <- list(covar_type = args[[1]], exposure_id = args[[2]])
  i <- 3
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) stop("Unexpected positional argument: ", key)
    key <- sub("^--", "", key)
    if (key %in% c("include-train", "no-include-train", "save-scores")) {
      out[[gsub("-", "_", key)]] <- TRUE
      i <- i + 1
    } else {
      if (i == length(args)) stop("Missing value for --", key)
      out[[gsub("-", "_", key)]] <- args[[i + 1]]
      i <- i + 2
    }
  }
  out
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

fit_zscore_params <- function(x) {
  list(center = mean(x, na.rm = TRUE), scale = stats::sd(x, na.rm = TRUE))
}

zscore_with <- function(x, center, scale) {
  if (!is.finite(scale) || scale == 0) return(rep(0, length(x)))
  (x - center) / scale
}

clip01 <- function(p, eps = 1e-15) pmin(pmax(p, eps), 1 - eps)

apply_median_imputer <- function(X, imp) {
  X2 <- as.matrix(X)
  idx <- which(is.na(X2), arr.ind = TRUE)
  if (nrow(idx) > 0) X2[idx] <- imp$median[idx[, 2]]
  X2
}

apply_outcome_encoder <- function(y_raw, enc) {
  if (enc$family == "gaussian") return(as.numeric(y_raw))
  if (enc$kind == "binary_numeric") {
    yy <- as.numeric(y_raw)
    out <- rep(NA_integer_, length(yy))
    out[yy == enc$values[1]] <- 0L
    out[yy == enc$values[2]] <- 1L
    return(out)
  }
  if (enc$family == "binomial") {
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

read_portability_map <- function(path, corr_col, threshold) {
  if (is.null(path) || !file.exists(path)) {
    message_ts("Portability file unavailable; portability-aware modes will produce empty selections: ", path)
    return(tibble())
  }

  raw <- readLines(path, warn = FALSE, n = 25)
  header_line <- which(grepl("(^|,)seqid(,|$)", tolower(raw)))[1]
  skip <- ifelse(is.na(header_line), 0L, header_line - 1L)
  dt <- fread(path, skip = skip, showProgress = FALSE)
  names(dt) <- tolower(gsub("[^A-Za-z0-9_]+", "_", names(dt)))

  corr_col2 <- tolower(gsub("[^A-Za-z0-9_]+", "_", corr_col))
  if (!corr_col2 %in% names(dt)) {
    stop("Portability correlation column not found: ", corr_col, ". Available columns: ", paste(names(dt), collapse = ", "))
  }

  as_tibble(dt) %>%
    mutate(
      portability_corr = suppressWarnings(as.numeric(.data[[corr_col2]])),
      portability_abs_corr = abs(portability_corr),
      portable = is.finite(portability_abs_corr) & portability_abs_corr >= threshold,
      across(any_of(c("seqid", "oid", "uniprot", "target", "target_name", "gene_name")), as.character)
    )
}

normalize_key <- function(x) toupper(gsub("[^A-Za-z0-9]+", "", x))

match_one_portability <- function(protein, port) {
  empty <- tibble(
    protein = protein,
    portability_corr = NA_real_,
    portability_abs_corr = NA_real_,
    portable = FALSE,
    portability_match = NA_character_,
    portability_match_method = "unmatched",
    portability_uniprot = NA_character_,
    portability_target = NA_character_,
    portability_gene_name = NA_character_
  )
  if (nrow(port) == 0) return(empty)

  pkey <- normalize_key(protein)
  candidates <- bind_rows(
    port %>% mutate(.match_key = normalize_key(gene_name %||% NA_character_), .method = "gene_name"),
    port %>% mutate(.match_key = normalize_key(target %||% NA_character_), .method = "target"),
    port %>% mutate(.match_key = normalize_key(uniprot %||% NA_character_), .method = "uniprot"),
    port %>% mutate(.match_key = normalize_key(seqid %||% NA_character_), .method = "seqid")
  ) %>%
    filter(!is.na(.match_key), nzchar(.match_key), .match_key == pkey)

  if (nrow(candidates) == 0 && grepl("_", protein, fixed = TRUE)) {
    tokens <- normalize_key(unlist(strsplit(protein, "_", fixed = TRUE)))
    candidates <- bind_rows(
      port %>% mutate(.match_key = normalize_key(gene_name %||% NA_character_), .method = "gene_name_token"),
      port %>% mutate(.match_key = normalize_key(target %||% NA_character_), .method = "target_token")
    ) %>%
      filter(!is.na(.match_key), nzchar(.match_key), .match_key %in% tokens)
  }

  if (nrow(candidates) == 0) return(empty)
  best <- candidates %>%
    arrange(desc(portability_abs_corr), desc(abs(portability_corr))) %>%
    slice(1)

  tibble(
    protein = protein,
    portability_corr = best$portability_corr,
    portability_abs_corr = best$portability_abs_corr,
      portable = is.finite(best$portability_abs_corr),
    portability_match = best$seqid %||% NA_character_,
    portability_match_method = best$.method,
    portability_uniprot = best$uniprot %||% NA_character_,
    portability_target = best$target %||% NA_character_,
    portability_gene_name = best$gene_name %||% NA_character_
  )
}

coef_table <- function(fit, family, lambda, prot_cols, model_type) {
  if (model_type == "full") {
    model_type_label <- "full_protein_component"
  } else {
    model_type_label <- "prot"
  }

  if (family != "multinomial") {
    b <- as.matrix(coef(fit, s = lambda))
    tibble(
      model_type = model_type_label,
      protein = rownames(b),
      class = NA_character_,
      beta = as.numeric(b[, 1])
    ) %>%
      filter(protein %in% prot_cols, beta != 0) %>%
      mutate(abs_beta = abs(beta))
  } else {
    coefs <- coef(fit, s = lambda)
    bind_rows(lapply(names(coefs), function(cls) {
      cm <- as.matrix(coefs[[cls]])
      tibble(
        model_type = model_type_label,
        protein = rownames(cm),
        class = cls,
        beta = as.numeric(cm[, 1])
      )
    })) %>%
      filter(protein %in% prot_cols, beta != 0) %>%
      group_by(model_type, protein) %>%
      summarise(
        beta = NA_real_,
        abs_beta = sqrt(sum(beta^2, na.rm = TRUE)),
        beta_multiclass = paste(paste(class, signif(beta, 6), sep = ":"), collapse = ";"),
        .groups = "drop"
      )
  }
}

rank_proteins <- function(coefs, port, modes, gamma, thresholds) {
  if (nrow(coefs) == 0) return(tibble())
  port_tbl <- bind_rows(lapply(coefs$protein, match_one_portability, port = port))
  base <- coefs %>%
    left_join(port_tbl, by = "protein") %>%
    mutate(
      portability_abs_corr = ifelse(is.finite(portability_abs_corr), portability_abs_corr, NA_real_),
      has_portability_match = is.finite(portability_abs_corr),
      portability_weight = ifelse(is.finite(portability_abs_corr), portability_abs_corr^gamma, 0),
      weighted_abs_beta = abs_beta * portability_weight
    )

  bind_rows(lapply(modes, function(mode) {
    if (mode == "lasso") {
      return(base %>%
        arrange(desc(abs_beta), protein) %>%
        mutate(
          selection_mode = mode,
          portability_threshold = NA_real_,
          portable = has_portability_match,
          rank = row_number()
        ) %>%
        relocate(selection_mode, portability_threshold, model_type, protein, rank))
    }

    if (!mode %in% c("portable", "portable_weighted")) stop("Unknown selection mode: ", mode)

    bind_rows(lapply(thresholds, function(thr) {
      ranked <- base %>%
        mutate(portable = has_portability_match & portability_abs_corr >= thr) %>%
        filter(portable)
      if (mode == "portable") {
        ranked <- ranked %>% arrange(desc(abs_beta), protein)
      } else {
        ranked <- ranked %>% arrange(desc(weighted_abs_beta), protein)
      }
      ranked %>%
        mutate(selection_mode = mode, portability_threshold = thr, rank = row_number()) %>%
        relocate(selection_mode, portability_threshold, model_type, protein, rank)
    }))
  }))
}

subset_coef_list <- function(fit, family, lambda, keep, prot_cols) {
  if (family != "multinomial") {
    b <- as.matrix(coef(fit, s = lambda))
    b <- b[c("(Intercept)", intersect(keep, rownames(b))), , drop = FALSE]
    return(b)
  }
  coefs <- coef(fit, s = lambda)
  lapply(coefs, function(cm) {
    m <- as.matrix(cm)
    m[c("(Intercept)", intersect(keep, rownames(m))), , drop = FALSE]
  })
}

predict_from_subset <- function(X, coef_obj, enc) {
  if (enc$family != "multinomial") {
    terms <- setdiff(rownames(coef_obj), "(Intercept)")
    lp <- rep(as.numeric(coef_obj["(Intercept)", 1]), nrow(X))
    if (length(terms) > 0) lp <- lp + as.numeric(X[, terms, drop = FALSE] %*% coef_obj[terms, 1])
    if (enc$family == "binomial") return(list(linear = lp, pred = plogis(lp)))
    return(list(linear = lp, pred = lp))
  }

  levels <- enc$levels
  eta <- matrix(0, nrow = nrow(X), ncol = length(levels), dimnames = list(NULL, levels))
  for (cls in levels) {
    cm <- coef_obj[[cls]]
    if (is.null(cm)) next
    terms <- setdiff(rownames(cm), "(Intercept)")
    eta[, cls] <- as.numeric(cm["(Intercept)", 1])
    if (length(terms) > 0) eta[, cls] <- eta[, cls] + as.numeric(X[, terms, drop = FALSE] %*% cm[terms, 1])
  }
  eta <- eta - apply(eta, 1, max)
  prob <- exp(eta)
  prob <- prob / rowSums(prob)
  pred <- as.numeric(prob %*% matrix(enc$codes[colnames(prob)], ncol = 1))
  list(linear = pred, pred = pred)
}

metric_row <- function(y_raw, enc, pred, score = pred) {
  if (enc$family == "gaussian") {
    y <- apply_outcome_encoder(y_raw, enc)
    ok <- is.finite(y) & is.finite(pred)
    y <- y[ok]
    p <- pred[ok]
    ss_res <- sum((y - p)^2)
    ss_tot <- sum((y - mean(y))^2)
    return(tibble(
      n_eval = length(y),
      r2 = ifelse(ss_tot == 0, NA_real_, 1 - ss_res / ss_tot),
      correlation = ifelse(length(y) > 2 && stats::sd(y) > 0 && stats::sd(p) > 0, stats::cor(y, p), NA_real_),
      rmse = sqrt(mean((y - p)^2)),
      auc = NA_real_,
      logloss = NA_real_,
      correlation_code = NA_real_,
      rmse_code = NA_real_
    ))
  }

  if (enc$family == "binomial") {
    y <- apply_outcome_encoder(y_raw, enc)
    ok <- !is.na(y) & is.finite(pred) & is.finite(score)
    y <- y[ok]
    p <- clip01(pred[ok])
    s <- score[ok]
    return(tibble(
      n_eval = length(y),
      r2 = NA_real_,
      correlation = NA_real_,
      rmse = NA_real_,
      auc = ifelse(length(unique(y)) == 2, as.numeric(pROC::auc(pROC::roc(y, s, quiet = TRUE))), NA_real_),
      logloss = -mean(y * log(p) + (1 - y) * log(1 - p)),
      correlation_code = NA_real_,
      rmse_code = NA_real_
    ))
  }

  y_code <- outcome_code(y_raw, enc)
  ok <- is.finite(y_code) & is.finite(pred)
  y_code <- y_code[ok]
  p <- pred[ok]
  tibble(
    n_eval = length(y_code),
    r2 = NA_real_,
    correlation = NA_real_,
    rmse = NA_real_,
    auc = NA_real_,
    logloss = NA_real_,
    correlation_code = ifelse(length(y_code) > 2 && stats::sd(y_code) > 0 && stats::sd(p) > 0, stats::cor(y_code, p), NA_real_),
    rmse_code = sqrt(mean((y_code - p)^2))
  )
}

assemble_rows <- function(pxs, exposure_id, artifact, include_train, instances) {
  pxs <- as_long_pxs(pxs)
  join_key <- c("eid", "instance")
  E_df <- pxs$Elist %>% purrr::reduce(full_join, by = join_key)
  if (!exposure_id %in% names(E_df)) stop("Exposure not found in longitudinal Elist: ", exposure_id)

  prot_df <- as.data.frame(pxs$UKBprot_df)
  split_df <- as.data.frame(pxs$split_df %||% data.frame(eid = unique(prot_df$eid)))
  missing_prot <- setdiff(artifact$prot_cols, names(prot_df))
  if (length(missing_prot) > 0) stop("Proteins from artifact missing in loader: ", paste(head(missing_prot, 20), collapse = ", "))

  df <- prot_df[, c(join_key, artifact$prot_cols), drop = FALSE] %>%
    inner_join(E_df[, c(join_key, exposure_id), drop = FALSE], by = join_key) %>%
    left_join(split_df, by = "eid")

  keep <- (df$analysis_set == "holdout_repeat_proteomics" & df$instance %in% instances)
  if (include_train) keep <- keep | (df$analysis_set == "train_baseline_only" & df$instance == 0L)
  df[keep, , drop = FALSE]
}

score_compact_panels <- function(df, artifact, ranked, ks, model_type) {
  enc <- artifact$encoder
  fit <- if (model_type == "prot") artifact$fit_prot else artifact$fit_full
  lambda <- if (model_type == "prot") artifact$lambda$prot else artifact$lambda$full
  X <- apply_median_imputer(as.matrix(df[, artifact$prot_cols, drop = FALSE]), artifact$protein_imputer)

  full_coef <- subset_coef_list(fit, enc$family, lambda, keep = artifact$prot_cols, prot_cols = artifact$prot_cols)
  full_pred <- predict_from_subset(X, full_coef, enc)
  full_z_params <- fit_zscore_params(full_pred$pred[df$analysis_set == "train_baseline_only" & df$instance == 0L])

  score_rows <- list()
  protein_rows <- list()
  ranked <- ranked %>%
    filter(model_type == ifelse(model_type == "prot", "prot", "full_protein_component"))
  rank_keys <- ranked %>%
    distinct(selection_mode, portability_threshold) %>%
    arrange(selection_mode, portability_threshold)

  for (ii in seq_len(nrow(rank_keys))) {
    mode <- rank_keys$selection_mode[ii]
    thr <- rank_keys$portability_threshold[ii]
    mode_ranked <- ranked %>%
      filter(selection_mode == mode) %>%
      filter(if (is.na(thr)) is.na(portability_threshold) else portability_threshold == thr) %>%
      arrange(rank)
    if (nrow(mode_ranked) == 0) next

    max_rank <- max(mode_ranked$rank, na.rm = TRUE)
    if (!is.finite(max_rank)) next

    for (req_k in unique(ks)) {
      actual_k <- min(ifelse(is.na(req_k), max_rank, req_k), max_rank)
      requested_k <- ifelse(is.na(req_k), "all", as.character(as.integer(req_k)))
      keep <- head(mode_ranked$protein, actual_k)
      coef_k <- subset_coef_list(fit, enc$family, lambda, keep = keep, prot_cols = artifact$prot_cols)
      pred_k <- predict_from_subset(X, coef_k, enc)
      train_idx <- df$analysis_set == "train_baseline_only" & df$instance == 0L
      z_k <- fit_zscore_params(pred_k$pred[train_idx])

      row_key <- paste(mode, ifelse(is.na(thr), "NA", thr), requested_k, sep = "_")
      score_rows[[row_key]] <- tibble(
        eid = df$eid,
        instance = df$instance,
        analysis_set = df$analysis_set,
        exposure_id = artifact$exposure_id,
        exposure_type = artifact$exposure_type,
        y_raw = df[[artifact$exposure_id]],
        model_type = ifelse(model_type == "prot", "prot", "full_protein_component"),
        selection_mode = mode,
        portability_threshold = thr,
        requested_k = requested_k,
        actual_k = as.integer(actual_k),
        k = as.integer(actual_k),
        compact_pred = pred_k$pred,
        compact_linear = pred_k$linear,
        compact_pes_z = zscore_with(pred_k$pred, z_k$center, z_k$scale),
        full_model_pred = full_pred$pred,
        full_model_pes_z = zscore_with(full_pred$pred, full_z_params$center, full_z_params$scale)
      )

      protein_rows[[row_key]] <- mode_ranked %>%
        slice_head(n = actual_k) %>%
        mutate(
          requested_k = requested_k,
          actual_k = as.integer(actual_k),
          k = as.integer(actual_k),
          selected_in_panel = TRUE
        )
    }
  }

  list(scores = bind_rows(score_rows), proteins = bind_rows(protein_rows))
}

cross_sectional_metrics <- function(scores, enc) {
  scores$.threshold_key <- ifelse(is.na(scores$portability_threshold), "none", as.character(scores$portability_threshold))
  split(scores, list(scores$model_type, scores$selection_mode, scores$.threshold_key, scores$requested_k, scores$actual_k, scores$analysis_set, scores$instance), drop = TRUE) %>%
    lapply(function(dd) {
      metric_row(dd$y_raw, enc, dd$compact_pred, dd$compact_linear) %>%
        mutate(
          model_type = unique(dd$model_type),
          selection_mode = unique(dd$selection_mode),
          portability_threshold = unique(dd$portability_threshold),
          requested_k = unique(dd$requested_k),
          actual_k = unique(dd$actual_k),
          k = unique(dd$k),
          analysis_set = unique(dd$analysis_set),
          instance = unique(dd$instance),
          score_cor_full_pes = ifelse(
            nrow(dd) > 2 && sd(dd$compact_pes_z, na.rm = TRUE) > 0 && sd(dd$full_model_pes_z, na.rm = TRUE) > 0,
            cor(dd$compact_pes_z, dd$full_model_pes_z, use = "pairwise.complete.obs"),
            NA_real_
          )
        )
    }) %>%
    bind_rows() %>%
    relocate(model_type, selection_mode, portability_threshold, requested_k, actual_k, k, analysis_set, instance, n_eval, score_cor_full_pes)
}

delta_metrics <- function(scores, enc) {
  hold <- scores %>% filter(analysis_set == "holdout_repeat_proteomics")
  if (nrow(hold) == 0) return(tibble())
  keys <- hold %>% distinct(model_type, selection_mode, portability_threshold, requested_k, actual_k, k)
  out <- list()
  for (ii in seq_len(nrow(keys))) {
    kk <- keys[ii, ]
    dd <- hold %>%
      filter(
        model_type == kk$model_type,
        selection_mode == kk$selection_mode,
        if (is.na(kk$portability_threshold)) is.na(portability_threshold) else portability_threshold == kk$portability_threshold,
        requested_k == kk$requested_k,
        actual_k == kk$actual_k
      )
    for (pair in list(c(0L, 2L), c(0L, 3L), c(2L, 3L))) {
      a <- pair[1]
      b <- pair[2]
      wide <- dd %>%
        filter(instance %in% c(a, b)) %>%
        select(eid, instance, y_raw, compact_pes_z, full_model_pes_z) %>%
        pivot_wider(
          names_from = instance,
          values_from = c(y_raw, compact_pes_z, full_model_pes_z),
          values_fn = list(y_raw = dplyr::first, compact_pes_z = dplyr::first, full_model_pes_z = dplyr::first)
        )
      ca <- paste0("compact_pes_z_", a)
      cb <- paste0("compact_pes_z_", b)
      fa <- paste0("full_model_pes_z_", a)
      fb <- paste0("full_model_pes_z_", b)
      ya <- paste0("y_raw_", a)
      yb <- paste0("y_raw_", b)
      if (!all(c(ca, cb, fa, fb, ya, yb) %in% names(wide))) next
      dy <- outcome_code(wide[[yb]], enc) - outcome_code(wide[[ya]], enc)
      dc <- wide[[cb]] - wide[[ca]]
      df <- wide[[fb]] - wide[[fa]]
      ok <- is.finite(dc) & is.finite(df)
      ok_y <- ok & is.finite(dy)
      out[[paste(ii, a, b, sep = "_")]] <- tibble(
        model_type = kk$model_type,
        selection_mode = kk$selection_mode,
        portability_threshold = kk$portability_threshold,
        requested_k = kk$requested_k,
        actual_k = kk$actual_k,
        k = kk$k,
        visit_pair = paste0("i", a, "_to_i", b),
        n_pair = sum(ok),
        delta_cor_full_pes = ifelse(sum(ok) > 2 && sd(dc[ok]) > 0 && sd(df[ok]) > 0, cor(dc[ok], df[ok]), NA_real_),
        delta_cor_exposure = ifelse(sum(ok_y) > 2 && sd(dc[ok_y]) > 0 && sd(dy[ok_y]) > 0, cor(dc[ok_y], dy[ok_y]), NA_real_),
        mean_abs_delta_compact = mean(abs(dc[ok]), na.rm = TRUE),
        mean_abs_delta_full = mean(abs(df[ok]), na.rm = TRUE)
      )
    }
  }
  bind_rows(out)
}

main <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))
  covar_type <- args$covar_type
  exposure_id <- canonicalize_exposure_id(args$exposure_id)

  result_dir <- args$result_dir %||% cfg$result_dir
  pxs_rds <- args$loader_rds %||% cfg$heap_rds
  portability_file <- args$portability_file %||% cfg$portability_file
  out_tag <- args$output_tag %||% cfg$out_tag
  ks_raw <- parse_csv(args$ks %||% paste(cfg$ks, collapse = ","))
  ks <- suppressWarnings(as.numeric(ifelse(tolower(ks_raw) == "all", NA_character_, ks_raw)))
  modes <- parse_csv(args$selection_modes %||% paste(cfg$selection_modes, collapse = ","))
  model_types <- parse_csv(args$model_types %||% paste(cfg$model_types, collapse = ","))
  instances <- as.integer(parse_csv(args$instances %||% paste(cfg$instances, collapse = ","), numeric = TRUE))
  corr_col <- args$portability_corr_col %||% cfg$portability_corr_col
  port_threshold <- as.numeric(args$portability_threshold %||% cfg$portability_threshold)
  port_thresholds <- parse_csv(args$portability_thresholds %||% "", numeric = TRUE)
  if (length(port_thresholds) == 0 || all(is.na(port_thresholds))) {
    port_thresholds <- port_threshold
  }
  port_thresholds <- sort(unique(port_thresholds[is.finite(port_thresholds)]))
  if (length(port_thresholds) == 0) stop("No valid portability thresholds supplied.")
  port_gamma <- as.numeric(args$portability_gamma %||% cfg$portability_gamma)
  include_train <- !isTRUE(args$no_include_train)

  covar_out_dir <- file.path(result_dir, covar_type)
  out_prefix <- file.path(covar_out_dir, paste0("PESlong_", covar_type, "_", safe_name(exposure_id)))
  artifact_file <- paste0(out_prefix, "_FinalModelArtifact.rds")
  if (!file.exists(artifact_file)) stop("Missing final artifact: ", artifact_file)
  if (!file.exists(pxs_rds)) stop("Missing HEAP.rds: ", pxs_rds, " (run run_HEAP_loader.sh first)")

  message_ts("Loading final PES artifact: ", artifact_file)
  artifact <- readRDS(artifact_file)
  if (!identical(artifact$exposure_id, exposure_id)) {
    stop("Artifact exposure_id mismatch: expected ", exposure_id, ", got ", artifact$exposure_id)
  }

  message_ts("Reading portability map")
  port <- read_portability_map(portability_file, corr_col = corr_col, threshold = min(port_thresholds))

  all_ranked <- list()
  for (mt in model_types) {
    fit <- if (mt == "prot") artifact$fit_prot else if (mt == "full") artifact$fit_full else stop("Unsupported model type: ", mt)
    lambda <- if (mt == "prot") artifact$lambda$prot else artifact$lambda$full
    coefs <- coef_table(fit, artifact$encoder$family, lambda, artifact$prot_cols, mt)
    all_ranked[[mt]] <- rank_proteins(coefs, port, modes = modes, gamma = port_gamma, thresholds = port_thresholds)
  }
  ranked <- bind_rows(all_ranked)
  if (nrow(ranked) == 0) stop("No non-zero proteins found after ranking.")

  message_ts("Loading HEAP.rds and assembling score rows")
  pxs <- as_pxs_longitudinal(readRDS(pxs_rds))
  df <- assemble_rows(pxs, exposure_id, artifact, include_train = include_train, instances = instances)
  message_ts("Rows to score: ", nrow(df), " | selected proteins available: ", length(unique(ranked$protein)))

  scored <- list()
  panels <- list()
  for (mt in model_types) {
    res <- score_compact_panels(df, artifact, ranked, ks = ks, model_type = mt)
    scored[[mt]] <- res$scores
    panels[[mt]] <- res$proteins
  }
  score_tbl <- bind_rows(scored) %>% mutate(exposure_id = exposure_id, exposure_type = artifact$exposure_type)
  panel_tbl <- bind_rows(panels) %>% mutate(exposure_id = exposure_id, exposure_type = artifact$exposure_type)
  metrics_tbl <- cross_sectional_metrics(score_tbl, artifact$encoder) %>%
    mutate(exposure_id = exposure_id, exposure_type = artifact$exposure_type) %>%
    relocate(exposure_id, exposure_type)
  delta_tbl <- delta_metrics(score_tbl, artifact$encoder) %>%
    mutate(exposure_id = exposure_id, exposure_type = artifact$exposure_type) %>%
    relocate(exposure_id, exposure_type)

  summary_tbl <- panel_tbl %>%
    group_by(exposure_id, exposure_type, model_type, selection_mode, portability_threshold, requested_k, actual_k, k) %>%
    summarise(
      n_proteins = n_distinct(protein),
      n_portable = sum(portable, na.rm = TRUE),
      prop_portable = n_portable / n_proteins,
      median_abs_beta = median(abs_beta, na.rm = TRUE),
      median_portability_abs_corr = median(portability_abs_corr, na.rm = TRUE),
      min_portability_abs_corr = suppressWarnings(min(portability_abs_corr, na.rm = TRUE)),
      cumulative_abs_beta = sum(abs_beta, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      median_portability_abs_corr = ifelse(is.infinite(median_portability_abs_corr), NA_real_, median_portability_abs_corr),
      min_portability_abs_corr = ifelse(is.infinite(min_portability_abs_corr), NA_real_, min_portability_abs_corr)
    )

  for (nm in names(score_tbl)) if (is.list(score_tbl[[nm]])) score_tbl[[nm]] <- as.character(score_tbl[[nm]])
  for (nm in names(panel_tbl)) if (is.list(panel_tbl[[nm]])) panel_tbl[[nm]] <- as.character(panel_tbl[[nm]])

  # Per-person compact PES scores are large (~0.8 GB/exposure). Off by default; the
  # manuscript figures use the Metrics/DeltaMetrics/PanelSummary summaries below, not
  # these. Pass --save-scores to write them.
  if (isTRUE(args$save_scores)) {
    fwrite(as.data.table(score_tbl), paste0(out_prefix, "_", out_tag, "_Scores.tsv"), sep = "\t")
    message_ts("Saved per-person Scores.tsv (--save-scores)")
  } else {
    message_ts("Skipping per-person Scores.tsv (pass --save-scores to write it)")
  }
  fwrite(as.data.table(metrics_tbl), paste0(out_prefix, "_", out_tag, "_Metrics.tsv"), sep = "\t")
  fwrite(as.data.table(delta_tbl), paste0(out_prefix, "_", out_tag, "_DeltaMetrics.tsv"), sep = "\t")
  fwrite(as.data.table(panel_tbl), paste0(out_prefix, "_", out_tag, "_ProteinPanels.tsv"), sep = "\t")
  fwrite(as.data.table(summary_tbl), paste0(out_prefix, "_", out_tag, "_PanelSummary.tsv"), sep = "\t")

  readme <- c(
    "Module6 compact PES deployability outputs",
    paste0("exposure_id: ", exposure_id),
    paste0("covar_type: ", covar_type),
    paste0("artifact: ", artifact_file),
    paste0("portability_file: ", portability_file),
    paste0("portability_corr_col: ", corr_col),
    paste0("portability_thresholds: ", paste(port_thresholds, collapse = ",")),
    "Selection modes:",
    "  lasso = top K proteins by absolute final LASSO coefficient",
    "  portable = top K among proteins with absolute SomaScan/Olink correlation above threshold",
    "  portable_weighted = top K by abs(beta) * abs(correlation)^gamma among portable proteins",
    "Scores are compressed from the frozen final PES model. They are not new CV refits.",
    "Training rows are apparent compact scores; heldout_repeat_proteomics rows remain out-of-sample for PES training."
  )
  writeLines(readme, paste0(out_prefix, "_", out_tag, "_README.txt"))

  message_ts("Saved compact PES outputs with prefix: ", paste0(out_prefix, "_", out_tag))
  message_ts("DONE")
}

main()
