############################################################
# HEAP Module 1 -- Predictive R2 Decomposition via Score Partition
#
# Core estimand:
#   Fit one full penalized model per protein/fold:
#     P ~ C + Gcis + Gtrans + E + Gcis:E + Gtrans:E
#
#   Derive learned out-of-fold component scores:
#     C_score, Gcis_score, Gtrans_score, G_score,
#     PXS_total, GIS_total,
#     PXS_<category>, GIS_<category>
#
#   Compute unique held-out predictive R2 by score-level conditional drop:
#     unique E   = R2(C + G + E + GxE scores) - R2(C + G + GxE scores)
#     unique G   = R2(C + G + E + GxE scores) - R2(C + E + GxE scores)
#     unique GxE = R2(C + G + E + GxE scores) - R2(C + G + E scores)
#
# Notes:
#   - This avoids repeated full-feature lasso refits for every block/category.
#   - This estimates unique held-out predictive R2 attributable to learned
#     component scores, conditional on the other learned scores.
#   - It does not allocate shared variance by Shapley.
############################################################

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

# Config helpers: load_sample_filter / apply_sample_filter (sample-filter axis) and
# the covariate-set loaders. Sourced from the same workflow dir as 00_paths.R.
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
  library(tidyr)
  library(purrr)
  library(glmnet)
  library(caret)
  library(Matrix)
  library(tibble)
})

############################################################
# 0) Config
############################################################

cfg <- list(
  seed = 123,
  miss_rate = 0.2,
  kfold = as.integer(Sys.getenv("HEAP_KFOLD", unset = "5")),

  # Main mode for this rewritten script.
  decomposition_mode = Sys.getenv("HEAP_DECOMPOSITION_MODE", unset = "score_partition"),

  # Inner CV folds for cv.glmnet. Use 5 for final, 3 for faster testing.
  glmnet_inner_kfold = as.integer(Sys.getenv("GLMNET_INNER_KFOLD", unset = "5")),

  # These are now cheap because they are score-level drops.
  run_genetic_subblocks   = isTRUE(as.logical(Sys.getenv("HEAP_RUN_GENETIC_SUBBLOCKS",   "TRUE"))),
  run_exposure_categories = isTRUE(as.logical(Sys.getenv("HEAP_RUN_EXPOSURE_CATEGORIES", "TRUE"))),
  run_gxe_categories      = isTRUE(as.logical(Sys.getenv("HEAP_RUN_GXE_CATEGORIES",      "TRUE"))),

  # Sensitivity: Gene × Covariate and Exposure × Covariate interaction terms.
  # Added to the single full-model fit; unique R2 derived by score-level drop.
  run_gxc_sensitivity = isTRUE(as.logical(Sys.getenv("HEAP_RUN_GXC_SENSITIVITY", "FALSE"))),
  run_exc_sensitivity = isTRUE(as.logical(Sys.getenv("HEAP_RUN_EXC_SENSITIVITY", "FALSE"))),

  gxc_covars = local({
    raw <- Sys.getenv("HEAP_GXC_COVARS", unset = "")
    if (nzchar(raw)) trimws(strsplit(raw, ",", fixed = TRUE)[[1]])
    else c("age_when_attended_assessment_centre_f21003_0_0", "sex_f31_0_0")
  }),

  exc_covars = local({
    raw <- Sys.getenv("HEAP_EXC_COVARS", unset = "")
    if (nzchar(raw)) trimws(strsplit(raw, ",", fixed = TRUE)[[1]])
    else c("age_when_attended_assessment_centre_f21003_0_0", "sex_f31_0_0")
  }),

  paths = list(
    heap_rds      = heap_loader_rds,
    omicspred_map = heap_omicspred_or_legacy(
      "UKB_Olink_multi_ancestry_models_val_results_portal.csv"
    ),
    gs_cis_dir    = igloo_path("UKB", "ProtGScis"),
    gs_tr_dir     = igloo_path("UKB", "ProtGStrans"),
    omic_list     = heap_omicspred_protein_list,
    out_root      = heap_project_output("module1_predictive_r2_score_partition")
  ),

  known_factor_vars = c(
    "uk_biobank_assessment_centre_f54_0_0",
    "sex_f31_0_0"
  ),

  model_families = list(
    lasso = list(alpha_grid = 1),
    ridge = list(alpha_grid = 0),
    enet  = list(alpha_grid = c(0.1, 0.3, 0.5, 0.7, 0.9))
  )
)

set.seed(cfg$seed)

prot_clean <- function(protID) gsub("-", "_", protID)
`%||%` <- function(x, y) if (!is.null(x)) x else y

############################################################
# 1) PXS adapter
############################################################

as_pxs <- function(x) {
  if (is.list(x) && !isS4(x)) return(x)
  if (!isS4(x)) stop("PXS object must be S4 or list")

  list(
    Elist       = x@Elist,
    Elist_names = x@Elist_names,
    Eid_cat     = x@Eid_cat,
    ordinalIDs  = x@ordinalIDs,
    UKBprot_df  = x@UKBprot_df,
    protIDs     = x@protIDs,
    covars_df   = x@covars_df,
    covars_list = x@covars_list
  )
}

PXScovarSpec <- function(pxs, covars_subset) {
  avail <- intersect(covars_subset, names(pxs$covars_df))
  dropped <- setdiff(covars_subset, avail)

  if (length(dropped) > 0) {
    message(
      "PXScovarSpec: ", length(dropped),
      " requested covariate(s) not in covars_df and will be skipped: ",
      paste(head(dropped, 5), collapse = ", "),
      if (length(dropped) > 5) paste0(" ... (", length(dropped) - 5, " more)") else ""
    )
  }

  pxs$covars_list <- avail
  pxs$covars_df <- pxs$covars_df %>%
    dplyr::select(all_of(c("eid", avail)))

  pxs
}

############################################################
# 2) OmicsPred and genetic score loading
############################################################

make_omicspred_resolver <- function(map_path) {
  cache <- NULL

  function(protID) {
    if (is.null(cache)) cache <<- data.table::fread(map_path)

    opID <- cache$OMICSPRED_ID[match(protID, cache$Gene)]

    if (length(opID) == 0 || is.na(opID) || opID == "") {
      stop("No OMICSPRED_ID found for protID=", protID)
    }

    opID
  }
}

read_sscore_or_zero <- function(fpath, out_col, eids, strict = FALSE) {
  make_zero <- function(eids) {
    dt0 <- data.table::data.table(eid = as.integer(eids))
    dt0[, (out_col) := 0]
    dt0
  }

  if (!file.exists(fpath)) {
    msg <- paste0("Missing score file: ", fpath, " -> using 0 for ", out_col)
    if (strict) stop(msg) else warning(msg)
    return(make_zero(eids))
  }

  dt <- data.table::fread(fpath)

  id_col <- intersect(c("eid", "IID", "#IID", "id", "ID"), names(dt))[1]
  if (is.na(id_col)) id_col <- names(dt)[1]

  score_col <- intersect(
    c("SCORE1_AVG", "SCORE1_SUM", "SCORE1", "score", "SCORE"),
    names(dt)
  )[1]
  if (is.na(score_col)) score_col <- names(dt)[ncol(dt)]

  out <- dt[, .(
    eid = as.integer(get(id_col)),
    score = as.numeric(get(score_col))
  )]

  data.table::setnames(out, c("eid", out_col))
  out
}

extract_protGS <- function(protID,
                           eids,
                           op_resolve,
                           which = c("cis", "trans"),
                           strict = FALSE,
                           cfg = cfg) {
  which <- match.arg(which)

  dir <- if (which == "cis") cfg$paths$gs_cis_dir else cfg$paths$gs_tr_dir

  opID <- op_resolve(protID)
  fpath <- file.path(dir, paste0(opID, ".sscore"))

  p <- prot_clean(protID)
  out_col <- paste0(p, if (which == "cis") "_GScis" else "_GStrans")

  read_sscore_or_zero(
    fpath = fpath,
    out_col = out_col,
    eids = eids,
    strict = strict
  )
}

GS_struct <- function(protID,
                      UKBprot_df,
                      op_resolve,
                      strict = FALSE,
                      cfg = cfg) {
  p <- prot_clean(protID)

  omic_orig <- UKBprot_df %>%
    dplyr::select(all_of(c("eid", p)))

  eids <- omic_orig$eid

  cis <- extract_protGS(
    protID = protID,
    eids = eids,
    op_resolve = op_resolve,
    which = "cis",
    strict = strict,
    cfg = cfg
  )

  trans <- extract_protGS(
    protID = protID,
    eids = eids,
    op_resolve = op_resolve,
    which = "trans",
    strict = strict,
    cfg = cfg
  )

  omic_PGS <- merge(cis, trans, by = "eid", all = TRUE)

  list(
    combo = na.omit(merge(omic_PGS, omic_orig, by = "eid")),
    solo  = na.omit(omic_PGS)
  )
}

############################################################
# 3) Preprocessing helpers
############################################################

missing_cols <- function(df, miss_rate = 0.2) {
  na_rate <- colMeans(is.na(df))
  names(na_rate[na_rate > miss_rate])
}

drop_missing_cols <- function(df, miss_rate = 0.2) {
  drop <- missing_cols(df, miss_rate)
  dplyr::select(df, !all_of(drop))
}

continuous_finder <- function(df) {
  if (ncol(df) == 0) return(character(0))

  max_cols <- apply(df, 2, max, na.rm = TRUE)
  unique_vals <- vapply(
    df,
    function(x) length(unique(x[!is.na(x)])),
    numeric(1)
  )

  names(max_cols[max_cols > 5 | unique_vals > 2])
}

coerce_known_categoricals <- function(df, as_factor = cfg$known_factor_vars) {
  for (v in intersect(as_factor, names(df))) {
    df[[v]] <- as.factor(df[[v]])
  }

  df
}

categorical_handler <- function(df,
                                ordinal_names,
                                ordinal_contrast = "treatment") {
  for (v in ordinal_names) {
    if (!v %in% names(df)) next

    df[[v]] <- factor(df[[v]], ordered = FALSE)
    lv <- levels(df[[v]])

    if (length(lv) <= 1L) next

    if (ordinal_contrast == "treatment") {
      # Reference level = "0" so every coefficient is "level k vs 0". HEAP ordinals
      # are coded from 0 upward (verified: all 31 have min == 0); warn + fall back
      # to the lowest level if a "0" level is ever absent.
      base_idx <- match("0", lv)
      if (is.na(base_idx)) {
        warning("ordinal '", v, "' has no '0' level; using lowest level '", lv[1],
                "' as the treatment reference")
        base_idx <- 1L
      }
      contrasts(df[[v]]) <- contr.treatment(length(lv), base = base_idx)
    } else if (ordinal_contrast == "sum") {
      contrasts(df[[v]]) <- contr.sum(length(lv))
    }
  }

  df
}

get_factor_schema <- function(df) {
  facs <- names(df)[vapply(df, is.factor, logical(1))]
  levels_map <- setNames(lapply(facs, function(v) levels(df[[v]])), facs)
  list(levels = levels_map)
}

apply_factor_schema <- function(df, schema) {
  if (is.null(schema) || is.null(schema$levels)) return(df)

  for (v in intersect(names(schema$levels), names(df))) {
    df[[v]] <- factor(df[[v]], levels = schema$levels[[v]])
  }

  df
}

fit_scaler <- function(train_df, numeric_cols) {
  numeric_cols <- intersect(numeric_cols, names(train_df))

  mu <- lapply(train_df[numeric_cols], function(x) mean(x, na.rm = TRUE))
  sd <- lapply(train_df[numeric_cols], function(x) stats::sd(x, na.rm = TRUE))

  list(mean = mu, sd = sd, cols = numeric_cols)
}

apply_scaler <- function(df, scaler) {
  for (v in scaler$cols) {
    if (!v %in% names(df)) next

    s <- scaler$sd[[v]]

    if (!is.finite(s) || s == 0) {
      df[[v]] <- 0
    } else {
      df[[v]] <- (df[[v]] - scaler$mean[[v]]) / s
    }
  }

  df
}

############################################################
# 4) Dataset assembly
############################################################

assemble_df_for_prot <- function(protID, pxs, op_resolve, cfg = cfg) {
  pxs <- as_pxs(pxs)
  p <- prot_clean(protID)

  gs <- GS_struct(
    protID = protID,
    UKBprot_df = pxs$UKBprot_df,
    op_resolve = op_resolve,
    strict = FALSE,
    cfg = cfg
  )

  if (!is.null(pxs$base_ec)) {
    df <- merge(gs$combo, pxs$base_ec, by = "eid")
    E_ids_all <- pxs$E_ids_all
  } else {
    E_df <- pxs$Elist %>%
      purrr::reduce(full_join, by = "eid")

    df <- merge(gs$combo, E_df, by = "eid")
    df <- merge(df, pxs$covars_df, by = "eid")

    E_ids_all <- setdiff(names(E_df), "eid")
  }

  df <- as.data.frame(df)

  drop <- missing_cols(df, cfg$miss_rate)
  df <- drop_missing_cols(df, cfg$miss_rate)

  E_ids <- setdiff(E_ids_all, drop)
  covars_used <- setdiff(pxs$covars_list, drop)
  covars_used <- intersect(covars_used, names(df))

  df <- coerce_known_categoricals(df, cfg$known_factor_vars)
  df <- categorical_handler(df, pxs$ordinalIDs, ordinal_contrast = "treatment")

  schema <- get_factor_schema(df)

  model_cols <- unique(c(
    "eid",
    p,
    paste0(p, "_GScis"),
    paste0(p, "_GStrans"),
    E_ids,
    covars_used
  ))

  model_cols <- intersect(model_cols, names(df))

  before_n <- nrow(df)
  df <- df[stats::complete.cases(df[, model_cols, drop = FALSE]), , drop = FALSE]
  after_n <- nrow(df)

  message("assemble_df_for_prot(): kept ", after_n, " / ", before_n, " rows for ", protID)

  list(
    df = df,
    meta = list(
      E_ids = intersect(E_ids, names(df)),
      covars_used = covars_used,
      ordinalIDs = pxs$ordinalIDs,
      factor_schema = schema,
      Elist_names = pxs$Elist_names,
      Eid_cat = pxs$Eid_cat
    )
  )
}

CV_split <- function(protID,
                     pxs,
                     op_resolve,
                     kfold = cfg$kfold,
                     cfg = cfg) {
  assembled <- assemble_df_for_prot(protID, pxs, op_resolve, cfg)

  df <- assembled$df
  meta <- assembled$meta

  folds <- createFolds(seq_len(nrow(df)), k = kfold, list = TRUE)

  train_list <- vector("list", kfold)
  test_list <- vector("list", kfold)

  for (i in seq_len(kfold)) {
    val_idx <- folds[[i]]
    train_list[[i]] <- df[-val_idx, , drop = FALSE]
    test_list[[i]] <- df[val_idx, , drop = FALSE]
  }

  fold_assign <- purrr::map2_dfr(
    test_list,
    seq_along(test_list),
    ~ tibble::tibble(eid = .x$eid, fold = .y)
  )

  list(
    train = train_list,
    test = test_list,
    fold_assign = fold_assign,
    meta = meta
  )
}

preprocess_fold <- function(protID,
                            train_df,
                            test_df,
                            covars_used,
                            ordinal_names) {
  p <- prot_clean(protID)
  G_cis <- paste0(p, "_GScis")
  G_tr <- paste0(p, "_GStrans")

  train_df <- as.data.frame(train_df)
  test_df <- as.data.frame(test_df)

  rel_cols <- c(p, G_cis, G_tr, covars_used)
  rel_cols <- intersect(rel_cols, names(train_df))

  numeric_cols <- vapply(
    train_df[, rel_cols, drop = FALSE],
    function(col) is.numeric(col) && !(all(col %in% c(0, 1))),
    logical(1)
  )

  DG_cont <- names(numeric_cols[numeric_cols])

  drop_cols <- unique(c("eid", rel_cols, ordinal_names))
  df_E <- train_df %>%
    dplyr::select(!all_of(intersect(drop_cols, names(train_df))))

  # Continuous EXPOSURES from the declared variable_type (analysis_exposures.tsv);
  # fall back to the value-range heuristic only if the config is unavailable. Like
  # ordinalIDs (now routed by declared type in as_pxs_baseline), this stops the
  # small-scale IMD scores being mis-scaled/mis-encoded.
  .E_cont_cfg <- heap_exposures_of_type("continuous", names(df_E))
  E_cont <- if (is.null(.E_cont_cfg)) continuous_finder(df_E) else .E_cont_cfg
  num_cols <- unique(c(DG_cont, E_cont))

  scaler <- fit_scaler(train_df, num_cols)

  list(
    train = apply_scaler(train_df, scaler),
    test = apply_scaler(test_df, scaler),
    scaler = scaler
  )
}

############################################################
# 5) Block construction
############################################################

cross_terms <- function(a, b) {
  if (length(a) == 0 || length(b) == 0) return(character(0))
  as.vector(outer(a, b, paste, sep = ":"))
}

sanitize_block_name <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

make_exposure_category_lookup <- function(Eid_cat) {
  out <- Eid_cat$Category
  names(out) <- Eid_cat$Eid
  out
}

get_exposure_categories <- function(E_ids, Eid_cat) {
  lut <- make_exposure_category_lookup(Eid_cat)
  cats <- unique(unname(lut[E_ids]))
  cats <- cats[!is.na(cats)]
  sort(unique(cats))
}

get_E_ids_for_category <- function(E_ids, Eid_cat, category) {
  lut <- make_exposure_category_lookup(Eid_cat)
  E_ids[lut[E_ids] == category]
}

build_formula_from_terms <- function(protID, rhs_terms) {
  p <- prot_clean(protID)

  rhs_terms <- unique(rhs_terms)
  rhs_terms <- rhs_terms[nzchar(rhs_terms)]

  if (length(rhs_terms) == 0) {
    as.formula(paste(p, "~ 1"))
  } else {
    as.formula(paste(p, "~", paste(rhs_terms, collapse = " + ")))
  }
}

make_terms <- function(protID,
                       covars,
                       E_ids,
                       include_c = TRUE,
                       include_gcis = FALSE,
                       include_gtrans = FALSE,
                       include_e_ids = character(0),
                       include_gxe_ids = character(0)) {
  p <- prot_clean(protID)
  Gcis <- paste0(p, "_GScis")
  Gtrans <- paste0(p, "_GStrans")

  terms <- character(0)

  if (include_c) {
    terms <- c(terms, covars)
  }

  if (include_gcis) {
    terms <- c(terms, Gcis)
  }

  if (include_gtrans) {
    terms <- c(terms, Gtrans)
  }

  if (length(include_e_ids) > 0) {
    terms <- c(terms, include_e_ids)
  }

  if (length(include_gxe_ids) > 0) {
    if (include_gcis) {
      terms <- c(terms, cross_terms(Gcis, include_gxe_ids))
    }
    if (include_gtrans) {
      terms <- c(terms, cross_terms(Gtrans, include_gxe_ids))
    }
  }

  unique(terms)
}

build_model_terms_by_name <- function(protID, meta, model_name) {
  E_all <- meta$E_ids
  C <- meta$covars_used

  # Include a genetic component unless it was explicitly flagged absent for this
  # protein (no cis/trans .sscore -> all-zero column dropped as zero-variance).
  # Default to present when the flag is unset, so callers without the flag keep
  # the original behaviour. Excluding an absent term keeps the dropped column out
  # of the model formula (the source of the `object '<p>_GScis' not found` crash).
  gcis_ok   <- !identical(meta$gcis_present, FALSE)
  gtrans_ok <- !identical(meta$gtrans_present, FALSE)

  if (model_name == "C") {
    return(make_terms(protID, C, E_all, include_c = TRUE))
  }

  if (model_name == "C_G") {
    return(make_terms(
      protID, C, E_all,
      include_c = TRUE,
      include_gcis = gcis_ok,
      include_gtrans = gtrans_ok
    ))
  }

  if (model_name == "C_G_E") {
    return(make_terms(
      protID, C, E_all,
      include_c = TRUE,
      include_gcis = gcis_ok,
      include_gtrans = gtrans_ok,
      include_e_ids = E_all
    ))
  }

  if (model_name == "C_G_E_GxE") {
    return(make_terms(
      protID, C, E_all,
      include_c = TRUE,
      include_gcis = gcis_ok,
      include_gtrans = gtrans_ok,
      include_e_ids = E_all,
      include_gxe_ids = E_all
    ))
  }

  stop("Unknown model_name: ", model_name)
}

extract_contrasts <- function(df) {
  out <- list()

  for (nm in names(df)) {
    if (is.factor(df[[nm]])) {
      out[[nm]] <- contrasts(df[[nm]], contrasts = TRUE)
    }
  }

  out
}

build_mats_with_mapping <- function(formula, train_df, test_df) {
  contr <- extract_contrasts(train_df)

  X_train <- Matrix::sparse.model.matrix(
    formula,
    data = train_df,
    contrasts.arg = contr
  )

  X_test <- Matrix::sparse.model.matrix(
    formula,
    data = test_df,
    contrasts.arg = contr
  )

  term_labels <- attr(terms(formula), "term.labels")
  col_assign <- attr(X_train, "assign")

  dm_to_term <- setNames(
    term_labels[col_assign[col_assign > 0L]],
    colnames(X_train)[col_assign > 0L]
  )

  term_to_dm <- split(names(dm_to_term), dm_to_term)

  list(
    X_train = X_train,
    X_test = X_test,
    dm_to_term = dm_to_term,
    term_to_dm = term_to_dm
  )
}

############################################################
# 6) R2 and penalized model fitting
############################################################

r2_score_trainmean <- function(y_test, yhat_test, y_train_mean) {
  ss_res <- sum((y_test - yhat_test)^2)
  ss_tot <- sum((y_test - y_train_mean)^2)

  if (!is.finite(ss_tot) || ss_tot == 0) return(NA_real_)

  1 - ss_res / ss_tot
}

r2_score_internal <- function(y, yhat) {
  ss_res <- sum((y - yhat)^2)
  ss_tot <- sum((y - mean(y))^2)

  if (!is.finite(ss_tot) || ss_tot == 0) return(NA_real_)

  1 - ss_res / ss_tot
}

safe_delta <- function(a, b) {
  out <- a - b
  if (!is.finite(out)) NA_real_ else out
}

glmnet_fit_predict <- function(y_train,
                               y_test,
                               X_train,
                               X_test,
                               family_type = c("lasso", "ridge", "enet"),
                               seed = cfg$seed,
                               cfg = cfg) {
  family_type <- match.arg(family_type)
  alpha_grid <- cfg$model_families[[family_type]]$alpha_grid

  inner_k <- cfg$glmnet_inner_kfold %||% 5L
  inner_k <- as.integer(inner_k)

  if (!is.finite(inner_k) || inner_k < 2L) {
    inner_k <- 5L
  }

  inner_k <- min(inner_k, length(y_train))

  set.seed(seed)
  foldid <- sample(rep(seq_len(inner_k), length.out = length(y_train)))

  cv_res <- lapply(alpha_grid, function(a) {
    t0 <- Sys.time()

    cvfit <- cv.glmnet(
      x = X_train,
      y = y_train,
      alpha = a,
      standardize = FALSE,
      intercept = TRUE,
      foldid = foldid
    )

    message(
      "cv.glmnet done: family=", family_type,
      " alpha=", a,
      " n=", length(y_train),
      " p=", ncol(X_train),
      " time_min=", round(difftime(Sys.time(), t0, units = "mins"), 3)
    )

    list(
      alpha = a,
      lambda = cvfit$lambda.min,
      cvm = min(cvfit$cvm),
      cvfit = cvfit
    )
  })

  best_idx <- which.min(vapply(cv_res, `[[`, numeric(1), "cvm"))
  best <- cv_res[[best_idx]]

  fit <- glmnet(
    x = X_train,
    y = y_train,
    alpha = best$alpha,
    lambda = best$lambda,
    standardize = FALSE,
    intercept = TRUE
  )

  yhat_train <- as.numeric(predict(fit, newx = X_train, s = best$lambda))
  yhat_test <- as.numeric(predict(fit, newx = X_test, s = best$lambda))

  b <- coef(fit, s = best$lambda)
  b_vec <- as.numeric(b)
  names(b_vec) <- rownames(b)

  nonzero_names <- names(b_vec)[b_vec != 0]
  n_nonzero <- sum(b_vec != 0) - as.integer("(Intercept)" %in% nonzero_names)

  list(
    pred_train = yhat_train,
    pred_test = yhat_test,
    alpha = best$alpha,
    lambda = best$lambda,
    cvm = best$cvm,
    n_nonzero = n_nonzero,
    coefs = b_vec[names(b_vec) != "(Intercept)"],
    intercept = b_vec["(Intercept)"] %||% 0
  )
}

############################################################
# 7) Component score construction
############################################################

compute_lasso_component_scores <- function(protID,
                                           meta,
                                           X_train,
                                           X_test,
                                           coefs,
                                           term_to_dm,
                                           train_df,
                                           test_df,
                                           sensitivity_terms = list()) {
  p <- prot_clean(protID)

  E_ids <- meta$E_ids
  C_ids <- meta$covars_used
  Eid_cat <- meta$Eid_cat

  Gcis_term <- paste0(p, "_GScis")
  Gtrans_term <- paste0(p, "_GStrans")

  coef_names <- names(coefs)
  nz_coef_names <- names(coefs)[is.finite(coefs) & coefs != 0]

  get_dm_for_terms <- function(terms) {
    terms <- intersect(terms, names(term_to_dm))
    unlist(term_to_dm[terms], use.names = FALSE)
  }

  score_from_terms <- function(X, terms) {
    dm_cols <- get_dm_for_terms(terms)
    dm_cols <- intersect(dm_cols, colnames(X))
    dm_cols <- intersect(dm_cols, coef_names)

    if (length(dm_cols) == 0) {
      return(rep(0, nrow(X)))
    }

    as.numeric(as.matrix(X[, dm_cols, drop = FALSE]) %*% coefs[dm_cols])
  }

  C_train <- score_from_terms(X_train, C_ids)
  C_test <- score_from_terms(X_test, C_ids)

  Gcis_train <- score_from_terms(X_train, Gcis_term)
  Gcis_test <- score_from_terms(X_test, Gcis_term)

  Gtrans_train <- score_from_terms(X_train, Gtrans_term)
  Gtrans_test <- score_from_terms(X_test, Gtrans_term)

  E_train <- score_from_terms(X_train, E_ids)
  E_test <- score_from_terms(X_test, E_ids)

  gxe_terms <- c(
    cross_terms(Gcis_term, E_ids),
    cross_terms(Gtrans_term, E_ids)
  )

  GIS_train <- score_from_terms(X_train, gxe_terms)
  GIS_test <- score_from_terms(X_test, gxe_terms)

  Gcis_raw_train <- if (Gcis_term %in% names(train_df)) train_df[[Gcis_term]] else rep(0, nrow(train_df))
  Gcis_raw_test  <- if (Gcis_term %in% names(test_df))  test_df[[Gcis_term]]  else rep(0, nrow(test_df))

  Gtrans_raw_train <- if (Gtrans_term %in% names(train_df)) train_df[[Gtrans_term]] else rep(0, nrow(train_df))
  Gtrans_raw_test  <- if (Gtrans_term %in% names(test_df))  test_df[[Gtrans_term]]  else rep(0, nrow(test_df))

  # Structural availability of the cis/trans genetic component for this protein
  # (FALSE = no .sscore / no variants -> the *_raw and *_score columns below are a
  # true 0, not an estimated ~0). Carried into mediation scores so Module 3 can
  # skip the degenerate instrument and label it, rather than fit on an all-zero PGS.
  gcis_present_flag   <- !identical(meta$gcis_present, FALSE)
  gtrans_present_flag <- !identical(meta$gtrans_present, FALSE)

  train_scores <- tibble::tibble(
    protein_value = train_df[[p]],
    C_score = C_train,
    Gcis_score = Gcis_train,
    Gtrans_score = Gtrans_train,
    G_score = Gcis_train + Gtrans_train,
    PXS_total = E_train,
    GIS_total = GIS_train,
    Gcis_raw = Gcis_raw_train,
    Gtrans_raw = Gtrans_raw_train,
    G_raw = Gcis_raw_train + Gtrans_raw_train,
    gcis_present = gcis_present_flag,
    gtrans_present = gtrans_present_flag
  )

  test_scores <- tibble::tibble(
    eid = test_df$eid,
    protein_value = test_df[[p]],
    C_score = C_test,
    Gcis_score = Gcis_test,
    Gtrans_score = Gtrans_test,
    G_score = Gcis_test + Gtrans_test,
    PXS_total = E_test,
    GIS_total = GIS_test,
    Gcis_raw = Gcis_raw_test,
    Gtrans_raw = Gtrans_raw_test,
    G_raw = Gcis_raw_test + Gtrans_raw_test,
    gcis_present = gcis_present_flag,
    gtrans_present = gtrans_present_flag
  )

  cats <- get_exposure_categories(E_ids, Eid_cat)

  audit_rows <- list()

  for (cat in cats) {
    cat_safe <- sanitize_block_name(cat)
    cat_E_ids <- get_E_ids_for_category(E_ids, Eid_cat, cat)

    pxs_col <- paste0("PXS_", cat_safe)
    gis_col <- paste0("GIS_", cat_safe)

    train_scores[[pxs_col]] <- score_from_terms(X_train, cat_E_ids)
    test_scores[[pxs_col]] <- score_from_terms(X_test, cat_E_ids)

    cat_gxe_terms <- c(
      cross_terms(Gcis_term, cat_E_ids),
      cross_terms(Gtrans_term, cat_E_ids)
    )

    train_scores[[gis_col]] <- score_from_terms(X_train, cat_gxe_terms)
    test_scores[[gis_col]] <- score_from_terms(X_test, cat_gxe_terms)

    for (eid_var in cat_E_ids) {
      e_dm <- get_dm_for_terms(eid_var)
      gxe_dm <- get_dm_for_terms(c(
        cross_terms(Gcis_term, eid_var),
        cross_terms(Gtrans_term, eid_var)
      ))

      audit_rows[[length(audit_rows) + 1L]] <- tibble::tibble(
        exposure_variable = eid_var,
        category = cat,
        n_E_design_cols = length(e_dm),
        n_GxE_design_cols = length(gxe_dm),
        E_design_cols = paste(e_dm, collapse = "|"),
        GxE_design_cols = paste(gxe_dm, collapse = "|"),
        E_nonzero = any(intersect(e_dm, nz_coef_names) %in% nz_coef_names),
        GxE_nonzero = any(intersect(gxe_dm, nz_coef_names) %in% nz_coef_names)
      )
    }
  }

  # Optional sensitivity scores (GxC, ExC, or any named term set)
  for (sens_name in names(sensitivity_terms)) {
    col_name <- paste0(sens_name, "_score")
    train_scores[[col_name]] <- score_from_terms(X_train, sensitivity_terms[[sens_name]])
    test_scores[[col_name]]  <- score_from_terms(X_test,  sensitivity_terms[[sens_name]])
  }

  audit_tbl <- if (length(audit_rows) > 0) {
    dplyr::bind_rows(audit_rows)
  } else {
    tibble::tibble()
  }

  mediation_scores <- test_scores %>%
    dplyr::mutate(
      protID = p,
      .before = 1L
    )

  list(
    train_scores = train_scores,
    test_scores = test_scores,
    mediation_scores = mediation_scores,
    audit = audit_tbl
  )
}

############################################################
# 8) Score-level unique R2
############################################################

fit_score_lm_predict <- function(y_train,
                                 y_test,
                                 train_scores,
                                 test_scores,
                                 score_cols) {
  score_cols <- intersect(score_cols, names(train_scores))
  score_cols <- intersect(score_cols, names(test_scores))

  if (length(score_cols) == 0) {
    pred_train <- rep(mean(y_train), length(y_train))
    pred_test <- rep(mean(y_train), length(y_test))

    return(list(
      pred_train = pred_train,
      pred_test = pred_test,
      train_r2 = r2_score_internal(y_train, pred_train),
      test_r2 = r2_score_trainmean(y_test, pred_test, mean(y_train)),
      n_score_cols = 0L
    ))
  }

  train_model_df <- data.frame(
    y = y_train,
    train_scores[, score_cols, drop = FALSE],
    check.names = FALSE
  )

  test_model_df <- data.frame(
    test_scores[, score_cols, drop = FALSE],
    check.names = FALSE
  )

  # Remove zero-variance score columns in training.
  keep <- vapply(
    train_model_df[, score_cols, drop = FALSE],
    function(x) {
      sx <- stats::sd(x, na.rm = TRUE)
      is.finite(sx) && sx > 0
    },
    logical(1)
  )

  score_cols <- score_cols[keep]

  if (length(score_cols) == 0) {
    pred_train <- rep(mean(y_train), length(y_train))
    pred_test <- rep(mean(y_train), length(y_test))

    return(list(
      pred_train = pred_train,
      pred_test = pred_test,
      train_r2 = r2_score_internal(y_train, pred_train),
      test_r2 = r2_score_trainmean(y_test, pred_test, mean(y_train)),
      n_score_cols = 0L
    ))
  }

  fml <- as.formula(paste("y ~", paste(score_cols, collapse = " + ")))

  fit <- stats::lm(fml, data = train_model_df)

  pred_train <- as.numeric(stats::predict(fit, newdata = train_model_df))
  pred_test <- as.numeric(stats::predict(fit, newdata = test_model_df))

  list(
    pred_train = pred_train,
    pred_test = pred_test,
    train_r2 = r2_score_internal(y_train, pred_train),
    test_r2 = r2_score_trainmean(y_test, pred_test, mean(y_train)),
    n_score_cols = length(score_cols)
  )
}

unique_r2_drop <- function(y_train,
                           y_test,
                           train_scores,
                           test_scores,
                           full_cols,
                           drop_cols) {
  full_cols <- intersect(full_cols, names(train_scores))
  full_cols <- intersect(full_cols, names(test_scores))
  drop_cols <- intersect(drop_cols, full_cols)

  reduced_cols <- setdiff(full_cols, drop_cols)

  full_fit <- fit_score_lm_predict(
    y_train = y_train,
    y_test = y_test,
    train_scores = train_scores,
    test_scores = test_scores,
    score_cols = full_cols
  )

  reduced_fit <- fit_score_lm_predict(
    y_train = y_train,
    y_test = y_test,
    train_scores = train_scores,
    test_scores = test_scores,
    score_cols = reduced_cols
  )

  tibble::tibble(
    r2 = safe_delta(full_fit$test_r2, reduced_fit$test_r2),
    r2_full = full_fit$test_r2,
    r2_reduced = reduced_fit$test_r2,
    # train-side counterparts (same drop, evaluated on the training fold) so a
    # component's generalisation can be assessed per block, not just whole-model
    r2_train = safe_delta(full_fit$train_r2, reduced_fit$train_r2),
    r2_full_train = full_fit$train_r2,
    r2_reduced_train = reduced_fit$train_r2,
    n_full_score_cols = full_fit$n_score_cols,
    n_reduced_score_cols = reduced_fit$n_score_cols
  )
}

compute_score_level_unique_r2 <- function(protID,
                                          fold,
                                          family_type,
                                          y_train,
                                          y_test,
                                          train_scores,
                                          test_scores,
                                          meta,
                                          cfg = cfg) {
  p <- prot_clean(protID)

  # Base primary columns; extend with any sensitivity scores that were computed
  primary_cols <- c("C_score", "G_score", "PXS_total", "GIS_total")
  sensitivity_score_cols <- grep("^(GxC|ExC)_score$", names(train_scores), value = TRUE)
  all_cols <- c(primary_cols, sensitivity_score_cols)

  fit_C <- fit_score_lm_predict(
    y_train = y_train,
    y_test = y_test,
    train_scores = train_scores,
    test_scores = test_scores,
    score_cols = "C_score"
  )

  fit_CG <- fit_score_lm_predict(
    y_train = y_train,
    y_test = y_test,
    train_scores = train_scores,
    test_scores = test_scores,
    score_cols = c("C_score", "G_score")
  )

  fit_CGE <- fit_score_lm_predict(
    y_train = y_train,
    y_test = y_test,
    train_scores = train_scores,
    test_scores = test_scores,
    score_cols = c("C_score", "G_score", "PXS_total")
  )

  fit_full <- fit_score_lm_predict(
    y_train = y_train,
    y_test = y_test,
    train_scores = train_scores,
    test_scores = test_scores,
    score_cols = all_cols
  )

  unique_C <- unique_r2_drop(
    y_train, y_test, train_scores, test_scores,
    full_cols = all_cols, drop_cols = "C_score"
  )

  unique_G <- unique_r2_drop(
    y_train, y_test, train_scores, test_scores,
    full_cols = all_cols, drop_cols = "G_score"
  )

  unique_E <- unique_r2_drop(
    y_train, y_test, train_scores, test_scores,
    full_cols = all_cols, drop_cols = "PXS_total"
  )

  unique_GxE <- unique_r2_drop(
    y_train, y_test, train_scores, test_scores,
    full_cols = all_cols, drop_cols = "GIS_total"
  )

  r2_coarse <- dplyr::bind_rows(
    tibble::tibble(
      omic = p,
      fold = fold,
      family = family_type,
      method = "score_model_total",
      block = c("C", "C+G", "C+G+E", "C+G+E+GxE"),
      r2 = c(
        fit_C$test_r2,
        fit_CG$test_r2,
        fit_CGE$test_r2,
        fit_full$test_r2
      ),
      r2_train = c(
        fit_C$train_r2,
        fit_CG$train_r2,
        fit_CGE$train_r2,
        fit_full$train_r2
      ),
      base_model = NA_character_,
      full_model = c("C", "C+G", "C+G+E", "C+G+E+GxE"),
      r2_full = NA_real_,
      r2_reduced = NA_real_,
      r2_full_train = NA_real_,
      r2_reduced_train = NA_real_
    ),
    tibble::tibble(
      omic = p,
      fold = fold,
      family = family_type,
      method = "score_unique_drop",
      block = c("Covars", "G", "E", "GxE"),
      r2 = c(
        unique_C$r2,
        unique_G$r2,
        unique_E$r2,
        unique_GxE$r2
      ),
      r2_train = c(
        unique_C$r2_train,
        unique_G$r2_train,
        unique_E$r2_train,
        unique_GxE$r2_train
      ),
      base_model = c(
        "G+E+GxE",
        "C+E+GxE",
        "C+G+GxE",
        "C+G+E"
      ),
      full_model = "C+G+E+GxE",
      r2_full = c(
        unique_C$r2_full,
        unique_G$r2_full,
        unique_E$r2_full,
        unique_GxE$r2_full
      ),
      r2_reduced = c(
        unique_C$r2_reduced,
        unique_G$r2_reduced,
        unique_E$r2_reduced,
        unique_GxE$r2_reduced
      ),
      r2_full_train = c(
        unique_C$r2_full_train,
        unique_G$r2_full_train,
        unique_E$r2_full_train,
        unique_GxE$r2_full_train
      ),
      r2_reduced_train = c(
        unique_C$r2_reduced_train,
        unique_G$r2_reduced_train,
        unique_E$r2_reduced_train,
        unique_GxE$r2_reduced_train
      )
    )
  )

  r2_genetic <- tibble::tibble()

  if (isTRUE(cfg$run_genetic_subblocks)) {
    genetic_cols <- c("C_score", "Gcis_score", "Gtrans_score", "PXS_total", "GIS_total")

    unique_Gcis <- unique_r2_drop(
      y_train, y_test, train_scores, test_scores,
      full_cols = genetic_cols,
      drop_cols = "Gcis_score"
    )

    unique_Gtrans <- unique_r2_drop(
      y_train, y_test, train_scores, test_scores,
      full_cols = genetic_cols,
      drop_cols = "Gtrans_score"
    )

    r2_genetic <- tibble::tibble(
      omic = p,
      fold = fold,
      family = family_type,
      method = "score_unique_drop",
      block = c("Gcis", "Gtrans"),
      r2 = c(unique_Gcis$r2, unique_Gtrans$r2),
      r2_train = c(unique_Gcis$r2_train, unique_Gtrans$r2_train),
      # FALSE = component structurally absent (no cis/trans variants): r2 is a
      # true 0, distinguishable from an estimated ~0.
      present = c(!identical(meta$gcis_present, FALSE),
                 !identical(meta$gtrans_present, FALSE)),
      base_model = c(
        "C+Gtrans+E+GxE",
        "C+Gcis+E+GxE"
      ),
      full_model = "C+Gcis+Gtrans+E+GxE",
      r2_full = c(unique_Gcis$r2_full, unique_Gtrans$r2_full),
      r2_reduced = c(unique_Gcis$r2_reduced, unique_Gtrans$r2_reduced),
      r2_full_train = c(unique_Gcis$r2_full_train, unique_Gtrans$r2_full_train),
      r2_reduced_train = c(unique_Gcis$r2_reduced_train, unique_Gtrans$r2_reduced_train)
    )
  }

  r2_ecat <- tibble::tibble()

  if (isTRUE(cfg$run_exposure_categories)) {
    pxs_cat_cols <- grep("^PXS_", names(train_scores), value = TRUE)
    pxs_cat_cols <- setdiff(pxs_cat_cols, "PXS_total")

    ecat_full_cols <- c("C_score", "G_score", pxs_cat_cols, "GIS_total")

    if (length(pxs_cat_cols) > 0) {
      for (col in pxs_cat_cols) {
        u <- unique_r2_drop(
          y_train, y_test,
          train_scores, test_scores,
          full_cols = ecat_full_cols,
          drop_cols = col
        )

        cat <- sub("^PXS_", "", col)

        r2_ecat <- dplyr::bind_rows(
          r2_ecat,
          tibble::tibble(
            omic = p,
            fold = fold,
            family = family_type,
            method = "score_unique_drop",
            category = cat,
            block = col,
            r2 = u$r2,
            r2_train = u$r2_train,
            base_model = paste0("all_scores_without_", col),
            full_model = "C+G+PXS_categories+GIS",
            r2_full = u$r2_full,
            r2_reduced = u$r2_reduced,
            r2_full_train = u$r2_full_train,
            r2_reduced_train = u$r2_reduced_train
          )
        )
      }
    }
  }

  r2_gxecat <- tibble::tibble()

  if (isTRUE(cfg$run_gxe_categories)) {
    gis_cat_cols <- grep("^GIS_", names(train_scores), value = TRUE)
    gis_cat_cols <- setdiff(gis_cat_cols, "GIS_total")

    gxe_full_cols <- c("C_score", "G_score", "PXS_total", gis_cat_cols)

    if (length(gis_cat_cols) > 0) {
      for (col in gis_cat_cols) {
        u <- unique_r2_drop(
          y_train, y_test,
          train_scores, test_scores,
          full_cols = gxe_full_cols,
          drop_cols = col
        )

        cat <- sub("^GIS_", "", col)

        r2_gxecat <- dplyr::bind_rows(
          r2_gxecat,
          tibble::tibble(
            omic = p,
            fold = fold,
            family = family_type,
            method = "score_unique_drop",
            category = cat,
            block = col,
            r2 = u$r2,
            r2_train = u$r2_train,
            base_model = paste0("all_scores_without_", col),
            full_model = "C+G+PXS+GIS_categories",
            r2_full = u$r2_full,
            r2_reduced = u$r2_reduced,
            r2_full_train = u$r2_full_train,
            r2_reduced_train = u$r2_reduced_train
          )
        )
      }
    }
  }

  # Sensitivity score unique drops (GxC_score, ExC_score when present)
  r2_sensitivity <- tibble::tibble()

  for (sens_col in sensitivity_score_cols) {
    u <- unique_r2_drop(
      y_train, y_test, train_scores, test_scores,
      full_cols = all_cols, drop_cols = sens_col
    )
    r2_sensitivity <- dplyr::bind_rows(
      r2_sensitivity,
      tibble::tibble(
        omic       = p,
        fold       = fold,
        family     = family_type,
        method     = "score_unique_drop",
        block      = sens_col,
        r2         = u$r2,
        r2_train   = u$r2_train,
        base_model = paste0("all_scores_without_", sens_col),
        full_model = paste(all_cols, collapse = "+"),
        r2_full    = u$r2_full,
        r2_reduced = u$r2_reduced,
        r2_full_train = u$r2_full_train,
        r2_reduced_train = u$r2_reduced_train
      )
    )
  }

  list(
    R2Coarse = r2_coarse,
    R2GeneticSubblocks = r2_genetic,
    R2ExposureCategories = r2_ecat,
    R2GxECategories = r2_gxecat,
    R2Sensitivity = r2_sensitivity
  )
}

############################################################
# 9) Fold-level score-partition runner
############################################################

run_fold_score_partition <- function(protID,
                                     fold_obj,
                                     idx,
                                     family_type = c("lasso", "ridge", "enet"),
                                     cfg = cfg) {
  family_type <- match.arg(family_type)

  p <- prot_clean(protID)
  meta <- fold_obj$meta

  train <- fold_obj$train[[idx]]
  test <- fold_obj$test[[idx]]

  train <- apply_factor_schema(train, meta$factor_schema)
  test <- apply_factor_schema(test, meta$factor_schema)

  pp <- preprocess_fold(
    protID = protID,
    train_df = train,
    test_df = test,
    covars_used = meta$covars_used,
    ordinal_names = meta$ordinalIDs
  )

  train <- pp$train
  test <- pp$test

  # Drop ZERO-VARIANCE features from this fold (never the response p): singleton
  # factors (< 2 non-empty levels) AND constant numeric columns (globally constant
  # exposures such as former_alcohol_drinker_f3731, or features that become
  # constant within a fold). Constant predictors contribute nothing and produce
  # degenerate design-matrix columns.
  zero_var_cols <- setdiff(heap_zero_variance_cols(train, names(train)), p)

  # Per-protein genetic-component availability for this fold. A protein with no
  # cis (or trans) variants has no .sscore file, so read_sscore_or_zero fills an
  # all-zero column that surfaces here as zero-variance. Recording presence lets
  # build_model_terms_by_name / sensitivity terms exclude the absent component
  # (a true 0-R2 instead of a formula crash on a dropped column) and flags it for
  # the R2 outputs and downstream mediation.
  Gcis_col   <- paste0(p, "_GScis")
  Gtrans_col <- paste0(p, "_GStrans")
  meta$gcis_present   <- (Gcis_col   %in% names(train)) && !(Gcis_col   %in% zero_var_cols)
  meta$gtrans_present <- (Gtrans_col %in% names(train)) && !(Gtrans_col %in% zero_var_cols)

  if (length(zero_var_cols) > 0L) {
    message(
      "fold ", idx,
      ": dropping ", length(zero_var_cols),
      " zero-variance feature(s): ",
      paste(zero_var_cols, collapse = ", ")
    )

    train <- train[, setdiff(names(train), zero_var_cols), drop = FALSE]
    test <- test[, setdiff(names(test), zero_var_cols), drop = FALSE]

    meta$covars_used <- setdiff(meta$covars_used, zero_var_cols)
    meta$E_ids <- setdiff(meta$E_ids, zero_var_cols)
    meta$ordinalIDs <- setdiff(meta$ordinalIDs, zero_var_cols)
  }

  y_train <- train[[p]]
  y_test <- test[[p]]

  # Build optional sensitivity term sets (GxC and/or ExC)
  sensitivity_terms <- list()

  if (isTRUE(cfg$run_gxc_sensitivity)) {
    gxc_avail <- intersect(cfg$gxc_covars, meta$covars_used)
    if (length(gxc_avail) > 0) {
      gxc_t <- c(
        if (!identical(meta$gcis_present, FALSE))   cross_terms(paste0(p, "_GScis"),   gxc_avail),
        if (!identical(meta$gtrans_present, FALSE)) cross_terms(paste0(p, "_GStrans"), gxc_avail)
      )
      if (length(gxc_t) > 0) {
        sensitivity_terms$GxC <- gxc_t
        message("fold ", idx, ": GxC sensitivity — ", length(gxc_t), " terms")
      }
    }
  }

  if (isTRUE(cfg$run_exc_sensitivity)) {
    exc_avail <- intersect(cfg$exc_covars, meta$covars_used)
    if (length(exc_avail) > 0) {
      exc_t <- cross_terms(meta$E_ids, exc_avail)
      if (length(exc_t) > 0) {
        sensitivity_terms$ExC <- exc_t
        message("fold ", idx, ": ExC sensitivity — ", length(exc_t), " terms")
      }
    }
  }

  full_terms <- unique(c(
    build_model_terms_by_name(protID, meta, "C_G_E_GxE"),
    unlist(sensitivity_terms, use.names = FALSE)
  ))

  fml <- build_formula_from_terms(protID, full_terms)

  message(
    "Fitting full model once: protein=", protID,
    " family=", family_type,
    " fold=", idx,
    " n_train=", length(y_train),
    " n_terms=", length(full_terms)
  )

  mats <- build_mats_with_mapping(fml, train, test)

  fitpred <- glmnet_fit_predict(
    y_train = y_train,
    y_test = y_test,
    X_train = mats$X_train,
    X_test = mats$X_test,
    family_type = family_type,
    seed = cfg$seed + idx,
    cfg = cfg
  )

  score_obj <- compute_lasso_component_scores(
    protID          = protID,
    meta            = meta,
    X_train         = mats$X_train,
    X_test          = mats$X_test,
    coefs           = fitpred$coefs,
    term_to_dm      = mats$term_to_dm,
    train_df        = train,
    test_df         = test,
    sensitivity_terms = sensitivity_terms
  )

  r2_obj <- compute_score_level_unique_r2(
    protID = protID,
    fold = idx,
    family_type = family_type,
    y_train = y_train,
    y_test = y_test,
    train_scores = score_obj$train_scores,
    test_scores = score_obj$test_scores,
    meta = meta,
    cfg = cfg
  )

  fit_summary <- tibble::tibble(
    omic = p,
    fold = idx,
    family = family_type,
    model_label = "C_G_E_GxE_full_lasso_once",
    model_class = "score_partition_full_fit",
    train_r2 = r2_score_internal(y_train, fitpred$pred_train),
    test_r2 = r2_score_trainmean(y_test, fitpred$pred_test, mean(y_train)),
    alpha = fitpred$alpha %||% NA_real_,
    lambda = fitpred$lambda %||% NA_real_,
    cvm = fitpred$cvm %||% NA_real_,
    n_nonzero = fitpred$n_nonzero %||% NA_integer_,
    n_design_cols = ncol(mats$X_train),
    n_train = length(y_train),
    n_test = length(y_test),
    inner_kfold = cfg$glmnet_inner_kfold
  )

  oof_pred <- tibble::tibble(
    eid = test$eid,
    protID = p,
    fold = idx,
    family = family_type,
    model_label = "C_G_E_GxE_full_lasso_once",
    obs = y_test,
    pred = fitpred$pred_test
  )

  scale_tbl <- tibble::tibble(
    protID = p,
    fold = idx,
    var = pp$scaler$cols,
    mean = unlist(pp$scaler$mean),
    sd = unlist(pp$scaler$sd)
  )

  score_train <- score_obj$train_scores %>%
    dplyr::mutate(
      protID = p,
      fold = idx,
      .before = 1L
    )

  score_test <- score_obj$test_scores %>%
    dplyr::mutate(
      protID = p,
      fold = idx,
      .before = 1L
    )

  mediation_audit <- if (nrow(score_obj$audit) > 0) {
    score_obj$audit %>%
      dplyr::mutate(
        protID = p,
        fold = idx,
        .before = 1L
      )
  } else {
    tibble::tibble()
  }

  list(
    FitSummary           = fit_summary,
    R2Coarse             = r2_obj$R2Coarse,
    R2GeneticSubblocks   = r2_obj$R2GeneticSubblocks,
    R2ExposureCategories = r2_obj$R2ExposureCategories,
    R2GxECategories      = r2_obj$R2GxECategories,
    R2Sensitivity        = r2_obj$R2Sensitivity,
    OOFPred              = oof_pred,
    ScaleFold            = scale_tbl,
    MediationScore       = score_obj$mediation_scores,
    MediationAudit       = mediation_audit,
    ScoreTrain           = score_train,
    ScoreTest            = score_test
  )
}

############################################################
# 10) Protein-level runner
############################################################

run_protein_predictive_r2 <- function(protID,
                                      pxs,
                                      op_resolve,
                                      family_type = c("lasso", "ridge", "enet"),
                                      k = cfg$kfold,
                                      cfg = cfg) {
  family_type <- match.arg(family_type)

  fold_obj <- CV_split(
    protID = protID,
    pxs = pxs,
    op_resolve = op_resolve,
    kfold = k,
    cfg = cfg
  )

  out <- vector("list", k)

  for (i in seq_len(k)) {
    message(
      "Protein ", protID,
      " [", family_type, "] fold ", i, "/", k,
      " decomposition=score_partition"
    )

    out[[i]] <- run_fold_score_partition(
      protID = protID,
      fold_obj = fold_obj,
      idx = i,
      family_type = family_type,
      cfg = cfg
    )
  }

  list(
    FitSummary           = dplyr::bind_rows(lapply(out, `[[`, "FitSummary")),
    R2Coarse             = dplyr::bind_rows(lapply(out, `[[`, "R2Coarse")),
    R2GeneticSubblocks   = dplyr::bind_rows(lapply(out, `[[`, "R2GeneticSubblocks")),
    R2ExposureCategories = dplyr::bind_rows(lapply(out, `[[`, "R2ExposureCategories")),
    R2GxECategories      = dplyr::bind_rows(lapply(out, `[[`, "R2GxECategories")),
    R2Sensitivity        = dplyr::bind_rows(lapply(out, `[[`, "R2Sensitivity")),
    OOFPred              = dplyr::bind_rows(lapply(out, `[[`, "OOFPred")),
    ScaleFold            = dplyr::bind_rows(lapply(out, `[[`, "ScaleFold")),
    MediationScore       = dplyr::bind_rows(lapply(out, `[[`, "MediationScore")),
    MediationAudit       = dplyr::bind_rows(lapply(out, `[[`, "MediationAudit")),
    ScoreTrain           = dplyr::bind_rows(lapply(out, `[[`, "ScoreTrain")),
    ScoreTest            = dplyr::bind_rows(lapply(out, `[[`, "ScoreTest")),
    FoldAssign           = fold_obj$fold_assign,
    Meta                 = fold_obj$meta
  )
}

############################################################
# 11) Chunk runner
############################################################

run_chunk_predictive_r2 <- function(protlist,
                                    pxs,
                                    folder_id,
                                    idx,
                                    covar_spec = "base",
                                    family_type = c("lasso", "ridge", "enet"),
                                    cfg = cfg,
                                    cat_group_map = NULL) {
  family_type <- match.arg(family_type)
  pxs <- as_pxs(pxs)

  out_dir <- file.path(cfg$paths$out_root, folder_id, family_type)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  op_resolve <- make_omicspred_resolver(cfg$paths$omicspred_map)

  if (is.null(pxs$base_ec)) {
    message("Precomputing E_df + covars base table once for this chunk...")

    E_df_all <- pxs$Elist %>%
      purrr::reduce(full_join, by = "eid")

    base_ec <- merge(E_df_all, pxs$covars_df, by = "eid")

    pxs$base_ec <- base_ec
    pxs$E_ids_all <- setdiff(names(E_df_all), "eid")
  }

  fit_all        <- list()
  coarse_all     <- list()
  genetic_all    <- list()
  ecat_all       <- list()
  gxecat_all     <- list()
  sensitivity_all <- list()
  pred_all       <- list()
  scale_all      <- list()
  med_score_all  <- list()
  med_audit_all  <- list()
  score_test_all <- list()

  for (prot in protlist) {
    message("Running protein: ", prot)

    run_obj <- run_protein_predictive_r2(
      protID = prot,
      pxs = pxs,
      op_resolve = op_resolve,
      family_type = family_type,
      k = cfg$kfold,
      cfg = cfg
    )

    fit_all[[prot]]         <- run_obj$FitSummary
    coarse_all[[prot]]      <- run_obj$R2Coarse
    genetic_all[[prot]]     <- run_obj$R2GeneticSubblocks
    ecat_all[[prot]]        <- run_obj$R2ExposureCategories
    gxecat_all[[prot]]      <- run_obj$R2GxECategories
    sensitivity_all[[prot]] <- run_obj$R2Sensitivity
    pred_all[[prot]]        <- run_obj$OOFPred
    scale_all[[prot]]       <- run_obj$ScaleFold

    if (nrow(run_obj$MediationScore) > 0) {
      med_score_all[[prot]] <- run_obj$MediationScore
    }

    if (nrow(run_obj$MediationAudit) > 0) {
      med_audit_all[[prot]] <- run_obj$MediationAudit
    }

    if (nrow(run_obj$ScoreTest) > 0) {
      score_test_all[[prot]] <- run_obj$ScoreTest
    }

    saveRDS(
      run_obj,
      file = file.path(
        out_dir,
        paste0("PredictiveR2_", prot, "_", covar_spec, "_", family_type, "_score_partition.rds")
      )
    )
  }

  fwrite(
    dplyr::bind_rows(fit_all),
    file = file.path(out_dir, paste0("fit_summary_", idx, ".txt")),
    sep = "\t"
  )

  fwrite(
    dplyr::bind_rows(coarse_all),
    file = file.path(out_dir, paste0("predictive_r2_coarse_", idx, ".txt")),
    sep = "\t"
  )

  fwrite(
    dplyr::bind_rows(genetic_all),
    file = file.path(out_dir, paste0("predictive_r2_genetic_subblocks_", idx, ".txt")),
    sep = "\t"
  )

  fwrite(
    dplyr::bind_rows(ecat_all),
    file = file.path(out_dir, paste0("predictive_r2_exposure_categories_", idx, ".txt")),
    sep = "\t"
  )

  fwrite(
    dplyr::bind_rows(gxecat_all),
    file = file.path(out_dir, paste0("predictive_r2_gxe_categories_", idx, ".txt")),
    sep = "\t"
  )

  if (length(sensitivity_all) > 0) {
    sens_df <- dplyr::bind_rows(sensitivity_all)
    if (nrow(sens_df) > 0) {
      fwrite(
        sens_df,
        file = file.path(out_dir, paste0("predictive_r2_sensitivity_", idx, ".txt")),
        sep = "\t"
      )
    }
  }

  fwrite(
    dplyr::bind_rows(pred_all),
    file = file.path(out_dir, paste0("oof_predictions_", idx, ".txt")),
    sep = "\t"
  )

  fwrite(
    dplyr::bind_rows(scale_all),
    file = file.path(out_dir, paste0("scale_", idx, ".txt")),
    sep = "\t"
  )

  fwrite(
    tibble::tibble(
      idx = idx,
      covar_spec = covar_spec,
      family_type = family_type,
      kfold = cfg$kfold,
      miss_rate = cfg$miss_rate,
      seed = cfg$seed,
      glmnet_inner_kfold = cfg$glmnet_inner_kfold,
      estimand = "unique_held_out_predictive_r2_of_learned_component_scores",
      decomposition_mode = "score_partition",
      primary_model = "full_lasso_once_per_outer_fold: C + Gcis + Gtrans + E + Gcis:E + Gtrans:E",
      primary_decomposition = "score-level conditional drop R2 over C_score, G_score, PXS_total, GIS_total",
      secondary_genetics = cfg$run_genetic_subblocks,
      secondary_exposure_categories = cfg$run_exposure_categories,
      secondary_gxe_categories = cfg$run_gxe_categories,
      sensitivity_gxc = cfg$run_gxc_sensitivity,
      sensitivity_exc = cfg$run_exc_sensitivity,
      gxc_covars = paste(cfg$gxc_covars, collapse = "|"),
      exc_covars = paste(cfg$exc_covars, collapse = "|")
    ),
    file = file.path(out_dir, paste0("run_config_", idx, ".txt")),
    sep = "\t"
  )

  ##########################################################
  # Mediation score files
  ##########################################################

  score_dir <- file.path(out_dir, "mediation_scores")
  dir.create(score_dir, showWarnings = FALSE, recursive = TRUE)

  if (length(med_score_all) > 0) {
    score_df <- dplyr::bind_rows(med_score_all)

    if (!is.null(cat_group_map) && nrow(cat_group_map) > 0) {
      grp_map <- cat_group_map
      colnames(grp_map) <- c("source_category", "analysis_group")

      all_groups <- unique(grp_map$analysis_group)

      for (grp in all_groups) {
        grp_safe <- sanitize_block_name(grp)
        src_cats <- grp_map$source_category[grp_map$analysis_group == grp]
        src_cols <- paste0("PXS_", sanitize_block_name(src_cats))
        avail <- intersect(src_cols, names(score_df))

        if (length(avail) > 0) {
          score_df[[paste0("PXSgrp_", grp_safe)]] <- rowSums(
            score_df[, avail, drop = FALSE],
            na.rm = TRUE
          )
        }
      }
    }

    fwrite(
      score_df,
      file = file.path(score_dir, paste0("mediation_scores_", idx, ".txt")),
      sep = "\t"
    )

    message("Saved mediation scores: ", nrow(score_df), " rows for idx=", idx)

    pxs_cat_cols <- grep("^PXS_", names(score_df), value = TRUE)
    pxs_cat_cols <- setdiff(pxs_cat_cols, "PXS_total")
    pxsgrp_cols <- grep("^PXSgrp_", names(score_df), value = TRUE)
    gis_cat_cols <- grep("^GIS_", names(score_df), value = TRUE)
    gis_cat_cols <- setdiff(gis_cat_cols, "GIS_total")

    fwrite(
      tibble::tibble(
        idx = idx,
        covar_spec = covar_spec,
        family_type = family_type,
        kfold = cfg$kfold,
        seed = cfg$seed,
        score_type = "out_of_fold",
        score_origin = "full_lasso_score_partition",
        uses_filtered_exposures = TRUE,
        n_exposure_categories = length(pxs_cat_cols),
        exposure_categories = paste(pxs_cat_cols, collapse = "|"),
        n_gxe_categories = length(gis_cat_cols),
        gxe_categories = paste(gis_cat_cols, collapse = "|"),
        n_exposure_groups = length(pxsgrp_cols),
        exposure_groups = paste(pxsgrp_cols, collapse = "|")
      ),
      file = file.path(score_dir, paste0("mediation_score_config_", idx, ".txt")),
      sep = "\t"
    )

    if (length(med_audit_all) > 0) {
      fwrite(
        dplyr::bind_rows(med_audit_all),
        file = file.path(score_dir, paste0("mediation_score_mapping_audit_", idx, ".txt")),
        sep = "\t"
      )
    }
  }

  if (length(score_test_all) > 0) {
    score_partition_dir <- file.path(out_dir, "score_partition_oof")
    dir.create(score_partition_dir, showWarnings = FALSE, recursive = TRUE)

    fwrite(
      dplyr::bind_rows(score_test_all),
      file = file.path(score_partition_dir, paste0("score_partition_oof_", idx, ".txt")),
      sep = "\t"
    )
  }

  invisible(TRUE)
}

############################################################
# 12) Covariate specs and command-line entrypoint
############################################################

heap <- readRDS(cfg$paths$heap_rds)

if (!exists("heap_filter_exposures", mode = "function")) {
  stop(
    "heap_filter_exposures() not found. ",
    "Make sure 00_paths.R is sourced through HEAP_PATHS_FILE or auto-discovery."
  )
}

heap <- heap_filter_exposures(heap)
pxs0 <- as_pxs_baseline(heap)

# ============================================================
# CLI: manifest-driven or legacy positional arguments
#
# Manifest mode (preferred):
#   Rscript Module1_suggested.R --manifest <path> --array-index <N>
#
# Legacy positional mode (backward compatibility):
#   Rscript Module1_suggested.R <idx> <split_num> <covarType> <family>
# ============================================================

args <- commandArgs(trailingOnly = TRUE)

.parse_flag <- function(args, flag, default = NULL) {
  idx <- which(args == flag)
  if (length(idx) == 0 || idx[1] >= length(args)) return(default)
  args[idx[1] + 1L]
}

.is_manifest_mode <- length(args) >= 2 && args[1] == "--manifest"

if (.is_manifest_mode) {
  .manifest_path <- .parse_flag(args, "--manifest")
  .array_idx_str <- .parse_flag(args, "--array-index")

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

  idx         <- as.integer(.row$chunk_id)
  split_num   <- as.integer(.row$n_chunks)
  covarType   <- as.character(.row$covariate_set)
  family_type <- as.character(.row$family)
  experiment_name <- as.character(.row$experiment_name %||% "unknown")
  .sample_filter  <- as.character(.row$sample_filter %||% "none")

  # Override cfg from manifest row
  if (!is.na(.row$seed))               cfg$seed               <- as.integer(.row$seed)
  if (!is.na(.row$kfold))              cfg$kfold              <- as.integer(.row$kfold)
  if (!is.na(.row$glmnet_inner_kfold)) cfg$glmnet_inner_kfold <- as.integer(.row$glmnet_inner_kfold)
  if (!is.na(.row$miss_rate))          cfg$miss_rate          <- as.numeric(.row$miss_rate)
  if (!is.na(.row$decomposition_mode)) cfg$decomposition_mode <- as.character(.row$decomposition_mode)

  cfg$run_genetic_subblocks   <- as.integer(.row$run_genetic_subblocks   %||% 1L) == 1L
  cfg$run_exposure_categories <- as.integer(.row$run_exposure_categories  %||% 1L) == 1L
  cfg$run_gxe_categories      <- as.integer(.row$run_gxe_categories       %||% 1L) == 1L
  cfg$run_gxc_sensitivity     <- as.integer(.row$run_gxc_sensitivity      %||% 0L) == 1L
  cfg$run_exc_sensitivity     <- as.integer(.row$run_exc_sensitivity      %||% 0L) == 1L

  if (!is.na(.row$gxc_covars) && nzchar(.row$gxc_covars))
    cfg$gxc_covars <- trimws(strsplit(.row$gxc_covars, ",", fixed = TRUE)[[1]])
  if (!is.na(.row$exc_covars) && nzchar(.row$exc_covars))
    cfg$exc_covars <- trimws(strsplit(.row$exc_covars, ",", fixed = TRUE)[[1]])

  # Output root: use IGLOO-rooted path from manifest
  if (!is.null(.row$output_path) && !is.na(.row$output_path) && nzchar(.row$output_path))
    cfg$paths$out_root <- as.character(.row$output_path)

  message(sprintf("[manifest] experiment=%s  array_index=%d  chunk=%d/%d  covar=%s  family=%s",
                  experiment_name, .arr_idx, idx, split_num, covarType, family_type))
  message("[manifest] out_root: ", cfg$paths$out_root)

} else {
  # ---- Legacy positional mode ----
  if (length(args) < 4) {
    stop(
      "Usage (manifest):   Rscript Module1_suggested.R --manifest <path> --array-index <N>\n",
      "Usage (positional): Rscript Module1_suggested.R <idx> <split_num> <covarType> <family>\n",
      "  family   : lasso | ridge | enet\n",
      "  covarType: base | base_bmi | base_draw | base_clinical | base_prevalent\n",
      "  (sample filters require manifest mode)"
    )
  }
  idx         <- as.integer(args[1])
  split_num   <- as.integer(args[2])
  covarType   <- as.character(args[3])
  family_type <- as.character(args[4])
  experiment_name <- paste0("positional_", covarType, "_", family_type)
  .sample_filter  <- "none"
}

# ============================================================
# Validate family
# ============================================================

allowed_families <- names(cfg$model_families)
if (!family_type %in% allowed_families)
  stop("family must be one of: ", paste(allowed_families, collapse = ", "))

# ============================================================
# Resolve covariate set
# Prefer centralized config (covariate_sets.yml); fall back to inline.
# ============================================================

.resolve_covartype <- function(covarType, pxs0) {
  .cfg_path <- heap_config("covariates", "covariate_sets.yml")
  if (file.exists(.cfg_path) && requireNamespace("yaml", quietly = TRUE)) {
    tryCatch({
      .sets <- yaml::read_yaml(.cfg_path)$covariate_sets
      if (covarType %in% names(.sets)) {
        .entry <- .sets[[covarType]]
        .covars <- .entry$covariates
        if (is.null(.covars) || (length(.covars) == 1L && is.na(.covars[[1L]]))) {
          # NULL sentinel = use full loader covars_list (Type5 pattern)
          return(pxs0$covars_list)
        }
        return(as.character(.covars))
      }
    }, error = function(e) NULL)
  }
  # Inline fallback (used only if covariate_sets.yml is unreadable). Mirrors the
  # descriptive sets in covariate_sets.yml EXACTLY. base does NOT include BMI or
  # fasting. base_ses resolves to the base list here (its deprivation E->C remap is
  # Module2-only and not applied in Module1).
  .pcs  <- paste0("genetic_principal_components_f22009_0_", 1:20)
  .base <- c(
    "age_when_attended_assessment_centre_f21003_0_0", "sex_f31_0_0",
    "age2", "age_sex", "age2_sex",
    "uk_biobank_assessment_centre_f54_0_0", .pcs
  )
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

covars_vec <- .resolve_covartype(covarType, pxs0)

# ------------------------------------------------------------
# Sample filter (the WHO-IS-IN axis), applied to the FULL covariate frame BEFORE
# any covariate-column narrowing (PXScovarSpec drops the filter column unless it is
# in the active set). Downstream complete-case/inner joins propagate the row
# restriction to the protein/exposure frames.
# ------------------------------------------------------------
.sf_spec <- load_sample_filter(.sample_filter)
if (!is.null(.sf_spec)) {
  .sf <- apply_sample_filter(pxs0$covars_df, .sf_spec)
  pxs0$covars_df <- .sf$df
  message(sprintf("[sample_filter] %s: %d -> %d (dropped %d)",
                  .sf_spec$name, .sf$n_before, .sf$n_after, .sf$n_dropped))
}

pxs_use    <- PXScovarSpec(pxs0, covars_vec)

# Guard: output must not be under scratch
if (exists("assert_no_canonical_scratch_output", mode = "function")) {
  assert_no_canonical_scratch_output(
    cfg$paths$out_root,
    context = paste0("Module1 out_root (", covarType, "/", family_type, ")")
  )
} else if (grepl(HEAP_PATHS$scratch_root, cfg$paths$out_root, fixed = TRUE)) {
  stop("REPRODUCIBILITY VIOLATION: Module1 out_root points to scratch: ",
       cfg$paths$out_root)
}

# ============================================================
# Protein list and chunk splitting
# ============================================================

omiclist  <- scan(file = cfg$paths$omic_list, what = character(), quiet = TRUE)

if (!is.finite(split_num) || split_num < 1)
  stop("split_num must be a positive integer")

split_vectors <- if (split_num >= length(omiclist)) {
  unname(as.list(omiclist))
} else {
  groups <- cut(seq_along(omiclist), breaks = split_num, labels = FALSE)
  split(omiclist, groups)
}

if (!is.finite(idx) || idx < 1 || idx > length(split_vectors))
  stop("idx must be between 1 and ", length(split_vectors), " for split_num=", split_num)

# ============================================================
# Category-group map (for grouped mediation mode)
# ============================================================

cat_group_map_path <- heap_config("exposure_sets", "analysis_exposure_category_groups.tsv")
cat_group_map <- if (file.exists(cat_group_map_path)) {
  message("Loading category-group mapping: ", cat_group_map_path)
  read.delim(cat_group_map_path, stringsAsFactors = FALSE, check.names = FALSE)
} else {
  NULL
}

# ============================================================
# Write resolved run config artifact
# ============================================================

.run_out_dir <- file.path(cfg$paths$out_root, covarType, family_type)
dir.create(.run_out_dir, recursive = TRUE, showWarnings = FALSE)

tryCatch({
  .rc <- list(
    experiment_name  = experiment_name,
    module           = "module1",
    array_index      = if (.is_manifest_mode) .arr_idx else NA_integer_,
    chunk_id         = idx,
    n_chunks         = split_num,
    covariate_set    = covarType,
    family           = family_type,
    seed             = cfg$seed,
    kfold            = cfg$kfold,
    glmnet_inner_kfold = cfg$glmnet_inner_kfold,
    miss_rate        = cfg$miss_rate,
    decomposition_mode = cfg$decomposition_mode,
    run_genetic_subblocks   = cfg$run_genetic_subblocks,
    run_exposure_categories = cfg$run_exposure_categories,
    run_gxe_categories      = cfg$run_gxe_categories,
    run_gxc_sensitivity     = cfg$run_gxc_sensitivity,
    run_exc_sensitivity     = cfg$run_exc_sensitivity,
    out_root         = cfg$paths$out_root,
    heap_rds         = cfg$paths$heap_rds,
    manifest_path    = if (.is_manifest_mode) .manifest_path else NA_character_,
    run_timestamp    = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    run_host         = Sys.info()[["nodename"]]
  )
  if (requireNamespace("yaml", quietly = TRUE)) {
    yaml::write_yaml(.rc, file.path(.run_out_dir, paste0("run_config_", idx, ".yml")))
  } else {
    saveRDS(.rc, file.path(.run_out_dir, paste0("run_config_", idx, ".rds")))
  }
}, error = function(e) {
  warning("Could not write run_config: ", conditionMessage(e))
})

# ============================================================
# Run
# ============================================================

run_chunk_predictive_r2(
  protlist      = split_vectors[[idx]],
  pxs           = pxs_use,
  folder_id     = covarType,
  idx           = idx,
  covar_spec    = covarType,
  family_type   = family_type,
  cfg           = cfg,
  cat_group_map = cat_group_map
)