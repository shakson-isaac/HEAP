############################################################
# PGS + PXS + Shapley R2 partitioning: Module 1
############################################################

# Libraries -------------------------------------------------------------

library(data.table)
library(tidyverse)  # dplyr, tidyr, purrr, tibble, etc.
library(ggplot2)
library(ggpmisc)
library(glmnet)
library(car)
library(caret)      # for createFolds
library(progress)

set.seed(123)

############################################################
# Class + Loader
############################################################

PXSconstruct <- setClass(
  "PXSconstruct",
  slots = c(
    Elist       = "list",      # E dataframes separated by category
    Elist_names = "character", # Category names
    Eid_cat     = "data.frame",# Dataframe with environmental IDs + category
    ordinalIDs  = "character", # List of ordinal environmental variable names
    
    UKBprot_df  = "data.frame",# Protein variables dataframe
    protIDs     = "character", # Protein IDs
    
    covars_df   = "data.frame",# Covariates dataframe
    covars_list = "character"  # Covariate names
  )
)

# Load pre-constructed object (paths as in your original script)
PXSloader <- readRDS(
  file = "/n/scratch/users/s/shi872/UKB_intermediate/UKB_PGS_PXS_load.rds"
)

############################################################
# Genetic score loading
############################################################

# ---- NEW: cache OmicsPred mapping once (avoid repeated fread for every protein) ----
.OMICSPRED_MAP <- NULL
get_omicspred_map <- function() {
  if (!is.null(.OMICSPRED_MAP)) return(.OMICSPRED_MAP)
  .OMICSPRED_MAP <<- data.table::fread(
    file = "/n/groups/patel/shakson_ukb/UK_Biobank/Data/OMICSPRED/UKB_Olink_multi_ancestry_models_val_results_portal.csv"
  )
  .OMICSPRED_MAP
}

read_sscore_or_zero <- function(fpath, out_col, eids, strict = FALSE) {
  
  # helper: make the fallback table (ALWAYS has eid + out_col)
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
  
  # --- identify ID column ---
  id_candidates <- c("eid", "IID", "#IID", "id", "ID")
  id_col <- id_candidates[id_candidates %in% names(dt)][1]
  if (is.na(id_col)) id_col <- names(dt)[1]  # fallback: first column
  
  # --- identify score column ---
  score_candidates <- c("SCORE1_AVG", "SCORE1_SUM", "SCORE1", "score", "SCORE")
  score_col <- score_candidates[score_candidates %in% names(dt)][1]
  if (is.na(score_col)) score_col <- names(dt)[ncol(dt)]  # fallback: last column
  
  # keep only eid + score
  out <- dt[, .(eid = as.integer(get(id_col)),
                score = as.numeric(get(score_col)))]
  
  data.table::setnames(out, c("eid", out_col))
  out
}

extract_protGScis <- function(protID, eids, strict = FALSE){
  omicpredIDs <- get_omicspred_map()
  dir <- "/n/groups/patel/IGLOO/UKB/ProtGScis/"
  
  opID <- omicpredIDs$OMICSPRED_ID[match(protID, omicpredIDs$Gene)]
  if (length(opID) == 0 || is.na(opID) || opID == "") {
    stop("No OMICSPRED_ID found for protID=", protID)
  }
  
  opFile <- paste0(opID, ".sscore")
  fpath <- file.path(dir, opFile)
  
  protID_clean <- gsub("-", "_", protID)
  out_col <- paste0(protID_clean, "_GScis")
  
  read_sscore_or_zero(fpath, out_col = out_col, eids = eids, strict = strict)
}

extract_protGStrans <- function(protID, eids, strict = FALSE){
  omicpredIDs <- get_omicspred_map()
  dir <- "/n/groups/patel/IGLOO/UKB/ProtGStrans/"
  
  opID <- omicpredIDs$OMICSPRED_ID[match(protID, omicpredIDs$Gene)]
  if (length(opID) == 0 || is.na(opID) || opID == "") {
    stop("No OMICSPRED_ID found for protID=", protID)
  }
  
  opFile <- paste0(opID, ".sscore")
  fpath <- file.path(dir, opFile)
  
  protID_clean <- gsub("-", "_", protID)
  out_col <- paste0(protID_clean, "_GStrans")
  
  read_sscore_or_zero(fpath, out_col = out_col, eids = eids, strict = strict)
}

GS_struct <- function(protID, UKBprot_df, strict = FALSE){
  protID_clean <- gsub("-", "_", protID)
  
  # eid universe comes from the protein measurement table
  omic_orig <- UKBprot_df %>% dplyr::select(all_of(c("eid", protID_clean)))
  eids <- omic_orig$eid
  
  # safe score loads (cis may be all zeros if missing)
  omic_PGScis  <- extract_protGScis(protID, eids = eids, strict = strict)
  omic_PGStrans <- extract_protGStrans(protID, eids = eids, strict = strict)
  
  # keep all eids (don’t accidentally inner-merge away participants)
  omic_PGS <- merge(omic_PGScis, omic_PGStrans, by = "eid", all = TRUE)
  
  # if trans exists but has holes, you can leave NA (so complete.cases drops them),
  # but for cis-missing we already set 0. Optionally:
  # omic_PGS[is.na(get(paste0(protID_clean,"_GScis"))), (paste0(protID_clean,"_GScis")) := 0]
  
  omicGS <- list()
  omicGS[["combo"]] <- na.omit(merge(omic_PGS, omic_orig, by = "eid"))
  omicGS[["solo"]]  <- na.omit(omic_PGS)
  
  omicGS
}


# extract_protGScis <- function(protID){
#   omicpredIDs <- get_omicspred_map()
#   dir <- "/n/groups/patel/IGLOO/UKB/ProtGScis/"
#   opID <- omicpredIDs$OMICSPRED_ID[match(protID, omicpredIDs$Gene)]
#   if (length(opID) == 0 || is.na(opID) || opID == "") {
#     stop("No OMICSPRED_ID found for protID=", protID)
#   }
# 
#   opFile <- paste0(opID, ".sscore")
#   print(opFile)
#   protGS <- data.table::fread(file = file.path(dir, opFile))
#   colnames(protGS) <- c("eid", paste0(protID,"_GScis"))
# 
#   # Convert dash to underscore
#   colnames(protGS) <- gsub("-", "_", colnames(protGS))
#   protGS
# }
# 
# extract_protGStrans <- function(protID){
#   omicpredIDs <- get_omicspred_map()
#   dir <- "/n/groups/patel/IGLOO/UKB/ProtGStrans/"
#   opID <- omicpredIDs$OMICSPRED_ID[match(protID, omicpredIDs$Gene)]
#   if (length(opID) == 0 || is.na(opID) || opID == "") {
#     stop("No OMICSPRED_ID found for protID=", protID)
#   }
# 
#   opFile <- paste0(opID, ".sscore")
#   protGS <- data.table::fread(file = file.path(dir, opFile))
#   colnames(protGS) <- c("eid", paste0(protID,"_GStrans"))
# 
#   colnames(protGS) <- gsub("-", "_", colnames(protGS))
#   protGS
# }
# 
# GS_struct <- function(protID, UKBprot_df){
#   # Load cis + trans polygenic scores
#   omic_PGScis <- extract_protGScis(protID)
#   omic_PGStrans <- extract_protGStrans(protID)
#   omic_PGS <- merge(omic_PGScis, omic_PGStrans, by = "eid")
#   
#   # Convert dash to underscore for HLA genes in protein ID
#   protID_clean <- gsub("-", "_", protID)
#   omic_orig <- UKBprot_df %>% dplyr::select(all_of(c("eid", protID_clean)))
#   
#   omicGS <- list()
#   omicGS[["combo"]] <- na.omit(merge(omic_PGS, omic_orig, by = "eid"))
#   omicGS[["solo"]] <- na.omit(omic_PGS)
#   
#   omicGS
# }

############################################################
# Preprocessing helpers
############################################################

continuous_finder <- function(df){
  max_cols <- apply(df, 2, max, na.rm = TRUE)
  unique_vals <- sapply(df, function(x) length(unique(x[!is.na(x)])))
  # continuous if not binary AND (max > 5 OR >2 unique values)
  continuous <- names(max_cols[max_cols > 5 | unique_vals > 2])
  continuous
}

remove_missing_cols <- function(df, miss_rate = 0.2){
  NAcols <- colMeans(is.na(df))
  removecols <- names(NAcols[NAcols > miss_rate])
  df %>% dplyr::select(!all_of(removecols))
}

remove_cols_name <- function(df, miss_rate = 0.2){
  NAcols <- colMeans(is.na(df))
  names(NAcols[NAcols > miss_rate])
}

categorical_handler <- function(df, ordinal_names, ordinal_contrast = "treatment"){
  for (i in ordinal_names) {
    if (!i %in% colnames(df)) next
    
    # Treatment coding => use unordered factor (ordered factors default to contr.poly => .L/.Q/.C)
    df[[i]] <- factor(df[[i]], ordered = FALSE)
    
    # If there is 0 or 1 level in this fold, skip setting contrasts
    nlev <- nlevels(df[[i]])
    if (nlev <= 1L) {
      # no contrasts needed / possible
      next
    }
    
    if (ordinal_contrast == "treatment") {
      contrasts(df[[i]]) <- contr.treatment(nlev)
    } else if (ordinal_contrast == "sum") {
      contrasts(df[[i]]) <- contr.sum(nlev)
    }
  }
  df
}

# --- UPDATED: factor schema capture + application (stable model.matrix columns) ---
# Store ONLY factor levels (not raw contrast matrices) for robustness.
# Contrasts are derived from the training fold and applied to test/new data via contrasts.arg.
get_factor_schema <- function(df) {
  facs <- names(df)[vapply(df, is.factor, logical(1))]
  levels_map <- lapply(facs, function(v) levels(df[[v]]))
  names(levels_map) <- facs
  list(levels = levels_map)
}

apply_factor_schema <- function(df, schema) {
  if (is.null(schema) || is.null(schema$levels)) return(df)

  for (v in names(schema$levels)) {
    if (!v %in% names(df)) next

    # enforce factor + training level set (stable dummy columns)
    df[[v]] <- factor(df[[v]], levels = schema$levels[[v]])
  }

  df
}

# --- NEW: coerce known integer-coded categoricals to factor BEFORE schema capture ---
# Note: user confirmed sex_f31_0_0 is a factor.
coerce_known_categoricals <- function(
    df,
    as_factor = c(
      "uk_biobank_assessment_centre_f54_0_0",
      "sex_f31_0_0"
    )
) {
  for (v in as_factor) {
    if (!v %in% names(df)) next
    df[[v]] <- as.factor(df[[v]])
  }
  df
}

build_formula <- function(protID, covars_list, E_ids){
  protID_clean <- gsub("-", "_", protID)
  G_cis <- paste0(protID_clean, "_GScis")
  G_trans <- paste0(protID_clean, "_GStrans")
  
  covariates <- covars_list
  
  # E and GxE components (cis + trans)
  pred_E <- c()
  for(i in E_ids){
    pred_Ei <- c(
      i,
      paste(G_cis, "*", i),
      paste(G_trans, "*", i)
    )
    pred_E <- c(pred_E, pred_Ei)
  }
  
  # IMPORTANT: keep intercept (remove "0")
  pred_var <- c(G_cis, G_trans, pred_E, covariates)
  as.formula(paste(protID_clean, paste(pred_var, collapse="+"), sep="~"))
}


preprocess_data <- function(protID, covariates,
                            trainData, testData,
                            ordinal_contrast = "treatment",
                            ordinal_names) {
  protID_clean <- gsub("-", "_", protID)
  G_cis <- paste0(protID_clean,"_GScis")
  G_trans <- paste0(protID_clean,"_GStrans")

  # NOTE: na removal now happens in CV_split() via complete.cases(),
  # so do NOT drop rows here (keeps fold membership fixed).
  trainData <- as.data.frame(trainData)
  testData  <- as.data.frame(testData)

  # continuous GS + protein + covariates
  rel_columns <- c(protID_clean, G_cis, G_trans, covariates)
  numeric_cols <- sapply(trainData[, rel_columns, drop = FALSE],
                         function(col) is.numeric(col) && !(all(col %in% c(0, 1))))
  DG_cont_names <- names(numeric_cols[numeric_cols])

  # continuous E features (not ordinal, not covars, not GS or omic)
  df_E <- trainData %>%
    dplyr::select(!all_of(c("eid", rel_columns, ordinal_names)))
  E_cont_names <- continuous_finder(df_E)

  numeric_col_names <- c(DG_cont_names, E_cont_names)

  mean_train <- lapply(trainData[numeric_col_names], function(x) mean(x, na.rm = TRUE))
  sd_train <- lapply(trainData[numeric_col_names], function(x) sd(x, na.rm = TRUE))

  for (i in numeric_col_names) {
    s <- sd_train[[i]]
    if (!is.finite(s) || s == 0) {
      # avoid Inf/NaN; treat as unscaled constant
      trainData[[i]] <- 0
      testData[[i]]  <- 0
    } else {
      trainData[[i]] <- (trainData[[i]] - mean_train[[i]]) / s
      testData[[i]]  <- (testData[[i]] - mean_train[[i]]) / s
    }
  }

  list(
    train = trainData,
    test = testData,
    numeric_col_names = numeric_col_names,
    mean_train = mean_train,
    sd_train = sd_train
  )
}


############################################################
# Lasso helpers
############################################################

extract_contrasts <- function(data) {
  contrasts_list <- list()
  for (colname in names(data)) {
    if (is.factor(data[[colname]])) {
      contrasts_list[[colname]] <- contrasts(data[[colname]], contrasts = TRUE)
    }
  }
  contrasts_list
}

lasso_fit <- function(protID, trainData, formula){
  protID_clean <- gsub("-", "_", protID)
  
  contrast_list <- extract_contrasts(trainData)
  x <- model.matrix(formula, data = trainData, contrasts.arg = contrast_list)
  y <- trainData[[protID_clean]]
  
  cv_lasso <- cv.glmnet(x, y, alpha = 1)
  best_lambda <- cv_lasso$lambda.min
  
  lasso_model <- glmnet(x, y, alpha = 1, lambda = best_lambda)
  lasso_coef <- coef(lasso_model)
  non_zero_coef <- which(lasso_coef != 0)[-1]  # exclude intercept
  
  selected_features <- rownames(lasso_coef)[non_zero_coef]
  selected_coefficients <- lasso_coef[non_zero_coef]
  
  list(
    lasso = lasso_model,
    lambda.min = best_lambda,
    select.features = selected_features,
    select.coefficients = selected_coefficients
  )
}

# R2 helper used both for lasso fit and Shapley
r2_score <- function(y, yhat) {
  ss_res <- sum((y - yhat)^2)
  ss_tot <- sum((y - mean(y))^2)
  if (!is.finite(ss_tot) || ss_tot == 0) return(NA_real_)
  1 - ss_res / ss_tot
}

glmnet_predict_fun <- function(model, newX, lambda) {
  as.numeric(predict(model, newx = as.matrix(newX), s = lambda))
}

# ---- NEW: compute sub-scores by component from coef_table + design matrix ----
score_component_contributions_from_coef_table <- function(X, coef_tbl) {
  # Returns per-row scores by component plus Intercept and PredTotal.
  # coef_tbl must include: term, beta, component.
  X <- as.matrix(X)

  intercept <- coef_tbl$beta[coef_tbl$term == "(Intercept)"][1]
  if (length(intercept) == 0 || is.na(intercept)) intercept <- 0

  ct <- coef_tbl %>% dplyr::filter(term != "(Intercept)")
  if (nrow(ct) == 0) {
    out <- tibble::tibble(Intercept = rep(intercept, nrow(X)), PredTotal = rep(intercept, nrow(X)))
    return(out)
  }

  common <- intersect(colnames(X), ct$term)
  if (length(common) == 0) {
    out <- tibble::tibble(Intercept = rep(intercept, nrow(X)), PredTotal = rep(intercept, nrow(X)))
    return(out)
  }

  beta <- ct$beta[match(common, ct$term)]
  comp <- ct$component[match(common, ct$term)]

  contrib_mat <- sweep(X[, common, drop = FALSE], 2, beta, `*`)

  comps <- sort(unique(comp))
  comp_scores <- sapply(comps, function(g) rowSums(contrib_mat[, comp == g, drop = FALSE]))
  comp_scores <- as.data.frame(comp_scores)

  comp_scores$Intercept <- intercept
  comp_scores$PredTotal <- intercept + rowSums(comp_scores[, setdiff(names(comp_scores), c("Intercept", "PredTotal")), drop = FALSE])

  tibble::as_tibble(comp_scores)
}

############################################################
# Shapley R2 machinery
############################################################

# Mask all features not in S by replacing with baseline
mask_X <- function(X, S, baseline) {
  X_masked <- X
  all_cols <- colnames(X)
  
  if (length(S) == 0) {
    inactive <- all_cols
  } else {
    inactive <- setdiff(all_cols, S)
  }
  
  for (col in inactive) {
    X_masked[[col]] <- baseline[[col]]
  }
  X_masked
}

shapley_r2_generic <- function(
    model,
    X,
    y,
    predict_fun,
    n_perm = 100,
    baseline = c("mean", "zero")
) {
  baseline <- match.arg(baseline)
  X <- as.data.frame(X)
  feature_names <- colnames(X)
  p <- length(feature_names)
  
  # Set baseline values for each feature
  base_vals <- sapply(X, function(col) {
    if (baseline == "mean") {
      mean(col, na.rm = TRUE)
    } else {
      0
    }
  })
  names(base_vals) <- feature_names
  
  # Full R^2
  yhat_full <- predict_fun(model, X)
  r2_full <- r2_score(y, yhat_full)
  
  if (!is.finite(r2_full)) {
    stop("Degenerate R^2 (full model).")
  }
  
  contrib <- setNames(numeric(p), feature_names)
  
  # Monte Carlo over permutations
  for (b in seq_len(n_perm)) {
    perm_idx <- sample(p)
    perm <- feature_names[perm_idx]
    S <- character(0)
    
    X_S <- mask_X(X, S, baseline = base_vals)
    yhat_S <- predict_fun(model, X_S)
    r2_S <- r2_score(y, yhat_S)
    
    for (j_name in perm) {
      S_new <- c(S, j_name)
      
      X_S_new <- mask_X(X, S_new, baseline = base_vals)
      yhat_S_new <- predict_fun(model, X_S_new)
      r2_S_new <- r2_score(y, yhat_S_new)
      
      contrib[j_name] <- contrib[j_name] + (r2_S_new - r2_S)
      
      S <- S_new
      r2_S <- r2_S_new
    }
  }
  
  contrib <- contrib / n_perm
  
  X_empty <- mask_X(X, character(0), baseline = base_vals)
  r2_empty <- r2_score(y, predict_fun(model, X_empty))
  
  list(
    shapley = contrib,
    r2_full = r2_full,
    r2_empty = r2_empty,
    sum_contrib = sum(contrib)
  )
}

# ---- UPDATED: main-effect Shapley only (interaction Shapley removed for deployment safety) ----
shapley_r2_full_decomp <- function(
    model,
    X,
    y,
    predict_fun,
    group_map = NULL,
    n_perm = 200,
    baseline = "mean"
) {
  shap_main <- shapley_r2_generic(
    model = model,
    X = X,
    y = y,
    predict_fun = predict_fun,
    n_perm = n_perm,
    baseline = baseline
  )

  phi <- shap_main$shapley

  main_tbl <- tibble::tibble(
    feature = names(phi),
    r2_main_pure = as.numeric(phi),
    r2_main_shapley = as.numeric(phi)
  )

  if (is.null(group_map)) {
    return(list(
      main = main_tbl,
      r2_full = shap_main$r2_full,
      r2_empty = shap_main$r2_empty
    ))
  }

  group_main <- main_tbl %>%
    dplyr::mutate(group = unname(group_map[feature])) %>%
    dplyr::group_by(group) %>%
    dplyr::summarise(
      r2_main = sum(r2_main_pure, na.rm = TRUE),
      .groups = "drop"
    )

  list(
    main = main_tbl,
    groups = group_main,
    r2_full = shap_main$r2_full,
    r2_empty = shap_main$r2_empty
  )
}

############################################################
# Group map: map model-matrix columns -> groups (G, E, GxE, categories, covars)
############################################################

build_group_map <- function(
    colnames_X,
    protID,
    E_ids,
    Eid_cat,
    Elist_names,
    covars_list
) {
  protID_clean <- gsub("-", "_", protID)
  G_cis <- paste0(protID_clean,"_GScis")
  G_trans <- paste0(protID_clean,"_GStrans")
  
  gm <- setNames(rep(NA_character_, length(colnames_X)), colnames_X)
  
  # Genetics main
  gm[colnames_X == G_cis] <- "Gcis"
  gm[colnames_X == G_trans] <- "Gtrans"
  
  # Covariates (exact-name columns, e.g. continuous covars)
  gm[colnames_X %in% covars_list] <- "Covars"
  
  # Helper: base name before ":" (for interactions)
  base_name <- function(x) sub("^(.*?)(\\:.*)?$", "\\1", x)
  
  # E main (continuous E columns that keep their name)
  is_E_main <- base_name(colnames_X) %in% E_ids
  gm[is_E_main & is.na(gm)] <- "E"
  
  # GxEcis / GxEtrans (global)
  gm[grepl(paste0("^", G_cis, ":"), colnames_X)] <- "GxEcis"
  gm[grepl(paste0("^", G_trans, ":"), colnames_X)] <- "GxEtrans"

  # Optional: category-level groups
  for (cat in Elist_names) {
    cat_E <- Eid_cat$Eid[Eid_cat$Category == cat]

    # Category E group (handle dummy-coded columns too)
    is_cat_E <- base_name(colnames_X) %in% cat_E |
      vapply(base_name(colnames_X), function(nm) any(startsWith(nm, cat_E)), logical(1))
    gm[is_cat_E] <- paste0("E_", cat)

    # Category GxEcis (interaction partner can also be dummy-coded)
    is_cat_GxEcis <- grepl(paste0("^", G_cis, ":"), colnames_X) &
      vapply(sapply(strsplit(colnames_X, ":"), `[`, 2), function(nm) any(startsWith(nm, cat_E)), logical(1))
    gm[is_cat_GxEcis] <- paste0("GxEcis_", cat)

    # Category GxEtrans
    is_cat_GxEtrans <- grepl(paste0("^", G_trans, ":"), colnames_X) &
      vapply(sapply(strsplit(colnames_X, ":"), `[`, 2), function(nm) any(startsWith(nm, cat_E)), logical(1))
    gm[is_cat_GxEtrans] <- paste0("GxEtrans_", cat)
  }
  
  ## ---- NEW: catch dummy-coded columns ----
  # Any column whose name starts with a covariate name -> Covars
  for (cv in covars_list) {
    gm[startsWith(colnames_X, cv) & is.na(gm)] <- "Covars"
  }
  
  # Any column whose *base* name starts with an E id -> E
  # (this pulls in stuff like E_Exercise_Freq2, etc.)
  for (e in E_ids) {
    gm[startsWith(base_name(colnames_X), e) & is.na(gm)] <- "E"
  }
  
  # Optional: put any remaining unlabelled columns into "Other"
  gm[is.na(gm)] <- "Other"
  
  gm
}

############################################################
# CV splitting (same as your original structure)
############################################################

CV_split <- function(protID, PXSdata, kfold){
  omicDS <- GS_struct(protID, PXSdata@UKBprot_df)
  protID_clean <- gsub("-", "_", protID)

  E_df <- PXSdata@Elist %>% purrr::reduce(full_join, by = "eid")

  df <- merge(omicDS$combo, E_df, by = "eid")
  df <- merge(df, PXSdata@covars_df, by = "eid")

  # NEW: ensure base data.frame so subsetting behaves consistently
  df <- as.data.frame(df)

  # Remove columns with too much missingness
  removecols <- remove_cols_name(df)
  df <- remove_missing_cols(df)
  print(paste0("Too many missing values: ", paste(removecols, collapse = ", ")))

  # Environmental ids
  E_ids <- colnames(E_df)[-1]
  E_ids <- E_ids[!(E_ids %in% removecols)]

  # ---- NEW: keep covariates consistent with columns actually retained after missingness filter ----
  covars_used <- PXSdata@covars_list[!(PXSdata@covars_list %in% removecols)]
  covars_used <- covars_used[covars_used %in% names(df)]

  # Coerce known integer-coded categoricals before creating folds
  df <- coerce_known_categoricals(df)

  # Handle ordinals before splitting
  df <- categorical_handler(df, PXSdata@ordinalIDs, ordinal_contrast = "treatment")

  # Capture factor schema once (stable columns across folds/deployment)
  factor_schema <- get_factor_schema(df)

  # ---- NEW: Drop incomplete rows ONCE before splitting into folds ----
  # Keep only columns that can appear in the modeling formula (outcome + predictors + eid)
  model_cols <- unique(c(
    "eid",
    protID_clean,
    paste0(protID_clean, "_GScis"),
    paste0(protID_clean, "_GStrans"),
    E_ids,
    covars_used
  ))
  model_cols <- intersect(model_cols, names(df))

  before_n <- nrow(df)
  df <- df[stats::complete.cases(df[, model_cols, drop = FALSE]), , drop = FALSE]
  after_n <- nrow(df)
  message("CV_split(): complete.cases filter kept ", after_n, " / ", before_n, " rows.")

  # Now create folds on the cleaned data
  folds <- createFolds(seq_len(nrow(df)), k = kfold, list = TRUE)

  CVdata <- list()
  for(i in seq_along(folds)){
    valIndex <- folds[[i]]
    trainData <- df[-valIndex, , drop = FALSE]
    testData  <- df[valIndex, , drop = FALSE]

    CVdata[["train"]][[i]] <- trainData
    CVdata[["test"]][[i]] <- testData
  }

  CVdata$Eids <- E_ids
  CVdata$covars_used <- covars_used
  CVdata$ordinalVar <- PXSdata@ordinalIDs
  CVdata$factor_schema <- factor_schema

  # store fold assignment for test rows = OOF mapping (now consistent)
  fold_assign <- purrr::map2_dfr(
    CVdata$test,
    seq_along(CVdata$test),
    ~ tibble::tibble(eid = .x$eid, fold = .y)
  )
  CVdata$fold_assign <- fold_assign

  CVdata
}

############################################################
# Per-fold run with Shapley R2
############################################################

CV_protPXSrun_shapley <- function(protID, folds, PXSdata, idx,
                                  n_perm = 200) {
  protID_clean <- gsub("-", "_", protID)
  G_cis <- paste0(protID_clean,"_GScis")
  G_trans <- paste0(protID_clean,"_GStrans")
  
  Elist_names <- PXSdata@Elist_names
  Eid_cat <- PXSdata@Eid_cat

  # Use the post-missingness covariate set captured in CV_split()
  covars_list <- folds[["covars_used"]]
  
  # Data for this fold
  protData <- list(
    train = folds[["train"]][[idx]],
    test = folds[["test"]][[idx]],
    Eids = folds[["Eids"]],
    covars_used = covars_list,
    ordinalVar = folds[["ordinalVar"]],
    factor_schema = folds[["factor_schema"]]
  )
  
  # Formula and data prep
  formula <- build_formula(protID, protData$covars_used, protData$Eids)
  
  pData <- preprocess_data(
    protID = protID,
    covariates = protData$covars_used,
    trainData = protData$train,
    testData = protData$test,
    ordinal_contrast = "treatment",
    ordinal_names = protData$ordinalVar
  )

  # NEW: enforce factor schema after na.omit (avoid level drift across folds)
  pData$train <- apply_factor_schema(pData$train, protData$factor_schema)
  pData$test <- apply_factor_schema(pData$test, protData$factor_schema)

  # Lasso fit
  lasso <- lasso_fit(protID, trainData = pData$train, formula = formula)

  # Model matrices (use training contrasts)
  contrast_list <- extract_contrasts(pData$train)
  X_train <- model.matrix(formula, data = pData$train, contrasts.arg = contrast_list)
  X_test <- model.matrix(formula, data = pData$test, contrasts.arg = contrast_list)
  
  y_train <- pData$train[[protID_clean]]
  y_test <- pData$test[[protID_clean]]
  
  # R2 of lasso fit
  ytrain_pred <- glmnet_predict_fun(lasso$lasso, X_train, lambda = lasso$lambda.min)
  ytest_pred <- glmnet_predict_fun(lasso$lasso, X_test, lambda = lasso$lambda.min)
  
  lasso_fit_metrics <- data.frame(
    omic = protID_clean,
    train_lasso = r2_score(y_train, ytrain_pred),
    test_lasso = r2_score(y_test, ytest_pred)
  )
  
  # Out-of-fold PXS
  test_PXS <- data.frame(
    eid = pData$test$eid,
    ID = protID_clean,
    PXS = ytest_pred
  )
  
  # Restrict X to selected, non-zero features for Shapley
  sel <- intersect(colnames(X_test), lasso$select.features)
  X_test_sel <- X_test[, sel, drop = FALSE]
  
  # Group map on selected features
  group_map <- build_group_map(
    colnames_X = colnames(X_test_sel),
    protID = protID,
    E_ids = protData$Eids,
    Eid_cat = Eid_cat,
    Elist_names = Elist_names,
    covars_list = protData$covars_used
  )
  
  # IMPORTANT: extract coefficients at lambda.min (includes intercept) and store intercept explicitly
  beta <- coef(lasso$lasso, s = lasso$lambda.min)
  beta_vec <- as.numeric(beta)
  names(beta_vec) <- rownames(beta)

  intercept <- unname(beta_vec["(Intercept)"])
  if (length(intercept) == 0 || is.na(intercept)) intercept <- 0

  coef_vec <- beta_vec[sel]
  coef_vec[is.na(coef_vec)] <- 0

  predict_fun_shap <- function(model, newX) {
    as.numeric(intercept + as.matrix(newX[, sel, drop = FALSE]) %*% coef_vec)
  }
  
  # Shapley R2 decomposition on test set (main effects only)
  shap_res <- shapley_r2_full_decomp(
    model = lasso$lasso,
    X = X_test_sel,
    y = y_test,
    predict_fun = predict_fun_shap,
    group_map = group_map,
    n_perm = n_perm,
    baseline = "mean"
  )
  
  group_df <- shap_res$groups
  group_df$omic <- protID_clean
  group_df$fold <- idx
  
  shap_main <- shap_res$main %>%
    dplyr::mutate(omic = protID_clean, fold = idx)
  
  # --- UPDATED: tidy coefficient table for this fold (store intercept as a row) ---
  base_name <- function(x) sub("^(.*?)(\\:.*)?$", "\\1", x)

  fold_coef_tbl <- tibble::tibble(
    protID = protID_clean,
    fold = idx,
    term = c("(Intercept)", sel),
    beta = c(intercept, as.numeric(coef_vec)),
    component = c("Intercept", unname(group_map[sel])),
    raw_var = c("(Intercept)", base_name(sel))
  ) %>%
    dplyr::mutate(
      side = dplyr::case_when(
        term == "(Intercept)" ~ "Intercept",
        component == "E" ~ "PXS_main",
        grepl("^GxEcis", component) ~ "GIS_cis",
        grepl("^GxEtrans", component) ~ "GIS_trans",
        component %in% c("Gcis", "Gtrans") ~ "G_main",
        component == "Covars" ~ "Covar",
        TRUE ~ "Other"
      )
    )
  
  # scaling table for this fold
  scale_table_fold <- tibble::tibble(
    protID = protID_clean,
    fold = idx,
    var = names(pData$mean_train),
    mean = unlist(pData$mean_train),
    sd = unlist(pData$sd_train)
  )
  
  list(
    R2.lasso.fit = lasso_fit_metrics,
    R2.groups = group_df,
    shap_main = shap_main,
    test.PXS = test_PXS,
    coef_table_fold = fold_coef_tbl,
    scale_table_fold = scale_table_fold,
    design_cols_fold = colnames(X_train),
    group_map_fold = group_map
  )
}

############################################################
# Run across folds for one protein
############################################################

prot_PXS_CVrun_shapley <- function(protID, PXSdata, k,
                                   n_perm = 200) {
  
  folds <- CV_split(protID, PXSdata, kfold = k)
  fold_assign <- folds$fold_assign
  factor_schema <- folds$factor_schema

  Lasso_all <- data.frame()
  R2groups_all <- data.frame()
  OOF_PXS_all <- data.frame()
  ShapMain_all <- data.frame()
  Coef_all <- data.frame()
  Scale_all <- data.frame()
  GroupMap_all <- vector("list", length = k)
  
  pb <- progress_bar$new(
    format = paste0("Protein ", protID, " [:bar] :percent | Fold :current/:total (:eta remaining)"),
    total = k,
    clear = FALSE,
    width = 60
  )
  
  design_cols_fold1 <- NULL

  for (i in seq_len(k)) {
    pb$tick()
    
    model <- CV_protPXSrun_shapley(
      protID = protID,
      folds = folds,
      PXSdata = PXSdata,
      idx = i,
      n_perm = n_perm
    )
    
    if (i == 1) {
      design_cols_fold1 <- model$design_cols_fold
    }

    Lasso_all <- rbind(Lasso_all, model$R2.lasso.fit)
    R2groups_all <- rbind(R2groups_all, model$R2.groups)
    OOF_PXS_all <- rbind(OOF_PXS_all, model$test.PXS)
    ShapMain_all <- rbind(ShapMain_all, model$shap_main)
    Coef_all <- rbind(Coef_all, model$coef_table_fold)
    Scale_all <- rbind(Scale_all, model$scale_table_fold)
    GroupMap_all[[i]] <- model$group_map_fold
  }
  
  list(
    Lasso = Lasso_all,
    R2groups = R2groups_all,
    OOF_PXS = OOF_PXS_all,
    ShapMain = ShapMain_all,
    CoefFold = Coef_all,
    ScaleFold = Scale_all,
    FoldAssign = fold_assign,
    FactorSchema = factor_schema,
    DesignColsFold1 = design_cols_fold1,
    EidsUsed = folds$Eids,
    CovarsUsed = folds$covars_used,
    GroupMap = GroupMap_all
  )
}

######################################
# Define PXSGIS model card for each protein
######################################

PXSGISModel <- setClass(
  "PXSGISModel",
  slots = c(
    protID = "character",
    covar_spec = "character",
    n_folds = "integer",
    coef_table = "data.frame",
    scale_table = "data.frame",
    meta = "list"
  )
)

make_PXSGISModel <- function(protID, covar_spec, k, run_obj, PXSdata) {
  PXSGISModel(
    protID = protID,
    covar_spec = covar_spec,
    n_folds = as.integer(k),
    coef_table = run_obj$CoefFold,
    scale_table = run_obj$ScaleFold,
    meta = list(
      Elist_names = PXSdata@Elist_names,
      Eid_cat = PXSdata@Eid_cat,
      covars_list = PXSdata@covars_list,
      covars_used = run_obj$CovarsUsed,
      E_ids_used = run_obj$EidsUsed,
      fold_assign = run_obj$FoldAssign,
      factor_schema = run_obj$FactorSchema,
      design_cols_fold1 = run_obj$DesignColsFold1,
      group_map_by_fold = run_obj$GroupMap 
    )
  )
}

############################################################
# Parallel runner for a chunk of proteins
############################################################

runProt_PGS_PXS_multi_shapley <- function(protlist, PXSdata, folder_id, idx,
                                          covar_spec = "Type3",
                                          n_perm = 200,
                                          save_oof_components = TRUE,
                                          save_oof_predtotal_only = FALSE) {
  
  out_dir <- file.path(
    "/n/groups/patel/shakson_ukb/UK_Biobank/Data/Parallel/Prot_PGSPXS_shapley",
    folder_id
  )
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  # --- helper: build the full (non-NA) dataset used for CV for this protein ---
  # This mirrors the CV data assembly: eid + protein + GS + E + covars, then complete.cases().
  build_full_df_for_prot <- function(protID, PXSdata) {
    protID_clean <- gsub("-", "_", protID)

    # Genetic predictors + (eid, protein)
    gs_list <- GS_struct(protID, PXSdata@UKBprot_df)
    df_base <- gs_list[["combo"]]  # contains eid + GScis + GStrans + protein

    # E data (already pre-split by category in PXSdata)
    if (!is.null(PXSdata@Elist) && length(PXSdata@Elist) > 0) {
      for (nm in names(PXSdata@Elist)) {
        df_base <- merge(df_base, PXSdata@Elist[[nm]], by = "eid", all.x = TRUE)
      }
    }

    # Covariates
    df_base <- merge(df_base, PXSdata@covars_df, by = "eid", all.x = TRUE)

    # Drop rows with any missing across the modeled columns (as CV_split() does)
    df_base <- df_base[stats::complete.cases(df_base), , drop = FALSE]
    df_base
  }
  
  # Summary collectors (optional)
  Lasso_all_list <- list()
  R2groups_all_list <- list()
  ShapMain_all_list <- list()
  
  for (prot in protlist) {
    message("Running protein: ", prot)
    
    # ---- run CV ----
    model <- prot_PXS_CVrun_shapley(
      protID = prot,
      PXSdata = PXSdata,
      k = 10,
      n_perm = n_perm
    )
    
    # ---- build the model card ----
    model_card <- make_PXSGISModel(
      protID = prot,
      covar_spec = covar_spec,
      k = 10L,
      run_obj = model,
      PXSdata = PXSdata
    )
    
    # ---- save it ----
    saveRDS(
      model_card,
      file = file.path(out_dir, paste0("PXSGISModel_", prot, "_", covar_spec, ".rds"))
    )

    # ---- OPTIONAL: save OOF component scores ----
    # NOTE: This can be large (N x #components). Use save_oof_predtotal_only to store only PredTotal.
    if (isTRUE(save_oof_components) || isTRUE(save_oof_predtotal_only)) {
      full_df <- build_full_df_for_prot(prot, PXSdata)
      oof_comp <- get_oof_components_from_card(model_card, full_df)

      if (isTRUE(save_oof_predtotal_only)) {
        oof_comp <- oof_comp %>% dplyr::select(eid, protID, fold, PredTotal)
      }

      saveRDS(
        oof_comp,
        file = file.path(out_dir, paste0("OOF_components_", prot, "_", covar_spec, ".rds"))
      )
    }
    
    # ---- OPTIONAL summary collectors ----
    Lasso_all_list[[prot]] <- model$Lasso
    R2groups_all_list[[prot]] <- model$R2groups
    ShapMain_all_list[[prot]] <- model$ShapMain
  }
  
  # Combine summaries (optional)
  Lasso_all <- dplyr::bind_rows(Lasso_all_list)
  R2groups_all <- dplyr::bind_rows(R2groups_all_list)
  ShapMain_all <- dplyr::bind_rows(ShapMain_all_list)
  
  # Write summaries (optional)
  write.table(Lasso_all,
              file = file.path(out_dir, paste0("lassofit_", idx, ".txt")),
              row.names = FALSE)
  
  write.table(R2groups_all,
              file = file.path(out_dir, paste0("R2groups_", idx, ".txt")),
              row.names = FALSE)
  
  write.table(ShapMain_all,
              file = file.path(out_dir, paste0("ShapMain_", idx, ".txt")),
              row.names = FALSE)
}

###### FOR DEPLOYMENT LATER (AS MODEL CARD)
get_deployment_params <- function(model_card) {
  coef_tbl <- model_card@coef_table

  # Build a complete grid of (fold x term) so missing terms are treated as 0 in the mean.
  folds <- sort(unique(coef_tbl$fold))

  term_meta <- coef_tbl %>%
    dplyr::group_by(term) %>%
    dplyr::arrange(dplyr::desc(abs(beta))) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::select(term, raw_var, component, side)

  full_grid <- tidyr::expand_grid(
    fold = folds,
    term = term_meta$term
  )

  coef_complete <- full_grid %>%
    dplyr::left_join(coef_tbl %>% dplyr::select(fold, term, beta), by = c("fold", "term")) %>%
    dplyr::mutate(beta = tidyr::replace_na(beta, 0)) %>%
    dplyr::left_join(term_meta, by = "term")

  coef_deploy <- coef_complete %>%
    dplyr::group_by(term, raw_var, component, side) %>%
    dplyr::summarise(beta = mean(beta), .groups = "drop")

  scale_deploy <- model_card@scale_table %>%
    dplyr::group_by(var) %>%
    dplyr::summarise(
      mean = mean(mean),
      sd = mean(sd),
      .groups = "drop"
    )

  list(
    coef = coef_deploy,
    scale = scale_deploy
  )
}

# helper: null-coalescing
`%||%` <- function(x, y) if (!is.null(x)) x else y

score_new_from_card <- function(model_card, new_df) {
  params <- get_deployment_params(model_card)
  coef_deploy <- params$coef
  scale_deploy <- params$scale

  protID <- model_card@protID
  protID_clean <- gsub("-", "_", protID)

  # Coerce known categoricals first (helps if deployment integers arrive as numeric)
  new_df <- coerce_known_categoricals(new_df)

  # Enforce factor level schema for stable dummy coding
  if (!is.null(model_card@meta$factor_schema)) {
    new_df <- apply_factor_schema(new_df, model_card@meta$factor_schema)
  }

  # 1) scale numeric vars as in training
  for (v in scale_deploy$var) {
    if (!v %in% names(new_df)) next
    m <- scale_deploy$mean[scale_deploy$var == v][1]
    s <- scale_deploy$sd[scale_deploy$var == v][1]
    if (!is.finite(s) || s == 0) next
    new_df[[v]] <- (new_df[[v]] - m) / s
  }

  # 2) rebuild model.matrix using same formula
  covars_list <- (model_card@meta$covars_used %||% model_card@meta$covars_list)
  E_ids <- (model_card@meta$E_ids_used %||% model_card@meta$Eid_cat$Eid)
  formula <- build_formula(protID, covars_list, E_ids)

  # Contrasts defined from the (schema-applied) new data
  contrast_list <- extract_contrasts(new_df)

  X_new <- model.matrix(formula, data = new_df, contrasts.arg = contrast_list)

  # 3) Apply averaged coefficients
  intercept_deploy <- coef_deploy$beta[coef_deploy$term == "(Intercept)"][1]
  if (length(intercept_deploy) == 0 || is.na(intercept_deploy)) intercept_deploy <- 0

  coef_terms <- coef_deploy %>% dplyr::filter(term != "(Intercept)")
  common_terms <- intersect(colnames(X_new), coef_terms$term)
  X_use <- X_new[, common_terms, drop = FALSE]
  beta <- coef_terms$beta[match(common_terms, coef_terms$term)]

  pred <- as.numeric(intercept_deploy + as.matrix(X_use) %*% beta)

  tibble::tibble(
    eid = new_df$eid,
    protID = protID_clean,
    predProt = pred
  )
}

get_oof_predictions_from_card <- function(model_card, full_df) {
  fold_assign <- model_card@meta$fold_assign
  protID <- model_card@protID
  protID_clean <- gsub("-", "_", protID)

  covars_list <- (model_card@meta$covars_used %||% model_card@meta$covars_list)
  E_ids <- (model_card@meta$E_ids_used %||% model_card@meta$Eid_cat$Eid)
  formula <- build_formula(protID, covars_list, E_ids)

  full_df2 <- full_df %>%
    dplyr::inner_join(fold_assign, by = "eid")

  res_list <- list()

  for (f in sort(unique(full_df2$fold))) {
    df_fold <- full_df2 %>% dplyr::filter(fold == f)

    # Coerce known categoricals then enforce stored factor levels
    df_fold <- coerce_known_categoricals(df_fold)
    if (!is.null(model_card@meta$factor_schema)) {
      df_fold <- apply_factor_schema(df_fold, model_card@meta$factor_schema)
    }

    coef_f <- model_card@coef_table %>% dplyr::filter(fold == f)
    scale_f <- model_card@scale_table %>% dplyr::filter(fold == f)

    # fold-specific scaling
    df_scaled <- df_fold
    for (v in scale_f$var) {
      if (!v %in% names(df_scaled)) next
      m <- scale_f$mean[scale_f$var == v][1]
      s <- scale_f$sd[scale_f$var == v][1]
      if (!is.finite(s) || s == 0) next
      df_scaled[[v]] <- (df_scaled[[v]] - m) / s
    }

    # Contrasts defined from the fold data AFTER schema enforcement
    contrast_list <- extract_contrasts(df_scaled)

    X <- model.matrix(formula, data = df_scaled, contrasts.arg = contrast_list)

    intercept_f <- coef_f$beta[coef_f$term == "(Intercept)"][1]
    if (length(intercept_f) == 0 || is.na(intercept_f)) intercept_f <- 0

    coef_terms_f <- coef_f %>% dplyr::filter(term != "(Intercept)")
    common <- intersect(colnames(X), coef_terms_f$term)
    X_use <- X[, common, drop = FALSE]
    beta <- coef_terms_f$beta[match(common, coef_terms_f$term)]

    pred <- as.numeric(intercept_f + as.matrix(X_use) %*% beta)

    res_list[[as.character(f)]] <- tibble::tibble(
      eid = df_fold$eid,
      fold = f,
      protID = protID_clean,
      predProt = pred
    )
  }

  dplyr::bind_rows(res_list)
}

############################################################
# Covariate specs + main entrypoint
############################################################

PXScovarSpec <- function(PXSdata, covars_subset){
  PXSdata@covars_list <- covars_subset
  PXSdata@covars_df <- PXSdata@covars_df %>%
    dplyr::select(all_of(c("eid", covars_subset)))
  PXSdata
}

CovarSpec <- list()
CovarSpec[["Type1"]] <- c(
  "age_when_attended_assessment_centre_f21003_0_0",
  "sex_f31_0_0"
)
CovarSpec[["Type2"]] <- c(
  "age_when_attended_assessment_centre_f21003_0_0",
  "sex_f31_0_0",
  "body_mass_index_bmi_f23104_0_0",
  "fasting_time_f74_0_0"
)
CovarSpec[["Type3"]] <- c(
  "age_when_attended_assessment_centre_f21003_0_0",
  "sex_f31_0_0",
  "age2","age_sex", "age2_sex",
  "body_mass_index_bmi_f23104_0_0",
  "fasting_time_f74_0_0",
  "uk_biobank_assessment_centre_f54_0_0",
  paste0("genetic_principal_components_f22009_0_",1:20)
)
CovarSpec[["Type4"]] <- c(
  "age_when_attended_assessment_centre_f21003_0_0",
  "sex_f31_0_0",
  "age2","age_sex", "age2_sex",
  "body_mass_index_bmi_f23104_0_0",
  "fasting_time_f74_0_0",
  "uk_biobank_assessment_centre_f54_0_0",
  paste0("genetic_principal_components_f22009_0_",1:20),
  "combined_Blood_pressure_medication",
  "combined_Hormone_replacement_therapy",
  "combined_Oral_contraceptive_pill_or_minipill",
  "combined_Insulin",
  "combined_Cholesterol_lowering_medication",
  "combined_Do_not_know",
  "combined_None_of_the_above",
  "combined_Prefer_not_to_answer"
)
CovarSpec[["Type5"]] <- PXSloader@covars_list

score_new_components_from_card <- function(model_card, new_df) {
  params <- get_deployment_params(model_card)
  coef_deploy <- params$coef
  scale_deploy <- params$scale

  protID <- model_card@protID
  protID_clean <- gsub("-", "_", protID)

  new_df <- coerce_known_categoricals(new_df)
  if (!is.null(model_card@meta$factor_schema)) {
    new_df <- apply_factor_schema(new_df, model_card@meta$factor_schema)
  }

  # scale as in deployment
  for (v in scale_deploy$var) {
    if (!v %in% names(new_df)) next
    m <- scale_deploy$mean[scale_deploy$var == v][1]
    s <- scale_deploy$sd[scale_deploy$var == v][1]
    if (!is.finite(s) || s == 0) next
    new_df[[v]] <- (new_df[[v]] - m) / s
  }

  covars_list <- (model_card@meta$covars_used %||% model_card@meta$covars_list)
  E_ids <- (model_card@meta$E_ids_used %||% model_card@meta$Eid_cat$Eid)
  formula <- build_formula(protID, covars_list, E_ids)

  contrast_list <- extract_contrasts(new_df)
  X_new <- model.matrix(formula, data = new_df, contrasts.arg = contrast_list)

  # component contributions (Intercept, Covars, E_*, GxE*, etc.) as defined in coef_deploy$component
  comp_scores <- score_component_contributions_from_coef_table(X_new, coef_deploy)

  dplyr::bind_cols(
    tibble::tibble(eid = new_df$eid, protID = protID_clean),
    comp_scores
  )
}

get_oof_components_from_card <- function(model_card, full_df) {
  fold_assign <- model_card@meta$fold_assign
  protID <- model_card@protID
  protID_clean <- gsub("-", "_", protID)

  covars_list <- (model_card@meta$covars_used %||% model_card@meta$covars_list)
  E_ids <- (model_card@meta$E_ids_used %||% model_card@meta$Eid_cat$Eid)
  formula <- build_formula(protID, covars_list, E_ids)

  full_df2 <- full_df %>% dplyr::inner_join(fold_assign, by = "eid")

  out <- vector("list", length = length(unique(full_df2$fold)))
  names(out) <- as.character(sort(unique(full_df2$fold)))

  for (f in sort(unique(full_df2$fold))) {
    df_fold <- full_df2 %>% dplyr::filter(fold == f)

    df_fold <- coerce_known_categoricals(df_fold)
    if (!is.null(model_card@meta$factor_schema)) {
      df_fold <- apply_factor_schema(df_fold, model_card@meta$factor_schema)
    }

    coef_f <- model_card@coef_table %>% dplyr::filter(fold == f)
    scale_f <- model_card@scale_table %>% dplyr::filter(fold == f)

    df_scaled <- df_fold
    for (v in scale_f$var) {
      if (!v %in% names(df_scaled)) next
      m <- scale_f$mean[scale_f$var == v][1]
      s <- scale_f$sd[scale_f$var == v][1]
      if (!is.finite(s) || s == 0) next
      df_scaled[[v]] <- (df_scaled[[v]] - m) / s
    }

    contrast_list <- extract_contrasts(df_scaled)
    X <- model.matrix(formula, data = df_scaled, contrasts.arg = contrast_list)

    comp_scores <- score_component_contributions_from_coef_table(X, coef_f)

    out[[as.character(f)]] <- dplyr::bind_cols(
      tibble::tibble(eid = df_fold$eid, fold = f, protID = protID_clean),
      comp_scores
    )
  }

  dplyr::bind_rows(out)
}

##### DEPLOYMENT/Parallel VERSION ####
# Protein list
omiclist <- scan(
  file = "/n/groups/patel/shakson_ukb/UK_Biobank/BScripts/ProtPGS_PXS/OMICPREDproteins.txt",
  what = character()
)

args <- commandArgs(trailingOnly = TRUE)
idx <- as.integer(args[1])
split_num <- as.integer(args[2])
covarType <- as.character(args[3])

idx = 787
split_num = 1000
covarType = "Type1"

# Split protein list into chunks
groups <- cut(seq_along(omiclist), breaks = split_num, labels = FALSE)
split_vectors <- split(omiclist, groups)

# Run for this chunk with chosen covariate spec
# runProt_PGS_PXS_multi_shapley(
#   protlist = split_vectors[[idx]],
#   PXSdata = PXScovarSpec(PXSloader, CovarSpec[[covarType]]),
#   folder_id = covarType,
#   idx = idx,
#   covar_spec = covarType,
#   n_perm = 10, #200,
#   save_oof_components = TRUE,
#   save_oof_predtotal_only = FALSE
# )

runProt_PGS_PXS_multi_shapley(
  protlist = split_vectors[[idx]],
  PXSdata = PXScovarSpec(PXSloader, CovarSpec[[covarType]]),
  folder_id = covarType,
  idx = idx,
  covar_spec = covarType,
  n_perm = 10, #200,
  save_oof_components = TRUE,
  save_oof_predtotal_only = FALSE
)

#RAD23B is a case where there is no cis-variant what do I do then?
xx1 <- readRDS("/n/groups/patel/shakson_ukb/UK_Biobank/Data/Parallel/Prot_PGSPXS_shapley/Type1/PXSGISModel_RAD23B_Type1.rds") 
xx2 <- readRDS("/n/groups/patel/shakson_ukb/UK_Biobank/Data/Parallel/Prot_PGSPXS_shapley/Type1/OOF_components_RAD23B_Type1.rds")
