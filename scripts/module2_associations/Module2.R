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

# Config helpers: load_sample_filter / apply_sample_filter (sample-filter axis) and
# load_covariate_set_remapping (base_ses deprivation E->C). Sourced from the same
# workflow dir as 00_paths.R.
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

# ============================================================
#  HEAP Univariate E + (Gcis, Gtrans) + (Gcis:E, Gtrans:E)
#  Optionally adds E x Covar and G x Covar interactions
#  (sensitivity analysis) via the int_covariates argument.
#
#  covariate sets base/base_bmi/base_draw/base_clinical/base_ses/base_prevalent
#                               -> HEAP/output/module2/
#  covar_variant sens_core      -> HEAP/output/module2_sens/  (core ExCov/GxCov)
#  covar_variant sens_extended  -> HEAP/output/module2_sens/  (extended ExCov/GxCov)
#  sample_filter exclude_prevalent -> healthy-at-baseline subset (reviewer exclusion)
#
#  Output per run: list(train, test), each a 4-element list:
#    [[1]] statE      : E main-effect coef table (factor-safe)
#    [[2]] statGxE    : GxE interaction coef table (cis + trans)
#    [[3]] statR2     : partitioned R2 (Type-1 SS; internal)
#    [[4]] statFblock : nested-model block F-tests
#             p_E_block   : drop E main effect (GxE kept); marginal E
#             p_GxE_joint : drop both GxE blocks
#             p_GcisxE    : drop cis interaction only
#             p_GtrxE     : drop trans interaction only
#             n_int_covars: number of ExCov/GxCov interaction covariates used
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(ggplot2)
  library(ggpmisc)
  library(car)
})

set.seed(123)

# ----------------------------
# Config: paths for cis/trans
# ----------------------------
CFG <- list(
  omicpred_map = heap_omicspred_or_legacy("UKB_Olink_multi_ancestry_models_val_results_portal.csv"),
  # Shared IGLOO genetics resources (already canonical — preserve)
  gs_cis_dir   = igloo_path("UKB", "ProtGScis"),
  gs_tr_dir    = igloo_path("UKB", "ProtGStrans"),
  # HEAP-specific outputs: IGLOO-rooted canonical location.
  out_root     = heap_project_output("module2"),
  out_root_sens = heap_project_output("module2_sens")
)

prot_clean <- function(protID) gsub("-", "_", protID)
`%||%` <- function(x, y) if (!is.null(x)) x else y

# ----------------------------
# Load data from HEAP.rds
# ----------------------------
if (!exists("heap_filter_exposures", mode = "function"))
  stop("heap_filter_exposures() not found. Ensure 00_paths.R is sourced correctly.")

HEAP <- readRDS(heap_loader_rds)
HEAP <- heap_filter_exposures(HEAP)
# as_pxs_baseline() now routes ordinalIDs by the DECLARED variable_type from
# analysis_exposures.tsv (see 00_paths.R), so small-scale continuous scores
# (income/employment/crime/health IMD scores) are no longer mis-factorised into
# hundreds of spurious polynomial terms. Exposure typing is config-driven, not
# value-range heuristics.
PXSloader <- as_pxs_baseline(HEAP)
rm(HEAP)

# ----------------------------
# Helpers
# ----------------------------
continuous_finder <- function(df){
  max_cols <- apply(df, 2, max, na.rm = TRUE)
  unique_vals <- sapply(df, function(x) length(unique(x[!is.na(x)])))
  names(max_cols[max_cols > 5 | unique_vals > 2])
}

categorical_handler <- function(df, ordinal_names, ordinal_contrast = "treatment"){
  for (i in ordinal_names) {
    if (!i %in% names(df)) next
    df[[i]] <- factor(df[[i]], ordered = TRUE)
    lv <- levels(df[[i]])
    if (length(lv) < 2) next

    if (ordinal_contrast == "treatment") {
      # Reference level = "0" so every coefficient is "level k vs 0". HEAP ordinals
      # are coded from 0 upward (verified: all 31 have min == 0); warn + fall back
      # to the lowest level if a "0" level is ever absent.
      base_idx <- match("0", lv)
      if (is.na(base_idx)) {
        warning("ordinal '", i, "' has no '0' level; using lowest level '", lv[1],
                "' as the treatment reference")
        base_idx <- 1L
      }
      contrasts(df[[i]]) <- contr.treatment(length(lv), base = base_idx)
    } else if (ordinal_contrast == "sum") {
      contrasts(df[[i]]) <- contr.sum(length(lv))
    }
  }
  df
}

# ------------------------------------------------------------
# OmicsPred resolver (cached)
# ------------------------------------------------------------
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
op_resolve <- make_omicspred_resolver(CFG$omicpred_map)

# ------------------------------------------------------------
# Read .sscore with fallback to zeros
# ------------------------------------------------------------
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

  score_col <- intersect(c("SCORE1_AVG", "SCORE1_SUM", "SCORE1", "score", "SCORE"), names(dt))[1]
  if (is.na(score_col)) score_col <- names(dt)[ncol(dt)]

  out <- dt[, .(eid = as.integer(get(id_col)), score = as.numeric(get(score_col)))]
  data.table::setnames(out, c("eid", out_col))
  out
}

extract_protGS_cis_trans <- function(protID, eids, op_resolve, strict = FALSE, cfg = CFG) {
  opID <- op_resolve(protID)
  p <- prot_clean(protID)

  cis_col <- paste0(p, "_GScis")
  tr_col  <- paste0(p, "_GStrans")

  cis_path <- file.path(cfg$gs_cis_dir, paste0(opID, ".sscore"))
  tr_path  <- file.path(cfg$gs_tr_dir,  paste0(opID, ".sscore"))

  cis_dt <- read_sscore_or_zero(cis_path, out_col = cis_col, eids = eids, strict = strict)
  tr_dt  <- read_sscore_or_zero(tr_path,  out_col = tr_col,  eids = eids, strict = strict)

  merge(cis_dt, tr_dt, by = "eid", all = TRUE)
}

GS_struct_cis_trans <- function(protID, UKBprot_df, op_resolve, strict = FALSE, cfg = CFG) {
  p <- prot_clean(protID)
  omic_orig <- UKBprot_df %>% dplyr::select(all_of(c("eid", p)))
  eids <- omic_orig$eid

  gs2 <- extract_protGS_cis_trans(protID, eids, op_resolve, strict = strict, cfg = cfg)

  list(
    combo = na.omit(merge(gs2, omic_orig, by = "eid")),
    solo  = na.omit(gs2)
  )
}

# ------------------------------------------------------------
# Partitioned R2 from ANOVA table (Type 1; order-dependent)
# ------------------------------------------------------------
calculate_partitioned_r2 <- function(model, anovaType = 1) {
  if (anovaType > 1) {
    anova_table <- car::Anova(model, type = anovaType)
  } else {
    anova_table <- anova(model)
  }

  anova_df <- as.data.frame(anova_table)

  ss_col <- intersect(c("Sum Sq", "Sum.Sq", "Sum Sq."), names(anova_df))[1]
  if (is.na(ss_col)) stop("Could not find sum-of-squares column in ANOVA table.")

  rn <- rownames(anova_df)
  ss_total <- sum(anova_df[[ss_col]], na.rm = TRUE)

  if ("Residuals" %in% rn) {
    ss_residual <- anova_df["Residuals", ss_col]
    ss_factors <- anova_df[rn != "Residuals", ss_col, drop = TRUE]
    fac_names <- rn[rn != "Residuals"]
  } else {
    ss_residual <- NA_real_
    ss_factors <- anova_df[[ss_col]]
    fac_names <- rn
  }

  r2_total <- if (is.finite(ss_residual)) 1 - ss_residual / ss_total else NA_real_
  r2_factors <- ss_factors / ss_total

  out <- data.frame(
    Factor = fac_names,
    Partitioned.R2 = as.numeric(r2_factors),
    stringsAsFactors = FALSE
  )
  out <- out %>% column_to_rownames("Factor")

  if (is.finite(r2_total)) {
    out <- rbind(out, data.frame(Partitioned.R2 = r2_total, row.names = "Total"))
  }

  out
}

# ------------------------------------------------------------
# 1) Load + merge + split
# ------------------------------------------------------------
create_data <- function(protID, PXSdata, split_ratio){
  omicDS <- GS_struct_cis_trans(protID, PXSdata$UKBprot_df, op_resolve = op_resolve, strict = FALSE)

  protID <- prot_clean(protID)
  E_df <- PXSdata$Elist %>% reduce(full_join, by = "eid")

  df <- merge(omicDS$combo, E_df, by = "eid")
  df <- merge(df, PXSdata$covars_df, by = "eid")

  E_ids <- setdiff(colnames(E_df), "eid")

  # Set factors/ordinals BEFORE split (consistent levels)
  df <- categorical_handler(df, PXSdata$ordinalIDs, ordinal_contrast = "treatment")

  trainIndex <- sample(seq_len(nrow(df)), size = split_ratio * nrow(df))
  list(
    train = df[trainIndex, ],
    test  = df[-trainIndex, ],
    Eids  = E_ids,
    ordinalVar = PXSdata$ordinalIDs
  )
}

# ------------------------------------------------------------
# 2) Scaling (omic + Gcis + Gtrans + covars + continuous E)
# ------------------------------------------------------------
preprocess_data <- function(protID, covariates,
                            trainData, testData,
                            ordinal_contrast = "treatment",
                            ordinal_names){

  protID <- prot_clean(protID)
  Gcis_var <- paste0(protID,"_GScis")
  Gtr_var  <- paste0(protID,"_GStrans")

  trainData <- as.data.frame(trainData)
  testData  <- as.data.frame(testData)

  rel_columns <- c(protID, Gcis_var, Gtr_var, covariates)
  rel_columns <- intersect(rel_columns, names(trainData))

  numeric_cols <- sapply(trainData[, rel_columns, drop = FALSE],
                         function(col) is.numeric(col) && !(all(col %in% c(0, 1))))
  DG_cont_names <- names(numeric_cols[numeric_cols])

  dfE <- trainData %>% dplyr::select(!all_of(intersect(c("eid", rel_columns, ordinal_names), names(trainData))))
  # Continuous EXPOSURES come from the declared variable_type (analysis_exposures.tsv,
  # via heap_exposures_of_type); fall back to the value-range heuristic only if the
  # config is unavailable. (binary exposures stay 0/1 and are not standardised.)
  .E_cont_cfg <- heap_exposures_of_type("continuous", names(dfE))
  E_cont_names <- if (is.null(.E_cont_cfg))
                    intersect(continuous_finder(dfE), names(trainData)) else .E_cont_cfg

  numeric_col_names <- unique(c(DG_cont_names, E_cont_names))

  mean_train <- lapply(trainData[numeric_col_names], function(x) mean(x, na.rm = TRUE))
  sd_train   <- lapply(trainData[numeric_col_names], function(x) sd(x, na.rm = TRUE))

  for (i in numeric_col_names) {
    s <- sd_train[[i]]
    if (!is.finite(s) || s == 0) {
      trainData[[i]] <- 0
      testData[[i]]  <- 0
    } else {
      trainData[[i]] <- (trainData[[i]] - mean_train[[i]]) / s
      testData[[i]]  <- (testData[[i]]  - mean_train[[i]]) / s
    }
  }

  list(train = trainData, test = testData)
}

# ------------------------------------------------------------
# helper: handle A:B vs B:A rowname differences
# ------------------------------------------------------------
match_term <- function(term, rn) {
  if (term %in% rn) return(term)
  if (grepl(":", term, fixed = TRUE)) {
    parts <- strsplit(term, ":", fixed = TRUE)[[1]]
    revt <- paste0(parts[2], ":", parts[1])
    if (revt %in% rn) return(revt)
  }
  NA_character_
}

# ------------------------------------------------------------
# helper: safely build cross-product interaction terms
# ------------------------------------------------------------
make_interaction_terms <- function(left_terms, right_terms) {
  if (length(left_terms) == 0 || length(right_terms) == 0) return(character(0))
  as.vector(outer(left_terms, right_terms, paste, sep = ":"))
}

# ------------------------------------------------------------
# 3) Univariate loop: per exposure E
#
# Full model:
#   P ~ Gcis + Gtrans + E + covars
#       + Gcis:E + Gtrans:E
#       [+ E:int_covariates + Gcis:int_covariates + Gtrans:int_covariates]
#
# Returns: list(statE, statGxE, statR2, statFblock)
# ------------------------------------------------------------
GxE_assoc_fin <- function(df, protID, covariates, Eids, int_covariates = NULL,
                          int_families = c("E", "G")){

  statE      <- list()
  statGxE    <- list()
  statR2     <- list()
  statFblock <- list()

  protID  <- gsub("-", "_", protID)
  res_var <- protID
  Gcis    <- paste0(protID, "_GScis")
  Gtr     <- paste0(protID, "_GStrans")

  fml <- function(rhs_terms) {
    rhs_terms <- rhs_terms[!is.na(rhs_terms) & nzchar(rhs_terms)]
    as.formula(paste(res_var, paste(rhs_terms, collapse = " + "), sep = " ~ "))
  }

  p_from_anova <- function(fit_reduced, fit_full) {
    if (is.null(fit_reduced) || is.null(fit_full)) return(NA_real_)
    out <- tryCatch(anova(fit_reduced, fit_full)$`Pr(>F)`[2], error = function(err) NA_real_)
    as.numeric(out)
  }

  int_covariates <- int_covariates %||% character(0)
  int_covariates <- intersect(int_covariates, names(df))

  # Drop zero-variance COVARIATES once for this split (constant -> aliased /
  # uninformative). `df` here is the train OR test split, so this also handles
  # covariates that become constant within a split.
  covariates <- setdiff(covariates, heap_zero_variance_cols(df, covariates))

  for (e in Eids) {

    # -----------------------------
    # Required columns
    # -----------------------------
    req <- unique(c(res_var, e, Gcis, Gtr, covariates, int_covariates))
    req <- intersect(req, names(df))
    if (!all(c(res_var, e, Gcis, Gtr) %in% req)) next

    df2 <- df %>% dplyr::select(all_of(req))

    # Skip zero-variance EXPOSURES: a constant E (globally constant, e.g.
    # former_alcohol_drinker_f3731, or constant within this split) has no
    # estimable main effect or interaction.
    if (length(heap_zero_variance_cols(df2, e))) next

    # -----------------------------
    # Build RHS terms
    # -----------------------------
    rhs_base <- c(e, Gcis, Gtr, covariates)

    GxE_terms     <- c(paste0(Gcis, ":", e), paste0(Gtr, ":", e))
    # Interaction-sensitivity covariate terms, gated by which families (E / G) are
    # requested: sens_ExC keeps only E x covar, sens_GxC keeps only G x covar, the
    # combined variant keeps both. Gating these vectors here propagates to rhs_full
    # AND every reduced model (they all reference the same vectors).
    ExCov_terms   <- if ("E" %in% int_families) make_interaction_terms(e,    int_covariates) else character(0)
    GcisCov_terms <- if ("G" %in% int_families) make_interaction_terms(Gcis, int_covariates) else character(0)
    GtrCov_terms  <- if ("G" %in% int_families) make_interaction_terms(Gtr,  int_covariates) else character(0)

    rhs_full <- c(rhs_base, ExCov_terms, GcisCov_terms, GtrCov_terms, GxE_terms)

    # -----------------------------
    # Fit FULL model
    # -----------------------------
    fit_full <- tryCatch(lm(fml(rhs_full), data = df2), error = function(err) NULL)
    if (is.null(fit_full)) next

    res <- summary(fit_full)

    # -----------------------------
    # Coeff table + R2/N
    # -----------------------------
    CI <- as.data.frame(res$coefficients)
    CI$R2 <- summary(fit_full)$r.squared
    CI$adj.R2 <- summary(fit_full)$adj.r.squared
    CI$samplesize <- length(fit_full$fitted.values)

    stat_tbl <- CI %>%
      rownames_to_column(var = "ID") %>%
      pivot_longer(cols = colnames(CI),
                   names_to = "stats",
                   values_to = "value")

    # -----------------------------
    # E terms (main effect only; factor-safe)
    # -----------------------------
    statE[[e]] <- stat_tbl %>%
      filter(ID == e | startsWith(ID, e) | grepl(e, ID, fixed = TRUE)) %>%
      filter(!grepl(":", ID, fixed = TRUE)) %>%
      filter(!grepl(Gcis, ID, fixed = TRUE)) %>%
      filter(!grepl(Gtr,  ID, fixed = TRUE))

    # -----------------------------
    # GxE terms (cis + trans), restricted to this exposure e
    # Parse interaction rownames to handle A:B vs B:A ordering
    # -----------------------------
    ids <- stat_tbl$ID
    is_int <- grepl(":", ids, fixed = TRUE)
    lhs <- sub(":.*$", "", ids)
    rhs <- sub("^.*:", "", ids)

    is_gcis_int <- is_int & (lhs == Gcis | rhs == Gcis)
    is_gtr_int  <- is_int & (lhs == Gtr  | rhs == Gtr)

    other <- ifelse(lhs %in% c(Gcis, Gtr), rhs, lhs)
    is_this_e <- (other == e) | startsWith(other, e)

    statGxE[[e]] <- stat_tbl %>%
      mutate(
        E_id = e,
        E_term = other,
        G_component = dplyr::case_when(
          is_gcis_int ~ "cis",
          is_gtr_int  ~ "trans",
          TRUE ~ NA_character_
        )
      ) %>%
      filter(is_int, is_this_e, is_gcis_int | is_gtr_int)

    # -----------------------------
    # Partitioned R2 (Type 1; internal use)
    # -----------------------------
    partR2 <- calculate_partitioned_r2(fit_full, anovaType = 1)
    rn <- rownames(partR2)

    cis_int <- paste0(Gcis, ":", e)
    tr_int  <- paste0(Gtr,  ":", e)

    e_key   <- match_term(e, rn)
    cis_key <- match_term(cis_int, rn)
    tr_key  <- match_term(tr_int, rn)

    statR2[[e]] <- data.frame(
      ID = c(e, cis_int, tr_int),
      R2 = c(
        if (!is.na(e_key))   partR2[e_key,   "Partitioned.R2"] else NA_real_,
        if (!is.na(cis_key)) partR2[cis_key, "Partitioned.R2"] else NA_real_,
        if (!is.na(tr_key))  partR2[tr_key,  "Partitioned.R2"] else NA_real_
      ),
      stringsAsFactors = FALSE
    )

    # -----------------------------
    # Nested-model (block) F-tests
    # All reduced models keep ExCov and GxCov terms (if present).
    # For standard runs (int_covariates = NULL) these degenerate to
    # the same tests as the original Module2 without sensitivity terms.
    # -----------------------------

    # p_E_block: drop E main effect; GxE + ExCov + GxCov remain
    rhs_noE <- unique(c(Gcis, Gtr, covariates,
                        ExCov_terms, GcisCov_terms, GtrCov_terms,
                        GxE_terms))
    fit_noE <- tryCatch(lm(fml(rhs_noE), data = df2), error = function(err) NULL)

    # p_GxE_joint: drop both GxE blocks; E main + ExCov + GxCov remain
    rhs_noGxE <- c(rhs_base, ExCov_terms, GcisCov_terms, GtrCov_terms)
    fit_noGxE  <- tryCatch(lm(fml(rhs_noGxE), data = df2), error = function(err) NULL)

    # p_GcisxE: drop cis interaction only
    fit_noCis <- tryCatch(
      lm(fml(c(rhs_noGxE, paste0(Gtr, ":", e))), data = df2),
      error = function(err) NULL
    )

    # p_GtrxE: drop trans interaction only
    fit_noTr <- tryCatch(
      lm(fml(c(rhs_noGxE, paste0(Gcis, ":", e))), data = df2),
      error = function(err) NULL
    )

    statFblock[[e]] <- data.frame(
      ID           = e,
      p_E_block    = p_from_anova(fit_noE,   fit_full),
      p_GxE_joint  = p_from_anova(fit_noGxE, fit_full),
      p_GcisxE     = p_from_anova(fit_noCis, fit_full),
      p_GtrxE      = p_from_anova(fit_noTr,  fit_full),
      R2           = summary(fit_full)$r.squared,
      adj.R2       = summary(fit_full)$adj.r.squared,
      samplesize   = nobs(fit_full),
      n_int_covars = length(int_covariates),
      stringsAsFactors = FALSE
    )
  }

  # Aggregate outputs
  statE <- if (length(statE)) {
    do.call(rbind, statE) %>% pivot_wider(names_from = stats, values_from = value)
  } else data.frame()

  statGxE <- if (length(statGxE)) {
    do.call(rbind, statGxE) %>% pivot_wider(names_from = stats, values_from = value)
  } else data.frame()

  statR2 <- if (length(statR2)) {
    out <- do.call(rbind, statR2)
    rownames(out) <- NULL
    out
  } else data.frame()

  statFblock <- if (length(statFblock)) do.call(rbind, statFblock) else data.frame()

  list(statE, statGxE, statR2, statFblock)
}

# ------------------------------------------------------------
# 4) Attach category labels to stats tables
# ------------------------------------------------------------
find_all_matches <- function(partial_id, unique_ids) {
  matches <- unique_ids[str_detect(unique_ids, fixed(partial_id))]
  if (length(matches) > 0) matches else NA
}

add_cat_info <- function(Eid_cat, stat_df){
  if (is.null(stat_df) || nrow(stat_df) == 0) return(stat_df)
  if (!("ID" %in% names(stat_df))) return(stat_df)

  df <- Eid_cat
  df$ID <- lapply(df$Eid, find_all_matches, unique_ids = stat_df$ID)

  df_long <- df %>% unnest(ID)
  out <- merge(df_long, stat_df, by = "ID")
  out
}

# ------------------------------------------------------------
# 5) Run per protein list
# ------------------------------------------------------------
runProt_univar <- function(protlist, PXSdata, int_covariates = NULL,
                           int_families = c("E", "G")){
  stat_batch_train <- NULL
  stat_batch_test  <- NULL

  for (prot in protlist) {
    message("Running: ", prot)

    protData <- create_data(protID = prot, PXSdata = PXSdata, split_ratio = 0.8)

    protPP <- preprocess_data(
      protID = prot,
      covariates = PXSdata$covars_list,
      trainData = protData$train,
      testData  = protData$test,
      ordinal_contrast = "treatment",
      ordinal_names = protData$ordinalVar
    )

    stat_train <- GxE_assoc_fin(
      df = protPP$train,
      protID = prot,
      covariates = PXSdata$covars_list,
      Eids = protData$Eids,
      int_covariates = int_covariates,
      int_families = int_families
    )

    stat_test <- GxE_assoc_fin(
      df = protPP$test,
      protID = prot,
      covariates = PXSdata$covars_list,
      Eids = protData$Eids,
      int_covariates = int_covariates,
      int_families = int_families
    )

    # Attach category info
    stat_train <- lapply(stat_train, function(x) add_cat_info(PXSdata$Eid_cat, x))
    stat_test  <- lapply(stat_test,  function(x) add_cat_info(PXSdata$Eid_cat, x))

    # Attach omic ID
    stat_train <- purrr::map(stat_train, ~ dplyr::mutate(.x, omicID = prot))
    stat_test  <- purrr::map(stat_test,  ~ dplyr::mutate(.x, omicID = prot))

    # Accumulate batch. `.x` is NULL on the first protein (slots initialised to
    # NULL), and dplyr::bind_rows(NULL, .y) == .y, so the first protein is kept
    # exactly once. (The previous `(.x %||% .y) %>% bind_rows(.y)` bound the
    # first protein's table to itself, double-counting the first protein of
    # every batch.)
    stat_batch_train <- purrr::map2(
      stat_batch_train %||% vector("list", length(stat_train)),
      stat_train,
      ~ dplyr::bind_rows(.x, .y)
    )

    stat_batch_test <- purrr::map2(
      stat_batch_test %||% vector("list", length(stat_test)),
      stat_test,
      ~ dplyr::bind_rows(.x, .y)
    )
  }

  list(train = stat_batch_train, test = stat_batch_test)
}

# ------------------------------------------------------------
# 6) Covariate spec utilities
# ------------------------------------------------------------
orgElistnames <- function(Elist_names){
  Elist_names <- gsub(" ","_", Elist_names)
  Elist_names <- gsub(paste(c("[(]", "[)]"), collapse = "|"),"",Elist_names)
  Elist_names
}

orgEids <- function(Elist){
  for(i in names(Elist)){
    Eids <- gsub(" ","_", colnames(Elist[[i]]))
    Eids <- make.names(Eids)
    colnames(Elist[[i]]) <- Eids
  }
  Elist
}

orgEtoCat <- function(Elist, Elist_names){
  Ecategory <- lapply(Elist, function(x) setdiff(colnames(x), "eid"))
  names(Ecategory) <- Elist_names
  Eid_cat <- stack(Ecategory)
  colnames(Eid_cat) <- c("Eid", "Category")
  Eid_cat$Category <- as.character(Eid_cat$Category)
  Eid_cat
}

covariate_to_E <- function(PXSnew, FeatIds, FeatCat){
  df <- PXSnew$covars_df %>% select(all_of(c("eid", FeatIds)))
  PXSnew$Elist[[FeatCat]] <- df
  PXSnew$Elist_names <- orgElistnames(names(PXSnew$Elist))
  PXSnew$Elist <- orgEids(PXSnew$Elist)
  PXSnew$Eid_cat <- orgEtoCat(PXSnew$Elist, PXSnew$Elist_names)

  PXSnew$covars_df <- PXSnew$covars_df %>% select(-all_of(c(FeatIds)))
  PXSnew$covars_list <- PXSnew$covars_list[!(PXSnew$covars_list %in% FeatIds)]
  PXSnew
}

E_to_Covariate <- function(PXSnew, Eids, Ecat){
  Edf <- PXSnew$Elist[[Ecat]] %>% select(all_of(c("eid", Eids)))
  PXSnew$covars_df <- merge(PXSnew$covars_df, Edf, by = "eid")
  PXSnew$covars_list <- c(PXSnew$covars_list, Eids)

  PXSnew$Elist[[Ecat]] <- NULL
  PXSnew$Elist_names <- PXSnew$Elist_names[!(PXSnew$Elist_names %in% Ecat)]
  PXSnew$Eid_cat <- orgEtoCat(PXSnew$Elist, PXSnew$Elist_names)
  PXSnew
}

PXScovarSpec <- function(PXSdata, covars_subset){
  PXSdata$covars_list <- covars_subset
  PXSdata$covars_df <- PXSdata$covars_df %>% select(all_of(c("eid", covars_subset)))
  PXSdata
}

# ------------------------------------------------------------
# 7) Covariate specs
# Type1-5: loaded from config/covariates/covariate_sets.yml (inline fallback).
# Type6/7: dynamically rearrange exposures as covariates — must remain inline
#          because they modify the PXS object at runtime, not just the covariate list.
# ------------------------------------------------------------

.resolve_covartype_m2 <- function(name, pxsloader) {
  .cfg_path <- heap_config("covariates", "covariate_sets.yml")
  if (file.exists(.cfg_path) && requireNamespace("yaml", quietly = TRUE)) {
    tryCatch({
      .sets <- yaml::read_yaml(.cfg_path)$covariate_sets
      if (name %in% names(.sets)) {
        .entry  <- .sets[[name]]
        .covars <- .entry$covariates
        if (is.null(.covars) || (length(.covars) == 1L && is.na(.covars[[1L]])))
          return(pxsloader$covars_list)
        return(as.character(.covars))
      }
    }, error = function(e) NULL)
  }
  # Inline fallback (used only if covariate_sets.yml is unreadable). Mirrors the
  # descriptive sets in covariate_sets.yml. base_ses resolves to the base list here;
  # its deprivation variables are added at runtime via E_to_Covariate.
  .pcs  <- paste0("genetic_principal_components_f22009_0_", 1:20)
  .base <- c("age_when_attended_assessment_centre_f21003_0_0", "sex_f31_0_0",
             "age2", "age_sex", "age2_sex",
             "uk_biobank_assessment_centre_f54_0_0", .pcs)
  switch(name,
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
    base_prevalent = c(.base, "prevalent_major_disease"),
    NULL
  )
}

CovarSpec <- list()
for (.ct in c("base", "base_bmi", "base_draw", "base_clinical",
              "base_ses", "base_prevalent")) {
  CovarSpec[[.ct]] <- .resolve_covartype_m2(.ct, PXSloader)
}
rm(.ct)

# ------------------------------------------------------------
# 7b) Sensitivity-analysis interaction covariates
# Only these covariates are allowed to interact with E and G.
# ------------------------------------------------------------
SensIntCovars <- list()

SensIntCovars[["core"]] <- c(
  "age_when_attended_assessment_centre_f21003_0_0",
  "sex_f31_0_0",
  "body_mass_index_bmi_f23104_0_0",
  "fasting_time_f74_0_0"
)

SensIntCovars[["extended"]] <- CovarSpec[["base"]]

# {age, sex} only -- mirrors the Module-1 interaction-term sensitivity (FigS13)
# so the GxE F-test "is it an age/sex-interaction artifact?" can be asked at the
# association level, separated into G x {age,sex} (sens_GxC) and E x {age,sex}
# (sens_ExC) and both (sens_GxC_ExC).
SensIntCovars[["agesex"]] <- c(
  "age_when_attended_assessment_centre_f21003_0_0",
  "sex_f31_0_0"
)

# ------------------------------------------------------------
# 8) Batch runner
# ------------------------------------------------------------
omiclist <- scan(
  file  = heap_omicspred_protein_list,
  what  = character(),
  quiet = TRUE
)

# ------------------------------------------------------------
# 8a) CLI: manifest-driven or legacy positional arguments
#
# Manifest mode (preferred):
#   Rscript Module2.R --manifest <path> --array-index <N>
#
# Legacy positional mode (backward compatibility):
#   Rscript Module2.R <idx> <split_num> <covarType>
# ------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

.parse_flag_m2 <- function(args, flag, default = NULL) {
  i <- which(args == flag)
  if (length(i) == 0 || i[1] >= length(args)) return(default)
  args[i[1] + 1L]
}

.is_manifest_mode_m2 <- length(args) >= 2 && args[1] == "--manifest"

if (.is_manifest_mode_m2) {
  .manifest_path <- .parse_flag_m2(args, "--manifest")
  .array_idx_str <- .parse_flag_m2(args, "--array-index")

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
  experiment_name <- as.character(.row$experiment_name %||% "unknown")
  .covar_variant  <- as.character(.row$covar_variant %||% "")
  .sample_filter  <- as.character(.row$sample_filter %||% "none")
  .manifest_out_root <- if (!is.null(.row$output_path) && !is.na(.row$output_path) &&
                             nzchar(.row$output_path))
    as.character(.row$output_path) else NULL

  message(sprintf("[manifest] experiment=%s  array_index=%d  chunk=%d/%d  covar=%s",
                  experiment_name, .arr_idx, idx, split_num, covarType))

} else {
  if (length(args) < 3) stop(
    "Usage (manifest):   Rscript Module2.R --manifest <path> --array-index <N>\n",
    "Usage (positional): Rscript Module2.R <idx> <split_num> <covarType>\n",
    "  covarType: base | base_bmi | base_draw | base_clinical | base_ses | base_prevalent\n",
    "  (interaction sensitivities and sample filters require manifest mode)"
  )
  idx            <- as.integer(args[1])
  split_num      <- as.integer(args[2])
  covarType      <- as.character(args[3])
  experiment_name <- paste0("positional_", covarType)
  .covar_variant  <- ""
  .sample_filter  <- "none"
  .manifest_out_root <- NULL
}

groups        <- cut(seq_along(omiclist), breaks = split_num, labels = FALSE)
split_vectors <- split(omiclist, groups)

# Resolve output root: manifest takes priority; fall back to CFG defaults
.m2_out_root <- if (!is.null(.manifest_out_root)) {
  .manifest_out_root
} else if (startsWith(.covar_variant, "sens")) {
  CFG$out_root_sens
} else {
  CFG$out_root
}
if (grepl(HEAP_PATHS$scratch_root, normalizePath(.m2_out_root, mustWork = FALSE), fixed = TRUE))
  stop("REPRODUCIBILITY VIOLATION: Module2 out_root points to scratch: ", .m2_out_root)

# Write run config before computation begins
.m2_job_dir <- file.path(.m2_out_root, covarType)
dir.create(.m2_job_dir, recursive = TRUE, showWarnings = FALSE)
tryCatch({
  .rc <- list(
    module          = "module2",
    experiment_name = experiment_name,
    chunk_id        = idx,
    n_chunks        = split_num,
    covariate_set   = covarType,
    covar_variant   = .covar_variant,
    out_root        = .m2_out_root,
    heap_rds        = heap_loader_rds,
    manifest_path   = if (.is_manifest_mode_m2) .manifest_path else NA_character_,
    array_index     = if (.is_manifest_mode_m2) .arr_idx else NA_integer_,
    run_timestamp   = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    run_host        = Sys.info()[["nodename"]]
  )
  if (requireNamespace("yaml", quietly = TRUE)) {
    yaml::write_yaml(.rc, file.path(.m2_job_dir, paste0("run_config_", idx, ".yml")))
  } else {
    saveRDS(.rc, file.path(.m2_job_dir, paste0("run_config_", idx, ".rds")))
  }
}, error = function(e) warning("Could not write run_config: ", conditionMessage(e)))

# ------------------------------------------------------------
# 8c) Sample filter (the WHO-IS-IN axis), applied to the FULL covariate frame
# before any covariate-column narrowing. The model frame is built by an INNER join
# on covars_df (see runProt: merge(df, PXSdata$covars_df, by = "eid")), so restricting
# covars_df rows here restricts the analysis sample for every covariate set / run.
# ------------------------------------------------------------
.sfilter_spec <- load_sample_filter(.sample_filter)
if (!is.null(.sfilter_spec)) {
  .sf <- apply_sample_filter(PXSloader$covars_df, .sfilter_spec)
  PXSloader$covars_df <- .sf$df
  message(sprintf(
    "[sample_filter] %s: covars %d -> %d participants (dropped %d); inner-joins propagate to the analysis sample.",
    .sfilter_spec$name, .sf$n_before, .sf$n_after, .sf$n_dropped))
}

runProt_univar_assoc <- function(protlist, PXSdata, folder_id, idx,
                                 out_root = .m2_out_root,
                                 int_covariates = NULL,
                                 int_families = c("E", "G")){
  out_dir <- file.path(out_root, folder_id)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  batch_run <- runProt_univar(
    protlist       = protlist,
    PXSdata        = PXSdata,
    int_covariates = int_covariates,
    int_families   = int_families
  )

  saveRDS(batch_run, file = file.path(out_dir, paste0("univar_assoc_", idx, ".rds")))
}

# ------------------------------------------------------------
# 9) Entrypoint
# ------------------------------------------------------------
.STANDARD_SETS <- c("base", "base_bmi", "base_draw", "base_clinical", "base_prevalent")

if (covarType %in% .STANDARD_SETS) {

  if (startsWith(.covar_variant, "sens")) {
    # Interaction (effect-modification) sensitivity: layer E x covar and/or G x covar
    # terms on the chosen covariate set. Each variant maps to (interaction covariates,
    # which families interact). sens_GxC / sens_ExC / sens_GxC_ExC use {age,sex} and
    # mirror the Module-1 interaction-term sensitivity (FigS13) at the association level.
    .sens_map <- list(
      sens_core      = list(cov = "core",     fam = c("E", "G")),
      sens_extended  = list(cov = "extended", fam = c("E", "G")),
      sens_GxC       = list(cov = "agesex",   fam = c("G")),
      sens_ExC       = list(cov = "agesex",   fam = c("E")),
      sens_GxC_ExC   = list(cov = "agesex",   fam = c("E", "G"))
    )
    .sm <- .sens_map[[.covar_variant]]
    if (is.null(.sm)) stop("Unknown covar_variant: ", .covar_variant)
    runProt_univar_assoc(
      split_vectors[[idx]],
      PXSdata        = PXScovarSpec(PXSloader, CovarSpec[[covarType]]),
      folder_id      = covarType,
      idx            = idx,
      out_root       = .m2_out_root,
      int_covariates = SensIntCovars[[.sm$cov]],
      int_families   = .sm$fam
    )
  } else {
    runProt_univar_assoc(
      split_vectors[[idx]],
      PXSdata   = PXScovarSpec(PXSloader, CovarSpec[[covarType]]),
      folder_id = covarType,
      idx       = idx
    )
  }

} else if (covarType == "base_ses") {

  # Socioeconomic deprivation (Module2 only): move the configured deprivation
  # variables from the Deprivation_Indices exposure into the covariate matrix
  # (E_to_Covariate), then run. The variable list + source category come from the
  # set's `remapping` block in covariate_sets.yml (England scores + household income).
  .remap <- load_covariate_set_remapping("base_ses")
  if (is.null(.remap) || is.null(.remap$variables))
    stop("base_ses requires a `remapping` block (variables, source_category) in covariate_sets.yml")
  PXS_SES <- E_to_Covariate(
    PXSnew = PXScovarSpec(PXSloader, CovarSpec[["base_ses"]]),
    Eids   = as.character(unlist(.remap$variables)),
    Ecat   = as.character(.remap$source_category)
  )
  runProt_univar_assoc(
    split_vectors[[idx]],
    PXSdata   = PXS_SES,
    folder_id = covarType,
    idx       = idx
  )

} else {
  stop("Unknown covarType: ", covarType,
       " (expected one of: ", paste(c(.STANDARD_SETS, "base_ses"), collapse = ", "), ")")
}
