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
options(warn = 1)

script_file <- grep("^--file=", commandArgs(), value = TRUE)
script_dir <- if (length(script_file) == 0L) getwd() else dirname(normalizePath(sub("^--file=", "", script_file[1L])))
source(file.path(script_dir, "common.R"))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3L) {
  stop(
    paste(
      "Usage: export_architecture_inputs.R <config.R> <run_id> <covar_spec>",
      "[--proteins=A,B] [--max-samples=N] [--seed=1]",
      "[--complete-case-exposures=true|false] [--protein-specific-prep=true|false]",
      "[--force=true]"
    ),
    call. = FALSE
  )
}

config_path <- args[1L]
run_id <- args[2L]
covar_spec_name <- args[3L]
opts <- parse_optional_args(args[-(1L:3L)])

cfg <- load_config(config_path)
if (!covar_spec_name %in% names(cfg$covariate_specs)) {
  stopf("Unknown covariate specification: %s", covar_spec_name)
}

seed <- as.integer(get_opt(opts, "seed", 1L))
max_samples <- as.integer(get_opt(opts, "max_samples", NA_integer_))
complete_case_exposures <- as_bool(get_opt(opts, "complete_case_exposures", FALSE), default = FALSE)
protein_specific_prep <- as_bool(get_opt(opts, "protein_specific_prep", FALSE), default = FALSE)
strict_scores <- as_bool(get_opt(opts, "strict_scores", FALSE), default = FALSE)
force <- as_bool(get_opt(opts, "force", FALSE))
paths <- resolve_run_paths(cfg, run_id, covar_spec_name, exposure_mode = get_opt(opts, "exposure_mode", "centered"))

sample_manifest_path <- file.path(paths$inputs, "sample_manifest.tsv")
if (file.exists(sample_manifest_path) && !force) {
  stopf("Input export already exists at %s. Re-run with --force=true to overwrite.", sample_manifest_path)
}

timestamp_msg("Loading staged PXS loader from", cfg$loader_rds)
pxs_loader <- normalize_loader(read_pxs_loader(cfg$loader_rds))
pxs <- as_pxs(pxs_loader)
proteins <- resolve_proteins(
  pxs_loader,
  protein_arg = get_opt(opts, "proteins", NULL),
  default_proteins = cfg$pilot_proteins %||% NULL
)
timestamp_msg("Selected proteins:", paste(proteins, collapse = ", "))
if (protein_specific_prep && length(proteins) != 1L) {
  stopf("Protein-specific prep requires exactly one protein, got %s.", length(proteins))
}

covar_spec <- cfg$covariate_specs[[covar_spec_name]]
all_covars <- unique(c(covar_spec$discrete, covar_spec$quantitative))

geno_ids <- read_psam(cfg$genotype_pfile)
geno_iid <- geno_ids$IID
prot_df <- as.data.frame(pxs$UKBprot_df[, c("eid", proteins), drop = FALSE])
covar_df <- as.data.frame(pxs$covars_df)
covar_df <- covar_df[, c("eid", all_covars), drop = FALSE]
covar_df$sex_numeric <- ifelse(covar_df$sex_f31_0_0 == "Male", 1, 0)

analysis_metrics <- list()
if (protein_specific_prep) {
  protein_setup <- build_protein_specific_analysis(
    protein_id = proteins[1L],
    pxs_loader = pxs_loader,
    covariate_names = all_covars,
    genotype_ids = geno_iid,
    cfg = cfg,
    strict_scores = strict_scores,
    use_scores = FALSE
  )
  prot_df <- protein_setup$protein_df
  covar_df <- protein_setup$covar_df
  covar_df$sex_numeric <- ifelse(covar_df$sex_f31_0_0 == "Male", 1, 0)
  analysis_ids <- protein_setup$analysis_ids
  analysis_metrics <- list(
    target_protein = proteins[1L],
    protein_nonmissing = as.integer(protein_setup$counts[["protein_nonmissing"]]),
    pre_covars_overlap = as.integer(protein_setup$counts[["pre_covars_overlap"]]),
    after_covars = as.integer(protein_setup$counts[["after_covars"]]),
    after_genotype = as.integer(protein_setup$counts[["after_genotype"]])
  )
  timestamp_msg(
    "Protein-specific prep for",
    proteins[1L],
    "retained",
    length(analysis_ids),
    "samples after phenotype/covariate/genotype overlap."
  )
} else {
  base_ids <- intersect(unique(prot_df$eid), geno_iid)
  base_ids <- intersect(base_ids, unique(covar_df$eid))
  covar_complete <- stats::complete.cases(covar_df[covar_df$eid %in% base_ids, all_covars, drop = FALSE])
  covar_keep <- covar_df$eid[covar_df$eid %in% base_ids][covar_complete]
  analysis_ids <- geno_iid[geno_iid %in% covar_keep]
  analysis_metrics <- list(
    target_protein = "",
    protein_nonmissing = "",
    pre_covars_overlap = "",
    after_covars = length(covar_keep),
    after_genotype = length(analysis_ids)
  )
}

if (length(analysis_ids) == 0L) {
  stopf("No overlapping analysis samples remained for %s.", covar_spec_name)
}

build_exposure_feature_audit <- function(exposure_df, pxs_obj, missing_threshold,
                                         selected_features = character(),
                                         prefix = "pre") {
  feature_names <- setdiff(names(exposure_df), "eid")
  if (length(feature_names) == 0L) {
    return(data.frame())
  }

  feature_frame <- exposure_df[, feature_names, drop = FALSE]
  nonmissing_n <- colSums(!is.na(feature_frame))
  missing_rate <- colMeans(is.na(feature_frame))
  category_lookup <- setNames(as.character(pxs_obj$Eid_cat$Category), pxs_obj$Eid_cat$Eid)

  data.frame(
    feature = feature_names,
    category = unname(category_lookup[feature_names]),
    nonmissing_n = as.integer(nonmissing_n),
    missing_rate = as.numeric(missing_rate),
    keep_under_missing_threshold = if (is.finite(missing_threshold)) missing_rate <= missing_threshold else NA,
    selected_for_export = feature_names %in% selected_features,
    stage = prefix,
    stringsAsFactors = FALSE
  )
}

summarize_exposure_categories <- function(feature_audit, ordered_categories) {
  if (nrow(feature_audit) == 0L) {
    return(data.frame())
  }

  out <- do.call(
    rbind,
    lapply(ordered_categories, function(cat_name) {
      idx <- feature_audit$category == cat_name
      sub <- feature_audit[idx, , drop = FALSE]
      if (nrow(sub) == 0L) {
        return(NULL)
      }
      data.frame(
        category = cat_name,
        n_features = nrow(sub),
        n_with_any_data = sum(sub$nonmissing_n > 0L),
        n_keep_under_missing_threshold = sum(sub$keep_under_missing_threshold %in% TRUE),
        n_selected_for_export = sum(sub$selected_for_export),
        median_missing_rate = stats::median(sub$missing_rate),
        stringsAsFactors = FALSE
      )
    })
  )
  if (is.null(out)) {
    data.frame()
  } else {
    out
  }
}

analysis_ids_pre_exposure <- analysis_ids
timestamp_msg("Merging exposures for", length(analysis_ids_pre_exposure), "analysis samples.")
exposure_raw <- merge_exposure_list(pxs$Elist, analysis_ids_pre_exposure)

# Align the environment kernel to the curated analysis exposure set
# (analysis_exposures.tsv, include == 1) -- the SAME exposures used by modules
# 1/2/3/6 and the exposure GWAS (single source of truth) -- rather than the full
# loader Elist. Intersect on exposure variable name; keep eid. Falls back to the
# full Elist (with a warning) only if the curated config is unreadable.
.curated_exposure_path <- tryCatch(heap_analysis_config(), error = function(e)
  tryCatch(heap_config("exposure_sets", "analysis_exposures.tsv"),
           error = function(e2) NA_character_))
.included_vars <- if (!is.na(.curated_exposure_path) && file.exists(.curated_exposure_path)) {
  .ec <- utils::read.delim(.curated_exposure_path, stringsAsFactors = FALSE, check.names = FALSE)
  .ec$variable[suppressWarnings(as.integer(.ec$include)) == 1L]
} else NULL
if (!is.null(.included_vars) && length(.included_vars) > 0L) {
  .n_before <- ncol(exposure_raw) - 1L
  exposure_raw <- exposure_raw[, c("eid", intersect(names(exposure_raw), .included_vars)), drop = FALSE]
  timestamp_msg(
    "Environment kernel restricted to curated analysis_exposures.tsv (include==1):",
    ncol(exposure_raw) - 1L, "of", .n_before, "exposure columns kept."
  )
} else {
  timestamp_msg("WARNING: curated analysis_exposures.tsv unavailable; environment kernel uses the full loader Elist.")
}

exposure_aligned_pre <- align_to_ids(exposure_raw, analysis_ids_pre_exposure)
exposure_features_pre <- exposure_aligned_pre[, setdiff(names(exposure_aligned_pre), "eid"), drop = FALSE]
exposure_missing_rate_max <- suppressWarnings(as.numeric(cfg$exposure_missing_rate_max %||% NA_real_))

n_exposures_complete_case <- NA_integer_
n_samples_after_exposure_complete <- NA_integer_
exposure_cols_kept <- setdiff(names(exposure_raw), "eid")
exposure_cols_dropped <- character()
if (complete_case_exposures) {
  if (is.finite(exposure_missing_rate_max)) {
    dropped_exposure_cols <- missing_cols(exposure_features_pre, miss_rate = exposure_missing_rate_max)
    if (length(dropped_exposure_cols) > 0L) {
      timestamp_msg(
        "Dropping",
        length(dropped_exposure_cols),
        "exposure columns with missingness >",
        exposure_missing_rate_max,
        "before complete-case sampling."
      )
      exposure_cols_dropped <- dropped_exposure_cols
      timestamp_msg(
        "First dropped exposure columns:",
        paste(utils::head(exposure_cols_dropped, 15L), collapse = ", ")
      )
    }
  }

  exposure_cols_kept <- setdiff(colnames(exposure_features_pre), exposure_cols_dropped)
  n_exposures_complete_case <- length(exposure_cols_kept)
  if (n_exposures_complete_case == 0L) {
    stopf("No exposure columns remained after complete-case pruning for %s.", run_id)
  }

  exposure_complete <- stats::complete.cases(exposure_features_pre[, exposure_cols_kept, drop = FALSE])
  analysis_ids <- analysis_ids_pre_exposure[exposure_complete]
  n_samples_after_exposure_complete <- length(analysis_ids)
  timestamp_msg(
    "Retained",
    n_samples_after_exposure_complete,
    "analysis samples after exposure complete-case filtering on",
    n_exposures_complete_case,
    "retained exposure columns."
  )
  timestamp_msg(
    "Exposure complete-case summary:",
    "total =", ncol(exposure_features_pre),
    "; kept =", n_exposures_complete_case,
    "; dropped =", length(exposure_cols_dropped)
  )
}

if (is.finite(max_samples) && !is.na(max_samples) && max_samples > 0L && length(analysis_ids) > max_samples) {
  set.seed(seed)
  sampled <- sort(sample(analysis_ids, size = max_samples, replace = FALSE))
  analysis_ids <- analysis_ids[analysis_ids %in% sampled]
}

if (length(analysis_ids) == 0L) {
  stopf("No analysis samples remained after applying the requested filters for %s.", covar_spec_name)
}

exposure_aligned_final <- align_to_ids(exposure_raw, analysis_ids)
exposure_features_final <- exposure_aligned_final[, intersect(exposure_cols_kept, names(exposure_aligned_final)), drop = FALSE]
feature_audit_pre <- build_exposure_feature_audit(
  exposure_df = exposure_aligned_pre,
  pxs_obj = pxs,
  missing_threshold = exposure_missing_rate_max,
  selected_features = exposure_cols_kept,
  prefix = "pre_complete_case"
)
feature_audit_final <- build_exposure_feature_audit(
  exposure_df = data.frame(eid = analysis_ids, exposure_features_final, check.names = FALSE),
  pxs_obj = pxs,
  missing_threshold = exposure_missing_rate_max,
  selected_features = exposure_cols_kept,
  prefix = "final_export"
)
feature_audit <- merge(
  feature_audit_pre[, c("feature", "category", "nonmissing_n", "missing_rate", "keep_under_missing_threshold", "selected_for_export"), drop = FALSE],
  feature_audit_final[, c("feature", "nonmissing_n", "missing_rate"), drop = FALSE],
  by = "feature",
  all = TRUE,
  suffixes = c("_pre", "_final"),
  sort = FALSE
)
if (!"category" %in% names(feature_audit)) {
  feature_audit$category <- NA_character_
}
feature_audit <- feature_audit[, c(
  "feature",
  "category",
  "nonmissing_n_pre",
  "missing_rate_pre",
  "keep_under_missing_threshold",
  "selected_for_export",
  "nonmissing_n_final",
  "missing_rate_final"
), drop = FALSE]
category_audit <- summarize_exposure_categories(feature_audit_pre, pxs$Elist_names)

sample_manifest <- data.frame(
  FID = analysis_ids,
  IID = analysis_ids,
  eid = analysis_ids,
  stringsAsFactors = FALSE
)

covar_aligned <- align_to_ids(covar_df, analysis_ids)
disc_names <- unique(covar_spec$discrete)
quant_names <- unique(covar_spec$quantitative)
disc_keep <- intersect(disc_names, names(covar_aligned))
quant_keep <- intersect(quant_names, names(covar_aligned))

covar_discrete <- data.frame(FID = analysis_ids, IID = analysis_ids, covar_aligned[, disc_keep, drop = FALSE], check.names = FALSE)
covar_quant <- data.frame(FID = analysis_ids, IID = analysis_ids, covar_aligned[, quant_keep, drop = FALSE], check.names = FALSE)
core_covars <- unique(covar_spec$kernel_quantitative)
core_covars <- core_covars[core_covars %in% names(covar_aligned)]
covar_core <- data.frame(FID = analysis_ids, IID = analysis_ids, covar_aligned[, core_covars, drop = FALSE], check.names = FALSE)

exposure_export <- data.frame(
  FID = analysis_ids,
  IID = analysis_ids,
  exposure_features_final,
  check.names = FALSE
)

prot_aligned <- align_to_ids(prot_df, analysis_ids)
protein_export <- data.frame(FID = analysis_ids, IID = analysis_ids, prot_aligned[, proteins, drop = FALSE], check.names = FALSE)

write_tsv(sample_manifest, sample_manifest_path)
write_tsv(data.frame(FID = analysis_ids, IID = analysis_ids, stringsAsFactors = FALSE), file.path(paths$inputs, "keep_ids.txt"), header = FALSE)
write_tsv(covar_discrete, file.path(paths$inputs, "covar_discrete.tsv"))
write_tsv(covar_quant, file.path(paths$inputs, "covar_quantitative.tsv"))
write_tsv(covar_core, file.path(paths$inputs, "covar_kernel_core.tsv"))
write_tsv(exposure_export, file.path(paths$inputs, "exposure_raw.tsv"))
write_tsv(feature_audit, file.path(paths$inputs, "exposure_feature_audit.tsv"))
write_tsv(category_audit, file.path(paths$inputs, "exposure_category_audit.tsv"))
write_tsv(protein_export, file.path(paths$inputs, "protein_matrix.tsv"))
write_tsv(
  data.frame(protein = proteins, stringsAsFactors = FALSE),
  file.path(paths$inputs, "protein_list.tsv")
)
write_tsv(
  rbind(
    data.frame(
      metric = c(
        "n_samples",
        "n_proteins",
        "n_exposures_raw",
        "n_covars_quant",
        "n_covars_discrete",
        "protein_specific_prep",
        "target_protein",
        "n_samples_before_exposure_complete",
        "n_protein_nonmissing",
        "n_overlap_pre_covars",
        "n_after_covars",
        "n_after_genotype_overlap"
      ),
      value = c(
        length(analysis_ids),
        length(proteins),
        ncol(exposure_export) - 2L,
        ncol(covar_quant) - 2L,
        ncol(covar_discrete) - 2L,
        as.integer(protein_specific_prep),
        analysis_metrics$target_protein,
        length(analysis_ids_pre_exposure),
        analysis_metrics$protein_nonmissing,
        analysis_metrics$pre_covars_overlap,
        analysis_metrics$after_covars,
        analysis_metrics$after_genotype
      ),
      stringsAsFactors = FALSE
    ),
    data.frame(
      metric = c(
        "complete_case_exposures",
        "n_samples_after_exposure_complete",
        "n_exposures_retained_complete_case",
        "n_exposures_dropped_pre_complete_case"
      ),
      value = c(
        as.integer(complete_case_exposures),
        if (is.na(n_samples_after_exposure_complete)) "" else as.character(n_samples_after_exposure_complete),
        if (is.na(n_exposures_complete_case)) "" else as.character(n_exposures_complete_case),
        if (length(exposure_cols_dropped) == 0L) "0" else as.character(length(exposure_cols_dropped))
      ),
      stringsAsFactors = FALSE
    )
  ),
  file.path(paths$inputs, "metadata.tsv")
)

if (nrow(category_audit) > 0L) {
  kept_strings <- paste0(
    category_audit$category,
    "=",
    category_audit$n_selected_for_export,
    "/",
    category_audit$n_features
  )
  timestamp_msg(
    "Exposure category retained counts:",
    paste(kept_strings, collapse = "; ")
  )
}
timestamp_msg(
  "Exposure audit files:",
  file.path(paths$inputs, "exposure_feature_audit.tsv"),
  "and",
  file.path(paths$inputs, "exposure_category_audit.tsv")
)

timestamp_msg("Export complete:", sample_manifest_path)
