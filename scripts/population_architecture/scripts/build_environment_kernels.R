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

script_file <- grep("^--file=", commandArgs(), value = TRUE)
script_dir <- if (length(script_file) == 0L) getwd() else dirname(normalizePath(sub("^--file=", "", script_file[1L])))
source(file.path(script_dir, "common.R"))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3L) {
  stop("Usage: build_environment_kernels.R <config.R> <run_id> <covar_spec> [--force=true] [--include-covar-kernels=true] [--center-exposures=true|false]", call. = FALSE)
}

config_path <- args[1L]
run_id <- args[2L]
covar_spec_name <- args[3L]
opts <- parse_optional_args(args[-(1L:3L)])

cfg <- load_config(config_path)
force <- as_bool(get_opt(opts, "force", FALSE))
include_covar_kernels <- as_bool(get_opt(opts, "include_covar_kernels", TRUE))
center_exposures <- as_bool(get_opt(opts, "center_exposures", TRUE), default = TRUE)
exposure_mode <- if (center_exposures) "centered" else "uncentered"
paths <- resolve_run_paths(cfg, run_id, covar_spec_name, exposure_mode = exposure_mode)

geno_prefix <- file.path(paths$kernels, "geno_ld_pruned")
if (!file.exists(paste0(geno_prefix, ".grm.bin"))) {
  stopf("Missing genotype GRM: %s. Run build_genotype_grm.R first.", paste0(geno_prefix, ".grm.bin"))
}

ids <- read_grm_ids(geno_prefix)
exposure_raw <- read_tsv(file.path(paths$inputs, "exposure_raw.tsv"), header = TRUE)
sample_manifest <- read_tsv(file.path(paths$inputs, "sample_manifest.tsv"), header = TRUE, col_classes = "character")
covar_core <- read_tsv(file.path(paths$inputs, "covar_kernel_core.tsv"), header = TRUE)
metadata_path <- file.path(paths$inputs, "metadata.tsv")
metadata_df <- if (file.exists(metadata_path)) read_tsv(metadata_path, header = TRUE) else NULL
pxs_loader <- normalize_loader(read_pxs_loader(cfg$loader_rds))

metadata_value <- function(metric_name, default = NA_character_) {
  if (is.null(metadata_df) || !"metric" %in% names(metadata_df) || !"value" %in% names(metadata_df)) {
    return(default)
  }
  idx <- which(metadata_df$metric == metric_name)
  if (length(idx) == 0L) {
    return(default)
  }
  as.character(metadata_df$value[idx[1L]])
}

complete_case_requested <- identical(metadata_value("complete_case_exposures", "0"), "1")
exposure_missing_strategy <- get_opt(
  opts,
  "exposure_missing_strategy",
  if (complete_case_requested) "complete_case" else "impute"
)
if (!exposure_missing_strategy %in% c("impute", "complete_case")) {
  stopf("Unsupported exposure_missing_strategy: %s", exposure_missing_strategy)
}

exposure_imputation <- if (identical(exposure_missing_strategy, "complete_case")) {
  "none"
} else {
  cfg$exposure_imputation %||% "mean"
}
categorical_missing_mode <- if (identical(exposure_missing_strategy, "complete_case")) {
  "error"
} else {
  "missing_level"
}

align_export <- function(df, ids_df) {
  merged <- merge(ids_df, df, by = c("FID", "IID"), all.x = TRUE, sort = FALSE)
  merged
}

exposure_aligned <- align_export(exposure_raw, ids)
exposure_features <- exposure_aligned[, -(1L:2L), drop = FALSE]
exposure_missing_rate_max <- suppressWarnings(as.numeric(cfg$exposure_missing_rate_max %||% NA_real_))
if (is.finite(exposure_missing_rate_max)) {
  dropped_exposure_cols <- missing_cols(exposure_features, miss_rate = exposure_missing_rate_max)
  if (length(dropped_exposure_cols) > 0L) {
    timestamp_msg(
      "Dropping",
      length(dropped_exposure_cols),
      "exposure columns with missingness >",
      exposure_missing_rate_max
    )
    exposure_features <- drop_missing_cols(exposure_features, miss_rate = exposure_missing_rate_max)
  }
}
# Encode ordinals by the DECLARED variable_type (analysis_exposures.tsv), NOT the
# loader's max<=5 heuristic ordinalIDs. The heuristic mis-flags small-scale
# CONTINUOUS scores (income/employment/crime/health IMD scores, pm2.5) as ordinal,
# and make_numeric_matrix one-hot encodes every ordinal id -> crime_score alone
# would add ~481 dummy columns, inflating the feature count (M) and distorting the
# E / GxE GRMs. Config-driven: continuous -> one numeric column; declared ordinals
# -> one-hot (per level); binary -> one numeric 0/1 column. Falls back to the
# loader heuristic only if the config is unavailable. (heap_exposures_of_type and
# 00_paths.R are sourced above.)
.ord_cfg <- if (exists("heap_exposures_of_type", mode = "function"))
              heap_exposures_of_type("ordinal", names(exposure_features)) else NULL
ordinal_ids_use <- if (is.null(.ord_cfg)) pxs_loader$ordinalIDs else .ord_cfg
timestamp_msg(sprintf(
  "Ordinal encoding ids: %d (declared ordinal) vs %d (loader heuristic ordinalIDs)",
  length(ordinal_ids_use), length(intersect(pxs_loader$ordinalIDs, names(exposure_features)))))
exposure_encoded <- make_numeric_matrix(
  exposure_features,
  ordinal_ids = ordinal_ids_use,
  known_factor_vars = cfg$known_factor_vars %||% character(),
  impute_method = exposure_imputation,
  center = center_exposures,
  scale_columns = TRUE,
  categorical_missing = categorical_missing_mode
)

e_prefix <- file.path(paths$kernels, paste0("E_", exposure_mode))
if (!file.exists(paste0(e_prefix, ".grm.bin")) || force) {
  timestamp_msg("Building exposure GRM:", e_prefix)
  write_gcta_grm_from_design(
    design = exposure_encoded$matrix,
    ids = ids,
    prefix = e_prefix,
    block_size = cfg$block_size %||% 500L,
    n_contributors = exposure_encoded$n_features
  )
  write_tsv(exposure_encoded$feature_meta, paste0(e_prefix, "_features.tsv"))
  write_tsv(
    data.frame(
      metric = c(
        "exposure_missing_strategy",
        "complete_case_exposures_requested",
        "centered",
        "n_samples",
        "n_encoded_features"
      ),
      value = c(
        exposure_missing_strategy,
        as.integer(complete_case_requested),
        as.integer(center_exposures),
        nrow(exposure_encoded$matrix),
        exposure_encoded$n_features
      ),
      stringsAsFactors = FALSE
    ),
    paste0(e_prefix, "_build_metadata.tsv")
  )
}

gxe_prefix <- file.path(paths$kernels, paste0("GxE_", exposure_mode))
if (!file.exists(paste0(gxe_prefix, ".grm.bin")) || force) {
  timestamp_msg("Building Hadamard GxE GRM:", gxe_prefix)
  stream_hadamard_grm(
    prefix_a = geno_prefix,
    prefix_b = e_prefix,
    prefix_out = gxe_prefix,
    contributor_count = exposure_encoded$n_features
  )
}

if (include_covar_kernels) {
  covar_aligned <- align_export(covar_core, ids)
  covar_features <- covar_aligned[, -(1L:2L), drop = FALSE]
  covar_encoded <- make_numeric_matrix(
    covar_features,
    ordinal_ids = character(),
    impute_method = "mean",
    center = TRUE,
    scale_columns = TRUE
  )

  c_prefix <- file.path(paths$kernels, "Covars_core")
  if (!file.exists(paste0(c_prefix, ".grm.bin")) || force) {
    timestamp_msg("Building covariate sensitivity GRM:", c_prefix)
    write_gcta_grm_from_design(
      design = covar_encoded$matrix,
      ids = ids,
      prefix = c_prefix,
      block_size = cfg$block_size %||% 500L,
      n_contributors = covar_encoded$n_features
    )
    write_tsv(covar_encoded$feature_meta, paste0(c_prefix, "_features.tsv"))
  }

  cg_prefix <- file.path(paths$kernels, "CovarsxG")
  ce_prefix <- file.path(paths$kernels, paste0("CovarsxE_", exposure_mode))

  if (!file.exists(paste0(cg_prefix, ".grm.bin")) || force) {
    timestamp_msg("Building CovarsxG sensitivity GRM:", cg_prefix)
    stream_hadamard_grm(geno_prefix, c_prefix, cg_prefix, contributor_count = covar_encoded$n_features)
  }
  if (!file.exists(paste0(ce_prefix, ".grm.bin")) || force) {
    timestamp_msg("Building CovarsxE sensitivity GRM:", ce_prefix)
    stream_hadamard_grm(c_prefix, e_prefix, ce_prefix, contributor_count = min(covar_encoded$n_features, exposure_encoded$n_features))
  }
}

write_tsv(
  data.frame(
    kernel = c("G", "E", "GxE", "Covars", "CovarsxG", "CovarsxE"),
    prefix = c(
      geno_prefix,
      e_prefix,
      gxe_prefix,
      file.path(paths$kernels, "Covars_core"),
      file.path(paths$kernels, "CovarsxG"),
      file.path(paths$kernels, paste0("CovarsxE_", exposure_mode))
    ),
    stringsAsFactors = FALSE
  ),
  file.path(paths$kernels, "kernel_manifest.tsv")
)

timestamp_msg("Environment and interaction kernels are ready in", paths$kernels)
