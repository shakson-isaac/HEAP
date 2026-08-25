#!/usr/bin/env Rscript

# HEAP shared config loading helpers.
#
# Source this file after workflow/00_paths.R to get all config loaders.
# All functions read from config/ files under HEAP root; none write.
#
# Quick reference:
#   load_exposure_config()           -> data.frame of analysis_exposures.tsv
#   load_exposure_groups()           -> data.frame of category -> broad_group mapping
#   load_covariate_sets()            -> full parsed YAML list
#   load_covariate_set("Type3")      -> character vector of covariate names
#                                       (or NULL for Type5)
#   load_module_experiments("module1") -> full parsed YAML experiment list
#   load_experiment_config("module1", "M1_Type3_lasso") -> single experiment list
#   resolve_manifest_row(manifest_path, array_index) -> single-row data.frame
#   write_run_config(cfg_list, output_dir) -> writes run_config.yml artifact

# ---------------------------------------------------------------------------
# Internal utilities
# ---------------------------------------------------------------------------

.heap_load_yaml <- function(path) {
  if (!requireNamespace("yaml", quietly = TRUE)) {
    stop(
      "Package 'yaml' is required for HEAP config loading.\n",
      "Install with: install.packages('yaml')"
    )
  }
  if (!file.exists(path)) stop("Config file not found: ", path)
  yaml::read_yaml(path)
}

.heap_require_paths <- function() {
  if (!exists("HEAP_PATHS", envir = .GlobalEnv)) {
    stop(
      "HEAP_PATHS not found. Source workflow/00_paths.R before config_helpers.R."
    )
  }
}

# ---------------------------------------------------------------------------
# Exposure config
# ---------------------------------------------------------------------------

#' Load the primary exposure inclusion/exclusion config.
#' Returns a data.frame with columns: variable, exposure_id, exposure_label,
#' category, analysis_group, data_type, missingness, include, reason_excluded, notes.
load_exposure_config <- function(path = NULL) {
  .heap_require_paths()
  if (is.null(path)) path <- heap_config("exposure_sets", "analysis_exposures.tsv")
  if (!file.exists(path)) stop("Exposure config not found: ", path)
  cfg <- utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
  cfg
}

#' Return only the included exposure variable names (include == 1).
included_exposure_vars <- function(config = NULL) {
  if (is.null(config)) config <- load_exposure_config()
  config$variable[config$include == 1L]
}

#' Load the broad category grouping file for grouped mediation (Module 3).
load_exposure_groups <- function(path = NULL) {
  .heap_require_paths()
  if (is.null(path)) {
    path <- heap_config("exposure_sets", "analysis_exposure_category_groups.tsv")
  }
  if (!file.exists(path)) stop("Exposure groups file not found: ", path)
  utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
}

# ---------------------------------------------------------------------------
# Covariate sets
# ---------------------------------------------------------------------------

#' Load the full covariate sets YAML as a list.
load_covariate_sets <- function(path = NULL) {
  .heap_require_paths()
  if (is.null(path)) path <- heap_config("covariates", "covariate_sets.yml")
  cfg <- .heap_load_yaml(path)
  if (is.null(cfg$covariate_sets)) stop("covariate_sets key missing in: ", path)
  cfg
}

#' Return the covariate names for a named set.
#'
#' Returns NULL for Type5 (full loader covariates resolved at runtime).
#' Returns a character vector for all other types.
#'
#' @param name  Covariate set name, e.g., "Type3"
#' @param path  Optional path override for covariate_sets.yml
load_covariate_set <- function(name, path = NULL) {
  cfg <- load_covariate_sets(path)
  sets <- cfg$covariate_sets
  if (!name %in% names(sets)) {
    stop(
      "Covariate set '", name, "' not found in covariate_sets.yml.\n",
      "Available sets: ", paste(names(sets), collapse = ", ")
    )
  }
  entry <- sets[[name]]
  covars <- entry$covariates

  if (is.null(covars) || (length(covars) == 1L && is.na(covars[[1]]))) {
    return(NULL)
  }

  as.character(covars)
}

#' Return the discrete covariate sub-list for a named set (population_architecture use).
load_covariate_set_discrete <- function(name, path = NULL) {
  cfg <- load_covariate_sets(path)
  entry <- cfg$covariate_sets[[name]]
  if (is.null(entry)) stop("Covariate set '", name, "' not found.")
  as.character(entry$discrete %||% character(0))
}

#' Return the kernel_quantitative sub-list (population_architecture use).
load_covariate_set_kernel <- function(name, path = NULL) {
  cfg <- load_covariate_sets(path)
  entry <- cfg$covariate_sets[[name]]
  if (is.null(entry)) stop("Covariate set '", name, "' not found.")
  as.character(entry$kernel_quantitative %||% character(0))
}

#' Return the quantitative sub-list for a named set (population_architecture use).
load_covariate_set_quantitative <- function(name, path = NULL) {
  cfg <- load_covariate_sets(path)
  entry <- cfg$covariate_sets[[name]]
  if (is.null(entry)) stop("Covariate set '", name, "' not found.")
  as.character(entry$quantitative %||% character(0))
}

#' Return the E->C remapping spec for a named covariate set, or NULL if none.
#'
#' base_ses carries a `remapping` block describing deprivation variables that live
#' in an exposure category (Deprivation_Indices) and must be moved into the
#' covariate matrix at runtime (Module2 only). Returns a list with $type,
#' $source_category, $variables, or NULL when the set has no remapping.
load_covariate_set_remapping <- function(name, path = NULL) {
  cfg <- load_covariate_sets(path)
  entry <- cfg$covariate_sets[[name]]
  if (is.null(entry)) stop("Covariate set '", name, "' not found.")
  entry$remapping
}

#' Resolve Type5 (NULL) to the loader's covars_list at runtime.
#'
#' @param covar_set_name  Named covariate set (e.g., "Type5")
#' @param pxs             PXS/HEAP object with $covars_list
resolve_covariate_set <- function(covar_set_name, pxs) {
  covars <- load_covariate_set(covar_set_name)
  if (is.null(covars)) {
    if (is.null(pxs$covars_list)) {
      stop("A NULL covariate set requires pxs$covars_list, but it is NULL.")
    }
    covars <- pxs$covars_list
  }
  covars
}

# ---------------------------------------------------------------------------
# Sample filters (the WHO-IS-IN axis, orthogonal to covariate adjustment)
# ---------------------------------------------------------------------------

#' Load the full sample_filters YAML as a list.
load_sample_filters <- function(path = NULL) {
  .heap_require_paths()
  if (is.null(path)) path <- heap_config("samples", "sample_filters.yml")
  cfg <- .heap_load_yaml(path)
  if (is.null(cfg$sample_filters)) stop("sample_filters key missing in: ", path)
  cfg
}

#' Return a single sample-filter spec by name, or NULL for no restriction.
#'
#' NULL (no-op) is returned for a NULL/NA/empty name, the literal "none", or a
#' filter entry whose `column` is null. Otherwise returns the spec list with an
#' added $name element: list(name, column, keep_values, drop_na, ...).
load_sample_filter <- function(name, path = NULL) {
  if (is.null(name) || (length(name) == 1L && (is.na(name) || !nzchar(name))) ||
      identical(name, "none")) {
    return(NULL)
  }
  cfg  <- load_sample_filters(path)
  filt <- cfg$sample_filters[[name]]
  if (is.null(filt)) {
    stop("Sample filter '", name, "' not found in sample_filters.yml.\n",
         "Available: ", paste(names(cfg$sample_filters), collapse = ", "))
  }
  if (is.null(filt$column)) return(NULL)   # e.g. the explicit "none" entry
  filt$name <- name
  filt
}

#' Apply a sample-filter spec to a covariate data.frame.
#'
#' Keeps rows where df[[spec$column]] is in spec$keep_values. Rows with NA in the
#' filter column are dropped iff spec$drop_na is TRUE (else retained). A NULL spec
#' is a no-op. MUST be called on the FULL covariate frame BEFORE the covariate set
#' narrows columns, or the filter column may already be gone.
#'
#' @return list(df, n_before, n_after, n_dropped, filter)
apply_sample_filter <- function(df, spec) {
  if (is.null(spec)) {
    return(list(df = df, n_before = nrow(df), n_after = nrow(df),
                n_dropped = 0L, filter = NA_character_))
  }
  col <- spec$column
  if (!col %in% names(df)) {
    stop("Sample filter '", spec$name %||% "?", "' needs column '", col,
         "', absent from the data frame. Was the loader rerun with the feature, ",
         "and is the filter applied before covariate narrowing?")
  }
  n_before <- nrow(df)
  x    <- df[[col]]
  keep <- x %in% spec$keep_values            # NA %in% ... is FALSE -> NA rows out
  if (!isTRUE(spec$drop_na)) keep <- keep | is.na(x)
  out <- df[keep, , drop = FALSE]
  rownames(out) <- NULL
  list(df = out, n_before = n_before, n_after = nrow(out),
       n_dropped = n_before - nrow(out), filter = spec$name %||% NA_character_)
}

# ---------------------------------------------------------------------------
# Module experiment configs
# ---------------------------------------------------------------------------

.module_yml_path <- function(module) {
  .heap_require_paths()
  module_clean <- sub("^module", "module", tolower(module))
  # Accept "module1", "1", "module1_variance_decomposition" etc.
  short <- gsub("^module_?", "", module_clean)
  short <- gsub("_.*", "", short)

  # Build candidate paths
  candidates <- c(
    heap_config("modules", paste0(module_clean, "_experiments.yml")),
    heap_config("modules", paste0("module", short, "_experiments.yml")),
    heap_config("modules", paste0(module, "_experiments.yml"))
  )
  hit <- candidates[file.exists(candidates)][1]
  if (is.na(hit)) {
    stop(
      "No experiment config found for module '", module, "'.\n",
      "Searched: ", paste(candidates, collapse = "\n  ")
    )
  }
  hit
}

#' Load all experiments for a module.
#' Returns the full parsed YAML list.
#'
#' @param module  Module identifier: "module1", "module3", "module5",
#'                "module6", or "population_architecture"
load_module_experiments <- function(module) {
  path <- .module_yml_path(module)
  cfg <- .heap_load_yaml(path)
  if (is.null(cfg$experiments)) stop("'experiments' key missing in: ", path)
  cfg
}

#' Load a single named experiment config for a module.
#'
#' @param module           Module identifier
#' @param experiment_name  Name matching a key in the experiments list
load_experiment_config <- function(module, experiment_name) {
  cfg <- load_module_experiments(module)
  exps <- cfg$experiments
  if (!experiment_name %in% names(exps)) {
    stop(
      "Experiment '", experiment_name, "' not found in module '", module, "'.\n",
      "Available experiments: ", paste(names(exps), collapse = ", ")
    )
  }

  # Merge defaults onto the experiment entry
  defaults <- cfg$defaults %||% list()
  exp <- exps[[experiment_name]]
  for (k in names(defaults)) {
    if (is.null(exp[[k]])) exp[[k]] <- defaults[[k]]
  }
  exp$experiment_name <- experiment_name
  exp$module <- module
  exp
}

# ---------------------------------------------------------------------------
# Manifest row resolution
# ---------------------------------------------------------------------------

#' Read a manifest TSV and return the row for a given Slurm array index.
#'
#' @param manifest_path  Path to the manifest TSV (must have array_index column).
#' @param array_index    Integer index (1-based, matching SLURM_ARRAY_TASK_ID).
resolve_manifest_row <- function(manifest_path, array_index) {
  if (!file.exists(manifest_path)) {
    stop(
      "Manifest not found: ", manifest_path, "\n",
      "Run generate_*_manifest() first."
    )
  }
  mf <- utils::read.delim(
    manifest_path, stringsAsFactors = FALSE, check.names = FALSE
  )
  idx <- as.integer(array_index)
  if (!"array_index" %in% names(mf)) {
    stop("Manifest missing 'array_index' column: ", manifest_path)
  }
  row <- mf[mf$array_index == idx, , drop = FALSE]
  if (nrow(row) == 0L) {
    stop(
      "No manifest row found for array_index=", idx,
      " in: ", manifest_path,
      "\n(manifest has ", nrow(mf), " rows, indices ",
      min(mf$array_index), "-", max(mf$array_index), ")"
    )
  }
  if (nrow(row) > 1L) {
    stop("Duplicate array_index=", idx, " in manifest: ", manifest_path)
  }
  as.list(row[1, ])
}

# ---------------------------------------------------------------------------
# Run config artifact writer
# ---------------------------------------------------------------------------

#' Write a resolved run config artifact to a job's output directory.
#'
#' Each module job should call this before doing any computation.
#' The artifact records exactly which experiment, covariate set, model
#' family, and paths were used for this specific run.
#'
#' @param cfg_list    Named list of all resolved parameters.
#' @param output_dir  Directory to write run_config.yml into.
#' @param filename    Filename (default: run_config.yml)
write_run_config <- function(cfg_list, output_dir, filename = "run_config.yml") {
  if (!requireNamespace("yaml", quietly = TRUE)) {
    warning("Package 'yaml' not available; writing run_config as R list dump instead.")
    out_path <- file.path(output_dir, sub("\\.yml$", ".rds", filename))
    saveRDS(cfg_list, out_path)
    return(invisible(out_path))
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  cfg_list$run_timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  cfg_list$run_host <- Sys.info()[["nodename"]]
  out_path <- file.path(output_dir, filename)
  yaml::write_yaml(cfg_list, out_path)
  message("[run_config] Written: ", out_path)
  invisible(out_path)
}

# ---------------------------------------------------------------------------
# Utility
# ---------------------------------------------------------------------------

`%||%` <- function(x, y) if (!is.null(x)) x else y
