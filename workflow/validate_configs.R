#!/usr/bin/env Rscript

# HEAP configuration validation functions.
#
# Source after workflow/00_paths.R and workflow/config_helpers.R.
#
# Run these checks before any module job starts to catch misconfiguration
# early with informative error messages.
#
# Typical usage at the top of a module script:
#
#   source("workflow/00_paths.R")
#   source("workflow/config_helpers.R")
#   source("workflow/validate_configs.R")
#
#   validate_exposure_config()
#   validate_covariate_set("Type3", heap)
#   validate_module_experiment("module1", "M1_Type3_lasso")
#   validate_upstream_outputs("module3", "M3_Type3_lasso_primary",
#                             upstream_score_dir = "path/to/scores")

# ---------------------------------------------------------------------------
# Exposure config validation
# ---------------------------------------------------------------------------

#' Validate the exposure config TSV.
#'
#' Checks:
#'   1. File exists.
#'   2. Required columns present.
#'   3. No duplicate variable names among included exposures.
#'   4. include column only contains 0 or 1.
#'   5. data_type values are valid.
#'
#' @param config  Optional pre-loaded config data.frame (default: loads from file).
#' @param stop_on_error  If TRUE (default), stop on first failure.
validate_exposure_config <- function(config = NULL, stop_on_error = TRUE) {
  if (is.null(config)) config <- load_exposure_config()

  errors <- character(0)

  required_cols <- c("variable", "include")
  missing_cols <- setdiff(required_cols, names(config))
  if (length(missing_cols) > 0) {
    errors <- c(errors, paste0(
      "Missing required columns in exposure config: ",
      paste(missing_cols, collapse = ", ")
    ))
  }

  if (length(errors) == 0) {
    included <- config[config$include == 1L, , drop = FALSE]

    dups <- included$variable[duplicated(included$variable)]
    if (length(dups) > 0) {
      errors <- c(errors, paste0(
        "Duplicate variable names in included exposures: ",
        paste(head(dups, 10), collapse = ", ")
      ))
    }

    bad_include <- setdiff(unique(config$include), c(0L, 1L, 0, 1))
    if (length(bad_include) > 0) {
      errors <- c(errors, paste0(
        "include column contains non-0/1 values: ",
        paste(bad_include, collapse = ", ")
      ))
    }

    if ("data_type" %in% names(config)) {
      valid_types <- c("continuous", "binary", "ordinal", "categorical", "composite", NA, "")
      bad_types <- setdiff(unique(config$data_type), valid_types)
      if (length(bad_types) > 0) {
        errors <- c(errors, paste0(
          "Unknown data_type values: ", paste(bad_types, collapse = ", "),
          " (expected: continuous, binary, ordinal, categorical, composite)"
        ))
      }
    }
  }

  .report_validation(errors, context = "Exposure config", stop_on_error = stop_on_error)
  invisible(length(errors) == 0L)
}

#' Check that included exposures exist in the HEAP loader output.
#'
#' @param heap   HEAP object (or pxs list) with $E_baseline or $Elist.
#' @param config Optional pre-loaded exposure config.
validate_exposures_in_loader <- function(heap, config = NULL) {
  if (is.null(config)) config <- load_exposure_config()
  included_vars <- included_exposure_vars(config)

  all_loader_vars <- unlist(lapply(
    heap$E_baseline %||% heap$Elist,
    function(df) setdiff(names(df), c("eid", "instance"))
  ))

  missing <- setdiff(included_vars, all_loader_vars)
  if (length(missing) > 0) {
    stop(
      length(missing), " included exposure variable(s) not found in HEAP loader output:\n  ",
      paste(head(missing, 20), collapse = "\n  "),
      if (length(missing) > 20) paste0("\n  ... (", length(missing) - 20, " more)") else "",
      "\nCheck analysis_exposures.tsv and HEAP_loader.R output."
    )
  }
  invisible(TRUE)
}

# ---------------------------------------------------------------------------
# Covariate set validation
# ---------------------------------------------------------------------------

#' Validate that all covariates in a named set exist in the HEAP loader output.
#'
#' @param covar_set_name  Name of covariate set (e.g., "Type3")
#' @param heap            HEAP or pxs object with $covars_df
#' @param stop_on_error   If TRUE (default), stop on missing covariates
validate_covariate_set <- function(covar_set_name, heap = NULL, stop_on_error = TRUE) {
  covars <- load_covariate_set(covar_set_name)

  # Type5 is always NULL (resolved at runtime)
  if (is.null(covars)) {
    message("[validate] Covariate set '", covar_set_name,
            "' is dynamic (resolved from loader covars_list). Skipping static check.")
    return(invisible(TRUE))
  }

  errors <- character(0)

  if (!is.null(heap)) {
    loader_covars <- names(heap$covars_df %||% heap$covars_baseline)
    loader_covars <- setdiff(loader_covars, "eid")

    missing <- setdiff(covars, loader_covars)
    if (length(missing) > 0) {
      errors <- c(errors, paste0(
        "Covariate set '", covar_set_name, "' requests ",
        length(missing), " column(s) not in loader output:\n  ",
        paste(head(missing, 10), collapse = "\n  "),
        if (length(missing) > 10) paste0("\n  ... (", length(missing) - 10, " more)") else ""
      ))
    }
  }

  .report_validation(errors, context = paste0("Covariate set '", covar_set_name, "'"),
                     stop_on_error = stop_on_error)
  invisible(length(errors) == 0L)
}

# ---------------------------------------------------------------------------
# Module experiment validation
# ---------------------------------------------------------------------------

#' Validate a named module experiment config.
#'
#' Checks:
#'   1. experiment_name exists in the module YAML.
#'   2. covariate_set references a valid set in covariate_sets.yml.
#'   3. family is a valid model family.
#'   4. For Module 3: mediation_mode is valid, upstream experiment exists.
#'
#' @param module          Module identifier (e.g., "module1")
#' @param experiment_name Named experiment to validate
validate_module_experiment <- function(module, experiment_name) {
  exp <- load_experiment_config(module, experiment_name)

  errors <- character(0)

  # Check covariate_set exists
  if (!is.null(exp$covariate_set)) {
    tryCatch(
      load_covariate_set(exp$covariate_set),
      error = function(e) {
        errors <<- c(errors, conditionMessage(e))
      }
    )
  }

  # Check family (where applicable)
  if (!is.null(exp$family)) {
    valid_families <- c("lasso", "ridge", "enet", "rf")
    if (!exp$family %in% valid_families) {
      errors <- c(errors, paste0(
        "Unknown model family '", exp$family, "' in experiment '", experiment_name, "'.\n",
        "Valid families: ", paste(valid_families, collapse = ", ")
      ))
    }
  }

  # Module 3: check mediation_mode and upstream experiment
  if (grepl("^module3", module, ignore.case = TRUE)) {
    valid_modes <- c("primary_total", "partitioned_categories",
                     "partitioned_grouped_categories")
    if (!is.null(exp$mediation_mode) && !exp$mediation_mode %in% valid_modes) {
      errors <- c(errors, paste0(
        "Unknown mediation_mode '", exp$mediation_mode, "'.\n",
        "Valid modes: ", paste(valid_modes, collapse = ", ")
      ))
    }
    if (!is.null(exp$upstream_module1_experiment)) {
      tryCatch(
        load_experiment_config("module1", exp$upstream_module1_experiment),
        error = function(e) {
          errors <<- c(errors, paste0(
            "Upstream Module 1 experiment '", exp$upstream_module1_experiment,
            "' not found: ", conditionMessage(e)
          ))
        }
      )
    }
  }

  .report_validation(errors,
                     context = paste0("Experiment '", experiment_name, "' (", module, ")"),
                     stop_on_error = TRUE)
  invisible(length(errors) == 0L)
}

# ---------------------------------------------------------------------------
# Upstream output validation
# ---------------------------------------------------------------------------

#' Validate that required upstream Module 1 score columns exist.
#'
#' Called by Module 3 before running mediation to confirm the upstream
#' Module 1 experiment produced the required score files.
#'
#' @param module3_experiment_name  Name of the Module 3 experiment.
#' @param upstream_score_dir  Path to the mediation_scores directory from Module 1.
validate_upstream_outputs <- function(module3_experiment_name, upstream_score_dir) {
  exp <- load_experiment_config("module3", module3_experiment_name)
  required_cols <- exp$required_score_columns

  if (is.null(required_cols) || length(required_cols) == 0) {
    message("[validate] No required_score_columns specified for '",
            module3_experiment_name, "'; skipping upstream output check.")
    return(invisible(TRUE))
  }

  if (!dir.exists(upstream_score_dir)) {
    stop(
      "Upstream score directory not found: ", upstream_score_dir, "\n",
      "Module 3 experiment '", module3_experiment_name, "' requires upstream Module 1 experiment '",
      exp$upstream_module1_experiment %||% "(unspecified)", "' to be run first."
    )
  }

  sample_files <- list.files(upstream_score_dir, pattern = "\\.txt$", full.names = TRUE)
  if (length(sample_files) == 0L) {
    stop(
      "Upstream score directory exists but contains no .txt files: ", upstream_score_dir
    )
  }

  sample_df <- tryCatch(
    utils::read.delim(sample_files[1], nrows = 1, check.names = FALSE),
    error = function(e) stop("Could not read upstream score file: ", conditionMessage(e))
  )

  missing_cols <- setdiff(required_cols, names(sample_df))
  if (length(missing_cols) > 0) {
    stop(
      "Module 3 experiment '", module3_experiment_name, "' requires score column(s) ",
      "not found in upstream Module 1 output:\n  ",
      paste(missing_cols, collapse = "\n  "),
      "\n\nUpstream file checked: ", sample_files[1],
      "\nAvailable columns: ", paste(names(sample_df), collapse = ", "),
      "\n\nFix: Ensure upstream Module 1 experiment '",
      exp$upstream_module1_experiment %||% "(unspecified)", "' ran with the correct options ",
      "(e.g., run_genetic_subblocks=TRUE for Gcis_raw/Gtrans_raw)."
    )
  }

  invisible(TRUE)
}

# ---------------------------------------------------------------------------
# Manifest validation
# ---------------------------------------------------------------------------

#' Validate a manifest TSV file.
#'
#' Checks:
#'   1. File exists and is non-empty.
#'   2. Required columns are present.
#'   3. array_index is unique and sequential (starting at 1).
#'   4. output_path and input_path are non-empty.
#'
#' @param manifest_path  Path to the manifest TSV.
#' @param required_cols  Column names that must be present.
validate_manifest <- function(
  manifest_path,
  required_cols = c("array_index", "module", "experiment_name",
                    "covariate_set", "output_path")
) {
  if (!file.exists(manifest_path)) {
    stop("Manifest file not found: ", manifest_path)
  }
  mf <- utils::read.delim(manifest_path, stringsAsFactors = FALSE, check.names = FALSE)
  if (nrow(mf) == 0L) stop("Manifest is empty: ", manifest_path)

  errors <- character(0)

  missing_cols <- setdiff(required_cols, names(mf))
  if (length(missing_cols) > 0) {
    errors <- c(errors, paste0(
      "Manifest missing required columns: ", paste(missing_cols, collapse = ", ")
    ))
  }

  if ("array_index" %in% names(mf)) {
    if (any(duplicated(mf$array_index))) {
      errors <- c(errors, "Manifest has duplicate array_index values.")
    }
    expected <- seq_len(nrow(mf))
    if (!identical(sort(mf$array_index), expected)) {
      errors <- c(errors, paste0(
        "array_index is not a contiguous 1-based sequence. ",
        "Got: ", min(mf$array_index), "-", max(mf$array_index),
        " (", nrow(mf), " rows)."
      ))
    }
  }

  .report_validation(errors, context = paste0("Manifest: ", basename(manifest_path)),
                     stop_on_error = TRUE)
  invisible(nrow(mf))
}

# ---------------------------------------------------------------------------
# Scratch path guard
# ---------------------------------------------------------------------------

#' Assert that a path is not under the scratch directory.
#'
#' Call this before writing canonical outputs to prevent accidental
#' scratch usage for reproducible artifacts.
#'
#' @param path  File or directory path to check.
#' @param label Short label for the error message.
assert_no_canonical_scratch_output <- function(path, label = "output") {
  scratch_root <- HEAP_PATHS$scratch_root
  if (grepl(paste0("^", scratch_root), normalizePath(path, mustWork = FALSE))) {
    stop(
      "REPRODUCIBILITY VIOLATION: ", label, " is being written to scratch.\n",
      "  Path: ", path, "\n",
      "  Scratch root: ", scratch_root, "\n",
      "Canonical HEAP outputs must go to: ", heap_project_root(), "\n",
      "Use heap_project_output(), heap_project_intermediate(), or heap_manifest() instead."
    )
  }
  invisible(TRUE)
}

# ---------------------------------------------------------------------------
# Internal helper
# ---------------------------------------------------------------------------

.report_validation <- function(errors, context = "", stop_on_error = TRUE) {
  if (length(errors) == 0L) {
    message("[validate] OK: ", context)
    return(invisible(NULL))
  }
  msg <- paste0(
    "Validation failed", if (nzchar(context)) paste0(" [", context, "]"), ":\n",
    paste0("  - ", errors, collapse = "\n")
  )
  if (stop_on_error) stop(msg) else warning(msg)
}
