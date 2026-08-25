#!/usr/bin/env Rscript

# HEAP manifest generation functions.
#
# Source after workflow/00_paths.R and workflow/config_helpers.R.
#
# Each generate_*_manifest() function reads a named experiment from the
# config YAML, expands it into one row per Slurm job, and writes a TSV
# manifest to heap_manifest(module, experiment_name.tsv).
#
# The manifest is the contract between:
#   - Human  : chose the experiment name and parameters
#   - Config : encoded those choices in a versioned YAML
#   - Script : generate_*_manifest() expanded them into concrete jobs
#   - Slurm  : schedules jobs by array index; reads exactly this file
#   - Module R script: reads its row; runs; writes outputs to output_path
#
# Quick reference:
#   generate_module1_manifest("M1_Type3_lasso")
#   generate_module3_manifest("M3_Type3_lasso_primary")
#   generate_module5_manifest("MR_UKB_primary")
#   generate_module6_manifest("M6_compact_Type5")
#   generate_population_architecture_manifest("PopArch_Type3_greml_total")
#   validate_manifest(heap_manifest("module1", "M1_Type3_lasso.tsv"))

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

.write_manifest <- function(df, module, experiment_name, overwrite = TRUE) {
  out_dir <- heap_manifest(module)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out_path <- file.path(out_dir, paste0(experiment_name, ".tsv"))
  if (file.exists(out_path) && !overwrite) {
    stop("Manifest already exists (use overwrite=TRUE): ", out_path)
  }
  utils::write.table(
    df, out_path,
    sep = "\t", row.names = FALSE, quote = FALSE
  )
  message("[manifest] Written: ", out_path, " (", nrow(df), " rows)")
  invisible(out_path)
}

.load_protein_list <- function(protein_set = "omicspred_2704") {
  .heap_require_paths()
  prot_path <- heap_config("protein_sets", "omicspred_proteins.txt")
  if (!file.exists(prot_path)) stop("Protein list not found: ", prot_path)
  prots <- readLines(prot_path)
  prots <- trimws(prots[nzchar(trimws(prots))])
  prots
}

# ---------------------------------------------------------------------------
# Module 1 manifest
# ---------------------------------------------------------------------------

#' Generate the Module 1 manifest for a named experiment.
#'
#' Each row represents one protein chunk × covariate_set × family job.
#' Slurm array index maps to a row.
#'
#' @param experiment_name  Name from config/modules/module1_experiments.yml
#' @param n_chunks         Number of array jobs (default: from config, fallback 400)
#' @param overwrite        Overwrite existing manifest if TRUE
generate_module1_manifest <- function(experiment_name,
                                      n_chunks = NULL,
                                      overwrite = TRUE) {
  exp <- load_experiment_config("module1", experiment_name)
  validate_module_experiment("module1", experiment_name)

  n_chunks <- n_chunks %||% exp$n_chunks %||% 400L
  n_chunks <- as.integer(n_chunks)

  out_root <- heap_project_output(
    exp$output_subdir %||% "module1_predictive_r2_score_partition",
    experiment_name
  )

  rows <- data.frame(
    array_index      = seq_len(n_chunks),
    module           = "module1",
    experiment_name  = experiment_name,
    covariate_set    = exp$covariate_set,
    family           = exp$family,
    sample_filter    = exp$sample_filter %||% "none",
    chunk_id         = seq_len(n_chunks),
    n_chunks         = n_chunks,
    seed             = exp$seed %||% 123L,
    kfold            = exp$kfold %||% 5L,
    glmnet_inner_kfold = exp$glmnet_inner_kfold %||% 5L,
    miss_rate        = exp$miss_rate %||% 0.20,
    decomposition_mode = exp$decomposition_mode %||% "score_partition",
    run_genetic_subblocks   = as.integer(
      exp$run_genetic_subblocks %||% TRUE
    ),
    run_exposure_categories = as.integer(
      exp$run_exposure_categories %||% TRUE
    ),
    run_gxe_categories      = as.integer(
      exp$run_gxe_categories %||% TRUE
    ),
    run_gxc_sensitivity     = as.integer(
      exp$run_gxc_sensitivity %||% FALSE
    ),
    run_exc_sensitivity     = as.integer(
      exp$run_exc_sensitivity %||% FALSE
    ),
    gxc_covars = paste(
      exp$gxc_covars %||% c(
        "age_when_attended_assessment_centre_f21003_0_0",
        "sex_f31_0_0"
      ),
      collapse = ","
    ),
    exc_covars = paste(
      exp$exc_covars %||% c(
        "age_when_attended_assessment_centre_f21003_0_0",
        "sex_f31_0_0"
      ),
      collapse = ","
    ),
    config_path = heap_config("modules", "module1_experiments.yml"),
    output_path = out_root,
    stringsAsFactors = FALSE
  )

  .write_manifest(rows, "module1", experiment_name, overwrite = overwrite)
}

# ---------------------------------------------------------------------------
# Module 2 manifest
# ---------------------------------------------------------------------------

#' Generate the Module 2 manifest for a named experiment.
#'
#' Each row represents one protein chunk for the given covariate_set run.
#'
#' @param experiment_name  Name from config/modules/module2_experiments.yml
#' @param n_chunks         Number of array jobs (default: from config, fallback 400)
#' @param overwrite        Overwrite existing manifest if TRUE
generate_module2_manifest <- function(experiment_name,
                                      n_chunks = NULL,
                                      overwrite = TRUE) {
  exp <- load_experiment_config("module2", experiment_name)

  n_chunks <- n_chunks %||% exp$n_chunks %||% 400L
  n_chunks <- as.integer(n_chunks)

  out_root <- heap_project_output(
    exp$output_subdir %||% "module2",
    experiment_name
  )

  rows <- data.frame(
    array_index      = seq_len(n_chunks),
    module           = "module2",
    experiment_name  = experiment_name,
    covariate_set    = exp$covariate_set,
    covar_variant    = exp$covar_variant %||% "",
    sample_filter    = exp$sample_filter %||% "none",
    chunk_id         = seq_len(n_chunks),
    n_chunks         = n_chunks,
    config_path      = heap_config("modules", "module2_experiments.yml"),
    output_path      = out_root,
    stringsAsFactors = FALSE
  )

  .write_manifest(rows, "module2", experiment_name, overwrite = overwrite)
}

# ---------------------------------------------------------------------------
# Module 3 manifest
# ---------------------------------------------------------------------------

#' Generate the Module 3 manifest for a named mediation experiment.
#'
#' Each row represents one protein chunk × mediation_mode job.
#'
#' @param experiment_name  Name from config/modules/module3_experiments.yml
#' @param n_chunks         Number of array jobs (default: from config, fallback 1000)
#' @param overwrite        Overwrite existing manifest if TRUE
generate_module3_manifest <- function(experiment_name,
                                      n_chunks = NULL,
                                      overwrite = TRUE) {
  exp <- load_experiment_config("module3", experiment_name)
  validate_module_experiment("module3", experiment_name)

  n_chunks <- n_chunks %||% exp$n_chunks %||% 1000L
  n_chunks <- as.integer(n_chunks)

  upstream_exp <- exp$upstream_module1_experiment

  # Upstream score dir: must point to the Module 1 experiment's mediation_scores output.
  # Module 1 writes:  out_root / covarType / family / mediation_scores
  # where out_root = heap_project_output(<m1_output_subdir>, <m1_experiment_name>)
  .m1_exp <- tryCatch(
    load_experiment_config("module1", upstream_exp),
    error = function(e) list(output_subdir = "module1_predictive_r2_score_partition")
  )
  upstream_score_dir <- heap_project_output(
    .m1_exp$output_subdir %||% "module1_predictive_r2_score_partition",
    upstream_exp,
    exp$covariate_set,
    exp$family,
    "mediation_scores"
  )

  out_root <- heap_project_output(
    exp$output_subdir %||% "module3",
    experiment_name
  )

  rows <- data.frame(
    array_index      = seq_len(n_chunks),
    module           = "module3",
    experiment_name  = experiment_name,
    covariate_set    = exp$covariate_set,
    family           = exp$family,
    mediation_mode   = exp$mediation_mode,
    sample_filter    = exp$sample_filter %||% "none",
    chunk_id         = seq_len(n_chunks),
    n_chunks         = n_chunks,
    upstream_module1_experiment = upstream_exp %||% "",
    upstream_score_dir = upstream_score_dir,
    delta_sd_total   = exp$delta_sd_total %||% 1.0,
    min_cases_per_disease = exp$min_cases_per_disease %||% 100L,
    min_eid_overlap  = exp$min_eid_overlap %||% 500L,
    disease_filter   = exp$disease_filter %||% "all_available",
    config_path      = heap_config("modules", "module3_experiments.yml"),
    output_path      = out_root,
    stringsAsFactors = FALSE
  )

  .write_manifest(rows, "module3", experiment_name, overwrite = overwrite)
}

# ---------------------------------------------------------------------------
# Module 5 manifest
# ---------------------------------------------------------------------------

#' Generate the Module 5 (MR) manifest for a named experiment.
#'
#' Each row represents one edge_type × chunk_id job. Unlike a fixed-n_chunks
#' grid, chunks are sized PER EDGE TYPE from the actual edge-list row counts
#' (e.g. PD/DP have thousands of pairs, EP/PE only a handful), so no array
#' task processes an empty chunk. The triad/edge lists must already exist:
#' run scripts/module5_mr/Module5_load.R first.
#'
#' The `runner` column tells the launcher which R script to dispatch
#' (Module5.R for the split-sample UKB arm, Module5_deCODE.R for deCODE).
#'
#' @param experiment_name  Name from config/modules/module5_experiments.yml
#' @param pairs_per_chunk  Target edge pairs per array task (default: from
#'                         config `pairs_per_chunk`, fallback 25)
#' @param overwrite        Overwrite existing manifest if TRUE
generate_module5_manifest <- function(experiment_name,
                                      pairs_per_chunk = NULL,
                                      overwrite = TRUE) {
  exp <- load_experiment_config("module5", experiment_name)

  edge_types <- exp$edge_types %||% c("EP", "ED", "PD", "PE", "DE", "DP")
  runner     <- exp$runner %||% "Module5.R"
  ppc        <- as.integer(pairs_per_chunk %||% exp$pairs_per_chunk %||% 25L)

  # HEAP_MR_EDGES_DIR override -> generate the manifest against the registry's
  # delta edge lists (only missing edges) for an incremental run.
  edges_dir <- Sys.getenv("HEAP_MR_EDGES_DIR", unset = heap_project_output("mr_edges", "global_edges"))

  rows_list <- lapply(edge_types, function(et) {
    ef <- file.path(edges_dir, paste0("edges_", et, ".tsv"))
    if (!file.exists(ef))
      stop("Edge file not found: ", ef,
           "\nRun scripts/module5_mr/Module5_load.R first to build the ",
           "triad/edge lists.")
    n_edges  <- nrow(data.table::fread(ef, showProgress = FALSE))
    n_chunks <- max(1L, as.integer(ceiling(n_edges / ppc)))
    data.frame(
      edge_type = et,
      chunk_id  = seq_len(n_chunks),
      n_chunks  = n_chunks,
      n_edges   = n_edges,
      stringsAsFactors = FALSE
    )
  })
  rows <- do.call(rbind, rows_list)
  rows$array_index     <- seq_len(nrow(rows))
  rows$module          <- "module5"
  rows$experiment_name <- experiment_name
  rows$runner          <- runner
  rows$p_threshold     <- exp$p_threshold %||% 5e-8
  rows$f_statistic_min <- exp$f_statistic_min %||% 10L
  rows$clump_r2        <- exp$clump_r2 %||% 0.001
  rows$clump_kb        <- exp$clump_kb %||% 10000L
  rows$cis_window_bp   <- exp$cis_window_bp %||% 1e6
  rows$output_path     <- heap_project_output(
    exp$output_subdir %||% "mr_edges",
    experiment_name
  )
  rows$config_path     <- heap_config("modules", "module5_experiments.yml")

  # Reorder
  rows <- rows[, c(
    "array_index", "module", "experiment_name", "runner", "edge_type",
    "chunk_id", "n_chunks", "n_edges", "p_threshold", "f_statistic_min",
    "clump_r2", "clump_kb", "cis_window_bp", "output_path", "config_path"
  )]

  .write_manifest(rows, "module5", experiment_name, overwrite = overwrite)
}

# ---------------------------------------------------------------------------
# Module 6 manifest
# ---------------------------------------------------------------------------

#' Generate the Module 6 manifest for a named experiment.
#'
#' For prod sub-workflow: one row per exposure in the exposure manifest.
#' For compact sub-workflow: one row per exposure × K × selection_mode.
#'
#' @param experiment_name  Name from config/modules/module6_experiments.yml
#' @param overwrite        Overwrite existing manifest if TRUE
generate_module6_manifest <- function(experiment_name, overwrite = TRUE) {
  exp <- load_experiment_config("module6", experiment_name)

  sub_wf <- exp$sub_workflow %||% "prod"

  # Load exposure list
  exp_manifest_path <- file.path(HEAP_PATHS$heap_root, exp$exposure_manifest)
  if (!file.exists(exp_manifest_path)) {
    # Try as absolute path
    exp_manifest_path <- exp$exposure_manifest
  }
  if (!file.exists(exp_manifest_path)) {
    stop("Exposure manifest not found: ", exp_manifest_path,
         "\nSet exposure_manifest in module6_experiments.yml.")
  }

  if (grepl("\\.txt$", exp_manifest_path)) {
    exposures <- readLines(exp_manifest_path)
    exposures <- trimws(exposures[nzchar(trimws(exposures))])
  } else {
    ecfg <- utils::read.delim(exp_manifest_path, stringsAsFactors = FALSE)
    exposures <- ecfg$variable[ecfg$include == 1L]
  }

  out_root <- heap_project_output(
    exp$output_subdir %||% "module6",
    experiment_name
  )

  if (sub_wf == "prod") {
    rows <- data.frame(
      array_index      = seq_along(exposures),
      module           = "module6",
      experiment_name  = experiment_name,
      sub_workflow     = "prod",
      covariate_set    = exp$covariate_set,
      sample_filter    = exp$sample_filter %||% "none",
      exposure_id      = exposures,
      seed             = exp$seed %||% 123L,
      kfold            = exp$kfold %||% 5L,
      miss_rate_prot   = exp$miss_rate_prot %||% 0.20,
      run_cross_sectional  = as.integer(exp$run_cross_sectional  %||% TRUE),
      run_longitudinal     = as.integer(exp$run_longitudinal     %||% TRUE),
      run_cox_validation   = as.integer(exp$run_cox_validation   %||% FALSE),
      output_path      = file.path(out_root, exposures),
      config_path      = heap_config("modules", "module6_experiments.yml"),
      stringsAsFactors = FALSE
    )

  } else if (sub_wf == "compact") {
    ks_raw <- exp$ks %||% c(10L, 25L, 50L, 100L, 200L, 500L)
    ks <- as.character(ks_raw)
    sel_modes <- exp$selection_modes %||% c("lasso", "portable", "portable_weighted")

    grid <- expand.grid(
      exposure_id    = exposures,
      k              = ks,
      selection_mode = sel_modes,
      stringsAsFactors = FALSE
    )
    grid$array_index      <- seq_len(nrow(grid))
    grid$module           <- "module6"
    grid$experiment_name  <- experiment_name
    grid$sub_workflow     <- "compact"
    grid$covariate_set    <- exp$covariate_set
    grid$sample_filter    <- exp$sample_filter %||% "none"
    grid$portability_threshold <- exp$portability_threshold %||% 0.70
    grid$portability_thresholds <- paste(
      exp$portability_thresholds %||% c(0.50, 0.60, 0.70, 0.80, 0.90),
      collapse = ","
    )
    grid$output_tag       <- exp$output_tag %||% "CompactPESThresholdSweep"
    grid$output_path      <- file.path(out_root, grid$exposure_id)
    grid$config_path      <- heap_config("modules", "module6_experiments.yml")

    rows <- grid[, c(
      "array_index", "module", "experiment_name", "sub_workflow",
      "covariate_set", "sample_filter", "exposure_id", "k", "selection_mode",
      "portability_threshold", "portability_thresholds",
      "output_tag", "output_path", "config_path"
    )]

  } else {
    rows <- data.frame(
      array_index      = seq_along(exposures),
      module           = "module6",
      experiment_name  = experiment_name,
      sub_workflow     = sub_wf,
      covariate_set    = exp$covariate_set,
      sample_filter    = exp$sample_filter %||% "none",
      exposure_id      = exposures,
      output_path      = file.path(out_root, exposures),
      config_path      = heap_config("modules", "module6_experiments.yml"),
      stringsAsFactors = FALSE
    )
  }

  .write_manifest(rows, "module6", experiment_name, overwrite = overwrite)
}

# ---------------------------------------------------------------------------
# Population architecture manifest
# ---------------------------------------------------------------------------

#' Generate the population architecture manifest for a named experiment.
#'
#' Each row represents one protein (or protein batch) × method job.
#'
#' @param experiment_name  Name from config/modules/population_architecture_experiments.yml
#' @param overwrite        Overwrite existing manifest if TRUE
generate_population_architecture_manifest <- function(experiment_name,
                                                       overwrite = TRUE) {
  exp <- load_experiment_config("population_architecture", experiment_name)

  proteins <- .load_protein_list(exp$protein_set %||% "omicspred_2704")

  out_root <- heap_project_output(
    exp$output_subdir %||% "population_architecture",
    experiment_name
  )

  rows <- data.frame(
    array_index          = seq_along(proteins),
    module               = "population_architecture",
    experiment_name      = experiment_name,
    protein_id           = proteins,
    covariate_set        = exp$covariate_set,
    method               = exp$method %||% "greml",
    grm_cutoff           = exp$grm_cutoff %||% 0.05,
    genetic_partition    = exp$genetic_partition %||% "total",
    exposure_kernel      = exp$exposure_kernel %||% "linear",
    gxe_mode             = exp$gxe_mode %||% "",
    sample_inclusion     = exp$sample_inclusion %||% "unrelated_european",
    cis_window_bp        = exp$cis_window_bp %||% 1000000L,
    min_protein_n        = exp$min_protein_n %||% 2000L,
    threads              = exp$threads %||% 8L,
    plink_memory_mb      = exp$plink_memory_mb %||% 96000L,
    exposure_missing_rate_max = exp$exposure_missing_rate_max %||% 0.20,
    config_path          = heap_config("modules", "population_architecture_experiments.yml"),
    output_path          = file.path(out_root, proteins),
    stringsAsFactors     = FALSE
  )

  .write_manifest(
    rows, "population_architecture", experiment_name, overwrite = overwrite
  )
}

# ---------------------------------------------------------------------------
# Convenience: generate all manifests for all ready experiments
# ---------------------------------------------------------------------------

#' Generate all manifests for experiments with status == "ready".
#'
#' Useful for a fresh run or after updating config files.
#'
#' @param modules  Which modules to generate for (default: all)
#' @param overwrite  Overwrite existing manifests
generate_all_ready_manifests <- function(
  modules = c("module1", "module2", "module3", "module5", "module6",
              "population_architecture"),
  overwrite = TRUE
) {
  results <- list()

  for (mod in modules) {
    tryCatch({
      cfg <- load_module_experiments(mod)
      exps <- cfg$experiments
      ready_names <- names(exps)[vapply(exps, function(e) {
        identical(e$status, "ready")
      }, logical(1))]

      if (length(ready_names) == 0L) {
        message("[manifest] No ready experiments for ", mod)
        next
      }

      gen_fn <- switch(
        mod,
        module1 = generate_module1_manifest,
        module2 = generate_module2_manifest,
        module3 = generate_module3_manifest,
        module5 = generate_module5_manifest,
        module6 = generate_module6_manifest,
        population_architecture = generate_population_architecture_manifest,
        NULL
      )

      if (is.null(gen_fn)) {
        message("[manifest] No generator for module: ", mod)
        next
      }

      for (nm in ready_names) {
        tryCatch(
          { path <- gen_fn(nm, overwrite = overwrite); results[[nm]] <- path },
          error = function(e) message("[manifest] FAILED ", nm, ": ", conditionMessage(e))
        )
      }
    }, error = function(e) {
      message("[manifest] Could not load experiments for ", mod, ": ", conditionMessage(e))
    })
  }

  invisible(results)
}
