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

build_mgrm_file <- function(prefixes, path) {
  writeLines(prefixes, con = path)
  path
}

model_prefixes <- function(paths, model_name, exposure_mode) {
  prefixes <- list(
    primary = c(
      file.path(paths$kernels, "geno_ld_pruned"),
      file.path(paths$kernels, paste0("E_", exposure_mode)),
      file.path(paths$kernels, paste0("GxE_", exposure_mode))
    ),
    sensitivity = c(
      file.path(paths$kernels, "geno_ld_pruned"),
      file.path(paths$kernels, paste0("E_", exposure_mode)),
      file.path(paths$kernels, paste0("GxE_", exposure_mode)),
      file.path(paths$kernels, "CovarsxG"),
      file.path(paths$kernels, paste0("CovarsxE_", exposure_mode))
    )
  )
  prefixes[[model_name]]
}

parse_result_row <- function(hsq_path, log_path, protein_id, covar_spec_name, prep_run_id,
                             model_name,
                             exposure_mode, n_complete, fixed_covar_share,
                             predictive_cmp, keep_artifact_paths = TRUE,
                             grm_cutoff = NA_real_) {
  hsq_df <- parse_hsq(hsq_path)
  variance_g <- extract_hsq_value(hsq_df, "V(G1)/Vp")
  variance_e <- extract_hsq_value(hsq_df, "V(G2)/Vp")
  variance_gxe <- extract_hsq_value(hsq_df, "V(G3)/Vp")
  se_g <- extract_hsq_value(hsq_df, "V(G1)/Vp", "SE")
  se_e <- extract_hsq_value(hsq_df, "V(G2)/Vp", "SE")
  se_gxe <- extract_hsq_value(hsq_df, "V(G3)/Vp", "SE")

  variance_cg <- if (model_name == "sensitivity") extract_hsq_value(hsq_df, "V(G4)/Vp") else NA_real_
  variance_ce <- if (model_name == "sensitivity") extract_hsq_value(hsq_df, "V(G5)/Vp") else NA_real_
  se_cg <- if (model_name == "sensitivity") extract_hsq_value(hsq_df, "V(G4)/Vp", "SE") else NA_real_
  se_ce <- if (model_name == "sensitivity") extract_hsq_value(hsq_df, "V(G5)/Vp", "SE") else NA_real_

  data.frame(
    protein = protein_id,
    n = n_complete,
    covariate_spec = covar_spec_name,
    prep_run_id = prep_run_id,
    grm_cutoff = grm_cutoff,
    method = paste0("GCTA_multi_kernel_REML_", model_name),
    exposure_kernel_mode = exposure_mode,
    variance_G = variance_g,
    variance_E = variance_e,
    variance_GxE = variance_gxe,
    variance_Covars_fixed = fixed_covar_share,
    variance_CovarsxG_sensitivity = variance_cg,
    variance_CovarsxE_sensitivity = variance_ce,
    se_G = se_g,
    se_E = se_e,
    se_GxE = se_gxe,
    se_CovarsxG_sensitivity = se_cg,
    se_CovarsxE_sensitivity = se_ce,
    predictive_gxe_base_r2 = predictive_cmp$base_r2,
    predictive_gxe_full_r2 = predictive_cmp$full_r2,
    predictive_gxe_delta_r2 = predictive_cmp$delta_r2,
    converged = detect_convergence(log_path),
    warnings = collect_log_warnings(log_path),
    hsq_path = if (keep_artifact_paths) hsq_path else NA_character_,
    log_path = if (keep_artifact_paths) log_path else NA_character_,
    stringsAsFactors = FALSE
  )
}

make_failure_row <- function(protein_id, covar_spec_name, prep_run_id, model_name, exposure_mode,
                             n_complete, fixed_covar_share, predictive_cmp,
                             log_path, hsq_path = NA_character_, keep_artifact_paths = TRUE,
                             grm_cutoff = NA_real_) {
  data.frame(
    protein = protein_id,
    n = n_complete,
    covariate_spec = covar_spec_name,
    prep_run_id = prep_run_id,
    grm_cutoff = grm_cutoff,
    method = paste0("GCTA_multi_kernel_REML_", model_name),
    exposure_kernel_mode = exposure_mode,
    variance_G = NA_real_,
    variance_E = NA_real_,
    variance_GxE = NA_real_,
    variance_Covars_fixed = fixed_covar_share,
    variance_CovarsxG_sensitivity = NA_real_,
    variance_CovarsxE_sensitivity = NA_real_,
    se_G = NA_real_,
    se_E = NA_real_,
    se_GxE = NA_real_,
    se_CovarsxG_sensitivity = NA_real_,
    se_CovarsxE_sensitivity = NA_real_,
    predictive_gxe_base_r2 = predictive_cmp$base_r2,
    predictive_gxe_full_r2 = predictive_cmp$full_r2,
    predictive_gxe_delta_r2 = predictive_cmp$delta_r2,
    converged = FALSE,
    warnings = collect_log_warnings(log_path),
    hsq_path = if (keep_artifact_paths) hsq_path else NA_character_,
    log_path = if (keep_artifact_paths) log_path else NA_character_,
    stringsAsFactors = FALSE
  )
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3L) {
  stop("Usage: run_population_architecture.R <config.R> <run_id> <covar_spec> [--prep-run-id=id] [--model=primary|sensitivity] [--proteins=A,B] [--force=true] [--continue-on-error=true] [--center-exposures=true|false] [--write-combined-summary=true|false] [--keep-artifact-paths=true|false]", call. = FALSE)
}

config_path <- args[1L]
run_id <- args[2L]
covar_spec_name <- args[3L]
opts <- parse_optional_args(args[-(1L:3L)])

cfg <- load_config(config_path)
model_name <- get_opt(opts, "model", "primary")
if (!model_name %in% c("primary", "sensitivity")) {
  stopf("Unknown model name: %s", model_name)
}

force <- as_bool(get_opt(opts, "force", FALSE))
continue_on_error <- as_bool(get_opt(opts, "continue_on_error", FALSE))
center_exposures <- as_bool(get_opt(opts, "center_exposures", TRUE), default = TRUE)
write_combined_summary <- as_bool(get_opt(opts, "write_combined_summary", TRUE), default = TRUE)
keep_artifact_paths <- as_bool(get_opt(opts, "keep_artifact_paths", TRUE), default = TRUE)
threads <- resolve_threads(cfg)
min_protein_n <- as.integer(get_opt(opts, "min_protein_n", cfg$min_protein_n %||% 2000L))
grm_cutoff <- suppressWarnings(as.numeric(get_opt(opts, "grm_cutoff", NA_real_)))
exposure_mode <- if (center_exposures) "centered" else "uncentered"
paths <- resolve_run_paths(cfg, run_id, covar_spec_name, exposure_mode = exposure_mode)
prep_run_id <- get_opt(opts, "prep_run_id", run_id)
prep_output_root <- get_opt(opts, "prep_output_root", "")
prep_paths <- if (identical(prep_run_id, run_id)) {
  paths
} else if (nzchar(prep_output_root)) {
  resolve_run_paths_from_root(prep_output_root, prep_run_id, covar_spec_name, exposure_mode = exposure_mode)
} else {
  resolve_run_paths(cfg, prep_run_id, covar_spec_name, exposure_mode = exposure_mode)
}

protein_matrix <- read_tsv(file.path(prep_paths$inputs, "protein_matrix.tsv"), header = TRUE)
metadata_path <- file.path(prep_paths$inputs, "metadata.tsv")
metadata_df <- if (file.exists(metadata_path)) read_tsv(metadata_path, header = TRUE) else NULL
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
protein_names <- setdiff(names(protein_matrix), c("FID", "IID"))
proteins <- read_protein_arg(get_opt(opts, "proteins", NULL), default_proteins = protein_names)
proteins <- unique(gsub("-", "_", proteins))
protein_specific_prep <- identical(metadata_value("protein_specific_prep", "0"), "1")
target_protein <- metadata_value("target_protein", "")
if (protein_specific_prep) {
  if (length(proteins) != 1L) {
    stopf("Protein-specific prep at %s can only be used with one protein at a time.", prep_paths$inputs)
  }
  if (nzchar(target_protein) && !identical(proteins[1L], target_protein)) {
    stopf(
      "Requested protein %s does not match the protein-specific prep target %s.",
      proteins[1L],
      target_protein
    )
  }
}

mgrm_prefixes <- model_prefixes(prep_paths, model_name, exposure_mode)
missing_kernel <- mgrm_prefixes[!file.exists(paste0(mgrm_prefixes, ".grm.bin"))]
if (length(missing_kernel) > 0L) {
  stopf("Missing kernel GRMs for model %s: %s", model_name, paste(missing_kernel, collapse = ", "))
}

model_dir <- if (model_name == "primary") paths$models_primary else paths$models_sensitivity
grm_ids <- read_grm_ids(mgrm_prefixes[1L])
protein_aligned <- merge(grm_ids, protein_matrix, by = c("FID", "IID"), all.x = TRUE, sort = FALSE)

covar_discrete <- read_tsv(file.path(prep_paths$inputs, "covar_discrete.tsv"), header = TRUE)
covar_quant <- read_tsv(file.path(prep_paths$inputs, "covar_quantitative.tsv"), header = TRUE)
covar_disc_aligned <- merge(grm_ids, covar_discrete, by = c("FID", "IID"), all.x = TRUE, sort = FALSE)
covar_quant_aligned <- merge(grm_ids, covar_quant, by = c("FID", "IID"), all.x = TRUE, sort = FALSE)

use_disc <- ncol(covar_disc_aligned) > 2L
use_quant <- ncol(covar_quant_aligned) > 2L

summary_rows <- list()

for (protein_id in proteins) {
  if (!protein_id %in% names(protein_aligned)) {
    stopf("Protein %s is not present in the exported protein matrix.", protein_id)
  }

  out_prefix <- file.path(model_dir, protein_id)
  hsq_path <- paste0(out_prefix, ".hsq")
  log_path <- paste0(out_prefix, ".run.log")
  row_path <- paste0(out_prefix, "_summary.tsv")

  if (file.exists(row_path) && !force) {
    timestamp_msg("Skipping existing result for", protein_id)
    summary_rows[[protein_id]] <- read_tsv(row_path, header = TRUE)
    next
  }

  y <- as.numeric(protein_aligned[[protein_id]])
  n_complete <- sum(is.finite(y))
  if (n_complete < min_protein_n) {
    timestamp_msg("Skipping", protein_id, "because n =", n_complete, "is below threshold.")
    next
  }

  pheno <- data.frame(
    FID = grm_ids$FID,
    IID = grm_ids$IID,
    phenotype = ifelse(is.finite(y), y, -9),
    stringsAsFactors = FALSE
  )
  pheno_path <- paste0(out_prefix, ".phen")
  write_tsv(pheno, pheno_path, header = FALSE)
  mgrm_path <- build_mgrm_file(mgrm_prefixes, file.path(model_dir, paste0(protein_id, "_mgrm.txt")))
  disc_path <- file.path(model_dir, paste0(protein_id, "_covar_discrete.txt"))
  quant_path <- file.path(model_dir, paste0(protein_id, "_covar_quantitative.txt"))
  write_tsv(covar_disc_aligned, disc_path, header = FALSE)
  write_tsv(covar_quant_aligned, quant_path, header = FALSE)
  fixed_share <- compute_fixed_covar_share(
    y = y,
    discrete_df = covar_disc_aligned,
    quantitative_df = covar_quant_aligned
  )
  predictive_cmp <- compute_predictive_gxe_delta(
    protein_id = protein_id,
    covar_spec = covar_spec_name,
    predictive_root = cfg$predictive_oof_root %||% "",
    phenotype_df = data.frame(eid = grm_ids$IID, protein_aligned[, protein_id, drop = FALSE], check.names = FALSE)
  )

  gcta_args <- c(
    "--reml",
    "--mgrm", mgrm_path,
    "--pheno", pheno_path,
    "--thread-num", as.character(threads),
    "--out", out_prefix
  )
  if (is.finite(grm_cutoff) && !is.na(grm_cutoff) && grm_cutoff > 0) {
    gcta_args <- c(gcta_args, "--grm-cutoff", as.character(grm_cutoff))
  }
  if (use_quant) {
    gcta_args <- c(gcta_args, "--qcovar", quant_path)
  }
  if (use_disc) {
    gcta_args <- c(gcta_args, "--covar", disc_path)
  }

  run_ok <- TRUE
  tryCatch(
    {
      timestamp_msg("Running", model_name, "REML for", protein_id, "with", threads, "thread(s).")
      run_command(cfg$gcta_bin, gcta_args, log_path = log_path)
    },
    error = function(err) {
      run_ok <<- FALSE
      if (!continue_on_error) {
        stop(err)
      }
      timestamp_msg("GCTA failed for", protein_id, ":", conditionMessage(err))
    }
  )

  if (!run_ok || !file.exists(hsq_path)) {
    result_row <- make_failure_row(
      protein_id = protein_id,
      covar_spec_name = covar_spec_name,
      prep_run_id = prep_run_id,
      model_name = model_name,
      exposure_mode = exposure_mode,
      n_complete = n_complete,
      fixed_covar_share = fixed_share,
      predictive_cmp = predictive_cmp,
      log_path = log_path,
      hsq_path = if (file.exists(hsq_path)) hsq_path else NA_character_,
      keep_artifact_paths = keep_artifact_paths,
      grm_cutoff = grm_cutoff
    )
    write_tsv(result_row, row_path)
    summary_rows[[protein_id]] <- result_row
    next
  }

  result_row <- parse_result_row(
    hsq_path = hsq_path,
    log_path = log_path,
    protein_id = protein_id,
    covar_spec_name = covar_spec_name,
    prep_run_id = prep_run_id,
    model_name = model_name,
    exposure_mode = exposure_mode,
    n_complete = n_complete,
    fixed_covar_share = fixed_share,
    predictive_cmp = predictive_cmp,
    keep_artifact_paths = keep_artifact_paths,
    grm_cutoff = grm_cutoff
  )
  write_tsv(result_row, row_path)
  summary_rows[[protein_id]] <- result_row
}

if (length(summary_rows) == 0L) {
  stopf("No protein results were produced for model %s.", model_name)
}

if (write_combined_summary) {
  summary_df <- do.call(rbind, summary_rows)
  write_tsv(summary_df, file.path(paths$summary, paste0("per_protein_summary_", model_name, ".tsv")))
  timestamp_msg("Saved per-protein summary for", model_name, "to", paths$summary)
} else {
  timestamp_msg("Per-protein row outputs written for", length(summary_rows), "protein(s); combined summary skipped by request.")
}
