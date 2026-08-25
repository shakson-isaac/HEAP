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

# Export collaborator-friendly compact lasso protein weights from Module 6
# compact PES outputs. Files are plain tab-delimited .txt tables with protein
# names and beta weights, plus useful imputation metadata when available.

suppressPackageStartupMessages({
  library(data.table)
  library(glmnet)
})

cfg <- list(
  result_dir = heap_project_output("module6_pes_longitudinal"),
  out_dir = heap_project_output("module6_pes_longitudinal", "deployable_lasso_weights"),
  covar_type = "base",
  output_tag = "CompactPESThresholdSweep",
  ks = c("10", "25", "50", "100", "200", "500", "all")
)

`%||%` <- function(x, y) if (!is.null(x)) x else y

safe_name <- function(x) gsub("[^A-Za-z0-9_.-]+", "_", x)

parse_csv <- function(x) {
  if (is.null(x) || !nzchar(x)) return(character())
  out <- trimws(unlist(strsplit(x, ",")))
  out[nzchar(out)]
}

parse_args <- function(args) {
  out <- list()
  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) stop("Unexpected positional argument: ", key)
    if (i == length(args)) stop("Missing value for ", key)
    out[[gsub("-", "_", sub("^--", "", key))]] <- args[[i + 1]]
    i <- i + 2
  }
  out
}

extract_exposure_id <- function(panel_file, covar_type, tag) {
  base <- basename(panel_file)
  sub(paste0("^PESlong_", covar_type, "_"), "", sub(paste0("_", tag, "_ProteinPanels\\.tsv$"), "", base))
}

read_artifact_metadata <- function(artifact_file) {
  if (!file.exists(artifact_file)) {
    return(list(impute = data.table(protein = character(), impute_median = numeric()), intercept = NA_real_))
  }
  art <- readRDS(artifact_file)
  imp <- data.table(
    protein = names(art$protein_imputer$median),
    impute_median = as.numeric(art$protein_imputer$median)
  )
  intercept <- NA_real_
  if (!is.null(art$fit_prot) && !is.null(art$lambda$prot)) {
    cf <- as.matrix(coef(art$fit_prot, s = art$lambda$prot))
    if ("(Intercept)" %in% rownames(cf)) intercept <- as.numeric(cf["(Intercept)", 1])
  }
  list(impute = imp, intercept = intercept)
}

write_one_panel <- function(panel, out_file) {
  dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
  fwrite(panel, out_file, sep = "\t", quote = FALSE, na = "NA")
}

main <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))
  result_dir <- args$result_dir %||% cfg$result_dir
  out_dir <- args$out_dir %||% cfg$out_dir
  covar_type <- args$covar_type %||% cfg$covar_type
  tag <- args$output_tag %||% cfg$output_tag
  ks <- parse_csv(args$ks %||% paste(cfg$ks, collapse = ","))

  in_dir <- file.path(result_dir, covar_type)
  panel_files <- list.files(
    in_dir,
    pattern = paste0("_", tag, "_ProteinPanels\\.tsv$"),
    full.names = TRUE
  )
  if (length(panel_files) == 0) stop("No protein panel files found in ", in_dir, " for tag ", tag)

  manifest <- list()
  for (panel_file in panel_files) {
    exposure_id <- extract_exposure_id(panel_file, covar_type, tag)
    panel <- fread(panel_file)
    panel <- panel[selection_mode == "lasso"]
    if (nrow(panel) == 0) next

    artifact_file <- file.path(in_dir, paste0("PESlong_", covar_type, "_", safe_name(exposure_id), "_FinalModelArtifact.rds"))
    meta <- read_artifact_metadata(artifact_file)
    panel <- merge(panel, meta$impute, by = "protein", all.x = TRUE, sort = FALSE)

    for (kk in ks) {
      pp <- panel[requested_k == kk]
      if (nrow(pp) == 0 && kk == "all") pp <- panel[actual_k == max(actual_k, na.rm = TRUE)]
      if (nrow(pp) == 0) next
      pp <- pp[order(rank)]

      out_subdir <- file.path(out_dir, covar_type, safe_name(exposure_id))
      out_file <- file.path(out_subdir, paste0(safe_name(exposure_id), "_lasso_k", kk, "_weights.txt"))
      export <- pp[, .(
        protein,
        beta,
        rank,
        requested_k,
        actual_k,
        impute_median,
        portability_abs_corr,
        portability_match_method,
        portability_uniprot,
        portability_target,
        portability_gene_name
      )]
      write_one_panel(export, out_file)

      meta_file <- file.path(out_subdir, paste0(safe_name(exposure_id), "_lasso_k", kk, "_metadata.txt"))
      writeLines(c(
        paste0("exposure_id\t", exposure_id),
        paste0("covar_type\t", covar_type),
        paste0("source_panel_file\t", panel_file),
        paste0("source_artifact_file\t", artifact_file),
        paste0("selection_mode\tlasso"),
        paste0("requested_k\t", kk),
        paste0("actual_k\t", unique(export$actual_k)[1]),
        paste0("intercept_prot_model\t", meta$intercept),
        "score_formula\tlinear_score = intercept + sum((protein_value_imputed_if_missing) * beta)",
        "note\tWeights are from the final frozen prot-only glmnet PES model on the original protein scale. Use impute_median for missing protein values if reproducing this pipeline."
      ), meta_file)

      manifest[[length(manifest) + 1]] <- data.table(
        exposure_id = exposure_id,
        covar_type = covar_type,
        requested_k = kk,
        actual_k = unique(export$actual_k)[1],
        n_weights = nrow(export),
        weights_file = out_file,
        metadata_file = meta_file
      )
    }
  }

  manifest_tbl <- rbindlist(manifest, fill = TRUE)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  fwrite(manifest_tbl, file.path(out_dir, paste0("deployable_lasso_weights_manifest_", covar_type, "_", tag, ".tsv")), sep = "\t")
  message("Wrote ", nrow(manifest_tbl), " deployable lasso weight files under ", out_dir)
}

main()
