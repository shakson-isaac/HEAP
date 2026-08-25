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
  stop(
    "Usage: filter_prep_inputs_by_keep.R <config.R> <run_id> <covar_spec> [--keep-file=path] [--center-exposures=true|false] [--cutoff=0.025]",
    call. = FALSE
  )
}

config_path <- args[1L]
run_id <- args[2L]
covar_spec_name <- args[3L]
opts <- parse_optional_args(args[-(1L:3L)])

cfg <- load_config(config_path)
center_exposures <- as_bool(get_opt(opts, "center_exposures", TRUE), default = TRUE)
exposure_mode <- if (center_exposures) "centered" else "uncentered"
keep_file <- get_opt(opts, "keep_file", "")
cutoff_value <- get_opt(opts, "cutoff", "")

if (!nzchar(keep_file) || !file.exists(keep_file)) {
  stopf("Missing keep file: %s", keep_file)
}

paths <- resolve_run_paths(cfg, run_id, covar_spec_name, exposure_mode = exposure_mode)
inputs_dir <- paths$inputs

read_keep <- function(path) {
  keep <- read_tsv(path, header = FALSE, col_classes = "character")
  if (ncol(keep) < 2L) {
    stopf("Keep file must have at least two columns (FID IID): %s", path)
  }
  keep <- keep[, 1:2, drop = FALSE]
  names(keep) <- c("FID", "IID")
  unique(keep)
}

filter_export <- function(path, keep_df) {
  df <- read_tsv(path, header = TRUE, col_classes = "character")
  if (!all(c("FID", "IID") %in% names(df))) {
    stopf("Expected FID and IID columns in %s.", path)
  }
  filtered <- merge(keep_df, df, by = c("FID", "IID"), all.x = FALSE, all.y = FALSE, sort = FALSE)
  write_tsv(filtered, path)
  nrow(filtered)
}

keep_df <- read_keep(keep_file)
sample_manifest_path <- file.path(inputs_dir, "sample_manifest.tsv")
sample_manifest <- read_tsv(sample_manifest_path, header = TRUE, col_classes = "character")
pre_n <- nrow(sample_manifest)

manifest_keep <- merge(sample_manifest[, c("FID", "IID", "eid"), drop = FALSE], keep_df, by = c("FID", "IID"), all.x = FALSE, all.y = FALSE, sort = FALSE)
write_tsv(manifest_keep, sample_manifest_path)
write_tsv(unique(manifest_keep[, c("FID", "IID"), drop = FALSE]), file.path(inputs_dir, "keep_ids.txt"), header = FALSE)

invisible(filter_export(file.path(inputs_dir, "covar_discrete.tsv"), keep_df))
invisible(filter_export(file.path(inputs_dir, "covar_quantitative.tsv"), keep_df))
invisible(filter_export(file.path(inputs_dir, "covar_kernel_core.tsv"), keep_df))
invisible(filter_export(file.path(inputs_dir, "exposure_raw.tsv"), keep_df))
post_n <- invisible(filter_export(file.path(inputs_dir, "protein_matrix.tsv"), keep_df))

metadata_path <- file.path(inputs_dir, "metadata.tsv")
if (file.exists(metadata_path)) {
  metadata_df <- read_tsv(metadata_path, header = TRUE)
  metadata_df$metric <- as.character(metadata_df$metric)
  metadata_df$value <- as.character(metadata_df$value)
  if ("n_samples" %in% metadata_df$metric) {
    metadata_df$value[metadata_df$metric == "n_samples"] <- as.character(post_n)
  } else {
    metadata_df <- rbind(metadata_df, data.frame(metric = "n_samples", value = as.character(post_n), stringsAsFactors = FALSE))
  }
  extra_metrics <- data.frame(
    metric = c("n_samples_before_relatedness_cutoff", "n_samples_after_relatedness_cutoff", "grm_relatedness_cutoff"),
    value = c(as.character(pre_n), as.character(post_n), as.character(cutoff_value)),
    stringsAsFactors = FALSE
  )
  metadata_df <- metadata_df[!metadata_df$metric %in% extra_metrics$metric, , drop = FALSE]
  metadata_df <- rbind(metadata_df, extra_metrics)
  write_tsv(metadata_df, metadata_path)
}

timestamp_msg("Filtered prep inputs from", pre_n, "to", post_n, "samples using", keep_file)
