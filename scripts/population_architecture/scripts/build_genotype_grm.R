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
  stop("Usage: build_genotype_grm.R <config.R> <run_id> <covar_spec> [--force=true] [--exposure-mode=centered]", call. = FALSE)
}

config_path <- args[1L]
run_id <- args[2L]
covar_spec_name <- args[3L]
opts <- parse_optional_args(args[-(1L:3L)])
cfg <- load_config(config_path)
force <- as_bool(get_opt(opts, "force", FALSE))
exposure_mode <- get_opt(opts, "exposure_mode", "centered")
threads <- resolve_threads(cfg)
plink_memory_mb <- suppressWarnings(as.integer(cfg$plink_memory_mb %||% NA_integer_))
paths <- resolve_run_paths(cfg, run_id, covar_spec_name, exposure_mode = exposure_mode)

keep_path <- file.path(paths$inputs, "keep_ids.txt")
if (!file.exists(keep_path)) {
  stopf("Missing keep file: %s. Run export_architecture_inputs.R first.", keep_path)
}

variant_ids_path <- file.path(paths$kernels, "ld_pruned_variant_ids.txt")
if (!file.exists(variant_ids_path) || force) {
  timestamp_msg("Reading LD-pruned variant IDs from", cfg$ld_pruned_pvar)
  pvar <- read_tsv(cfg$ld_pruned_pvar, header = TRUE, col_classes = "character")
  if (!"ID" %in% names(pvar)) {
    stopf("Expected an ID column in %s.", cfg$ld_pruned_pvar)
  }
  write.table(pvar$ID, file = variant_ids_path, row.names = FALSE, col.names = FALSE, quote = FALSE)
}

subset_prefix <- file.path(paths$kernels, "geno_ld_pruned_subset")
bed_target <- paste0(subset_prefix, ".bed")
if (!file.exists(bed_target) || force) {
  plink_log <- file.path(paths$logs, "plink_subset.log")
  if (is.finite(plink_memory_mb) && plink_memory_mb > 0L) {
    timestamp_msg(
      "Building LD-pruned genotype subset with",
      threads,
      "thread(s) and a PLINK2 memory cap of",
      plink_memory_mb,
      "MB."
    )
  } else {
    timestamp_msg("Building LD-pruned genotype subset with", threads, "thread(s).")
  }
  plink_args <- c(
    "--pfile", cfg$genotype_pfile,
    "--keep", keep_path,
    "--extract", variant_ids_path,
    "--make-bed",
    "--out", subset_prefix,
    "--threads", as.character(threads)
  )
  if (is.finite(plink_memory_mb) && plink_memory_mb > 0L) {
    plink_args <- c(plink_args, "--memory", as.character(plink_memory_mb))
  }
  run_command(
    cfg$plink2_bin,
    plink_args,
    log_path = plink_log
  )
}

grm_prefix <- file.path(paths$kernels, "geno_ld_pruned")
grm_target <- paste0(grm_prefix, ".grm.bin")
if (!file.exists(grm_target) || force) {
  gcta_log <- file.path(paths$logs, "gcta_make_grm.log")
  timestamp_msg("Constructing genotype GRM with", threads, "thread(s).")
  run_command(
    cfg$gcta_bin,
    c(
      "--bfile", subset_prefix,
      "--make-grm",
      "--thread-num", as.character(threads),
      "--out", grm_prefix
    ),
    log_path = gcta_log
  )
}

ids <- read_grm_ids(grm_prefix)
write_tsv(ids, file.path(paths$kernels, "geno_ld_pruned_ids.tsv"))
timestamp_msg("Genotype GRM ready:", grm_prefix)
