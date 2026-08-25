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
if (length(args) < 2L) {
  stop("Usage: write_protein_set_from_loader.R <config.R> <out.txt>", call. = FALSE)
}

config_path <- args[1L]
out_path <- args[2L]

cfg <- load_config(config_path)
pxs_loader <- normalize_loader(read_pxs_loader(cfg$loader_rds))
protein_ids <- sort(unique(as.character(pxs_loader$protIDs)))
protein_ids <- protein_ids[nzchar(protein_ids)]

dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
writeLines(protein_ids, con = out_path)
message("Wrote ", length(protein_ids), " protein IDs to ", out_path)
