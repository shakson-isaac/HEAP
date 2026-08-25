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
  stop("Usage: write_protein_resource_tiers.R <config.R> <out_dir> [--small-threshold=N] [--medium-threshold=N]", call. = FALSE)
}

config_path <- args[1L]
out_dir <- args[2L]
opts <- parse_optional_args(args[-(1L:2L)])

small_threshold <- as.integer(get_opt(opts, "small_threshold", 47000L))
medium_threshold <- as.integer(get_opt(opts, "medium_threshold", 50000L))
if (small_threshold >= medium_threshold) {
  stopf("small_threshold (%s) must be smaller than medium_threshold (%s).", small_threshold, medium_threshold)
}

cfg <- load_config(config_path)
pxs_loader <- normalize_loader(read_pxs_loader(cfg$loader_rds))
prot_df <- as.data.frame(pxs_loader$UKBprot_df)
protein_ids <- sort(unique(as.character(pxs_loader$protIDs)))
protein_ids <- protein_ids[nzchar(protein_ids)]

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

nonmissing_counts <- vapply(
  protein_ids,
  function(protein_id) {
    if (!protein_id %in% names(prot_df)) {
      return(NA_integer_)
    }
    sum(is.finite(as.numeric(prot_df[[protein_id]])))
  },
  integer(1L)
)

resource_plan <- data.frame(
  protein = protein_ids,
  protein_nonmissing_n = as.integer(nonmissing_counts),
  tier = ifelse(
    nonmissing_counts < small_threshold,
    "small",
    ifelse(nonmissing_counts < medium_threshold, "medium", "large")
  ),
  stringsAsFactors = FALSE
)

write_tsv(resource_plan, file.path(out_dir, "resource_plan.tsv"))

for (tier_name in c("small", "medium", "large")) {
  tier_proteins <- resource_plan$protein[resource_plan$tier == tier_name]
  writeLines(tier_proteins, con = file.path(out_dir, paste0(tier_name, ".txt")))
}

message("Wrote resource plan for ", nrow(resource_plan), " proteins to ", out_dir)
