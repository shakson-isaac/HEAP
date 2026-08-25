#!/usr/bin/env Rscript
# Reads analysis_plan.tsv and returns parameters for a given experiment_id.
# Used by SLURM scripts to avoid hardcoding run parameters.
#
# Usage (from a SLURM script or CLI):
#   Rscript read_analysis_plan.R M1_Type3_lasso
#   Rscript read_analysis_plan.R --list              # print all experiments
#   Rscript read_analysis_plan.R --list main          # print main-priority only

local({
  candidates <- c(
    Sys.getenv("HEAP_PATHS_FILE", unset = ""),
    file.path(getwd(), "workflow", "00_paths.R"),
    file.path(getwd(), "..", "workflow", "00_paths.R"),
    file.path(getwd(), "..", "..", "workflow", "00_paths.R")
  )
  candidates <- candidates[nzchar(candidates)]
  hit <- candidates[file.exists(candidates)][1]
  if (!is.na(hit)) source(hit)
})

plan_path <- heap_config("analysis_plan.tsv")

read_plan <- function() {
  if (!file.exists(plan_path))
    stop("analysis_plan.tsv not found at: ", plan_path)
  read.delim(plan_path, stringsAsFactors = FALSE, check.names = FALSE)
}

get_experiment <- function(experiment_id) {
  plan <- read_plan()
  row  <- plan[plan$experiment_id == experiment_id, ]
  if (nrow(row) == 0)
    stop("Unknown experiment_id: ", experiment_id,
         "\nAvailable: ", paste(plan$experiment_id, collapse = ", "))
  as.list(row[1, ])
}

list_experiments <- function(priority = NULL) {
  plan <- read_plan()
  if (!is.null(priority))
    plan <- plan[plan$priority == priority, ]
  plan
}

# ---- CLI entrypoint ----
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  cat("Usage: Rscript read_analysis_plan.R <experiment_id>\n")
  cat("       Rscript read_analysis_plan.R --list [priority]\n")
  quit(status = 0)
}

if (args[1] == "--list") {
  priority <- if (length(args) >= 2) args[2] else NULL
  plan <- list_experiments(priority)
  write.table(plan, stdout(), sep = "\t", quote = FALSE, row.names = FALSE)
} else {
  exp <- get_experiment(args[1])
  # Print as KEY=VALUE pairs suitable for shell eval
  for (nm in names(exp)) {
    cat(sprintf('%s="%s"\n', nm, exp[[nm]]))
  }
}
