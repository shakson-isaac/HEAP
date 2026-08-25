#!/usr/bin/env Rscript
# prewarm_instruments.R — pre-compute (clump) each UNIQUE MR instrument exactly
# once and cache it to disk, so the main Module 5 array starts with warm caches
# and avoids redundant concurrent clumping (the same protein otherwise gets
# clumped by ~8 tasks at once on a cold run).
#
# It sources the relevant runner in WARM_ONLY mode (HEAP_MR_WARM_ONLY=1), which
# defines CFG + the get_*_instruments() functions but skips the edge run — so the
# cache it writes is byte-identical to what the runners read (no logic drift).
#
# Usage: Rscript prewarm_instruments.R <kind> <slice_idx> <n_slices>
#   kind:      ukb_protein | decode_protein | exposure | disease
#   slice_idx: 1..n_slices ; this task handles entities where
#              ((i-1) %% n_slices) + 1 == slice_idx
#
# Caches written (shared with the runners):
#   ukb_protein    -> output/mr/protein_inst/        (cis + trans)
#   decode_protein -> output/mr/protein_inst_decode/  (cis + trans)
#   exposure       -> output/mr/clumps/               (shared by both arms)
#   disease        -> output/mr/disease_inst/         (shared by both arms)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3)
  stop("Usage: Rscript prewarm_instruments.R <ukb_protein|decode_protein|exposure|disease> <slice_idx> <n_slices>")
kind      <- args[1]
slice_idx <- as.integer(args[2])
n_slices  <- as.integer(args[3])

# Source the runner for its functions ONLY (no edge run).
Sys.setenv(HEAP_MR_WARM_ONLY = "1")
heap_root <- Sys.getenv("HEAP_ROOT", unset = "/n/groups/patel/shakson_ukb/HEAP")
runner <- if (kind == "decode_protein")
  file.path(heap_root, "scripts", "module5_mr", "Module5_deCODE.R") else
  file.path(heap_root, "scripts", "module5_mr", "Module5.R")
if (!file.exists(runner)) stop("Runner not found: ", runner)
source(runner)

suppressPackageStartupMessages(library(data.table))

ge <- function(f) file.path(CFG$edges_dir, f)
rd <- function(f) if (file.exists(ge(f))) fread(ge(f), showProgress = FALSE) else data.table()

entities <- switch(kind,
  ukb_protein    = c(rd("edges_PD.tsv")$Protein,  rd("edges_PE.tsv")$Protein),
  decode_protein = c(rd("edges_PD.tsv")$Protein,  rd("edges_PE.tsv")$Protein),
  exposure       = c(rd("edges_EP.tsv")$Exposure, rd("edges_ED.tsv")$Exposure),
  disease        = c(rd("edges_DP.tsv")$Disease,  rd("edges_DE.tsv")$Disease),
  stop("unknown kind: ", kind)
)
entities <- sort(unique(entities[!is.na(entities) & nzchar(entities)]))
mine <- entities[(((seq_along(entities) - 1L) %% n_slices) + 1L) == slice_idx]

message(sprintf("[prewarm] kind=%s slice=%d/%d : %d of %d unique entities",
                kind, slice_idx, n_slices, length(mine), length(entities)))

warm_one <- function(e) {
  if (kind %in% c("ukb_protein", "decode_protein")) {
    get_protein_instruments(e, "cis")
    get_protein_instruments(e, "trans")
  } else if (kind == "exposure") {
    get_exposure_instruments(e)
  } else {
    get_disease_instruments(e)
  }
  "OK"
}

ok <- 0L; skip <- 0L
for (e in mine) {
  res <- tryCatch(warm_one(e), error = function(err) paste("SKIP:", conditionMessage(err)))
  if (identical(res, "OK")) ok <- ok + 1L
  else { skip <- skip + 1L; message("  ", e, " -> ", res) }
  gc(FALSE)  # release the (large, esp. ~900 MB deCODE genome-wide) GWAS before the next entity
}
message(sprintf("[prewarm] kind=%s slice=%d/%d DONE: %d cached, %d skipped",
                kind, slice_idx, n_slices, ok, skip))
