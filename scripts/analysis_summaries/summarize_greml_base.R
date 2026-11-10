#!/usr/bin/env Rscript
# ============================================================================
# summarize_greml_base.R
# ----------------------------------------------------------------------------
# Aggregate the per-protein GREML summaries from the BASE covariate specification
# into one table for Supplementary Table S3.
#
# Source: population_architecture/base/grm_cutoff_0p025/primary/<PROT>_summary.tsv
#         (2,918 proteins attempted, covariate_spec = "base", grm_cutoff = 0.025)
#
# This is the run Figure 1b reads (via concordance_greml_vs_heap_r2.tsv) and the
# one the Methods describe ("Of 2,918 proteins, 2,238 converged and were
# retained"). The older aggregate at UK_Biobank/Output/GREMLcutoff/ was fit under
# the retired `Type3` spec -- which is NOT base (base drops BMI + fasting) -- and
# has 2,846 rows, matching neither the figure nor the Methods. Do not use it.
#
# All attempted proteins are kept, with `converged` as a column, so the 680
# non-converged proteins are visible rather than silently dropped
# (SUPPLEMENT_DATA_STANDARD.md §4).
#
#   Run:  Rscript scripts/analysis_summaries/summarize_greml_base.R
# ============================================================================
suppressPackageStartupMessages(library(data.table))
options(scipen = 999)

local({
  cand <- c(Sys.getenv("HEAP_PATHS_FILE", unset = ""),
            "/n/groups/patel/shakson_ukb/HEAP/workflow/00_paths.R")
  hit <- cand[nzchar(cand) & file.exists(cand)][1]; if (!is.na(hit)) source(hit)
})

SRC <- heap_project_output("population_architecture", "base", "grm_cutoff_0p025", "primary")
OUT <- file.path(heap_root, "docs", "manuscript_stats", "module1_variance")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

f <- list.files(SRC, pattern = "_summary\\.tsv$", full.names = TRUE)
if (!length(f)) stop("no per-protein GREML summaries under ", SRC)
message("Reading ", length(f), " per-protein GREML summaries (base spec) ...")
d <- rbindlist(lapply(f, fread), fill = TRUE)

## drop columns that are empty across every protein, and internal O2 paths
## (SUPPLEMENT_DATA_STANDARD.md §3 -- these must not ship)
internal <- c("hsq_path", "log_path", "prep_run_id")
allna <- names(d)[vapply(d, function(x)
  all(is.na(x) | trimws(as.character(x)) %in% c("", "NA")), logical(1))]
drop <- intersect(names(d), union(internal, allna))
if (length(drop)) message("  dropping: ", paste(drop, collapse = ", "))
d[, (drop) := NULL]

setorder(d, protein)
p <- file.path(OUT, "greml_per_protein_base.tsv")
fwrite(d, p, sep = "\t")
message(sprintf("  %s  %d proteins x %d cols  (%d converged, %d not)",
                basename(p), nrow(d), ncol(d),
                sum(d$converged %in% c(TRUE, "TRUE")),
                sum(!(d$converged %in% c(TRUE, "TRUE")))))
message("  covariate_spec: ", paste(unique(d$covariate_spec), collapse = ", "),
        " | grm_cutoff: ", paste(unique(d$grm_cutoff), collapse = ", "))
