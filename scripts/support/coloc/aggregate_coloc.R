#!/usr/bin/env Rscript
# ============================================================================
# support/coloc/aggregate_coloc.R
# ----------------------------------------------------------------------------
# Collate all per-locus coloc.abf summaries from run_coloc_locus.R into a single
# results table, applying the PP.H4 >= 0.8 hard gate that replaces the Module-5
# "pending" coloc_status with a verdict:
#   confirmed  (PP.H4 >= 0.8)            -> shared causal variant; MR edge stands
#   distinct   (PP.H3 >= 0.8)            -> distinct causal variants; LD artifact
#   ambiguous  (neither passes)          -> underpowered / inconclusive
#   failed     (no coloc; status != ok)  -> too_few_snps / missing inputs etc.
#
# Reads:  output/support/coloc/per_locus/*_coloc_summary.tsv
# Writes: output/support/coloc/coloc_results.tsv   (one row per locus)
#         output/support/coloc/coloc_index.tsv     (locus -> PP.H4, for figures)
#
# NOTE: this script does NOT modify build_mr_tables.R or the MR sensitivity
# tables; it only produces the coloc verdict table that a downstream join can
# use to flip coloc_status from "pending".
#
# Usage:
#   module load gcc/14.2.0 R/4.4.2
#   HEAP_PATHS_FILE=.../workflow/00_paths.R \
#     Rscript scripts/support/coloc/aggregate_coloc.R
# ============================================================================
local({
  cand <- c(Sys.getenv("HEAP_PATHS_FILE", unset = ""),
            file.path(getwd(), "workflow", "00_paths.R"),
            "/n/groups/patel/shakson_ukb/HEAP/workflow/00_paths.R")
  hit <- cand[nzchar(cand) & file.exists(cand)][1]
  if (is.na(hit)) stop("Could not locate workflow/00_paths.R (set HEAP_PATHS_FILE).")
  source(hit)
})
suppressPackageStartupMessages({ library(data.table) })

PP4_CONFIRM <- as.numeric(Sys.getenv("COLOC_PP4", "0.8"))

out_dir   <- heap_project_output("support", "coloc")
locus_dir <- file.path(out_dir, "per_locus")
if (!dir.exists(locus_dir)) stop("per-locus dir not found (run run_coloc_locus.R first): ", locus_dir)

files <- list.files(locus_dir, pattern = "_coloc_summary\\.tsv$", full.names = TRUE)
if (!length(files)) stop("No *_coloc_summary.tsv under ", locus_dir)
message("Collating ", length(files), " per-locus coloc summaries.")

rows <- rbindlist(lapply(files, function(f) {
  s <- tryCatch(fread(f), error = function(e) NULL)
  if (is.null(s) || !nrow(s)) return(NULL)
  s <- s[1]
  getc <- function(col, default = NA) if (col %in% names(s)) s[[col]][1] else default
  data.table(
    arm      = getc("arm"),
    protID   = getc("protID"),
    disease  = getc("disease"),
    edge_dir = getc("edge_dir"),
    lead_snp = getc("lead_snp"),
    chr      = getc("chr"),
    pos      = getc("pos"),
    nsnps    = as.integer(getc("nsnps", NA_integer_)),
    PP.H3    = as.numeric(getc("PP.H3", NA_real_)),
    PP.H4    = as.numeric(getc("PP.H4", NA_real_)),
    run_status = as.character(getc("status", NA_character_))
  )
}), use.names = TRUE, fill = TRUE)

# verdict / hard gate
classify <- function(pp4, pp3, run_status) {
  if (is.na(pp4) || (!is.na(run_status) && run_status != "ok")) return("failed")
  if (pp4 >= PP4_CONFIRM) return("confirmed")
  if (!is.na(pp3) && pp3 >= PP4_CONFIRM) return("distinct")
  "ambiguous"
}
rows[, status := mapply(classify, PP.H4, PP.H3, run_status)]
setorder(rows, -PP.H4, na.last = TRUE)

results_fp <- file.path(out_dir, "coloc_results.tsv")
fwrite(rows[, .(arm, protID, disease, lead_snp, nsnps, PP.H3, PP.H4, status)], results_fp, sep = "\t")

# refresh coloc_index.tsv (locus -> PP.H4) for the figure layer; merge with any
# previously-staged legacy index rows that are not in this fresh batch.
idx_new <- rows[!is.na(PP.H4),
                .(locus = paste(arm, protID, disease, sep = "__"),
                  lead_snp, chr, pos, nsnps, PP.H3, PP.H4)]
idx_fp <- file.path(out_dir, "coloc_index.tsv")
if (file.exists(idx_fp)) {
  old <- tryCatch(fread(idx_fp), error = function(e) NULL)
  if (!is.null(old) && "locus" %in% names(old)) {
    keep_old <- old[!(locus %in% idx_new$locus)]
    idx_new <- rbind(idx_new, keep_old, fill = TRUE)
  }
}
fwrite(idx_new, idx_fp, sep = "\t")

# ---- report ----------------------------------------------------------------
cat("\n================ COLOC RESULTS ================\n")
cat("Loci collated:", nrow(rows), "  ->  ", results_fp, "\n\n")
brk <- rows[, .N, by = status][order(-N)]
print(brk)
cat(sprintf("\nconfirmed (PP.H4 >= %.2f): %d / %d\n",
            PP4_CONFIRM, sum(rows$status == "confirmed"), nrow(rows)))
cat("\nTop loci by PP.H4:\n")
print(head(rows[, .(arm, protID, disease, nsnps,
                    PP.H3 = round(PP.H3, 3), PP.H4 = round(PP.H4, 3), status)], 15))
cat("\nWrote:\n  ", results_fp, "\n  ", idx_fp, "\nDONE.\n", sep = "")
