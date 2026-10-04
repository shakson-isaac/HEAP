#!/usr/bin/env Rscript
# ============================================================================
# support/coloc/run_coloc.R
# ----------------------------------------------------------------------------
# Stage colocalization results into the canonical IGLOO support/coloc/ location
# so the figure layer (fig_mr_coloc) can read them with the usual path discipline,
# and write a coloc index (locus -> PP.H4) for prioritization.
#
# The heavy coloc analysis itself (coloc.abf / SuSiE fine-mapping + LocusZoom LD
# from the harmonized pQTL x disease GWAS) lives in the legacy
#   ModuleMR/COLOC/{Susie_Coloc.R, LocusZoom.R}
# which already produced per-locus result sets under
#   UK_Biobank/Output/Coloc/Results/<lead>_<protein>_<disease>_{plot_table,coloc_summary,
#                                    susie_susie_coloc_summary,harmonized_snps}.tsv
# This script copies those (the plotted-data tables) to the canonical tree; a full
# in-repo re-run of the SuSiE/LocusZoom port (PLINK EUR LD + coloc) is the
# regeneration path and remains TODO (needs the LD ref + GWAS staged to IGLOO).
#
# Usage:
#   module load gcc/14.2.0 R/4.4.2
#   HEAP_PATHS_FILE=.../workflow/00_paths.R \
#     Rscript scripts/support/coloc/run_coloc.R [<legacy_results_dir>]
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

args     <- commandArgs(trailingOnly = TRUE)
src_dir  <- if (length(args) >= 1) args[1] else
  "/n/groups/patel/shakson_ukb/UK_Biobank/Output/Coloc/Results"
out_dir  <- heap_project_output("support", "coloc")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(src_dir)) stop("legacy coloc Results dir not found: ", src_dir)

# loci = stems with a *_plot_table.tsv (the plotted-data table the figure needs)
pt <- list.files(src_dir, pattern = "_plot_table\\.tsv$", full.names = TRUE)
if (!length(pt)) stop("No *_plot_table.tsv under ", src_dir)
loci <- sub("_plot_table\\.tsv$", "", basename(pt))
message("Found ", length(loci), " coloc locus result set(s): ", paste(loci, collapse = ", "))

WANT <- c("plot_table", "coloc_summary", "susie_susie_coloc_summary", "harmonized_snps")
idx  <- rbindlist(lapply(loci, function(locus) {
  for (w in WANT) {
    f <- file.path(src_dir, paste0(locus, "_", w, ".tsv"))
    if (file.exists(f)) file.copy(f, file.path(out_dir, basename(f)), overwrite = TRUE)
  }
  cs <- file.path(src_dir, paste0(locus, "_coloc_summary.tsv"))
  row <- data.table(locus = locus)
  if (file.exists(cs)) {
    s <- fread(cs)
    for (cc in intersect(c("lead_snp","chr","pos","nsnps","PP.H3","PP.H4"), names(s)))
      row[[cc]] <- s[[cc]][1]
  }
  row
}), fill = TRUE)
fwrite(idx, file.path(out_dir, "coloc_index.tsv"), sep = "\t")

cat("\nStaged coloc results -> ", out_dir, "\n", sep = "")
print(idx)
cat("DONE.\n")
