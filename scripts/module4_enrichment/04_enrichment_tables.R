#!/usr/bin/env Rscript
# ============================================================================
# 04_enrichment_tables.R
#
# Purpose : Flatten the GSEA result object (HEAPgsea.qs) into tidy figure-input
#           CSVs — one row per (exposure x tissue/pathway) with the enrichment
#           statistics used by the Module 4 figures.
#
# Inputs  : heap_project_output("module4_enrichment","HEAPgsea.qs")
#             S4 HEAPpathconstruct: HEAPtgsea (tissue GSEA @result per exposure),
#             HEAPpgsea (Reactome GSEA @result per exposure).
#
# Outputs :
#   - heap_project_output("module4_enrichment","tissue_enrichment.csv")
#       (legacy EassocTissueEnrichment.csv shape: ID, setSize, NES, p.adjust, cID)
#   - heap_project_output("module4_enrichment","pathway_enrichment.csv")
#       (legacy EassocPathwayEnrichment.csv shape: Description, setSize, NES,
#        p.adjust, cID)
#   where cID = exposure id.
#
# Run     : Rscript 04_enrichment_tables.R
#
# Ported-from: scripts/visualizations/Visualizations/Module2/HEAPassoc_pathwaytables.R
#   (map2_dfr over HEAPtgsea/HEAPpgsea selecting setSize/NES/p.adjust + cID).
# ============================================================================

local({
  cand <- c(Sys.getenv("HEAP_PATHS_FILE", unset = ""),
            file.path(getwd(), "workflow", "00_paths.R"),
            "/n/groups/patel/shakson_ukb/HEAP/workflow/00_paths.R")
  hit <- cand[nzchar(cand) & file.exists(cand)][1]
  if (!is.na(hit)) source(hit)
})
if (!exists("heap_project_output"))
  stop("04_enrichment_tables.R: could not load workflow/00_paths.R (set HEAP_PATHS_FILE).")

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(purrr)
  library(qs)
})

OUT <- heap_project_output("module4_enrichment")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# See 03_run_gsea.R: base keeps the unsuffixed names because tissue_enrichment
# and pathway_enrichment ship as cited supplementary tables.
.args <- commandArgs(trailingOnly = TRUE)
.suffix <- {
  i <- which(.args == "--out-suffix")
  if (length(i) && length(.args) >= i + 1L) .args[i + 1L] else ""
}
gsea_file <- file.path(OUT, paste0("HEAPgsea", .suffix, ".qs"))
if (!file.exists(gsea_file))
  stop("[module4_enrichment] missing GSEA object:\n  ", gsea_file,
       "\nRun 03_run_gsea.R first.", call. = FALSE)

HEAPht <- qread(gsea_file)

# Flatten one named list of @result data.frames -> tidy data.frame, keeping
# `keep_cols` and tagging the exposure id as cID. Skips NULL/empty/missing-col
# entries (an exposure can yield no enriched terms).
flatten_enrichment <- function(res_list, keep_cols) {
  rows <- imap(res_list, function(df, exposure_id) {
    if (is.null(df) || nrow(df) == 0L) return(NULL)
    if (!all(keep_cols %in% names(df))) return(NULL)
    df %>%
      select(all_of(keep_cols)) %>%
      mutate(cID = exposure_id)
  })
  out <- bind_rows(Filter(Negate(is.null), rows))
  rownames(out) <- NULL
  out
}

# Both tables now keep the gene-set ID *and* its Description, plus the raw
# `pvalue` alongside the adjusted one. The legacy selection kept only one
# identifier each -- tissue had a slug ID with no label, pathway had a label with
# no Reactome ID, so neither joined cleanly -- and dropped the raw p entirely,
# which the supplement standard requires reported alongside any adjusted value.
# Additive: downstream figures select columns by name and are unaffected.
KEEP <- c("ID", "Description", "setSize", "NES", "pvalue", "p.adjust")

# Tissue effects (GTEx term key = ID).
TissueEffects <- flatten_enrichment(HEAPht@HEAPtgsea, KEEP)

# Pathway effects (Reactome).
PathEffects <- flatten_enrichment(HEAPht@HEAPpgsea, KEEP)

tissue_out  <- file.path(OUT, paste0("tissue_enrichment", .suffix, ".csv"))
pathway_out <- file.path(OUT, paste0("pathway_enrichment", .suffix, ".csv"))
fwrite(TissueEffects, tissue_out)
fwrite(PathEffects,  pathway_out)
message("[tables] wrote ", tissue_out,  " (", nrow(TissueEffects), " rows)")
message("[tables] wrote ", pathway_out, " (", nrow(PathEffects),  " rows)")
