#!/usr/bin/env Rscript
# ============================================================================
# 00b_prepare_gsea_input.R
# ----------------------------------------------------------------------------
# Build the GSEA input table for Module 2 exposure->protein associations from
# the clean (experiment-nested) Module 2 output. GSEA ranks the FULL proteome
# per exposure, so we keep every protein; multi-level (treatment-coded) exposures
# are collapsed to ONE representative association per (base exposure, protein) =
# the term with the largest |t| (the protein's strongest association with that
# exposure), so the downstream heatmap is per EXPOSURE, not per level. Restricted
# to exposures with >=1 replicated (train+test Bonferroni) E-block association.
#
# Gene names collapsed to first symbol (complexes) and de-duplicated per exposure
# (keep max |t|) so clusterProfiler::GSEA gets a tie-free ranked vector.
#
# Output: module4_enrichment/gsea_input_<experiment>.csv
#         (cols: ID = base exposure, omicID = gene symbol, `t value`, Estimate, `Pr(>|t|)`)
#
# Run: HEAP_PATHS_FILE=.../00_paths.R Rscript 00b_prepare_gsea_input.R [covarType] [experiment]
# ============================================================================
suppressPackageStartupMessages(library(data.table))
local({ cand <- c(Sys.getenv("HEAP_PATHS_FILE", unset = ""),
                  "/n/groups/patel/shakson_ukb/HEAP/workflow/00_paths.R")
        source(cand[file.exists(cand)][1]) })
CM <- "/n/groups/patel/shakson_ukb/HEAP/scripts/visualizations/common"
source(file.path(CM, "figure_paths.R")); source(file.path(CM, "load_heap_results.R"))

a <- commandArgs(trailingOnly = TRUE)
covarType  <- if (length(a) >= 1) a[1] else "base"
experiment <- if (length(a) >= 2) a[2] else "M2_base_main"

message("Loading Module 2 statE (train+test) for ", experiment, "/", covarType, " ...")
te <- load_module2_results(covarType, "test",  experiment = experiment)$statE
# active exposures = >=1 replicated (both-split Bonferroni) E-block association
rep <- load_module2_replicated(covarType, experiment = experiment)
active <- unique(rep[replicated == TRUE]$Eid)
message("active exposures (>=1 replicated): ", length(active), " of ", uniqueN(te$Eid))

st <- te[Eid %in% active & is.finite(`t value`),
         .(Eid, omicID, t = `t value`, Estimate, p = `Pr(>|t|)`)]
st[, symbol := sub("_.*", "", omicID)]                       # complexes -> first symbol
st[, at := abs(t)]
# per (exposure, gene symbol): keep the strongest-|t| association
setorder(st, Eid, symbol, -at)
st <- st[, .SD[1], by = .(Eid, symbol)]
out <- st[, .(ID = Eid, omicID = symbol, `t value` = t, Estimate = Estimate, `Pr(>|t|)` = p)]

outdir <- heap_project_output("module4_enrichment"); dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
outf <- file.path(outdir, paste0("gsea_input_", experiment, ".csv"))
fwrite(out, outf)
cat("wrote", outf, "|", nrow(out), "rows |", uniqueN(out$ID), "exposures |",
    uniqueN(out$omicID), "genes | median genes/exposure:",
    round(median(out[, .N, by = ID]$N)), "\n")
