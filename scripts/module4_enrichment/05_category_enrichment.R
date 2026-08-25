#!/usr/bin/env Rscript
# ============================================================================
# 05_category_enrichment.R
# ----------------------------------------------------------------------------
# CATEGORY-level (directionless) GSEA for the exposures-by-tissue / by-pathway
# main figure. For each exposure CATEGORY, rank the proteome by the protein's
# STRONGEST |t| association across that category's exposures (UNSIGNED) and run
# tissue (GTEx) + Reactome GSEA. Ranking by |t| makes the enrichment
# DIRECTIONLESS: a positive NES = the category's strongly-associated proteins are
# over-represented in the tissue/pathway (no misleading up/down sign that would
# otherwise be an artifact of pooling opposite-direction exposures/levels).
#
# Input  : module4_enrichment/gsea_input_<experiment>.csv (per-exposure signed t,
#          from 00b_prepare_gsea_input.R) + analysis_exposures.tsv (category map)
# Output : module4_enrichment/{tissue,pathway}_enrichment_category.csv
#          (cols ID/Description, setSize, NES, p.adjust, cID = category)
#
# Run: HEAP_PATHS_FILE=.../00_paths.R Rscript 05_category_enrichment.R [experiment]
# ============================================================================
suppressPackageStartupMessages({ library(data.table) })
local({ cand <- c(Sys.getenv("HEAP_PATHS_FILE", unset = ""),
                  "/n/groups/patel/shakson_ukb/HEAP/workflow/00_paths.R")
        source(cand[file.exists(cand)][1]) })
.dir <- "/n/groups/patel/shakson_ukb/HEAP/scripts/module4_enrichment"
source(file.path(.dir, "02_enrichment_core.R"))
suppressPackageStartupMessages({ library(pbapply) })
if (requireNamespace("BiocParallel", quietly = TRUE))
  suppressWarnings(BiocParallel::register(BiocParallel::SerialParam()))

a <- commandArgs(trailingOnly = TRUE)
experiment <- if (length(a) >= 1) a[1] else "M2_base_main"
OUT <- heap_project_output("module4_enrichment")

# --- build category-level ranked input (|t| aggregated) ----------------------
gin <- fread(file.path(OUT, paste0("gsea_input_", experiment, ".csv")))
ae  <- fread(heap_config("exposure_sets", "analysis_exposures.tsv"))
gin[, category := setNames(ae$category, ae$variable)[ID]]
gin <- gin[!is.na(category)]
gin[, at := abs(`t value`)]
setorder(gin, category, omicID, -at)
catin <- gin[, .SD[1], by = .(category, omicID)]            # per (category, gene): peak |t|
catin <- catin[, .(ID = category, omicID, `t value` = at, Estimate = at, `Pr(>|t|)` = 1)]
cats <- sort(unique(catin$ID))
message("category-level GSEA over ", length(cats), " categories")

df_entrez <- load_entrez_map()
TI <- pblapply(cats, function(x) suppressMessages(suppressWarnings({
  gl <- CreateGenelists(id = x, assoc_df = catin, rank_stat = "t value", prefiltered = TRUE)
  tryCatch(PSEA_tissue(gl), error = function(e) NULL) })))
PA <- pblapply(cats, function(x) suppressMessages(suppressWarnings({
  gl <- CreateGenelists(id = x, assoc_df = catin, rank_stat = "t value", prefiltered = TRUE)
  gl <- convertEntrez(gl, df_entrez = df_entrez)
  tryCatch(PSEA_paths(gl), error = function(e) NULL) })))
names(TI) <- names(PA) <- cats

flat <- function(lst, cols) rbindlist(lapply(names(lst), function(cc) {
  d <- lst[[cc]]; if (is.null(d) || nrow(d) == 0L) return(NULL)
  d <- as.data.table(d); if (!all(cols %in% names(d))) return(NULL)
  out <- d[, ..cols]; out[, cID := cc]; out }), fill = TRUE)

fwrite(flat(TI, c("ID", "setSize", "NES", "p.adjust")),
       file.path(OUT, "tissue_enrichment_category.csv"))
fwrite(flat(PA, c("Description", "setSize", "NES", "p.adjust")),
       file.path(OUT, "pathway_enrichment_category.csv"))
cat("wrote tissue_enrichment_category.csv + pathway_enrichment_category.csv |",
    "tissue rows:", nrow(flat(TI, c("ID","setSize","NES","p.adjust"))),
    "| pathway rows:", nrow(flat(PA, c("Description","setSize","NES","p.adjust"))), "\n")
