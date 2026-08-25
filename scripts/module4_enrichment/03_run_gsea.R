#!/usr/bin/env Rscript
# ============================================================================
# 03_run_gsea.R  (PRIMARY / manuscript path)
#
# Purpose : Gene-Set Enrichment Analysis of Module 2 exposure->protein
#           associations. For every significant (replicated) exposure, rank its
#           proteins by t value and run:
#             - tissue GSEA  : clusterProfiler::GSEA over GTEx tissue Term2Gene
#             - pathway GSEA : ReactomePA::gsePathway over Entrez-ranked list
#           Results are bundled into a HEAPpathconstruct S4 object and saved.
#
# Inputs  :
#   - heap_project_output("module2","ReplicatedEassoc.csv")
#       (cols: ID, omicID, t value_test/Estimate_test, ...) — canonical sigBOTH.
#       Fallback (--covar-type): module2/<covarType>/univar_assoc_*.rds statE.
#   - heap_project_output("module4_enrichment","genesets","GTEX_tissue.txt")
#   - heap_project_output("module4_enrichment","OlinkEntrezConv.txt")  (Entrez map)
#
# Outputs :
#   - heap_project_output("module4_enrichment","HEAPgsea.qs")
#       S4 HEAPpathconstruct: slots HEAPtgsea (tissue GSEA @result per exposure),
#       HEAPpgsea (Reactome GSEA @result per exposure).
#
# Run     : Rscript 03_run_gsea.R [--covar-type Type3]
#
# Ported-from: scripts/visualizations/Visualizations/Module2/HEAPassoc_pathway.R
#   (PSEA_tissue, PSEA_paths, TissueEnrichGSEA, PathwayEnrichGSEA; tissue GSEA
#    over GTEx Term2Gene + Reactome gsePathway). Legacy GSEA params preserved:
#    minGSSize=10, maxGSSize=500, pvalueCutoff=0.05, pAdjustMethod="BH",
#    TERM2NAME=NA (tissue), eps=0 (Reactome). Ranking statistic = t value.
# ============================================================================

local({
  cand <- c(Sys.getenv("HEAP_PATHS_FILE", unset = ""),
            file.path(getwd(), "workflow", "00_paths.R"),
            "/n/groups/patel/shakson_ukb/HEAP/workflow/00_paths.R")
  hit <- cand[nzchar(cand) & file.exists(cand)][1]
  if (!is.na(hit)) source(hit)
})
if (!exists("heap_project_output"))
  stop("03_run_gsea.R: could not load workflow/00_paths.R (set HEAP_PATHS_FILE).")

# --- source the function library (robust script-dir resolution) --------------
.this_script_dir <- local({
  args <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", args[grepl("^--file=", args)])
  if (length(f)) dirname(normalizePath(f[1])) else getwd()
})
.core <- file.path(.this_script_dir, "02_enrichment_core.R")
if (!file.exists(.core))
  .core <- heap_script("module4_enrichment", "02_enrichment_core.R")
source(.core)
if (!isTRUE(.MODULE4_ENRICHMENT_CORE_LOADED))
  stop("03_run_gsea.R: failed to load 02_enrichment_core.R")

suppressPackageStartupMessages({
  library(qs)
  library(pbapply)
})

# Parallelise the per-exposure loop across forks (HEAP_GSEA_CORES, default 1) and
# force the inner clusterProfiler/fgsea to run SERIALLY so we don't nest forks.
gsea_cores <- max(1L, suppressWarnings(as.integer(Sys.getenv("HEAP_GSEA_CORES", "1"))))
gsea_cl <- if (gsea_cores > 1L) gsea_cores else NULL
if (requireNamespace("BiocParallel", quietly = TRUE))
  suppressWarnings(BiocParallel::register(BiocParallel::SerialParam()))
message("[gsea] per-exposure parallelism: ", gsea_cores, " core(s)")

OUT <- heap_project_output("module4_enrichment")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# --- parse optional --covar-type (selects per-batch rds fallback) ------------
.args <- commandArgs(trailingOnly = TRUE)
covar_type <- {
  i <- which(.args == "--covar-type")
  if (length(i) && length(.args) >= i + 1L) .args[i + 1L] else NA_character_
}

# --- load Module 2 associations ----------------------------------------------
# Explicit --input <csv> (e.g. the full-proteome ranked table from
# 00b_prepare_gsea_input.R) takes precedence; treated as prefiltered (GSEA ranks
# the full list, so the ranked vector is the whole proteome per exposure).
input_csv <- { i <- which(.args == "--input"); if (length(i) && length(.args) >= i + 1L) .args[i + 1L] else NA_character_ }
repl_csv <- heap_project_output("module2", "ReplicatedEassoc.csv")
if (!is.na(input_csv)) {
  message("[gsea] using explicit --input: ", input_csv)
  assoc <- load_replicated_assoc(input_csv)
  prefiltered <- TRUE
} else if (file.exists(repl_csv)) {
  message("[gsea] using canonical replicated associations: ", repl_csv)
  assoc <- load_replicated_assoc(repl_csv)
  prefiltered <- TRUE
} else if (!is.na(covar_type)) {
  message("[gsea] ReplicatedEassoc.csv missing; falling back to per-batch rds for covarType=", covar_type)
  assoc <- load_univar_statE(covar_type, split = "test")
  prefiltered <- FALSE
} else {
  stop("[module4_enrichment] No Module 2 input available.\n",
       "  Expected: ", repl_csv, "\n",
       "  Or pass --covar-type <Type> to use per-batch univar_assoc_*.rds.",
       call. = FALSE)
}

# Pre-load cached Entrez map once (shared by all pathway GSEA conversions).
df_entrez <- load_entrez_map()

# Exposures with at least one significant hit (sorted, unique).
Eid <- sort(unique(assoc$ID))
if (length(Eid) == 0L)
  stop("[module4_enrichment] no exposure IDs in the Module 2 table.", call. = FALSE)
message("[gsea] running GSEA across ", length(Eid), " exposures")

# --- tissue GSEA per exposure ------------------------------------------------
HEAPtissue <- pblapply(Eid, function(x) {
  suppressMessages(suppressWarnings({
    gl <- CreateGenelists(id = x, assoc_df = assoc, rank_stat = "t value",
                          prefiltered = prefiltered)
    tryCatch(PSEA_tissue(gl), error = function(e) {
      message("[gsea/tissue] ", x, ": ", conditionMessage(e)); NULL
    })
  }))
}, cl = gsea_cl)
names(HEAPtissue) <- Eid

# --- pathway (Reactome) GSEA per exposure ------------------------------------
HEAPpathway <- pblapply(Eid, function(x) {
  suppressMessages(suppressWarnings({
    gl <- CreateGenelists(id = x, assoc_df = assoc, rank_stat = "t value",
                          prefiltered = prefiltered)
    gl <- convertEntrez(gl, df_entrez = df_entrez)
    tryCatch(PSEA_paths(gl), error = function(e) {
      message("[gsea/pathway] ", x, ": ", conditionMessage(e)); NULL
    })
  }))
}, cl = gsea_cl)
names(HEAPpathway) <- Eid

# --- bundle into S4 object (matches legacy HEAPpathconstruct) -----------------
HEAPpathconstruct <- setClass(
  "HEAPpathconstruct",
  slots = c(
    HEAPtgsea = "list",   # tissue GSEA @result per exposure
    HEAPpgsea = "list"    # Reactome pathway GSEA @result per exposure
  )
)
HEAPgsea <- HEAPpathconstruct(HEAPtgsea = HEAPtissue, HEAPpgsea = HEAPpathway)

gc()
# Spec-suffixed output. The BASE spec deliberately keeps the original filename,
# because HEAPgsea.qs is the input from which tissue_enrichment.csv and
# pathway_enrichment.csv are derived -- and both of those are LIVE supplementary
# workbook tables in HEAP_manuscript/config/supp_tables.tsv (keys tissue_enrich
# and pathway_enrich, each cited once in the results). A second specification
# writing to the unsuffixed name would silently rewrite two published tables.
# Same convention Module 6 uses: base keeps the original names, everything else
# is suffixed.
.suffix <- {
  i <- which(.args == "--out-suffix")
  if (length(i) && length(.args) >= i + 1L) .args[i + 1L] else ""
}
out_file <- file.path(OUT, paste0("HEAPgsea", .suffix, ".qs"))
qsave(HEAPgsea, file = out_file)
message("[gsea] wrote ", out_file)
