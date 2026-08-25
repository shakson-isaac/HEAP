#!/usr/bin/env Rscript
# ============================================================================
# 03b_run_ora.R  (SEPARATE / optional functionality — not the figure default)
#
# Purpose : Over-Representation Analysis (ORA) of Module 2 significant
#           exposure->protein associations. For every exposure, take its
#           significant protein set (split complexes by underscore) and test:
#             - tissue ORA  : clusterProfiler::enricher over GTEx Term2Gene
#             - pathway ORA : KEGG enricher + enrichGO(BP/CC/MF) +
#                             enrichPathway (Reactome) + enrichDO
#           Universe = OmicsPred proteins surveyed in UK Biobank.
#
# Inputs  :
#   - heap_project_output("module2","ReplicatedEassoc.csv")  (canonical sigBOTH)
#       Fallback (--covar-type): module2/<covarType>/univar_assoc_*.rds statE.
#   - genesets/GTEX_tissue.txt, genesets/KEGG_T2G.txt, genesets/KEGG_T2N.txt
#   - OlinkEntrezConv.txt  (cached Entrez map)
#
# Outputs :
#   - heap_project_output("module4_enrichment","ora_tissue.csv")
#   - heap_project_output("module4_enrichment","ora_pathway.csv")
#
# Run     : Rscript 03b_run_ora.R [--covar-type Type3] [--direction all|up|down]
#
# Ported-from: UK_Biobank/RScripts/Pure_StatGen/TissueSpec/TissueSpec_RunExample.R
#   (ORA_tissue, ORA_paths). Legacy ORA params preserved: pvalueCutoff=0.05,
#   pAdjustMethod="BH", qvalueCutoff=0.2, minGSSize=10, maxGSSize=500, universe =
#   OmicsPred (symbols for tissue, Entrez for pathways), per-enrichment tryCatch.
# ============================================================================

local({
  cand <- c(Sys.getenv("HEAP_PATHS_FILE", unset = ""),
            file.path(getwd(), "workflow", "00_paths.R"),
            "/n/groups/patel/shakson_ukb/HEAP/workflow/00_paths.R")
  hit <- cand[nzchar(cand) & file.exists(cand)][1]
  if (!is.na(hit)) source(hit)
})
if (!exists("heap_project_output"))
  stop("03b_run_ora.R: could not load workflow/00_paths.R (set HEAP_PATHS_FILE).")

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
  stop("03b_run_ora.R: failed to load 02_enrichment_core.R")

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(pbapply)
})

OUT <- heap_project_output("module4_enrichment")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# --- args --------------------------------------------------------------------
.args <- commandArgs(trailingOnly = TRUE)
.opt <- function(flag, default = NA_character_) {
  i <- which(.args == flag)
  if (length(i) && length(.args) >= i + 1L) .args[i + 1L] else default
}
covar_type <- .opt("--covar-type")
direction  <- .opt("--direction", "all")        # all | up | down
ora_set     <- paste0("ORA_", direction)         # ORA_all / ORA_up / ORA_down
ora_set_ez  <- paste0(ora_set, "_entrez")

# --- load Module 2 significant associations ----------------------------------
repl_csv <- heap_project_output("module2", "ReplicatedEassoc.csv")
if (file.exists(repl_csv)) {
  message("[ora] using canonical replicated associations: ", repl_csv)
  assoc <- load_replicated_assoc(repl_csv)
  prefiltered <- TRUE
} else if (!is.na(covar_type)) {
  message("[ora] ReplicatedEassoc.csv missing; per-batch rds for covarType=", covar_type)
  assoc <- load_univar_statE(covar_type, split = "test")
  prefiltered <- FALSE
} else {
  stop("[module4_enrichment] No Module 2 input available.\n",
       "  Expected: ", repl_csv, "\n",
       "  Or pass --covar-type <Type> to use per-batch univar_assoc_*.rds.",
       call. = FALSE)
}

df_entrez <- load_entrez_map()
Eid <- sort(unique(assoc$ID))
if (length(Eid) == 0L)
  stop("[module4_enrichment] no exposure IDs in the Module 2 table.", call. = FALSE)
message("[ora] running ORA (", ora_set, ") across ", length(Eid), " exposures")

# --- tissue ORA per exposure -------------------------------------------------
tissue_rows <- pblapply(Eid, function(x) {
  suppressMessages(suppressWarnings({
    gl <- CreateGenelists(id = x, assoc_df = assoc, prefiltered = prefiltered)
    if (length(gl[[ora_set]]) == 0L) return(NULL)
    res <- tryCatch(ORA_tissue(gl, Type = ora_set), error = function(e) {
      message("[ora/tissue] ", x, ": ", conditionMessage(e)); NULL
    })
    if (is.null(res) || nrow(res) == 0L) return(NULL)
    res$ID_exposure <- x
    res$direction   <- direction
    res
  }))
})
tissue_df <- data.table::rbindlist(Filter(Negate(is.null), tissue_rows), fill = TRUE)

# --- pathway ORA per exposure (KEGG/GO/Reactome/DO) --------------------------
pathway_rows <- pblapply(Eid, function(x) {
  suppressMessages(suppressWarnings({
    gl <- CreateGenelists(id = x, assoc_df = assoc, prefiltered = prefiltered)
    gl <- convertEntrez(gl, df_entrez = df_entrez)
    if (length(gl[[ora_set_ez]]) == 0L || all(is.na(gl[[ora_set_ez]]))) return(NULL)
    enr <- tryCatch(ORA_paths(gl, Type = ora_set_ez, df_entrez = df_entrez),
                    error = function(e) { message("[ora/path] ", x, ": ", conditionMessage(e)); NULL })
    if (is.null(enr)) return(NULL)
    # Flatten the named list of @result data.frames, tagging the source database.
    parts <- lapply(names(enr), function(db) {
      r <- enr[[db]]
      if (is.null(r) || nrow(r) == 0L) return(NULL)
      r$database    <- db
      r$ID_exposure <- x
      r$direction   <- direction
      r
    })
    data.table::rbindlist(Filter(Negate(is.null), parts), fill = TRUE)
  }))
})
pathway_df <- data.table::rbindlist(Filter(Negate(is.null), pathway_rows), fill = TRUE)

# --- write -------------------------------------------------------------------
tissue_out  <- file.path(OUT, "ora_tissue.csv")
pathway_out <- file.path(OUT, "ora_pathway.csv")
fwrite(tissue_df,  tissue_out)
fwrite(pathway_df, pathway_out)
message("[ora] wrote ", tissue_out, " (", nrow(tissue_df), " rows)")
message("[ora] wrote ", pathway_out, " (", nrow(pathway_df), " rows)")
