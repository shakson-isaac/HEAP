#!/usr/bin/env Rscript
# ============================================================================
# 02_enrichment_core.R
# Module 4 enrichment — sourced function library.
#
# Purpose : Shared helpers for tissue/pathway enrichment of Module 2
#           exposure->protein associations. Provides:
#             - CreateGenelists()      build per-exposure ranked GSEA vector + ORA gene sets
#             - load_entrez_map()      read cached BioMart map (OlinkEntrezConv.txt)
#             - convertEntrez()        attach Entrez IDs to gene lists (GSEA + ORA)
#             - convertEntreztoSymbol()lookup symbol(s) for an Entrez ID
#             - ukb_universe()         OmicsPred protein universe (symbols + entrez)
#             - PSEA_tissue()          GSEA over GTEx tissue Term2Gene  (primary)
#             - PSEA_paths()           Reactome gsePathway over Entrez ranks (primary)
#             - ORA_tissue()           enricher over GTEx tissue Term2Gene (separate)
#             - ORA_paths()            KEGG/GO/Reactome/DO over-representation (separate)
#
# Inputs (resolved via workflow/00_paths.R helpers; consumed by 03_/03b_/04_):
#   - Module 2 significant assoc : heap_project_output("module2","ReplicatedEassoc.csv")
#                                  (cols: ID, omicID, <stat>_train, <stat>_test, AssocID)
#                                  fallback: per-batch module2/<covarType>/univar_assoc_*.rds
#   - Gene sets                  : heap_project_output("module4_enrichment","genesets",*.txt)
#                                  GTEX_tissue.txt, KEGG_T2G.txt, KEGG_T2N.txt, HPA_*.txt
#   - Entrez map                 : heap_project_output("module4_enrichment","OlinkEntrezConv.txt")
#   - Protein universe           : heap_omicspred_or_legacy(<OmicsPred portal csv>)
#
# Outputs : none (function library)
#
# Ported-from:
#   scripts/visualizations/Visualizations/Module2/HEAPassoc_pathway.R
#       (CreateGenelists, PSEA_tissue, PSEA_paths, convertEntrez, entrezMap)
#   UK_Biobank/RScripts/Pure_StatGen/TissueSpec/TissueSpec_RunExample.R
#       (ORA_tissue, ORA_paths)
#
# Notes on the legacy->new port:
#   * Legacy read `HEAPassoc@HEAPlist[[CovarType]][[Split]][[1]]` (statE, all assoc)
#     then Bonferroni-filtered per exposure. The new canonical input
#     ReplicatedEassoc.csv is already the replicated-in-both (sigBOTH) table, so
#     every row is significant: the per-exposure Bonferroni filter is no longer
#     re-applied (it was applied upstream in summarize_replicated_associations.R).
#   * Legacy ran enrichment on the TEST split; we keep that by ranking on the
#     test-split t value / Estimate ("<stat>_test") when present, else the bare
#     column name (per-batch rds path or legacy single-split table).
#   * Legacy live BioMart `entrezMap()` is replaced by the cached
#     OlinkEntrezConv.txt produced by 01_entrez_map.R (one BioMart call).
# ============================================================================

# --- resolve path helpers ----------------------------------------------------
local({
  cand <- c(Sys.getenv("HEAP_PATHS_FILE", unset = ""),
            file.path(getwd(), "workflow", "00_paths.R"),
            "/n/groups/patel/shakson_ukb/HEAP/workflow/00_paths.R")
  hit <- cand[nzchar(cand) & file.exists(cand)][1]
  if (!is.na(hit)) source(hit)
})

if (!exists("heap_project_output"))
  stop("02_enrichment_core.R: could not load workflow/00_paths.R ",
       "(set HEAP_PATHS_FILE or run from inside the HEAP tree).")

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(stringr)
})

# Canonical Module 4 output dir + sub-paths
.m4_out      <- function(...) heap_project_output("module4_enrichment", ...)
.m4_geneset  <- function(name) .m4_out("genesets", name)
.m4_entrez   <- function() .m4_out("OlinkEntrezConv.txt")

# OmicsPred universe portal CSV (canonical IGLOO, else legacy full path).
.omicspred_universe_file <- function() {
  heap_omicspred_or_legacy("UKB_Olink_multi_ancestry_models_val_results_portal.csv")
}

# Informative existence check used throughout the module.
.require_input <- function(path, what) {
  if (!file.exists(path))
    stop(sprintf(
      "[module4_enrichment] missing %s:\n  %s\nRun the upstream step that produces it (see README.md).",
      what, path), call. = FALSE)
  invisible(path)
}

# ----------------------------------------------------------------------------
# Module 2 input loading
#
# Canonical: heap_project_output("module2","ReplicatedEassoc.csv") — the
# replicated-in-both (sigBOTH) exposure->protein associations. Columns:
#   ID, omicID, <stat>_train, <stat>_test, AssocID
# Fallback : per-batch module2/<covarType>/univar_assoc_*.rds (list(train,test)
#   each a list whose [[1]] element is statE with columns
#   ID, Eid, Category, Estimate, Std. Error, t value, Pr(>|t|), ..., omicID).
# ----------------------------------------------------------------------------

# Resolve the per-exposure stat columns for a given metric, preferring the
# test-split copy (legacy ran enrichment on TEST), then the plain column.
.pick_col <- function(df, base) {
  nm <- names(df)
  cand <- c(paste0(base, "_test"), base, paste0(base, "_train"))
  hit <- cand[cand %in% nm][1]
  if (is.na(hit))
    stop(sprintf("[module4_enrichment] expected a '%s' column (or %s_test) in the Module 2 table; found: %s",
                 base, base, paste(nm, collapse = ", ")), call. = FALSE)
  hit
}

# Resolve a Module 2 output path IGLOO-first, then local HEAP/output staging
# (consistent with summarize_replicated_associations.R / run_intervention_compare.R).
.m2_resolve <- function(...) {
  ig <- heap_project_output("module2", ...)
  if (file.exists(ig) || dir.exists(ig)) return(ig)
  loc <- heap_output("module2", ...)
  if (file.exists(loc) || dir.exists(loc)) return(loc)
  ig                                     # canonical path for the error message
}

# Load the canonical replicated-association table (or a provided data.frame).
load_replicated_assoc <- function(path = .m2_resolve("ReplicatedEassoc.csv")) {
  .require_input(path, "Module 2 replicated associations (ReplicatedEassoc.csv)")
  df <- as.data.frame(fread(path), stringsAsFactors = FALSE, check.names = FALSE)
  if (!all(c("ID", "omicID") %in% names(df)))
    stop("[module4_enrichment] ReplicatedEassoc.csv must contain 'ID' and 'omicID' columns; found: ",
         paste(names(df), collapse = ", "), call. = FALSE)
  df
}

# Fallback loader: aggregate ALL per-batch univar_assoc_*.rds statE (component
# [[1]]) for one covarType + split. Files are named by protein index (e.g.
# univar_assoc_2113.rds), so we glob rather than assume a 1..N index range.
load_univar_statE <- function(covarType, split = "test", batch_dir = NULL) {
  if (is.null(batch_dir)) batch_dir <- .m2_resolve(covarType)
  if (!dir.exists(batch_dir))
    stop("[module4_enrichment] univar batch dir not found: ", batch_dir,
         " (run Module 2 for ", covarType, ").", call. = FALSE)
  files <- list.files(batch_dir, pattern = "^univar_assoc_.*\\.rds$", full.names = TRUE)
  parts <- list()
  for (f in files) {
    obj <- tryCatch(readRDS(f), error = function(e) NULL)
    if (is.null(obj) || is.null(obj[[split]])) next
    statE <- obj[[split]][[1]]            # component [[1]] = statE
    if (!is.null(statE) && nrow(statE) > 0) parts[[length(parts) + 1L]] <- statE
  }
  if (length(parts) == 0L)
    stop("[module4_enrichment] no univar_assoc statE rows found under ", batch_dir,
         " (split=", split, ").", call. = FALSE)
  data.table::rbindlist(parts, fill = TRUE) |> as.data.frame(check.names = FALSE)
}

# ----------------------------------------------------------------------------
# CreateGenelists()
#
# Build, for one exposure id, the gene lists used downstream:
#   GSEA   : named numeric vector (names = gene symbol, value = ranking stat),
#            sorted decreasing; complexes split to first symbol to avoid ties.
#   Betas  : named numeric vector of effect sizes (same naming).
#   ORA_all/ORA_up/ORA_down : character vectors of significant gene symbols.
#
# `assoc_df` is the Module 2 table (canonical ReplicatedEassoc.csv = all rows
# already significant, or a per-exposure subset). `rank_stat` selects the GSEA
# ranking metric ("t value" in legacy). When `prefiltered = TRUE` (canonical)
# every row is treated as significant; otherwise a Bonferroni filter on the
# p-value column is applied as in the legacy CreateGenelists().
# ----------------------------------------------------------------------------
CreateGenelists <- function(id, assoc_df, rank_stat = "t value",
                            prefiltered = TRUE) {
  GeneLists <- list()

  Eassoc_df <- assoc_df %>% filter(.data[["ID"]] == id)
  if (nrow(Eassoc_df) == 0L)
    warning("[module4_enrichment] no associations for exposure id: ", id)

  estimate_col <- .pick_col(Eassoc_df, "Estimate")
  rank_col     <- .pick_col(Eassoc_df, rank_stat)

  # Significant rows for ORA. Canonical input is already sigBOTH (prefiltered);
  # otherwise reproduce the legacy per-exposure Bonferroni filter.
  if (prefiltered) {
    sig_df <- Eassoc_df
  } else {
    p_col <- .pick_col(Eassoc_df, "Pr(>|t|)")
    bonf  <- 0.05 / nrow(Eassoc_df)
    sig_df <- Eassoc_df %>% filter(.data[[p_col]] < bonf)
  }

  split_complex <- function(x) unlist(strsplit(x, "_"))  # split complexes by protein

  ORA_all  <- sig_df %>% pull(omicID) %>% split_complex()
  ORA_up   <- sig_df %>% filter(.data[[estimate_col]] > 0) %>% pull(omicID) %>% split_complex()
  ORA_down <- sig_df %>% filter(.data[[estimate_col]] < 0) %>% pull(omicID) %>% split_complex()

  # GSEA ranked vector: use the ranking statistic (t value) — significance+direction.
  GSEAgenelist <- Eassoc_df[[rank_col]]
  names(GSEAgenelist) <- Eassoc_df$omicID
  names(GSEAgenelist) <- gsub("_.*", "", names(GSEAgenelist))  # first symbol only (avoid ties)
  GSEAgenelist <- sort(GSEAgenelist, decreasing = TRUE)

  # Sorted Beta estimates (used for cnetplot foldChange in viz).
  Betas <- Eassoc_df[[estimate_col]]
  names(Betas) <- Eassoc_df$omicID
  names(Betas) <- gsub("_.*", "", names(Betas))
  Betas <- sort(Betas, decreasing = TRUE)

  GeneLists[["ORA_all"]]  <- ORA_all
  GeneLists[["ORA_up"]]   <- ORA_up
  GeneLists[["ORA_down"]] <- ORA_down
  GeneLists[["GSEA"]]     <- GSEAgenelist
  GeneLists[["Betas"]]    <- Betas

  GeneLists
}

# ----------------------------------------------------------------------------
# Entrez mapping
#
# Replaces the legacy live-BioMart entrezMap()/df_entrez with the cached
# OlinkEntrezConv.txt produced by 01_entrez_map.R. That file is the OmicsPred
# portal CSV left-joined to BioMart hgnc_symbol -> entrezgene_id, so it carries
# at least: Gene (Olink protein/complex id), genes (split symbol), entrezgene_id.
# ----------------------------------------------------------------------------
load_entrez_map <- function(path = .m4_entrez()) {
  .require_input(path, "cached BioMart Entrez map (OlinkEntrezConv.txt)")
  df <- as.data.frame(fread(path), stringsAsFactors = FALSE, check.names = FALSE)
  needed <- c("Gene", "entrezgene_id")
  if (!all(needed %in% names(df)))
    stop("[module4_enrichment] OlinkEntrezConv.txt must contain columns Gene and entrezgene_id; found: ",
         paste(names(df), collapse = ", "), call. = FALSE)
  df
}

# Attach Entrez IDs to each element of a gene list. For ORA character vectors
# the result is an (unnamed) Entrez vector; for named GSEA vectors the *names*
# are remapped to Entrez while values (ranks) are preserved. Mirrors legacy
# convertEntrez(): adds a "<name>_entrez" entry for every input entry.
convertEntrez <- function(GeneList, df_entrez = load_entrez_map()) {
  for (i in names(GeneList)) {
    name <- paste0(i, "_entrez")
    if (is.null(names(GeneList[[i]]))) {          # ORA: unnamed symbol vector
      GeneList[[name]] <- df_entrez$entrezgene_id[match(GeneList[[i]], df_entrez$Gene)]
    } else {                                       # GSEA: named ranked vector
      GeneList[[name]] <- GeneList[[i]]
      names(GeneList[[name]]) <- df_entrez$entrezgene_id[match(names(GeneList[[i]]), df_entrez$Gene)]
    }
  }
  GeneList
}

# Convert Entrez ID(s) back to gene symbol(s) using the cached map.
convertEntreztoSymbol <- function(entrez_ids, df_entrez = load_entrez_map()) {
  df_entrez$Gene[match(as.character(entrez_ids), as.character(df_entrez$entrezgene_id))]
}

# OmicsPred protein universe (symbols, complexes split) and Entrez equivalents.
ukb_universe <- function(df_entrez = NULL) {
  f <- .omicspred_universe_file()
  .require_input(f, "OmicsPred protein universe (portal csv)")
  omicpredIDs <- fread(f)
  syms <- unlist(strsplit(omicpredIDs$Gene, "_"))
  out <- list(symbol = syms)
  if (!is.null(df_entrez)) {
    ez <- df_entrez$entrezgene_id[match(syms, df_entrez$Gene)]
    out$entrez <- as.character(ez)
  }
  out
}

# ----------------------------------------------------------------------------
# GSEA (primary, manuscript path)
# ----------------------------------------------------------------------------

# Tissue GSEA over GTEx Term2Gene. Legacy params: minGSSize=10, maxGSSize=500,
# pvalueCutoff=0.05, pAdjustMethod="BH", TERM2NAME=NA. Returns @result data.frame.
PSEA_tissue <- function(GeneList,
                        gtex_file = .m4_geneset("GTEX_tissue.txt")) {
  if (!requireNamespace("clusterProfiler", quietly = TRUE))
    stop("[module4_enrichment] package 'clusterProfiler' is required for GSEA.", call. = FALSE)
  .require_input(gtex_file, "GTEx tissue Term2Gene (genesets/GTEX_tissue.txt)")
  GTEX <- fread(gtex_file)

  gsea <- clusterProfiler::GSEA(
    geneList     = GeneList[["GSEA"]],
    minGSSize    = 10,
    maxGSSize    = 500,
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    TERM2GENE    = GTEX,
    TERM2NAME    = NA)
  gsea@result
}

# Reactome pathway GSEA over Entrez-ranked list. Legacy params: minGSSize=10,
# maxGSSize=500, pvalueCutoff=0.05, pAdjustMethod="BH", eps=0. @result data.frame.
PSEA_paths <- function(GeneList) {
  if (!requireNamespace("ReactomePA", quietly = TRUE))
    stop("[module4_enrichment] package 'ReactomePA' is required for pathway GSEA.", call. = FALSE)
  ranked <- GeneList[["GSEA_entrez"]]
  ranked <- ranked[!is.na(names(ranked))]          # drop genes without Entrez map
  Reactome <- ReactomePA::gsePathway(
    ranked,
    minGSSize    = 10,
    maxGSSize    = 500,
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    verbose      = TRUE,
    eps          = 0)
  Reactome@result
}

# ----------------------------------------------------------------------------
# ORA (separate, optional functionality)
# ----------------------------------------------------------------------------

# Tissue over-representation over GTEx Term2Gene. `Type` selects which gene set
# (ORA_all/ORA_up/ORA_down). Legacy params: pvalueCutoff=0.05, BH, minGSSize=10,
# maxGSSize=500, qvalueCutoff=0.2, universe = OmicsPred symbols, TERM2NAME=NA.
ORA_tissue <- function(GeneList, Type = "ORA_all",
                       gtex_file = .m4_geneset("GTEX_tissue.txt")) {
  if (!requireNamespace("clusterProfiler", quietly = TRUE))
    stop("[module4_enrichment] package 'clusterProfiler' is required for ORA.", call. = FALSE)
  .require_input(gtex_file, "GTEx tissue Term2Gene (genesets/GTEX_tissue.txt)")
  GTEX <- fread(gtex_file)
  universe_syms <- ukb_universe()$symbol

  ora <- clusterProfiler::enricher(
    gene         = GeneList[[Type]],
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    universe     = universe_syms,
    minGSSize    = 10,
    maxGSSize    = 500,
    qvalueCutoff = 0.2,
    TERM2GENE    = GTEX,
    TERM2NAME    = NA)
  if (is.null(ora)) return(NULL)
  ora@result
}

# Pathway over-representation (KEGG, GO BP/CC/MF, Reactome, DO) over Entrez sets.
# Legacy params preserved: pvalueCutoff=0.05, BH, qvalueCutoff=0.2,
# minGSSize=10, maxGSSize=500, universe = OmicsPred Entrez. Each enrichment is
# wrapped in tryCatch -> NULL (legacy ORA_paths v2 behavior). Returns a named
# list of @result data.frames (or NULL for failed/empty enrichments).
ORA_paths <- function(GeneList, Type = "ORA_all_entrez",
                      df_entrez = load_entrez_map()) {
  if (!requireNamespace("clusterProfiler", quietly = TRUE))
    stop("[module4_enrichment] package 'clusterProfiler' is required for ORA.", call. = FALSE)
  EnrichORA <- list()

  # KEGG gene-set tables. KEGG_T2N columns are stored flipped (name,id) in the
  # legacy generator -> reorder to (id,name) as TERM2NAME expects.
  kegg_t2g_f <- .m4_geneset("KEGG_T2G.txt")
  kegg_t2n_f <- .m4_geneset("KEGG_T2N.txt")
  KEGG_T2G <- if (file.exists(kegg_t2g_f)) fread(kegg_t2g_f) else NULL
  KEGG_T2N <- if (file.exists(kegg_t2n_f)) { tn <- fread(kegg_t2n_f); tn[, c(2, 1)] } else NULL

  universe_entrez <- ukb_universe(df_entrez)$entrez
  genes <- GeneList[[Type]]
  genes <- as.character(genes[!is.na(genes)])

  EnrichORA[["KEGG"]] <- if (!is.null(KEGG_T2G)) tryCatch(
    clusterProfiler::enricher(
      gene = genes, pvalueCutoff = 0.05, pAdjustMethod = "BH",
      universe = universe_entrez, qvalueCutoff = 0.2,
      minGSSize = 10, maxGSSize = 500,
      TERM2GENE = KEGG_T2G, TERM2NAME = KEGG_T2N)@result,
    error = function(e) NULL) else NULL

  go_one <- function(ont) {
    if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) return(NULL)
    tryCatch(
      clusterProfiler::enrichGO(
        gene = genes, OrgDb = org.Hs.eg.db::org.Hs.eg.db, keyType = "ENTREZID",
        ont = ont, pvalueCutoff = 0.05, pAdjustMethod = "BH",
        universe = universe_entrez, qvalueCutoff = 0.2,
        minGSSize = 10, maxGSSize = 500, readable = TRUE)@result,
      error = function(e) NULL)
  }
  EnrichORA[["GO_BP"]] <- go_one("BP")
  EnrichORA[["GO_CC"]] <- go_one("CC")
  EnrichORA[["GO_MF"]] <- go_one("MF")

  EnrichORA[["Reactome"]] <- if (requireNamespace("ReactomePA", quietly = TRUE)) tryCatch(
    ReactomePA::enrichPathway(
      gene = genes, pvalueCutoff = 0.05, readable = TRUE,
      universe = universe_entrez, qvalueCutoff = 0.2,
      minGSSize = 10, maxGSSize = 500)@result,
    error = function(e) NULL) else NULL

  EnrichORA[["DO"]] <- if (requireNamespace("DOSE", quietly = TRUE)) tryCatch(
    DOSE::enrichDO(
      gene = genes, ont = "DO", pvalueCutoff = 0.05, pAdjustMethod = "BH",
      universe = universe_entrez, qvalueCutoff = 0.2,
      minGSSize = 10, maxGSSize = 500, readable = TRUE)@result,
    error = function(e) NULL) else NULL

  EnrichORA
}

# Marker so sourcing scripts can assert the library loaded.
.MODULE4_ENRICHMENT_CORE_LOADED <- TRUE
