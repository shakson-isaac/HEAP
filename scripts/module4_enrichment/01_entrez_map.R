#!/usr/bin/env Rscript
# ============================================================================
# 01_entrez_map.R
# Module 4 — Olink protein symbol -> Entrez ID map (BioMart, cached).
#
# Purpose:
#   Build the HGNC gene-symbol -> Entrez gene-id mapping for the Olink protein
#   universe (OmicsPred) so the enrichment scripts can convert protein lists to
#   Entrez IDs for pathway analysis. The OmicsPred "Gene" field can encode
#   protein complexes as underscore-joined symbols (e.g. "A_B"); these are split
#   so each constituent symbol is mapped individually, then re-joined back onto
#   the original protein rows.
#
# Cache-first behavior:
#   If OlinkEntrezConv.txt already exists it is loaded and returned WITHOUT a
#   BioMart query (BioMart needs network and is slow/flaky). Set the env var
#   HEAP_FORCE_BIOMART=1 to force a refresh.
#
# Inputs (IGLOO-canonical via workflow/00_paths.R helpers):
#   heap_omicspred_or_legacy("UKB_Olink_multi_ancestry_models_val_results_portal.csv")
#       -> Olink protein universe (column "Gene")
#   biomaRt -> Ensembl hsapiens_gene_ensembl (network)
#
# Outputs (under heap_project_output("module4_enrichment", ...)):
#   OlinkEntrezConv.txt   (original OmicsPred columns + genes, entrezgene_id)
#
# Ported from:
#   scripts/visualizations/Visualizations/ModuleExt/ObtainEntrezMapping.R
#
# Deviations from legacy:
#   * Legacy hard-coded UK_Biobank paths replaced by IGLOO helpers.
#   * Added cache-first load so re-runs reuse OlinkEntrezConv.txt instead of
#     re-querying BioMart (per Module 4 I/O contract).
# ============================================================================

local({
  cand <- c(Sys.getenv("HEAP_PATHS_FILE", unset = ""),
            file.path(getwd(), "workflow", "00_paths.R"),
            "/n/groups/patel/shakson_ukb/HEAP/workflow/00_paths.R")
  hit <- cand[nzchar(cand) & file.exists(cand)][1]
  if (!is.na(hit)) source(hit)
})

suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
})

OUT <- heap_project_output("module4_enrichment")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

CACHE <- file.path(OUT, "OlinkEntrezConv.txt")
force_refresh <- nzchar(Sys.getenv("HEAP_FORCE_BIOMART", unset = ""))

# Function: Obtain Olink protein symbols -> Entrez IDs map via BioMart.
entrezMap <- function() {
  suppressPackageStartupMessages({
    library(org.Hs.eg.db)
    library(biomaRt)
  })

  omics_file <- heap_omicspred_or_legacy(
    "UKB_Olink_multi_ancestry_models_val_results_portal.csv")
  omicpredIDs <- fread(file = omics_file)

  # Connect to the Ensembl database.
  ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

  map_prot_to_entrez <- function(hgnc_ids) {
    getBM(attributes = c("hgnc_symbol", "entrezgene_id"),
          filters = "hgnc_symbol",
          values = hgnc_ids,
          mart = ensembl)
  }

  # Split protein-complex symbols ("A_B") into separate genes.
  df_split <- omicpredIDs %>%
    mutate(genes = str_split(Gene, "_")) %>%
    unnest(genes)

  entrez_results <- map_prot_to_entrez(df_split$genes)

  df_entrez <- df_split %>%
    left_join(entrez_results, by = c("genes" = "hgnc_symbol"))

  # Report unmapped symbols (informational; matches legacy intent).
  OlinkIDs_noconvert <- unique(df_entrez$Gene[is.na(df_entrez$entrezgene_id)])
  if (length(OlinkIDs_noconvert) > 0L)
    message("[biomaRt] ", length(OlinkIDs_noconvert),
            " Olink Gene entries had >=1 unmapped symbol.")

  df_entrez
}

if (file.exists(CACHE) && !force_refresh) {
  message("[cache] OlinkEntrezConv.txt exists; loading instead of querying BioMart: ",
          CACHE)
  UKBentrez <- fread(file = CACHE)
} else {
  message("[biomaRt] querying Ensembl for HGNC -> Entrez map ...")
  UKBentrez <- entrezMap()
  fwrite(UKBentrez, file = CACHE)
  message("[biomaRt] wrote ", CACHE, " (", nrow(UKBentrez), " rows)")
}

# Utility across other scripts:
#   df_entrez <- fread(file = file.path(
#     heap_project_output("module4_enrichment"), "OlinkEntrezConv.txt"))
invisible(UKBentrez)
