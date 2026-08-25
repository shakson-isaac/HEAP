#!/usr/bin/env Rscript
# ============================================================================
# 00_build_genesets.R
# Module 4 — Tissue / Pathway enrichment gene-set construction.
#
# Purpose:
#   Build the canonical TERM2GENE gene-set tables used by the enrichment
#   (GSEA/ORA) scripts:
#     * GTEx tissue specificity:
#         - per-gene Tau score + max expression (tissue specificity metric)
#         - per-tissue gene set defined by fold-change > 4 (a gene is assigned
#           to a tissue if its expression there is > 4x the mean of the other
#           tissues), emitted as a two-column Term2Gene table (Tissue, Gene).
#     * HPA tissue specificity (3 nested levels) + secretome:
#         - Specific  : IHC Level == "High"
#         - Enriched  : IHC Level in {High, Medium}
#         - Expressed : IHC Level in {High, Medium, Low, Ascending, Descending}
#         - Secretome : subcellular "Extracellular location" == "Predicted to be secreted"
#     * KEGG pathway tables (gene->pathway T2G, pathway->name T2N) via KEGGREST.
#
# Inputs (IGLOO-canonical via workflow/00_paths.R helpers):
#   heap_gtex_rna_dir()                                  -> 54 gene_tpm_v10_<tissue>.gct files
#   heap_hpa_or_legacy("normal_ihc_data.tsv")            -> HPA IHC tissue levels
#   heap_hpa_or_legacy("subcellular_location.tsv")       -> HPA subcellular / secretome
#   heap_omicspred_or_legacy("UKB_Olink_multi_ancestry_models_val_results_portal.csv")
#                                                        -> Olink protein universe (Gene)
#   KEGGREST (network) for KEGG T2G/T2N (optional; skipped with a warning on failure).
#
# Outputs (under heap_project_output("module4_enrichment", ...)):
#   genesets/GTEX_tissue.txt      (Tissue, Gene)
#   genesets/HPA_specific.txt     (Tissue, Gene name)
#   genesets/HPA_enriched.txt     (Tissue, Gene name)
#   genesets/HPA_expressed.txt    (Tissue, Gene name)
#   genesets/HPA_secretome.txt    (Pathway, Gene name)
#   genesets/KEGG_T2G.txt         (pathway_id, gene_id)   [if KEGG available]
#   genesets/KEGG_T2N.txt         (pathway_id, pathway_name) [if KEGG available]
#   GTEX_tau_scores.csv           (Name, Description, tau_score, max_expr)
#
# Ported from:
#   scripts/visualizations/Visualizations/ModuleExt/GTEX_ident.R  (GTEx Tau + FC>4)
#   scripts/visualizations/Visualizations/ModuleExt/HPA_ident.R   (HPA + KEGG)
#
# Notes / deviations from legacy (see file footer for details):
#   * Legacy hard-coded /n/groups/patel/shakson_ukb/UK_Biobank/... paths are
#     replaced by IGLOO helpers.
#   * The exploratory ggplot Tau-density plot in GTEX_ident.R is dropped here
#     (plotting belongs in 05_enrichment_figures.R); the Tau scores it consumed
#     are written to GTEX_tau_scores.csv instead.
#   * KEGG (network) is isolated in tryCatch so the file-based GTEx/HPA outputs
#     still complete if KEGGREST/network is unavailable.
# ============================================================================

local({
  cand <- c(Sys.getenv("HEAP_PATHS_FILE", unset = ""),
            file.path(getwd(), "workflow", "00_paths.R"),
            "/n/groups/patel/shakson_ukb/HEAP/workflow/00_paths.R")
  hit <- cand[nzchar(cand) & file.exists(cand)][1]
  if (!is.na(hit)) source(hit)
})

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
})

OUT <- heap_project_output("module4_enrichment")
GENESETS <- file.path(OUT, "genesets")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(GENESETS, recursive = TRUE, showWarnings = FALSE)

# ============================================================================
# 1) GTEx tissue specificity (ported from GTEX_ident.R)
# ============================================================================
message("[GTEx] reading median-TPM GCTs from ", heap_gtex_rna_dir())

file_list <- list.files(path = heap_gtex_rna_dir(), full.names = TRUE,
                        pattern = "\\.gct$")
if (length(file_list) == 0L)
  stop("No GTEx .gct files found in ", heap_gtex_rna_dir())

# fread auto-skips the 2-line GCT header (#1.2 + dims) and reads the table
# starting at the Name/Description column row.
data_list <- lapply(file_list, fread)

tissue_ids <- basename(file_list)
tissue_ids <- gsub("gene_tpm_v10_", "", tissue_ids)
tissue_ids <- gsub("\\.gct", "", tissue_ids)
names(data_list) <- tissue_ids

# General transcriptome: collapse each tissue's per-sample columns to the
# per-gene median expression, renamed to the tissue id.
median_expr <- function(df, tissueID) {
  df2 <- df %>%
    pivot_longer(-c(Name, Description)) %>%
    group_by(Name) %>%
    summarize(median_value = median(value), .groups = "drop") %>%
    left_join(df, by = "Name") %>%
    select(Name, Description, median_value)

  colnames(df2) <- c("Name", "Description", tissueID)
  df2
}

data_list <- lapply(names(data_list), function(id) {
  median_expr(data_list[[id]], id)
})

# Main GTEx data frame: one column per tissue (median TPM).
GTEX_df <- data_list %>% reduce(full_join, by = c("Name", "Description"))

# Tau tissue-specificity score (and max expression) per gene.
# Tau = sum(1 - x_i/max(x)) / (n - 1), computed on log1p(TPM) so all values
# are non-negative (legacy formula, unchanged).
GTEX_spec <- GTEX_df %>%
  pivot_longer(-c(Name, Description)) %>%
  group_by(Name) %>%
  mutate(value = log1p(value)) %>%
  summarize(tau_score = sum(1 - (value / max(value))) / (length(value) - 1),
            max_expr = max(value),
            .groups = "drop") %>%
  left_join(GTEX_df, by = "Name") %>%
  select(Name, Description, tau_score, max_expr)

# Fold change per tissue: x_i / mean(x_-i) = x_i*(n-1) / (sum(x) - x_i).
expression_data <- GTEX_df %>%
  select(-Name, -Description)

n_tissue <- ncol(expression_data)

fold_change_results <- apply(expression_data, 1, function(x) {
  total_sum <- sum(x)
  x / ((total_sum - x) / (length(x) - 1))
})

fold_change_df <- as.data.frame(t(fold_change_results))
fold_change_df$Name <- GTEX_df$Name
fold_change_df$Description <- GTEX_df$Description

# Per-tissue gene set: genes with fold change > 4 in that tissue.
tissue_cols <- colnames(fold_change_df)[seq_len(n_tissue)]
GTEx_GeneSet <- lapply(tissue_cols, function(x) {
  fold_change_df %>%
    filter(.[[x]] > 4) %>%
    pull(Description)
})
names(GTEx_GeneSet) <- tissue_cols

# Build the two-column Term2Gene (Tissue, Gene) table.
Term2Gene <- data.frame(Term = character(0), Gene = character(0))
for (term in names(GTEx_GeneSet)) {
  genes <- GTEx_GeneSet[[term]]
  if (length(genes) == 0L) next
  Term2Gene <- rbind(Term2Gene,
                     data.frame(Term = rep(term, length(genes)), Gene = genes))
}

write.table(Term2Gene, file = file.path(GENESETS, "GTEX_tissue.txt"),
            sep = "\t", row.names = FALSE)
fwrite(GTEX_spec, file = file.path(OUT, "GTEX_tau_scores.csv"))
message("[GTEx] wrote GTEX_tissue.txt (", nrow(Term2Gene), " rows) and GTEX_tau_scores.csv")

# ============================================================================
# 2) HPA tissue specificity + secretome (ported from HPA_ident.R)
# ============================================================================
hpa_ihc  <- heap_hpa_or_legacy("normal_ihc_data.tsv")
hpa_subc <- heap_hpa_or_legacy("subcellular_location.tsv")
message("[HPA] reading ", hpa_ihc, " and ", hpa_subc)

tissueinfo  <- fread(file = hpa_ihc)
secreteinfo <- fread(file = hpa_subc)

# Nested IHC-level heuristics.
SpecificSet <- tissueinfo %>%
  filter(Level %in% c("High")) %>%
  select(c("Gene", "Gene name", "Tissue"))
EnrichedSet <- tissueinfo %>%
  filter(Level %in% c("High", "Medium")) %>%
  select(c("Gene", "Gene name", "Tissue"))
ExpressedSet <- tissueinfo %>%
  filter(Level %in% c("High", "Medium", "Low", "Ascending", "Descending")) %>%
  select(c("Gene", "Gene name", "Tissue"))
Secretome <- secreteinfo %>%
  filter(`Extracellular location` %in% c("Predicted to be secreted")) %>%
  mutate(Pathway = "Secretome") %>%
  select(c("Gene", "Gene name", "Pathway"))

write.table(SpecificSet[, c("Tissue", "Gene name")],
            file = file.path(GENESETS, "HPA_specific.txt"),
            sep = "\t", row.names = FALSE)
write.table(EnrichedSet[, c("Tissue", "Gene name")],
            file = file.path(GENESETS, "HPA_enriched.txt"),
            sep = "\t", row.names = FALSE)
write.table(ExpressedSet[, c("Tissue", "Gene name")],
            file = file.path(GENESETS, "HPA_expressed.txt"),
            sep = "\t", row.names = FALSE)
write.table(Secretome[, c("Pathway", "Gene name")],
            file = file.path(GENESETS, "HPA_secretome.txt"),
            sep = "\t", row.names = FALSE)
message("[HPA] wrote HPA_specific/enriched/expressed/secretome.txt")

# ============================================================================
# 3) KEGG pathway tables (ported from HPA_ident.R) — network, isolated.
# ============================================================================
kegg_ok <- tryCatch({
  suppressPackageStartupMessages(library(KEGGREST))

  gene2pathway <- keggLink("pathway", "hsa")
  df3 <- data.frame(
    gene_id    = gsub("hsa:", "", names(gene2pathway)),
    pathway_id = gsub("path:", "", gene2pathway)
  )

  pathway2name <- keggList("pathway", "hsa")
  df4 <- data.frame(
    gene_id    = names(pathway2name),
    pathway_id = gsub(" - Homo sapiens \\(human\\)", "", pathway2name)
  )

  write.table(df3[, c("pathway_id", "gene_id")],
              file = file.path(GENESETS, "KEGG_T2G.txt"),
              sep = "\t", row.names = FALSE)
  write.table(df4[, c("pathway_id", "gene_id")],
              file = file.path(GENESETS, "KEGG_T2N.txt"),
              sep = "\t", row.names = FALSE)
  message("[KEGG] wrote KEGG_T2G.txt and KEGG_T2N.txt")
  TRUE
}, error = function(e) {
  warning("[KEGG] skipped (KEGGREST/network unavailable): ",
          conditionMessage(e), call. = FALSE)
  FALSE
})

message("Done. Outputs under ", OUT,
        if (!kegg_ok) "  (KEGG tables skipped)" else "")
