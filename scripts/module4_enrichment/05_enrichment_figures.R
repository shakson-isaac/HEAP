#!/usr/bin/env Rscript
# ============================================================================
# 05_enrichment_figures.R  [scaffold — not yet implemented]
# Module 4 enrichment. Port from: Module2/HEAPassoc_pathwayviz.R
# Purpose: PLOTTING ONLY: ComplexHeatmap tissue/pathway + GTEx Tau density
# Inputs : tissue_enrichment.csv; pathway_enrichment.csv (+ figures dir)
# Outputs: figures/main|supplement (via export_helpers)  (under heap_project_output("module4_enrichment", ...))
# See docs/SUPPORT_ANALYSIS_PLAN.md §2.
# ============================================================================
local({
  cand <- c(Sys.getenv("HEAP_PATHS_FILE", unset = ""),
            file.path(getwd(), "workflow", "00_paths.R"),
            "/n/groups/patel/shakson_ukb/HEAP/workflow/00_paths.R")
  hit <- cand[nzchar(cand) & file.exists(cand)][1]
  if (!is.na(hit)) source(hit)
})
OUT <- heap_project_output("module4_enrichment")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

stop("module4_enrichment/05_enrichment_figures.R is a scaffold stub. Implement by porting from: Module2/HEAPassoc_pathwayviz.R")
