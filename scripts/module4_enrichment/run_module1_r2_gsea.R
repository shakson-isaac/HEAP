#!/usr/bin/env Rscript

# ============================================================================
# run_module1_r2_gsea.R
# ----------------------------------------------------------------------------
# Gene-Set Enrichment Analysis (GSEA) of the Module 1 variance decomposition:
# rank the measured proteome by each component's UNIQUE out-of-fold predictive
# R2 (Genetic R2_G, Exposomic R2_E, GxE R2_GxE; score_unique_drop, mean over
# folds) and ask which GO terms / Reactome pathways are enriched among the
# proteins each component explains best. This is manuscript Figure 2E.
#
# Engine: clusterProfiler::gseGO (keyType = SYMBOL, OrgDb = org.Hs.eg.db,
#   ont = ALL -> BP/MF/CC) and ReactomePA::gsePathway (Entrez-ranked). The
#   universe is the ranked list itself = the 2,686 measured proteins (correct
#   proteomic background). Params mirror the HEAP enrichment scaffold:
#   minGSSize=10, maxGSSize=500, pvalueCutoff=0.05, pAdjustMethod="BH", eps=0,
#   seed=TRUE (reproducible permutations).
#
# Output: docs/manuscript_stats/module1_enrichment/
#   SUMMARY.md                          headline + top terms per component
#   gsea_GO_<component>.tsv             all significant GO terms (p.adjust<0.05)
#   gsea_Reactome_<component>.tsv       all significant Reactome pathways
#   gsea_all_significant.tsv            long table across components/databases
#   + the full gseaResult objects in module4_enrichment/module1_r2_gsea.qs
#
# Run (a few minutes):
#   HEAP_PATHS_FILE=.../workflow/00_paths.R \
#     Rscript scripts/module4_enrichment/run_module1_r2_gsea.R [covarType] [method]
#   Defaults: base / lasso / M1_base_lasso.
# ============================================================================

local({
  cand <- c(Sys.getenv("HEAP_PATHS_FILE", unset = ""),
            file.path(getwd(), "workflow", "00_paths.R"),
            "/n/groups/patel/shakson_ukb/HEAP/workflow/00_paths.R")
  hit <- cand[nzchar(cand) & file.exists(cand)][1]
  if (!is.na(hit)) source(hit)
})
if (!exists("heap_project_output"))
  stop("could not load workflow/00_paths.R (set HEAP_PATHS_FILE).")
source(file.path(heap_path(), "scripts", "visualizations", "common", "figure_paths.R"))
source(file.path(heap_path(), "scripts", "visualizations", "common", "load_heap_results.R"))

suppressPackageStartupMessages({
  library(data.table)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ReactomePA)
})

a <- commandArgs(trailingOnly = TRUE); a <- a[!startsWith(a, "--")]
covarType  <- if (length(a) >= 1) a[1] else "base"
method     <- if (length(a) >= 2) a[2] else "lasso"
experiment <- if (length(a) >= 3) a[3] else "M1_base_lasso"

# ---- per-protein unique R2 by component -------------------------------------
co <- load_module1_predictive_r2(covarType, method, level = "coarse", experiment = experiment)
ud <- co[get("method") == "score_unique_drop" & block %in% c("G", "E", "GxE"),
         .(r2 = mean(r2)), by = .(omic, block)]
W  <- dcast(ud, omic ~ block, value.var = "r2")

components <- list(Genetic = "G", Exposomic = "E", GxE = "GxE")

# build a decreasing named (SYMBOL -> R2) vector for one component
make_ranklist <- function(col) {
  d <- W[is.finite(get(col)), .(omic, r2 = get(col))]
  d <- d[, .(r2 = max(r2)), by = omic]          # 1 value per symbol
  v <- setNames(d$r2, d$omic)
  sort(v, decreasing = TRUE)
}

run_go <- function(gl) {
  tryCatch(
    clusterProfiler::gseGO(geneList = gl, OrgDb = org.Hs.eg.db, keyType = "SYMBOL",
                           ont = "ALL", minGSSize = 10, maxGSSize = 500,
                           pvalueCutoff = 0.05, pAdjustMethod = "BH",
                           eps = 0, seed = TRUE, verbose = FALSE),
    error = function(e) { message("  gseGO failed: ", conditionMessage(e)); NULL })
}
run_reactome <- function(gl) {
  sym <- names(gl)
  map <- suppressWarnings(clusterProfiler::bitr(sym, "SYMBOL", "ENTREZID", org.Hs.eg.db))
  map <- map[!duplicated(map$SYMBOL), ]
  gl2 <- gl[map$SYMBOL]; names(gl2) <- map$ENTREZID
  gl2 <- sort(gl2[!is.na(names(gl2))], decreasing = TRUE)
  tryCatch(
    ReactomePA::gsePathway(geneList = gl2, organism = "human", minGSSize = 10,
                           maxGSSize = 500, pvalueCutoff = 0.05,
                           pAdjustMethod = "BH", eps = 0, seed = TRUE, verbose = FALSE),
    error = function(e) { message("  gsePathway failed: ", conditionMessage(e)); NULL })
}

res <- list(); all_sig <- list()
for (comp in names(components)) {
  message("== ", comp, " ==")
  gl <- make_ranklist(components[[comp]])
  message("  ranked ", length(gl), " proteins (top: ", names(gl)[1], " R2=", round(gl[1], 3), ")")
  go  <- run_go(gl)
  re  <- run_reactome(gl)
  res[[comp]] <- list(GO = go, Reactome = re)
  if (!is.null(go) && nrow(as.data.frame(go)))
    all_sig[[paste0(comp, "_GO")]] <-
      data.table(component = comp, database = "GO", as.data.table(as.data.frame(go)))
  if (!is.null(re) && nrow(as.data.frame(re)))
    all_sig[[paste0(comp, "_Reactome")]] <-
      data.table(component = comp, database = "Reactome", ONTOLOGY = NA_character_,
                 as.data.table(as.data.frame(re)))
  message("  GO sig: ", if (is.null(go)) 0 else nrow(as.data.frame(go)),
          " | Reactome sig: ", if (is.null(re)) 0 else nrow(as.data.frame(re)))
}

# ============================================================================
# Write outputs
# ============================================================================
outdir <- file.path(heap_path(), "docs", "manuscript_stats", "module1_enrichment")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

slim <- function(dt) {
  keep <- intersect(c("component", "database", "ONTOLOGY", "ID", "Description",
                      "setSize", "NES", "enrichmentScore", "pvalue", "p.adjust",
                      "qvalue", "core_enrichment"), names(dt))
  dt[, ..keep]
}
combined <- if (length(all_sig)) rbindlist(all_sig, fill = TRUE) else data.table()
if (nrow(combined)) {
  combined <- slim(combined)
  setorder(combined, component, p.adjust)
  fwrite(combined, file.path(outdir, "gsea_all_significant.tsv"), sep = "\t")
  for (comp in names(components)) {
    g <- combined[component == comp & database == "GO"]
    r <- combined[component == comp & database == "Reactome"]
    if (nrow(g)) fwrite(g, file.path(outdir, paste0("gsea_GO_", comp, ".tsv")), sep = "\t")
    if (nrow(r)) fwrite(r, file.path(outdir, paste0("gsea_Reactome_", comp, ".tsv")), sep = "\t")
  }
}
# full objects for any downstream figure
if (requireNamespace("qs", quietly = TRUE))
  qs::qsave(res, file.path(heap_project_output("module4_enrichment"), "module1_r2_gsea.qs"))

# ---- SUMMARY.md -------------------------------------------------------------
fmt <- function(x, d = 2) formatC(x, format = "f", digits = d)
sci <- function(x) formatC(x, format = "e", digits = 1)
mk <- character(0); add <- function(...) mk <<- c(mk, paste0(...))
add("# Module 1 — GSEA of variance-decomposition R2 rankings")
add("")
add("_Generated by `scripts/module4_enrichment/run_module1_r2_gsea.R` from `",
    experiment, "/", covarType, "/", method, "`. Re-run to refresh._")
add("")
add("GSEA (clusterProfiler::gseGO ont=ALL + ReactomePA::gsePathway) on the measured ",
    "proteome ranked by each component's **unique out-of-fold R2**. A positive NES ",
    "means the gene set is concentrated among the proteins that component explains ",
    "best. Significance: BH p.adjust < 0.05. Universe = the ", nrow(W), " ranked proteins.")
add("")
if (!nrow(combined)) {
  add("**No gene sets reached p.adjust < 0.05 for any component.**")
} else {
  add("## Counts of significant gene sets (p.adjust < 0.05)")
  add("")
  add("| component | GO | Reactome |")
  add("|---|--:|--:|")
  for (comp in names(components))
    add("| ", comp, " | ", nrow(combined[component == comp & database == "GO"]),
        " | ", nrow(combined[component == comp & database == "Reactome"]), " |")
  add("")
  for (comp in names(components)) {
    cc <- combined[component == comp][order(p.adjust)]
    if (!nrow(cc)) { add("### ", comp, ": no significant terms"); add(""); next }
    add("### ", comp, " — top terms (by p.adjust)")
    add("")
    add("| database | term | setSize | NES | p.adjust |")
    add("|---|---|--:|--:|--:|")
    for (i in seq_len(min(12, nrow(cc))))
      add("| ", cc$database[i], if (!is.na(cc$ONTOLOGY[i])) paste0(" (", cc$ONTOLOGY[i], ")") else "",
          " | ", cc$Description[i], " | ", cc$setSize[i], " | ", fmt(cc$NES[i]),
          " | ", sci(cc$p.adjust[i]), " |")
    add("")
  }
}
add("_Full tables: gsea_all_significant.tsv + per-component gsea_GO_*/gsea_Reactome_*.tsv. ",
    "gseaResult objects: module4_enrichment/module1_r2_gsea.qs._")
writeLines(mk, file.path(outdir, "SUMMARY.md"))

message("\nWrote enrichment summary to: ", outdir)
cat(sprintf("Significant gene sets: %s\n",
            if (nrow(combined)) paste(sapply(names(components), function(c)
              sprintf("%s=%d", c, nrow(combined[component == c]))), collapse = " | ") else "none"))
