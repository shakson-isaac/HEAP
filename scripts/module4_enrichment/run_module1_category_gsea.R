#!/usr/bin/env Rscript

# ============================================================================
# run_module1_category_gsea.R
# ----------------------------------------------------------------------------
# PER-CATEGORY GSEA of the Module 1 exposomic decomposition. Instead of ranking
# proteins by the AGGREGATE exposomic R2 (which is dominated by diet/exercise and
# thus confounded), rank the measured proteome by EACH exposure category's unique
# out-of-fold R2 (score_unique_drop, level = exposure_categories, mean over folds)
# and ask which GO terms / Reactome pathways each category's proteins enrich.
#
# Engine mirrors run_module1_r2_gsea.R: clusterProfiler::gseGO (SYMBOL, ont=ALL)
# + ReactomePA::gsePathway (Entrez). minGSSize=10, maxGSSize=500, BH<0.05, eps=0,
# seed=TRUE. Universe = the ranked proteome.
#
# Only categories with REACH >= MIN_REACH proteins (R2 > 0.005) are run (sparse
# categories cannot give a meaningful enrichment). MIN_REACH via env (default 10).
#
# Output: docs/manuscript_stats/module1_enrichment_bycategory/
#   gsea_all_significant_bycategory.tsv   long table category x database x term
#   SUMMARY.md                            top terms per category + reach counts
#   module4_enrichment/module1_category_gsea.qs   full gseaResult objects
#
# Run (~10-20 min):
#   HEAP_PATHS_FILE=.../workflow/00_paths.R \
#     Rscript scripts/module4_enrichment/run_module1_category_gsea.R [covarType] [method]
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
source(file.path(heap_path(), "scripts", "visualizations", "common", "label_helpers.R"))

suppressPackageStartupMessages({
  library(data.table)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ReactomePA)
})

a <- commandArgs(trailingOnly = TRUE); a <- a[!startsWith(a, "--")]
covarType  <- if (length(a) >= 1) a[1] else "base"
method     <- if (length(a) >= 2) a[2] else "lasso"
MIN_REACH  <- as.integer(Sys.getenv("HEAP_MIN_REACH", "10"))
THR        <- 0.005

# ---- per-protein unique R2 by exposure CATEGORY -----------------------------
ec <- load_module1_predictive_r2(covarType, method, level = "exposure_categories")
if ("method" %in% names(ec)) ec <- ec[get("method") == "score_unique_drop"]
pp <- ec[, .(r2 = mean(r2)), by = .(omic, category)]
reach <- pp[, .(reach = sum(r2 > THR)), by = category][order(-reach)]
cats  <- reach[reach >= MIN_REACH, as.character(category)]
message("categories by reach (R2>", THR, "):")
print(reach)
message("\nrunning GSEA on ", length(cats), " categories with reach >= ", MIN_REACH,
        ": ", paste(cats, collapse = ", "))

make_ranklist <- function(cat) {
  d <- pp[as.character(category) == cat & is.finite(r2), .(omic, r2)]
  d <- d[, .(r2 = max(r2)), by = omic]            # 1 value per symbol
  v <- setNames(d$r2, d$omic)
  sort(v, decreasing = TRUE)
}
run_go <- function(gl) tryCatch(
  clusterProfiler::gseGO(geneList = gl, OrgDb = org.Hs.eg.db, keyType = "SYMBOL",
                         ont = "ALL", minGSSize = 10, maxGSSize = 500,
                         pvalueCutoff = 0.05, pAdjustMethod = "BH", eps = 0,
                         seed = TRUE, verbose = FALSE),
  error = function(e) { message("  gseGO failed: ", conditionMessage(e)); NULL })
run_reactome <- function(gl) {
  sym <- names(gl)
  map <- suppressWarnings(clusterProfiler::bitr(sym, "SYMBOL", "ENTREZID", org.Hs.eg.db))
  map <- map[!duplicated(map$SYMBOL), ]
  gl2 <- gl[map$SYMBOL]; names(gl2) <- map$ENTREZID
  gl2 <- sort(gl2[!is.na(names(gl2))], decreasing = TRUE)
  tryCatch(
    ReactomePA::gsePathway(geneList = gl2, organism = "human", minGSSize = 10,
                           maxGSSize = 500, pvalueCutoff = 0.05, pAdjustMethod = "BH",
                           eps = 0, seed = TRUE, verbose = FALSE),
    error = function(e) { message("  gsePathway failed: ", conditionMessage(e)); NULL })
}

res <- list(); all_sig <- list()
for (cat in cats) {
  message("== ", cat, " ==")
  gl <- make_ranklist(cat)
  message("  ranked ", length(gl), " proteins (top: ", names(gl)[1], " R2=", round(gl[1], 3), ")")
  go <- run_go(gl); re <- run_reactome(gl)
  res[[cat]] <- list(GO = go, Reactome = re)
  if (!is.null(go) && nrow(as.data.frame(go)))
    all_sig[[paste0(cat, "_GO")]] <-
      data.table(category = cat, database = "GO", as.data.table(as.data.frame(go)))
  if (!is.null(re) && nrow(as.data.frame(re)))
    all_sig[[paste0(cat, "_Reactome")]] <-
      data.table(category = cat, database = "Reactome", ONTOLOGY = NA_character_,
                 as.data.table(as.data.frame(re)))
  message("  GO sig: ", if (is.null(go)) 0 else nrow(as.data.frame(go)),
          " | Reactome sig: ", if (is.null(re)) 0 else nrow(as.data.frame(re)))
}

outdir <- file.path(heap_path(), "docs", "manuscript_stats", "module1_enrichment_bycategory")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
slim <- function(dt) {
  keep <- intersect(c("category","database","ONTOLOGY","ID","Description","setSize",
                      "NES","enrichmentScore","pvalue","p.adjust","qvalue","core_enrichment"), names(dt))
  dt[, ..keep]
}
combined <- if (length(all_sig)) rbindlist(all_sig, fill = TRUE) else data.table()
if (nrow(combined)) {
  combined <- slim(combined); setorder(combined, category, p.adjust)
  fwrite(combined, file.path(outdir, "gsea_all_significant_bycategory.tsv"), sep = "\t")
}
if (requireNamespace("qs", quietly = TRUE))
  qs::qsave(res, file.path(heap_project_output("module4_enrichment"), "module1_category_gsea.qs"))

# ---- SUMMARY.md -------------------------------------------------------------
fmt <- function(x, d = 2) formatC(x, format = "f", digits = d)
sci <- function(x) formatC(x, format = "e", digits = 1)
mk <- character(0); add <- function(...) mk <<- c(mk, paste0(...))
add("# Module 1 — PER-CATEGORY GSEA of exposomic R2 rankings"); add("")
add("_`scripts/module4_enrichment/run_module1_category_gsea.R` (", covarType, "/", method,
    "). Proteome ranked by EACH category's unique out-of-fold R2; GSEA (gseGO ont=ALL + ",
    "gsePathway). BH<0.05. Categories with reach>=", MIN_REACH, " (R2>", THR, ") only._"); add("")
add("## Reach per category (proteins R2>", THR, ")"); add("")
add("| category | reach | run? | GO sig | Reactome sig |"); add("|---|--:|:--:|--:|--:|")
for (i in seq_len(nrow(reach))) {
  ct <- as.character(reach$category[i]); rn <- ct %in% cats
  ng <- if (nrow(combined)) nrow(combined[category == ct & database == "GO"]) else 0
  nr <- if (nrow(combined)) nrow(combined[category == ct & database == "Reactome"]) else 0
  add("| ", ct, " | ", reach$reach[i], " | ", if (rn) "yes" else "no", " | ",
      if (rn) ng else "-", " | ", if (rn) nr else "-", " |")
}
add("")
for (cat in cats) {
  cc <- if (nrow(combined)) combined[category == cat][order(p.adjust)] else data.table()
  add("### ", cat, if (nrow(cc)) "" else " — no significant terms"); add("")
  if (nrow(cc)) {
    add("| database | term | setSize | NES | p.adjust |"); add("|---|---|--:|--:|--:|")
    for (i in seq_len(min(12, nrow(cc))))
      add("| ", cc$database[i], if (!is.na(cc$ONTOLOGY[i])) paste0(" (", cc$ONTOLOGY[i], ")") else "",
          " | ", cc$Description[i], " | ", cc$setSize[i], " | ", fmt(cc$NES[i]), " | ", sci(cc$p.adjust[i]), " |")
    add("")
  }
}
writeLines(mk, file.path(outdir, "SUMMARY.md"))
message("\nWrote per-category enrichment to: ", outdir)
cat(sprintf("Categories run: %d | total significant terms: %d\n",
            length(cats), if (nrow(combined)) nrow(combined) else 0))
