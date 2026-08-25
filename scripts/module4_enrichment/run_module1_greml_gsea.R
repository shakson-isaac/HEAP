#!/usr/bin/env Rscript
# ============================================================================
# run_module1_greml_gsea.R
# ----------------------------------------------------------------------------
# Pathway-level concordance of the Module-1 G/E/GxE decomposition between the
# two estimators. Over the SAME proteins (the converged-GREML ∩ Module-1 set),
# rank by (i) GREML variance component and (ii) HEAP unique predictive R2, run
# the identical GSEA engine, and compare which pathways each ranking enriches.
# Answers: "do GREML and HEAP rankings recover the same biology, per component?"
#
# Input  : population_architecture/<covar>/grm_cutoff_<cut>/concordance_greml_vs_heap_r2.tsv
#          (protein, component, greml, r2 -- built by support/module1_greml_vs_r2.R)
# Engine : clusterProfiler::gseGO (SYMBOL, ont=ALL) + ReactomePA::gsePathway
#          (mirror run_module1_r2_gsea.R params).
# Output : docs/manuscript_stats/module1_enrichment_greml/
#            gsea_all_significant.tsv     (metric x component x database, long)
#            COMPARISON.md                (overlap + NES concordance, GREML vs HEAP)
# ============================================================================
local({
  cand <- c(Sys.getenv("HEAP_PATHS_FILE", unset = ""),
            "/n/groups/patel/shakson_ukb/HEAP/workflow/00_paths.R")
  hit <- cand[nzchar(cand) & file.exists(cand)][1]; if (!is.na(hit)) source(hit)
})
source(file.path(heap_path(), "scripts", "visualizations", "common", "figure_paths.R"))
suppressPackageStartupMessages({
  library(data.table); library(clusterProfiler); library(org.Hs.eg.db); library(ReactomePA)
})

covarType <- "base"; grm_cut <- "0p025"
pa_dir <- file.path(heap_project_output("population_architecture"), covarType, paste0("grm_cutoff_", grm_cut))
d <- fread(file.path(pa_dir, "concordance_greml_vs_heap_r2.tsv"))   # protein, component, greml, r2
comps <- c("G", "E", "GxE")

run_go <- function(gl) tryCatch(
  clusterProfiler::gseGO(geneList = gl, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "ALL",
    minGSSize = 10, maxGSSize = 500, pvalueCutoff = 0.05, pAdjustMethod = "BH", eps = 0,
    seed = TRUE, verbose = FALSE), error = function(e) NULL)
run_re <- function(gl) {
  map <- suppressWarnings(clusterProfiler::bitr(names(gl), "SYMBOL", "ENTREZID", org.Hs.eg.db))
  map <- map[!duplicated(map$SYMBOL), ]; gl2 <- gl[map$SYMBOL]; names(gl2) <- map$ENTREZID
  gl2 <- sort(gl2[!is.na(names(gl2))], decreasing = TRUE)
  tryCatch(ReactomePA::gsePathway(geneList = gl2, organism = "human", minGSSize = 10,
    maxGSSize = 500, pvalueCutoff = 0.05, pAdjustMethod = "BH", eps = 0, seed = TRUE,
    verbose = FALSE), error = function(e) NULL)
}
as_dt <- function(g) if (is.null(g) || !nrow(as.data.frame(g))) NULL else as.data.table(as.data.frame(g))

all_sig <- list()
for (metric in c("greml", "r2")) for (comp in comps) {
  sub <- d[component == comp & is.finite(get(metric))]
  gl <- sort(setNames(sub[[metric]], sub$protein), decreasing = TRUE)
  message("== ", metric, " / ", comp, " (", length(gl), " proteins) ==")
  go <- as_dt(run_go(gl)); re <- as_dt(run_re(gl))
  if (!is.null(go)) all_sig[[paste(metric, comp, "GO")]]  <- data.table(metric, component = comp, database = "GO", go)
  if (!is.null(re)) all_sig[[paste(metric, comp, "RE")]]  <- data.table(metric, component = comp, database = "Reactome", ONTOLOGY = NA_character_, re)
  message("   GO sig=", if (is.null(go)) 0 else nrow(go), " | Reactome sig=", if (is.null(re)) 0 else nrow(re))
}
S <- rbindlist(all_sig, fill = TRUE)
outdir <- file.path(heap_path(), "docs", "manuscript_stats", "module1_enrichment_greml")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
keep <- intersect(c("metric","component","database","ONTOLOGY","ID","Description","setSize","NES","pvalue","p.adjust"), names(S))
fwrite(S[, ..keep][order(component, database, metric, p.adjust)], file.path(outdir, "gsea_all_significant.tsv"), sep = "\t")

# ---- concordance: overlap of significant terms + NES correlation -----------
mk <- c("# Module 1 — GREML vs HEAP ranking: pathway concordance", "",
        "_Same proteins (converged-GREML ∩ Module 1). GSEA on each ranking; compares which",
        "gene sets each estimator enriches per component (BH p.adjust<0.05)._", "",
        "| component | database | sig(GREML) | sig(HEAP) | shared | Jaccard | NES r (shared) |",
        "|---|---|--:|--:|--:|--:|--:|")
sumrows <- list()
for (comp in comps) for (db in c("GO", "Reactome")) {
  g <- S[metric == "greml" & component == comp & database == db]
  h <- S[metric == "r2"    & component == comp & database == db]
  ids_g <- unique(g$ID); ids_h <- unique(h$ID); sh <- intersect(ids_g, ids_h)
  sumrows[[paste(comp, db)]] <- data.table(component = comp, database = db,
    n_greml = length(ids_g), n_heap = length(ids_h), n_shared = length(sh))
  jac <- if (length(union(ids_g, ids_h))) length(sh) / length(union(ids_g, ids_h)) else NA_real_
  nesr <- NA_real_
  if (length(sh) >= 3) {
    m <- merge(g[ID %in% sh, .(ID, NES_g = NES)], h[ID %in% sh, .(ID, NES_h = NES)], by = "ID")
    nesr <- suppressWarnings(cor(m$NES_g, m$NES_h))
  }
  mk <- c(mk, sprintf("| %s | %s | %d | %d | %d | %s | %s |", comp, db,
                      length(ids_g), length(ids_h), length(sh),
                      ifelse(is.na(jac), "—", formatC(jac, format="f", digits=2)),
                      ifelse(is.na(nesr), "—", formatC(nesr, format="f", digits=2))))
}
writeLines(mk, file.path(outdir, "COMPARISON.md"))
fwrite(rbindlist(sumrows), file.path(pa_dir, "pathway_concordance_summary.tsv"), sep = "\t")  # thin-plotter input
cat(paste(mk, collapse = "\n"), "\n")
message("\nWrote ", outdir)
