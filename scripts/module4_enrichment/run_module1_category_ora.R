#!/usr/bin/env Rscript
# ============================================================================
# run_module1_category_ora.R
# ----------------------------------------------------------------------------
# PER-CATEGORY over-representation analysis (ORA) of the Module 1 exposomic
# decomposition. For each exposure CATEGORY, take the proteins it shapes
# (unique out-of-fold R2 > THR; score_unique_drop, level=exposure_categories,
# mean over folds) and test which KEGG pathways / GTEx tissues are over-
# represented vs the measured-proteome background (hypergeometric / phyper).
#
# WHY ORA not GSEA: the per-protein clusterProfiler/ReactomePA env is not
# installed here (compute nodes have no package network); ORA needs only the
# CACHED gene-set Term2Gene files + base R, and answers the same question
# without the diet/exercise aggregation confound of run_module1_r2_gsea.R.
#
# Inputs (cached): output/module4_enrichment/{OlinkEntrezConv.txt,
#   genesets/KEGG_T2G.txt, genesets/KEGG_T2N.txt, genesets/GTEX_tissue.txt}
# Output: docs/manuscript_stats/module1_enrichment_bycategory/
#   ora_kegg_bycategory.tsv, ora_gtex_bycategory.tsv, SUMMARY.md
#
# Run (seconds, login node):
#   HEAP_PATHS_FILE=.../workflow/00_paths.R \
#     Rscript scripts/module4_enrichment/run_module1_category_ora.R [covarType] [method]
# ============================================================================
local({
  cand <- c(Sys.getenv("HEAP_PATHS_FILE", unset = ""),
            file.path(getwd(), "workflow", "00_paths.R"),
            "/n/groups/patel/shakson_ukb/HEAP/workflow/00_paths.R")
  hit <- cand[nzchar(cand) & file.exists(cand)][1]; if (!is.na(hit)) source(hit)
})
if (!exists("heap_project_output")) stop("could not load workflow/00_paths.R (set HEAP_PATHS_FILE).")
source(file.path(heap_path(), "scripts", "visualizations", "common", "figure_paths.R"))
source(file.path(heap_path(), "scripts", "visualizations", "common", "load_heap_results.R"))
suppressPackageStartupMessages(library(data.table))

a <- commandArgs(trailingOnly = TRUE); a <- a[!startsWith(a, "--")]
covarType <- if (length(a) >= 1) a[1] else "base"
method    <- if (length(a) >= 2) a[2] else "lasso"
THR       <- as.numeric(Sys.getenv("HEAP_R2_THR", "0.005"))
MIN_REACH <- as.integer(Sys.getenv("HEAP_MIN_REACH", "10"))
MIN_K     <- 5L     # min gene-set size within universe
MIN_OV    <- 3L     # min foreground overlap to test

M4 <- heap_project_output("module4_enrichment")
GS <- file.path(M4, "genesets")

# ---- per-category R2 + universe ---------------------------------------------
ec <- load_module1_predictive_r2(covarType, method, level = "exposure_categories")
if ("method" %in% names(ec)) ec <- ec[get("method") == "score_unique_drop"]
pp <- ec[, .(r2 = mean(r2)), by = .(omic, category)]
universe_sym <- unique(pp$omic)
reach <- pp[, .(reach = sum(r2 > THR)), by = category][order(-reach)]
cats  <- reach[reach >= MIN_REACH, as.character(category)]

# ---- symbol -> entrez (for KEGG) --------------------------------------------
conv <- fread(file.path(M4, "OlinkEntrezConv.txt"))
sym2ent <- unique(conv[!is.na(entrezgene_id) & nzchar(Gene),
                       .(Gene, entrez = as.character(entrezgene_id))], by = "Gene")

# ---- gene sets restricted to the universe -----------------------------------
kegg  <- fread(file.path(GS, "KEGG_T2G.txt"));  setnames(kegg, c("setid", "id"))
kegg[, id := as.character(id)]
keggN <- fread(file.path(GS, "KEGG_T2N.txt"));  setnames(keggN, c("pname", "setid"))
gtex  <- fread(file.path(GS, "GTEX_tissue.txt")); setnames(gtex, c("setid", "id"))

uni_ent <- unique(sym2ent[Gene %in% universe_sym, entrez])
uni_sym <- intersect(universe_sym, unique(gtex$id))
kegg_u  <- unique(kegg[id %in% uni_ent], by = c("setid", "id"))
gtex_u  <- unique(gtex[id %in% uni_sym], by = c("setid", "id"))

# generic hypergeometric ORA: foreground ids vs sets within a universe
ora <- function(fg, set_dt, universe_ids) {
  N <- length(universe_ids); n <- length(fg)
  sizes <- set_dt[, .(K = .N), by = setid]
  ov <- set_dt[id %in% fg, .(k = .N), by = setid]
  r <- merge(sizes, ov, by = "setid", all.x = TRUE); r[is.na(k), k := 0L]
  r <- r[K >= MIN_K & k >= MIN_OV]
  if (!nrow(r)) return(r[0])
  r[, p := phyper(k - 1L, K, N - K, n, lower.tail = FALSE)]
  r[, padj := p.adjust(p, "BH")]
  r[, `:=`(n_fg = n, N_uni = N)]
  r[order(padj)]
}

kegg_all <- list(); gtex_all <- list()
for (cat in cats) {
  fg_sym <- pp[as.character(category) == cat & r2 > THR, omic]
  fg_ent <- unique(sym2ent[Gene %in% fg_sym, entrez])
  rk <- ora(fg_ent, kegg_u, uni_ent)
  if (nrow(rk)) { rk[keggN, pname := i.pname, on = "setid"]; rk[, category := cat]; kegg_all[[cat]] <- rk }
  rg <- ora(intersect(fg_sym, uni_sym), gtex_u, uni_sym)
  if (nrow(rg)) { rg[, `:=`(pname = setid, category = cat)]; gtex_all[[cat]] <- rg }
  message(sprintf("%-22s reach=%4d | KEGG sig(padj<.05)=%2d | GTEx sig=%2d",
                  cat, reach[as.character(category) == cat, reach],
                  if (nrow(rk)) sum(rk$padj < 0.05) else 0,
                  if (nrow(rg)) sum(rg$padj < 0.05) else 0))
}
KEGG <- if (length(kegg_all)) rbindlist(kegg_all, fill = TRUE) else data.table()
GTEX <- if (length(gtex_all)) rbindlist(gtex_all, fill = TRUE) else data.table()

outdir <- file.path(heap_path(), "docs", "manuscript_stats", "module1_enrichment_bycategory")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
ordc <- function(dt) if (nrow(dt)) dt[, .(category, set = pname, setid, K, k, n_fg, N_uni, p, padj)][order(category, padj)] else dt
if (nrow(KEGG)) fwrite(ordc(KEGG), file.path(outdir, "ora_kegg_bycategory.tsv"), sep = "\t")
if (nrow(GTEX)) fwrite(ordc(GTEX), file.path(outdir, "ora_gtex_bycategory.tsv"), sep = "\t")

# ---- SUMMARY.md -------------------------------------------------------------
sci <- function(x) formatC(x, format = "e", digits = 1)
mk <- character(0); add <- function(...) mk <<- c(mk, paste0(...))
add("# Module 1 — PER-CATEGORY ORA (hypergeometric) of exposure-responsive proteins"); add("")
add("_`run_module1_category_ora.R` (", covarType, "/", method, "). Foreground = proteins with ",
    "unique OOF R2 > ", THR, " per category; background = ", length(universe_sym),
    " measured proteins. Hypergeometric over cached KEGG / GTEx gene sets; BH per category. ",
    "Categories with reach >= ", MIN_REACH, " only. (ORA, not GSEA: clusterProfiler env unavailable.)_"); add("")
add("## Reach + significant sets (padj<0.05)"); add("")
add("| category | reach | KEGG sig | GTEx sig |"); add("|---|--:|--:|--:|")
for (i in seq_len(nrow(reach))) {
  ct <- as.character(reach$category[i]); run <- ct %in% cats
  nk <- if (nrow(KEGG)) sum(KEGG$category == ct & KEGG$padj < 0.05) else 0
  ng <- if (nrow(GTEX)) sum(GTEX$category == ct & GTEX$padj < 0.05) else 0
  add("| ", ct, " | ", reach$reach[i], " | ", if (run) nk else "-", " | ", if (run) ng else "-", " |")
}
add("")
for (cat in cats) {
  add("### ", cat); add("")
  kk <- if (nrow(KEGG)) KEGG[category == cat & padj < 0.05][order(padj)] else data.table()
  add("**KEGG pathways** (top 8):", if (!nrow(kk)) " none" else ""); add("")
  if (nrow(kk)) { add("| pathway | k/K | padj |"); add("|---|--:|--:|")
    for (i in seq_len(min(8, nrow(kk)))) add("| ", kk$pname[i], " | ", kk$k[i], "/", kk$K[i], " | ", sci(kk$padj[i]), " |"); add("") }
  gg <- if (nrow(GTEX)) GTEX[category == cat & padj < 0.05][order(padj)] else data.table()
  add("**GTEx tissues** (top 6):", if (!nrow(gg)) " none" else ""); add("")
  if (nrow(gg)) { add("| tissue | k/K | padj |"); add("|---|--:|--:|")
    for (i in seq_len(min(6, nrow(gg)))) add("| ", gg$pname[i], " | ", gg$k[i], "/", gg$K[i], " | ", sci(gg$padj[i]), " |"); add("") }
}
writeLines(mk, file.path(outdir, "SUMMARY.md"))
message("\nWrote per-category ORA to: ", outdir)
