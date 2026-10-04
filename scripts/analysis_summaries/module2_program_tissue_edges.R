#!/usr/bin/env Rscript
# ============================================================================
# module2_program_tissue_edges.R   (ANALYSIS -- run with the GSEA group library)
# ----------------------------------------------------------------------------
# Curate the exposure -> biological program -> tissue associations behind the
# Module-2 tripartite panel (Fig 2b). Reads module4_enrichment/HEAPgsea.qs
# (per-exposure pathway + tissue GSEA), maps pathways to coarse program clusters
# (heap_program_cluster) and tissues to organ systems, links a program to a
# tissue when their leading-edge proteins overlap (>=3 shared, same NES sign),
# and writes three tables consumed by the plotter:
#   exposure_program_edges.tsv  -- exemplar exposure -> program (npath, net dir)
#   program_tissue_edges.tsv    -- program -> tissue (n_exp, up/down) across all
#   cluster_membership.tsv      -- pathway -> cluster audit (prevalence, dir)
#
# Run: module load gcc/14.2.0 R/4.4.2; export HEAP_PATHS_FILE=.../00_paths.R
#      R_LIBS=/n/groups/patel/shakson_ukb/Rlib/gsea_env \
#        Rscript scripts/analysis_summaries/module2_program_tissue_edges.R
# ============================================================================
suppressPackageStartupMessages({ library(qs); library(data.table) }); source(Sys.getenv("HEAP_PATHS_FILE"))
local({ c <- "/n/groups/patel/shakson_ukb/HEAP/scripts/visualizations/common"
  for (f in c("figure_paths","label_helpers","program_clusters")) source(file.path(c, paste0(f, ".R"))) })

ED  <- heap_project_output("module4_enrichment")
OUT <- file.path(ED, "program_tissue"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
g <- qread(file.path(ED, "HEAPgsea.qs")); P <- g@HEAPpgsea; Tt <- g@HEAPtgsea
conv <- fread(file.path(ED, "OlinkEntrezConv.txt"))[, .(entrez = as.character(entrezgene_id), Gene)][!is.na(entrez) & entrez != ""]
e2s <- setNames(conv$Gene, conv$entrez)
t2o <- setNames(rep(names(HEAP_PROGRAM_TISSUES), lengths(HEAP_PROGRAM_TISSUES)), unlist(HEAP_PROGRAM_TISSUES))
le  <- function(x) unlist(strsplit(x, "/"))
exps <- names(P)[!sapply(P, is.null)]; nE <- length(exps)

## ---- cluster membership audit ----
PATH <- rbindlist(lapply(exps, function(e){ pe <- as.data.table(as.data.frame(P[[e]]))[p.adjust < 0.05]
  if (!nrow(pe)) return(NULL); pe[, .(e = e, Description, NES)] }))
PATH[, clust := heap_program_cluster(Description)]
mem <- PATH[, .(n_exp = uniqueN(e), pct = round(100 * uniqueN(e)/nE),
                n_up = uniqueN(e[NES > 0]), n_dn = uniqueN(e[NES < 0])), by = .(clust, Description)][order(clust, -n_exp)]
fwrite(mem, file.path(OUT, "cluster_membership.tsv"), sep = "\t")

## ---- program -> tissue (curated across all exposures) ----
LINK <- rbindlist(lapply(exps, function(e){
  pe <- as.data.table(as.data.frame(P[[e]]))[p.adjust < 0.05]
  tt <- as.data.table(as.data.frame(Tt[[e]]))[p.adjust < 0.05 & ID %in% names(t2o)]
  if (!nrow(pe) || !nrow(tt)) return(NULL); pe[, cl := heap_program_cluster(Description)]; rows <- list()
  for (i in seq_len(nrow(pe))) for (j in seq_len(nrow(tt))) {
    if (sign(pe$NES[i]) != sign(tt$NES[j])) next
    sh <- length(intersect(unname(na.omit(e2s[le(pe$core_enrichment[i])])), le(tt$core_enrichment[j]))); if (sh < 3) next
    rows[[length(rows)+1]] <- data.table(e = e, clust = pe$cl[i], organ = t2o[tt$ID[j]], dir = sign(pe$NES[i])) }
  if (length(rows)) rbindlist(rows) }))
PT <- LINK[, .(n_exp = uniqueN(e), n_up = uniqueN(e[dir > 0]), n_dn = uniqueN(e[dir < 0])), by = .(clust, organ)][order(-n_exp)]
fwrite(PT, file.path(OUT, "program_tissue_edges.tsv"), sep = "\t")

## ---- exemplar exposure -> program (one true direction each) ----
EX <- as.data.table(HEAP_TRIPARTITE_EXEMPLARS)
EP <- rbindlist(lapply(seq_len(nrow(EX)), function(k){ e <- EX$e[k]; if (is.null(P[[e]])) return(NULL)
  pe <- as.data.table(as.data.frame(P[[e]]))[p.adjust < 0.05]; if (!nrow(pe)) return(NULL)
  pe[, cl := heap_program_cluster(Description)]
  s <- pe[cl %in% HEAP_PROGRAM_LEVELS, .(npath = .N, net = sum(sign(NES) * abs(NES))), by = cl]; if (!nrow(s)) return(NULL)
  s[, `:=`(exposure = e, lab = EX$lab[k], grp = EX$grp[k], dir = ifelse(net >= 0, "up", "down"))]; s }))
setnames(EP, "cl", "clust"); EP <- EP[npath >= 1]
fwrite(EP, file.path(OUT, "exposure_program_edges.tsv"), sep = "\t")

cat(sprintf("exposures=%d | clusters=%d | program-tissue edges=%d | exposure-program edges=%d\n",
            nE, uniqueN(PATH$clust), nrow(PT), nrow(EP)))
cat("missing exemplars:", paste(setdiff(EX$lab, unique(EP$lab)), collapse = ", "), "\n")
cat("wrote", OUT, "/{exposure_program,program_tissue}_edges.tsv + cluster_membership.tsv\n")
