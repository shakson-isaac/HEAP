#!/usr/bin/env Rscript
# ============================================================================
# build_mr_supp_tables.R
# ----------------------------------------------------------------------------
# The Module-5 MR tables that SHIP: two Supplementary Tables (triad-level) and
# three Supplementary Data files (triad- and edge-level). Replaces a set of
# Python prototypes in HEAPweb -- anything that writes into HEAP/ is R, so the
# supplement builds itself and does not depend on the website repo.
#
# Outputs (docs/manuscript_stats/module5/):
#   mr_triads_wide.tsv    18,780 x 72   Supp Table  S_mr_triads   (motif-carrying)
#   mr_triads_all.tsv    199,646 x 72   Supp Data                 (all tested)
#   mr_edges_ukb.tsv      77,573 x 48   Supp Data   mr_edges      (Olink panel)
#   mr_edges_decode.tsv   58,593 x 48   Supp Data   mr_edges_decode (SomaScan)
#
# THREE THINGS THESE ADD over what shipped before:
#   1. the deCODE panel. Tier 1+ is DEFINED as replication across both pQTL
#      panels, and the old triad table carried only UK Biobank, so the rung
#      could not be checked from it.
#   2. mr_tier_final per edge. The presence flags (pEP ... pDE) are TRUE for
#      Tier 1 OR Tier 1+, collapsing the strongest rung into a boolean.
#   3. effect sizes in the edge table. What shipped as `mr_edges`
#      (mr_tiered_edges.tsv) had 14 columns of tiers and flags and NO beta, SE
#      or p-value -- a reader could see an edge was Tier 1 but not its estimate.
#
# Two naming decisions, applied at this boundary:
#   * edge_class "polygenic" -> "genome-wide". The source word reads as
#     "polygenic score", which these are not, and it conflated exposure- and
#     disease-origin instruments. `instrument_type` names the origin.
#   * coloc_status "pending" is a WORK QUEUE MARKER, not a result
#     (aggregate_coloc.R exists to replace it with a verdict). The verdict is
#     joined from support/coloc/coloc_results.tsv under the PP.H4 >= 0.8 gate.
#
#   Run:  module load gcc/14.2.0 R/4.4.2
#         Rscript scripts/analysis_summaries/build_mr_supp_tables.R
# ============================================================================
suppressPackageStartupMessages(library(data.table))
options(scipen = 999)

local({
  cand <- c(Sys.getenv("HEAP_PATHS_FILE", unset = ""),
            "/n/groups/patel/shakson_ukb/HEAP/workflow/00_paths.R")
  hit <- cand[nzchar(cand) & file.exists(cand)][1]; if (!is.na(hit)) source(hit)
})
SUMD <- file.path(heap_project_output("mr_edges"), "summary")
# Large generated tables go to the IGLOO output root, not into the repo: the
# biggest file here is ~180 MB and everything in docs/manuscript_stats/ is ~1 MB.
# config/supp_tables.tsv already resolves absolute paths, which is how the
# existing mr_edges and coloc_summary rows point at IGLOO.
OUT  <- file.path(SUMD, "supp")
STATS <- file.path(heap_root, "docs", "manuscript_stats", "module5")
COL  <- file.path(heap_project_output("support"), "coloc", "coloc_results.tsv")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
ts <- function(...) message(format(Sys.time(), "%H:%M:%S "), ...)

## ---- edge -> triad mapping -------------------------------------------------
## suffix in the wide layout, its direction, instrument class, and the
## endpoints of the edge within the (E, P, D) triad.
EDGE <- data.table(
  suf   = c("EP","PDcis","PDtrans","ED","PEcis","PEtrans","DP","DE"),
  tdir  = c("E_to_P","Pcis_to_D","Ptrans_to_D","E_to_D",
            "Pcis_to_E","Ptrans_to_E","D_to_P","D_to_E"),
  cls   = c("genome-wide","cis","trans","genome-wide",
            "cis","trans","genome-wide","genome-wide"),
  #        EP  PDcis PDtrans ED  PEcis PEtrans DP  DE
  from  = c("E","P",  "P",    "E","P",  "P",    "D","D"),
  to    = c("P","D",  "D",    "D","E",  "E",    "P","E"))
## deCODE is an alternative PROTEIN panel, so it can recompute any edge with the
## protein at either end. E<->D has no protein and is identical across panels.
EDGE[, protein_edge := suf %in% c("EP","PDcis","PDtrans","PEcis","PEtrans","DP")]

## ---- colocalization verdicts, replacing "pending" --------------------------
coloc <- if (file.exists(COL)) fread(COL) else data.table()
if (nrow(coloc)) {
  coloc[, panel := fifelse(arm == "UKB", "UKB", "dCODE")]
  coloc[, verdict := fcase(is.finite(`PP.H4`) & `PP.H4` >= 0.8, "colocalized",
                           is.finite(`PP.H3`) & `PP.H3` >= 0.8, "distinct_variants",
                           default = "ambiguous")]
  ts("coloc verdicts: ", nrow(coloc))
}

te <- fread(file.path(SUMD, "mr_tiered_edges.tsv"))
te[, panel := fifelse(dataset == "UKB", "UKB", "dCODE")]
setkey(te, panel, edge_dir, src_id, tgt_id)
ts("tiered edges: ", format(nrow(te), big.mark = ","))

## ============================================================================
## 1-2. triad tables (wide), all tested and motif-carrying
## ============================================================================
mot <- fread(file.path(STATS, "mr_triad_motifs.tsv"))   # Tier-1 motif per triad
setkey(mot, Exposure, Protein, Disease)

wide <- fread(file.path(SUMD, "MRmotifs.tsv"))[dataset == "UKB"]
dec  <- fread(file.path(SUMD, "DECODE", "MRmotifs.tsv"))
setkey(dec, Exposure, Protein, Disease)
ts("triads: ", format(nrow(wide), big.mark = ","))

out <- wide[, .(triad_id = paste(Exposure, Protein, Disease, sep = "|"),
                Exposure, Exposure_category = ExposureCategory, Protein,
                Disease, Disease_UKB, ICD10)]
out[, motif := mot[.(wide$Exposure, wide$Protein, wide$Disease), motif]]
out[is.na(motif), motif := ""]

for (i in seq_len(nrow(EDGE))) {
  e <- EDGE[i]
  for (tag in if (e$protein_edge) c("UKB", "dCODE") else "UKB") {
    src <- if (tag == "UKB") wide else dec[.(wide$Exposure, wide$Protein, wide$Disease)]
    for (stat in c("beta", "se", "padj"))
      set(out, j = sprintf("%s_%s_%s", stat, e$suf, tag),
          value = src[[sprintf("%s_%s", stat, e$suf)]])
    nodes <- list(E = wide$Exposure, P = wide$Protein, D = wide$Disease)
    set(out, j = sprintf("tier_%s_%s", e$suf, tag),
        value = te[.(tag, e$tdir, nodes[[e$from]], nodes[[e$to]]), mr_tier_final])
    if (e$suf %in% c("PDcis", "PEcis") && nrow(coloc)) {
      cv <- coloc[.(tag, e$tdir, wide$Protein, nodes[[e$to]]),
                  on = .(panel, edge_dir, protID, target)]
      set(out, j = sprintf("coloc_%s_%s", e$suf, tag),     value = cv$verdict)
      set(out, j = sprintf("colocPPH4_%s_%s", e$suf, tag), value = cv$`PP.H4`)
    }
  }
}
fwrite(out, file.path(OUT, "mr_triads_all.tsv"), sep = "\t", na = "")
ts("wrote mr_triads_all.tsv  ", format(nrow(out), big.mark = ","), " x ", ncol(out))

fwrite(out[motif != ""], file.path(OUT, "mr_triads_wide.tsv"), sep = "\t", na = "")
ts("wrote mr_triads_wide.tsv ", format(nrow(out[motif != ""]), big.mark = ","),
   " x ", ncol(out))

## ============================================================================
## 3-4. per-panel edge tables
## ============================================================================
## mr_sensitivity_long.tsv has the estimate and every diagnostic but not the
## final tier; mr_tiered_edges.tsv has the final tier and the cross-panel
## replication flag but no estimate. Same rows. Join them.
for (p in list(list(tag = "UKB", dir = ".", label = "UK Biobank",
                    plat = "Olink", file = "mr_edges_ukb.tsv"),
               list(tag = "dCODE", dir = "DECODE", label = "deCODE",
                    plat = "SomaScan", file = "mr_edges_decode.tsv"))) {
  s <- fread(file.path(SUMD, p$dir, "mr_sensitivity_long.tsv"))
  v <- te[panel == p$tag, .(edge_dir, src_id, tgt_id, mr_tier_final, replicated)]
  s <- merge(s, v, by = c("edge_dir", "src_id", "tgt_id"), all.x = TRUE, sort = FALSE)
  s[, `:=`(panel = p$label, platform = p$plat,
           instrument_type = fcase(edge_class == "cis",   "cis-pQTL",
                                   edge_class == "trans", "trans-pQTL",
                                   startsWith(edge_dir, "E_"), "exposure GWAS",
                                   default = "disease GWAS"))]
  setcolorder(s, c("panel", "platform", "instrument_type"))
  fwrite(s, file.path(OUT, p$file), sep = "\t", na = "")
  ts("wrote ", p$file, "  ", format(nrow(s), big.mark = ","), " x ", ncol(s))
}
ts("done")
