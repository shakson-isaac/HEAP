#!/usr/bin/env Rscript

# ============================================================================
# annotate_mr.R  —  HEAP support analysis (intervention <-> MR integration)
# ----------------------------------------------------------------------------
# Annotates the intervention-comparison scatter (HEAP exposure->protein effects
# vs HERITAGE / GLP1 protein shifts) with Mendelian-randomization causal-edge
# evidence, so the Module-4 figures can colour each protein by the significant
# causal edge and shape it by UKB<->DECODE replication.
#
# Ported from the analysis core of the legacy
#   Visualizations/ModuleMR/COMBO/compare_GLP1v2.R
# but: (1) operates on the canonical intervention_scatter.tsv (per exposure-term
# x protein) instead of the dead HEAPintv2.qs S4 object; (2) keeps the analysis
# OUT of the plotter (tables are emitted; the figure reads them).
#
# Two grains are emitted because the figures use both:
#   * intervention_mr_edges.tsv   DISEASE-RESOLVED (mr_key x protein x disease):
#       the per-triplet edge call used by the disease-specific scatter panels.
#       Replication ("Both") here means the SAME edge type is significant in
#       BOTH UKB and DECODE for the SAME (exposure, protein, disease) triplet
#       (the legacy trip_key rule).
#   * intervention_scatter_mr.tsv PER (exposure-term x protein): the strongest
#       edge across diseases (any-disease overview) + the driving disease.
#
# MR edge types (protein->disease causal direction is PDcis / PDtrans):
#   PDcis    protein->disease, cis-instrumented        (most specific)
#   PDtrans  protein->disease, trans-instrumented
#   DP       disease->protein (reverse)                (weakest / most common)
# Per-point colour priority: PDcis > PDtrans > DP > None.
#
# Inputs:
#   support/intervention_compare/intervention_scatter.tsv   (run_intervention_compare.R)
#     needs cols: covarType, exposure_id, Eid, Category, protein, beta_HEAP,
#                 se_HEAP, HERITAGE_effect, GLP1_effect1, GLP1_effect2, olink_soma_r
#   MRmotifs (UKB)     HEAP_MRMOTIFS_UKB    or legacy summary/MRmotifs.csv
#   MRmotifs (DECODE)  HEAP_MRMOTIFS_DECODE or legacy summary/DECODE/MRmotifs.csv
#     needs cols: Exposure, Protein, Disease, padj_{PDcis,PDtrans,DP},
#                 beta_{PDcis,PDtrans,DP}, any_sig
#
# Outputs -> support/intervention_compare/ (+ figures/data/):
#   intervention_mr_edges.tsv    mr_key, protein, disease, mr_edge_sig,
#                                mr_support, padj_edge, beta_edge
#   intervention_scatter_mr.tsv  covarType, exposure_id, Eid, Category, mr_key,
#                                protein, beta_HEAP, se_HEAP, HERITAGE_effect,
#                                GLP1_effect1, GLP1_effect2, olink_soma_r,
#                                mr_edge_sig, mr_support, best_disease,
#                                padj_edge, beta_edge, n_dz_edge
#
# Usage:
#   module load gcc/14.2.0 R/4.4.2
#   HEAP_PATHS_FILE=.../workflow/00_paths.R \
#     Rscript scripts/support/intervention_compare/annotate_mr.R
# ============================================================================

local({
  cand <- c(Sys.getenv("HEAP_PATHS_FILE", unset = ""),
            file.path(getwd(), "workflow", "00_paths.R"),
            "/n/groups/patel/shakson_ukb/HEAP/workflow/00_paths.R")
  hit <- cand[nzchar(cand) & file.exists(cand)][1]
  if (is.na(hit)) stop("Could not locate workflow/00_paths.R (set HEAP_PATHS_FILE).")
  source(hit)
})
suppressPackageStartupMessages({ library(data.table) })

MR_ALPHA <- 0.05
EDGES    <- c("PDcis", "PDtrans", "DP")   # colour priority is reverse of this
EDGE_RANK <- c(PDcis = 3L, PDtrans = 2L, DP = 1L, None = 0L)

# --- locate inputs ----------------------------------------------------------
sc_path <- heap_project_output("support", "intervention_compare",
                               "intervention_scatter.tsv")
if (!file.exists(sc_path))
  stop("Missing ", sc_path, " — run run_intervention_compare.R first.", call. = FALSE)

# canonical Module 5 motif tables are now produced by
# scripts/support/mr_tables/build_mr_tables.R from the fresh per-edge tree, so
# default to those (override via env). Falls back to the legacy CSVs if absent.
.mr_default <- function(...) {
  fresh <- heap_project_output("mr_edges", "summary", ...)
  legacy <- file.path("/n/groups/patel/shakson_ukb/UK_Biobank/Output/MR_edges/summary",
                      sub("\\.tsv$", ".csv", file.path(...)))
  if (file.exists(fresh)) fresh else legacy
}
mr_ukb_fp <- Sys.getenv("HEAP_MRMOTIFS_UKB", unset = .mr_default("MRmotifs.tsv"))
mr_dec_fp <- Sys.getenv("HEAP_MRMOTIFS_DECODE", unset = .mr_default("DECODE", "MRmotifs.tsv"))
for (p in c(mr_ukb_fp, mr_dec_fp))
  if (!file.exists(p)) stop("Missing MRmotifs input: ", p, call. = FALSE)

sc <- fread(sc_path)
message("intervention_scatter: ", nrow(sc), " rows, ",
        uniqueN(sc$exposure_id), " exposure-terms, ", uniqueN(sc$protein), " proteins.")
sc_proteins <- unique(sc$protein)

# --- read one MRmotifs file -> per-triplet edge significance flags -----------
read_edges <- function(fp, tag) {
  keep <- c("Exposure", "Protein", "Disease", "any_sig",
            paste0("padj_", EDGES), paste0("beta_", EDGES))
  x <- fread(fp, select = keep)
  x <- x[any_sig == TRUE]
  x[, any_sig := NULL]
  x <- x[Protein %in% sc_proteins]               # lean: only proteins we plot
  for (e in EDGES)
    x[, (paste0("sig_", e)) := is.finite(get(paste0("padj_", e))) &
                               get(paste0("padj_", e)) < MR_ALPHA]
  vcols <- c(paste0("padj_", EDGES), paste0("beta_", EDGES), paste0("sig_", EDGES))
  setnames(x, vcols, paste0(vcols, "_", tag))
  x
}

# universe of MR exposures (unfiltered) for the ID-vs-Eid key choice
mr_exposures <- union(fread(mr_ukb_fp, select = "Exposure")$Exposure,
                      fread(mr_dec_fp, select = "Exposure")$Exposure)

message("Reading UKB MRmotifs ...")
ukb <- read_edges(mr_ukb_fp, "ukb")
message("Reading DECODE MRmotifs ...")
dec <- read_edges(mr_dec_fp, "dec")

# --- disease-resolved combine (trip_key = Exposure x Protein x Disease) ------
edges <- merge(ukb, dec, by = c("Exposure", "Protein", "Disease"), all = TRUE)
for (e in EDGES) for (tg in c("ukb", "dec")) {
  col <- paste0("sig_", e, "_", tg)
  edges[is.na(get(col)), (col) := FALSE]
}
# combined significance + strongest edge (priority PDcis > PDtrans > DP)
for (e in EDGES)
  edges[, (paste0("sigc_", e)) := get(paste0("sig_", e, "_ukb")) |
                                  get(paste0("sig_", e, "_dec"))]
edges[, mr_edge_sig := fifelse(sigc_PDcis, "PDcis",
                        fifelse(sigc_PDtrans, "PDtrans",
                         fifelse(sigc_DP, "DP", "None")))]
edges <- edges[mr_edge_sig != "None"]            # any_sig may flag EP/ED-only rows
edges[, mr_edge_sig := factor(mr_edge_sig, levels = c("None", "DP", "PDtrans", "PDcis"))]

# replication of the CHOSEN edge type at this triplet -> shape
edges[, sig_u := fifelse(mr_edge_sig == "PDcis",   sig_PDcis_ukb,
                  fifelse(mr_edge_sig == "PDtrans", sig_PDtrans_ukb, sig_DP_ukb))]
edges[, sig_d := fifelse(mr_edge_sig == "PDcis",   sig_PDcis_dec,
                  fifelse(mr_edge_sig == "PDtrans", sig_PDtrans_dec, sig_DP_dec))]
edges[, mr_support := fifelse(sig_u & sig_d, "Both",
                       fifelse(sig_u, "UKB only",
                        fifelse(sig_d, "DECODE only", "None")))]
edges[, mr_support := factor(mr_support, levels = c("None", "UKB only", "DECODE only", "Both"))]

# padj / beta for the chosen edge type (dataset with the smaller padj) ---------
edges[, `:=`(padj_edge = NA_real_, beta_edge = NA_real_)]
for (e in EDGES) {
  sel <- which(edges$mr_edge_sig == e)
  if (!length(sel)) next
  pu <- edges[[paste0("padj_", e, "_ukb")]][sel]; pd <- edges[[paste0("padj_", e, "_dec")]][sel]
  bu <- edges[[paste0("beta_", e, "_ukb")]][sel]; bd <- edges[[paste0("beta_", e, "_dec")]][sel]
  use_u <- is.finite(pu) & (!is.finite(pd) | pu <= pd)
  set(edges, sel, "padj_edge", ifelse(use_u, pu, pd))
  set(edges, sel, "beta_edge", ifelse(use_u, bu, bd))
}

edges_out <- edges[, .(mr_key = Exposure, protein = Protein, disease = Disease,
                       mr_edge_sig, mr_support, padj_edge, beta_edge)]

# --- any-disease aggregate per (mr_key, protein): strongest edge, best disease
edges_out[, erank := EDGE_RANK[as.character(mr_edge_sig)]]
agg <- edges_out[order(-erank, padj_edge)][, .SD[1L],
                 by = .(mr_key, protein)][
                 , .(mr_key, protein, mr_edge_sig, mr_support,
                     best_disease = disease, padj_edge, beta_edge)]
# n_dz_edge: how many diseases carry the strongest edge type for this protein
ndz <- edges_out[, .(n_dz_edge = .N),
                 by = .(mr_key, protein, mr_edge_sig)]
agg <- merge(agg, ndz, by = c("mr_key", "protein", "mr_edge_sig"), all.x = TRUE)

# --- annotate the scatter (any-disease overview) ----------------------------
sc[, mr_key := fifelse(exposure_id %in% mr_exposures, exposure_id,
                fifelse(Eid %in% mr_exposures, Eid, exposure_id))]
sc_mr <- merge(sc, agg, by = c("mr_key", "protein"), all.x = TRUE)
sc_mr[is.na(mr_edge_sig),  mr_edge_sig  := "None"]
sc_mr[is.na(mr_support),   mr_support   := "None"]
sc_mr[is.na(n_dz_edge),    n_dz_edge    := 0L]
sc_mr[, mr_edge_sig := factor(mr_edge_sig, levels = c("None", "DP", "PDtrans", "PDcis"))]
sc_mr[, mr_support  := factor(mr_support,  levels = c("None", "UKB only", "DECODE only", "Both"))]

out_cols <- c("covarType", "exposure_id", "Eid", "Category", "mr_key", "protein",
              "beta_HEAP", "se_HEAP", "HERITAGE_effect", "GLP1_effect1",
              "GLP1_effect2", "olink_soma_r",
              "mr_edge_sig", "mr_support", "best_disease",
              "padj_edge", "beta_edge", "n_dz_edge")
out_cols <- out_cols[out_cols %in% names(sc_mr)]
sc_mr <- sc_mr[, ..out_cols]

# --- write ------------------------------------------------------------------
out_dir  <- heap_project_output("support", "intervention_compare")
fig_data <- file.path(heap_project_root("figures"), "data")
dir.create(out_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(fig_data, recursive = TRUE, showWarnings = FALSE)
fwrite(edges_out[, .(mr_key, protein, disease, mr_edge_sig, mr_support, padj_edge, beta_edge)],
       file.path(out_dir,  "intervention_mr_edges.tsv"), sep = "\t")
fwrite(edges_out[, .(mr_key, protein, disease, mr_edge_sig, mr_support, padj_edge, beta_edge)],
       file.path(fig_data, "intervention_mr_edges.tsv"), sep = "\t")
fwrite(sc_mr, file.path(out_dir,  "intervention_scatter_mr.tsv"), sep = "\t")
fwrite(sc_mr, file.path(fig_data, "intervention_scatter_mr.tsv"), sep = "\t")

# --- summary ----------------------------------------------------------------
cat("\n================ MR ANNOTATION SUMMARY ================\n")
cat(sprintf("disease-resolved edges (proteins in scatter): %d triplets, %d (exposure,protein,disease)\n",
            nrow(edges_out), uniqueN(edges_out[, .(mr_key, protein, disease)])))
cat("\nDisease-resolved edge calls (mr_edge_sig x mr_support):\n")
print(dcast(edges_out, mr_edge_sig ~ mr_support, value.var = "protein",
            fun.aggregate = length))
cat("\nPer-(exposure,protein) overview (strongest edge):\n")
print(sc_mr[, .N, by = mr_edge_sig][order(-N)])
cat("\nReplicated (Both) protein->disease triplets by disease (top 15):\n")
print(head(edges_out[mr_edge_sig %in% c("PDcis", "PDtrans") & mr_support == "Both",
                     .N, by = disease][order(-N)], 15))
cat("\nWrote:\n  ", file.path(out_dir, "intervention_mr_edges.tsv"),
    "\n  ", file.path(out_dir, "intervention_scatter_mr.tsv"), "\nDONE.\n", sep = "")
