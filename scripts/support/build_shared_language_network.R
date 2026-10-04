#!/usr/bin/env Rscript
# ============================================================================
# build_shared_language_network.R
# ----------------------------------------------------------------------------
# Node + edge tables for the Fig5 "shared language" tripartite network:
#   observational lifestyle exposures + interventional trials  (sources, left)
#     -> shared plasma proteins, split into                    (hubs, middle)
#          CAUSAL INTERMEDIATES  (forward protein->disease MR)  and
#          DISEASE REPORTERS     (reverse disease->protein, a marker)
#     -> cardiometabolic disease                               (hub, right)
#
# Thesis (see project_heap_narrative): lifestyle and trials move a shared set of
# proteins; a minority are genetically causal for disease (colocalized / MR),
# the majority are downstream reporters that carry the record of disease.
#
# Disease set = the SAME four panel b (fig_m4_panel_b) uses: type-2 diabetes,
# obesity, lipids, hypertension (so the two panels are consistent). Disease
# SEEDING guarantees each of the four appears (>=1 linked protein).
#
# GENETIC EDGES come from the DISEASE-RESOLVED intervention_mr_edges.tsv (the
# full protein x disease MR edge list) NOT the scatter's single best_disease, so
# a protein shows ALL its causal disease links and diseases like lipids are not
# lost to a stronger competing disease. Re-run annotate_mr.R first if the MR
# tables are stale.
#
# Inputs (support/intervention_compare/):
#   intervention_scatter_mr.tsv   obs beta_HEAP, HERITAGE/GLP1 effects, breadth
#   intervention_mr_edges.tsv     protein x disease x mr_edge_sig x beta_edge
# Output: support/intervention_compare/shared_language_network_{nodes,edges}.tsv
# Run   : HEAP_PATHS_FILE=.../00_paths.R \
#           Rscript scripts/support/build_shared_language_network.R
# ============================================================================
local({
  cand <- c(file.path(getwd(), "scripts", "visualizations", "common"),
            "/n/groups/patel/shakson_ukb/HEAP/scripts/visualizations/common")
  cm <- cand[dir.exists(cand)][1]
  for (f in c("figure_paths.R", "load_heap_results.R")) source(file.path(cm, f))
})
suppressPackageStartupMessages({ library(data.table) })

outdir <- heap_project_output("support", "intervention_compare")
d <- fread(file.path(outdir, "intervention_scatter_mr.tsv"))
d <- d[!is.na(protein) & protein != ""]
e  <- fread(file.path(outdir, "intervention_mr_edges.tsv"))    # disease-resolved edges (used for DP reporters)
tt <- fread(file.path(outdir, "mr_pd_tiered.tsv"))             # CANONICAL Fig-4 tier table (forward causal)

## ---- disease set = panel b's four; robust to finngen codes OR readable names -
dzc <- function(x) { xl <- tolower(x); data.table::fcase(
  grepl("t2d|type.?2.?diab|non.?insulin|\\be11\\b", xl),                 "Type-2 diabetes",
  grepl("obes|\\be66\\b", xl),                                            "Obesity",
  grepl("lipoprot|lipidaem|lipidem|hyperchol|\\be78\\b", xl),            "Lipid disorder",
  grepl("hyptens|hypertens|\\bi10\\b", xl),                              "Hypertension",
  default = NA_character_) }
dcanon <- c("Type-2 diabetes", "Obesity", "Lipid disorder", "Hypertension")

## ---- curated observational exemplars (label -> id regex) --------------------
selx <- c("Strenuous sports"  = "strenuous",
          "Processed meat"     = "processed_meat",
          "Usual walking pace" = "usual_walking_pace",
          "Current smoking"    = "current_tobacco_smoking")
d[, exlab := NA_character_]
for (nm in names(selx)) d[grepl(selx[nm], exposure_id, ignore.case = TRUE), exlab := nm]
# multi-level (one-hot) exposures give several rows per (protein, exposure); collapse
# to ONE edge at the protein's STRONGEST level (largest |beta|) -> no double-counting.
obs  <- d[!is.na(exlab) & is.finite(beta_HEAP), .(protein, exlab, ob = beta_HEAP)]
obs  <- obs[order(protein, exlab, -abs(ob))][, .SD[1], by = .(protein, exlab)]
intv <- unique(d[, .(protein, H = HERITAGE_effect, G = GLP1_effect1)])[!is.na(H) | !is.na(G)]
br   <- d[, .(breadth = uniqueN(exposure_id)), by = protein]

## ---- genetic edges ----------------------------------------------------------
# FORWARD causal edges = Tier 1 or above (cis-pQTL) from the CANONICAL tier table
# (identical gate to Fig 4): tier_rank >= 4 AND edge_class == "cis". coloc-
# confirmed (PP.H4>=0.8) -> "colocalized" (solid); cis coloc-pending -> "cis
# (coloc pending)" (dashed). This EXCLUDES trans-pQTL (Suggestive), LD-confounded
# cis (Tier2) and non-significant cis (Null) from the causal set.
tt[, dz := dzc(disease)]
fwd <- tt[!is.na(dz) & !is.na(beta) & tier_rank >= 4L & edge_class == "cis"]
fwd <- fwd[order(protein, dz, -coloc_confirmed, -tier_rank)][, .SD[1], by = .(protein, dz)]
fwd[, tier := fifelse(coloc_confirmed == TRUE, "colocalized", "cis (coloc pending)")]
fwd <- fwd[, .(protein, dz, tier, sign = as.integer(sign(beta)), weight = abs(beta))]
# REVERSE edges = disease->protein (DP), strongest per protein-disease
e[, dz := dzc(disease)]
rev_all <- e[!is.na(dz) & !is.na(beta_edge) & mr_edge_sig == "DP"][order(protein, dz, -abs(beta_edge))][, .SD[1], by = .(protein, dz)]
rev_all <- rev_all[, .(protein, dz, sign = as.integer(sign(beta_edge)), weight = abs(beta_edge))]

## ---- corroborated core: ALL Tier1 causal + top reporters by breadth ---------
readmoved <- Reduce(intersect, list(obs$protein, intv$protein))
causal_prot <- intersect(unique(fwd$protein), readmoved)
fwd <- fwd[protein %in% causal_prot]
# reporters = non-causal proteins with a DP edge, that are read + moved by interventions
rev <- rev_all[protein %in% readmoved & !protein %in% causal_prot]
rep_pool <- merge(unique(rev[, .(protein)]), br, by = "protein")[order(-breadth)]
selp <- head(unique(c(causal_prot, rep_pool$protein)), 18)     # keep ALL Tier1 causal, then fill reporters

sel <- merge(data.table(protein = selp, class = fifelse(selp %in% causal_prot, "causal", "reporter")), br, by = "protein")
# per-protein exposome R2 (Module 1 base/lasso variance decomposition) -> node colour
r2fp <- heap_project_output("module1_predictive_r2_score_partition", "M1_base_lasso", "base", "lasso", "module1_component_r2.tsv")
r2 <- if (file.exists(r2fp)) fread(r2fp)[, .(protein, R2_E)] else data.table(protein = character(), R2_E = numeric())
sel <- merge(sel, r2, by = "protein", all.x = TRUE)
# causal proteins: ALL Tier1 forward edges (solid, disease-SPECIFIC) + their OWN
# reverse edges (rev_cau, drawn FAINT) -- they too are broad obesity/T2D markers;
# the distinguishing feature is the specific causal edge.
# reporters: ALL reverse edges (dotted), broad multi-disease markers.
fwd_cau <- fwd[protein %in% sel[class == "causal"]$protein]
rev_rep <- rev[protein %in% sel[class == "reporter"]$protein]
rev_cau <- rev_all[protein %in% sel[class == "causal"]$protein]

## ---- NODES ------------------------------------------------------------------
src <- data.table(id = c(names(selx), "HERITAGE", "GLP1-RA"),
                  kind = c(rep("exp_obs", 4), "exp_rct", "exp_rct"))
dzset <- intersect(dcanon, Reduce(union, list(fwd_cau$dz, rev_rep$dz, rev_cau$dz)))
nodes <- rbindlist(list(
  src[, .(id, kind, label = id, class = NA_character_, breadth = NA_integer_, R2_E = NA_real_)],
  sel[, .(id = protein, kind = "protein", label = protein, class, breadth, R2_E)],
  data.table(id = dzset, kind = "disease", label = dzset, class = NA_character_, breadth = NA_integer_, R2_E = NA_real_)))

## ---- EDGES ------------------------------------------------------------------
sgn <- function(x) fifelse(is.na(x) | x == 0, NA_integer_, as.integer(sign(x)))
edges <- rbindlist(list(
  obs[protein %in% selp, .(from = exlab, to = protein, etype = "obs", tier = NA_character_, sign = sgn(ob), weight = abs(ob))],
  intv[protein %in% selp & !is.na(H), .(from = "HERITAGE", to = protein, etype = "interv", tier = NA_character_, sign = sgn(H), weight = abs(H))],
  intv[protein %in% selp & !is.na(G), .(from = "GLP1-RA",  to = protein, etype = "interv", tier = NA_character_, sign = sgn(G), weight = abs(G))],
  fwd_cau[, .(from = protein, to = dz, etype = "gen_fwd", tier, sign, weight)],
  rev_rep[, .(from = dz, to = protein, etype = "gen_rev", tier = "reporter (reverse)", sign, weight)],
  rev_cau[, .(from = dz, to = protein, etype = "gen_rev_causal", tier = "reporter (reverse)", sign, weight)]),
  use.names = TRUE)

fwrite(nodes, file.path(outdir, "shared_language_network_nodes.tsv"), sep = "\t")
fwrite(edges, file.path(outdir, "shared_language_network_edges.tsv"), sep = "\t")
cat(sprintf("shared-language network: %d proteins (%d causal, %d reporter), %d diseases [%s], %d sources\n",
            nrow(sel), sum(sel$class == "causal"), sum(sel$class == "reporter"),
            length(dzset), paste(dzset, collapse = ", "), nrow(src)))
cat(sprintf("edges obs=%d interv=%d gen_fwd=%d gen_rev=%d\n",
            nrow(edges[etype == "obs"]), nrow(edges[etype == "interv"]),
            nrow(edges[etype == "gen_fwd"]), nrow(edges[etype == "gen_rev"])))
