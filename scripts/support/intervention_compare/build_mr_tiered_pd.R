#!/usr/bin/env Rscript
# ============================================================================
# build_mr_tiered_pd.R  — refined protein->disease MR edge table (tier + coloc)
# ----------------------------------------------------------------------------
# Replaces the flat "significant or not" protein->disease calls (intervention_mr_
# edges.tsv) with the new Module-5 TIERED + COLOC output. For each (protein,
# disease) it keeps the BEST protein->disease edge across arms (UKB/deCODE) and
# instrument class (cis/trans), recording confidence:
#   tier  Tier1plus > Tier1 > Tier2 > Suggestive > Reverse > Null
#   coloc cis-pQTL colocalization (PP.H4 >= 0.8 = confirmed)  [gold standard]
# Direction (beta sign) is joined from MRmotifs (beta_PDcis/PDtrans).
#
# Inputs : mr_edges/summary/mr_tiered_edges.tsv  (edge_dir Pcis_to_D / Ptrans_to_D)
#          mr_edges/summary/MRmotifs.tsv         (beta_PDcis / beta_PDtrans)
# Output : support/intervention_compare/mr_pd_tiered.tsv
#          protein, disease, mr_tier, tier_rank, edge_class, coloc_confirmed,
#          coloc_pph4, n_arms_qualified, replicated, beta, sign
# Run    : module load gcc/14.2.0 R/4.4.2; HEAP_PATHS_FILE=.../00_paths.R Rscript .../build_mr_tiered_pd.R
# ============================================================================
local({
  cand <- c(Sys.getenv("HEAP_PATHS_FILE", unset = ""),
            "/n/groups/patel/shakson_ukb/HEAP/workflow/00_paths.R")
  source(cand[nzchar(cand) & file.exists(cand)][1])
})
suppressPackageStartupMessages(library(data.table))

te <- fread(heap_project_output("mr_edges", "summary", "mr_tiered_edges.tsv"))
pd <- te[edge_dir %in% c("Pcis_to_D", "Ptrans_to_D")]
pd[, eclass := fifelse(edge_dir == "Pcis_to_D", "cis", "trans")]
TR <- c(Tier1plus = 5L, Tier1 = 4L, Tier2 = 3L, Suggestive = 2L, Reverse = 1L, Null = 0L)
pd[, trank := TR[mr_tier_final]]; pd[is.na(trank), trank := 0L]
# coloc: mr_tiered_edges' coloc_status is INCOMPLETE (mostly pending/unavailable),
# so use the AUTHORITATIVE systematic coloc table (build_coloc_table.R) — e.g. it
# correctly marks ASGR1→Lipoprotein (PP.H4=0.998) which mr_tiered_edges left "pending".
cf <- file.path(heap_project_output("support","intervention_compare"), "coloc_pph4.tsv")
cz <- if (file.exists(cf)) fread(cf) else data.table(protein=character(), disease=character(), pph4=numeric())
pd[, coloc_status := NULL]
pd <- merge(pd, cz[, .(src_id = protein, tgt_id = disease, coloc_pph4_auth = pph4)],
            by = c("src_id","tgt_id"), all.x = TRUE)
pd[, coloc_pph4 := fifelse(is.finite(coloc_pph4_auth), coloc_pph4_auth, as.numeric(coloc_pph4))]
pd[, colocok := is.finite(coloc_pph4) & coloc_pph4 >= 0.8]

# best edge per (protein, disease): highest tier, then coloc-confirmed, then cis, then replicated
ord <- pd[order(-trank, -colocok, eclass != "cis", -as.integer(replicated))]
best <- ord[, .SD[1], by = .(protein = src_id, disease = tgt_id)]
# how many arms reach >= Tier2 for this protein->disease (cross-cohort support)
qual <- pd[trank >= 3L, .(n_arms_qualified = uniqueN(dataset)), by = .(protein = src_id, disease = tgt_id)]
best <- merge(best, qual, by = c("protein","disease"), all.x = TRUE)
best[is.na(n_arms_qualified), n_arms_qualified := 0L]

# direction from MRmotifs (beta of the chosen instrument class)
mm <- unique(fread(heap_project_output("mr_edges","summary","MRmotifs.tsv"),
                   select = c("Protein","Disease","beta_PDcis","beta_PDtrans")))
best <- merge(best, mm, by.x = c("protein","disease"), by.y = c("Protein","Disease"), all.x = TRUE)
best[, beta := fifelse(eclass == "cis", beta_PDcis, beta_PDtrans)]
best[, sign := sign(beta)]

out <- best[, .(protein, disease, mr_tier = mr_tier_final, tier_rank = trank,
                edge_class = eclass, coloc_confirmed = colocok, coloc_pph4,
                n_arms_qualified, replicated, beta = round(beta, 4), sign)]
setorder(out, -tier_rank, -coloc_confirmed)
f <- file.path(heap_project_output("support","intervention_compare"), "mr_pd_tiered.tsv")
fwrite(out, f, sep = "\t")
message(sprintf("Wrote %s\n  %d protein-disease edges | Tier1+:%d Tier1:%d Tier2:%d Suggestive:%d | coloc-confirmed:%d",
                f, nrow(out), sum(out$mr_tier=="Tier1plus"), sum(out$mr_tier=="Tier1"),
                sum(out$mr_tier=="Tier2"), sum(out$mr_tier=="Suggestive"), sum(out$coloc_confirmed)))
cat("\n=== T2D edges for the Fig6 cast ===\n")
cast <- c("ICAM1","FABP4","FURIN","ADIPOQ","ADH1B","SULT2A1")
print(out[disease == "finngen_R12_T2D" & protein %in% cast][order(-tier_rank),
          .(protein, mr_tier, edge_class, coloc_confirmed, pph4 = round(coloc_pph4,2), n_arms_qualified, beta)])
