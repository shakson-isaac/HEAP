#!/usr/bin/env Rscript
# ============================================================================
# build_mr_protein_class.R  — unified MR classification for Fig5 b & c
# NOTE ON WORDING: this script keeps Tier1 OR Tier1plus, i.e. "Tier 1 or above".
# Do NOT call that set Tier 1+: in the manuscript Tier 1+ means replicated across
# BOTH pQTL platforms, and only ASGR1 and PCSK9 qualify. The other causal
# proteins here (FURIN, ICAM1, ALCAM, PRSS8, ADM, SOST) are Tier 1 only.
# ----------------------------------------------------------------------------
# ONE MR source for both panels (replaces the flat annotate_mr coloring in c and
# the ad-hoc gold count in b). From mr_tiered_edges, keep only **Tier1 or above**
# edges, restricted to a UNIFIED cardiometabolic disease set, and split by
# direction: protein→disease = CAUSAL, disease→protein = REPORTER. Counts are
# kept PER DISEASE (not lumped).
#
# Outputs (support/intervention_compare/):
#   mr_protein_class_t2d.tsv   per protein, T2D-SPECIFIC class (causal/reporter)
#                              + edge_class/coloc/replicated  -> PANEL c
#   exposure_mr_validation.tsv per exposure, # Tier-1-or-above CAUSAL proteins
#                              PER DISEASE (nT2D/nObesity/nLipids/nHTN) + total
#                              + which proteins -> PANEL b
# Run: module load gcc/14.2.0 R/4.4.2; HEAP_PATHS_FILE=.../00_paths.R Rscript ...
# ============================================================================
local({
  cand <- c(Sys.getenv("HEAP_PATHS_FILE", unset=""), "/n/groups/patel/shakson_ukb/HEAP/workflow/00_paths.R")
  source(cand[nzchar(cand) & file.exists(cand)][1])
})
suppressPackageStartupMessages(library(data.table))
icd <- heap_project_output("support","intervention_compare")

# unified cardiometabolic disease set (the GLP1/exercise target space)
CM <- c(finngen_R12_T2D="T2D", finngen_R12_T2D_WIDE="T2D",
        finngen_R12_E4_OBESITY="Obesity", finngen_R12_E4_OBESITYNAS="Obesity", finngen_R12_E4_OBESITYCAL="Obesity",
        finngen_R12_E4_LIPOPROT="Lipids", finngen_R12_I9_HYPTENSESS="Hypertension")
T2DSET <- c("finngen_R12_T2D","finngen_R12_T2D_WIDE")
DZLEV  <- c("T2D","Obesity","Lipids","Hypertension")

te <- fread(heap_project_output("mr_edges","summary","mr_tiered_edges.tsv"))
t1 <- te[mr_tier_final %in% c("Tier1","Tier1plus")]              # Tier1 or above ONLY
repl <- function(x) { x <- tolower(as.character(x)); x %in% c("true","t","1","yes","both") }

# ---- CAUSAL (protein->disease) Tier 1 or above, cardiometabolic -----------
caus <- t1[edge_dir %in% c("Pcis_to_D","Ptrans_to_D") & tgt_id %in% names(CM)]
caus[, `:=`(protein=src_id, dz=CM[tgt_id], edge_class=fifelse(edge_dir=="Pcis_to_D","cis","trans"),
            coloc=tolower(coloc_status)=="confirmed", replicated=repl(replicated))]
caus_best <- caus[order(-coloc, edge_class!="cis", -replicated)][, .SD[1], by=.(protein, dz)]

# ---- REPORTER (disease->protein) Tier 1 or above, cardiometabolic ---------
rep <- t1[edge_dir=="D_to_P" & src_id %in% names(CM)]
rep[, `:=`(protein=tgt_id, dz=CM[src_id], replicated=repl(replicated))]
rep_best <- rep[, .(replicated=any(replicated)), by=.(protein, dz)]

# ---- (1) T2D-specific per-protein class (panel c) -------------------------
caus_t2d <- caus[tgt_id %in% T2DSET, .(protein, edge_class, coloc, replicated)][order(-coloc, edge_class!="cis")][, .SD[1], by=protein]
rep_t2d  <- rep[src_id %in% T2DSET, .(replicated=any(replicated)), by=protein]
cls <- rbind(
  caus_t2d[, .(protein, t2d_class="causal", edge_class, coloc, replicated)],
  rep_t2d[!protein %in% caus_t2d$protein, .(protein, t2d_class="reporter", edge_class=NA_character_, coloc=FALSE, replicated)]
)
fwrite(cls, file.path(icd, "mr_protein_class_t2d.tsv"), sep="\t")

# ---- (2) per-exposure Tier-1-or-above CAUSAL counts PER DISEASE (panel b) -
sc <- fread(file.path(icd,"intervention_scatter_mr.tsv"))
expp <- unique(sc[, .(exposure_id, Eid, Category, protein)])
caus_map <- caus_best[, .(protein, dz)]                          # protein -> causal disease(s)
ec <- merge(expp, caus_map, by="protein", allow.cartesian=TRUE)  # exposure x its causal proteins
perdz <- dcast(ec[, .(n=uniqueN(protein)), by=.(exposure_id, dz)], exposure_id ~ dz, value.var="n", fill=0)
for (d in setdiff(DZLEV, names(perdz))) perdz[[d]] <- 0L
prot <- ec[, .(causal_proteins=paste(sort(unique(paste0(protein,"(",dz,")"))), collapse=";"),
               n_causal=uniqueN(protein)), by=exposure_id]
val <- merge(unique(expp[, .(exposure_id, Eid, Category)]), perdz, by="exposure_id", all.x=TRUE)
val <- merge(val, prot, by="exposure_id", all.x=TRUE)
for (d in DZLEV) val[is.na(get(d)), (d):=0L]
val[is.na(n_causal), n_causal:=0L]
val[, mr_anchored := n_causal >= 1L]
setcolorder(val, c("exposure_id","Eid","Category", DZLEV, "n_causal","mr_anchored","causal_proteins"))
setorder(val, -n_causal)
fwrite(val, file.path(icd, "exposure_mr_validation.tsv"), sep="\t")

message(sprintf("mr_protein_class_t2d.tsv: %d causal + %d reporter (Tier 1 or above, T2D)",
                cls[t2d_class=="causal",.N], cls[t2d_class=="reporter",.N]))
message(sprintf("CAUSAL proteins per disease (Tier 1 or above): %s",
                paste(sprintf("%s=%s", DZLEV, sapply(DZLEV, function(d) caus_best[dz==d, uniqueN(protein)])), collapse="  ")))
message(sprintf("exposure_mr_validation.tsv: %d exposures, %d MR-anchored (>=1 Tier-1-or-above causal protein)", nrow(val), sum(val$mr_anchored)))
cat("\n=== causal proteins (Tier 1 or above, cardiometabolic) ===\n"); print(caus_best[, .(protein, dz, edge_class, coloc)][order(dz)])