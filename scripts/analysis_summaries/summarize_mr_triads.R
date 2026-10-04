#!/usr/bin/env Rscript
# ============================================================================
# summarize_mr_triads.R
# ----------------------------------------------------------------------------
# Supplementary tables for the Module-5 MR triads (exposure -> protein ->
# disease) under the CANONICAL Tier-1 motif definition.
#
# CRITICAL -- two motif definitions exist and they are NOT nested:
#   * MRmotifs.tsv carries precomputed motif_* flags defined on nominal
#     significance (padj). Under that rule the mediator motif has 84 triads.
#   * The manuscript (main Fig 4b, and fig_mr_motif_overview.R) recomputes every
#     motif with all six edges required at Tier 1 / Tier1plus. Under that rule
#     the mediator motif has 7 triads / 3 proteins.
# The main text quotes the Tier-1 numbers, so this script applies the Tier-1 rule
# and IGNORES the precomputed motif_* columns. Using them would contradict the
# paper by an order of magnitude. FURIN is the canonical example of the
# non-nesting: a Tier-1 mediator that the significance rule disqualifies.
# Motif rules are copied verbatim from fig_mr_motif_overview.R:80-85.
#
# Outputs (manuscript_stats/module5/):
#   mr_triad_motifs.tsv    one row per motif-carrying triad, with the per-edge
#                          effect estimates for the edges the motif asserts
#   mr_motif_counts.tsv    triad and protein counts per motif, under BOTH rules,
#                          so the difference is documented rather than hidden
#
#   ARM SCOPE (settled 2026-08-29). An edge that touches the protein is measured on
#   one pQTL platform, so it is evaluated WITHIN that platform: UKB Olink or deCODE
#   SomaScan. A motif is therefore assembled within a platform -- its four
#   protein-involving edges must come from the same panel -- and is Tier 1 if EITHER
#   platform supports it, Tier 1+ if BOTH do. The exposure->disease and
#   disease->exposure legs carry no protein and are common to the two.
#
#   Evaluating instead on a pooled edge set would let a deCODE reverse edge veto a
#   UKB triad, mixing platforms inside one motif; that drops FURIN from the mediator
#   set and is NOT the rule.
#
#   Run:  Rscript scripts/analysis_summaries/summarize_mr_triads.R
#         (both arms are evaluated; no argument)
# ============================================================================
suppressPackageStartupMessages(library(data.table))
options(scipen = 999)

local({
  cand <- c(Sys.getenv("HEAP_PATHS_FILE", unset = ""),
            "/n/groups/patel/shakson_ukb/HEAP/workflow/00_paths.R")
  hit <- cand[nzchar(cand) & file.exists(cand)][1]; if (!is.na(hit)) source(hit)
})
SUMD <- file.path(heap_project_output("mr_edges"), "summary")
OUT  <- file.path(heap_root, "docs", "manuscript_stats", "module5")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

ARMS <- c("UKB", "DECODE")
# the deCODE triad table lives in its own subdirectory, not as a `dataset` partition
# of the top-level file -- reading MRmotifs.tsv and filtering dataset=="DECODE"
# silently returns zero rows.
motif_src <- function(arm) {
  if (arm == "UKB") file.path(SUMD, "MRmotifs.tsv")
  else              file.path(SUMD, "DECODE", "MRmotifs.tsv")
}
TE_ALL <- fread(file.path(SUMD, "mr_tiered_edges.tsv"))

## ---- evaluate the motif signature WITHIN each arm --------------------------
NAME <- c(A = "A Mediator (E->P->D)", B = "B Biomarker", C = "C Exposure-marker",
          D = "D Reverse (P->E)",     E = "E Disease-liability (D->P)")

eval_arm <- function(arm) {
  m  <- fread(motif_src(arm))
  T1 <- TE_ALL[dataset == arm & mr_tier_final %in% c("Tier1", "Tier1plus")]
  k  <- function(d) T1[edge_dir %in% d, paste(src_id, tgt_id)]
  EP <- k("E_to_P"); PD <- k(c("Pcis_to_D", "Ptrans_to_D")); ED <- k("E_to_D")
  PE <- k(c("Pcis_to_E", "Ptrans_to_E")); DP <- k("D_to_P"); DE <- k("D_to_E")
  t <- m[, .(Exposure, Protein, Disease, ExposureCategory, ICD10, Disease_UKB)]
  t[, `:=`(pEP = paste(Exposure, Protein) %in% EP,
           pPD = paste(Protein,  Disease)  %in% PD,
           pED = paste(Exposure, Disease)  %in% ED,
           pPE = paste(Protein,  Exposure) %in% PE,
           pDP = paste(Disease,  Protein)  %in% DP,
           pDE = paste(Disease,  Exposure) %in% DE)]
  t[, A := pEP &  pPD & pED & !pPE & !pDP & !pDE]
  t[, B := pEP & !pPD & pED & !pPE &  pDP & !pDE]
  t[, C := pEP & !pPD & pED & !pPE & !pDP]
  t[, D := !pEP & pPE]
  t[, E := pDP & pDE]
  message(sprintf("  %-7s %s triads | Tier-1 edges %s | mediator triads %d",
                  arm, format(nrow(m), big.mark = ","),
                  format(nrow(T1), big.mark = ","), sum(t$A)))
  list(flags = t, raw = m)
}
message("evaluating each pQTL platform separately:")
AR <- setNames(lapply(ARMS, eval_arm), ARMS)

# a triad carries a motif at Tier 1 if EITHER platform supports it; Tier 1+ if BOTH.
key3 <- function(d) paste(d$Exposure, d$Protein, d$Disease)
sets <- lapply(names(NAME), function(g)
  setNames(lapply(ARMS, function(a) key3(AR[[a]]$flags[get(g) == TRUE])), ARMS))
names(sets) <- names(NAME)

## ---- counts: per arm, union (Tier 1) and both-arm (Tier 1+) ----------------
# `nominal_*` keeps the precomputed motif_* flags for contrast; those are defined on
# NOMINAL significance in the UKB table and are NOT nested with the Tier-1 rule.
FLAG <- c(A = "motif_A_mediator", B = "motif_B_biomarker", C = "motif_C_exposure_marker",
          D = "motif_D_P_to_E",   E = "motif_E_disease_liability")
tru  <- function(x) x %in% c(TRUE, "TRUE", 1, "1")
mUKB <- AR[["UKB"]]$raw
prot_of <- function(keys) uniqueN(sub("^\\S+ (\\S+) .*$", "\\1", keys))

cnt <- rbindlist(lapply(names(NAME), function(g) {
  u <- sets[[g]][["UKB"]]; d <- sets[[g]][["DECODE"]]
  un <- union(u, d); bo <- intersect(u, d)
  data.table(
    motif            = NAME[[g]],
    ukb_triads       = length(u),
    decode_triads    = length(d),
    tier1_triads     = length(un),          # EITHER platform
    tier1_proteins   = prot_of(un),
    tier1plus_triads = length(bo),          # BOTH platforms
    tier1plus_proteins = prot_of(bo),
    nominal_triads   = sum(tru(mUKB[[FLAG[[g]]]])),
    nominal_proteins = uniqueN(mUKB$Protein[tru(mUKB[[FLAG[[g]]]])]))
}))
fwrite(cnt, file.path(OUT, "mr_motif_counts.tsv"), sep = "\t")
print(cnt)

## ---- one row per motif-carrying triad --------------------------------------
# Union of the two platforms. `arms` records which supported the motif, so a reader
# can see at a glance that the mediators are UKB-only while the disease-liability
# set draws on both. Effect estimates come from UKB where it supports the triad
# (the primary panel), otherwise from deCODE.
label_motifs <- function(t) {
  t[, motif := {
    lab <- character(.N)
    for (g in names(NAME)) lab <- ifelse(get(g), ifelse(nzchar(lab),
                                          paste(lab, NAME[[g]], sep = "; "), NAME[[g]]), lab)
    lab }]
  t[nzchar(motif)]
}
triU <- label_motifs(copy(AR[["UKB"]]$flags))
triD <- label_motifs(copy(AR[["DECODE"]]$flags))
triU[, arms := "UKB"]; triD[, arms := "DECODE"]
kU <- key3(triU); kD <- key3(triD)
triD <- triD[!(kD %in% kU)]                       # UKB wins on overlap
triU[key3(triU) %in% kD, arms := "UKB+DECODE"]
tri  <- rbind(triU, triD, fill = TRUE)

num <- function(x) suppressWarnings(as.numeric(x))
add <- function(dst, src, pfx) {
  for (s in c("beta", "se", "padj")) {
    cc <- paste0(s, "_", pfx)
    if (cc %in% names(src)) set(dst, j = cc, value = signif(num(src[[cc]]), 4))
  }
  dst
}
# pull each triad's estimates from the arm that actually supported it
for (a in ARMS) {
  ra  <- AR[[a]]$raw
  idx <- which(tri$arms == a | (a == "UKB" & tri$arms == "UKB+DECODE"))
  if (!length(idx)) next
  i <- match(paste(tri$Exposure[idx], tri$Protein[idx], tri$Disease[idx]), key3(ra))
  for (pf in c("EP", "PDcis", "PDtrans", "ED", "PEcis", "PEtrans", "DP", "DE"))
    for (st in c("beta", "se", "padj")) {
      cc <- paste0(st, "_", pf)
      if (cc %in% names(ra)) {
        if (!cc %in% names(tri)) set(tri, j = cc, value = NA_real_)
        set(tri, i = idx, j = cc, value = signif(num(ra[[cc]][i]), 4))
      }
    }
}

setnames(tri, c("ExposureCategory", "ICD10"), c("Exposure category", "ICD10"))
tri[, c("A", "B", "C", "D", "E", "pEP", "pPD", "pED", "pPE", "pDP", "pDE") := NULL]
setcolorder(tri, c("Exposure", "Exposure category", "Protein", "Disease", "ICD10",
                   "Disease_UKB", "motif", "arms"))
setorder(tri, motif, Exposure, Protein, Disease)
f <- file.path(OUT, "mr_triad_motifs.tsv")
fwrite(tri, f, sep = "\t")
message(sprintf("\n  %s  %s triads x %d cols  (%s distinct proteins)",
                basename(f), format(nrow(tri), big.mark = ","), ncol(tri),
                format(uniqueN(tri$Protein), big.mark = ",")))
message(sprintf("  arms: %s", paste(sprintf("%s=%s", names(table(tri$arms)),
                                            format(as.integer(table(tri$arms)), big.mark=",")), collapse=" | ")))
