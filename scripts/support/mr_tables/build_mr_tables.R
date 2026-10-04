#!/usr/bin/env Rscript
# ============================================================================
# support/mr_tables/build_mr_tables.R   <cohort>      cohort in {UKB, DECODE}
# ----------------------------------------------------------------------------
# Assemble the canonical Module-5 MR summary tables for one instrument arm by
# walking the FRESH per-edge output tree directly:
#
#   mr_edges/MR_UKB_primary/<edge_dir>/<src>/<tgt>/<edge_dir>_<suffix>.tsv     (UKB)
#   mr_edges_decode/MR_deCODE_replication/<edge_dir>/<src>/<tgt>/...           (DECODE)
#
# This REPLACES the dead chain MRviz.R -> summary/<EDGE>res.csv -> MRmotifs.csv,
# which read a stale Jan-2025 aggregate. The per-edge tree (8 directed edges) is
# the source of truth; this script aggregates the primary IVW/Wald estimate plus
# the heterogeneity (IVW Cochran Q), pleiotropy (Egger intercept) and Steiger
# directionality diagnostics into:
#
#   <out>/mr_sensitivity_long.tsv     one row per edge: primary estimate + BH q +
#                                     het/pleio/steiger flags + hit_after_sens
#   <out>/sensitivity_by_edgedir.tsv  per-edge_dir hit counts + retention %
#   <out>/sensitivity_overall.tsv     per-dataset roll-up
#   <out>/MRmotifs.tsv                wide one-row-per-(E,P,D) triad table with
#                                     beta_/se_/padj_ for all 8 edges + motif
#                                     A-E classification (drop-in for annotate_mr.R
#                                     and the shared-vs-unique triad figure)
#
# out = mr_edges/summary           (UKB)
#       mr_edges/summary/DECODE    (DECODE)
#
# Motif logic ported verbatim from ModuleMR/MRviz.R; sensitivity logic from
# ModuleMR/COMBO/MRviolation_tagging.R. Missing diagnostic = pass (same rule).
#
# Usage:
#   module load gcc/14.2.0 R/4.4.2
#   HEAP_PATHS_FILE=.../workflow/00_paths.R CPUS=16 \
#     Rscript scripts/support/mr_tables/build_mr_tables.R UKB
# Env:
#   CPUS               file-read parallelism (default detectCores()-1)
#   HEAP_MR_TABLES_MAXEDGES   cap edges per edge_dir (smoke test; default 0 = all)
# ============================================================================

local({
  cand <- c(Sys.getenv("HEAP_PATHS_FILE", unset = ""),
            file.path(getwd(), "workflow", "00_paths.R"),
            "/n/groups/patel/shakson_ukb/HEAP/workflow/00_paths.R")
  hit <- cand[nzchar(cand) & file.exists(cand)][1]
  if (is.na(hit)) stop("Could not locate workflow/00_paths.R (set HEAP_PATHS_FILE).")
  source(hit)
})
suppressPackageStartupMessages({ library(data.table); library(parallel) })

`%||%` <- function(a, b) if (!is.null(a) && length(a) && !is.na(a)) a else b
ts <- function(...) message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...)))

# ---- args / config ---------------------------------------------------------
args   <- commandArgs(trailingOnly = TRUE)
cohort <- toupper(args[1] %||% "UKB")

ALPHA_Q    <- 0.05   # BH q hit threshold (within edge_dir)
ALPHA_SENS <- 0.05   # heterogeneity Q_pval / Egger intercept p threshold
ALPHA_DIR  <- 0.05   # Steiger directionality p; direction must be ESTABLISHED, not merely signed
MAXEDGES   <- as.integer(Sys.getenv("HEAP_MR_TABLES_MAXEDGES", "0"))
nthreads   <- suppressWarnings(as.integer(Sys.getenv("CPUS", "")))
if (is.na(nthreads) || nthreads < 1) nthreads <- max(1L, detectCores() - 1L)
setDTthreads(1L)  # we parallelise at the file level

COHORTS <- list(
  UKB    = list(tag = "UKB",
                per_edge = heap_project_output("mr_edges", "MR_UKB_primary"),
                out      = heap_project_output("mr_edges", "summary")),
  DECODE = list(tag = "DECODE",
                per_edge = heap_project_output("mr_edges_decode", "MR_deCODE_replication"),
                out      = heap_project_output("mr_edges", "summary", "DECODE"))
)
if (!cohort %in% names(COHORTS)) stop("cohort must be UKB or DECODE; got ", cohort)
CFG       <- COHORTS[[cohort]]
out_override <- Sys.getenv("HEAP_MR_TABLES_OUT", "")   # smoke-test / alt-location override
if (nzchar(out_override)) CFG$out <- out_override
edges_dir <- heap_project_output("mr_edges", "global_edges")
dir.create(CFG$out, recursive = TRUE, showWarnings = FALSE)

ts("cohort=", cohort, " per_edge=", CFG$per_edge)
ts("nthreads=", nthreads, "  out=", CFG$out)
if (!dir.exists(CFG$per_edge))
  stop("per-edge tree not found: ", CFG$per_edge)

# edge_dir -> (edge list file, src column, tgt column, wide-table category) -----
EDGE_SPEC <- list(
  E_to_P      = list(list_file = "edges_EP.tsv", src = "Exposure", tgt = "Protein",  cat = "EP"),
  Pcis_to_E   = list(list_file = "edges_PE.tsv", src = "Protein",  tgt = "Exposure", cat = "PEcis"),
  Ptrans_to_E = list(list_file = "edges_PE.tsv", src = "Protein",  tgt = "Exposure", cat = "PEtrans"),
  Pcis_to_D   = list(list_file = "edges_PD.tsv", src = "Protein",  tgt = "Disease",  cat = "PDcis"),
  Ptrans_to_D = list(list_file = "edges_PD.tsv", src = "Protein",  tgt = "Disease",  cat = "PDtrans"),
  E_to_D      = list(list_file = "edges_ED.tsv", src = "Exposure", tgt = "Disease",  cat = "ED"),
  D_to_E      = list(list_file = "edges_DE.tsv", src = "Disease",  tgt = "Exposure", cat = "DE"),
  D_to_P      = list(list_file = "edges_DP.tsv", src = "Disease",  tgt = "Protein",  cat = "DP")
)

# ---- build the (edge_dir, src, tgt) universe from the edge lists ------------
meta <- rbindlist(lapply(names(EDGE_SPEC), function(ed) {
  sp <- EDGE_SPEC[[ed]]
  el <- fread(file.path(edges_dir, sp$list_file), showProgress = FALSE)
  dt <- unique(data.table(edge_dir = ed, cat = sp$cat,
                          src = el[[sp$src]], tgt = el[[sp$tgt]]))
  if (MAXEDGES > 0 && nrow(dt) > MAXEDGES) dt <- dt[1:MAXEDGES]
  dt
}), fill = TRUE)
meta[, .row := .I]
ts("edge universe: ", nrow(meta), " directed edges across ", uniqueN(meta$edge_dir), " edge_dirs")

# ---- file readers ----------------------------------------------------------
safe_fread <- function(p) {
  if (!file.exists(p)) return(NULL)
  s <- file.info(p)$size
  if (is.na(s) || s == 0) return(NULL)
  tryCatch(fread(p, sep = "\t", showProgress = FALSE), error = function(e) NULL)
}

# pull one tidy row of values from each per-edge file (NULL => drop)
pick_summary <- function(d) {                       # <edge>_summary.tsv
  if (is.null(d) || !nrow(d) || !"method" %in% names(d)) return(NULL)
  k <- d[method %in% c("Inverse variance weighted", "Wald ratio")]
  if (!nrow(k)) k <- d
  k <- k[1L]
  data.table(method = k$method, nsnp = as.integer(k$nsnp),
             b = as.numeric(k$b), se = as.numeric(k$se), pval = as.numeric(k$pval))
}
pick_het <- function(d) {                           # <edge>_heterogeneity.tsv (IVW Cochran Q)
  if (is.null(d) || !nrow(d) || !"Q_pval" %in% names(d)) return(NULL)
  k <- if ("method" %in% names(d)) d[method == "Inverse variance weighted"] else d
  if (!nrow(k)) return(NULL)
  k <- k[1L]
  data.table(het_Q = as.numeric(k$Q), het_Qdf = as.numeric(k$Q_df), het_pval = as.numeric(k$Q_pval))
}
pick_ple <- function(d) {                           # <edge>_pleiotropy.tsv (Egger intercept)
  if (is.null(d) || !nrow(d) || !"pval" %in% names(d)) return(NULL)
  d <- d[1L]
  data.table(egger_intercept = as.numeric(d$egger_intercept), egger_pval = as.numeric(d$pval))
}
pick_stg <- function(d) {                           # <edge>_steiger.tsv (directionality)
  if (is.null(d) || !nrow(d) || !"correct_causal_direction" %in% names(d)) return(NULL)
  d <- d[1L]
  data.table(steiger_dir = as.logical(d$correct_causal_direction),
             steiger_pval = as.numeric(d$steiger_pval))
}
pick_presso <- function(d) {                        # <edge>_presso.tsv (outlier correction)
  if (is.null(d) || !nrow(d) || !"presso_global_pval" %in% names(d)) return(NULL)
  gp <- suppressWarnings(as.numeric(sub("^<", "", as.character(d$presso_global_pval[1]))))
  oc <- if ("MR Analysis" %in% names(d)) d[`MR Analysis` == "Outlier-corrected"] else d[0]
  data.table(presso_global_p = gp,
             presso_corr_b = if (nrow(oc)) suppressWarnings(as.numeric(oc$`Causal Estimate`[1])) else NA_real_,
             presso_corr_p = if (nrow(oc)) suppressWarnings(as.numeric(oc$`P-value`[1])) else NA_real_)
}
pick_methods <- function(d) {                       # <edge>_mr_methods.tsv (robust estimators)
  if (is.null(d) || !nrow(d) || !"method" %in% names(d)) return(NULL)
  wm <- d[method == "Weighted median"]; eg <- d[method == "MR Egger"]
  data.table(
    wm_b = if (nrow(wm)) as.numeric(wm$b[1]) else NA_real_,
    wm_p = if (nrow(wm)) as.numeric(wm$pval[1]) else NA_real_,
    egslope_b = if (nrow(eg)) as.numeric(eg$b[1]) else NA_real_,
    egslope_p = if (nrow(eg)) as.numeric(eg$pval[1]) else NA_real_)
}

gather_suffix <- function(meta, per_edge, suffix, picker, nthreads) {
  files <- file.path(per_edge, meta$edge_dir, meta$src, meta$tgt,
                     paste0(meta$edge_dir, "_", suffix, ".tsv"))
  ix <- which(file.exists(files))
  ts("  ", suffix, ": ", length(ix), "/", length(files), " files present")
  if (!length(ix)) return(data.table(.row = integer(0)))
  parts <- mclapply(ix, function(i) {
    v <- picker(safe_fread(files[i]))
    if (is.null(v) || !nrow(v)) return(NULL)
    v[, .row := meta$.row[i]]
    v
  }, mc.cores = nthreads)
  rbindlist(Filter(Negate(is.null), parts), fill = TRUE)
}

# ---- gather all four suffixes ---------------------------------------------
ts("Reading per-edge files (", nthreads, " cores) ...")
prim <- gather_suffix(meta, CFG$per_edge, "summary",       pick_summary, nthreads)
het  <- gather_suffix(meta, CFG$per_edge, "heterogeneity", pick_het,     nthreads)
ple  <- gather_suffix(meta, CFG$per_edge, "pleiotropy",    pick_ple,     nthreads)
stg  <- gather_suffix(meta, CFG$per_edge, "steiger",       pick_stg,     nthreads)

if (!nrow(prim)) stop("No primary MR estimates found under ", CFG$per_edge,
                      " — has the MR run produced output yet?")

# ---- assemble the long per-edge sensitivity table --------------------------
L <- merge(meta, prim, by = ".row", all = FALSE)            # keep edges with a primary estimate
for (d in list(het, ple, stg)) if (nrow(d)) L <- merge(L, d, by = ".row", all.x = TRUE)
# ensure all diagnostic cols exist even if a whole suffix was empty
for (cc in c("het_Q","het_Qdf","het_pval","egger_intercept","egger_pval","steiger_pval"))
  if (!cc %in% names(L)) L[, (cc) := NA_real_]
if (!"steiger_dir" %in% names(L)) L[, steiger_dir := NA]

safe_padj <- function(p) if (!length(p) || all(is.na(p))) rep(NA_real_, length(p)) else p.adjust(p, "BH")
L[, pval_adj := safe_padj(pval), by = edge_dir]            # BH within edge_dir (== MRviz)
L[, mr_hit := !is.na(pval_adj) & pval_adj < ALPHA_Q]
L[, het_flag    := !is.na(het_pval)   & het_pval   < ALPHA_SENS]
L[, pleio_flag  := !is.na(egger_pval) & egger_pval < ALPHA_SENS]
L[, steiger_flag:= !is.na(steiger_dir) & steiger_dir == FALSE]    # wrong causal direction
L[, sens_pass := (is.na(het_pval)   | het_pval   >= ALPHA_SENS) &
                 (is.na(egger_pval) | egger_pval >= ALPHA_SENS)]
L[, hit_after_sens := mr_hit & sens_pass]

# ============================================================================
# CONFIDENCE TIERING (see mr_schematics/mr_hit_flowchart.tex)
#  trunk: significant -> instrument sufficiency (>=3 SNP, or cis low-SNP) ->
#         robustness (Q/Egger clean, or rescued by MR-PRESSO/weighted median) ->
#         Steiger direction -> qualified; then class-specific tail.
#  Tier 1+ (cross-arm replication) is added later by compare_arms.R.
# ============================================================================
EDGE_CLASS <- c(Pcis_to_D = "cis", Pcis_to_E = "cis", Ptrans_to_D = "trans",
                Ptrans_to_E = "trans", E_to_P = "polygenic", D_to_P = "polygenic",
                E_to_D = "polygenic", D_to_E = "polygenic")
EDGE_LANE  <- c(Pcis_to_D = "pQTL", Pcis_to_E = "pQTL", Ptrans_to_D = "pQTL",
                Ptrans_to_E = "pQTL", E_to_P = "protein_outcome", D_to_P = "protein_outcome",
                E_to_D = "protein_free", D_to_E = "protein_free")
L[, edge_class := EDGE_CLASS[edge_dir]]
L[, lane := EDGE_LANE[edge_dir]]
L[, het_clean   := is.na(het_pval)   | het_pval   >= ALPHA_SENS]
L[, pleio_clean := is.na(egger_pval) | egger_pval >= ALPHA_SENS]
L[, clean := het_clean & pleio_clean]

# pass 2: read MR-PRESSO + robust-estimator files ONLY for sig, >=3 SNP, not-clean
# edges (the rescue candidates) — keeps IO small.
resc_rows <- L[mr_hit == TRUE & nsnp >= 3 & clean == FALSE, .row]
ts("rescue candidates (sig, >=3 SNP, not clean): ", length(resc_rows))
if (length(resc_rows)) {
  meta_resc <- meta[.row %in% resc_rows]
  pr <- gather_suffix(meta_resc, CFG$per_edge, "presso",     pick_presso, nthreads)
  mm <- gather_suffix(meta_resc, CFG$per_edge, "mr_methods", pick_methods, nthreads)
  if (nrow(pr)) L <- merge(L, pr, by = ".row", all.x = TRUE)
  if (nrow(mm)) L <- merge(L, mm, by = ".row", all.x = TRUE)
}
for (cc in c("presso_global_p","presso_corr_b","presso_corr_p","wm_b","wm_p","egslope_b","egslope_p"))
  if (!cc %in% names(L)) L[, (cc) := NA_real_]

L[, rescued_presso := !is.na(presso_global_p) & presso_global_p < ALPHA_SENS &
                      is.finite(presso_corr_p) & presso_corr_p < ALPHA_Q &
                      is.finite(presso_corr_b) & sign(presso_corr_b) == sign(b)]
L[, rescued_median := is.finite(wm_p) & wm_p < ALPHA_Q & is.finite(wm_b) & sign(wm_b) == sign(b)]
L[is.na(rescued_presso), rescued_presso := FALSE]
L[is.na(rescued_median), rescued_median := FALSE]
L[, robust_pass := clean | rescued_presso | rescued_median]
L[, het_status := fcase(
    nsnp < 3 & edge_class == "cis", "na_cis_lowsnp",
    nsnp < 3,                       "na_lowsnp",
    clean,                          "homogeneous",
    rescued_presso,                 "rescued_presso",
    rescued_median,                 "rescued_median",
    default =                       "heterogeneous")]
# ---- causal direction: ESTABLISHED, not merely signed ----------------------
# directionality_test() returns two things from the SAME pair of numbers:
#   steiger_dir  == (snp_r2.exposure > snp_r2.outcome)   -- the SIGN of the gap
#   steiger_pval == whether those two R^2 differ at all  -- the CONFIDENCE
# Reading the sign alone treats a 1.05x gap (p~0.9, indistinguishable from noise)
# exactly like a 10x gap (p~1e-20). Because a weakly instrumented exposure has a
# small snp_r2.exposure by construction, that penalised exposures in proportion to
# how poorly they were measured: the weakest-instrumented third of exposures lost
# 48% of eligible E->P edges to "Reverse", the strongest third 1% (Spearman
# rho = -0.82). It also cut both ways -- 199 of 219 E->P demotions rested on a
# non-significant test, while 440 Tier-1 edges were credited with a "correct
# causal direction" that was equally unestablished.
# So: require the test to be significant before acting on its sign, in EITHER
# direction. Non-significant -> direction unresolved -> Suggestive (see below).
# steiger_dir missing (3 edges) keeps the previous pass behaviour.
L[, direction_established := !is.na(steiger_pval) & steiger_pval < ALPHA_DIR]
L[, direction_unresolved  := !is.na(steiger_dir) & !direction_established]
L[, steiger_ok := is.na(steiger_dir) | (steiger_dir == TRUE & direction_established)]

L[, mr_tier := fcase(
    mr_hit == FALSE,                                  "Null",
    nsnp < 3 & edge_class != "cis",                   "Suggestive",   # insufficient instruments
    nsnp >= 3 & robust_pass == FALSE,                 "Suggestive",   # heterogeneous/pleiotropic
    direction_unresolved,                             "Suggestive",   # Steiger not significant
    steiger_ok == FALSE,                              "Reverse",      # established, and backwards
    lane == "pQTL" & edge_class == "trans",           "Tier2",        # trans-only
    lane == "pQTL" & edge_class == "cis",             "Tier1",        # cis (coloc-pending)
    lane == "protein_outcome",                        "Tier1",
    lane == "protein_free",                           "Tier1",        # single-source ceiling
    default =                                         "Null")]
L[, tier_reason := fcase(
    mr_tier == "Null",                                            "not_significant",
    mr_tier == "Suggestive" & nsnp < 3 & edge_class != "cis",     "insufficient_instruments",
    mr_tier == "Suggestive" & nsnp >= 3 & robust_pass == FALSE,   "heterogeneous_pleiotropic",
    mr_tier == "Suggestive",                                      "direction_unresolved",
    mr_tier == "Reverse",                                         "reverse_direction",
    mr_tier == "Tier2",                                           "trans_only",
    mr_tier == "Tier1" & lane == "pQTL",                          "cis_robust_directional",
    mr_tier == "Tier1" & lane == "protein_outcome",               "robust_directional",
    mr_tier == "Tier1" & lane == "protein_free",                  "robust_directional_singlesource",
    default =                                                     "")]
# --- colocalization HARD GATE (PP.H4 >= 0.8) for cis-pQTL edges --------------
# Systematic coloc.abf results from support/coloc/run_coloc_systematic.R. A
# colocalized cis edge keeps Tier1 (-> Tier1+ if it also replicates cross-arm);
# a cis edge whose pQTL & disease/exposure signals do NOT share a causal variant
# (PP.H4 < 0.8 -> LD confounding) is DOWNGRADED to Tier2. cis edges with no coloc
# result yet keep Tier1 but are tagged coloc_unavailable (not penalised for an
# un-run test). Non-cis / non-pQTL edges are coloc-N/A.
L[, PP_H4 := NA_real_]
coloc_fp <- heap_project_output("support", "coloc", "coloc_results.tsv")
if (file.exists(coloc_fp)) {
  CO <- fread(coloc_fp)
  if (all(c("arm","protID","target","PP.H4") %in% names(CO))) {
    CO <- CO[arm == CFG$tag & is.finite(PP.H4), .(PP_H4 = max(PP.H4)), by = .(protID, target)]
    if (nrow(CO)) L[CO, on = c(src = "protID", tgt = "target"), PP_H4 := i.PP_H4]
  }
}
L[, coloc_status := fcase(
    !(lane == "pQTL" & edge_class == "cis"),  "na",
    !is.finite(PP_H4),                        "coloc_unavailable",
    PP_H4 >= 0.8,                             "colocalized",
    default =                                 "not_colocalized")]
# Hard gate: cis-pQTL Tier1 edge that fails coloc -> Tier2.
L[mr_tier == "Tier1" & lane == "pQTL" & edge_class == "cis" & coloc_status == "not_colocalized",
  `:=`(mr_tier = "Tier2", tier_reason = "cis_not_colocalized")]
L[mr_tier == "Tier1" & coloc_status == "colocalized", tier_reason := "cis_colocalized"]

L[, dataset := CFG$tag]
setnames(L, c("src", "tgt"), c("src_id", "tgt_id"))

long_cols <- c("dataset","edge_dir","cat","edge_class","lane","src_id","tgt_id","method","nsnp",
               "b","se","pval","pval_adj","mr_hit",
               "het_Q","het_Qdf","het_pval","het_flag","het_clean",
               "egger_intercept","egger_pval","pleio_flag","pleio_clean",
               "steiger_dir","steiger_pval","steiger_flag","steiger_ok",
               "direction_established","direction_unresolved",
               "presso_global_p","presso_corr_p","wm_p",
               "clean","rescued_presso","rescued_median","robust_pass","het_status",
               "sens_pass","hit_after_sens","mr_tier","tier_reason","coloc_status","PP_H4")
fwrite(L[, ..long_cols], file.path(CFG$out, "mr_sensitivity_long.tsv"), sep = "\t")
ts("wrote mr_sensitivity_long.tsv (", nrow(L), " edges)")

# ---- sensitivity summaries (per edge_dir + overall) ------------------------
mk_summary <- function(dt, by) dt[, {
  nh  <- sum(mr_hit, na.rm = TRUE)
  nha <- sum(hit_after_sens, na.rm = TRUE)
  nwh <- sum(mr_hit & !is.na(het_pval),   na.rm = TRUE)
  nwp <- sum(mr_hit & !is.na(egger_pval), na.rm = TRUE)
  .(n_edges = .N, n_hits = nh, n_hits_after = nha,
    n_hit_with_het = nwh, n_hit_with_ple = nwp,
    n_het_flag = sum(mr_hit & !is.na(het_pval)   & het_flag,   na.rm = TRUE),
    n_ple_flag = sum(mr_hit & !is.na(egger_pval) & pleio_flag, na.rm = TRUE),
    hit_retention_pct = if (nh > 0) 100 * nha / nh else NA_real_)
}, by = by]
by_edge <- mk_summary(L, c("dataset", "edge_dir"))
by_edge[, het_flag_rate_pct := fifelse(n_hit_with_het > 0, 100 * n_het_flag / n_hit_with_het, NA_real_)]
by_edge[, ple_flag_rate_pct := fifelse(n_hit_with_ple > 0, 100 * n_ple_flag / n_hit_with_ple, NA_real_)]
fwrite(by_edge, file.path(CFG$out, "sensitivity_by_edgedir.tsv"), sep = "\t")
fwrite(mk_summary(L, "dataset"), file.path(CFG$out, "sensitivity_overall.tsv"), sep = "\t")
ts("wrote sensitivity_by_edgedir.tsv + sensitivity_overall.tsv")

# ============================================================================
# WIDE MOTIF TABLE  (one row per E-P-D triad; ports MRviz.R motif logic)
# ============================================================================
triads <- fread(file.path(edges_dir, "mr_triads.tsv"), showProgress = FALSE)
keep_meta <- intersect(c("ExposureCategory","Disease_UKB","ICD10"), names(triads))
trip <- triads[, .SD[1L], by = .(Exposure, Protein, Disease), .SDcols = keep_meta]
ts("triad master: ", nrow(trip), " (E,P,D) triads")

# per-category wide slice from L (primary estimate keyed to triad coordinates)
sl <- function(catname, src_name, tgt_name, tag) {
  x <- L[cat == catname, .(src_id, tgt_id, b, se, pval_adj)]
  setnames(x, c("src_id","tgt_id","b","se","pval_adj"),
           c(src_name, tgt_name, paste0("beta_", tag), paste0("se_", tag), paste0("padj_", tag)))
  unique(x)
}
m <- trip
m <- merge(m, sl("EP",      "Exposure","Protein", "EP"),      by = c("Exposure","Protein"), all.x = TRUE)
m <- merge(m, sl("PDcis",   "Protein","Disease",  "PDcis"),   by = c("Protein","Disease"),  all.x = TRUE)
m <- merge(m, sl("PDtrans", "Protein","Disease",  "PDtrans"), by = c("Protein","Disease"),  all.x = TRUE)
m <- merge(m, sl("ED",      "Exposure","Disease", "ED"),      by = c("Exposure","Disease"), all.x = TRUE)
m <- merge(m, sl("PEcis",   "Protein","Exposure", "PEcis"),   by = c("Protein","Exposure"), all.x = TRUE)
m <- merge(m, sl("PEtrans", "Protein","Exposure", "PEtrans"), by = c("Protein","Exposure"), all.x = TRUE)
m <- merge(m, sl("DP",      "Disease","Protein",  "DP"),      by = c("Disease","Protein"),  all.x = TRUE)
m <- merge(m, sl("DE",      "Disease","Exposure", "DE"),      by = c("Disease","Exposure"), all.x = TRUE)

# edge state: NA -> "0"; sig & beta>0 -> "+"; sig & beta<0 -> "-"; else "0"
edge_state <- function(beta, padj) {
  out <- rep("0", length(beta))
  ok <- !is.na(beta) & !is.na(padj) & padj < ALPHA_Q
  out[ok & beta > 0] <- "+"
  out[ok & beta < 0] <- "-"
  out
}
EDGES8 <- c("EP","PDcis","PDtrans","ED","PEcis","PEtrans","DP","DE")
for (e in EDGES8)
  m[, (paste0("state_", e)) := edge_state(get(paste0("beta_", e)), get(paste0("padj_", e)))]

is_present <- function(x) !is.na(x) & x != "0"
collapse_any <- function(a, b) fifelse(a == "+" | b == "+", "+", fifelse(a == "-" | b == "-", "-", "0"))
m[, state_PDany := collapse_any(state_PDcis, state_PDtrans)]
m[, state_PEany := collapse_any(state_PEcis, state_PEtrans)]
m[, pres_EP := is_present(state_EP)]
m[, pres_PD := is_present(state_PDany)]
m[, pres_ED := is_present(state_ED)]
m[, pres_PE := is_present(state_PEany)]
m[, pres_DP := is_present(state_DP)]
m[, pres_DE := is_present(state_DE)]

# motif definitions (verbatim from MRviz.R)
m[, motif_A_mediator        :=  pres_EP &  pres_PD &  pres_ED & !pres_PE & !pres_DP & !pres_DE]
m[, motif_B_biomarker       :=  pres_EP & !pres_PD &  pres_ED & !pres_PE &  pres_DP & !pres_DE]
m[, motif_C_exposure_marker :=  pres_EP & !pres_PD &  pres_ED & !pres_PE & !pres_DP]
m[, motif_D_P_to_E          := !pres_EP &  pres_PE]
m[, motif_E_disease_liability :=  pres_DP &  pres_DE]

m[, signature8 := paste0(state_EP, state_PDcis, state_PDtrans, state_ED,
                         state_PEcis, state_PEtrans, state_DP, state_DE)]
state_cols <- paste0("state_", EDGES8)
m[, any_sig := Reduce(`|`, lapply(state_cols, function(c) m[[c]] %in% c("+","-")))]
m[, n_motifs := motif_A_mediator + motif_B_biomarker + motif_C_exposure_marker +
                motif_D_P_to_E + motif_E_disease_liability]
m[, motif_label := fcase(
  !any_sig,                  "Null (no MR hits)",
  n_motifs > 1,              "Multiple motifs",
  motif_A_mediator,          "A",
  motif_B_biomarker,         "B",
  motif_C_exposure_marker,   "C",
  motif_D_P_to_E,            "D",
  motif_E_disease_liability, "E",
  default =                  "Other (has hits)")]
m[, triplet := paste(Exposure, Protein, Disease, sep = " | ")]
m[, dataset := CFG$tag]

fwrite(m, file.path(CFG$out, "MRmotifs.tsv"), sep = "\t")
ts("wrote MRmotifs.tsv (", nrow(m), " triads)")

# ---- console roll-up -------------------------------------------------------
cat("\n================ ", cohort, " MR TABLE SUMMARY ================\n", sep = "")
cat("Edges with a primary estimate, by edge_dir:\n")
print(by_edge[, .(edge_dir, n_edges, n_hits, n_hits_after, hit_retention_pct = round(hit_retention_pct, 1))])
cat("\nTier distribution by lane (per-arm; Tier1+ added cross-arm):\n")
print(dcast(L, lane ~ factor(mr_tier, levels = c("Null","Suggestive","Reverse","Tier2","Tier1")),
            value.var = "src_id", fun.aggregate = length))
cat("\nMotif label distribution:\n")
print(m[, .N, by = motif_label][order(-N)])
cat("\nWrote to: ", CFG$out, "\nDONE.\n", sep = "")
