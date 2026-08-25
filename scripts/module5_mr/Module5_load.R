#!/usr/bin/env Rscript
# Module5_load.R -- Build MR priority edge lists from Module 3 mediation results.
#
# Reads:
#   analysis_plan.tsv          -> which Module 3 experiments to load
#   Module 3 MDres files       -> partitioned_categories long-format results
#   ReplicatedEassoc.csv       -> replicated exposure-protein associations (Module 2)
#   UKBFinnGenDisease.csv      -> UKB first-occurrence ICD10 -> FinnGen disease ID map
#
# Produces edge files under heap_project_output("mr_edges", "global_edges"):
#   MR_priority_table.tsv      -- full protein x disease priority table
#   edges_PD.tsv               -- all P->D pairs (any significant mediation or protein_p)
#   edges_PD_Gcis.tsv          -- P->D pairs where cis-genetic NIE is significant
#   edges_PD_Gtrans.tsv        -- P->D pairs where trans-genetic NIE is significant
#   edges_PD_anyG.tsv          -- P->D pairs where either Gcis or Gtrans NIE is significant
#   edges_EPD_categories.tsv   -- E(category) x P x D triples: PXS category NIE significant
#   edges_EP.tsv               -- E->P pairs from assoc, filtered to proteins in PD set
#   edges_ED.tsv               -- E->D pairs derived from EPD category edges
#   HEAPres.tsv                -- full expanded table with FinnGen IDs (legacy format)

local({
  candidates <- c(
    Sys.getenv("HEAP_PATHS_FILE", unset = ""),
    file.path(getwd(), "workflow", "00_paths.R"),
    file.path(getwd(), "..", "workflow", "00_paths.R"),
    file.path(getwd(), "..", "..", "workflow", "00_paths.R"),
    file.path(getwd(), "..", "..", "..", "workflow", "00_paths.R")
  )
  candidates <- candidates[nzchar(candidates)]
  hit <- candidates[file.exists(candidates)][1]
  if (!is.na(hit)) source(hit)
})

# config_helpers.R lives alongside 00_paths.R; needed for load_module_experiments()
local({
  candidates <- c(
    if (exists("HEAP_PATHS") && !is.null(HEAP_PATHS$heap_root))
      file.path(HEAP_PATHS$heap_root, "workflow", "config_helpers.R") else "",
    file.path(getwd(), "workflow", "config_helpers.R"),
    file.path(getwd(), "..", "workflow", "config_helpers.R"),
    file.path(getwd(), "..", "..", "workflow", "config_helpers.R"),
    file.path(getwd(), "..", "..", "..", "workflow", "config_helpers.R")
  )
  candidates <- candidates[nzchar(candidates)]
  hit <- candidates[file.exists(candidates)][1]
  if (!is.na(hit)) source(hit)
})

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
})

ALPHA <- 0.05

# ============================================================
# 1) Resolve the main partitioned Module 3 experiments from config
#
# Source of truth is config/modules/module3_experiments.yml. The old
# analysis_plan.tsv (Type1-7 scheme) is retired after the covariate-set
# restructure (covariate sets are now descriptive: base, base_clinical, ...).
# We take every experiment whose priority == "main", status == "ready", and
# mediation_mode == "partitioned_categories", then union their results below.
# ============================================================

m3_cfg <- load_module_experiments("module3")$experiments

m3_keep <- Filter(function(e) {
  identical(e$priority, "main") &&
  identical(e$status, "ready") &&
  identical(e$mediation_mode, "partitioned_categories")
}, m3_cfg)

if (length(m3_keep) == 0)
  stop("No main partitioned_categories Module 3 experiments (priority=main, ",
       "status=ready) in config/modules/module3_experiments.yml")

m3_exps <- rbindlist(lapply(names(m3_keep), function(nm) {
  e <- m3_keep[[nm]]
  data.table(
    experiment_id  = nm,
    covarType      = e$covariate_set,
    family         = e$family,
    mediation_mode = e$mediation_mode
  )
}))

message("Module 3 experiments (config-resolved, main partitioned):")
print(m3_exps)

# ============================================================
# 2) Load and combine all MDres files
# ============================================================

load_mdres <- function(experiment_id, covarType, family, mediation_mode) {
  # Module 3 (manifest mode) writes results under:
  #   heap_project_output("module3", <experiment_id>, <covarType>, <family>, <mediation_mode>)
  # IGLOO canonical, fall back to local staging if present.
  m3_root <- heap_project_output("module3")
  if (!dir.exists(m3_root) && dir.exists(heap_output("module3")))
    m3_root <- heap_output("module3")
  score_dir <- file.path(m3_root, experiment_id, covarType, family, mediation_mode)
  files <- list.files(score_dir, pattern = "^MDres_.*\\.txt$", full.names = TRUE)
  if (length(files) == 0) {
    warning("No MDres files found — skipping: ", score_dir)
    return(NULL)
  }
  message("  Loading ", length(files), " MDres file(s) from ", score_dir)
  rbindlist(lapply(files, fread), fill = TRUE)
}

md_list <- lapply(seq_len(nrow(m3_exps)), function(i) {
  exp <- m3_exps[i]
  dt  <- load_mdres(exp$experiment_id, exp$covarType, exp$family, exp$mediation_mode)
  if (is.null(dt)) return(NULL)
  dt[, experiment_id := exp$experiment_id]
  dt[, covarType     := exp$covarType]
  dt
})
md_list <- Filter(Negate(is.null), md_list)

if (length(md_list) == 0)
  stop("No Module 3 results found for any ready experiment. Run Module3.R first.")

md_all <- rbindlist(md_list, fill = TRUE)

message("Combined: ", nrow(md_all), " rows across ",
        length(unique(md_all$experiment_id)), " experiment(s)")

# ============================================================
# 2b) PXS_total primary mediation -> drives MR TRIAD SELECTION
#
# The partitioned (per-category) results above are kept as DESCRIPTIVE
# provenance (priority table + edges_EPD_categories) to explain WHICH lifestyle
# category contributes to a mediation. But the MR triad set is selected on the
# TOTAL-exposome NIE (PXS_total): the hypothesis "does the exposome, in
# aggregate, mediate P->D". Exposures are then resolved by the replicated E->P
# associations (any category). Source = main primary_total Module 3 experiments.
# ============================================================
m3_primary <- Filter(function(e) {
  identical(e$priority, "main") && identical(e$status, "ready") &&
  identical(e$mediation_mode, "primary_total")
}, m3_cfg)
if (length(m3_primary) == 0)
  stop("No main primary_total Module 3 experiment (PXS_total) in config — ",
       "needed for MR triad selection.")
md_total <- rbindlist(lapply(names(m3_primary), function(nm) {
  e  <- m3_primary[[nm]]
  dt <- load_mdres(nm, e$covariate_set, e$family, e$mediation_mode)
  if (is.null(dt)) return(NULL)
  dt[, experiment_id := nm][]
}), fill = TRUE)
if (is.null(md_total) || !nrow(md_total))
  stop("No primary_total MDres found (PXS_total) — run the M3 primary experiment first.")
nie_tot <- md_total[effect_type == "NIE" & predictor == "PXS_total" & !is.na(delta_p)]
nie_tot[, bonf := ALPHA / uniqueN(paste(protID, DZ_ID)), by = experiment_id]
nie_tot[, sig := delta_p < bonf]
# union of significant PXS_total NIE pairs across the primary experiments
pxs_total_sig <- unique(nie_tot[sig == TRUE, .(protID, DZ_ID)])
message("PXS_total NIE-significant P-D pairs (MR triad anchor): ", nrow(pxs_total_sig),
        " (", uniqueN(pxs_total_sig$protID), " proteins, ",
        uniqueN(pxs_total_sig$DZ_ID), " diseases)")

# ============================================================
# 3) Significance filtering (NIE only)
#
# Bonferroni correction per predictor per experiment:
#   threshold = ALPHA / n_unique(protID x DZ_ID) for each predictor
#
# A protein-disease pair is significant if it passes the threshold in
# ANY experiment (union across covariate specifications).
# ============================================================

nie_df <- md_all[effect_type == "NIE"]

# Compute per-predictor Bonferroni threshold within each experiment
nie_df[, n_tests := uniqueN(paste(protID, DZ_ID)), by = .(experiment_id, predictor)]
nie_df[, bonf_threshold := ALPHA / n_tests]
nie_df[, sig := !is.na(delta_p) & delta_p < bonf_threshold]

# Union significant hits across experiments: sig in ANY experiment
sig_any <- nie_df[sig == TRUE,
                  .(experiment_id, covarType, protID, DZ_ID, predictor,
                    predictor_class, effect_logHR, effect_HR,
                    delta_se, delta_l95, delta_u95, delta_p,
                    protein_HR, protein_HR_l95, protein_HR_u95, protein_p,
                    mediator_adjR2, cox_cindex, n, n_cases)]

message("Significant NIE hits (union across experiments): ",
        uniqueN(paste(sig_any$protID, sig_any$DZ_ID)),
        " unique protein-disease pairs, ",
        nrow(sig_any), " predictor-specific rows")

# Also flag protein_p significant pairs (direct protein-disease association)
prot_sig <- md_all[effect_type == "NIE",
                   .(protID, DZ_ID, protein_HR, protein_p,
                     n_tests_pd = uniqueN(paste(protID, DZ_ID))),
                   by = .(experiment_id)][
  , bonf_pd := ALPHA / n_tests_pd][
  , prot_sig := !is.na(protein_p) & protein_p < bonf_pd]

prot_sig_pairs <- unique(prot_sig[prot_sig == TRUE, .(protID, DZ_ID,
                                                       protein_HR, protein_p)])

# ============================================================
# 4) Build the priority summary table
#    One row per protein x disease; columns flag which predictors are significant
# ============================================================

# Pivot significant predictor flags wide
pred_flags <- dcast(
  sig_any[, .(protID, DZ_ID, predictor, sig = TRUE)],
  protID + DZ_ID ~ predictor,
  value.var = "sig",
  fill = FALSE
)

# Join protein_p info — aggregate the direct protein->disease Cox stats PER
# (protein, disease). The `by` is essential: without it, mean()/max() collapse to
# a single global scalar that data.table recycles across every row, making
# protein_HR/protein_p constant for all 451k pairs (the bug that flat-lined
# fig_mr_priority). The NIE rows repeat the same direct stats across predictors,
# so the per-pair mean just recovers that pair's value.
prot_info <- md_all[effect_type == "NIE",
                    .(protein_HR = mean(protein_HR, na.rm = TRUE),
                      protein_p = mean(protein_p, na.rm = TRUE),
                      mediator_adjR2 = mean(mediator_adjR2, na.rm = TRUE),
                      cox_cindex = mean(cox_cindex, na.rm = TRUE),
                      n = max(n), n_cases = max(n_cases)),
                    by = .(protID, DZ_ID)]

priority_tbl <- merge(pred_flags, prot_info, by = c("protID", "DZ_ID"), all = TRUE)
priority_tbl[is.na(priority_tbl)] <- FALSE

# sig_any_predictor: TRUE if any NIE predictor significant
pred_cols <- setdiff(names(priority_tbl),
                     c("protID", "DZ_ID", "protein_HR", "protein_p",
                       "mediator_adjR2", "cox_cindex", "n", "n_cases"))
priority_tbl[, sig_any_predictor := rowSums(.SD) > 0, .SDcols = pred_cols]
# Direct protein-disease association significant for THIS pair (pairwise match;
# the previous `protID %in% ... & DZ_ID %in% ...` form was a cartesian over-match
# that flagged any protein-in-set x disease-in-set combination).
prot_sig_keys <- paste(prot_sig_pairs$protID, prot_sig_pairs$DZ_ID)
priority_tbl[, sig_protein_direct := paste(protID, DZ_ID) %in% prot_sig_keys]
priority_tbl[, mr_priority := sig_any_predictor | sig_protein_direct]

message("Priority table: ", nrow(priority_tbl), " protein-disease pairs assessed, ",
        sum(priority_tbl$mr_priority), " flagged for MR")

# ============================================================
# 5) Build edge lists
# ============================================================

# P->D: any significant mediation OR direct protein association
edges_PD <- unique(priority_tbl[mr_priority == TRUE, .(Protein = protID, Disease = DZ_ID)])

# P->D stratified by genetic instrument type
g_cols_cis   <- grep("^Gcis_raw$",   names(priority_tbl), value = TRUE)
g_cols_trans <- grep("^Gtrans_raw$", names(priority_tbl), value = TRUE)

edges_PD_Gcis <- if (length(g_cols_cis) > 0 && any(priority_tbl[[g_cols_cis]])) {
  unique(priority_tbl[get(g_cols_cis) == TRUE, .(Protein = protID, Disease = DZ_ID)])
} else data.table(Protein = character(), Disease = character())

edges_PD_Gtrans <- if (length(g_cols_trans) > 0 && any(priority_tbl[[g_cols_trans]])) {
  unique(priority_tbl[get(g_cols_trans) == TRUE, .(Protein = protID, Disease = DZ_ID)])
} else data.table(Protein = character(), Disease = character())

edges_PD_anyG <- unique(rbind(edges_PD_Gcis, edges_PD_Gtrans))

# E(category)->P->D: PXS_<category> NIE significant
pxs_cat_cols <- grep("^PXS_", names(priority_tbl), value = TRUE)

epd_rows <- lapply(pxs_cat_cols, function(col) {
  cat_name <- sub("^PXS_", "", col)
  dt <- priority_tbl[get(col) == TRUE,
                     .(ExposureCategory = cat_name,
                       Protein = protID,
                       Disease = DZ_ID)]
  dt
})
edges_EPD_categories <- rbindlist(epd_rows, fill = TRUE)

# E->P: from ReplicatedEassoc, restricted to proteins appearing in PD edges.
# NOTE (migration): ReplicatedEassoc.csv is a HEAP-DERIVED table (replicated
# exposure->protein associations), not raw data. Canonical home is an IGLOO HEAP
# output; the producing step still needs to be wired (open decision #5 in
# docs/NON_VISUALIZATION_DEPENDENCY_AUDIT.md). Prefer the IGLOO copy, fall back
# to the legacy App/Tables copy until the producer is migrated.
assoc_igloo <- heap_project_output("module2", "ReplicatedEassoc.csv")
assoc_path  <- if (file.exists(assoc_igloo)) assoc_igloo else
  legacy_ukb_path("Output", "App", "Tables", "ReplicatedEassoc.csv")

if (file.exists(assoc_path)) {
  assoc <- fread(assoc_path)
  # Filter to: proteins in the MR-priority PD set
  # AND: exposure categories that have at least one significant EPD edge
  sig_cats <- unique(edges_EPD_categories$ExposureCategory)
  sig_prots <- unique(edges_PD$Protein)

  edges_EP <- unique(assoc[
    omicID %in% sig_prots & Category_train %in% sig_cats,
    .(Exposure = Eid_train, ExposureCategory = Category_train, Protein = omicID)
  ])
  message("E->P edges (assoc, filtered to MR-priority proteins + significant categories): ",
          nrow(edges_EP))
  # Ungated replicated E->P (ALL categories) — used for PXS_total triad selection.
  edges_EP_all <- unique(assoc[, .(Exposure = Eid_train,
                                   ExposureCategory = Category_train, Protein = omicID)])
} else {
  warning("ReplicatedEassoc.csv not found at: ", assoc_path)
  edges_EP <- data.table(Exposure = character(), ExposureCategory = character(),
                          Protein = character())
  edges_EP_all <- copy(edges_EP)
}

# E->D: derived from EPD edges (collapse protein dimension)
edges_ED <- unique(edges_EPD_categories[, .(ExposureCategory, Disease)])

# ============================================================
# 6) Map UKB disease IDs -> FinnGen IDs
# ============================================================

finngen_path <- igloo_path("FinnGen", "UKBFinnGenDisease.csv")
if (file.exists(finngen_path)) {
  fg_map <- fread(finngen_path)
  setnames(fg_map, c("Disease", "regex_1", "FinnGen"), c("Disease_UKB", "ICD10", "FinnGen"))

  # Some UKB diseases map to MULTIPLE FinnGen phenotypes, stored ";"-separated
  # (e.g. "T2D; T2D_WIDE", "E4_OBESITY; E4_OBESITYCAL; E4_OBESITYNAS"). Split
  # into one row per FinnGen id so each becomes a separate, file-resolvable
  # outcome (<id>.gz) instead of an invalid concatenated id.
  fg_map <- fg_map[, .(FinnGen = trimws(unlist(strsplit(as.character(FinnGen), ";", fixed = TRUE)))),
                   by = .(Disease_UKB, ICD10)]
  fg_map <- unique(fg_map[nzchar(FinnGen)])

  map_to_finngen <- function(dt, disease_col = "Disease") {
    dt <- copy(dt)
    setnames(dt, disease_col, "Disease_UKB")
    # 1-to-many: a UKB disease can map to several FinnGen phenotypes (split above).
    dt <- merge(dt, fg_map, by = "Disease_UKB", all.x = TRUE, allow.cartesian = TRUE)
    dt <- dt[!is.na(FinnGen) & FinnGen != ""]
    dt[, FinnGen := paste0("finngen_R12_", FinnGen)]
    dt
  }

  edges_PD_fg      <- map_to_finngen(edges_PD,      "Disease")
  edges_PD_Gcis_fg <- map_to_finngen(edges_PD_Gcis, "Disease")
  edges_PD_Gtrans_fg <- map_to_finngen(edges_PD_Gtrans, "Disease")
  edges_PD_anyG_fg <- map_to_finngen(edges_PD_anyG, "Disease")
  edges_EPD_fg     <- map_to_finngen(edges_EPD_categories, "Disease")
  edges_ED_fg      <- map_to_finngen(edges_ED,       "Disease")

  message("FinnGen mapping: ",
          nrow(edges_PD_fg), " P-D pairs, ",
          nrow(edges_EPD_fg), " E-P-D triples with FinnGen IDs")
} else {
  warning("UKBFinnGenDisease.csv not found; FinnGen mapping skipped")
  edges_PD_fg <- edges_PD_Gcis_fg <- edges_PD_Gtrans_fg <-
    edges_PD_anyG_fg <- edges_EPD_fg <- edges_ED_fg <- NULL
}

# ============================================================
# 7) Write outputs
# ============================================================

# IGLOO canonical — matches where Module5.R / Module5_deCODE.R read edges from.
outdir <- heap_project_output("mr_edges", "global_edges")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ---- 7a) Analytic / provenance tables (NOT consumed directly by the runners) --
fwrite(priority_tbl,
       file.path(outdir, "MR_priority_table.tsv"), sep = "\t")

# Full priority P->D union (UKB disease names) — reference only.
fwrite(edges_PD,
       file.path(outdir, "edges_PD_all_priority_ukb.tsv"), sep = "\t")
fwrite(edges_PD_Gcis,
       file.path(outdir, "edges_PD_Gcis.tsv"), sep = "\t")
fwrite(edges_PD_Gtrans,
       file.path(outdir, "edges_PD_Gtrans.tsv"), sep = "\t")
fwrite(edges_PD_anyG,
       file.path(outdir, "edges_PD_anyG.tsv"), sep = "\t")
fwrite(edges_EPD_categories,
       file.path(outdir, "edges_EPD_categories.tsv"), sep = "\t")
fwrite(edges_ED,
       file.path(outdir, "edges_ED_categories.tsv"), sep = "\t")

# FinnGen-mapped provenance variants
if (!is.null(edges_PD_fg)) {
  fwrite(edges_PD_fg,
         file.path(outdir, "edges_PD_finngen.tsv"), sep = "\t")
  fwrite(edges_PD_Gcis_fg,
         file.path(outdir, "edges_PD_Gcis_finngen.tsv"), sep = "\t")
  fwrite(edges_PD_Gtrans_fg,
         file.path(outdir, "edges_PD_Gtrans_finngen.tsv"), sep = "\t")
  fwrite(edges_PD_anyG_fg,
         file.path(outdir, "edges_PD_anyG_finngen.tsv"), sep = "\t")
  fwrite(edges_EPD_fg,
         file.path(outdir, "edges_EPD_categories_finngen.tsv"), sep = "\t")
  fwrite(edges_ED_fg,
         file.path(outdir, "edges_ED_finngen.tsv"), sep = "\t")

  # HEAPres: Exposure + Protein + Disease(UKB/ICD10/FinnGen) triple provenance.
  heap_res <- merge(
    edges_EPD_fg[, .(ExposureCategory, Protein, Disease_UKB, ICD10, FinnGen)],
    edges_EP[, .(ExposureCategory, Exposure, Protein)],
    by = c("ExposureCategory", "Protein"), all.x = TRUE, allow.cartesian = TRUE
  )
  heap_res <- na.omit(heap_res)
  heap_res[, FinnGen := sub("finngen_R12_", "", FinnGen)]
  fwrite(heap_res, file.path(outdir, "HEAPres.tsv"), sep = "\t")
}

# ---- 7b) CANONICAL RUNNER-FACING EDGE FILES -------------------------------
# These are what Module5.R / Module5_deCODE.R read as edges_<EDGE_TYPE>.tsv.
# Schema the runners require:
#   E<->P : columns Exposure, Protein         (exposure = REGENIE exposure id)
#   P<->D : columns Protein,  Disease         (disease  = FinnGen id, <id>.gz)
#   E<->D : columns Exposure, Disease         (FinnGen id)
# Each undirected pair is written under BOTH directional names so the
# bidirectional network (EP/PE, PD/DP, ED/DE) is fully covered.
if (is.null(edges_EPD_fg)) {
  warning("No FinnGen mapping — canonical runner-facing edge files NOT written.")
} else {
  # ---- MR-SELECTION: TOTAL-exposome mediation (PXS_total NIE) ---------------
  # P<->D pairs are selected on PXS_total NIE significance (Bonferroni; from the
  # primary_total Module 3 experiment, computed in section 2b) — the hypothesis
  # "does the exposome, in aggregate, mediate P->D". The per-category EPD edges
  # and the genetic Gcis/Gtrans edges remain DESCRIPTIVE provenance
  # (edges_EPD_categories*, edges_PD_*G*); they no longer drive MR selection.
  pd_total_fg <- map_to_finngen(pxs_total_sig[, .(protID, Disease = DZ_ID)], "Disease")
  run_PD <- unique(pd_total_fg[, .(Protein = protID, Disease = FinnGen)])

  # E<->P: ALL replicated exposure->protein associations (any category) for the
  # PXS_total-prioritized proteins — the exposome-aggregate hypothesis lets MR
  # test each replicated exposure that acts on the protein.
  prots_total <- unique(pd_total_fg$protID)
  run_EP <- unique(edges_EP_all[Protein %in% prots_total,
                                .(Exposure, ExposureCategory, Protein)])
  # Drop exposures with NO REGENIE GWAS — they cannot be MR-instrumented (E->P) or
  # used as an outcome (P->E), and the runner stop()s on a missing GWAS, crashing
  # the whole chunk. (e.g. the degenerate activity one-hot multi-level.) Removing
  # them here keeps un-runnable E<->P / E<->D edges out of the run entirely.
  .expdir <- heap_project_output("gwas", "regenie_step2")
  .n_ep0 <- nrow(run_EP)
  run_EP <- run_EP[file.exists(file.path(.expdir, Exposure, paste0(Exposure, ".regenie")))]
  message("E->P edges after dropping GWAS-less exposures: ", nrow(run_EP),
          " (removed ", .n_ep0 - nrow(run_EP), "; exposures kept: ", uniqueN(run_EP$Exposure), ")")

  # Fully-specified E->P->D triads = PXS_total-significant (P,D) x replicated
  # E->P for that protein (joined on Protein; NOT category-gated). THE documented
  # "what triads the MR tests" table.
  mr_triads <- if (nrow(run_EP) > 0) {
    tri <- merge(
      pd_total_fg[, .(Protein = protID, Disease = FinnGen, Disease_UKB, ICD10)],
      run_EP[, .(Exposure, ExposureCategory, Protein)],
      by = "Protein", allow.cartesian = TRUE)
    unique(tri[!is.na(Exposure),
               .(Exposure, ExposureCategory, Protein, Disease, Disease_UKB, ICD10)])
  } else data.table()

  # E<->D: collapse the protein dimension of the triads.
  run_ED <- if (nrow(mr_triads) > 0)
    unique(mr_triads[, .(Exposure, Disease)]) else
    data.table(Exposure = character(), Disease = character())

  fwrite(run_PD, file.path(outdir, "edges_PD.tsv"), sep = "\t")
  fwrite(run_PD, file.path(outdir, "edges_DP.tsv"), sep = "\t")
  fwrite(run_EP, file.path(outdir, "edges_EP.tsv"), sep = "\t")
  fwrite(run_EP, file.path(outdir, "edges_PE.tsv"), sep = "\t")
  fwrite(run_ED, file.path(outdir, "edges_ED.tsv"), sep = "\t")
  fwrite(run_ED, file.path(outdir, "edges_DE.tsv"), sep = "\t")
  fwrite(mr_triads, file.path(outdir, "mr_triads.tsv"), sep = "\t")

  cat("\nCanonical runner-facing edge files (FinnGen disease ids):\n",
      "  edges_PD / edges_DP :", nrow(run_PD), "protein<->disease pairs (exposome E->P->D mediation anchored)\n",
      "  edges_EP / edges_PE :", nrow(run_EP), "exposure<->protein pairs (replicated)\n",
      "  edges_ED / edges_DE :", nrow(run_ED), "exposure<->disease pairs\n",
      "  mr_triads           :", nrow(mr_triads), "fully-specified E->P->D triads (connection inventory)\n",
      sep = "")
}

cat("\nAnalytic edge summary (provenance):\n",
    "  MR_priority_table        :", nrow(priority_tbl), "protein-disease pairs assessed\n",
    "  edges_PD_all_priority_ukb:", nrow(edges_PD), "(any sig mediation OR direct, UKB names)\n",
    "  edges_PD_Gcis            :", nrow(edges_PD_Gcis), "(cis-genetic NIE significant)\n",
    "  edges_PD_Gtrans          :", nrow(edges_PD_Gtrans), "(trans-genetic NIE significant)\n",
    "  edges_PD_anyG            :", nrow(edges_PD_anyG), "(Gcis OR Gtrans NIE significant)\n",
    "  edges_EPD_categories     :", nrow(edges_EPD_categories), "(category-level E->P->D triples)\n",
    sep = "")
