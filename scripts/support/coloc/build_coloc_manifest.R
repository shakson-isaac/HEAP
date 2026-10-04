#!/usr/bin/env Rscript
# ============================================================================
# support/coloc/build_coloc_manifest.R
# ----------------------------------------------------------------------------
# Assemble the colocalization MANIFEST: one row per cis-pQTL locus whose
# coloc_status is currently "pending" in the Module-5 MR sensitivity tables, so
# a systematic coloc.abf pass (run_coloc_locus.R) can replace "pending" with a
# PP.H4 >= 0.8 hard gate.
#
# Source of truth for the pending loci:
#   output/mr_edges/summary/mr_sensitivity_long.tsv         (UKB arm)
#   output/mr_edges/summary/DECODE/mr_sensitivity_long.tsv  (deCODE arm)
# Filter: coloc_status == "pending" on the cis-pQTL lane (lane==pQTL &
# edge_class==cis). Loci are src_id (protein) x tgt_id (FinnGen disease for
# P->D, UKB exposure for P->E). ~57 UKB + ~16 deCODE = ~73.
#
# Each row resolves the four inputs the runner needs and records WHY a locus is
# (un)runnable in resolve_status:
#   OK              -> all inputs present, runnable
#   missing_pqtl    -> no UKB parquet / deCODE somascan file for the protein
#   missing_gwas    -> no FinnGen .gz (P->D) or REGENIE .regenie (P->E) outcome
#   missing_lead    -> no cis instrument clump for the protein (lead SNP/window)
#   not_casecontrol -> P->E exposure-outcome locus has no FinnGen file; coloc.abf
#                      still runs (type="quant") but flagged so callers can skip
#                      the disease cc-gate if desired. (Reported, still runnable.)
#
# Output -> output/support/coloc/coloc_manifest.tsv
#
# Usage:
#   module load gcc/14.2.0 R/4.4.2
#   HEAP_PATHS_FILE=.../workflow/00_paths.R \
#     Rscript scripts/support/coloc/build_coloc_manifest.R
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

WINDOW_KB <- as.numeric(Sys.getenv("COLOC_WINDOW_KB", "500"))

summ_dir  <- heap_project_output("mr_edges", "summary")
out_dir   <- heap_project_output("support", "coloc")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

inst_ukb  <- heap_project_output("mr", "protein_inst")
inst_dec  <- heap_project_output("mr", "protein_inst_decode")
UKB_PQTL  <- "/n/groups/patel/IGLOO/UKB/pQTL"
DEC_PQTL  <- "/n/groups/patel/IGLOO/DECODE/pQTL/final_somascan_smp"
DEC_MAP   <- fread("/n/groups/patel/IGLOO/DECODE/pQTLmetadata/somascan_protein_map.tsv")
FINNGEN   <- "/n/groups/patel/IGLOO/FinnGen/SummaryStats"
FG_MAN    <- fread("/n/groups/patel/IGLOO/FinnGen/finngen_R12_manifest.tsv")
EXPO_GWAS <- heap_gwas("regenie_step2")

norm_chr <- function(x) gsub("^chr", "", as.character(x))

# ---- collect pending cis loci from both arms -------------------------------
read_pending <- function(fp) {
  if (!file.exists(fp)) { message("WARN: sensitivity file not found: ", fp); return(NULL) }
  d <- fread(fp)
  need <- c("dataset", "edge_dir", "lane", "edge_class", "src_id", "tgt_id", "coloc_status")
  miss <- setdiff(need, names(d))
  if (length(miss)) stop("sensitivity file missing cols: ", paste(miss, collapse = ", "))
  d <- d[lane == "pQTL" & edge_class == "cis" & coloc_status == "pending",
         .(arm = dataset, edge_dir, protID = src_id, disease_or_exposure_id = tgt_id)]
  unique(d)
}

pend <- rbindlist(list(
  read_pending(file.path(summ_dir, "mr_sensitivity_long.tsv")),
  read_pending(file.path(summ_dir, "DECODE", "mr_sensitivity_long.tsv"))
), use.names = TRUE)
pend <- unique(pend)
if (!nrow(pend)) stop("No pending cis-pQTL loci found in either sensitivity table.")
message(sprintf("Pending cis-pQTL loci: %d (UKB=%d, DECODE=%d)",
                nrow(pend), sum(pend$arm == "UKB"), sum(pend$arm == "DECODE")))

# ---- resolvers (mirror run_coloc_shortlist.R) ------------------------------
# cis lead SNP + window from the clumped instrument cache (smallest p)
resolve_lead <- function(protID, arm) {
  f <- file.path(if (arm == "UKB") inst_ukb else inst_dec,
                 paste0("protein_", protID, "_cis_clumped.tsv"))
  if (!file.exists(f)) return(list(lead = NA_character_, chr = NA_character_, pos = NA_integer_, ok = FALSE))
  cl <- fread(f)
  if (!nrow(cl) || !("pval.exposure" %in% names(cl))) return(list(lead = NA_character_, chr = NA_character_, pos = NA_integer_, ok = FALSE))
  cl <- cl[order(pval.exposure)][1]
  lead <- if ("SNP" %in% names(cl)) as.character(cl$SNP[1]) else NA_character_
  list(lead = lead, chr = norm_chr(cl$chr.exposure[1]), pos = as.integer(cl$pos.exposure[1]),
       ok = !is.na(lead) && nzchar(lead))
}

# pQTL summary-stats file for the protein
resolve_pqtl <- function(protID, arm) {
  if (arm == "UKB") {
    hits <- list.files(UKB_PQTL, pattern = paste0("^", protID, "_.*\\.parquet$"), full.names = TRUE)
    if (length(hits)) hits[1] else NA_character_
  } else {
    mp <- DEC_MAP[EntrezGeneSymbol == protID]
    if (!nrow(mp)) return(NA_character_)
    f <- file.path(DEC_PQTL, mp$file[1])
    if (file.exists(f)) f else NA_character_
  }
}

# outcome GWAS file + case fraction (s) + N. FinnGen disease => cc; UKB exposure => quant.
resolve_gwas <- function(id) {
  if (startsWith(id, "finngen_R12_")) {
    f <- file.path(FINNGEN, paste0(id, ".gz"))
    ph <- sub("^finngen_R12_", "", id); man <- FG_MAN[phenocode == ph]
    nca <- if (nrow(man)) as.numeric(man$num_cases[1]) else NA_real_
    nco <- if (nrow(man)) as.numeric(man$num_controls[1]) else NA_real_
    list(path = if (file.exists(f)) f else NA_character_,
         is_cc = TRUE,
         s = if (is.finite(nca) && is.finite(nco) && (nca + nco) > 0) nca / (nca + nco) else NA_real_,
         n = if (is.finite(nca) && is.finite(nco)) nca + nco else NA_real_)
  } else {
    f <- file.path(EXPO_GWAS, id, paste0(id, ".regenie"))
    list(path = if (file.exists(f)) f else NA_character_,
         is_cc = FALSE, s = NA_real_, n = NA_real_)   # N inferred per-locus from sumstats
  }
}

# ---- build the manifest ----------------------------------------------------
man <- rbindlist(lapply(seq_len(nrow(pend)), function(i) {
  r <- pend[i]
  arm <- r$arm; protID <- r$protID; id <- r$disease_or_exposure_id
  lead <- resolve_lead(protID, arm)
  pq   <- resolve_pqtl(protID, arm)
  gw   <- resolve_gwas(id)

  # resolve_status precedence: missing inputs block the run; not_casecontrol is
  # an INFO flag for runnable P->E exposure loci (no FinnGen, coloc still runs).
  status <- "OK"
  if (is.na(pq))         status <- "missing_pqtl"
  else if (is.na(gw$path)) status <- "missing_gwas"
  else if (!lead$ok)    status <- "missing_lead"
  else if (!gw$is_cc)   status <- "not_casecontrol"

  data.table(
    arm = arm, protID = protID, disease_or_exposure_id = id, edge_dir = r$edge_dir,
    pqtl_path = pq, gwas_path = gw$path,
    lead_snp = lead$lead, chr = lead$chr, pos = lead$pos, window_kb = WINDOW_KB,
    outcome_type = if (gw$is_cc) "cc" else "quant",
    s_casefrac = gw$s, n_gwas = gw$n,
    resolve_status = status
  )
}), use.names = TRUE)

setorder(man, arm, edge_dir, protID, disease_or_exposure_id)
man[, row := .I]
setcolorder(man, "row")

out_fp <- file.path(out_dir, "coloc_manifest.tsv")
fwrite(man, out_fp, sep = "\t")

# ---- report ----------------------------------------------------------------
cat("\n================ COLOC MANIFEST ================\n")
cat("Rows:", nrow(man), "  ->  ", out_fp, "\n\n")
cat("resolve_status breakdown:\n")
brk <- man[, .N, by = resolve_status][order(-N)]
print(brk)
# runnable = anything not blocked by a missing input (OK + not_casecontrol)
runnable <- man[resolve_status %in% c("OK", "not_casecontrol")]
cat(sprintf("\nRunnable loci (OK + not_casecontrol): %d / %d\n", nrow(runnable), nrow(man)))
cat(sprintf("  - cc disease (OK):            %d\n", sum(man$resolve_status == "OK")))
cat(sprintf("  - quant exposure (not_cc):    %d\n", sum(man$resolve_status == "not_casecontrol")))
blocked <- man[!resolve_status %in% c("OK", "not_casecontrol")]
if (nrow(blocked)) {
  cat("\nBlocked loci:\n")
  print(blocked[, .(arm, protID, disease_or_exposure_id, resolve_status)])
}
cat("\nDONE.\n")
