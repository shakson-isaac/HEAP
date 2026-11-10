#!/usr/bin/env Rscript
# ============================================================================
# summarize_module3_mediation.R
# ----------------------------------------------------------------------------
# Rebuilds the Module 3 (GEM) mediation supplementary outputs FROM THE CANONICAL
# RUN DIRECTORIES. Replaces scripts/support/build_mediation_supp_tables.R, which
# read its main input from a temporary session scratch directory under /tmp that no
# longer exists -- so it could not be re-run and its outputs were orphaned.
#
# Everything here derives from the IGLOO output tree only.
#
# BASE specification (the main reporting spec, per SUPPLEMENT_DATA_STANDARD §2.2):
#   M3_base_lasso_primary/base/lasso/primary_total
#       predictors: PXS_total (exposure_total), G_raw (genetic_total)
#   M3_base_lasso_partitioned/base/lasso/partitioned_categories
#       predictors: Gcis_raw, Gtrans_raw, PXS_<13 exposure categories>
#
# Each (protein, disease, predictor) has an NDE and an NIE row; they are pivoted
# into one link. TE(logHR) = NIE + NDE; PM = logNIE/logTE when sign-consistent;
# q = BH FDR of the NIE p within each driver class. A link is a significant
# mediator when q < 0.05 and 0 <= PM <= 1.
#
# Outputs
#   supp_deposit/med_exposure_all.tsv   all exposure/PXS-category links   (Data)
#   supp_deposit/med_genetic_all.tsv    all genetic links (cis/trans/total) (Data)
#   manuscript_stats/module3/med_disease_summary.tsv   per-disease   (S12, Tier 1)
#   manuscript_stats/module3/med_protein_summary.tsv   per-protein
#   manuscript_stats/module3/med_dominance.tsv         exposure- vs genetic-dominant
#
#   Run:  Rscript scripts/analysis_summaries/summarize_module3_mediation.R
#         Rscript scripts/analysis_summaries/summarize_module3_mediation.R base_bmi
#
#   No argument = BASE, writing the deposit filenames above. Those are cited
#   supplementary data, so that path is unchanged.
#
#   An argument names one of base_bmi / base_clinical / base_draw /
#   base_exclprev and writes the same three Data files under
#     supp_deposit/med_specs_partitioned/<spec>/
#   with the per-disease and per-protein summaries under
#     docs/manuscript_stats/module3/specs/<spec>/
#   Nothing else differs: the same FDR rule (BH within predictor_class), the same
#   PM clamping, the same estimable/structural-zero handling. That is the point --
#   a specification is only comparable to base if it was reduced the same way.
# ============================================================================
suppressPackageStartupMessages(library(data.table))
setDTthreads(0)

local({
  cand <- c(Sys.getenv("HEAP_PATHS_FILE", unset = ""),
            "/n/groups/patel/shakson_ukb/HEAP/workflow/00_paths.R")
  hit <- cand[nzchar(cand) & file.exists(cand)][1]; if (!is.na(hit)) source(hit)
})
## ---- which specification -------------------------------------------------
## No argument reproduces the previous behaviour exactly: the BASE specification,
## written to the same deposit filenames it has always used. Those files are
## cited supplementary data, so the default path had to stay byte-identical and
## nothing below branches on SPEC except the run directories and the destination.
##
## With an argument, the same code reads that specification's runs and writes
## under supp_deposit/med_specs_partitioned/<spec>/. The FDR rule is untouched --
## BH within predictor_class -- because computing it any other way for one
## specification would make it silently incomparable to the others, which is the
## whole reason these are produced by this script rather than ad hoc.
##
## exclprev is a SAMPLE specification: its covariate directory is `base`, not
## `base_exclprev`. Read from the tree, not guessed from the experiment name.
SPEC_DIRS <- list(
  base          = list(exp = "M3_base_lasso",          cov = "base"),
  base_bmi      = list(exp = "M3_base_bmi_lasso",      cov = "base_bmi"),
  base_clinical = list(exp = "M3_base_clinical_lasso", cov = "base_clinical"),
  base_draw     = list(exp = "M3_base_draw_lasso",     cov = "base_draw"),
  base_exclprev = list(exp = "M3_base_exclprev_lasso", cov = "base")
)
.args <- commandArgs(trailingOnly = TRUE)
SPEC  <- if (length(.args) && nzchar(.args[1])) .args[1] else "base"
if (!SPEC %in% names(SPEC_DIRS))
  stop("unknown specification '", SPEC, "'. One of: ", paste(names(SPEC_DIRS), collapse = ", "))
SD <- SPEC_DIRS[[SPEC]]

M3   <- heap_project_output("module3")
DEP  <- heap_project_output("supp_deposit")
if (!identical(SPEC, "base")) DEP <- file.path(DEP, "med_specs_partitioned", SPEC)
STAT <- file.path(heap_root, "docs", "manuscript_stats", "module3")
# med_disease_summary.tsv is Supplementary Table S12. A sensitivity run must not
# land on top of it, so every non-base specification gets its own subdirectory.
if (!identical(SPEC, "base")) STAT <- file.path(STAT, "specs", SPEC)
dir.create(DEP,  recursive = TRUE, showWarnings = FALSE)
dir.create(STAT, recursive = TRUE, showWarnings = FALSE)

KEEP <- c("protID","DZ_ID","predictor","predictor_class","effect_type",
          "effect_logHR","effect_HR","delta_se","delta_l95_HR","delta_u95_HR",
          "delta_p","protein_HR","protein_p","n","n_cases","instrument_present")

## DZ_ID -> ICD-10 code + readable disease name
parse_dz <- function(x) {
  s <- sub("^age_", "", x); s <- sub("_f[0-9]+_[0-9]+_[0-9]+$", "", s)
  code <- toupper(sub("^([a-z][0-9]+)_.*$", "\\1", s))
  nm   <- sub("^[a-z][0-9]+_", "", s); nm <- sub("^first_reported_", "", nm)
  list(icd = code, disease = gsub("_", " ", nm))
}

## read one run directory -> one row per (protein, disease, predictor)
read_run <- function(dir) {
  f <- list.files(dir, "^MDres_.*\\.txt$", full.names = TRUE)
  if (!length(f)) stop("no MDres files under ", dir)
  message("  ", basename(dirname(dirname(dirname(dir)))), "/", basename(dir),
          ": ", length(f), " files")
  acc <- vector("list", length(f))
  for (i in seq_along(f)) {
    d <- fread(f[i], select = KEEP, showProgress = FALSE)
    nie <- d[effect_type == "NIE",
             .(protID, DZ_ID, predictor, predictor_class, n, n_cases,
               protein_HR, protein_p, driver_present = instrument_present,
               nie_logHR = effect_logHR, nie_HR = effect_HR, nie_se = delta_se,
               nie_l95 = delta_l95_HR, nie_u95 = delta_u95_HR, nie_p = delta_p)]
    # the NDE row carries its OWN delta-method CI and p -- keep them, not just
    # the point estimate
    nde <- d[effect_type == "NDE",
             .(protID, DZ_ID, predictor, nde_logHR = effect_logHR, nde_HR = effect_HR,
               nde_l95 = delta_l95_HR, nde_u95 = delta_u95_HR, nde_p = delta_p)]
    acc[[i]] <- merge(nie, nde, by = c("protID", "DZ_ID", "predictor"))
    if (i %% 100 == 0) message("    ", i, "/", length(f))
  }
  rbindlist(acc, fill = TRUE)
}

# assembling the links means reading ~6 GB, so cache it: a downstream change
# should not force a re-read. Delete the cache to force a rebuild.
# cache lives in intermediate/, NOT in supp_deposit/ -- that directory holds only
# files destined for the public deposit and gets uploaded wholesale.
CACHE <- heap_project_intermediate(sprintf(".med_links_%s_cache.rds", SPEC))
if (file.exists(CACHE) && !nzchar(Sys.getenv("MED_REBUILD"))) {
  message("Reading cached links (", CACHE, ") -- set MED_REBUILD=1 to re-read")
  L <- readRDS(CACHE)
} else {
  message("Reading ", toupper(SPEC), " mediation runs ...")
  d_prim <- file.path(M3, paste0(SD$exp, "_primary"),     SD$cov, "lasso", "primary_total")
  d_part <- file.path(M3, paste0(SD$exp, "_partitioned"), SD$cov, "lasso", "partitioned_categories")
  message("  primary     : ", d_prim)
  message("  partitioned : ", d_part)
  # Fail on a partial run rather than summarising it. A partitioned array that
  # lost tasks still leaves hundreds of MDres files, and the resulting deposit
  # would look complete while quietly missing proteins -- worse than an error.
  for (.d in c(d_prim, d_part)) {
    if (!dir.exists(.d)) stop("missing run directory: ", .d)
    .nf <- length(list.files(.d, "^MDres_.*\\.txt$"))
    if (.nf < 400) stop(sprintf("%s has %d of 400 MDres files -- run incomplete", .d, .nf))
  }
  prim <- read_run(d_prim)
  part <- read_run(d_part)
  L <- rbindlist(list(prim, part), fill = TRUE); rm(prim, part); invisible(gc())
  # one stray unparsable value in any of the 800 files makes rbindlist coerce a
  # whole column to character -- force the numeric columns back.
  numcols <- c("nie_logHR","nie_HR","nie_l95","nie_u95","nie_p",
               "nde_logHR","nde_HR","protein_HR","protein_p","n","n_cases")
  for (cc in intersect(numcols, names(L))) if (!is.numeric(L[[cc]])) {
    bad <- sum(is.na(suppressWarnings(as.numeric(L[[cc]]))) & !is.na(L[[cc]]))
    message("  coercing '", cc, "' to numeric (", bad, " unparsable -> NA)")
    set(L, j = cc, value = suppressWarnings(as.numeric(L[[cc]])))
  }
  saveRDS(L, CACHE)
}

## total effect, proportion mediated, BH-q within driver class
L[, total_logHR := nie_logHR + nde_logHR]
L[, total_HR    := exp(total_logHR)]
# Semantics replicated EXACTLY from the retired build_mediation_supp_tables.R so
# that the only thing changing here is provenance, not the numbers. Note its
# header comment said "sign-consistent (0<=PM<=1)" but the code CLAMPS PM>1 to 1
# rather than excluding it, and drops only PM<0. Keep the code's behaviour.
L[, prop_mediated := fifelse(sign(nie_logHR) == sign(total_logHR) & total_logHR != 0,
                             nie_logHR / total_logHR, NA_real_)]
L[prop_mediated < 0, prop_mediated := NA_real_]
# PM > 1 is clamped to 1 (inherited behaviour, kept for continuity), but flagged
# so a clamped 1.0 is distinguishable from a genuine one. A ratio above 1 means
# |NIE| exceeded |TE|, which signals an unstable total effect rather than
# complete mediation.
L[, pm_clamped := !is.na(prop_mediated) & prop_mediated > 1]
L[prop_mediated > 1, prop_mediated := 1]
L[, q := NA_real_]
L[!is.na(nie_p), q := p.adjust(nie_p, "BH"), by = predictor_class]
L[, sig := !is.na(q) & q < 0.05 & !is.na(prop_mediated)]

# ---- distinguish the two kinds of exactly-zero effect ----------------------
# Module3.R emits logHR = 0 / HR = 1 in two situations that mean different things:
#  (a) instrument_present == FALSE -- the protein has no cis/trans variants, so
#      there is genuinely no genetic pathway. Zero is CORRECT and is kept
#      (mirrors Module 1's "R2 = 0 + present = FALSE" convention).
#  (b) zscore_numeric() collapsed a zero-variance predictor to a constant 0
#      within THAT disease's analysis sample (Module3.R:118-130). Nothing was
#      estimated, so reporting 0/1 asserts a null we never measured -> set NA.
# Both already have p = NA, hence q = NA and sig = FALSE, so this changes no
# count; it only stops an unestimable cell being read as a precise zero.
L[, estimable := !(is.na(nie_se) & driver_present %in% c(TRUE, "TRUE"))]
unest <- which(!L$estimable)
for (cc in c("nie_HR","nie_l95","nie_u95","nde_HR","nde_l95","nde_u95",
             "total_HR","prop_mediated","nie_logHR"))
  if (cc %in% names(L)) set(L, i = unest, j = cc, value = NA_real_)
message(sprintf("unestimable links set to NA: %s (%.2f%%) | structural zeros kept: %s",
                format(length(unest), big.mark = ","), 100 * length(unest) / nrow(L),
                format(sum(!(L$driver_present %in% c(TRUE, "TRUE"))), big.mark = ",")))
dz <- parse_dz(L$DZ_ID); L[, `:=`(icd = dz$icd, disease = dz$disease)]
message(sprintf("links: %s | significant: %s", format(nrow(L), big.mark = ","),
                format(sum(L$sig), big.mark = ",")))

## ---- Tier-2 Data: all links, split exposure vs genetic ---------------------
# Column order follows SUPPLEMENT_DATA_STANDARD §3:
#   identifiers | n | estimates (each with its own CI) | inference | flags
COLS <- function(lab) c(
  lab, "protID", "icd", "disease", "DZ_ID",          # identifiers
  "n_cases",                                          # sample
  "nie_HR", "nie_l95", "nie_u95",                     # indirect effect + CI
  "nde_HR", "nde_l95", "nde_u95",                     # direct effect + CI
  "total_HR", "prop_mediated",                        # total, proportion mediated
  "protein_HR", "protein_p",                          # protein -> disease
  "nie_p", "nde_p", "q",                              # inference
  "sig", "estimable", "driver_present", "pm_clamped") # flags
rnd <- function(d) { for (c in c("nie_HR","nie_l95","nie_u95","nde_HR","total_HR",
                                 "prop_mediated","protein_HR"))
                       if (c %in% names(d)) set(d, j = c, value = signif(d[[c]], 4)); d }

# The aggregate PXS_total and the 13 per-category scores are shipped SEPARATELY.
# Proportion mediated is always relative to its own driver's total effect, so a
# PM against one category is not on the same scale as a PM against the whole
# exposome (median 0.159 vs 0.127 among significant links). Splitting the files
# makes pooling them impossible rather than merely discouraged.
cols_ex <- COLS("exposure_category")

tot <- L[predictor_class == "exposure_total"]
tot[, exposure_category := "PXS_total"]
fwrite(rnd(tot)[, ..cols_ex], file.path(DEP, "med_exposure_total.tsv"), sep = "\t")
message("  med_exposure_total.tsv       ", format(nrow(tot), big.mark = ","),
        " links (aggregate PXS_total)")

cat_ <- L[predictor_class == "exposure_category"]
cat_[, exposure_category := sub("^PXS_", "", predictor)]
fwrite(rnd(cat_)[, ..cols_ex], file.path(DEP, "med_exposure_categories.tsv"), sep = "\t")
message("  med_exposure_categories.tsv  ", format(nrow(cat_), big.mark = ","),
        " links (", uniqueN(cat_$exposure_category), " categories)")

# genetic: cis and trans are components, PGS_total is the combined score
gn <- L[predictor_class %in% c("genetic_total","genetic_cis","genetic_trans")]
gn[, component := fifelse(predictor_class == "genetic_total", "PGS_total",
                          sub("^genetic_", "", predictor_class))]
cols_gn <- COLS("component")
fwrite(rnd(gn)[, ..cols_gn], file.path(DEP, "med_genetic_all.tsv"), sep = "\t")
message("  med_genetic_all.tsv          ", format(nrow(gn), big.mark = ","), " links")

## ---- Tier-1 summaries (base spec) ------------------------------------------
# The original counts exposure mediation from `exposure_category` ONLY (not the
# exposure_total driver), and emits a row for EVERY disease present in the link
# set -- including diseases with no significant mediator, which is why it has 181
# rows rather than only the diseases that reach significance.
sigL <- L[sig == TRUE]
EXPCLS <- "exposure_category"
GENCLS <- c("genetic_total", "genetic_cis", "genetic_trans")
allDz  <- unique(L[, .(DZ_ID, icd, disease)])
agg <- sigL[, .(n_exposure_mediated_proteins = uniqueN(protID[predictor_class == EXPCLS]),
                n_genetic_mediated_proteins  = uniqueN(protID[predictor_class %in% GENCLS]),
                n_cis_proteins   = uniqueN(protID[predictor_class == "genetic_cis"]),
                n_trans_proteins = uniqueN(protID[predictor_class == "genetic_trans"]),
                n_both_proteins  = uniqueN(protID[predictor_class == "genetic_total"]),
                n_shared_proteins = length(intersect(protID[predictor_class == EXPCLS],
                                                     protID[predictor_class %in% GENCLS]))),
            by = DZ_ID]
dis <- merge(allDz, agg, by = "DZ_ID", all.x = TRUE)
for (cc in setdiff(names(dis), c("DZ_ID", "icd", "disease")))
  set(dis, i = which(is.na(dis[[cc]])), j = cc, value = 0L)
setorder(dis, -n_exposure_mediated_proteins, -n_genetic_mediated_proteins)
fwrite(dis, file.path(STAT, "med_disease_summary.tsv"), sep = "\t")
message("  med_disease_summary.tsv  ", nrow(dis), " diseases")

prot <- sigL[, .(n_diseases = uniqueN(DZ_ID),
                 n_exposure_categories = uniqueN(predictor[predictor_class == "exposure_category"]),
                 strongest_exposure_NIE_HR = signif(nie_HR[which.max(abs(nie_logHR * (predictor_class %in%
                   c("exposure_total","exposure_category"))))], 4),
                 strongest_genetic_NIE_HR = signif(nie_HR[which.max(abs(nie_logHR * (predictor_class %in%
                   c("genetic_total","genetic_cis","genetic_trans"))))], 4)), by = protID]
setorder(prot, -n_diseases)
fwrite(prot, file.path(STAT, "med_protein_summary.tsv"), sep = "\t")
message("  med_protein_summary.tsv  ", nrow(prot), " proteins")

dom <- sigL[, .(max_exposure_logNIE = suppressWarnings(max(abs(nie_logHR[predictor_class %in%
                  c("exposure_total","exposure_category")]), na.rm = TRUE)),
                max_genetic_logNIE  = suppressWarnings(max(abs(nie_logHR[predictor_class %in%
                  c("genetic_total","genetic_cis","genetic_trans")]), na.rm = TRUE))), by = protID]
dom <- dom[is.finite(max_exposure_logNIE) & is.finite(max_genetic_logNIE)]
dom[, gem_balance := signif(max_genetic_logNIE - max_exposure_logNIE, 4)]
dom[, dominance := fifelse(gem_balance > 0, "genetic-dominant", "exposure-dominant")]
setorder(dom, gem_balance)
fwrite(dom, file.path(STAT, "med_dominance.tsv"), sep = "\t")
message("  med_dominance.tsv  ", nrow(dom), " proteins")
