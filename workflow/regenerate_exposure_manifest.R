#!/usr/bin/env Rscript

# regenerate_exposure_manifest.R
#
# Rebuilds config/exposure_sets/analysis_exposures.tsv from HEAP.rds so the
# exposure manifest stays in sync with the loader output.
#
# IMPORTANT: HEAP_loader.R does NOT write this file. It writes HEAP.rds and the
# *_audit_*.tsv sidecars. The manifest is a separate curated config that modules
# read (include == 1 AND miss_rate_prot_i0 < 0.20). Run THIS script after
# regenerating HEAP.rds whenever exposure features are added/removed/renamed
# (e.g., the income field 738 row and the air-pollution *_mean columns).
#
# For every variable in heap$E_baseline (instance 0) it computes:
#   miss_rate_prot_i0 = missingness within the proteomics instance-0 subset
# and sets, by the standard rule:
#   include = 1 if miss_rate_prot_i0 < MISS_THRESH (0.20) else 0
#   notes   = "miss_rate_prot_i0 > 0.20" when excluded by missingness
#
# Sync semantics vs the existing manifest:
#   - variables now in HEAP but not in the manifest  -> ADDED
#   - variables in the manifest but no longer in HEAP -> REMOVED (e.g. the old
#                                                        *_nearest_assessment rows)
#   - variable_type is PRESERVED from the existing manifest where present
#     (only newly-added variables are auto-classified), to avoid churn.
#   - MANUAL include/notes overrides are PRESERVED: if a variable's existing
#     include flag disagrees with the rule applied to its *existing* recorded
#     miss-rate, that decision (include + notes) is carried forward and flagged.
#   - SUPERSEDED air-pollution per-year columns are forced to include=0: NO2
#     (2005/06/07/10) and PM2.5 (2010) are excluded in favour of no2_mean /
#     pm2_5_mean (see SUPERSEDED_AIR below). PM10's per-year columns are kept
#     (no mean; 2007 vs 2010 disagree at r ~ 0.39).
#   - NON-RESPONSE one-hot indicator columns are forced to include=0:
#     "*_Prefer_not_to_answer" and "*_Do_not_know" (see NONRESPONSE_PATTERNS).
#     "None of the above" is intentionally kept (meaningful for multi-select).
#
# Usage:
#   module load gcc/14.2.0 R/4.4.2
#   Rscript workflow/regenerate_exposure_manifest.R            # dry-run: prints diff only
#   Rscript workflow/regenerate_exposure_manifest.R --write    # writes the manifest

# ---------------------------------------------------------------------------
# Bootstrap paths/helpers (same pattern as the other workflow scripts)
# ---------------------------------------------------------------------------
local({
  candidates <- c(
    Sys.getenv("HEAP_PATHS_FILE", unset = ""),
    file.path(Sys.getenv("HEAP_ROOT", unset = ""), "workflow", "00_paths.R"),
    file.path(getwd(), "workflow", "00_paths.R"),
    file.path(getwd(), "..", "workflow", "00_paths.R"),
    file.path(getwd(), "00_paths.R")
  )
  candidates <- candidates[nzchar(candidates)]
  hit <- candidates[file.exists(candidates)][1]
  if (is.na(hit)) stop("Could not locate workflow/00_paths.R")
  source(hit)
})

suppressPackageStartupMessages(library(stringr))

args      <- commandArgs(trailingOnly = TRUE)
DO_WRITE  <- "--write" %in% args
MISS_THRESH <- 0.20
MISS_NOTE   <- "miss_rate_prot_i0 > 0.20"
# Monomorphic (zero-variance) exposures carry no information: a constant predictor
# is uninformative in Modules 1/2/GREML and a constant phenotype is unmodellable in
# the exposure GWAS / Module 6. Variance is assessed on the proteomics-i0 analysis
# sample (the rows the modules actually use).
MONO_NOTE <- "monomorphic (zero variance in proteomics-i0 sample)"

# Air-pollution per-year fields that are SUPERSEDED by a mean-across-years
# feature: their raw per-year columns are forced to include=0 because the
# {pollutant}_mean feature is used instead (avoids using redundant repeated
# measures). Mirrors cfg$pollution_year_map in HEAP_loader.R. The mean feature
# itself stays include=1 (by the standard rule). PM10 is intentionally NOT here:
# its 2007 vs 2010 surfaces disagree (r ~ 0.39), so they are kept as two
# separate include=1 features. The supersede rule is gated on the mean column
# actually existing in HEAP (so the raw is never dropped without a replacement).
AIR_CATEGORY   <- "Residential_Air_Pollution"
SUPERSEDED_AIR <- list(
  no2_mean   = c(24016L, 24017L, 24018L, 24003L),  # NO2 2005/06/07/10 -> no2_mean
  pm2_5_mean = c(24006L)                           # PM2.5 2010        -> pm2_5_mean
)

# One-hot NON-RESPONSE indicator columns are not informative exposures, so they
# are forced to include=0. Matched (case-insensitively) as a suffix on the
# canonicalized variable name; "[._]" tolerates either separator. "None of the
# above" is deliberately EXCLUDED from this list: for multi-select fields it is a
# meaningful answer (= none of the listed options).
NONRESPONSE_PATTERNS <- c(
  "prefer[._]not[._]to[._]answer$",   # "..._Prefer_not_to_answer"
  "do[._]not[._]know$"                # "..._Do_not_know"
)
NONRESPONSE_RE   <- paste(NONRESPONSE_PATTERNS, collapse = "|")
NONRESPONSE_NOTE <- "non-response indicator (one-hot); not an informative exposure"
is_nonresponse   <- function(v) grepl(NONRESPONSE_RE, v, ignore.case = TRUE)

ts <- function(...) message(format(Sys.time(), "[%H:%M:%S] "), ...)

HEAP_RDS    <- heap_loader_rds
MANIFEST    <- heap_analysis_config()

if (!file.exists(HEAP_RDS))
  stop("HEAP.rds not found at: ", HEAP_RDS,
       "\nRun HEAP_loader.R first (or set HEAP_LOADER_RDS).")

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

# Variable-type classification, consistent with exposure_coding_check.R.
classify_var <- function(x) {
  if (is.character(x) || is.factor(x)) return("categorical")
  ux <- sort(unique(stats::na.omit(as.numeric(x))))
  if (length(ux) == 0L)                       return("continuous")
  if (all(ux %in% c(0, 1)))                   return("binary")
  if (length(ux) <= 10L && all(ux == floor(ux))) return("ordinal")
  "continuous"
}

# UKB field id embedded in a canonicalized variable name (_f1558_ -> 1558), else NA.
extract_field_id <- function(v) {
  m <- regmatches(v, regexpr("_f([0-9]+)_", v))
  if (length(m) == 0L) return(NA_integer_)
  as.integer(gsub("[^0-9]", "", m))
}

# ---------------------------------------------------------------------------
# Load inputs
# ---------------------------------------------------------------------------
ts("Loading HEAP.rds: ", HEAP_RDS)
heap <- readRDS(HEAP_RDS)

prot_eids <- unique(heap$prot_baseline$eid)
ts("Proteomics instance-0 participants: ", length(prot_eids))

# Resolve which air-pollution per-year field IDs are superseded by a *_mean
# feature, gated on that mean column existing in HEAP. -> map "field_id" -> mean.
air_vars <- names(heap$E_baseline[[AIR_CATEGORY]])
superseded_fid_to_mean <- character(0)
for (mean_feat in names(SUPERSEDED_AIR)) {
  if (mean_feat %in% air_vars) {
    superseded_fid_to_mean[as.character(SUPERSEDED_AIR[[mean_feat]])] <- mean_feat
  } else {
    ts("NOTE: ", mean_feat, " absent from ", AIR_CATEGORY,
       "; its per-year fields will NOT be superseded (kept by the standard rule).")
  }
}
is_superseded <- function(cat_name, fid) {
  !is.na(fid) && identical(cat_name, AIR_CATEGORY) &&
    as.character(fid) %in% names(superseded_fid_to_mean)
}
if (length(superseded_fid_to_mean))
  ts("Superseded air per-year fields (-> include=0): ",
     paste(names(superseded_fid_to_mean), collapse = ", "))

old <- if (file.exists(MANIFEST)) {
  utils::read.delim(MANIFEST, stringsAsFactors = FALSE, colClasses = "character")
} else {
  ts("No existing manifest; creating fresh.")
  data.frame(category = character(), field_id = character(), variable = character(),
             variable_type = character(), miss_rate_prot_i0 = character(),
             include = character(), notes = character(), stringsAsFactors = FALSE)
}

# Variables whose existing include disagrees with the rule on their OWN recorded
# miss-rate are treated as manual overrides and preserved verbatim -- EXCEPT
# rule-driven exclusions (superseded air fields, non-response indicators), so
# those rules stay authoritative and idempotent rather than being mistaken for a
# manual edit.
old_superseded <- if (nrow(old) == 0L) logical(0) else
  mapply(is_superseded, old$category,
         vapply(old$variable, extract_field_id, integer(1)), USE.NAMES = FALSE)
old_nonresponse <- if (nrow(old) == 0L) logical(0) else is_nonresponse(old$variable)
old_rule       <- ifelse(suppressWarnings(as.numeric(old$miss_rate_prot_i0)) < MISS_THRESH, 1L, 0L)
override_mask  <- !is.na(old_rule) & as.integer(old$include) != old_rule &
                  !old_superseded & !old_nonresponse
overrides      <- setNames(
  Map(function(inc, nt) list(include = inc, notes = nt),
      old$include[override_mask], old$notes[override_mask]),
  old$variable[override_mask]
)
old_type       <- setNames(old$variable_type, old$variable)
if (length(overrides))
  ts("Preserving ", length(overrides), " manual include/notes override(s).")

# ---------------------------------------------------------------------------
# Recompute one row per E_baseline variable on the proteomics-i0 subset
# ---------------------------------------------------------------------------
rows <- list()
for (cat_name in names(heap$E_baseline)) {
  edf <- heap$E_baseline[[cat_name]]
  if (is.null(edf) || nrow(edf) == 0L) next
  sub   <- edf[edf$eid %in% prot_eids, , drop = FALSE]
  vcols <- setdiff(names(edf), "eid")

  for (v in vcols) {
    miss <- round(mean(is.na(sub[[v]])), 4)
    fid  <- extract_field_id(v)
    vtype <- if (!is.na(old_type[v])) old_type[[v]] else classify_var(edf[[v]])
    n_distinct <- length(unique(stats::na.omit(sub[[v]])))   # variance on analysis sample

    if (n_distinct < 2L) {                    # monomorphic / zero variance -> excluded
      inc  <- 0L                              # highest precedence: a constant is useless
      note <- MONO_NOTE                       # even if a manual override said include=1
    } else if (is_superseded(cat_name, fid)) {  # air per-year col -> excluded for the *_mean
      inc  <- 0L
      note <- paste0("superseded by ", superseded_fid_to_mean[[as.character(fid)]],
                     " (mean-across-years feature used instead)")
    } else if (is_nonresponse(v)) {           # one-hot non-response indicator -> excluded
      inc  <- 0L
      note <- NONRESPONSE_NOTE
    } else if (!is.null(overrides[[v]])) {    # preserve manual decision
      inc  <- as.integer(overrides[[v]]$include)
      note <- overrides[[v]]$notes
    } else {                                  # standard rule
      inc  <- if (miss < MISS_THRESH) 1L else 0L
      note <- if (inc == 0L) MISS_NOTE else ""
    }

    rows[[length(rows) + 1L]] <- data.frame(
      category = cat_name, field_id = if (is.na(fid)) NA_character_ else as.character(fid),
      variable = v, variable_type = vtype,
      miss_rate_prot_i0 = formatC(miss, format = "f", digits = 4),
      include = inc, notes = note, stringsAsFactors = FALSE
    )
  }
}
new <- do.call(rbind, rows)

# Order: category, then numeric field_id (derived NA-id rows last), then variable.
fid_num <- suppressWarnings(as.integer(new$field_id))
new <- new[order(new$category, is.na(fid_num), fid_num, new$variable), ]
rownames(new) <- NULL

# ---------------------------------------------------------------------------
# Diff report
# ---------------------------------------------------------------------------
added   <- setdiff(new$variable, old$variable)
removed <- setdiff(old$variable, new$variable)
common  <- intersect(new$variable, old$variable)
inc_old <- setNames(as.integer(old$include), old$variable)
inc_new <- setNames(new$include, new$variable)
flipped <- common[inc_old[common] != inc_new[common]]

cat("\n================ exposure manifest diff ================\n")
cat(sprintf("HEAP.rds : %s\n", HEAP_RDS))
cat(sprintf("manifest : %s\n", MANIFEST))
cat(sprintf("rows: %d existing -> %d regenerated\n", nrow(old), nrow(new)))
cat(sprintf("\nADDED (%d):\n", length(added)))
if (length(added)) cat(paste0("  + ", added, collapse = "\n"), "\n") else cat("  (none)\n")
cat(sprintf("\nREMOVED (%d):\n", length(removed)))
if (length(removed)) cat(paste0("  - ", removed, collapse = "\n"), "\n") else cat("  (none)\n")
cat(sprintf("\nINCLUDE FLIPS (%d):\n", length(flipped)))
if (length(flipped)) {
  for (v in flipped) cat(sprintf("  ~ %s: include %d -> %d\n", v, inc_old[v], inc_new[v]))
} else cat("  (none)\n")
cat(sprintf("\ninclude=1: %d | include=0: %d\n", sum(new$include == 1L), sum(new$include == 0L)))

# ---------------------------------------------------------------------------
# Write (only with --write)
# ---------------------------------------------------------------------------
if (DO_WRITE) {
  dir.create(dirname(MANIFEST), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(new, MANIFEST, sep = "\t", row.names = FALSE, quote = FALSE)
  ts("WROTE manifest: ", MANIFEST, " (", nrow(new), " rows)")
} else {
  cat("\n[dry-run] No file written. Re-run with --write to update the manifest.\n")
}
