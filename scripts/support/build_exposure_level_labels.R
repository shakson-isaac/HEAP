#!/usr/bin/env Rscript
# ============================================================================
# build_exposure_level_labels.R
# ----------------------------------------------------------------------------
# A codebook the project did not have: for every ordinal exposure, the readable
# meaning of each dummy-term suffix that Module 2 emits.
#
# WHY THIS IS NEEDED. The loader encodes ordinals to 0-indexed integers and
# discards the label strings, so the levels are NOT recoverable from HEAP.rds --
# an earlier attempt to read them back out of the 9 GB object returned nothing.
# Figures were therefore hardcoding label maps by hand (fig_dose_response.R
# carried one for alcohol), which is unverifiable and goes stale silently.
#
# THE MAPPING, derived and then VALIDATED against that hardcoded alcohol map:
#   1. UKB codings are ranked ascending -> 0-indexed  (`scaled`)
#   2. fields in HEAP_loader.R `ordinal_reversals` are flipped, new = max - x,
#      so 0 is always the LOWEST exposure                (`recoded`)
#   3. Module 2 dummy-codes the factor and appends the 1-based level index, so
#      the term suffix is recoded + 1; level 1 is the reference and is not emitted
# Applying this to field 1558 reproduces the hand-written alcohol labels exactly
# (suffix 2 = "Special occasions only" ... suffix 6 = "Daily or almost daily"),
# which is the check that the rule is right rather than merely plausible.
#
# Fields whose UKB coding holds only negative special values (-1 "Do not know",
# -3 "Prefer not to answer") are plain counts -- the days-per-week activity
# fields -- so their level IS the number, and they are labelled as such.
#
# Output: docs/manuscript_stats/exposure_level_labels.tsv
# ============================================================================
suppressPackageStartupMessages(library(data.table))
local({
  cand <- c(Sys.getenv("HEAP_PATHS_FILE", unset = ""),
            "/n/groups/patel/shakson_ukb/HEAP/workflow/00_paths.R")
  hit <- cand[nzchar(cand) & file.exists(cand)][1]; if (!is.na(hit)) source(hit)
})

COD_F <- "/n/groups/patel/IGLOO/UKB/RAW/codings/Codings.csv"
DIC_F <- "/n/groups/patel/uk_biobank/Data_Dictionary_Showcase.csv"
COD <- fread(COD_F, quote = "\"", showProgress = FALSE)
DIC <- fread(DIC_F, showProgress = FALSE)

# from HEAP_loader.R `ordinal_reversals` -- keep in sync if that list grows
REV <- c(`1558` = 5L, `1628` = 2L, `1249` = 3L, `3506` = 2L)

field_levels <- function(fid) {
  row <- DIC[FieldID == fid]
  if (!nrow(row)) return(NULL)
  cid <- row$Coding[1]; nm <- row$Field[1]
  cd  <- COD[Coding == cid & suppressWarnings(as.integer(Value)) >= 0]
  if (!nrow(cd)) {                       # plain count field (days per week)
    return(data.table(field_id = fid, field = nm,
                      suffix = 2:8, recoded = 1:7,
                      level_label = paste(1:7, ifelse(1:7 == 1, "day", "days"))))
  }
  cd <- cd[order(as.integer(Value))]
  cd[, scaled := seq_len(.N) - 1L]
  mx <- max(cd$scaled)
  cd[, recoded := if (as.character(fid) %in% names(REV)) mx - scaled else scaled]
  data.table(field_id = fid, field = nm, suffix = cd$recoded + 1L,
             recoded = cd$recoded, level_label = cd$Meaning)[order(suffix)]
}

FIELDS <- as.integer(unique(na.omit(DIC$FieldID)))
ORD <- fread(file.path(heap_root, "config", "exposure_sets", "analysis_exposures.tsv"))
ORD <- ORD[include == 1 & variable_type == "ordinal"]
ORD[, fid := suppressWarnings(as.integer(sub("^.*_f([0-9]+)_.*$", "\\1", variable)))]
use <- sort(unique(na.omit(ORD$fid)))
message(sprintf("ordinal analysis exposures with a field id: %d", length(use)))

out <- rbindlist(lapply(use, field_levels), fill = TRUE)
out <- out[!is.na(level_label)]
f <- file.path(heap_root, "docs", "manuscript_stats", "exposure_level_labels.tsv")
fwrite(out, f, sep = "\t")
message(sprintf("wrote %s: %d fields, %d level rows", f, uniqueN(out$field_id), nrow(out)))

# the validation that makes this trustworthy
alc <- out[field_id == 1558][order(suffix)]
EXPECT <- c(`2` = "Special occasions only", `6` = "Daily or almost daily")
ok <- all(sapply(names(EXPECT), function(k)
  identical(alc[suffix == as.integer(k)]$level_label, unname(EXPECT[k]))))
message(sprintf("alcohol cross-check against the hand-written map: %s", if (ok) "PASS" else "FAIL"))
if (!ok) stop("level derivation no longer reproduces the verified alcohol labels")
