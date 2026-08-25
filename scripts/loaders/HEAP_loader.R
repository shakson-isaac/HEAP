#!/usr/bin/env Rscript
# HEAP_loader.R -- Unified UK Biobank data loader for the HEAP workflow.
#
# Replaces three separate loaders:
#   loaderProt_PGS_PXSv2.R             (baseline proteomics + E + covariates)
#   loaderProt_PGS_PXS_mediationv2.R   (disease T2E outcomes)
#   loaderProt_PGS_PXS_longitudinal.R  (multi-visit proteomics + E + covariates)
#
# Output: HEAP.rds saved to heap_loader_rds (IGLOO-rooted canonical path).
#   Default: /n/groups/patel/IGLOO/UKB/HEAP/intermediate/HEAP.rds
#   Override: set HEAP_LOADER_RDS env var (scratch only for temporary testing).
#
# Structure of HEAP.rds:
#   $meta              -- provenance and field catalog
#   $protIDs           -- protein identifiers (shared)
#   $Elist_names       -- exposure category names (shared)
#   $Eid_cat           -- Eid -> Category lookup (shared)
#   $ordinalIDs        -- ordinal exposure variable names (shared)
#   $covars_list       -- full covariate column names (shared)
#   $prot_long         -- (eid, instance) x proteins data.frame
#   $E_long            -- list of (eid, instance) x exposure data.frames by category
#   $covars_long       -- (eid, instance) x covariates data.frame (includes assessment timing)
#   $split_df          -- eid-level train/holdout assignments
#   $prot_baseline     -- eid x proteins (instance 0 only, no instance col)
#   $E_baseline        -- list of eid x exposure data.frames (instance 0 only)
#   $covars_baseline   -- eid x covariates (instance 0 only, no instance col)
#   $disease           -- list(DZ_df, DZ_ids) first-occurrence T2E outcomes
#   $prevalent_disease -- eid x baseline prevalent-disease flags (major chronic
#                         disease binary + multimorbidity count); also merged into
#                         covars_long/covars_baseline. Cols listed in meta$prevalent_disease_cols
#   $audit             -- data quality tables written to sidecar TSVs
#
# Env vars:
#   HEAP_SKIP_DISEASE=TRUE   skip first-occurrence disease loading (dev runs)
#
# Downstream compatibility:
#   as_pxs_baseline(heap)     -> PXSconstruct-compatible list for Modules 1/2/3/5
#   as_pxs_longitudinal(heap) -> longitudinal-compatible list for Module 6
# These functions are defined at the bottom of this file and in 00_paths.R
# after HEAP.rds is loaded by a module script.

local({
  candidates <- c(
    Sys.getenv("HEAP_PATHS_FILE", unset = ""),
    file.path(Sys.getenv("HEAP_ROOT", unset = ""), "workflow", "00_paths.R"),
    file.path(getwd(), "workflow", "00_paths.R"),
    file.path(getwd(), "..", "workflow", "00_paths.R"),
    file.path(getwd(), "..", "..", "workflow", "00_paths.R"),
    file.path(getwd(), "..", "..", "..", "workflow", "00_paths.R")
  )
  candidates <- candidates[nzchar(candidates)]
  hit <- candidates[file.exists(candidates)][1]
  if (!is.na(hit)) source(hit)
})

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(stringr)
  library(lubridate)
  library(future)
  library(arrow)
})

############################################################
# 0) Configuration
############################################################

cfg <- list(
  projID             = 52887L,
  instances          = c(0L, 2L, 3L),
  canonical_instance = 0L,
  skip_disease       = isTRUE(as.logical(Sys.getenv("HEAP_SKIP_DISEASE", "FALSE"))),
  dataloader_dir     = heap_script("loaders"),

  paths = list(
    # Raw Olink + protein-id map: IGLOO RAW canonical, shared UKB project tree fallback.
    olink       = heap_raw_or_legacy(c("proteomics", "olink_data_52887.txt"),
                    "/n/groups/patel/uk_biobank/olink_22881_52887/olink_data_52887.txt"),
    protein_id  = heap_raw_or_legacy(c("proteomics", "protein_id_conv.txt"),
                    "/n/groups/patel/uk_biobank/project_22881_672185/protein_id_conv.txt"),
    # allpaths dictionary + data codings: IGLOO RAW canonical, legacy fallback.
    path_all    = heap_raw_or_legacy("allpaths_52887.txt",
                    legacy_ukb_path("Data", "Paths", "52887", "allpaths.txt")),
    codings     = heap_raw_or_legacy(c("codings", "Codings.csv"),
                    legacy_ukb_path("RScripts", "Extract_Raw", "Finalized", "Codings.csv")),
    out_rds     = heap_loader_rds
  ),

  # Lifestyle exposure categories loaded via load_lifestyle_category() (path IDs)
  lifestyle_path_ids = c(
    Alcohol        = 100051L,
    Diet_Weekly    = 100052L,
    Smoking        = 100058L,
    Exercise_MET   = 54L,
    Exercise_Freq  = 100054L,
    Internet_Usage = 100053L,
    Sleep          = 100057L,
    Sun_Exposure   = 100055L,
    Sexual_Factors = 100056L
  ),

  # UKBdict category IDs for field-based exposure loading (section 2)
  deprivation_category   = 76L,
  geographical_category  = 711L,
  geographical_exclude   = c(54L),   # field 54 = assessment centre (already a covariate)

  # Individual household income (field 738; Category 100066 Sociodemographics >
  # Household). NOT in deprivation category 76, so it is loaded explicitly and
  # merged into Deprivation_Indices. The raw data stores label strings, which are
  # encoded to a 0-indexed ordinal (0 = lowest bracket ... 4 = highest income).
  # "Do not know" / "Prefer not to answer" (and any unmapped value) -> NA.
  income_field     = 738L,
  income_label_map = c(
    "Less than 18,000"     = 0L,
    "18,000 to 30,999"     = 1L,
    "31,000 to 51,999"     = 2L,
    "52,000 to 100,000"    = 3L,
    "Greater than 100,000" = 4L
  ),

  # Vitamins loaded via specific field IDs (multi-array one-hot)
  vitamin_fields = c(6155L, 6179L),

  # Residential air pollution fields (section 3)
  # NO2 (multiple years), NOx, PM10, PM2.5 variants, traffic exposure
  air_pollution_fields = c(
    24016L, 24017L, 24018L, 24003L,  # NO2: 2005, 2006, 2007, 2010
    24004L,                           # NOx: 2010
    24019L, 24005L,                   # PM10: 2007, 2010
    24006L, 24007L, 24008L,           # PM2.5, PM2.5 absorbance, PM2.5-10: 2010
    24009L, 24011L, 24013L, 24015L   # Traffic intensity / road load
  ),

  # Residential noise pollution fields (section 3)
  noise_pollution_fields = c(
    24020L,  # Average daytime sound level
    24021L,  # Average evening sound level
    24022L   # Average night-time sound level
  ),

  # Multi-year pollution field -> calendar year map (section 4)
  # Groups the per-year UKB fields whose estimates are collapsed into a single
  # per-individual MEAN across model years (na.rm = TRUE), NOT a temporal
  # "nearest-assessment" pick. The UKB air pollution fields are land-use-
  # regression spatial surfaces (not annual personal measurements), so averaging
  # across years yields a robust spatial exposure with no past/future temporal-
  # matching assumption. Only pollutants whose model years actually agree belong
  # here:
  #   no2   -- 4 years (2005/06/07/10), cross-year r >= 0.79 -> averaged.
  #   pm2_5 -- single 2010 surface, so its "mean" is just that estimate.
  # PM10 is intentionally EXCLUDED: its 2007 and 2010 surfaces only correlate
  # r ~ 0.39 (different model vintages, means 22.1 vs 16.2 ug/m3), so they are
  # kept as two SEPARATE features (raw columns) rather than averaged.
  # The per-year raw columns are retained for all pollutants so cross-year
  # stability can be verified (see exposure_coding_check.R).
  pollution_year_map = list(
    no2   = c("24016" = 2005L, "24017" = 2006L, "24018" = 2007L, "24003" = 2010L),
    pm2_5 = c("24006" = 2010L)
  ),

  # Covariates loaded via specific field IDs.
  # Field 53 = date of attending assessment centre (needed for timing derivation).
  covar_fields      = c(21003L, 31L, 23104L, 74L, 54L, 22009L, 53L),
  medication_fields = c(6177L, 6153L),

  # First-occurrence disease categories
  disease_categories = c(
    2401L, 2403L, 2404L, 2405L,
    2406L, 2407L, 2408L, 2409L,
    2410L, 2411L, 2412L, 2413L,
    2414L, 2415L, 2416L, 2417L
  ),
  disease_min_cases = 100L,

  # Baseline additional fields for disease T2E
  baseline_fields_t2e = c(40007L, 40000L, 54L, 21003L, 53L, 34L, 52L, 191L, 190L),

  # Prevalent disease at baseline (covariate derivation; see derive_prevalent_disease).
  # A first-occurrence disease is "prevalent at baseline" for a participant when its
  # first-occurrence age <= the participant's age at baseline assessment
  # (recode_age_of_assessment_0_0) -- the mirror image of the incident-case test used
  # for T2E outcomes (incident => disease_age - assessment_age > 0; load_disease()).
  #
  # NOTE on definition: a flag over ALL retained first-occurrence diseases is ~91%
  # positive in this cohort (it includes minor/acute conditions: back pain, acute
  # URTIs, chickenpox, gastritis, ...), so it has near-zero variance and is a poor
  # adjustment covariate / exclusion set. The curated map below restricts to major
  # chronic conditions (~15% positive), which is usable for both adjustment and a
  # healthy-at-baseline exclusion sensitivity analysis. Keyed by 3-char ICD-10 prefix.
  # Caveat: DZ_df only retains diseases with >= disease_min_cases INCIDENT cases, so a
  # disease that is common-prevalent but rare-incident could be absent; the major
  # chronic conditions below all clear that threshold.
  prevalent_major_disease = list(
    diabetes        = c("E10", "E11", "E12", "E13", "E14"),
    ischemic_heart  = c("I20", "I21", "I22", "I23", "I24", "I25"),
    heart_failure   = c("I50"),
    stroke_cvd      = c("I60", "I61", "I62", "I63", "I64"),
    copd            = c("J43", "J44"),
    ckd             = c("N18"),
    liver_cirrhosis = c("K70", "K74"),
    dementia        = c("F00", "F01", "F03", "G30")
  ),

  # Smoking structural-missingness: fields where NA is structural for never-smokers
  # and 0 is a valid imputed value (verified against UKB encoding before adding).
  #
  # Field 20161 (pack years of smoking): continuous, no special encoding.
  # Range is 0+ for smokers. 0 is a valid extension of the scale for never-smokers
  # (0 pack-years = never smoked). Verified: encoding is unconstrained continuous,
  # not an integer lookup table, so 0 does not collide with a special code.
  smoking_structural_fields = c(20161L),

  # Smoking fields flagged for manual review -- do NOT auto-recode.
  # These require checking whether 0 is a valid imputed value before recoding:
  #   20162: pack years as proportion of life span -- continuous like 20161, probably OK
  #   3436:  age started smoking (current smokers) -- encoding 100291; 0 is NOT valid (0 age)
  #   2867:  age started smoking (former smokers)  -- encoding 100291; 0 is NOT valid
  #   6183:  number of cigarettes currently smoked  -- encoding 100353 has -10=<1/day; 0 plausible
  #   2887:  number of cigarettes previously smoked -- encoding 100353; 0 plausible
  smoking_manual_review_fields = c(
    20162L,  # Pack years as proportion of life span (continuous; likely OK to set NA->0)
    3436L,   # Age started smoking in current smokers (DO NOT set to 0 -- age 0 invalid)
    2867L,   # Age started smoking in former smokers  (DO NOT set to 0 -- age 0 invalid)
    6183L,   # Cigarettes per day currently smoked (0 plausible but verify)
    2887L    # Cigarettes per day previously smoked (0 plausible but verify)
  ),

  # Known ordinal-encoding issues for audit (field_id as character -> note)
  # All ordinals are 0-indexed after recoding (min=0 = lowest exposure).
  known_ordinal_issues = c(
    "1558" = "Recoded: was 1=Daily, 6=Never; reversed to 0=Never, 5=Daily (verified)",
    "1628" = "Recoded: was 1=More, 3=Less; reversed to 0=Less, 2=More (verified)",
    "1249" = "Recoded: was 1=Most days, 4=Never smoked; reversed to 0=Never, 3=Most days (verified)",
    "3506" = "Recoded: was 1=More nowadays, 3=Less; reversed to 0=Less, 2=More (verified). NOTE: field 2644 was initially assumed to be this field but is actually a binary 0/1 field (smoked>=100 cigarettes lifetime) and is NOT recoded.",
    "1239" = "Recoded: irregular values 0/1/2 where 1>2 in intensity; remapped to 0=No, 1=Occasional, 2=Most days (verified)",
    "1031" = "Recoded: was 1=Almost daily, 6=Never; reversed to 0=Never, 5=Almost daily; value 7 (no friends) set to NA (verified)",
    "6160" = "Leisure/social activities: verify ordering direction"
  ),

  # Simple ordinal direction reversals: formula is new = max_val - x (0-indexed output).
  # IMPORTANT: preprocess_UKB_df / scale_ordinal_columns 0-indexes these fields before
  # recode_ordinal_directions runs, so max_val here is (n_levels - 1), NOT the UKB original max.
  # Add a field here only after manually verifying the UKB Data-Coding direction.
  ordinal_reversals = list(
    # Field 1558: Alcohol intake frequency.
    # UKB: 1=Daily...6=Never. Post-scale: 0=Daily...5=Never (6 levels).
    # Reversed so 0=Never (lowest), 5=Daily (highest).
    "1558" = 5L,
    # Field 1628: Alcohol intake versus 10 years previously.
    # UKB: 1=More, 2=Same, 3=Less. Post-scale: 0=More, 1=Same, 2=Less (3 levels).
    # Reversed so 0=Less, 2=More.
    "1628" = 2L,
    # Field 1249: Past tobacco smoking.
    # UKB: 1=Most days, 2=Occasionally, 3=Just tried, 4=Never. Post-scale: 0..3.
    # Reversed so 0=Never (lowest), 3=Most days (highest).
    "1249" = 3L,
    # Field 3506: Smoking compared to 10 years previous.
    # UKB encoding 100360: 1=More, 2=Same, 3=Less. Post-scale: 0..2.
    # Reversed so 0=Less, 2=More.
    # NOTE: field 2644 (encoding 100349, 0=No/1=Yes) is a different field
    # ("smoked >=100 cigarettes lifetime") and is NOT recoded.
    "3506" = 2L
  ),

  # Custom value remaps for fields where the coding is not a clean reversal.
  # Each entry is a named character vector: old_value_string -> new_numeric_value.
  # Values not listed are left as NA (treated as missing).
  # All remaps produce 0-indexed output (0 = lowest exposure).
  # Add a field here only after manually verifying the UKB Data-Coding direction.
  ordinal_custom_remaps = list(
    # Field 1239: Current tobacco smoking.
    # UKB: 0=No, 1=Yes on most/all days, 2=Only occasionally.
    # Intensity order: No(0) < Occasionally(2) < Most days(1).
    # Remapped to 0-indexed: 0->0 (none), 2->1 (occasional), 1->2 (daily).
    "1239" = c("0" = 0L, "1" = 2L, "2" = 1L),
    # Field 1031: Frequency of friend/family visits.
    # UKB: 1=Almost daily, 2=2-4x/week, 3=~once/week, 4=~once/month,
    #       5=Once every few months, 6=Never or almost never,
    #       7=No friends/family outside household (structural, not a frequency).
    # Reversed to 0-indexed: 6=Never->0, 1=Almost daily->5; 7 set to NA.
    "1031" = c("1"=5L, "2"=4L, "3"=3L, "4"=2L, "5"=1L, "6"=0L, "7"=NA_integer_)
  ),

  # Sun Exposure field IDs to drop after loading (phenotypic traits, not exposures).
  sun_exposure_exclude_fields = c(1717L, 1757L, 1747L),

  # Fields where negative values are UKB special codes (Do not know / Prefer not to
  # answer) that should be recoded to NA. Applied after E_long is assembled.
  negative_to_na_fields = c(
    2139L,  # Age first had sexual intercourse: -1=Do not know, -3=Prefer not to answer
    864L    # Days/week walked 10+ min: -2=Unable to walk (treat as missing, not 0 days)
  ),

  # Ordinal fields to convert to one-hot after loading. Each entry:
  #   field_id (int), label_map (named char vector: value_str -> level_label).
  # Negative UKB special codes are set to NA before encoding.
  ordinal_to_onehot_fields = list(
    list(
      field_id  = 1180L,
      # Field 1180: Morning/evening person (chronotype).
      # UKB: 1=Definitely morning, 2=More morning, 3=More evening, 4=Definitely evening.
      label_map = c("1" = "Definitely_morning", "2" = "More_morning",
                    "3" = "More_evening",        "4" = "Definitely_evening")
    )
  )
)

`%||%` <- function(x, y) if (!is.null(x)) x else y

ts_msg <- function(...) message(format(Sys.time(), "[%H:%M:%S]"), " ", ...)

############################################################
# 1) Dataloader helper bootstrap
############################################################

source_helpers <- function() {
  projID <<- cfg$projID
  old <- getwd()
  on.exit(setwd(old), add = TRUE)
  setwd(cfg$dataloader_dir)
  source(file.path(cfg$dataloader_dir, "dataloader_functions_parallel.R"))
  load_project(cfg$projID)
  invisible(TRUE)
}

############################################################
# 2) UKBdict field-discovery helpers
############################################################

# Returns field IDs from UKBdict for a given category, optionally excluding some.
get_category_fields <- function(UKBdict, category_id, exclude_field_ids = integer()) {
  UKBdict %>%
    filter(Category == category_id) %>%
    pull(FieldID) %>%
    unique() %>%
    setdiff(exclude_field_ids)
}

# Finds the first column in df whose name contains the pattern _f{field_id}_.
find_field_col <- function(df, field_id) {
  pat  <- paste0("_f", field_id, "_")
  cols <- grep(pat, names(df), value = TRUE, fixed = TRUE)
  if (length(cols) == 0L) NULL else cols[[1L]]
}

# Extracts the UKB field ID embedded in a canonicalized column name (e.g. _f1558_).
extract_field_id <- function(varname) {
  m <- regmatches(varname, regexpr("_f([0-9]+)_", varname))
  if (length(m) == 0L) return(NA_integer_)
  as.integer(sub("_f", "", sub("_$", "", sub("^_f", "", m))))
}

############################################################
# 3) Shared preprocessing utilities
############################################################

clean_names <- function(df) {
  nm <- gsub(" ", "_", names(df))
  nm <- make.names(nm)
  names(df) <- nm
  df
}

instance_pat <- function(inst) paste0("_", inst, "_")

canonicalize_names <- function(nm, inst, canon = 0L) {
  gsub(paste0("(_f[0-9]+)_", inst, "_"), paste0("\\1_", canon, "_"), nm)
}

canonicalize_visit <- function(df, inst, canon = 0L) {
  df <- as.data.frame(df)
  names(df) <- canonicalize_names(names(df), inst, canon)
  df <- clean_names(df)
  df$instance <- as.integer(inst)
  df %>% relocate(eid, instance)
}

empty_visit <- function(inst) data.frame(eid = integer(), instance = integer())

is_binary_like <- function(x) {
  ux <- unique(stats::na.omit(x))
  length(ux) > 0 && all(ux %in% c(0, 1, FALSE, TRUE))
}

fill_binary_zeros <- function(df) {
  vcols <- setdiff(names(df), c("eid", "instance"))
  bcols <- vcols[vapply(df[vcols], is_binary_like, logical(1))]
  for (cc in bcols) {
    for (ii in sort(unique(df$instance))) {
      idx <- df$instance == ii
      if (any(idx) && all(is.na(df[[cc]][idx])) && any(!is.na(df[[cc]][!idx]))) {
        df[[cc]][idx] <- 0
      }
    }
  }
  df
}

bind_visits <- function(dfs) {
  dfs <- dfs[vapply(dfs, nrow, integer(1)) > 0]
  if (length(dfs) == 0) return(data.frame(eid = integer(), instance = integer()))
  fill_binary_zeros(bind_rows(dfs))
}

ordinal_finder <- function(df) {
  vcols <- setdiff(names(df), c("eid", "instance"))
  vcols[vapply(vcols, function(v) {
    x <- suppressWarnings(as.numeric(df[[v]]))
    ux <- unique(stats::na.omit(x))
    length(ux) > 1 && max(ux, na.rm = TRUE) <= 5
  }, logical(1))]
}

orgEtoCat <- function(Elist, Elist_names) {
  ecats <- lapply(Elist, function(x) setdiff(names(x), c("eid", "instance")))
  names(ecats) <- Elist_names
  Eid_cat <- stack(ecats)
  colnames(Eid_cat) <- c("Eid", "Category")
  Eid_cat$Category <- as.character(Eid_cat$Category)
  as.data.frame(Eid_cat)
}

drop_instance_col <- function(df, inst = cfg$canonical_instance) {
  df <- df[df$instance == inst, setdiff(names(df), "instance"), drop = FALSE]
  rownames(df) <- NULL
  df
}

############################################################
# 4) Assessment timing derivation (covariates, section 5)
############################################################

# Adds assessment_year, assessment_month, assessment_season to covars_df.
# Uses field 53 (date of attending assessment centre). No-ops with a warning
# if the column is absent.
derive_assessment_timing <- function(covars_df) {
  date_col <- find_field_col(covars_df, 53L)
  if (is.null(date_col)) {
    warning("Field 53 (assessment date) not found in covariates; ",
            "assessment timing variables not added. ",
            "Check that 53L is in cfg$covar_fields.")
    return(covars_df)
  }
  dates <- suppressWarnings(as.Date(covars_df[[date_col]]))
  covars_df$assessment_year   <- as.integer(format(dates, "%Y"))
  covars_df$assessment_month  <- as.integer(format(dates, "%m"))
  covars_df$assessment_season <- dplyr::case_when(
    covars_df$assessment_month %in% c(12L, 1L, 2L) ~ "Winter",
    covars_df$assessment_month %in% c(3L, 4L, 5L)  ~ "Spring",
    covars_df$assessment_month %in% c(6L, 7L, 8L)  ~ "Summer",
    covars_df$assessment_month %in% c(9L, 10L, 11L) ~ "Autumn",
    TRUE ~ NA_character_
  )
  covars_df
}

############################################################
# 5) Mean-across-years pollution feature derivation (section 4)
############################################################

# For pollutants whose model years agree (NO2: 2005/06/07/10), collapse the
# per-year estimates into a single feature equal to the per-individual MEAN
# across all available model years (na.rm = TRUE). Single-year pollutants
# (PM2.5: 2010 only) reduce to that one estimate. Pollutants NOT listed in
# year_map_list (e.g. PM10, whose 2007/2010 surfaces disagree at r ~ 0.39) are
# left untouched as their separate per-year raw columns.
#
# Rationale: the UKB air pollution fields are ESCAPE land-use-regression spatial
# surfaces, not annual personal measurements; their spatial contrast is highly
# stable across model years. Averaging across years uses all the spatial
# information and avoids any past/future temporal-matching assumption (so there
# is no dependence on assessment_year). The per-year raw columns are RETAINED in
# air_df, so cross-year stability can be verified separately
# (see exposure_coding_check.R). NO2, PM10 and PM2.5 stay as separate features.
#
# Args:
#   air_df       : E_long[["Residential_Air_Pollution"]] data.frame
#   year_map_list: cfg$pollution_year_map (named list: pollutant -> named int vector
#                  of field_id_str -> calendar_year)
#
# Returns list(air_df = <updated, with {pollutant}_mean columns>, audit = <df>)
derive_mean_pollution_features <- function(air_df, year_map_list) {
  if (is.null(air_df) || nrow(air_df) == 0) {
    ts_msg("Residential_Air_Pollution is empty; skipping mean-across-years derivation")
    return(list(air_df = air_df, audit = data.frame()))
  }

  audit_rows <- list()

  for (pollutant in names(year_map_list)) {
    year_map  <- year_map_list[[pollutant]]   # named: "field_id" -> calendar_year
    field_ids <- as.integer(names(year_map))

    # Locate the per-year columns present in air_df
    col_map <- lapply(field_ids, function(fid) find_field_col(air_df, fid))
    names(col_map) <- as.character(field_ids)
    available <- field_ids[!vapply(col_map, is.null, logical(1))]

    if (length(available) == 0L) {
      warning("No columns found for pollutant '", pollutant, "'; skipping.")
      next
    }

    avail_cols  <- unlist(col_map[as.character(available)], use.names = FALSE)
    avail_years <- unname(year_map[as.character(available)])

    # Per-individual mean across available model years (na.rm = TRUE).
    # cbind() keeps this a matrix even when only one year is available.
    mat       <- do.call(cbind, lapply(avail_cols, function(cc) as.numeric(air_df[[cc]])))
    n_present <- rowSums(!is.na(mat))
    mean_vec  <- rowMeans(mat, na.rm = TRUE)
    mean_vec[n_present == 0L] <- NA_real_       # all years missing -> NA (not NaN)

    derived_col <- paste0(pollutant, "_mean")
    air_df[[derived_col]] <- mean_vec

    audit_rows[[length(audit_rows) + 1L]] <- data.frame(
      derived_variable    = derived_col,
      source_fields       = paste(available,   collapse = ","),
      source_years        = paste(avail_years, collapse = ","),
      n_years_available   = length(available),
      n_participants      = nrow(air_df),
      n_with_any_year     = sum(n_present > 0L),
      n_missing_all_years = sum(n_present == 0L),
      derivation_rule     = "mean_across_available_model_years"
    )
  }

  list(
    air_df = as.data.frame(air_df),
    audit  = if (length(audit_rows) > 0L) bind_rows(audit_rows) else data.frame()
  )
}

############################################################
# 6) Ordinal encoding audit (section 6)
############################################################

# Loads Codings.csv (UKB encoding_id -> value -> label; latin-1 encoded).
# Returns a keyed data.table, or NULL if the file is unavailable.
load_codings_table <- function(path) {
  if (!file.exists(path)) {
    warning("Codings.csv not found at: ", path, ". Audit will omit coding labels.")
    return(NULL)
  }
  ct <- tryCatch(
    fread(path, encoding = "Latin-1", colClasses = "character",
          col.names = c("encoding_id", "value", "meaning")),
    error = function(e) { warning("Failed to load Codings.csv: ", e$message); NULL }
  )
  if (!is.null(ct)) setkey(ct, encoding_id, value)
  ct
}

# Returns "value=label; ..." for positive codes of a given encoding_id.
coding_labels_string <- function(codings_dt, encoding_id_str) {
  if (is.null(codings_dt) || is.na(encoding_id_str) || !nzchar(encoding_id_str))
    return(NA_character_)
  rows <- codings_dt[encoding_id == encoding_id_str]
  rows <- rows[suppressWarnings(as.integer(value) > 0)]
  if (nrow(rows) == 0L) return(NA_character_)
  rows <- rows[order(suppressWarnings(as.integer(value)))]
  paste(paste0(rows$value, "=", rows$meaning), collapse = "; ")
}

# Scans all Elist categories for ordinal-looking variables, looks up UKBdict
# metadata and actual coding labels from Codings.csv, and flags known or
# suspected direction problems for manual review.
# Does NOT silently recode any variable.
audit_ordinal_encoding <- function(Elist, UKBdict, codings_dt = NULL) {
  rows <- list()

  for (cat_name in names(Elist)) {
    cat_df <- Elist[[cat_name]]
    vcols  <- setdiff(names(cat_df), c("eid", "instance"))

    for (vv in vcols) {
      x  <- suppressWarnings(as.numeric(cat_df[[vv]]))
      ux <- sort(unique(stats::na.omit(x)))
      # Audit integer-range variables that look ordinal (2–10 levels)
      if (length(ux) < 2L || max(ux, na.rm = TRUE) > 10L) next

      fid_int <- extract_field_id(vv)
      fid_str <- as.character(fid_int)

      dict_row <- if (!is.na(fid_int) && !is.null(UKBdict)) {
        UKBdict[UKBdict$FieldID == fid_int, ][1L, ]
      } else NULL

      field_desc    <- tryCatch(as.character(dict_row$Field),     error = function(e) NA_character_)
      value_type    <- tryCatch(as.character(dict_row$ValueType), error = function(e) NA_character_)
      coding_id     <- tryCatch(as.character(dict_row$Coding),    error = function(e) NA_character_)
      coding_labels <- coding_labels_string(codings_dt, coding_id)

      known_issue <- if (!is.na(fid_str) && fid_str %in% names(cfg$known_ordinal_issues))
                       cfg$known_ordinal_issues[[fid_str]]
                     else NA_character_

      rec_action <- if (!is.na(known_issue)) "manual_review" else "keep"

      # Heuristic: frequency/how-often fields starting at 1 are common reversal candidates
      if (is.na(known_issue) && !is.na(field_desc) &&
          grepl("frequen|how often", tolower(field_desc)) && min(ux) == 1L) {
        known_issue <- "Frequency field starting at 1 — verify whether 1 = most or least frequent"
        rec_action  <- "manual_review"
      }

      rows[[length(rows) + 1L]] <- data.frame(
        category           = cat_name,
        field_id           = fid_int,
        variable_name      = vv,
        field_description  = field_desc    %||% NA_character_,
        value_type         = value_type    %||% NA_character_,
        data_coding_id     = coding_id     %||% NA_character_,
        coding_labels      = coding_labels %||% NA_character_,
        observed_min       = min(ux),
        observed_max       = max(ux),
        n_unique_values    = length(ux),
        observed_values    = paste(ux, collapse = ","),
        known_issue        = known_issue,
        observed_ordering  = ifelse(!is.na(known_issue), "check_direction", "increasing"),
        recommended_action = rec_action
      )
    }
  }

  if (length(rows) == 0L) return(data.frame())
  bind_rows(rows)
}

############################################################
# 7) Structural missingness handler (section 7)
############################################################

# Applies documented structural-missingness rules to the Smoking Elist entry.
# Only recodes when smoking status == Never (field 20116 = 0 or "Never"):
#   - pack_years_of_smoking (field 20161): NA -> 0
# Additional smoking fields are flagged in the audit as manual_review without
# any automatic recoding.
#
# Returns list(Elist = <updated>, audit = <data.frame>)
apply_structural_missingness_rules <- function(Elist, UKBdict) {
  audit_rows <- list()

  if (!"Smoking" %in% names(Elist)) {
    warning("Smoking category not found in Elist; structural missingness rules not applied.")
    return(list(Elist = Elist, audit = data.frame()))
  }

  smoking_df  <- Elist[["Smoking"]]
  status_col  <- find_field_col(smoking_df, 20116L)

  # Helper: determine "is never smoker" mask using numeric or character coding
  is_never_smoker <- function(df, col) {
    if (is.null(col)) return(rep(FALSE, nrow(df)))
    sv      <- df[[col]]
    sv_num  <- suppressWarnings(as.numeric(as.character(sv)))
    sv_chr  <- as.character(sv)
    (!is.na(sv_num) & sv_num == 0L) | grepl("^[Nn]ever$", sv_chr)
  }

  never_mask <- is_never_smoker(smoking_df, status_col)

  # Auto-recode: pack years -> 0 for never-smokers
  for (fid in cfg$smoking_structural_fields) {
    col <- find_field_col(smoking_df, fid)
    if (is.null(col)) {
      warning("Field ", fid, " not found in Smoking category; skipping auto-recode.")
      next
    }
    n_before  <- sum(is.na(smoking_df[[col]]))
    recode_ok <- never_mask & is.na(smoking_df[[col]])
    smoking_df[[col]][recode_ok] <- 0
    n_recoded <- sum(recode_ok)
    n_after   <- sum(is.na(smoking_df[[col]]))

    audit_rows[[length(audit_rows) + 1L]] <- data.frame(
      category                   = "Smoking",
      variable                   = col,
      field_id                   = fid,
      n_missing_before           = n_before,
      n_recoded_to_zero          = n_recoded,
      n_missing_after            = n_after,
      rule_used                  = paste0("smoking_status(f20116)==Never -> ",
                                          col, " NA->0"),
      automatic_or_manual_review = "automatic"
    )
  }

  # Manual review: flag candidates without recoding
  for (fid in cfg$smoking_manual_review_fields) {
    col <- find_field_col(smoking_df, fid)
    if (is.null(col)) next
    n_missing <- sum(is.na(smoking_df[[col]]))
    audit_rows[[length(audit_rows) + 1L]] <- data.frame(
      category                   = "Smoking",
      variable                   = col,
      field_id                   = fid,
      n_missing_before           = n_missing,
      n_recoded_to_zero          = 0L,
      n_missing_after            = n_missing,
      rule_used                  = NA_character_,
      automatic_or_manual_review = "manual_review"
    )
  }

  Elist[["Smoking"]] <- smoking_df
  list(
    Elist = Elist,
    audit = if (length(audit_rows) > 0L) bind_rows(audit_rows) else data.frame()
  )
}

############################################################
# 8) Ordinal direction recoding (verified fields only)
############################################################

# Applies verified ordinal direction fixes from cfg$ordinal_reversals and
# cfg$ordinal_custom_remaps. Only fields manually verified against UKB
# Data-Coding are listed in those config entries.
#
# Simple reversal formula: new = max_val - x  (0-indexed; 0 = lowest exposure)
# Custom remap: value-by-value lookup; values not in the map become NA.
recode_ordinal_directions <- function(Elist) {
  for (cat_name in names(Elist)) {
    df <- Elist[[cat_name]]

    # Simple reversals
    for (fid_str in names(cfg$ordinal_reversals)) {
      col <- find_field_col(df, as.integer(fid_str))
      if (is.null(col)) next
      max_val <- cfg$ordinal_reversals[[fid_str]]
      x <- suppressWarnings(as.numeric(df[[col]]))
      df[[col]] <- ifelse(is.na(x), NA_real_, max_val - x)
      ts_msg("Reversed ordinal: ", col, " (field ", fid_str, ", max=", max_val, ")")
    }

    # Custom value remaps
    for (fid_str in names(cfg$ordinal_custom_remaps)) {
      col <- find_field_col(df, as.integer(fid_str))
      if (is.null(col)) next
      remap <- cfg$ordinal_custom_remaps[[fid_str]]
      x_chr <- as.character(suppressWarnings(as.numeric(df[[col]])))
      new_x <- remap[x_chr]            # NA for any value not in the map
      df[[col]] <- as.numeric(new_x)
      ts_msg("Custom remap ordinal: ", col, " (field ", fid_str, ")")
    }

    Elist[[cat_name]] <- df
  }
  Elist
}

############################################################
# 9) Post-load exposure cleanups
############################################################

# Drop phenotypic (non-exposure) columns from Sun_Exposure by field ID.
drop_sun_exposure_phenotypic_cols <- function(Elist) {
  df <- Elist[["Sun_Exposure"]]
  if (is.null(df)) return(Elist)
  for (fid in cfg$sun_exposure_exclude_fields) {
    col <- find_field_col(df, fid)
    if (!is.null(col)) {
      drop_cols <- grep(paste0("_f", fid, "_"), names(df), value = TRUE)
      df <- df[, setdiff(names(df), drop_cols), drop = FALSE]
      ts_msg("Dropped Sun_Exposure phenotypic field ", fid, ": ", paste(drop_cols, collapse = ", "))
    }
  }
  Elist[["Sun_Exposure"]] <- df
  Elist
}

# Recode negative values to NA for specified fields (UKB special codes).
recode_negatives_to_na <- function(Elist) {
  for (cat_name in names(Elist)) {
    df <- Elist[[cat_name]]
    for (fid in cfg$negative_to_na_fields) {
      col <- find_field_col(df, fid)
      if (is.null(col)) next
      x <- suppressWarnings(as.numeric(df[[col]]))
      n_recoded <- sum(!is.na(x) & x < 0)
      if (n_recoded > 0) {
        df[[col]] <- ifelse(!is.na(x) & x < 0, NA_real_, x)
        ts_msg("Recoded ", n_recoded, " negative values to NA: ", col)
      }
    }
    Elist[[cat_name]] <- df
  }
  Elist
}

# Convert specified ordinal fields to one-hot dummy columns.
# The original column is replaced by binary indicator columns.
convert_ordinal_to_onehot <- function(Elist) {
  for (spec in cfg$ordinal_to_onehot_fields) {
    fid       <- spec$field_id
    label_map <- spec$label_map
    for (cat_name in names(Elist)) {
      df  <- Elist[[cat_name]]
      col <- find_field_col(df, fid)
      if (is.null(col)) next
      x_chr <- as.character(suppressWarnings(as.numeric(df[[col]])))
      # Derive column name prefix from the existing column (drop trailing suffix after _fXXXX_)
      prefix <- sub("(_f[0-9]+_[0-9]+_[0-9]+).*$", "\\1", col)
      for (val_str in names(label_map)) {
        new_col <- paste0(prefix, "_", label_map[[val_str]])
        df[[new_col]] <- as.integer(!is.na(x_chr) & x_chr == val_str)
        df[[new_col]][is.na(x_chr) | x_chr == "NA"] <- NA_integer_
      }
      df[[col]] <- NULL
      ts_msg("Converted ordinal to one-hot: ", col, " -> ", length(label_map), " columns")
      Elist[[cat_name]] <- df
    }
  }
  Elist
}

############################################################
# 10) Audit helpers
############################################################

audit_df <- function(df, category) {
  if (!all(c("eid", "instance") %in% names(df))) {
    return(data.table(category = character(), variable = character()))
  }
  vcols <- setdiff(names(df), c("eid", "instance"))
  rbindlist(lapply(vcols, function(vv) {
    rbindlist(lapply(sort(unique(df$instance)), function(ii) {
      idx <- df$instance == ii
      data.table(
        category = category, variable = vv, instance = ii,
        n = sum(idx),
        n_nonmiss = sum(!is.na(df[[vv]][idx])),
        miss_rate = mean(is.na(df[[vv]][idx]))
      )
    }))
  }), fill = TRUE)
}

write_audits <- function(audit, out_rds) {
  prefix <- sub("\\.rds$", "", out_rds)
  for (nm in names(audit)) {
    tbl <- audit[[nm]]
    if (is.null(tbl) || (is.data.frame(tbl) && nrow(tbl) == 0L)) next
    fwrite(as.data.table(tbl),
           file = paste0(prefix, "_audit_", nm, ".tsv"), sep = "\t")
  }
}

############################################################
# 9) Proteomics loader (all instances)
############################################################

load_proteomics <- function() {
  ts_msg("Loading Olink proteomics (instances: ", paste(cfg$instances, collapse = ","), ")")

  prot_raw <- fread(
    cfg$paths$olink,
    select = c("eid", "ins_index", "protein_id", "result"),
    showProgress = FALSE
  )
  prot_raw <- prot_raw[ins_index %in% cfg$instances]

  prot_id_map <- fread(cfg$paths$protein_id, fill = TRUE, quote = "", encoding = "UTF-8")
  prot_id_map[, c("prot_id", "prot_name") := tstrsplit(meaning, ";", fixed = TRUE, keep = 1:2)]
  prot_id_map[, prot_id := gsub("-", "_", prot_id)]

  ts_msg("Casting proteomics to wide format")
  prot_wide  <- dcast.data.table(prot_raw, eid + ins_index ~ protein_id, value.var = "result")
  coding_map <- setNames(prot_id_map$prot_id, as.character(prot_id_map$coding))
  old_nm <- setdiff(names(prot_wide), c("eid", "ins_index"))
  new_nm <- coding_map[old_nm]
  new_nm[is.na(new_nm)] <- old_nm[is.na(new_nm)]
  setnames(prot_wide, old_nm, new_nm)
  setnames(prot_wide, "ins_index", "instance")
  prot_wide[, instance := as.integer(instance)]

  prot_ids <- setdiff(names(prot_wide), c("eid", "instance"))

  vm       <- unique(prot_wide[, .(eid, instance)])
  split_df <- vm[, .(
    has_i0    = any(instance == 0L),
    has_i2    = any(instance == 2L),
    has_i3    = any(instance == 3L),
    n_visits  = uniqueN(instance),
    instances = paste(sort(unique(instance)), collapse = ",")
  ), by = eid]
  split_df[, analysis_set := fifelse(
    has_i0 & !has_i2 & !has_i3, "train_baseline_only",
    fifelse(has_i2 | has_i3, "holdout_repeat_proteomics", "other")
  )]

  list(
    prot_long = as.data.frame(prot_wide),
    protIDs   = prot_ids,
    split_df  = as.data.frame(split_df)
  )
}

############################################################
# 10) Exposure loader helpers
############################################################

load_lifestyle_category <- function(path_id, category_name, raw_df) {
  visit_dfs <- lapply(cfg$instances, function(inst) {
    out <- tryCatch(
      preprocess_UKB_df(
        path_id = path_id, missingness = 1,
        timepoint = instance_pat(inst),
        df = raw_df, feature_engineer = TRUE
      ),
      error = function(e) {
        warning("path ", path_id, " instance ", inst, ": ", conditionMessage(e))
        NULL
      }
    )
    if (is.null(out) || ncol(out) <= 1L) return(empty_visit(inst))
    canonicalize_visit(out, inst, cfg$canonical_instance)
  })
  bind_visits(visit_dfs)
}

load_field_instances <- function(raw_df, label, process_fn = NULL) {
  visit_dfs <- lapply(cfg$instances, function(inst) {
    out <- tryCatch(
      UKB_instances(raw_df, instance_pat(inst)),
      error = function(e) {
        warning(label, " instance ", inst, ": ", conditionMessage(e))
        NULL
      }
    )
    if (is.null(out) || ncol(out) <= 1L) return(empty_visit(inst))
    if (!is.null(process_fn)) out <- process_fn(out, inst)
    if (is.null(out) || ncol(out) <= 1L) return(empty_visit(inst))
    canonicalize_visit(out, inst, cfg$canonical_instance)
  })
  bind_visits(visit_dfs)
}

process_multi_onehot <- function(df, inst) clean_names(UKB_onehot_handle(UKB_multiarray_handle(df)))

# Load individual household income (field 738) and encode its label strings into
# a 0-indexed ordinal column (0 = lowest bracket ... 4 = highest income). This
# field is NOT in the deprivation category and is stored as label strings, so it
# bypasses the ordinal feature-engineering pipeline used for lifestyle paths and
# must be encoded here. "Do not know"/"Prefer not to answer" and any value not in
# cfg$income_label_map become NA. Returns an (eid, instance) x income data.frame
# (instances per cfg$instances), or NULL if the field is unavailable.
load_income_ordinal <- function(UKBdict) {
  if (!cfg$income_field %in% UKBdict$FieldID) {
    warning("Income field ", cfg$income_field, " not in UKBdict; income not loaded.")
    return(NULL)
  }
  raw_inc <- fast_dataloader_viafield(UKBdict, cfg$income_field, directoryInfo)
  inc_df  <- load_field_instances(raw_inc, "Income")
  col <- find_field_col(inc_df, cfg$income_field)
  if (is.null(col)) {
    warning("Income column (field ", cfg$income_field, ") not found after loading.")
    return(NULL)
  }
  x <- as.character(inc_df[[col]])
  inc_df[[col]] <- unname(cfg$income_label_map[x])   # NA for unmapped / special codes
  ts_msg("Encoded income (field ", cfg$income_field, ") to ordinal 0-",
         max(cfg$income_label_map), ": ", col,
         " (n non-missing = ", sum(!is.na(inc_df[[col]])), ")")
  inc_df
}

############################################################
# 11) Exposure loader (all categories, all instances)
############################################################

load_exposures <- function(UKBdict) {
  Elist <- list()

  # -- Lifestyle categories (loaded via path ID) --
  for (cat_name in names(cfg$lifestyle_path_ids)) {
    path_id <- unname(cfg$lifestyle_path_ids[[cat_name]])
    ts_msg("Loading E: ", cat_name, " (path ", path_id, ")")
    load_path(cfg$projID, path_id)
    raw_df <- fast_dataloader(path, directoryInfo)
    Elist[[cat_name]] <- load_lifestyle_category(path_id, cat_name, raw_df)
    gc()
  }

  # -- Deprivation Indices (UKB category 76: Indices of Multiple Deprivation) --
  ts_msg("Loading E: Deprivation_Indices (category 76)")
  deprivation_fields <- get_category_fields(UKBdict, cfg$deprivation_category)
  if (length(deprivation_fields) > 0L) {
    raw_dep <- fast_dataloader_viafield(UKBdict, deprivation_fields, directoryInfo)
    Elist[["Deprivation_Indices"]] <- load_field_instances(raw_dep, "Deprivation_Indices")
  } else {
    warning("No fields found for UKB category ", cfg$deprivation_category,
            "; Deprivation_Indices not loaded.")
  }

  # Append individual household income (field 738) as an ordinal feature.
  # (Category 100066, not in deprivation category 76; encoded to 0-4 ordinal.)
  ts_msg("Loading E: individual income (field ", cfg$income_field,
         ") -> Deprivation_Indices")
  inc_df <- load_income_ordinal(UKBdict)
  if (!is.null(inc_df) && ncol(inc_df) > 2L) {
    if (!is.null(Elist[["Deprivation_Indices"]])) {
      Elist[["Deprivation_Indices"]] <- full_join(
        Elist[["Deprivation_Indices"]], inc_df, by = c("eid", "instance")
      )
    } else {
      Elist[["Deprivation_Indices"]] <- inc_df
    }
  }
  gc()

  # -- Geographical Measures (UKB category 711, excluding assessment centre field 54) --
  ts_msg("Loading E: Geographical_Measures (category 711, excluding field 54)")
  geographical_fields <- get_category_fields(
    UKBdict, cfg$geographical_category, exclude_field_ids = cfg$geographical_exclude
  )
  if (length(geographical_fields) > 0L) {
    raw_geo <- fast_dataloader_viafield(UKBdict, geographical_fields, directoryInfo)
    Elist[["Geographical_Measures"]] <- load_field_instances(raw_geo, "Geographical_Measures")
  } else {
    warning("No fields found for UKB category ", cfg$geographical_category,
            " after exclusions; Geographical_Measures not loaded.")
  }
  gc()

  # -- Residential Air Pollution (NO2, NOx, PM variants, traffic; section 3) --
  ts_msg("Loading E: Residential_Air_Pollution")
  air_fids_present <- cfg$air_pollution_fields[
    cfg$air_pollution_fields %in% UKBdict$FieldID
  ]
  if (length(air_fids_present) > 0L) {
    raw_air <- fast_dataloader_viafield(UKBdict, air_fids_present, directoryInfo)
    Elist[["Residential_Air_Pollution"]] <- load_field_instances(
      raw_air, "Residential_Air_Pollution"
    )
    if (length(air_fids_present) < length(cfg$air_pollution_fields)) {
      missing_fids <- setdiff(cfg$air_pollution_fields, air_fids_present)
      warning("Air pollution fields not found in UKBdict: ",
              paste(missing_fids, collapse = ", "))
    }
  } else {
    warning("No air pollution fields found in UKBdict; ",
            "Residential_Air_Pollution not loaded.")
  }
  gc()

  # -- Residential Noise Pollution (section 3) --
  ts_msg("Loading E: Residential_Noise_Pollution")
  noise_fids_present <- cfg$noise_pollution_fields[
    cfg$noise_pollution_fields %in% UKBdict$FieldID
  ]
  if (length(noise_fids_present) > 0L) {
    raw_noise <- fast_dataloader_viafield(UKBdict, noise_fids_present, directoryInfo)
    Elist[["Residential_Noise_Pollution"]] <- load_field_instances(
      raw_noise, "Residential_Noise_Pollution"
    )
  } else {
    warning("No noise pollution fields found in UKBdict; ",
            "Residential_Noise_Pollution not loaded.")
  }
  gc()

  # -- Vitamins (multi-array one-hot encoding) --
  ts_msg("Loading E: Vitamins")
  raw_vit <- fast_dataloader_viafield(UKBdict, cfg$vitamin_fields, directoryInfo)
  Elist[["Vitamins"]] <- load_field_instances(raw_vit, "Vitamins", process_multi_onehot)
  gc()

  Elist
}

############################################################
# 12) Covariate loader (all instances, includes assessment timing)
############################################################

load_covariates <- function(UKBdict) {
  ts_msg("Loading covariates (fields: ", paste(cfg$covar_fields, collapse = ","), ")")
  raw_cov   <- fast_dataloader_viafield(UKBdict, cfg$covar_fields, directoryInfo)
  covars_df <- load_field_instances(raw_cov, "Covariates")

  # Propagate static fields (sex, genetic PCs) across instances via fill
  static_cols <- grep("^sex_f31_0_0$|^genetic_principal_components_f22009_0_",
                      names(covars_df), value = TRUE)
  if (length(static_cols) > 0L) {
    covars_df <- covars_df %>%
      arrange(eid, instance) %>%
      group_by(eid) %>%
      fill(all_of(static_cols), .direction = "downup") %>%
      ungroup()
  }

  # Derived age interaction terms
  age_col <- "age_when_attended_assessment_centre_f21003_0_0"
  sex_col <- "sex_f31_0_0"
  if (all(c(age_col, sex_col) %in% names(covars_df))) {
    male_ind <- ifelse(as.character(covars_df[[sex_col]]) == "Male", 1L, 0L)
    covars_df$age2     <- as.numeric(covars_df[[age_col]])^2
    covars_df$age_sex  <- as.numeric(covars_df[[age_col]]) * male_ind
    covars_df$age2_sex <- as.numeric(covars_df[[age_col]])^2 * male_ind
  }

  # Assessment centre as factor
  ctr_col <- "uk_biobank_assessment_centre_f54_0_0"
  if (ctr_col %in% names(covars_df)) covars_df[[ctr_col]] <- as.factor(covars_df[[ctr_col]])

  # Assessment timing variables (year, month, season) from field 53 (section 5)
  ts_msg("Deriving assessment timing variables from field 53")
  covars_df <- derive_assessment_timing(covars_df)

  # Medications (one-hot, joined to covariates)
  ts_msg("Loading medications")
  raw_med <- fast_dataloader_viafield(UKBdict, cfg$medication_fields, directoryInfo)

  process_medications <- function(df, inst) {
    out <- UKB_multiarray_handle(df) %>%
      pivot_longer(-eid, names_to = "variable", values_to = "value") %>%
      group_by(eid) %>%
      summarise(
        combined = paste(unique(unlist(strsplit(value[!is.na(value)], ";"))),
                         collapse = ";"),
        .groups = "drop"
      )
    out$combined[out$combined == ""] <- NA
    clean_names(UKB_onehot_handle(out))
  }

  medi_df <- load_field_instances(raw_med, "Medications", process_medications)
  if (ncol(medi_df) > 2L) {
    covars_df <- full_join(covars_df, medi_df, by = c("eid", "instance"))
  }

  covars_df <- covars_df %>% arrange(eid, instance)
  list(
    covars_df   = as.data.frame(covars_df),
    covars_list = setdiff(names(covars_df), c("eid", "instance"))
  )
}

############################################################
# 13) Disease T2E loader
############################################################

load_disease <- function(UKBdict) {
  if (cfg$skip_disease) {
    ts_msg("Skipping disease loading (HEAP_SKIP_DISEASE=TRUE)")
    return(list(DZ_df = data.frame(eid = integer()), DZ_ids = character()))
  }

  ts_msg("Loading first-occurrence disease fields (takes ~10 min)")

  raw_base    <- fast_dataloader_viafield(UKBdict, cfg$baseline_fields_t2e, directoryInfo)
  baseline_df <- UKB_instances(raw_base, "_0_")
  baseline_df <- baseline_df %>%
    mutate(month_of_birth_f52_0_0 = recode(
      month_of_birth_f52_0_0,
      January = 1, February = 2, March = 3, April = 4,
      May = 5, June = 6, July = 7, August = 8,
      September = 9, October = 10, November = 11, December = 12
    ))
  baseline_df$birth_date <- as.Date(with(baseline_df,
    paste(year_of_birth_f34_0_0, month_of_birth_f52_0_0, "15", sep = "-")
  ), "%Y-%m-%d")

  First_Occur_IDs <- UKBdict %>%
    filter(Category %in% cfg$disease_categories) %>%
    pull(FieldID)
  DZdf   <- fast_dataloader_viafield(UKBdict, First_Occur_IDs, directoryInfo)
  T2E_df <- list(DZdf, baseline_df) %>% reduce(full_join, by = "eid")
  T2E_df <- time2event_ages(T2E_df)

  DZid_date <- UKBdict %>%
    filter(Category %in% cfg$disease_categories) %>%
    filter(grepl("date", descriptive_colnames)) %>%
    pull(descriptive_colnames)
  DZid_age <- gsub("date", "age", DZid_date)

  age_of_DZ <- function(date_col, age_col, df) {
    df[[date_col]] <- as.Date(df[[date_col]])
    df[[age_col]]  <- round(as.numeric((df[[date_col]] - df$birth_date) / 365), 1)
    df
  }
  for (i in seq_along(DZid_date)) T2E_df <- age_of_DZ(DZid_date[i], DZid_age[i], T2E_df)

  keep_cols <- c("eid", DZid_age,
                 "recode_age_of_assessment_0_0", "recode_age_of_death_0_0",
                 "age_of_removal_0_0", "age_of_lastfollowup")
  DZ_age <- T2E_df %>% select(all_of(intersect(keep_cols, names(T2E_df))))

  DZcount <- vapply(DZid_age, function(i) {
    if (!i %in% names(DZ_age)) return(0L)
    as.integer(sum(DZ_age[[i]] - DZ_age$recode_age_of_assessment_0_0 > 0, na.rm = TRUE))
  }, integer(1))

  DZ_ids_keep <- names(DZcount[DZcount >= cfg$disease_min_cases])
  ts_msg("Retaining ", length(DZ_ids_keep), " disease outcomes (>= ", cfg$disease_min_cases, " cases)")

  DZ_df <- DZ_age %>% select(all_of(c(
    "eid", DZ_ids_keep,
    "recode_age_of_assessment_0_0", "recode_age_of_death_0_0",
    "age_of_removal_0_0", "age_of_lastfollowup"
  )))

  list(DZ_df = as.data.frame(DZ_df), DZ_ids = DZ_ids_keep)
}

############################################################
# 13b) Prevalent disease at baseline (covariate derivation)
############################################################

# Derives baseline prevalent-disease features from the disease T2E object.
# A retained first-occurrence disease is "prevalent at baseline" for a participant
# when its first-occurrence age <= the participant's age at baseline assessment
# (recode_age_of_assessment_0_0). This is the exact mirror of the incident-case test
# at the end of load_disease() (incident => disease_age - assessment_age > 0).
#
# Produces, over the curated major-chronic-disease set (cfg$prevalent_major_disease):
#   prevalent_<domain>          -- 0/1 per disease domain (diabetes, ischemic_heart, ...)
#   prevalent_major_disease     -- 0/1, any major chronic disease present at baseline
#   prevalent_major_disease_n   -- count of major chronic disease DOMAINS at baseline
#                                  (multimorbidity burden; a domain counts once even if
#                                   several of its ICD codes are present)
# and, for reference only (NOT recommended as a covariate -- ~91% positive):
#   prevalent_any_disease       -- 0/1 over all retained first-occurrence diseases
#
# Prevalence is defined only for participants with a baseline assessment age; rows
# without one are dropped. Returns an (eid x features) data.frame keyed by eid, or
# NULL if disease was skipped / DZ_df is empty.
derive_prevalent_disease <- function(disease_obj, major_map = cfg$prevalent_major_disease) {
  dz     <- disease_obj$DZ_df
  dz_ids <- intersect(disease_obj$DZ_ids, names(dz))
  asmt   <- "recode_age_of_assessment_0_0"
  if (is.null(dz) || nrow(dz) == 0L || length(dz_ids) == 0L || !asmt %in% names(dz)) {
    warning("derive_prevalent_disease: disease data unavailable; prevalent flags skipped.")
    return(NULL)
  }

  dz         <- as.data.frame(dz)
  age_assess <- suppressWarnings(as.numeric(dz[[asmt]]))

  # first-occurrence age - assessment age; <= 0 (and non-NA) => prevalent at baseline
  M    <- as.matrix(dz[, dz_ids, drop = FALSE]); storage.mode(M) <- "double"
  diff <- M - age_assess
  prev <- (diff <= 0) & !is.na(diff)

  # ICD-10 3-char prefix embedded in each column name:
  #   age_<letter><digits>_first_reported_<...>_fNNNNN_0_0
  icd <- toupper(sub("^age_([a-z][0-9]+)_.*$", "\\1", dz_ids))

  out          <- data.frame(eid = dz$eid)
  domain_flags <- matrix(0L, nrow = nrow(dz), ncol = 0L)
  for (domain in names(major_map)) {
    cols <- dz_ids[icd %in% major_map[[domain]]]
    flag <- if (length(cols)) as.integer(rowSums(prev[, cols, drop = FALSE]) > 0) else 0L
    if (length(cols) == 0L)
      warning("derive_prevalent_disease: no retained diseases matched domain '", domain, "'.")
    out[[paste0("prevalent_", domain)]] <- flag
    domain_flags <- cbind(domain_flags, flag)
  }
  out$prevalent_major_disease   <- as.integer(rowSums(domain_flags) > 0)
  out$prevalent_major_disease_n <- as.integer(rowSums(domain_flags))
  out$prevalent_any_disease     <- as.integer(rowSums(prev) > 0)

  out <- out[!is.na(age_assess), , drop = FALSE]
  rownames(out) <- NULL

  n_major <- sum(out$prevalent_major_disease)
  ts_msg("Prevalent disease at baseline: ", n_major, " / ", nrow(out),
         sprintf(" (%.1f%%)", 100 * n_major / nrow(out)),
         " with >=1 major chronic disease")
  out
}

############################################################
# 14) Main
############################################################

main <- function() {
  ts_msg("HEAP_loader starting")

  source_helpers()
  source(file.path(cfg$dataloader_dir, "time2event_functions.R"))

  UKBdict    <- fread(cfg$paths$path_all)
  codings_dt <- load_codings_table(cfg$paths$codings)
  plan(multisession, workers = 1)

  # Proteomics
  prot_obj <- load_proteomics()

  # Exposures (lifestyle + deprivation + geo + air pollution + noise + vitamins)
  E_long <- load_exposures(UKBdict)

  # Drop phenotypic (non-exposure) columns from Sun_Exposure
  ts_msg("Dropping Sun_Exposure phenotypic fields (skin colour, facial ageing, hair colour)")
  E_long <- drop_sun_exposure_phenotypic_cols(E_long)

  # Apply verified ordinal direction reversals (e.g. alcohol frequency field 1558)
  ts_msg("Applying ordinal direction reversals")
  E_long <- recode_ordinal_directions(E_long)

  # Recode negative UKB special codes to NA (e.g. age at first sexual intercourse)
  ts_msg("Recoding negative special codes to NA")
  E_long <- recode_negatives_to_na(E_long)

  # Convert specified ordinal fields to one-hot (e.g. morning/evening chronotype)
  ts_msg("Converting ordinal fields to one-hot")
  E_long <- convert_ordinal_to_onehot(E_long)

  # Covariates (includes assessment timing: assessment_year, _month, _season)
  covar_obj <- load_covariates(UKBdict)

  # Structural missingness: recode smoking pack-years for never-smokers
  ts_msg("Applying structural missingness rules (smoking)")
  sm_result <- apply_structural_missingness_rules(E_long, UKBdict)
  E_long    <- sm_result$Elist

  # Multi-year pollution derivation: per-individual MEAN across model years.
  # (No temporal matching; per-year raw columns are retained for cross-year
  # verification in exposure_coding_check.R. NO2/PM10/PM2.5 stay separate.)
  ts_msg("Deriving mean-across-years pollution features (NO2, PM10, PM2.5)")
  poll_result <- derive_mean_pollution_features(
    air_df        = E_long[["Residential_Air_Pollution"]],
    year_map_list = cfg$pollution_year_map
  )
  E_long[["Residential_Air_Pollution"]] <- poll_result$air_df

  # Disease
  disease_obj <- load_disease(UKBdict)

  # Prevalent disease at baseline (covariate / exclusion feature; section 13b).
  # Merged into the covariate frames by eid -- a time-fixed baseline attribute, so the
  # same value is broadcast across a participant's visits. Participants without a
  # baseline assessment age (or absent from the disease data) get NA. Deliberately NOT
  # added to covars_list: the DEFAULT model covariate set is left unchanged. Reference
  # the prevalent_* columns explicitly to adjust, or filter on prevalent_major_disease
  # (==0) for a healthy-at-baseline exclusion sensitivity analysis.
  prevalent_df <- derive_prevalent_disease(disease_obj)
  if (!is.null(prevalent_df)) {
    prevalent_cols      <- setdiff(names(prevalent_df), "eid")
    covar_obj$covars_df <- left_join(covar_obj$covars_df, prevalent_df, by = "eid")
  } else {
    prevalent_cols <- character()
  }

  # Shared lookups (computed after all E_long modifications)
  Elist_names <- names(E_long)
  Eid_cat     <- orgEtoCat(E_long, Elist_names)
  # NOTE: ordinal_finder() classifies "ordinal" as max(value) <= 5, which
  # mis-flags small-scale CONTINUOUS scores (England IMD income/employment/crime/
  # health, pm2.5 absorbance, all in [0,5]). This stored ordinalIDs is now
  # SUPERSEDED at module-load time: as_pxs_baseline()/as_pxs_longitudinal() in
  # workflow/00_paths.R re-route ordinalIDs by the DECLARED variable_type in
  # config/exposure_sets/analysis_exposures.tsv. Retained only as a fallback for
  # objects/contexts where that config is unavailable.
  ordinalIDs  <- unique(unlist(lapply(E_long, ordinal_finder), use.names = FALSE))

  # Baseline views (instance 0, no instance column)
  ts_msg("Deriving baseline (instance-0) views")
  prot_baseline   <- drop_instance_col(prot_obj$prot_long)
  E_baseline      <- lapply(E_long, drop_instance_col)
  covars_baseline <- drop_instance_col(covar_obj$covars_df)

  # Audits
  ts_msg("Building audit tables")
  prot_audit <- as.data.table(prot_obj$prot_long)[, .(
    n_eid = uniqueN(eid), n_rows = .N
  ), by = instance][order(instance)]

  exposure_audit <- rbindlist(
    Map(audit_df, E_long, names(E_long)), fill = TRUE
  )
  covar_audit    <- audit_df(covar_obj$covars_df, "Covariates")
  ordinal_audit  <- audit_ordinal_encoding(E_long, UKBdict, codings_dt)

  audit <- list(
    proteomics                = as.data.frame(prot_audit),
    exposures                 = as.data.frame(exposure_audit),
    covariates                = as.data.frame(covar_audit),
    ordinal_encoding          = ordinal_audit,
    structural_missingness    = sm_result$audit,
    residential_air_pollution = poll_result$audit
  )

  # Assemble HEAP
  HEAP <- list(
    meta = list(
      object_version     = "1.1",
      created_at         = as.character(Sys.time()),
      instances          = cfg$instances,
      canonical_instance = cfg$canonical_instance,
      projID             = cfg$projID,
      n_proteins         = length(prot_obj$protIDs),
      n_eid_baseline     = nrow(prot_baseline),
      n_exposure_cats    = length(Elist_names),
      exposure_cats      = Elist_names,
      skip_disease       = cfg$skip_disease,
      # Baseline prevalent-disease columns merged into covars_long / covars_baseline
      # (NOT in covars_list; opt-in for adjustment / exclusion sensitivity analyses).
      prevalent_disease_cols = prevalent_cols
    ),

    # Shared lookups
    protIDs     = prot_obj$protIDs,
    Elist_names = Elist_names,
    Eid_cat     = Eid_cat,
    ordinalIDs  = ordinalIDs,
    covars_list = covar_obj$covars_list,

    # Longitudinal (eid x instance)
    prot_long   = prot_obj$prot_long,
    E_long      = E_long,
    covars_long = covar_obj$covars_df,
    split_df    = prot_obj$split_df,

    # Baseline (instance 0, no instance column)
    prot_baseline   = prot_baseline,
    E_baseline      = E_baseline,
    covars_baseline = covars_baseline,

    # Disease outcomes
    disease = disease_obj,

    # Baseline prevalent-disease features (eid x flags); also merged into the
    # covariate frames. NULL when disease loading is skipped.
    prevalent_disease = prevalent_df,

    # Audit
    audit = audit
  )

  dir.create(dirname(cfg$paths$out_rds), recursive = TRUE, showWarnings = FALSE)
  saveRDS(HEAP, file = cfg$paths$out_rds)
  write_audits(HEAP$audit, cfg$paths$out_rds)

  ts_msg("Saved HEAP.rds: ", cfg$paths$out_rds)
  ts_msg("Proteins: ", HEAP$meta$n_proteins)
  ts_msg("Baseline participants: ", HEAP$meta$n_eid_baseline)
  ts_msg("Exposure categories: ", paste(Elist_names, collapse = ", "))
  ts_msg("Disease outcomes retained: ", length(HEAP$disease$DZ_ids))
  ts_msg("Proteomics by instance:")
  print(prot_audit)

  invisible(HEAP)
}

if (identical(environment(), globalenv()) && Sys.getenv("HEAP_LOADER_SOURCE_ONLY") != "1") {
  main()
}

############################################################
# Compatibility accessors -- also used by downstream module scripts.
# Source this file with HEAP_LOADER_SOURCE_ONLY=1 to get these without
# re-running main(), or simply copy them into your module scripts.
############################################################

as_pxs_baseline <- function(heap) {
  list(
    Elist       = heap$E_baseline,
    Elist_names = heap$Elist_names,
    Eid_cat     = heap$Eid_cat,
    # Route ordinalIDs by the DECLARED variable_type (mirrors the authoritative
    # definition in workflow/00_paths.R, sourced at the top of this file) so this
    # accessor can never silently revert the small-scale-continuous fix. Falls back
    # to the loader heuristic only if 00_paths.R was not sourced.
    ordinalIDs  = if (exists("heap_resolve_ordinal_ids", mode = "function"))
                    heap_resolve_ordinal_ids(heap, .heap_present_exposures(heap$E_baseline))
                  else heap$ordinalIDs,
    UKBprot_df  = heap$prot_baseline,
    protIDs     = heap$protIDs,
    covars_df   = heap$covars_baseline,
    covars_list = heap$covars_list
  )
}

as_pxs_longitudinal <- function(heap) {
  list(
    Elist              = heap$E_long,
    Elist_names        = heap$Elist_names,
    Eid_cat            = heap$Eid_cat,
    # Declared-type routing (see as_pxs_baseline note above); heuristic fallback only
    # if 00_paths.R was not sourced.
    ordinalIDs         = if (exists("heap_resolve_ordinal_ids", mode = "function"))
                           heap_resolve_ordinal_ids(heap, .heap_present_exposures(heap$E_long))
                         else heap$ordinalIDs,
    UKBprot_df         = heap$prot_long,
    protIDs            = heap$protIDs,
    covars_df          = heap$covars_long,
    covars_list        = heap$covars_list,
    split_df           = heap$split_df,
    instances          = heap$meta$instances,
    canonical_instance = heap$meta$canonical_instance
  )
}
