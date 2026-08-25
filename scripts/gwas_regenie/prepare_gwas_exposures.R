#!/usr/bin/env Rscript
# prepare_gwas_exposures.R
#
# Prepares per-exposure phenotype and covariate files for REGENIE exposure GWAS.
# Reads from the unified HEAP.rds loader and the canonical analysis_exposures.tsv.
#
# Two-sample MR design:
#   The exposure GWAS sample is defined as UKB participants who have
#   genotype data AND exposure data, EXCLUDING the proteomics cohort.
#   This ensures that exposure GWAS summary statistics used in downstream
#   MR analyses come from a non-overlapping sample relative to the
#   proteomics-based protein association analyses.
#
# Exposure eligibility (from analysis_exposures.tsv):
#   include == 1  AND  miss_rate_prot_i0 < 0.20
#
# REGENIE options applied:
#   Quantitative (continuous + ordinal):  --apply-rint in BOTH step 1 and step 2
#   Binary:                               --bt (steps 1 and 2) + Firth in step 2
#
# Outputs:
#   Single-exposure mode (--exposure-id, invoked per task by the REGENIE array jobs):
#     <GWAS_EXPOSURE_INPUT_DIR>/<exposure>/pheno.txt   -- FID IID <exposure>
#     <GWAS_EXPOSURE_INPUT_DIR>/<exposure>/covar.txt   -- FID IID <covars>
#   Batch mode (default, no --exposure-id):
#     <HEAP_SLURM_GWAS>/evars_continuous_heap.txt      -- one exposure per line
#     <HEAP_SLURM_GWAS>/evars_binary_heap.txt
#     <HEAP_OUTPUT_GWAS>/sample_counts.tsv
#     <HEAP_OUTPUT_GWAS>/exposure_gwas_qc.tsv
#   Batch mode does NOT stage per-exposure pheno/covar files; each array job
#   regenerates its own in single-exposure mode.

############################################################
# 0) Bootstrap paths
############################################################

# Group-writable outputs for hpc_patel team runs: umask 0002 -> files 664, dirs
# 775. Set here (not just in the caller) so every file/dir this script creates --
# evar lists, QC tables, pheno/covar -- is writable by any group member, no matter
# how it is launched (submit_gwas_exposures.sh srun, a direct Rscript, or the
# per-job single-exposure call inside the REGENIE array scripts).
Sys.umask("0002")

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
  if (is.na(hit)) stop("Cannot find 00_paths.R. Set HEAP_ROOT or run from inside the HEAP tree.")
  source(hit)
})

# Config helpers: load_covariate_set() so the GWAS adjustment uses the canonical
# `base` covariate set (single source of truth) rather than a hardcoded list.
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
  library(purrr)
  library(readr)
  library(tidyr)
})

ts_msg <- function(...) message(format(Sys.time(), "[%H:%M:%S]"), " ", ...)

############################################################
# 1) Configuration and CLI
#
# Two modes:
#
#   Batch mode (default, no args):
#     Generates pheno/covar files for ALL eligible exposures.
#     Also writes evars_*.txt lists to slurm/gwas_regenie/.
#     Historically used GWAS_EXPOSURE_INPUT_DIR on scratch; now
#     accepts --output-dir to put files in a caller-supplied directory.
#
#   Single-exposure mode (used by regenie Slurm jobs):
#     Rscript prepare_gwas_exposures.R \
#       --exposure-id <name> \
#       --output-dir <per-job-tmp-dir>
#     Generates pheno/covar only for the named exposure.
#     Does NOT write evars list files (batch responsibility only).
#
############################################################

.parse_cli_args <- function() {
  raw <- commandArgs(trailingOnly = TRUE)
  out <- list(
    exposure_id  = NULL,   # single-exposure mode if non-NULL
    output_dir   = NULL,   # explicit output directory
    batch        = TRUE    # TRUE = generate all; FALSE = single exposure
  )
  i <- 1L
  while (i <= length(raw)) {
    flag <- raw[[i]]
    val  <- if (i < length(raw)) raw[[i + 1L]] else NA_character_
    if (flag == "--exposure-id") {
      out$exposure_id <- val; out$batch <- FALSE; i <- i + 2L
    } else if (flag == "--output-dir") {
      out$output_dir <- val; i <- i + 2L
    } else {
      i <- i + 1L
    }
  }
  out
}
CLI <- .parse_cli_args()

# Directory where per-exposure pheno/covar files are written.
# Single-exposure mode:  caller supplies --output-dir (job-local tmp directory).
# Batch mode:            use GWAS_EXPOSURE_INPUT_DIR env var or scratch default.
GWAS_EXPOSURE_INPUT_DIR <- if (!is.null(CLI$output_dir)) {
  CLI$output_dir
} else {
  Sys.getenv("GWAS_EXPOSURE_INPUT_DIR",
             unset = scratch_path("GWAS_exposure_input"))
}

# HEAP output directory for QC tables and logs (batch mode only) — IGLOO canonical.
HEAP_OUTPUT_GWAS <- heap_project_output("gwas_regenie")

# HEAP slurm directory for evar list files (batch mode only)
HEAP_SLURM_GWAS <- heap_path("slurm", "gwas_regenie")

# Missingness threshold: exposures with miss_rate >= this value are excluded
MISS_RATE_THRESHOLD <- 0.20

# GWAS covariates: use the canonical `base` covariate set from
# config/covariates/covariate_sets.yml as the single source of truth, so the exposure
# GWAS adjustment matches the module `base` adjustment.
#   base = age + sex + age^2 + age*sex + age^2*sex + assessment centre + 20 PCs.
# sex is hand-encoded 0/1 (identical to its dummy coding) and assessment centre (a
# 22-level factor) is one-hot encoded in section 5, so REGENIE can treat EVERY covar.txt
# column as quantitative -- no --catCovarList needed (consistent encoding for all
# categoricals). The final column list (with centre dummies) is built in section 5.
GWAS_COVAR_SET  <- Sys.getenv("GWAS_COVAR_SET", unset = "base")
GWAS_CENTRE_COL <- "uk_biobank_assessment_centre_f54_0_0"

# Only single-exposure mode writes under GWAS_EXPOSURE_INPUT_DIR; in batch mode this
# would just leave an empty scratch tree, so gate it. HEAP_OUTPUT_GWAS holds the QC
# tables written in batch mode, so it is always created.
if (!CLI$batch) dir.create(GWAS_EXPOSURE_INPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(HEAP_OUTPUT_GWAS,        recursive = TRUE, showWarnings = FALSE)

############################################################
# 2) Load HEAP.rds
############################################################

ts_msg("Loading HEAP.rds from: ", heap_loader_rds)
if (!file.exists(heap_loader_rds)) {
  stop("HEAP.rds not found at: ", heap_loader_rds,
       "\nRun HEAP_loader.R first to generate it.")
}
heap <- readRDS(heap_loader_rds)
ts_msg("HEAP.rds loaded. Baseline participants: ",
       nrow(heap$prot_baseline), " (proteomics cohort)")

############################################################
# 3) Load and filter analysis_exposures.tsv
############################################################

exposure_config_path <- heap_analysis_config()
ts_msg("Reading exposure config: ", exposure_config_path)
if (!file.exists(exposure_config_path)) {
  stop("analysis_exposures.tsv not found at: ", exposure_config_path)
}
exp_cfg <- read_tsv(exposure_config_path, col_types = cols(.default = "c")) %>%
  mutate(
    include       = as.integer(include),
    miss_rate_prot_i0 = as.numeric(miss_rate_prot_i0)
  )

# Apply eligibility filters
eligible <- exp_cfg %>%
  filter(include == 1L, miss_rate_prot_i0 < MISS_RATE_THRESHOLD)

ts_msg("Eligible exposures after include==1 and miss_rate < ", MISS_RATE_THRESHOLD, ": ",
       nrow(eligible), " (out of ", nrow(exp_cfg), " total in TSV)")

if (nrow(eligible) == 0) {
  stop("No eligible exposures remain after filtering. ",
       "Check analysis_exposures.tsv and the MISS_RATE_THRESHOLD.")
}

# Single-exposure mode: restrict to the requested exposure
if (!CLI$batch && !is.null(CLI$exposure_id)) {
  target <- CLI$exposure_id
  eligible <- eligible %>% filter(variable == target)
  if (nrow(eligible) == 0) {
    stop("Exposure '", target, "' not found in eligible exposures.\n",
         "Check that include==1 and miss_rate < ", MISS_RATE_THRESHOLD,
         " in analysis_exposures.tsv for this exposure.")
  }
  ts_msg("Single-exposure mode: ", target)
}

############################################################
# 4) Collect exposure data from HEAP baseline
############################################################

ts_msg("Stacking HEAP baseline exposure categories")
E_all <- purrr::reduce(heap$E_baseline, dplyr::full_join, by = "eid")
ts_msg("Total UKB participants with any exposure data: ", nrow(E_all))

# Verify that all eligible exposures are present in the stacked data
eligible_vars  <- eligible$variable
present_vars   <- intersect(eligible_vars, names(E_all))
missing_vars   <- setdiff(eligible_vars, names(E_all))

if (length(missing_vars) > 0) {
  missing_msg <- paste(
    "The following exposures in analysis_exposures.tsv are not present in HEAP.rds:",
    paste(missing_vars, collapse = ", ")
  )
  # Write a warning file and stop -- do not silently proceed with a partial exposure list
  warn_file <- file.path(HEAP_OUTPUT_GWAS, "MISSING_EXPOSURES_ERROR.txt")
  writeLines(c(
    paste("ERROR generated at:", Sys.time()),
    "Exposures listed in analysis_exposures.tsv (include==1, miss_rate<0.20)",
    "but NOT found in HEAP.rds E_baseline:",
    missing_vars
  ), warn_file)
  stop(missing_msg,
       "\nCheck that HEAP.rds was built from the same HEAP_loader.R version.",
       "\nMissing variable list written to: ", warn_file)
}
ts_msg("All ", length(eligible_vars), " eligible exposure variables found in HEAP.rds")

############################################################
# 5) Covariates and sample selection
############################################################

covars_raw <- heap$covars_baseline

# Recode sex to 0/1 (Male=1, Female=0) for REGENIE
sex_col <- "sex_f31_0_0"
if (!sex_col %in% names(covars_raw)) {
  stop("sex_f31_0_0 not found in covars_baseline. Check HEAP.rds.")
}
covars_raw[[sex_col]] <- ifelse(as.character(covars_raw[[sex_col]]) == "Male", 1L, 0L)

# Resolve the GWAS covariates from the canonical `base` set.
base_covars <- load_covariate_set(GWAS_COVAR_SET)
if (is.null(base_covars) || length(base_covars) == 0L)
  stop("Covariate set '", GWAS_COVAR_SET, "' did not resolve to any columns.")
numeric_base    <- setdiff(base_covars, GWAS_CENTRE_COL)   # all-numeric base covars (sex is 0/1)
missing_numeric <- setdiff(numeric_base, names(covars_raw))
if (length(missing_numeric) > 0)
  stop("Required base covariate columns not found in HEAP.rds covars_baseline:\n",
       paste(missing_numeric, collapse = ", "))

# One-hot encode assessment centre (22-level factor) into dummy columns, dropping one
# reference level to avoid collinearity with REGENIE's intercept. NA centre -> NA dummies
# (those rows are dropped by the complete-case filter below). This mirrors sex being
# hand-encoded 0/1, so REGENIE treats every covar.txt column as quantitative (no
# --catCovarList). Reference level is encoded as all-zero dummies.
centre_dummy_cols <- character(0)
if (GWAS_CENTRE_COL %in% base_covars) {
  if (!GWAS_CENTRE_COL %in% names(covars_raw))
    stop("base includes ", GWAS_CENTRE_COL, " but it is absent from covars_baseline.")
  cen_chr    <- as.character(covars_raw[[GWAS_CENTRE_COL]])
  cen_levels <- sort(unique(cen_chr[!is.na(cen_chr)]))
  if (length(cen_levels) < 2L)
    stop("Assessment centre has <2 non-missing levels; cannot one-hot encode.")
  ref_level    <- cen_levels[1L]
  dummy_levels <- cen_levels[-1L]
  for (lv in dummy_levels) {
    col <- paste0("centre_", make.names(lv))
    covars_raw[[col]] <- ifelse(is.na(cen_chr), NA_integer_, as.integer(cen_chr == lv))
    centre_dummy_cols <- c(centre_dummy_cols, col)
  }
  ts_msg("One-hot encoded assessment centre: ", length(dummy_levels),
         " dummies (reference = ", ref_level, ")")
}

# Final GWAS covariate column list: numeric base covars + centre dummies.
GWAS_COVAR_COLS <- c(numeric_base, centre_dummy_cols)
ts_msg("GWAS covariate set '", GWAS_COVAR_SET, "': ", length(GWAS_COVAR_COLS),
       " columns (", length(numeric_base), " numeric + ",
       length(centre_dummy_cols), " centre dummies)")

# Check all required covariate columns are present
missing_covars <- setdiff(GWAS_COVAR_COLS, names(covars_raw))
if (length(missing_covars) > 0) {
  stop("Required GWAS covariate columns not found in HEAP.rds covars_baseline:\n",
       paste(missing_covars, collapse = ", "))
}

# Keep only required covariate columns + eid
covars_sel <- covars_raw %>% select(eid, all_of(GWAS_COVAR_COLS))
covars <- covars_sel[rowSums(is.na(covars_sel)) == 0, ]   # REGENIE requires no missing covariates

ts_msg("Participants with complete covariates: ", nrow(covars))

# Proteomics cohort IDs (excluded for two-sample MR)
prot_ids <- heap$prot_baseline$eid
ts_msg("Proteomics cohort size (to be excluded): ", length(prot_ids))

# Non-proteomics participants with complete covariates
gwas_base <- covars %>%
  filter(!(eid %in% prot_ids))
ts_msg("Non-proteomics participants with complete covariates (GWAS base): ", nrow(gwas_base))

# Global sample counts (written to log)
sample_counts <- data.frame(
  metric = c(
    "total_ukb_with_any_exposure_data",
    "participants_with_complete_covariates",
    "proteomics_cohort_excluded",
    "gwas_base_sample_non_proteomics"
  ),
  n = c(
    nrow(E_all),
    nrow(covars),
    length(prot_ids),
    nrow(gwas_base)
  )
)

############################################################
# 6) Per-exposure file generation
############################################################

ts_msg("Generating per-exposure phenotype and covariate files")

qc_rows <- list()
evars_continuous <- character(0)
evars_binary     <- character(0)

for (i in seq_len(nrow(eligible))) {
  evar       <- eligible$variable[i]
  etype      <- eligible$variable_type[i]   # "continuous", "ordinal", or "binary"
  miss_rate  <- eligible$miss_rate_prot_i0[i]
  is_quant   <- etype %in% c("continuous", "ordinal")
  apply_rint <- is_quant

  # Extract this exposure for the GWAS base sample
  edata <- E_all %>%
    select(eid, all_of(evar)) %>%
    filter(eid %in% gwas_base$eid)

  n_with_phenotype  <- sum(!is.na(edata[[evar]]))
  n_missing_pheno   <- sum(is.na(edata[[evar]]))

  if (n_with_phenotype == 0) {
    ts_msg("WARNING: ", evar, " has 0 non-missing values in GWAS base sample -- skipping")
    next
  }

  # Complete cases only (drop participants with NA phenotype for this exposure)
  subset <- edata %>%
    filter(!is.na(.data[[evar]])) %>%
    left_join(gwas_base, by = "eid")

  n_final <- nrow(subset)
  n_cases   <- NA_integer_
  n_controls <- NA_integer_
  if (!is_quant) {
    n_cases    <- sum(subset[[evar]] == 1, na.rm = TRUE)
    n_controls <- sum(subset[[evar]] == 0, na.rm = TRUE)
  }

  # Per-exposure pheno/covar files are written ONLY in single-exposure mode.
  # The REGENIE array jobs invoke this script with --exposure-id and regenerate their
  # own pheno/covar fresh per task, so writing all ~172 pairs in batch mode is dead
  # weight (and heavy scratch I/O). The counting above (n_with_phenotype, n_final,
  # n_cases/n_controls) and the evar-list/QC bookkeeping below still run for every
  # exposure, so the evar lists and QC tables are unchanged by this gate.
  if (!CLI$batch) {
    # Output directory for this exposure
    out_dir <- file.path(GWAS_EXPOSURE_INPUT_DIR, evar)
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    # Phenotype file: FID IID <trait>
    pheno_df <- subset %>%
      transmute(
        FID = eid,
        IID = eid,
        !!evar := .data[[evar]]
      )
    fwrite(pheno_df, file = file.path(out_dir, "pheno.txt"),
           sep = "\t", na = "NA", quote = FALSE)

    # Covariate file: FID IID <covars>
    covar_df <- subset %>%
      select(eid, all_of(GWAS_COVAR_COLS)) %>%
      mutate(FID = eid, IID = eid) %>%
      select(FID, IID, all_of(GWAS_COVAR_COLS)) %>%
      select(-any_of("eid"))
    fwrite(covar_df, file = file.path(out_dir, "covar.txt"),
           sep = "\t", na = "NA", quote = FALSE)
  }

  # Track evar lists
  if (is_quant) {
    evars_continuous <- c(evars_continuous, evar)
  } else {
    evars_binary <- c(evars_binary, evar)
  }

  # QC row
  qc_rows[[length(qc_rows) + 1L]] <- data.frame(
    exposure_id              = evar,
    category                 = eligible$category[i],
    exposure_type            = etype,
    miss_rate_from_tsv       = miss_rate,
    n_gwas_base              = nrow(gwas_base),
    n_with_phenotype         = n_with_phenotype,
    n_missing_phenotype      = n_missing_pheno,
    n_proteomics_excluded    = length(prot_ids),
    n_final_gwas             = n_final,
    apply_rint               = apply_rint,
    n_cases_binary           = n_cases,
    n_controls_binary        = n_controls,
    stringsAsFactors         = FALSE
  )

  if (i %% 10 == 0) ts_msg("  Progress: ", i, " / ", nrow(eligible))
}

ts_msg("Done: ", length(evars_continuous), " quantitative exposures, ",
       length(evars_binary), " binary exposures")

############################################################
# 7) Write evar list files for SLURM arrays (batch mode only)
############################################################

if (CLI$batch) {
  evars_cont_path <- file.path(HEAP_SLURM_GWAS, "evars_continuous_heap.txt")
  evars_bin_path  <- file.path(HEAP_SLURM_GWAS, "evars_binary_heap.txt")

  writeLines(evars_continuous, evars_cont_path)
  writeLines(evars_binary,     evars_bin_path)

  ts_msg("Written: ", evars_cont_path, " (", length(evars_continuous), " exposures)")
  ts_msg("Written: ", evars_bin_path,  " (", length(evars_binary),     " exposures)")
} else {
  ts_msg("Single-exposure mode: skipping evars list file update.")
}

############################################################
# 8) Write QC and sample count tables (batch mode only)
############################################################

if (CLI$batch) {
  qc_table <- bind_rows(qc_rows)
  qc_path  <- file.path(HEAP_OUTPUT_GWAS, "exposure_gwas_qc.tsv")
  write_tsv(qc_table, qc_path)
  ts_msg("Written QC table: ", qc_path)

  sc_path <- file.path(HEAP_OUTPUT_GWAS, "sample_counts.tsv")
  write_tsv(sample_counts, sc_path)
  ts_msg("Written sample counts: ", sc_path)
}

############################################################
# 9) Summary
############################################################

ts_msg("=== GWAS EXPOSURE PREP SUMMARY ===")
ts_msg("Mode:                           ", if (CLI$batch) "batch (all exposures)" else paste0("single-exposure (", CLI$exposure_id, ")"))
ts_msg("HEAP.rds path:                  ", heap_loader_rds)
ts_msg("analysis_exposures.tsv:         ", exposure_config_path)
ts_msg("Missingness threshold:          <", MISS_RATE_THRESHOLD)
ts_msg("Exposures processed:            ", nrow(eligible))
ts_msg("  Quantitative (cont+ordinal):  ", length(evars_continuous))
ts_msg("  Binary:                       ", length(evars_binary))
ts_msg("Proteomics IDs excluded:        ", length(prot_ids))
ts_msg("GWAS base sample size:          ", nrow(gwas_base))
if (!CLI$batch) ts_msg("Pheno/covar output dir:         ", GWAS_EXPOSURE_INPUT_DIR)
if (CLI$batch) {
  ts_msg("SLURM array sizes:")
  ts_msg("  --array=1-", length(evars_continuous), "  (continuous script)")
  ts_msg("  --array=1-", length(evars_binary),     "  (binary script)")
  ts_msg("QC table:                       ", file.path(HEAP_OUTPUT_GWAS, "exposure_gwas_qc.tsv"))
  ts_msg("Sample counts:                  ", file.path(HEAP_OUTPUT_GWAS, "sample_counts.tsv"))
}
ts_msg("=== DONE ===")
