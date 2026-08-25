#!/usr/bin/env Rscript

# Central path configuration for the HEAP reproducible workflow.
# Scripts copied into HEAP should source this file and write new analysis
# outputs under HEAP/output while reading protected legacy inputs in place.

heap_find_root <- function(start = getwd()) {
  # Explicit override wins FIRST. Required for git worktrees: their leaf dir is
  # named after the agent (not "HEAP"), so the auto-walk below would otherwise
  # climb past them to the canonical main checkout and ignore the worktree's code.
  env_root <- Sys.getenv("HEAP_ROOT", unset = "")
  if (nzchar(env_root)) return(normalizePath(env_root, mustWork = FALSE))
  path <- normalizePath(start, mustWork = FALSE)
  repeat {
    if (basename(path) == "HEAP" && dir.exists(file.path(path, "workflow"))) {
      return(path)
    }
    candidate <- file.path(path, "HEAP")
    if (dir.exists(file.path(candidate, "workflow"))) {
      return(normalizePath(candidate, mustWork = FALSE))
    }
    parent <- dirname(path)
    if (identical(parent, path)) break
    path <- parent
  }
  stop("Could not locate HEAP root. Set HEAP_ROOT or run from inside the HEAP workflow.")
}

heap_root <- heap_find_root()
workspace_root <- normalizePath(file.path(heap_root, ".."), mustWork = FALSE)

# Shared-group outputs: create files/dirs group-writable so any hpc_patel member
# can co-write the canonical IGLOO HEAP outputs. umask 002 => dirs rwxrwxr-x,
# files rw-rw-r--. Honors a stricter umask already set in the environment via
# HEAP_UMASK (e.g. "077") if a user wants private outputs.
Sys.umask(Sys.getenv("HEAP_UMASK", unset = "002"))

HEAP_PATHS <- list(
  heap_root = heap_root,
  workflow = file.path(heap_root, "workflow"),
  scripts = file.path(heap_root, "scripts"),
  slurm = file.path(heap_root, "slurm"),
  config = file.path(heap_root, "config"),
  output = file.path(heap_root, "output"),
  logs = file.path(heap_root, "logs"),
  legacy_ukb_root = Sys.getenv(
    "HEAP_LEGACY_UKB_ROOT",
    unset = file.path(workspace_root, "UK_Biobank")
  ),
  # Per-user scratch (O2 convention: /n/scratch/users/<first-letter>/<user>).
  # Derived from $USER so any group member's jobs stage to their own scratch;
  # override with HEAP_SCRATCH_ROOT.
  scratch_root = Sys.getenv(
    "HEAP_SCRATCH_ROOT",
    unset = local({
      u <- Sys.getenv("USER", unset = "unknown")
      file.path("/n/scratch/users", substr(u, 1, 1), u)
    })
  ),
  igloo_root = Sys.getenv(
    "HEAP_IGLOO_ROOT",
    unset = "/n/groups/patel/IGLOO"
  )
)

# ---------------------------------------------------------------------------
# Shared group R library
#
# So any hpc_patel member can run the HEAP R scripts without installing the
# package set themselves. When present, the shared IGLOO R library (versioned by
# R major.minor.patch) is prepended to .libPaths(); base/system libs still apply.
# Build/refresh it with scripts/setup/install_r_packages.R. Override with
# HEAP_RLIB, or set HEAP_RLIB="" to disable and use only your personal library.
# ---------------------------------------------------------------------------
local({
  rlib <- Sys.getenv(
    "HEAP_RLIB",
    unset = file.path(HEAP_PATHS$igloo_root, "Rlib", as.character(getRversion()))
  )
  if (nzchar(rlib) && dir.exists(rlib)) .libPaths(c(rlib, .libPaths()))
})

heap_path <- function(...) file.path(HEAP_PATHS$heap_root, ...)
heap_config <- function(...) file.path(HEAP_PATHS$config, ...)
heap_output <- function(...) file.path(HEAP_PATHS$output, ...)
heap_script <- function(...) file.path(HEAP_PATHS$scripts, ...)
heap_log <- function(...) file.path(HEAP_PATHS$logs, ...)
legacy_ukb_path <- function(...) file.path(HEAP_PATHS$legacy_ukb_root, ...)
scratch_path <- function(...) file.path(HEAP_PATHS$scratch_root, ...)
igloo_path <- function(...) file.path(HEAP_PATHS$igloo_root, ...)

# Protein list: prefer versioned copy in HEAP/config, fall back to legacy location.
heap_omicspred_protein_list <- local({
  config_copy <- heap_config("protein_sets", "omicspred_proteins.txt")
  if (file.exists(config_copy)) config_copy else legacy_ukb_path("BScripts", "ProtPGS_PXS", "OMICPREDproteins.txt")
})

# heap_loader_rds is defined below, after the IGLOO-rooted helpers are set up.
# Do not use scratch for heap_loader_rds in new scripts.
# The canonical path is heap_project_intermediate("HEAP.rds") (IGLOO-rooted).
# See lines below for the final definition.

# ---------------------------------------------------------------------------
# Exposure variable TYPES — single source of truth.
#
# The declared type of every exposure lives in the `variable_type` column of
# config/exposure_sets/analysis_exposures.tsv (binary | ordinal | continuous).
# Modules must route encoding by this DECLARED type, NOT by runtime value-range
# heuristics. The loader's ordinal_finder() ("ordinal if max value <= 5") and the
# per-module continuous_finder() ("continuous if max>5 | unique>2") mis-classify
# small-scale CONTINUOUS scores — e.g. the England IMD income/employment/health/
# crime scores and pm2.5 absorbance (all in [0, 5]) — as ordinal, so they get
# factorised into one term per distinct value (crime_score: ~471 spurious
# polynomial-contrast terms with ~1e13 coefficients). HEAP.rds stores every
# exposure as plain numeric; only the type ROUTING is wrong.
#
# heap_exposure_type_map() returns a named char vector variable -> type (or NULL
# when the config is unavailable, so callers can fall back to legacy behaviour).
# heap_exposures_of_type() filters it, optionally to columns actually present.
# ---------------------------------------------------------------------------
heap_exposure_type_map <- function() {
  f <- tryCatch(heap_config("exposure_sets", "analysis_exposures.tsv"),
                error = function(e) "")
  if (!nzchar(f) || !file.exists(f)) return(NULL)
  ae <- tryCatch(utils::read.delim(f, stringsAsFactors = FALSE, check.names = FALSE),
                 error = function(e) NULL)
  if (is.null(ae) || !all(c("variable", "variable_type") %in% names(ae))) return(NULL)
  if ("include" %in% names(ae))
    ae <- ae[is.na(ae$include) | ae$include == 1, , drop = FALSE]
  stats::setNames(as.character(ae$variable_type), as.character(ae$variable))
}

#' Declared exposures of the given type(s); NULL if the config is unavailable.
#' @param types one or more of "binary","ordinal","continuous"
#' @param present optional vector of column names to intersect with
heap_exposures_of_type <- function(types, present = NULL) {
  m <- heap_exposure_type_map()
  if (is.null(m)) return(NULL)
  v <- names(m)[m %in% types]
  if (!is.null(present)) v <- intersect(v, present)
  v
}

#' ordinalIDs derived from the declared types (the variables Module 1/2/3 should
#' factorise). Falls back to the loader's heuristic ordinalIDs when the config is
#' unavailable. `present` = exposure columns actually in the loader object.
heap_resolve_ordinal_ids <- function(heap, present) {
  cfg <- heap_exposures_of_type("ordinal", present)
  if (is.null(cfg)) heap$ordinalIDs else cfg
}

.heap_present_exposures <- function(Elist)
  unique(unlist(lapply(Elist, function(d) setdiff(names(d), "eid")), use.names = FALSE))

#' Names of ZERO-VARIANCE ("constant") columns among `cols` — features that
#' contribute nothing to a model and must be dropped: numeric columns that are
#' all-NA or have < 2 distinct non-NA values, and factors with < 2 non-empty
#' levels. Used by Module 1 & Module 2 to exclude constant exposures/covariates
#' (globally constant ones, e.g. former_alcohol_drinker_f3731, and any that
#' become constant within a train/test split).
heap_zero_variance_cols <- function(df, cols = names(df)) {
  cols <- intersect(cols, names(df))
  if (!length(cols)) return(character(0))
  is_zv <- function(x) {
    xx <- x[!is.na(x)]
    if (length(xx) == 0L) return(TRUE)
    if (is.factor(x)) return(nlevels(droplevels(factor(xx))) < 2L)
    length(unique(xx)) < 2L
  }
  cols[vapply(cols, function(c) is_zv(df[[c]]), logical(1))]
}

# Compatibility accessors for downstream module scripts.
# Returns a PXSconstruct-compatible list for cross-sectional modules (1, 2, 3, 5).
# ordinalIDs is now routed by the DECLARED variable_type (see above), so the
# small-scale continuous scores are no longer mis-factorised.
as_pxs_baseline <- function(heap) {
  list(
    Elist       = heap$E_baseline,
    Elist_names = heap$Elist_names,
    Eid_cat     = heap$Eid_cat,
    ordinalIDs  = heap_resolve_ordinal_ids(heap, .heap_present_exposures(heap$E_baseline)),
    UKBprot_df  = heap$prot_baseline,
    protIDs     = heap$protIDs,
    covars_df   = heap$covars_baseline,
    covars_list = heap$covars_list
  )
}

# Returns a longitudinal-compatible list for Module 6.
as_pxs_longitudinal <- function(heap) {
  list(
    Elist              = heap$E_long,
    Elist_names        = heap$Elist_names,
    Eid_cat            = heap$Eid_cat,
    ordinalIDs         = heap_resolve_ordinal_ids(heap, .heap_present_exposures(heap$E_long)),
    UKBprot_df         = heap$prot_long,
    protIDs            = heap$protIDs,
    covars_df          = heap$covars_long,
    covars_list        = heap$covars_list,
    split_df           = heap$split_df,
    instances          = heap$meta$instances,
    canonical_instance = heap$meta$canonical_instance
  )
}

dir.create(HEAP_PATHS$output, recursive = TRUE, showWarnings = FALSE)
dir.create(HEAP_PATHS$logs, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# IGLOO-rooted HEAP project paths
#
# Canonical generated HEAP outputs, manifests, and run artifacts live under:
#   /n/groups/patel/IGLOO/UKB/HEAP   (heap_project_root)
#
# These helpers are the authoritative paths for all large generated files that
# need to be reproducible and auditable across runs. The local HEAP/output
# folder is reserved for lightweight, git-ignored staging; large canonical
# outputs must go through heap_project_output() / heap_project_intermediate().
# ---------------------------------------------------------------------------

HEAP_PATHS$project_root <- Sys.getenv(
  "HEAP_PROJECT_ROOT",
  unset = file.path(HEAP_PATHS$igloo_root, "UKB", "HEAP")
)

heap_project_root <- function(...) file.path(HEAP_PATHS$project_root, ...)

# Generated analysis results (module outputs, score files, mediation tables)
heap_project_output <- function(...) heap_project_root("output", ...)

# Manifest TSV files (one per experiment; Slurm reads these)
heap_manifest <- function(...) heap_project_root("manifests", ...)

# Large intermediate files that are reproducible but not results
# (e.g., HEAP.rds, longitudinal RDS files)
heap_project_intermediate <- function(...) heap_project_root("intermediate", ...)

# Slurm logs for jobs writing to IGLOO-rooted outputs
heap_project_logs <- function(...) heap_project_root("logs", ...)

# Run configuration artifacts written by each job before it runs
heap_run_config <- function(...) heap_project_root("run_configs", ...)

# Documentation/reports generated by analysis runs
heap_project_docs <- function(...) heap_project_root("docs", ...)

# HEAP-specific GWAS outputs (exposure GWAS, regenie step2 final summstats).
# Canonical: /n/groups/patel/IGLOO/UKB/HEAP/output/gwas/...
# Distinct from shared genetics resources under igloo_path("UKB", "gwas").
heap_gwas <- function(...) heap_project_output("gwas", ...)

# Shared IGLOO genetics resources (read-only references, not generated by HEAP).
# igloo_path("UKB", "gwas")         -> shared genotype pfiles, LD reference
# igloo_path("UKB", "ProtGScis")    -> cis protein genetic scores (IGLOO canonical)
# igloo_path("UKB", "ProtGStrans")  -> trans protein genetic scores (IGLOO canonical)
# igloo_path("UKB", "pQTL")         -> UKB protein GWAS summary stats
# igloo_path("UKB", "pQTLmetadata") -> pQTL metadata files
# igloo_path("FinnGen", "SummaryStats") -> FinnGen disease GWAS summary stats
# igloo_path("DECODE", "pQTL", ...) -> deCODE pQTL summary stats

# ---------------------------------------------------------------------------
# IGLOO RAW + shared-resource helpers (migration targets)
#
# These centralize the canonical IGLOO locations for raw UKB inputs and shared
# reference resources, with a transparent fall-back to the legacy UK_Biobank
# tree so scripts keep working before/after each data copy. Once a file exists
# at the IGLOO location it is used automatically; until then the legacy copy is
# used. This makes the legacy->IGLOO migration safe and incremental.
# ---------------------------------------------------------------------------

# Converted raw UKB parquet/feather + path indices + codings.
# Canonical: /n/groups/patel/IGLOO/UKB/RAW
heap_raw <- function(...) igloo_path("UKB", "RAW", ...)

# Resolve a raw input preferring IGLOO RAW, else a legacy full path.
# igloo_rel: character vector of path parts under RAW (e.g. c("codings","Codings.csv")).
# legacy_full: absolute legacy path used if the IGLOO copy is not present yet.
# Returns the IGLOO path when neither exists, so errors point at the canonical home.
heap_raw_or_legacy <- function(igloo_rel, legacy_full) {
  ig <- do.call(heap_raw, as.list(igloo_rel))
  if (file.exists(ig)) ig else if (file.exists(legacy_full)) legacy_full else ig
}

# Shared IGLOO MR LD reference (1000G plink bfiles). Canonical: /n/groups/patel/IGLOO/LDref
# Used as a plink bfile *prefix* (e.g. EUR -> EUR.bed/bim/fam).
heap_ldref <- function(...) igloo_path("LDref", ...)
heap_ldref_prefix <- function(pop = "EUR",
                              legacy_full = legacy_ukb_path("RScripts", "Pure_StatGen",
                                                            "MR", "LDref", pop)) {
  ig <- heap_ldref(pop)
  if (file.exists(paste0(ig, ".bed"))) ig
  else if (file.exists(paste0(legacy_full, ".bed"))) legacy_full
  else ig
}

# Shared IGLOO GCTA install. Canonical: /n/groups/patel/IGLOO/GCTA
heap_gcta <- function(...) igloo_path("GCTA", ...)
heap_gcta_bin <- function(
    legacy_full = legacy_ukb_path("bin", "gcta-1.94.1-linux-kernel-3-x86_64", "gcta64")) {
  ig <- heap_gcta("gcta64")
  if (file.exists(ig)) ig else if (file.exists(legacy_full)) legacy_full else ig
}

# Shared IGLOO OmicsPred reference files. Canonical: /n/groups/patel/IGLOO/UKB/OMICSPRED
heap_omicspred <- function(...) igloo_path("UKB", "OMICSPRED", ...)
heap_omicspred_or_legacy <- function(rel) {
  ig <- do.call(heap_omicspred, as.list(rel))
  leg <- do.call(function(...) legacy_ukb_path("Data", "OMICSPRED", ...), as.list(rel))
  if (file.exists(ig)) ig else if (file.exists(leg)) leg else ig
}

# External intervention proteomics (GLP1 STEP1/STEP2, HERITAGE) used by the
# intervention-comparison support analysis. Canonical: /n/groups/patel/IGLOO/UKB/Interventions
#   GLP1_proteomics.xlsx  (STEP1=sheet S2_tx_STEP1, STEP2=sheet S3_tx_STEP2)
#   jciinsight_prot.xlsx  (HERITAGE; sheet 1, skip=2)
heap_interventions <- function(...) igloo_path("UKB", "Interventions", ...)
heap_interventions_or_legacy <- function(name) {
  ig <- heap_interventions(name)
  leg <- file.path("/n/groups/patel/shakson_ukb/Motrpac/Related_Data", name)
  if (file.exists(ig)) ig else if (file.exists(leg)) leg else ig
}
# Olink<->SomaScan cross-platform reliability (already IGLOO-canonical).
heap_olinksoma <- function(name = "OlinkSoma.csv") igloo_path("UKB", "OlinkSoma", name)

# GTEx v10 median-TPM GCT matrices (tissue gene-set construction).
# Canonical: /n/groups/patel/IGLOO/UKB/GTEX/RNA
heap_gtex <- function(...) igloo_path("UKB", "GTEX", ...)
heap_gtex_rna_dir <- function() {
  ig <- heap_gtex("RNA")
  if (dir.exists(ig)) ig else legacy_ukb_path("Data", "GTEX", "RNA")
}

# Human Protein Atlas tissue/subcellular tables. Canonical: /n/groups/patel/IGLOO/UKB/HPA
heap_hpa <- function(...) igloo_path("UKB", "HPA", ...)
heap_hpa_or_legacy <- function(name) {
  ig <- heap_hpa(name)
  leg <- legacy_ukb_path("Data", "HPA", name)
  if (file.exists(ig)) ig else if (file.exists(leg)) leg else ig
}

# Cold-storage genotype bgen (large, not backed up). Canonical:
# /n/no_backup2/patel/IGLOO/UKB/Genetics
heap_no_backup_igloo_root <- Sys.getenv(
  "HEAP_NO_BACKUP_IGLOO_ROOT", unset = "/n/no_backup2/patel/IGLOO")
heap_genetics_bgen <- function(name = "UKBallchr.bgen")
  file.path(heap_no_backup_igloo_root, "UKB", "Genetics", name)

# Canonical HEAP loader RDS — IGLOO-rooted.
# Production default: /n/groups/patel/IGLOO/UKB/HEAP/intermediate/HEAP.rds
# Override via HEAP_LOADER_RDS env var (e.g., to point at scratch during a
# transition period or for local testing only).
heap_loader_rds_canonical <- heap_project_intermediate("HEAP.rds")

heap_loader_rds <- local({
  env_override <- Sys.getenv("HEAP_LOADER_RDS", unset = "")
  if (nzchar(env_override)) {
    normalizePath(env_override, mustWork = FALSE)
  } else {
    heap_loader_rds_canonical
  }
})

# Create IGLOO-rooted directories on load so scripts can write immediately
invisible(lapply(
  list(
    HEAP_PATHS$project_root,
    heap_project_output(),
    heap_manifest(),
    heap_project_intermediate(),
    heap_project_logs(),
    heap_run_config()
  ),
  dir.create, recursive = TRUE, showWarnings = FALSE
))

# ---------------------------------------------------------------------------
# Exposure feature-selection config
#
# heap_filter_exposures(heap, config_path) applies analysis_exposures.tsv
# to a HEAP object, returning a modified copy with only include=1 variables
# retained in E_long, E_baseline, Eid_cat, and ordinalIDs.
#
# Usage in a module script:
#   heap <- readRDS(heap_loader_rds)
#   heap <- heap_filter_exposures(heap)      # uses default config path
# ---------------------------------------------------------------------------

heap_analysis_config <- function()
  heap_config("exposure_sets", "analysis_exposures.tsv")

heap_filter_exposures <- function(heap,
                                  config_path = heap_analysis_config()) {
  if (!file.exists(config_path))
    stop("Exposure config not found: ", config_path,
         "\nRun the config generator or check HEAP/config/exposure_sets/")

  cfg <- read.delim(config_path, stringsAsFactors = FALSE, check.names = FALSE)
  keep_vars <- cfg$variable[cfg$include == 1L]

  # Filter E_long and E_baseline for each category
  for (cat_name in names(heap$E_long)) {
    keep_cat <- keep_vars[keep_vars %in% names(heap$E_long[[cat_name]])]
    heap$E_long[[cat_name]]     <- heap$E_long[[cat_name]][,
      c("eid", "instance", keep_cat), drop = FALSE]
    if (!is.null(heap$E_baseline[[cat_name]])) {
      heap$E_baseline[[cat_name]] <- heap$E_baseline[[cat_name]][,
        c("eid", keep_cat), drop = FALSE]
    }
  }

  # Drop categories that end up with no exposure variables
  empty <- vapply(heap$E_long, function(df)
    length(setdiff(names(df), c("eid", "instance"))) == 0L, logical(1))
  heap$E_long     <- heap$E_long[!empty]
  heap$E_baseline <- heap$E_baseline[!empty]

  # Rebuild Eid_cat and ordinalIDs from filtered data
  if (!is.null(heap$Eid_cat)) {
    heap$Eid_cat <- heap$Eid_cat[heap$Eid_cat$Eid %in% keep_vars, ]
  }
  heap$ordinalIDs <- heap$ordinalIDs[heap$ordinalIDs %in% keep_vars]
  heap$Elist_names <- names(heap$E_long)

  heap
}
