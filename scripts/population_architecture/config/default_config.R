
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

# Load centralized covariate config helpers (reads covariate_sets.yml)
local({
  helpers_path <- file.path(HEAP_PATHS$workflow, "config_helpers.R")
  if (file.exists(helpers_path)) source(helpers_path)
})

# Helper: build population_architecture covariate_spec entry from covariate_sets.yml
.pa_covar_spec <- function(name) {
  tryCatch({
    list(
      discrete           = load_covariate_set_discrete(name),
      quantitative       = load_covariate_set_kernel(name),  # uses 'quantitative' sub-list
      kernel_quantitative = load_covariate_set_kernel(name)
    )
  }, error = function(e) {
    message("population_architecture: could not load covariate set '", name,
            "' from YAML (", conditionMessage(e), "); using inline fallback.")
    NULL
  })
}

# For population_architecture, the 'quantitative' and 'kernel_quantitative' sub-lists
# differ: quantitative includes all numeric covariates (with PCs), kernel_quantitative
# excludes PCs (used for environment kernel construction). Load both from YAML.
.pa_load_covar_spec <- function(name) {
  .cfg_path <- heap_config("covariates", "covariate_sets.yml")
  if (!file.exists(.cfg_path) || !requireNamespace("yaml", quietly = TRUE)) {
    return(NULL)
  }
  tryCatch({
    .sets <- yaml::read_yaml(.cfg_path)$covariate_sets
    if (!name %in% names(.sets)) return(NULL)
    .entry <- .sets[[name]]
    list(
      discrete            = as.character(.entry$discrete            %||% character(0)),
      quantitative        = as.character(.entry$quantitative        %||% character(0)),
      kernel_quantitative = as.character(.entry$kernel_quantitative %||% character(0))
    )
  }, error = function(e) NULL)
}

`%||%` <- function(x, y) if (!is.null(x)) x else y
cfg <- list(
  project_root = HEAP_PATHS$heap_root,
  loader_rds = heap_loader_rds,
  omicpred_map = heap_omicspred_or_legacy("UKB_Olink_multi_ancestry_models_val_results_portal.csv"),
  gs_cis_dir = igloo_path("UKB", "ProtGScis"),
  gs_tr_dir = igloo_path("UKB", "ProtGStrans"),
  genotype_pfile = igloo_path("UKB", "gwas", "UKBallchr"),
  ld_pruned_pvar = igloo_path("UKB", "gwas", "ukb_nonimputed_snps.pvar"),
  gcta_bin = heap_gcta_bin(),   # IGLOO/GCTA/gcta64 (legacy fallback)
  plink2_bin = igloo_path("UKB", "gwas", "plink2"),
  # Module 1 OOF scores: read from IGLOO canonical location (written by manifest runs)
  # or fall back to local heap_output for legacy runs.
  predictive_oof_root = Sys.getenv(
    "HEAP_MODULE1_OOF_ROOT",
    unset = heap_project_output("module1_predictive_r2_score_partition")
  ),
  # Population architecture outputs: IGLOO canonical location.
  # Override via POPARCH_OUTPUT_ROOT env var (used by ProtGremlGrmCutoff.sh).
  output_root = Sys.getenv(
    "POPARCH_OUTPUT_ROOT",
    unset = heap_project_output("population_architecture")
  ),
  threads = 8L,
  plink_memory_mb = 96000L,
  block_size = 500L,
  min_protein_n = 2000L,
  exposure_imputation = "mean",
  exposure_missing_rate_max = 0.2,
  known_factor_vars = c(
    "uk_biobank_assessment_centre_f54_0_0",
    "sex_f31_0_0"
  ),
  pilot_proteins = c("DKKL1", "TLR3", "LILRB5"),
  # Covariate specs: loaded from config/covariates/covariate_sets.yml.
  # Inline fallback used if YAML is unavailable.
  covariate_specs = local({
    .base          <- .pa_load_covar_spec("base")
    .base_clinical <- .pa_load_covar_spec("base_clinical")

    # Inline fallback (used if yaml package is absent or covariate_sets.yml missing).
    # These MUST mirror covariate_sets.yml exactly.
    # NOTE: base drops BMI + fasting (unlike old Type3).
    .inline_base <- list(
      discrete = c("sex_f31_0_0", "uk_biobank_assessment_centre_f54_0_0"),
      quantitative = c(
        "age_when_attended_assessment_centre_f21003_0_0",
        "age2", "age_sex", "age2_sex",
        paste0("genetic_principal_components_f22009_0_", 1:20)
      ),
      kernel_quantitative = c(
        "age_when_attended_assessment_centre_f21003_0_0",
        "sex_numeric", "age2", "age_sex", "age2_sex"
      )
    )
    .inline_base_clinical <- list(
      discrete = c(
        "sex_f31_0_0", "uk_biobank_assessment_centre_f54_0_0",
        "assessment_season"
      ),
      quantitative = c(
        "age_when_attended_assessment_centre_f21003_0_0",
        "age2", "age_sex", "age2_sex",
        paste0("genetic_principal_components_f22009_0_", 1:20),
        "body_mass_index_bmi_f23104_0_0", "fasting_time_f74_0_0",
        "combined_Blood_pressure_medication", "combined_Hormone_replacement_therapy",
        "combined_Oral_contraceptive_pill_or_minipill", "combined_Insulin",
        "combined_Cholesterol_lowering_medication"
      ),
      kernel_quantitative = c(
        "age_when_attended_assessment_centre_f21003_0_0",
        "sex_numeric", "age2", "age_sex", "age2_sex",
        "body_mass_index_bmi_f23104_0_0", "fasting_time_f74_0_0"
      )
    )

    list(
      base          = .base          %||% .inline_base,
      base_clinical = .base_clinical %||% .inline_base_clinical
    )
  })
)
