
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
library(data.table)
library(tidyverse)

suppressPackageStartupMessages({
  library(data.table)
})

# Single source of truth: derive the module6 exposure list AND per-exposure type
# from config/exposure_sets/analysis_exposures.tsv (include == 1), so module6 is
# consistent with module6 prod, modules 1-3, and the exposure GWAS. The previous
# version read the stale standalone evars_binary.txt / evars_continuous.txt files.
#
# Type mapping mirrors the GWAS is_quant logic (prepare_gwas_exposures.R):
#   continuous + ordinal -> "continuous"  (quantitative; gaussian PES outcome)
#   binary               -> "binary"      (binomial PES outcome)
# Multi-category exposures are already one-hot-expanded into binary rows upstream,
# so there is no separate categorical/multinomial type to handle here.
exposure_cfg_path <- heap_analysis_config()
out_path          <- heap_project_output("module6_pes_test", "exposure_specs.tsv")  # IGLOO canonical

ecfg <- fread(exposure_cfg_path, sep = "\t")
ecfg <- ecfg[as.integer(include) == 1L]
ecfg[, exposure_id := trimws(as.character(variable))]
ecfg <- ecfg[nzchar(exposure_id)]

ecfg[, exposure_type := fifelse(
  variable_type %in% c("continuous", "ordinal"), "continuous",
  fifelse(variable_type == "binary", "binary", NA_character_)
)]
if (anyNA(ecfg$exposure_type)) {
  bad <- sort(unique(ecfg$variable_type[is.na(ecfg$exposure_type)]))
  stop("Unmapped variable_type(s) in analysis_exposures.tsv: ", paste(bad, collapse = ", "))
}

manifest <- unique(ecfg[, .(exposure_id, exposure_type)])
if (anyDuplicated(manifest$exposure_id)) {
  dups <- manifest$exposure_id[duplicated(manifest$exposure_id)]
  stop("ERROR: duplicate exposure_ids in manifest:\n", paste(unique(dups), collapse = "\n"))
}
setorder(manifest, exposure_id, exposure_type)

# write exposure_specs.tsv (exposure_id + exposure_type; read by the module6 runtime)
fwrite(manifest, out_path, sep = "\t")
cat("Wrote manifest:", out_path, "\n",
    "Source:", exposure_cfg_path, "(include==1)\n",
    "N total:",      nrow(manifest), "\n",
    "N binary:",     sum(manifest$exposure_type == "binary"), "\n",
    "N continuous:", sum(manifest$exposure_type == "continuous"),
    "(continuous + ordinal)\n")

# write the module6 exposure list (one exposure_id per line; used by compact/cox manifests)
fwrite(as.data.frame(manifest$exposure_id),
       heap_path("config", "exposure_sets", "module6_exposures.txt"),
       col.names = FALSE)







#' #exposure_id = "summed_met_minutes_per_week_for_all_activity_f22040_0_0" #"alcohol_intake_frequency_f1558_0_0" #"fresh_fruit_intake_f1309_0_0" #"types_of_physical_activity_in_last_4_weeks_f6164_0_0.multi_Strenuous_sports" #"smoking_status_f20116_0_0_Current"
#' #disease_age_col = "age_e11_first_reported_non_insulin_dependent_diabetes_mellitus_f130708_0_0" #"age_n18_first_reported_chronic_renal_failure_f132032_0_0" #"age_e11_first_reported_non_insulin_dependent_diabetes_mellitus_f130708_0_0" #"age_e11_first_reported_non_insulin_dependent_diabetes_mellitus_f130708_0_0" #"age_j43_first_reported_emphysema_f131490_0_0"
#' 
#' p1 <- readRDS(scratch_path("UKB_intermediate", "UKB_PGS_PXS_load.rds"))
#' #p2 <- readRDS("/n/groups/patel/IGLOO/UKB/Mediation/Data/UKB_MDstore_Type5.rds")
#' 
#' #Get exposure and disease lists:
#' 
#' # Exposures
#' as_pxs <- function(x) {
#'   if (is.list(x) && !isS4(x)) return(x)
#'   if (!isS4(x)) stop("PXS object must be S4 or list")
#'   list(
#'     Elist       = x@Elist,
#'     Elist_names = x@Elist_names,
#'     Eid_cat     = x@Eid_cat,
#'     ordinalIDs  = x@ordinalIDs,
#'     UKBprot_df  = x@UKBprot_df,
#'     protIDs     = x@protIDs,
#'     covars_df   = x@covars_df,
#'     covars_list = x@covars_list
#'   )
#' }
#' p1_x <- as_pxs(p1)
#' exposures <- p1_x$Eid_cat$Eid
#' View(as.data.frame(exposures))
#' 
#' #'*Make the manifest:*




# # Disease
# p1_dz <- p2@DZ_df
# disease <- colnames(p1_dz)[2:274]
# View(as.data.frame(disease))
# 
# 
# exposure_v1 <- c("alcohol_intake_frequency_f1558_0_0",
#                  "fresh_fruit_intake_f1309_0_0",
#                  "summed_met_minutes_per_week_for_all_activity_f22040_0_0",
#                  "types_of_physical_activity_in_last_4_weeks_f6164_0_0.multi_Strenuous_sports",
#                  "smoking_status_f20116_0_0_Current",
#                  "time_spent_watching_television_tv_f1070_0_0")
# disease_v1 <- c("age_e11_first_reported_non_insulin_dependent_diabetes_mellitus_f130708_0_0",
#                 "age_j43_first_reported_emphysema_f131490_0_0",
#                 "age_j44_first_reported_other_chronic_obstructive_pulmonary_disease_f131492_0_0",
#                 "age_n18_first_reported_chronic_renal_failure_f132032_0_0",
#                 "age_i10_first_reported_essential_primary_hypertension_f131286_0_0")
# 
# combo_df <- expand.grid(
#   exposure = exposure_v1,
#   disease  = disease_v1,
#   stringsAsFactors = FALSE
# )
# 
# combo_df
# 
# fwrite(combo_df,"/n/groups/patel/shakson_ukb/UK_Biobank/BScripts/Module6/Production/module6_config.csv")
# 
# 




