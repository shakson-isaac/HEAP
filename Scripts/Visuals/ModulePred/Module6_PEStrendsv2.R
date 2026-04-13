#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
})

# ============================================================
# Config (same as before)
# ============================================================
cfg <- list(
  out_dir = "/n/groups/patel/shakson_ukb/UK_Biobank/Data/Parallel/PES_test",
  paths = list(
    pxs_rds = "/n/scratch/users/s/shi872/UKB_intermediate/UKB_PGS_PXS_load.rds",
    t2e_rds_prefix = "/n/groups/patel/IGLOO/UKB/Mediation/Data/UKB_MDstore_"
  ),
  CovarSpec = list(
    Type1 = c("age_when_attended_assessment_centre_f21003_0_0", "sex_f31_0_0"),
    Type2 = c("age_when_attended_assessment_centre_f21003_0_0", "sex_f31_0_0",
              "body_mass_index_bmi_f23104_0_0", "fasting_time_f74_0_0"),
    Type3 = c("age_when_attended_assessment_centre_f21003_0_0", "sex_f31_0_0",
              "age2","age_sex","age2_sex",
              "body_mass_index_bmi_f23104_0_0", "fasting_time_f74_0_0",
              "uk_biobank_assessment_centre_f54_0_0",
              paste0("genetic_principal_components_f22009_0_",1:20)),
    Type4 = c("age_when_attended_assessment_centre_f21003_0_0", "sex_f31_0_0",
              "age2","age_sex","age2_sex",
              "body_mass_index_bmi_f23104_0_0", "fasting_time_f74_0_0",
              "uk_biobank_assessment_centre_f54_0_0",
              paste0("genetic_principal_components_f22009_0_",1:20),
              "combined_Blood_pressure_medication",
              "combined_Hormone_replacement_therapy",
              "combined_Oral_contraceptive_pill_or_minipill",
              "combined_Insulin",
              "combined_Cholesterol_lowering_medication",
              "combined_Do_not_know",
              "combined_None_of_the_above",
              "combined_Prefer_not_to_answer"),
    Type5 = NULL
  )
)

as_pxs <- function(x) {
  if (is.list(x) && !isS4(x)) return(x)
  if (!isS4(x)) stop("PXS object must be S4 or list")
  list(
    covars_df   = x@covars_df,
    covars_list = x@covars_list
  )
}

# ============================================================
# NEW: survival-time recode that *keeps* prevalent cases
# - does NOT filter out baseline>=event_age
# - still computes DZ_status and DZ_survtime as before
# ============================================================
survival_time_keep_prevalent <- function(Time2Event_df, event_age_col,
                                         recode_status = "DZ_status",
                                         recode_survtime = "DZ_survtime") {
  
  T2E_df <- as.data.frame(Time2Event_df)
  
  if (!"recode_age_of_assessment_0_0" %in% names(T2E_df)) {
    stop("Time2Event_df missing recode_age_of_assessment_0_0")
  }
  T2E_df <- T2E_df[!is.na(T2E_df$recode_age_of_assessment_0_0), ]
  
  censor_age <- pmin(
    T2E_df$recode_age_of_death_0_0,
    T2E_df$age_of_removal_0_0,
    T2E_df$age_of_lastfollowup,
    na.rm = TRUE
  )
  
  # event occurred at some point before censoring (can be prevalent or incident)
  T2E_df[[recode_status]] <- as.integer(
    !is.na(T2E_df[[event_age_col]]) & (T2E_df[[event_age_col]] <= censor_age)
  )
  
  # for status==1, DZ_survtime is the event age (even if < baseline age)
  T2E_df[[recode_survtime]] <- ifelse(
    T2E_df[[recode_status]] == 1,
    T2E_df[[event_age_col]],
    censor_age
  )
  
  T2E_df
}

# ============================================================
# Args
# ============================================================
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop(paste0(
    "Usage:\n",
    "  Rscript PES_plot_time_to_diagnosis_prevalent_plus_incident.R <covarType> <exposure_id> <disease_age_col> [pes_col]\n\n",
    "Example:\n",
    "  Rscript PES_plot_time_to_diagnosis_prevalent_plus_incident.R Type5 smoking_status_f20116_0_0_Current age_j43_first_reported_emphysema_f131490_0_0 pes_prot_z\n"
  ))
}

covarType       <- as.character(args[1])
exposure_id     <- as.character(args[2])
disease_age_col <- as.character(args[3])
pes_col         <- if (length(args) >= 4) as.character(args[4]) else "pes_prot_z"


covarType <- "Type5"
exposure_id <-  "smoking_status_f20116_0_0_Current" #"usual_walking_pace_f924_0_0"
disease_age_col <- "age_j44_first_reported_other_chronic_obstructive_pulmonary_disease_f131492_0_0" #"age_j43_first_reported_emphysema_f131490_0_0" # "age_e11_first_reported_non_insulin_dependent_diabetes_mellitus_f130708_0_0"


# plot window (edit as you like)
XMIN <- -20
XMAX <-  10

# ============================================================
# Load existing OOF/PES (no refit)
# ============================================================
covar_out_dir <- file.path(cfg$out_dir, covarType)
out_prefix    <- file.path(covar_out_dir, paste0("PES_", covarType, "_", exposure_id))
oof_rds       <- paste0(out_prefix, "_OOF.rds")
if (!file.exists(oof_rds)) stop("Missing OOF RDS: ", oof_rds)

oof_tbl <- readRDS(oof_rds)
if (!pes_col %in% names(oof_tbl)) stop("pes_col not in oof_tbl: ", pes_col)

# covars + t2e
pxs0 <- as_pxs(readRDS(cfg$paths$pxs_rds))
covars_used <- if (covarType == "Type5") pxs0$covars_list else cfg$CovarSpec[[covarType]]
covars_df <- as.data.frame(pxs0$covars_df)

MDloader_path <- paste0(cfg$paths$t2e_rds_prefix, covarType, ".rds")
MDloader <- readRDS(MDloader_path)
t2e_df <- as.data.frame(MDloader@DZ_df)

needed_t2e <- c("eid","recode_age_of_assessment_0_0","recode_age_of_death_0_0","age_of_removal_0_0","age_of_lastfollowup")
stopifnot(all(needed_t2e %in% names(t2e_df)))
if (!disease_age_col %in% names(t2e_df)) stop("disease_age_col not found: ", disease_age_col)

# join
base <- oof_tbl %>%
  dplyr::select(eid, exposure_id, exposure_type, y_raw, all_of(pes_col)) %>%
  dplyr::rename(pes_z = all_of(pes_col)) %>%
  dplyr::inner_join(covars_df[, c("eid", covars_used), drop = FALSE], by = "eid") %>%
  dplyr::inner_join(t2e_df[, c(needed_t2e, disease_age_col), drop = FALSE], by = "eid")

# keep prevalent + incident
icd <- survival_time_keep_prevalent(base, event_age_col = disease_age_col,
                                    recode_status = "DZ_status", recode_survtime = "DZ_survtime")

# time-to-dx in years (event age minus baseline age)
icd$event_age  <- icd[[disease_age_col]]
icd$time_to_dx <- icd$event_age - icd$recode_age_of_assessment_0_0

# diagnosis-aligned x:
#   incident => time_to_dx > 0 => x < 0 (before dx)
#   prevalent => time_to_dx < 0 => x > 0 (after dx)
icd$x_years_rel_dx0 <- -icd$time_to_dx

# label case type among cases
icd$case_type <- dplyr::case_when(
  icd$DZ_status != 1 ~ "control/censored",
  is.finite(icd$time_to_dx) & icd$time_to_dx > 0 ~ "incident",
  is.finite(icd$time_to_dx) & icd$time_to_dx < 0 ~ "prevalent",
  TRUE ~ "case_unknown"
)

# plot data: cases only (both types)
plot_df <- icd %>%
  filter(DZ_status == 1) %>%
  filter(is.finite(pes_z), is.finite(x_years_rel_dx0)) %>%
  filter(x_years_rel_dx0 >= XMIN, x_years_rel_dx0 <= XMAX)

message("Cases plotted: n=", nrow(plot_df),
        " | incident=", sum(plot_df$case_type == "incident"),
        " | prevalent=", sum(plot_df$case_type == "prevalent"))

if (nrow(plot_df) < 20) {
  stop("Too few cases after windowing/filtering. Increase XMIN/XMAX or pick another disease.")
}

# plot
p <- ggplot(plot_df, aes(x = x_years_rel_dx0, y = pes_z)) +
  geom_point(aes(shape = case_type), alpha = 0.18, size = 0.7) +
  geom_smooth(method = "loess", se = TRUE, span = 0.9) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(
    title = paste0("PES vs time-to-diagnosis: ", disease_age_col),
    subtitle = paste0("Exposure: ", exposure_id,
                      " | CovarType: ", covarType,
                      " | PES: ", pes_col,
                      " | Incident + prevalent cases"),
    x = "Years relative to diagnosis (0 = diagnosis; negative = before; positive = after)",
    y = paste0(pes_col, " (z)"),
    shape = "Case type"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )
p

out_png <- file.path(covar_out_dir,
                     paste0("PES_vs_TimeToDx__INCplusPREV__", covarType, "__",
                            exposure_id, "__", disease_age_col, "__", pes_col, ".png"))
out_pdf <- sub("\\.png$", ".pdf", out_png)

#ggsave(out_png, p, width = 10, height = 5, dpi = 300)
#ggsave(out_pdf, p, width = 10, height = 5)

# save data used for plotting
out_tsv <- sub("\\.png$", ".plotdata.tsv", out_png)
#fwrite(as.data.table(plot_df[, c("eid","pes_z","case_type","x_years_rel_dx0","time_to_dx")]),
#       out_tsv, sep = "\t")

message("Saved:\n  ", out_png, "\n  ", out_pdf, "\n  ", out_tsv)
message("DONE.")