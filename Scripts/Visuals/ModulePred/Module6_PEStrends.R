#!/usr/bin/env Rscript

# ============================================================
# Plot PES vs time_to_diagnosis (aligned to dx = 0)
# Re-uses existing OOF/PES outputs (does NOT rerun PES models)
#
# Usage:
#   Rscript PES_plot_time_to_diagnosis.R <covarType> <exposure_id> <disease_age_col> [pes_col]
#
# Example:
#   Rscript PES_plot_time_to_diagnosis.R Type5 summed_met_minutes_per_week_for_all_activity_f22040_0_0 age_E11 pes_prot_z
#
# Notes:
# - disease_age_col must be a column in MDstore@DZ_df (e.g. "age_E11", etc.)
# - For incident cases, x-axis is negative years before diagnosis:
#     x = -(age_at_event - age_at_assessment)
#   so diagnosis occurs at x = 0.
# - By default plots incident cases only (recommended).
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(purrr)
})

# -----------------------------
# Config (match your pipeline)
# -----------------------------
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

# -----------------------------
# Helpers (copied logic)
# -----------------------------
as_pxs <- function(x) {
  if (is.list(x) && !isS4(x)) return(x)
  if (!isS4(x)) stop("PXS object must be S4 or list")
  list(
    covars_df   = x@covars_df,
    covars_list = x@covars_list
  )
}

survival_time <- function(Time2Event_df, event_age_col,
                          recode_status = "DZ_status",
                          recode_survtime = "DZ_survtime") {
  
  T2E_df <- as.data.frame(Time2Event_df)
  
  if (!"recode_age_of_assessment_0_0" %in% names(T2E_df)) {
    stop("Time2Event_df missing recode_age_of_assessment_0_0")
  }
  
  T2E_df <- T2E_df[!is.na(T2E_df$recode_age_of_assessment_0_0), ]
  
  # exclude prevalent cases at baseline
  T2E_df <- T2E_df[
    is.na(T2E_df[[event_age_col]]) |
      (T2E_df$recode_age_of_assessment_0_0 < T2E_df[[event_age_col]]),
  ]
  
  censor_age <- pmin(
    T2E_df$recode_age_of_death_0_0,
    T2E_df$age_of_removal_0_0,
    T2E_df$age_of_lastfollowup,
    na.rm = TRUE
  )
  
  T2E_df[[recode_status]] <- as.integer(
    !is.na(T2E_df[[event_age_col]]) & (T2E_df[[event_age_col]] <= censor_age)
  )
  
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
    "  Rscript PES_plot_time_to_diagnosis.R <covarType> <exposure_id> <disease_age_col> [pes_col]\n\n",
    "Example:\n",
    "  Rscript PES_plot_time_to_diagnosis.R Type5 summed_met_minutes_per_week_for_all_activity_f22040_0_0 age_E11 pes_prot_z\n"
  ))
}

covarType      <- as.character(args[1])
exposure_id    <- as.character(args[2])
disease_age_col<- as.character(args[3])
pes_col        <- if (length(args) >= 4) as.character(args[4]) else "pes_prot_z"

covarType <- "Type5"
exposure_id <- "smoking_status_f20116_0_0_Current"
disease_age_col <- "age_j43_first_reported_emphysema_f131490_0_0"


if (!covarType %in% names(cfg$CovarSpec) && covarType != "Type5") stop("Unknown covarType: ", covarType)

# ============================================================
# Locate existing OOF/PES
# ============================================================
covar_out_dir <- file.path(cfg$out_dir, covarType)
out_prefix    <- file.path(covar_out_dir, paste0("PES_", covarType, "_", exposure_id))
oof_rds       <- paste0(out_prefix, "_OOF.rds")

if (!file.exists(oof_rds)) {
  stop("Missing OOF/PES RDS (run your main pipeline first):\n  ", oof_rds)
}

message("Reading OOF/PES from:\n  ", oof_rds)
oof_tbl <- readRDS(oof_rds)
if (!pes_col %in% names(oof_tbl)) stop("pes_col not found in oof_tbl: ", pes_col)

# ============================================================
# Load covariates + disease time-to-event
# ============================================================
pxs0 <- as_pxs(readRDS(cfg$paths$pxs_rds))
covars_used <- if (covarType == "Type5") pxs0$covars_list else cfg$CovarSpec[[covarType]]
if (is.null(covars_used) || length(covars_used) == 0) stop("No covariates found for covarType: ", covarType)

covars_df <- as.data.frame(pxs0$covars_df)
if (!"eid" %in% names(covars_df)) stop("covars_df must include eid")

MDloader_path <- paste0(cfg$paths$t2e_rds_prefix, covarType, ".rds")
if (!file.exists(MDloader_path)) stop("Missing MDstore file:\n  ", MDloader_path)

MDloader <- readRDS(MDloader_path)
t2e_df <- as.data.frame(MDloader@DZ_df)
if (!"eid" %in% names(t2e_df)) stop("t2e_df must include eid")

needed_t2e <- c("eid","recode_age_of_assessment_0_0","recode_age_of_death_0_0","age_of_removal_0_0","age_of_lastfollowup")
miss_t2e <- setdiff(needed_t2e, names(t2e_df))
if (length(miss_t2e) > 0) stop("t2e_df missing columns: ", paste(miss_t2e, collapse=", "))

if (!disease_age_col %in% names(t2e_df)) {
  stop("disease_age_col not found in t2e_df: ", disease_age_col,
       "\nTip: list disease cols with: grep('^age_', names(t2e_df), value=TRUE)")
}

# ============================================================
# Build the same joined dataframe as Cox stage
# ============================================================
base <- oof_tbl %>%
  dplyr::select(eid, exposure_id, exposure_type, y_raw, all_of(pes_col)) %>%
  dplyr::rename(pes_z = all_of(pes_col)) %>%
  dplyr::inner_join(covars_df[, c("eid", covars_used), drop = FALSE], by = "eid") %>%
  dplyr::inner_join(t2e_df[, c(needed_t2e, disease_age_col), drop = FALSE], by = "eid")

icd <- survival_time(base, event_age_col = disease_age_col,
                     recode_status = "DZ_status", recode_survtime = "DZ_survtime")

# time from baseline assessment to event/censor (years)
icd$time_from_baseline <- icd$DZ_survtime - icd$recode_age_of_assessment_0_0

# For cases, event time uses the disease age column directly (more explicit)
icd$event_age <- icd[[disease_age_col]]
icd$time_to_dx <- icd$event_age - icd$recode_age_of_assessment_0_0

# diagnosis-aligned x: negative years before dx, dx at 0
icd$x_years_to_dx0 <- -icd$time_to_dx

# Basic filtering
icd <- icd %>%
  filter(is.finite(pes_z)) %>%
  filter(is.finite(recode_age_of_assessment_0_0))

# Recommended: cases-only for diagnosis-aligned plot
plot_df <- icd %>%
  filter(DZ_status == 1) %>%
  filter(is.finite(x_years_to_dx0)) %>%
  # keep plausible window (edit as you like)
  filter(x_years_to_dx0 <= 0, x_years_to_dx0 >= -20)

if (nrow(plot_df) < 50) {
  stop("Too few incident cases after filtering (n=", nrow(plot_df), "). Try widening the window or pick another disease.")
}

message("Plotting n=", nrow(plot_df), " incident cases for ", disease_age_col,
        " | window: [-20, 0] years relative to diagnosis.")

# ============================================================
# Plot (shaded phases + smooth)
# ============================================================
phase_risk  <- c(-20, -5)
phase_trans <- c(-5, 0)

p <- ggplot(plot_df, aes(x = x_years_to_dx0, y = pes_z)) +
  annotate("rect", xmin = phase_risk[1], xmax = phase_risk[2], ymin = -Inf, ymax = Inf, alpha = 0.08) +
  annotate("rect", xmin = phase_trans[1], xmax = phase_trans[2], ymin = -Inf, ymax = Inf, alpha = 0.12) +
  geom_point(alpha = 0.15, size = 0.7) +
  geom_smooth(method = "loess", se = TRUE, span = 0.9) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(
    title = paste0("PES vs time-to-diagnosis: ", disease_age_col),
    subtitle = paste0("Exposure: ", exposure_id, " | CovarType: ", covarType,
                      " | PES: ", pes_col, " | Incident cases only"),
    x = "Years relative to diagnosis (0 = diagnosis; negative = years before)",
    y = paste0(pes_col, " (z)")
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

# ============================================================
# Save
# ============================================================
out_png <- file.path(covar_out_dir,
                     paste0("PES_vs_TimeToDx__", covarType, "__", exposure_id, "__", disease_age_col, "__", pes_col, ".png"))
out_pdf <- sub("\\.png$", ".pdf", out_png)

ggsave(out_png, p, width = 10, height = 5, dpi = 300)
ggsave(out_pdf, p, width = 10, height = 5)

message("Saved:\n  ", out_png, "\n  ", out_pdf)

# Also write the plotted data (so you can reuse it in other figures)
out_tsv <- sub("\\.png$", ".plotdata.tsv", out_png)
fwrite(as.data.table(plot_df[, c("eid","pes_z","x_years_to_dx0","time_to_dx")]),
       out_tsv, sep = "\t")
message("Saved plot data:\n  ", out_tsv)

message("DONE.")


