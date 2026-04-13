#!/usr/bin/env Rscript

# ============================================================
# T2D diagnosis-aligned trajectories for ONE protein:
#   Panel A: Observed protein (z)
#   Panel B: Genetic component (Gcis + Gtrans, z) from Module1 OOF_components
#   Panel C: Non-genetic component (E + GxE, z) from Module1 OOF_components
#
# Includes incident + prevalent cases, aligned to diagnosis:
#   x = -(event_age - baseline_age)
#     incident -> x < 0  (years before diagnosis)
#     prevalent -> x > 0 (years after diagnosis; baseline after dx)
#
# Also outputs spline-trend p-values (cases only).
#
# Usage:
#   Rscript T2D_oneProtein_trajectory_decomp.R <covarType> <protID> <disease_age_col> [xmin] [xmax]
#
# Example:
#   Rscript T2D_oneProtein_trajectory_decomp.R Type5 LEP age_e11_first_reported_non_insulin_dependent_diabetes_f131294_0_0 -15 10
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(splines)
  library(gridExtra)
})

# -----------------------------
# Config (EDIT if needed)
# -----------------------------
cfg <- list(
  # PXS object with UKBprot_df (observed proteins)
  pxs_rds = "/n/scratch/users/s/shi872/UKB_intermediate/UKB_PGS_PXS_load.rds",
  
  # MDstore for disease timing + covariates (same as Module3)
  mdstore_prefix = "/n/groups/patel/IGLOO/UKB/Mediation/Data/UKB_MDstore_",
  
  # Module1 output root
  module1_root = "/n/groups/patel/shakson_ukb/UK_Biobank/Data/Parallel/Module1",
  
  # where to save figures/tables
  out_root = "/n/groups/patel/shakson_ukb/UK_Biobank/Output/Biomarker/Trajectories",
  
  # minimal covariates for adjustment if you want them in trend tests
  # (if Type5, we’ll use MDstore@covars_list by default)
  fallback_covars = c("age_when_attended_assessment_centre_f21003_0_0", "sex_f31_0_0")
)

prot_clean <- function(protID) gsub("-", "_", protID)

# -----------------------------
# survival recode (KEEP prevalent)
# -----------------------------
survival_time_keep_prevalent <- function(df, event_age_col,
                                         recode_status = "DZ_status",
                                         recode_survtime = "DZ_survtime") {
  
  df <- as.data.frame(df)
  df <- df[!is.na(df$recode_age_of_assessment_0_0), ]
  
  censor_age <- pmin(
    df$recode_age_of_death_0_0,
    df$age_of_removal_0_0,
    df$age_of_lastfollowup,
    na.rm = TRUE
  )
  
  df[[recode_status]] <- as.integer(
    !is.na(df[[event_age_col]]) & (df[[event_age_col]] <= censor_age)
  )
  
  df[[recode_survtime]] <- ifelse(
    df[[recode_status]] == 1,
    df[[event_age_col]],
    censor_age
  )
  
  df
}

zscore_vec <- function(x) {
  s <- sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(0, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}

# ------------------------------------------------------------
# Args
# ------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop(paste0(
    "Usage:\n",
    "  Rscript T2D_oneProtein_trajectory_decomp.R <covarType> <protID> <disease_age_col> [xmin] [xmax]\n\n",
    "Example:\n",
    "  Rscript T2D_oneProtein_trajectory_decomp.R Type5 LEP age_e11_first_reported_non_insulin_dependent_diabetes_f131294_0_0 -15 10\n"
  ))
}

covarType       <- as.character(args[1])
protID_raw      <- as.character(args[2])
disease_age_col <- as.character(args[3])
xmin <- if (length(args) >= 4) as.numeric(args[4]) else -20
xmax <- if (length(args) >= 5) as.numeric(args[5]) else  10

covarType = "Type5"
protID_raw =  "GFAP" #"LEP" 
disease_age_col = "age_g30_first_reported_alzheimers_disease_f131036_0_0" #"age_e66_first_reported_obesity_f130792_0_0"

#"CXCL17" 
#"ASGR1" 
#"LEP"

#"age_j43_first_reported_emphysema_f131490_0_0" 
#"age_e78_first_reported_disorders_of_lipoprotein_metabolism_and_other_lipidaemias_f130814_0_0" 
#"age_j44_first_reported_other_chronic_obstructive_pulmonary_disease_f131492_0_0" 
#"age_j43_first_reported_emphysema_f131490_0_0" 
#"age_e11_first_reported_non_insulin_dependent_diabetes_mellitus_f130708_0_0"
#"age_e66_first_reported_obesity_f130792_0_0"
#"age_e78_first_reported_disorders_of_lipoprotein_metabolism_and_other_lipidaemias_f130814_0_0"
#"age_g30_first_reported_alzheimers_disease_f131036_0_0"
#"age_n18_first_reported_chronic_renal_failure_f132032_0_0"
#"age_f10_first_reported_mental_and_behavioural_disorders_due_to_use_of_alcohol_f130854_0_0"
#"age_f00_first_reported_dementia_in_alzheimers_disease_f130836_0_0"


protID <- prot_clean(protID_raw)

# ------------------------------------------------------------
# Load Module1 OOF component contributions
# ------------------------------------------------------------
in_dir <- file.path(cfg$module1_root, covarType)

oof_comp_fp <- file.path(in_dir, paste0("OOF_components_", protID_raw, "_", covarType, ".rds"))
if (!file.exists(oof_comp_fp)) {
  # try cleaned name fallback
  oof_comp_fp2 <- file.path(in_dir, paste0("OOF_components_", protID, "_", covarType, ".rds"))
  if (file.exists(oof_comp_fp2)) oof_comp_fp <- oof_comp_fp2
}

if (!file.exists(oof_comp_fp)) {
  stop("Missing OOF_components file. Tried:\n  ", oof_comp_fp, "\n",
       "Tip: check the exact saved filename in: ", in_dir)
}

message("Reading: ", oof_comp_fp)
oof_comp <- readRDS(oof_comp_fp) %>% as.data.frame()
stopifnot(all(c("eid","PredTotal") %in% names(oof_comp)))

# Ensure component columns exist (some may be missing for some proteins)
if (!"Gcis" %in% names(oof_comp)) oof_comp$Gcis <- 0
if (!"Gtrans" %in% names(oof_comp)) oof_comp$Gtrans <- 0

# Non-genetic component columns
E_cols    <- grep("^E_", names(oof_comp), value = TRUE)
GxEcis_c  <- grep("^GxEcis_", names(oof_comp), value = TRUE)
GxEtr_c   <- grep("^GxEtrans_", names(oof_comp), value = TRUE)

# Summaries
oof_comp <- oof_comp %>%
  mutate(
    comp_G = Gcis + Gtrans,
    comp_E = if (length(E_cols) > 0) rowSums(across(all_of(E_cols)), na.rm = TRUE) else 0,
    comp_GxE = 0
  )

if (length(GxEcis_c) > 0) oof_comp$comp_GxE <- oof_comp$comp_GxE + rowSums(oof_comp[, GxEcis_c, drop = FALSE], na.rm = TRUE)
if (length(GxEtr_c)  > 0) oof_comp$comp_GxE <- oof_comp$comp_GxE + rowSums(oof_comp[, GxEtr_c,  drop = FALSE], na.rm = TRUE)

oof_comp <- oof_comp %>%
  mutate(comp_nonG = comp_E + comp_GxE)

# ------------------------------------------------------------
# Load observed protein from pxs_rds
# ------------------------------------------------------------
pxs0 <- readRDS(cfg$pxs_rds)

# Handle either S4 or list
UKBprot_df <- NULL
if (isS4(pxs0)) UKBprot_df <- as.data.frame(pxs0@UKBprot_df)
if (is.list(pxs0) && !is.null(pxs0$UKBprot_df)) UKBprot_df <- as.data.frame(pxs0$UKBprot_df)

if (is.null(UKBprot_df)) stop("Could not locate UKBprot_df in pxs object.")

if (!"eid" %in% names(UKBprot_df)) stop("UKBprot_df missing eid.")
if (!protID %in% names(UKBprot_df)) {
  stop("Protein column not found in UKBprot_df: ", protID,
       "\nTip: check if the protein is saved with '-' vs '_' and pass protID accordingly.")
}

obsP <- UKBprot_df %>% select(eid, !!protID) %>% rename(obs_protein = !!protID)

# ------------------------------------------------------------
# Load MDstore (disease timing + covariates)
# ------------------------------------------------------------
md_fp <- paste0(cfg$mdstore_prefix, covarType, ".rds")
if (!file.exists(md_fp)) stop("Missing MDstore: ", md_fp)

MD <- readRDS(md_fp)

DZ_df <- as.data.frame(MD@DZ_df)
cov_df <- as.data.frame(MD@covars_df)
covars_used <- MD@covars_list
if (is.null(covars_used) || length(covars_used) == 0) covars_used <- cfg$fallback_covars

needed_t2e <- c("eid","recode_age_of_assessment_0_0","recode_age_of_death_0_0","age_of_removal_0_0","age_of_lastfollowup")
miss <- setdiff(needed_t2e, names(DZ_df))
if (length(miss) > 0) stop("DZ_df missing columns: ", paste(miss, collapse = ", "))
if (!disease_age_col %in% names(DZ_df)) stop("disease_age_col not found in DZ_df: ", disease_age_col)

# ------------------------------------------------------------
# Build analysis dataframe
# ------------------------------------------------------------
df0 <- oof_comp %>%
  inner_join(obsP, by = "eid") %>%
  inner_join(cov_df[, c("eid", intersect(covars_used, names(cov_df))), drop = FALSE], by = "eid") %>%
  inner_join(DZ_df[, c(needed_t2e, disease_age_col), drop = FALSE], by = "eid")

df1 <- survival_time_keep_prevalent(df0, event_age_col = disease_age_col)

df1$event_age  <- df1[[disease_age_col]]
df1$time_to_dx <- df1$event_age - df1$recode_age_of_assessment_0_0
df1$x_years_rel_dx0 <- -df1$time_to_dx

df1$case_type <- dplyr::case_when(
  df1$DZ_status != 1 ~ "control/censored",
  is.finite(df1$time_to_dx) & df1$time_to_dx > 0 ~ "incident",
  is.finite(df1$time_to_dx) & df1$time_to_dx < 0 ~ "prevalent",
  TRUE ~ "case_unknown"
)

# Use cases only for diagnosis-aligned plot
plot_df <- df1 %>%
  filter(DZ_status == 1) %>%
  filter(is.finite(x_years_rel_dx0)) %>%
  filter(x_years_rel_dx0 >= xmin, x_years_rel_dx0 <= xmax) %>%
  mutate(
    obsP_z      = zscore_vec(obs_protein),
    compG_z     = zscore_vec(comp_G),
    compNonG_z  = zscore_vec(comp_nonG)
  )

message("Cases plotted: n=", nrow(plot_df),
        " | incident=", sum(plot_df$case_type == "incident", na.rm = TRUE),
        " | prevalent=", sum(plot_df$case_type == "prevalent", na.rm = TRUE))

if (nrow(plot_df) < 200) {
  warning("Low number of cases after filtering. Consider widening xmin/xmax.")
}

# ------------------------------------------------------------
# Trend significance tests (cases only)
#   Y ~ ns(x, df=3) + covars
#   p = LRT comparing spline vs no-spline
# ------------------------------------------------------------
do_spline_test <- function(df, ycol, xcol = "x_years_rel_dx0", covars = character()) {
  df <- df %>% filter(is.finite(.data[[ycol]]), is.finite(.data[[xcol]]))
  covars <- intersect(covars, names(df))
  base_rhs <- if (length(covars) > 0) paste(covars, collapse = " + ") else "1"
  
  f0 <- as.formula(paste0(ycol, " ~ ", base_rhs))
  f1 <- as.formula(paste0(ycol, " ~ ns(", xcol, ", df=3) + ", base_rhs))
  
  m0 <- lm(f0, data = df)
  m1 <- lm(f1, data = df)
  a <- anova(m0, m1)
  p <- a$`Pr(>F)`[2]
  tibble::tibble(y = ycol, p_spline = as.numeric(p), n = nrow(df))
}

trend_tbl <- dplyr::bind_rows(
  do_spline_test(plot_df, "obsP_z",     covars = covars_used),
  do_spline_test(plot_df, "compG_z",    covars = covars_used),
  do_spline_test(plot_df, "compNonG_z", covars = covars_used)
)

# ------------------------------------------------------------
# Plot panels
# ------------------------------------------------------------
base_theme <- theme_bw() +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

make_panel <- function(df, ycol, ylab, title) {
  ggplot(df, aes(x = x_years_rel_dx0, y = .data[[ycol]])) +
    geom_point(aes(shape = case_type), alpha = 0.18, size = 0.7) +
    geom_smooth(method = "loess", se = TRUE, span = 0.9) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    labs(title = title,
         x = "Years relative to diagnosis (0 = diagnosis; negative = before; positive = after)",
         y = ylab,
         shape = "Case type") +
    base_theme
}

pA <- make_panel(plot_df, "obsP_z",
                 ylab = paste0(protID_raw, " (observed, z)"),
                 title = "A  Observed protein")

pB <- make_panel(plot_df, "compG_z",
                 ylab = "Genetic component (Gcis + Gtrans, z)",
                 title = "B  Genetic component")

pC <- make_panel(plot_df, "compNonG_z",
                 ylab = "Non-genetic component (E + GxE, z)",
                 title = "C  Non-genetic component")

# Arrange (3 rows)
g <- gridExtra::arrangeGrob(pA, pB, pC, ncol = 1)


# ------------------------------------------------------------
# Save outputs
# ------------------------------------------------------------
out_dir <- file.path(cfg$out_root, covarType, protID_raw)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

stub <- paste0("T2Dtraj_", covarType, "__", protID_raw, "__", disease_age_col,
               "__", xmin, "_to_", xmax)

fig_png <- file.path(out_dir, paste0(stub, ".png"))
fig_pdf <- file.path(out_dir, paste0(stub, ".pdf"))
tbl_tsv <- file.path(out_dir, paste0(stub, "_trendTests.tsv"))
dat_tsv <- file.path(out_dir, paste0(stub, "_plotdata.tsv"))

ggsave(fig_png, g, width = 10, height = 12, dpi = 300)
ggsave(fig_pdf, g, width = 10, height = 12)

data.table::fwrite(as.data.table(trend_tbl), tbl_tsv, sep = "\t")
data.table::fwrite(as.data.table(plot_df), dat_tsv, sep = "\t")

message("Saved figure:\n  ", fig_png, "\n  ", fig_pdf)
message("Saved trend tests:\n  ", tbl_tsv)
message("Saved plot data:\n  ", dat_tsv)
message("DONE.")