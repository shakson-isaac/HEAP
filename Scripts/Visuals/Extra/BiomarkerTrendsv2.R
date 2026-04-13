#!/usr/bin/env Rscript

# ============================================================
# Nested case-control (risk-set) trajectories for T2D
# One protein example with decomposition:
#   A) Observed protein (z)
#   B) Genetic component: Gcis + Gtrans (z)   [from Module1 OOF_components]
#   C) Non-genetic component: sum(E_*) + sum(GxE*) (z) [from Module1 OOF_components]
#
# Matches each INCIDENT case to k controls on:
#   - baseline age (recode_age_of_assessment_0_0) within age_caliper_years
#   - sex (sex_f31_0_0) exact
#   - BMI (body_mass_index_bmi_f23104_0_0) within bmi_caliper_units
#
# Risk-set control eligibility at case index age:
#   - not diagnosed before index (event age is NA or > index_age)
#   - still under observation at index (censor_age >= index_age)
#
# Plots binned mean ± SE vs years-to-diagnosis (negative), incident cases only.
#
# Usage:
#   Rscript T2D_nested_case_control_trajectory_oneProtein.R <covarType> <protID> <t2d_age_col> [k_controls] [xmin] [xmax]
#
# Example:
#   Rscript T2D_nested_case_control_trajectory_oneProtein.R Type5 LEP age_e11_first_reported_non_insulin_dependent_diabetes_f131294_0_0 5 -15 -2
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(gridExtra)
})

# -----------------------------
# Config (EDIT if needed)
# -----------------------------
cfg <- list(
  # PXS object with UKBprot_df (observed proteins)
  pxs_rds = "/n/scratch/users/s/shi872/UKB_intermediate/UKB_PGS_PXS_load.rds",
  
  # MDstore for disease timing + covariates
  mdstore_prefix = "/n/groups/patel/IGLOO/UKB/Mediation/Data/UKB_MDstore_",
  
  # Module1 output root
  module1_root = "/n/groups/patel/shakson_ukb/UK_Biobank/Data/Parallel/Module1",
  
  # output root for this analysis
  out_root = "/n/groups/patel/shakson_ukb/UK_Biobank/Output/Biomarker/Trajectories",
  
  # matching calipers
  age_caliper_years = 2,
  bmi_caliper_units = 2,   # kg/m^2
  # if insufficient controls, progressively relax BMI then age
  bmi_relax_seq = c(2, 3, 5, 8),
  age_relax_seq = c(2, 3, 5),
  
  # binning (yearly bins)
  bin_width_years = 1
)

prot_clean <- function(protID) gsub("-", "_", protID)

zscore_vec <- function(x) {
  s <- sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(0, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}

# -----------------------------
# Incident-only survival helper (EXCLUDES prevalent)
# -----------------------------
survival_time_incident_only <- function(df, event_age_col,
                                        recode_status = "DZ_status",
                                        recode_survtime = "DZ_survtime") {
  df <- as.data.frame(df)
  
  stopifnot("recode_age_of_assessment_0_0" %in% names(df))
  df <- df[!is.na(df$recode_age_of_assessment_0_0), ]
  
  # exclude prevalent cases at baseline
  df <- df[
    is.na(df[[event_age_col]]) | (df$recode_age_of_assessment_0_0 < df[[event_age_col]]),
  ]
  
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
  
  df$censor_age <- censor_age
  df
}

# -----------------------------
# Risk-set matching
# -----------------------------
match_controls_for_case <- function(case_row, candidates, k,
                                    age_relax_seq, bmi_relax_seq) {
  # case_row: 1-row data.frame
  # candidates: eligible controls already filtered for risk-set at this index age
  
  for (ac in age_relax_seq) {
    for (bc in bmi_relax_seq) {
      pool <- candidates %>%
        filter(abs(recode_age_of_assessment_0_0 - case_row$recode_age_of_assessment_0_0) <= ac) %>%
        filter(abs(body_mass_index_bmi_f23104_0_0 - case_row$body_mass_index_bmi_f23104_0_0) <= bc)
      
      if (nrow(pool) >= k) {
        sel <- pool %>% slice_sample(n = k, replace = FALSE)
        sel$age_caliper_used <- ac
        sel$bmi_caliper_used <- bc
        return(sel)
      }
    }
  }
  
  # if still insufficient, take whatever exists at max relax (or none)
  pool <- candidates %>%
    filter(abs(recode_age_of_assessment_0_0 - case_row$recode_age_of_assessment_0_0) <= max(age_relax_seq)) %>%
    filter(abs(body_mass_index_bmi_f23104_0_0 - case_row$body_mass_index_bmi_f23104_0_0) <= max(bmi_relax_seq))
  
  if (nrow(pool) == 0) return(NULL)
  sel <- pool %>% slice_sample(n = min(k, nrow(pool)), replace = FALSE)
  sel$age_caliper_used <- max(age_relax_seq)
  sel$bmi_caliper_used <- max(bmi_relax_seq)
  sel$k_shortfall <- k - nrow(sel)
  sel
}

# -----------------------------
# CLI args
# -----------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop(paste0(
    "Usage:\n",
    "  Rscript T2D_nested_case_control_trajectory_oneProtein.R <covarType> <protID> <t2d_age_col> [k_controls] [xmin] [xmax]\n\n",
    "Example:\n",
    "  Rscript T2D_nested_case_control_trajectory_oneProtein.R Type5 LEP age_e11_first_reported_non_insulin_dependent_diabetes_f130708_0_0 5 -15 -2\n"
  ))
}

covarType <- as.character(args[1])
protID_raw <- as.character(args[2])
t2d_age_col <- as.character(args[3])

covarType = "Type5"
#protID_raw = "LEP"
protID_raw = "LEP"
t2d_age_col = "age_e66_first_reported_obesity_f130792_0_0"
#"age_e78_first_reported_disorders_of_lipoprotein_metabolism_and_other_lipidaemias_f130814_0_0" 
#"age_e11_first_reported_non_insulin_dependent_diabetes_mellitus_f130708_0_0"
  
#"CXCL17" 
#"ASGR1" 
#"LEP"
#"IGFBP2"
  #"NCAN"

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

k_controls <- if (length(args) >= 4) as.integer(args[4]) else 5
xmin <- if (length(args) >= 5) as.numeric(args[5]) else -15
xmax <- if (length(args) >= 6) as.numeric(args[6]) else -1

protID <- prot_clean(protID_raw)

# -----------------------------
# Load Module1 OOF components
# -----------------------------
in_dir <- file.path(cfg$module1_root, covarType)

oof_fp1 <- file.path(in_dir, paste0("OOF_components_", protID_raw, "_", covarType, ".rds"))
oof_fp2 <- file.path(in_dir, paste0("OOF_components_", protID, "_", covarType, ".rds"))
oof_comp_fp <- if (file.exists(oof_fp1)) oof_fp1 else oof_fp2
if (!file.exists(oof_comp_fp)) {
  stop("Missing OOF_components file. Tried:\n  ", oof_fp1, "\n  ", oof_fp2)
}

message("Reading OOF components: ", oof_comp_fp)
oof_comp <- readRDS(oof_comp_fp) %>% as.data.frame()
stopifnot(all(c("eid","PredTotal") %in% names(oof_comp)))

if (!"Gcis" %in% names(oof_comp)) oof_comp$Gcis <- 0
if (!"Gtrans" %in% names(oof_comp)) oof_comp$Gtrans <- 0

E_cols   <- grep("^E_", names(oof_comp), value = TRUE)
GxEcis_c <- grep("^GxEcis_", names(oof_comp), value = TRUE)
GxEtr_c  <- grep("^GxEtrans_", names(oof_comp), value = TRUE)

oof_comp <- oof_comp %>%
  mutate(
    comp_G = Gcis + Gtrans,
    comp_E = if (length(E_cols) > 0) rowSums(across(all_of(E_cols)), na.rm = TRUE) else 0
  )

gx <- rep(0, nrow(oof_comp))
if (length(GxEcis_c) > 0) gx <- gx + rowSums(oof_comp[, GxEcis_c, drop = FALSE], na.rm = TRUE)
if (length(GxEtr_c)  > 0) gx <- gx + rowSums(oof_comp[, GxEtr_c,  drop = FALSE], na.rm = TRUE)
oof_comp$comp_GxE <- gx
oof_comp$comp_nonG <- oof_comp$comp_E + oof_comp$comp_GxE

# -----------------------------
# Load observed protein from pxs_rds
# -----------------------------
pxs0 <- readRDS(cfg$pxs_rds)

UKBprot_df <- NULL
if (isS4(pxs0)) UKBprot_df <- as.data.frame(pxs0@UKBprot_df)
if (is.list(pxs0) && !is.null(pxs0$UKBprot_df)) UKBprot_df <- as.data.frame(pxs0$UKBprot_df)
if (is.null(UKBprot_df)) stop("Could not locate UKBprot_df in pxs object.")

if (!protID %in% names(UKBprot_df)) {
  stop("Protein column not found in UKBprot_df: ", protID,
       "\n(try passing the exact protein ID used in the object)")
}

obsP <- UKBprot_df %>% select(eid, !!protID) %>% rename(obs_protein = !!protID)

# -----------------------------
# Load MDstore (timing + covariates)
# -----------------------------
md_fp <- paste0(cfg$mdstore_prefix, covarType, ".rds")
if (!file.exists(md_fp)) stop("Missing MDstore: ", md_fp)
MD <- readRDS(md_fp)

DZ_df <- as.data.frame(MD@DZ_df)
cov_df <- as.data.frame(MD@covars_df)

needed_t2e <- c("eid","recode_age_of_assessment_0_0","recode_age_of_death_0_0","age_of_removal_0_0","age_of_lastfollowup")
miss <- setdiff(needed_t2e, names(DZ_df))
if (length(miss) > 0) stop("DZ_df missing: ", paste(miss, collapse = ", "))
if (!t2d_age_col %in% names(DZ_df)) stop("T2D age column not found in DZ_df: ", t2d_age_col)

# pull only needed covariates for matching
match_covs <- c("eid", "sex_f31_0_0", "body_mass_index_bmi_f23104_0_0")
miss2 <- setdiff(match_covs, names(cov_df))
if (length(miss2) > 0) stop("covars_df missing: ", paste(miss2, collapse = ", "))

# -----------------------------
# Assemble analysis base
# -----------------------------
df0 <- oof_comp %>%
  inner_join(obsP, by = "eid") %>%
  inner_join(cov_df[, match_covs, drop = FALSE], by = "eid") %>%
  inner_join(DZ_df[, c(needed_t2e, t2d_age_col), drop = FALSE], by = "eid")

# incident-only filtering + censor_age
df1 <- survival_time_incident_only(df0, event_age_col = t2d_age_col)

# define incident case set
cases <- df1 %>%
  filter(DZ_status == 1) %>%
  mutate(index_age = .data[[t2d_age_col]]) %>%
  filter(is.finite(index_age))

message("Incident cases available: ", nrow(cases))

if (nrow(cases) < 200) warning("Low number of incident cases. Matching may be unstable.")

# -----------------------------
# Nested case-control sampling (risk-set)
# -----------------------------
set.seed(1)

matched_list <- vector("list", nrow(cases))
match_meta <- vector("list", nrow(cases))

# Pre-split by sex to reduce work
df_by_sex <- split(df1, df1$sex_f31_0_0)

for (i in seq_len(nrow(cases))) {
  cse <- cases[i, , drop = FALSE]
  idx_age <- cse$index_age
  
  # risk-set eligible controls:
  # same sex, not yet T2D by idx_age, and still observed at idx_age
  pool0 <- df_by_sex[[as.character(cse$sex_f31_0_0)]]
  if (is.null(pool0) || nrow(pool0) == 0) next
  
  candidates <- pool0 %>%
    filter(eid != cse$eid) %>%
    filter(is.finite(body_mass_index_bmi_f23104_0_0)) %>%
    filter(is.finite(recode_age_of_assessment_0_0)) %>%
    filter(censor_age >= idx_age) %>%
    filter(is.na(.data[[t2d_age_col]]) | .data[[t2d_age_col]] > idx_age)
  
  if (nrow(candidates) == 0) next
  
  sel <- match_controls_for_case(
    case_row = cse,
    candidates = candidates,
    k = k_controls,
    age_relax_seq = cfg$age_relax_seq,
    bmi_relax_seq = cfg$bmi_relax_seq
  )
  if (is.null(sel)) next
  
  # create matched set
  set_id <- paste0("set_", i)
  case_rec <- cse %>%
    mutate(set_id = set_id, role = "case", index_age = idx_age,
           age_caliper_used = NA_real_, bmi_caliper_used = NA_real_, k_shortfall = 0)
  
  ctrl_rec <- sel %>%
    mutate(set_id = set_id, role = "control", index_age = idx_age)
  
  matched_list[[i]] <- bind_rows(case_rec, ctrl_rec)
  
  match_meta[[i]] <- tibble(
    set_id = set_id,
    case_eid = cse$eid,
    index_age = idx_age,
    n_controls = sum(ctrl_rec$role == "control"),
    age_caliper_used = ctrl_rec$age_caliper_used[1],
    bmi_caliper_used = ctrl_rec$bmi_caliper_used[1],
    k_shortfall = if ("k_shortfall" %in% names(ctrl_rec)) max(ctrl_rec$k_shortfall, na.rm = TRUE) else 0
  )
}

matched_df <- bind_rows(matched_list)
meta_df <- bind_rows(match_meta)

if (nrow(matched_df) == 0) stop("No matched sets were created. Check t2d_age_col and covariates availability.")

message("Matched sets created: ", n_distinct(matched_df$set_id))
message("Total rows in matched dataset: ", nrow(matched_df))

# -----------------------------
# Compute time axis and restrict to pre-diagnosis window
# (both cases and controls are aligned to case index date)
# -----------------------------
matched_df <- matched_df %>%
  mutate(
    years_to_index = index_age - recode_age_of_assessment_0_0,
    x_years = -years_to_index  # negative = before diagnosis
  ) %>%
  filter(is.finite(x_years)) %>%
  filter(x_years >= xmin, x_years <= xmax)

# z-score outcomes within matched analytic sample (keeps case/control comparable)
matched_df <- matched_df %>%
  mutate(
    obsP_z = zscore_vec(obs_protein),
    compG_z = zscore_vec(comp_G),
    compNonG_z = zscore_vec(comp_nonG)
  )

# -----------------------------
# Bin and summarize: mean ± SE by role and time bin
# -----------------------------
bin_width <- cfg$bin_width_years
matched_df <- matched_df %>%
  mutate(bin = floor(x_years / bin_width) * bin_width) %>%
  # ensure bins are ordered numeric
  mutate(bin = as.numeric(bin))

summ_bin <- function(df, ycol) {
  df %>%
    group_by(role, bin) %>%
    summarise(
      n = n(),
      mean = mean(.data[[ycol]], na.rm = TRUE),
      se = sd(.data[[ycol]], na.rm = TRUE) / sqrt(n),
      .groups = "drop"
    ) %>%
    mutate(y = ycol)
}

binned <- bind_rows(
  summ_bin(matched_df, "obsP_z"),
  summ_bin(matched_df, "compG_z"),
  summ_bin(matched_df, "compNonG_z")
)

# -----------------------------
# Simple trend test (like “slope difference” idea)
# Fit mean ~ bin for each role; compare slopes (interaction)
# (Uses binned means as response; lightweight and reviewer-readable)
# -----------------------------
trend_test <- function(binned_df, yname) {
  dd <- binned_df %>% filter(y == yname) %>% filter(is.finite(mean), is.finite(bin))
  # keep bins that exist for both roles
  keep_bins <- dd %>% count(bin) %>% filter(n == 2) %>% pull(bin)
  dd <- dd %>% filter(bin %in% keep_bins)
  
  if (nrow(dd) < 8) {
    return(tibble(y = yname, p_slope_diff = NA_real_, note = "too_few_bins"))
  }
  
  # weighted by n (or 1/se^2). We'll use n to avoid infinities.
  fit <- lm(mean ~ bin * role, data = dd, weights = n)
  a <- anova(fit)
  # interaction term tests slope difference
  # term name can vary depending on role encoding; grab the last row
  p <- a$`Pr(>F)`[nrow(a)]
  tibble(y = yname, p_slope_diff = as.numeric(p), note = "bin_mean_lm_bin_by_role")
}

trend_tbl <- bind_rows(
  trend_test(binned, "obsP_z"),
  trend_test(binned, "compG_z"),
  trend_test(binned, "compNonG_z")
)

# -----------------------------
# Plot function (mean ± SE with lines)
# -----------------------------
plot_panel <- function(binned_df, yname, title, ylab) {
  dd <- binned_df %>% filter(y == yname)
  
  ggplot(dd, aes(x = bin, y = mean, group = role, linetype = role)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 1.8) +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.25, linewidth = 0.6) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    scale_x_continuous(breaks = seq(xmin, xmax, by = 2)) +
    labs(
      title = title,
      x = "Time to diagnosis (years; negative = before diagnosis)",
      y = ylab,
      linetype = NULL
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      legend.position = "right"
    )
}

pA <- plot_panel(binned, "obsP_z",
                 title = paste0("A  Observed protein (", protID_raw, ")"),
                 ylab = "Observed protein (z)")

pB <- plot_panel(binned, "compG_z",
                 title = "B  Genetic component (Gcis + Gtrans)",
                 ylab = "Genetic component (z)")

pC <- plot_panel(binned, "compNonG_z",
                 title = "C  Non-genetic component (E + GxE)",
                 ylab = "Non-genetic component (z)")

fig <- gridExtra::arrangeGrob(pA, pB, pC, ncol = 1)

# -----------------------------
# Save outputs
# -----------------------------
out_dir <- file.path(cfg$out_root, covarType, protID_raw)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

stub <- paste0("NCC_", covarType, "__", protID_raw, "__", t2d_age_col,
               "__k", k_controls, "__", xmin, "_to_", xmax)

fig_png <- file.path(out_dir, paste0(stub, ".png"))
fig_pdf <- file.path(out_dir, paste0(stub, ".pdf"))
matched_tsv <- file.path(out_dir, paste0(stub, ".matched_sets.tsv"))
meta_tsv <- file.path(out_dir, paste0(stub, ".match_meta.tsv"))
binned_tsv <- file.path(out_dir, paste0(stub, ".binned.tsv"))
trend_tsv <- file.path(out_dir, paste0(stub, ".trend_tests.tsv"))

ggsave(fig_png, fig, width = 10, height = 12, dpi = 300)
ggsave(fig_pdf, fig, width = 10, height = 12)

fwrite(as.data.table(matched_df), matched_tsv, sep = "\t")
fwrite(as.data.table(meta_df), meta_tsv, sep = "\t")
fwrite(as.data.table(binned), binned_tsv, sep = "\t")
fwrite(as.data.table(trend_tbl), trend_tsv, sep = "\t")

message("Saved:\n  ", fig_png, "\n  ", fig_pdf,
        "\n  ", matched_tsv,
        "\n  ", meta_tsv,
        "\n  ", binned_tsv,
        "\n  ", trend_tsv)
message("DONE.")