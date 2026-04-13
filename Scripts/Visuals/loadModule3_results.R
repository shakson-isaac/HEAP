#Libraries:
library(data.table)
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(ggpmisc)
library(pbapply)
library(dplyr)
library(purrr)
library(qs) #For faster save/load of rds object.
setwd("/n/groups/patel/shakson_ukb/UK_Biobank/")

#Figures:
# General overall averaged result of mediation
# Variation of GEM statistic
# Highlight Specific Diseases (HR and c-index)
# Highlight E v G partitioning (Plot) - Bootstrapping Procedure (Script #2: mediationboot.R)
# Highlight variability across covariate specifications

#Tables:
#Save significant summary stats as excel files for tables.

#Tips:
#'HAVE geom_point come after geom_ribbon to make sure the points POP OUT!!!
#'Next usage: Use p-value instead of HR for each protein - counting associations*


#'*Load Files of Mediation Results*
loadMDAssoc <- function(covarType){
  error_files <- list()  # Track files w/ errors
  stat_all <- list()
  
  # Use pblapply for progress bar and iteration
  for(i in 1:1000){
    tryCatch({
      load <- fread(file=paste0("/n/groups/patel/shakson_ukb/UK_Biobank/Data/Parallel/Module3/",
                                covarType,"/MDres_",i,".txt"))
      
      stat_all[[i]] <- load
      
    }, error = function(e) {
      message(paste("Error reading file:", i))
      error_files <<- append(error_files, i)
    })
  }
  
  
  stat_all <- do.call("rbind", stat_all)
  # stat_all <- stat_all %>%
  #   mutate(E_HRi = exp(`Exposure Indirect Effect`),
  #          G_HRi = exp(`Genetic Indirect Effect`),
  #          delta_HRi = exp(abs(`Exposure Indirect Effect` - `Genetic Indirect Effect`)),
  #          total_HRi = exp(`Exposure Indirect Effect` + `Genetic Indirect Effect`),
  #          Eprop_mediated = abs(`Exposure Indirect Effect`)/(abs(`Exposure Indirect Effect`) + abs(`Exposure Direct Effect`)),
  #          Gprop_mediated = abs(`Genetic Indirect Effect`)/(abs(`Genetic Indirect Effect`) + abs(`Genetic Direct Effect`)))
  # stat_all$CovarSpec <- covarType
  # stat_all$HR_upper95 <- as.numeric(stat_all$HR_upper95)
  
  
  return(stat_all)
}

# Load results:
Type1 <- loadMDAssoc(covarType = "Type1")
Type2 <- loadMDAssoc(covarType = "Type2")
Type3 <- loadMDAssoc(covarType = "Type3")
Type4 <- loadMDAssoc(covarType = "Type4")
Type5 <- loadMDAssoc(covarType = "Type5")

fwrite(Type1, file = "/n/groups/patel/shakson_ukb/UK_Biobank/Output/Mediation/Results/Type1.csv")
fwrite(Type2, file = "/n/groups/patel/shakson_ukb/UK_Biobank/Output/Mediation/Results/Type2.csv")
fwrite(Type3, file = "/n/groups/patel/shakson_ukb/UK_Biobank/Output/Mediation/Results/Type3.csv")
fwrite(Type4, file = "/n/groups/patel/shakson_ukb/UK_Biobank/Output/Mediation/Results/Type4.csv")
fwrite(Type5, file = "/n/groups/patel/shakson_ukb/UK_Biobank/Output/Mediation/Results/Type5.csv")



colnames(Type5)
head(Type5)

#'*Obtain OLS slope across Specifications:*
#'* Obtain OLS slope across Specifications (logHR version)
#'   - Regress exposure-mediated indirect effect on genetic (cis) indirect effect
#'   - Uses: pxs_NIE_logHR ~ gcis_NIE_logHR
oSlope <- function(MDres) {
  req <- c("prot_term", "pxs_NIE_logHR", "gcis_NIE_logHR", "protein_p", "DZ_ID")
  miss <- setdiff(req, names(MDres))
  if (length(miss) > 0) stop("Missing columns in MDres: ", paste(miss, collapse = ", "))
  
  protContext <- pbapply::pblapply(unique(MDres$prot_term), function(p) {
    
    df0 <- dplyr::filter(MDres, prot_term == p)
    
    # Complete cases for regression
    df_fit <- df0 %>%
      dplyr::filter(is.finite(pxs_NIE_logHR), is.finite(gcis_NIE_logHR))
    
    # If too few points or no variation in X, slope isn't estimable
    if (nrow(df_fit) < 2 || isTRUE(stats::sd(df_fit$gcis_NIE_logHR) == 0)) {
      df_sig <- df0 %>% dplyr::filter(protein_p < 0.05 / nrow(MDres))
      num_diseases <- length(unique(df_sig$DZ_ID))
      
      max_gcis <- if (nrow(df_fit) > 0) max(abs(df_fit$gcis_NIE_logHR), na.rm = TRUE) else NA_real_
      max_pxs  <- if (nrow(df_fit) > 0) max(abs(df_fit$pxs_NIE_logHR),  na.rm = TRUE) else NA_real_
      
      return(data.frame(
        ID = p,
        Estimate = NA_real_,
        Std.Error = NA_real_,
        Pvalue = NA_real_,
        max_gcis = max_gcis,
        max_pxs = max_pxs,
        Estimate_Norm = NA_real_,
        NumDiseases = num_diseases,
        N_fit = nrow(df_fit)
      ))
    }
    
    fit <- stats::lm(pxs_NIE_logHR ~ gcis_NIE_logHR, data = df_fit)
    stats_mat <- summary(fit)$coefficients
    
    # Guard: sometimes still only intercept row after aliasing
    if (nrow(stats_mat) < 2) {
      slope <- se <- pval <- NA_real_
    } else {
      slope <- stats_mat[2, "Estimate"]
      se    <- stats_mat[2, "Std. Error"]
      pval  <- stats_mat[2, "Pr(>|t|)"]
    }
    
    # Max magnitudes on the same sample used for regression
    max_gcis <- max(abs(df_fit$gcis_NIE_logHR), na.rm = TRUE)
    max_pxs  <- max(abs(df_fit$pxs_NIE_logHR),  na.rm = TRUE)
    
    # Significant disease count (kept as you wrote it)
    df_sig <- df0 %>% dplyr::filter(protein_p < 0.05 / nrow(MDres))
    num_diseases <- length(unique(df_sig$DZ_ID))
    
    # Normalized slope: scale to compare typical ranges (dimensionless-ish)
    est_norm <- if (is.finite(slope) && is.finite(max_gcis) && is.finite(max_pxs) && max_pxs > 0) {
      slope * (max_gcis / max_pxs)
    } else {
      NA_real_
    }
    
    data.frame(
      ID = p,
      Estimate = slope,
      Std.Error = se,
      Pvalue = pval,
      max_gcis = max_gcis,
      max_pxs = max_pxs,
      Estimate_Norm = est_norm,
      NumDiseases = num_diseases,
      N_fit = nrow(df_fit)
    )
  })
  
  do.call(rbind, protContext)
}

Type1_o <- oSlope(Type1)
Type2_o <- oSlope(Type2)
Type3_o <- oSlope(Type3)
Type4_o <- oSlope(Type4)
Type5_o <- oSlope(Type5)


Type5_T2D <- Type5 %>%
                filter(pxs_NIE_logHR_delta_p < 0.05/nrow(Type5) |
                       gcis_NIE_logHR_delta_p < 0.05/nrow(Type5) |
                       gtrans_NIE_logHR_delta_p < 0.05/nrow(Type5)) %>%
                filter(DZ_ID == "age_e11_first_reported_non_insulin_dependent_diabetes_mellitus_f130708_0_0")



topCisEffects <- Type5_T2D %>% 
                    arrange(gcis_NIE_logHR_delta_p) %>%
                    select(c("prot_term", "DZ_ID",
                             "gcis_NIE_HR_delta_l95",
                             "gcis_NIE_HR_delta_u95",
                             "gcis_NIE_logHR_delta_p",
                             "gtrans_NIE_HR_delta_l95",
                             "gtrans_NIE_HR_delta_u95",
                             "gtrans_NIE_logHR_delta_p",
                             "pxs_NIE_HR_delta_l95",
                             "pxs_NIE_HR_delta_u95",
                             "pxs_NIE_logHR_delta_p"))
colnames(topCisEffects)

topTransEffects <- Type5_T2D %>% 
  arrange(gtrans_NIE_logHR_delta_p) %>%
  select(c("prot_term", "DZ_ID",
           "gcis_NIE_HR_delta_l95",
           "gcis_NIE_HR_delta_u95",
           "gcis_NIE_logHR_delta_p",
           "gtrans_NIE_HR_delta_l95",
           "gtrans_NIE_HR_delta_u95",
           "gtrans_NIE_logHR_delta_p",
           "pxs_NIE_HR_delta_l95",
           "pxs_NIE_HR_delta_u95",
           "pxs_NIE_logHR_delta_p"))

topPXSEffects <- Type5_T2D %>% 
  arrange(pxs_NIE_logHR_delta_p) %>%
  select(c("prot_term", "DZ_ID",
           "gcis_NIE_HR_delta_l95",
           "gcis_NIE_HR_delta_u95",
           "gcis_NIE_logHR_delta_p",
           "gtrans_NIE_HR_delta_l95",
           "gtrans_NIE_HR_delta_u95",
           "gtrans_NIE_logHR_delta_p",
           "pxs_NIE_HR_delta_l95",
           "pxs_NIE_HR_delta_u95",
           "pxs_NIE_logHR_delta_p"))




#Quick check:
Type1_SigE_NIE <- Type1 %>% filter(pxs_NIE_logHR_delta_p < 0.05/nrow(Type1))
Type1_SigGcis_NIE <- Type1 %>% filter(gcis_NIE_logHR_delta_p < 0.05/nrow(Type1))
Type1_SigGtrans_NIE <- Type1 %>% filter(gtrans_NIE_logHR_delta_p < 0.05/nrow(Type1))


#Quick check (v2):
Type5_SigE_NIE <- Type5 %>% filter(pxs_NIE_logHR_delta_p < 0.05/nrow(Type5))
Type5_SigGcis_NIE <- Type5 %>% filter(gcis_NIE_logHR_delta_p < 0.05/nrow(Type5))
Type5_SigGtrans_NIE <- Type5 %>% filter(gtrans_NIE_logHR_delta_p < 0.05/nrow(Type5))


colnames(Type1)

plot(Type1$gcis_NIE_HR, Type1$pxs_NIE_HR)
Type1$gcis_NIE_logHR


# I HAVE A GOOD IDEA:
# MR results compared with the 'level of covariate specification adjustment'
# Another principled way of understanding associations

#Heritable 'partioned regions' of exposure gwas to 'validate' GTEX?




####### OTHER VISUALIZATIONS ####
library(data.table)

DT <- as.data.table(Type5)

alpha <- 0.05 / nrow(DT)

DTf <- DT[
  pxs_NIE_logHR_delta_p   < alpha |
    gcis_NIE_logHR_delta_p  < alpha |
    gtrans_NIE_logHR_delta_p< alpha
]

# OPTIONAL: require a minimum mediated magnitude so plots don't fill with numerical dust
DTf <- DTf[ abs(pxs_NIE_logHR) + abs(gcis_NIE_logHR) + abs(gtrans_NIE_logHR) > 0.01 ]

DTf[, `:=`(
  NIE_E     = pxs_NIE_logHR,
  NIE_cis   = gcis_NIE_logHR,
  NIE_trans = gtrans_NIE_logHR
)]

DTf[, denom := abs(NIE_E) + abs(NIE_cis) + abs(NIE_trans)]
DTf <- DTf[denom > 0]

DTf[, `:=`(
  wE   = abs(NIE_E)/denom,
  wCIS = abs(NIE_cis)/denom,
  wTR  = abs(NIE_trans)/denom,
  
  # what direction is the TOTAL mediated effect?
  nie_dir = fifelse((NIE_E + NIE_cis + NIE_trans) >= 0, "Mediated risk ↑", "Mediated risk ↓"),
  
  # evidence strength for plotting (use the best p among components you filtered on)
  p_min = pmin(pxs_NIE_logHR_delta_p, gcis_NIE_logHR_delta_p, gtrans_NIE_logHR_delta_p),
  strength = -log10(pmax(p_min, 1e-300))
)]

# "dominant driver" for coloring / annotation
DTf[, driver := c("E-dominant","Gcis-dominant","Gtrans-dominant")[
  max.col(cbind(abs(NIE_E), abs(NIE_cis), abs(NIE_trans)), ties.method="first")
]]

###### ANOTHER Ternary Plot Version #####
suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(ggtern)     # ternary plots
  library(ggalluvial) # alluvial / sankey-style
  library(scales)
})
#install.packages('ggtern')
#install.packages('ggalluvial')

DF <- Type5 %>% as_tibble()

alpha <- 0.05 / nrow(DF)

DFf <- DF %>%
  # carry forward if ANY component NIE is significant
  filter(
    pxs_NIE_logHR_delta_p    < alpha |
      gcis_NIE_logHR_delta_p   < alpha |
      gtrans_NIE_logHR_delta_p < alpha
  ) %>%
  # OPTIONAL: avoid numerical dust in the ternary
  mutate(denom_abs = abs(pxs_NIE_logHR) + abs(gcis_NIE_logHR) + abs(gtrans_NIE_logHR)) %>%
  filter(denom_abs > 0.01) %>%
  mutate(
    # weights (composition)
    wE   = abs(pxs_NIE_logHR)   / denom_abs,
    wCIS = abs(gcis_NIE_logHR)  / denom_abs,
    wTR  = abs(gtrans_NIE_logHR)/ denom_abs,
    
    # direction of TOTAL mediated effect (sum of components)
    total_NIE_sum = pxs_NIE_logHR + gcis_NIE_logHR + gtrans_NIE_logHR,
    nie_dir = if_else(total_NIE_sum >= 0, "Mediated risk ↑", "Mediated risk ↓"),
    
    # evidence strength (use the best p among components)
    p_min = pmin(pxs_NIE_logHR_delta_p, gcis_NIE_logHR_delta_p, gtrans_NIE_logHR_delta_p, na.rm = TRUE),
    strength = -log10(pmax(p_min, 1e-300)),
    
    # dominant driver (largest absolute component)
    driver = case_when(
      abs(pxs_NIE_logHR) >= abs(gcis_NIE_logHR) & abs(pxs_NIE_logHR) >= abs(gtrans_NIE_logHR) ~ "E-dominant",
      abs(gcis_NIE_logHR) >= abs(pxs_NIE_logHR) & abs(gcis_NIE_logHR) >= abs(gtrans_NIE_logHR) ~ "Gcis-dominant",
      TRUE ~ "Gtrans-dominant"
    )
  )

DFf <- DFf %>%
  mutate(
    icd_letter = str_to_upper(str_match(DZ_ID, "^age_([a-z])")[,2]),
    dz_cat = case_when(
      icd_letter %in% c("A","B") ~ "Infectious",
      icd_letter == "C" ~ "Neoplasms",
      icd_letter == "D" ~ "Blood/Immune",
      icd_letter == "E" ~ "Endocrine/Metabolic",
      icd_letter == "F" ~ "Mental",
      icd_letter == "G" ~ "Nervous",
      icd_letter == "I" ~ "Circulatory",
      icd_letter == "J" ~ "Respiratory",
      icd_letter == "K" ~ "Digestive",
      icd_letter == "M" ~ "Musculoskeletal",
      TRUE ~ "Other"
    )
  )

p_tern_driver <- ggtern(data = DFf, aes(x = wE, y = wCIS, z = wTR)) +
  geom_point(aes(color = driver, shape = nie_dir, alpha = strength), size = 1.6) +
  scale_alpha_continuous(range = c(0.15, 0.9), guide = "none") +
  labs(
    title = "Composition of mediated effects across E, Gcis, and Gtrans",
    T = "Gtrans (|NIE| share)", L = "E (|NIE| share)", R = "Gcis (|NIE| share)",
    color = "Dominant driver",
    shape = "Direction (total NIE)"
  ) +
  theme_bw() +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold")
  )

p_tern_driver

p_tern_facet <- ggtern(data = DFf, aes(x = wE, y = wCIS, z = wTR)) +
  geom_point(aes(color = driver, alpha = strength), size = 1.2) +
  scale_alpha_continuous(range = c(0.15, 0.9), guide = "none") +
  facet_wrap(~ dz_cat, ncol = 3) +
  labs(
    title = "Mediated effect composition differs by disease category",
    T = "Gtrans share", L = "E share", R = "Gcis share",
    color = "Dominant driver"
  ) +
  theme_bw() +
  theme(plot.title = element_text(face = "bold"))

p_tern_facet





##### SANKEY PLOT
DFa <- DFf %>%
  group_by(DZ_ID) %>%
  arrange(desc(strength), .by_group = TRUE) %>%
  slice_head(n = 10) %>%
  ungroup() %>%
  mutate(weight = 1)

top_prots <- DFf %>%
  count(prot_term, sort = TRUE) %>%
  slice_head(n = 30) %>%
  pull(prot_term)

DFa <- DFf %>%
  filter(prot_term %in% top_prots) %>%
  mutate(weight = 1)

p_alluvial <- ggplot(DFa,
                     aes(axis1 = dz_cat, axis2 = prot_term, axis3 = driver, y = weight)
) +
  geom_alluvium(aes(fill = driver), alpha = 0.7, width = 1/12) +
  geom_stratum(width = 1/10, alpha = 0.9) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
  scale_x_discrete(limits = c("Disease category", "Mediator protein", "Dominant driver"),
                   expand = c(.05, .05)) +
  labs(title = "Disease → mediator proteins → dominant driver of mediation") +
  theme_bw() +
  theme(
    axis.title = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold")
  )

p_alluvial

length(unique(DTf$prot_term))
length(unique(DTf$DZ_ID))

unique(str_to_upper(str_match(DTf$DZ_ID, "^age_([a-z])")[,2]))

#### REDO Ternary and Alluvial Plot ####

# ============================================================
# HEAP Mediation viz (dplyr) — FINAL FIGURE VERSION
# Panel 1: Ternary hexbin (optionally faceted by NIE direction)
# Panel 2: Alluvial (Disease category → Mediator (top N + Other) → Driver)
#
# Inputs:
#   - Type5 (data.frame / tibble / data.table)
# Required columns:
#   DZ_ID, prot_term,
#   pxs_NIE_logHR, gcis_NIE_logHR, gtrans_NIE_logHR,
#   pxs_NIE_logHR_delta_p, gcis_NIE_logHR_delta_p, gtrans_NIE_logHR_delta_p
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(ggtern)
  library(ggalluvial)
})

# ----------------------------
# 0) Config
# ----------------------------
CFG <- list(
  alpha_global = 0.05,     # your Bonferroni familywise alpha
  denom_min = 0.01,        # remove "numerical dust" in ternary
  hex_bins_main = 35,      # ternary hex resolution
  hex_bins_facet = 22,     # if you facet ternary
  alluvial_topN_proteins = 25,   # show top N proteins + "Other proteins"
  alluvial_topK_diseases = 10,   # OPTIONAL: limit to top K disease categories by weight (set NULL to keep all)
  alluvial_min_weight = 1        # OPTIONAL: drop extremely tiny flows if using weight != 1
)

# ----------------------------
# 1) Start from Type5
# ----------------------------
DF <- Type5 %>% as_tibble()

alpha <- CFG$alpha_global / nrow(DF)

# ----------------------------
# 2) Carry-forward filter + composition + meta
# ----------------------------
DFf <- DF %>%
  # carry forward if ANY component NIE is significant
  filter(
    pxs_NIE_logHR_delta_p    < alpha |
      gcis_NIE_logHR_delta_p   < alpha |
      gtrans_NIE_logHR_delta_p < alpha
  ) %>%
  mutate(
    # absolute-sum for composition
    denom_abs = abs(pxs_NIE_logHR) + abs(gcis_NIE_logHR) + abs(gtrans_NIE_logHR)
  ) %>%
  # remove numerical dust
  filter(denom_abs > CFG$denom_min) %>%
  mutate(
    # ternary weights (absolute contribution shares)
    wE   = abs(pxs_NIE_logHR)   / denom_abs,
    wCIS = abs(gcis_NIE_logHR)  / denom_abs,
    wTR  = abs(gtrans_NIE_logHR)/ denom_abs,
    
    # direction of total mediated effect
    total_NIE_sum = pxs_NIE_logHR + gcis_NIE_logHR + gtrans_NIE_logHR,
    nie_dir = if_else(total_NIE_sum >= 0, "Mediated risk ↑", "Mediated risk ↓"),
    
    # evidence strength: best p among components (faithful to your rule)
    p_min = pmin(pxs_NIE_logHR_delta_p, gcis_NIE_logHR_delta_p, gtrans_NIE_logHR_delta_p, na.rm = TRUE),
    strength = -log10(pmax(p_min, 1e-300)),
    
    # dominant driver: which component has largest absolute NIE
    driver = case_when(
      abs(pxs_NIE_logHR)   >= abs(gcis_NIE_logHR) & abs(pxs_NIE_logHR)   >= abs(gtrans_NIE_logHR) ~ "E-dominant",
      abs(gcis_NIE_logHR)  >= abs(pxs_NIE_logHR)  & abs(gcis_NIE_logHR)  >= abs(gtrans_NIE_logHR) ~ "Gcis-dominant",
      TRUE ~ "Gtrans-dominant"
    )
  ) %>%
  # ICD10 letter -> coarse chapter categories
  mutate(
    icd_letter = str_to_upper(str_match(DZ_ID, "^age_([a-z])")[,2]),
    dz_cat = case_when(
      icd_letter %in% c("A","B") ~ "Infectious",
      icd_letter == "D" ~ "Blood/Immune",
      icd_letter == "E" ~ "Endocrine/Metab",
      icd_letter == "F" ~ "Mental",
      icd_letter == "G" ~ "Nervous",
      icd_letter == "H" ~ "Eye/Ear",
      icd_letter == "I" ~ "Circulatory",
      icd_letter == "J" ~ "Respiratory",
      icd_letter == "K" ~ "Digestive",
      icd_letter == "L" ~ "Skin",
      icd_letter == "M" ~ "MSK",
      icd_letter == "N" ~ "Genitourinary",
      TRUE ~ "Other"
    )
  )

message("After filtering: n = ", nrow(DFf),
        " rows; proteins = ", dplyr::n_distinct(DFf$prot_term),
        "; diseases = ", dplyr::n_distinct(DFf$DZ_ID))

# ============================================================
# Panel 1A: Ternary HEXBIN (single panel)
# ============================================================
p_tern_hex <- ggtern(DFf, aes(x = wE, y = wCIS, z = wTR)) +
  geom_hex_tern(bins = CFG$hex_bins_main) +
  labs(
    title = "Composition of mediated effects across E, Gcis, and Gtrans",
    L = "E (|NIE| share)",
    R = "Gcis (|NIE| share)",
    T = "Gtrans (|NIE| share)"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "right"
  )

print(p_tern_hex)

# ============================================================
# Panel 1B: Ternary HEXBIN faceted by NIE direction (optional)
# ============================================================
p_tern_hex_dir <- ggtern(DFf, aes(x = wE, y = wCIS, z = wTR)) +
  geom_hex_tern(bins = CFG$hex_bins_facet) +
  facet_wrap(~nie_dir) +
  labs(
    title = "Composition of mediated effects by direction of total NIE",
    L = "E (|NIE| share)",
    R = "Gcis (|NIE| share)",
    T = "Gtrans (|NIE| share)"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "none"
  )

print(p_tern_hex_dir)

# ============================================================
# (SUPPLEMENT) Faceted ternary by disease category (show only 4–6 big categories)
# - Choose categories with enough points
# ============================================================
top_cats <- DFf %>%
  count(dz_cat, sort = TRUE) %>%
  slice_head(n = 6) %>%
  pull(dz_cat)

DFf_supp <- DFf %>%
  mutate(dz_cat_supp = if_else(dz_cat %in% top_cats, dz_cat, "Other")) %>%
  mutate(dz_cat_supp = factor(dz_cat_supp, levels = c(top_cats, "Other")))

p_tern_hex_cat_supp <- ggtern(DFf_supp, aes(x = wE, y = wCIS, z = wTR)) +
  geom_hex_tern(bins = CFG$hex_bins_facet) +
  facet_wrap(~dz_cat_supp, ncol = 3) +
  labs(
    title = "Mediated effect composition by disease category (top categories; others collapsed)",
    L = "E share", R = "Gcis share", T = "Gtrans share"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "none",
    # reduce clutter in small facets
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )

print(p_tern_hex_cat_supp)

# ============================================================
# Panel 2: Alluvial (Disease category → Mediator (top N + Other) → Driver)
# Strategy:
#   - choose top N proteins overall by number of links (or by summed strength)
#   - map all other proteins to "Other proteins"
#   - aggregate counts (or weights)
# ============================================================

# Pick top proteins by # links (robust + easy to interpret)
top_prots <- DFf %>%
  count(prot_term, sort = TRUE) %>%
  slice_head(n = CFG$alluvial_topN_proteins) %>%
  pull(prot_term)

# Build aggregated flow table
DFa <- DFf %>%
  mutate(mediator = if_else(prot_term %in% top_prots, prot_term, "Other proteins")) %>%
  # weight choices:
  #   - weight = 1 counts links
  #   - OR weight = abs(total_NIE_sum) to emphasize magnitude
  mutate(weight = 1L) %>%
  count(dz_cat, mediator, driver, wt = weight, name = "weight") %>%
  filter(weight >= CFG$alluvial_min_weight)

# OPTIONAL: limit to top K disease categories by total flow, if you want less spaghetti
if (!is.null(CFG$alluvial_topK_diseases)) {
  keep_cats <- DFa %>%
    group_by(dz_cat) %>%
    summarise(tw = sum(weight), .groups = "drop") %>%
    arrange(desc(tw)) %>%
    slice_head(n = CFG$alluvial_topK_diseases) %>%
    pull(dz_cat)
  
  DFa <- DFa %>%
    mutate(dz_cat2 = if_else(dz_cat %in% keep_cats, dz_cat, "Other")) %>%
    select(-dz_cat) %>%
    rename(dz_cat = dz_cat2)
}

# Ensure factors look nice: order disease categories by total weight
dz_levels <- DFa %>%
  group_by(dz_cat) %>%
  summarise(tw = sum(weight), .groups = "drop") %>%
  arrange(desc(tw)) %>%
  pull(dz_cat)

DFa <- DFa %>%
  mutate(
    dz_cat = factor(dz_cat, levels = dz_levels),
    driver = factor(driver, levels = c("E-dominant", "Gcis-dominant", "Gtrans-dominant"))
  )

# Alluvial plot
p_alluvial <- ggplot(DFa,
                     aes(axis1 = dz_cat, axis2 = mediator, axis3 = driver, y = weight)) +
  geom_alluvium(aes(fill = driver), alpha = 0.75, width = 1/14) +
  geom_stratum(width = 1/10, fill = "white", color = "black") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
  scale_x_discrete(limits = c("Disease category", "Mediator", "Dominant driver"),
                   expand = c(.03, .03)) +
  labs(
    title = "Disease category → mediator proteins → dominant driver of mediation",
    fill = "Driver"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold"),
    axis.title = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank()
  )

print(p_alluvial)

# ============================================================
# OPTIONAL: Save outputs
# ============================================================
# ggsave("panel1_ternary_hex.png", p_tern_hex, width = 8.5, height = 6, dpi = 300)
# ggsave("panel1_ternary_hex_by_dir.png", p_tern_hex_dir, width = 8.5, height = 4, dpi = 300)
# ggsave("supp_ternary_hex_by_cat.png", p_tern_hex_cat_supp, width = 10, height = 7, dpi = 300)
# ggsave("panel2_alluvial.png", p_alluvial, width = 12, height = 6.5, dpi = 300)


#### Alluvial Plot Version ####
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(ggalluvial)
})

# ----------------------------
# 0) Settings
# ----------------------------
CFG <- list(
  alpha_global = 0.05,
  denom_min = 0.01,      # remove tiny composition rows
  K_modules = 12,        # try 8–20; 12 is a good start for readability
  seed = 1,
  profile_weight = c("count", "strength")[1],  # "count" or "strength"
  alluvial_weight = c("count", "strength")[1], # "count" or "strength"
  keep_top_dz_cat = 10   # set NULL to keep all categories
)

# ----------------------------
# 1) Build DFf (filtered links) from Type5
# ----------------------------
DF <- Type5 %>% as_tibble()
alpha <- CFG$alpha_global / nrow(DF)

DFf <- DF %>%
  filter(
    pxs_NIE_logHR_delta_p    < alpha |
      gcis_NIE_logHR_delta_p   < alpha |
      gtrans_NIE_logHR_delta_p < alpha
  ) %>%
  mutate(
    denom_abs = abs(pxs_NIE_logHR) + abs(gcis_NIE_logHR) + abs(gtrans_NIE_logHR)
  ) %>%
  filter(denom_abs > CFG$denom_min) %>%
  mutate(
    # dominant driver for this link
    driver = case_when(
      abs(pxs_NIE_logHR)   >= abs(gcis_NIE_logHR) & abs(pxs_NIE_logHR)   >= abs(gtrans_NIE_logHR) ~ "E-dominant",
      abs(gcis_NIE_logHR)  >= abs(pxs_NIE_logHR)  & abs(gcis_NIE_logHR)  >= abs(gtrans_NIE_logHR) ~ "Gcis-dominant",
      TRUE ~ "Gtrans-dominant"
    ),
    # evidence strength (for optional weighting)
    p_min = pmin(pxs_NIE_logHR_delta_p, gcis_NIE_logHR_delta_p, gtrans_NIE_logHR_delta_p, na.rm = TRUE),
    strength = -log10(pmax(p_min, 1e-300)),
    # total NIE direction (optional later)
    total_NIE_sum = pxs_NIE_logHR + gcis_NIE_logHR + gtrans_NIE_logHR,
    nie_dir = if_else(total_NIE_sum >= 0, "Mediated risk ↑", "Mediated risk ↓"),
    # ICD10 letter -> chapter
    icd_letter = str_to_upper(str_match(DZ_ID, "^age_([a-z])")[,2]),
    dz_cat = case_when(
      icd_letter %in% c("A","B") ~ "Infectious",
      icd_letter == "D" ~ "Blood/Immune",
      icd_letter == "E" ~ "Endocrine/Metab",
      icd_letter == "F" ~ "Mental",
      icd_letter == "G" ~ "Nervous",
      icd_letter == "H" ~ "Eye/Ear",
      icd_letter == "I" ~ "Circulatory",
      icd_letter == "J" ~ "Respiratory",
      icd_letter == "K" ~ "Digestive",
      icd_letter == "L" ~ "Skin",
      icd_letter == "M" ~ "MSK",
      icd_letter == "N" ~ "Genitourinary",
      TRUE ~ "Other"
    )
  )

message("Filtered links: ", nrow(DFf),
        " | proteins: ", n_distinct(DFf$prot_term),
        " | diseases: ", n_distinct(DFf$DZ_ID),
        " | dz_cats: ", n_distinct(DFf$dz_cat))

# ----------------------------
# 2) Build a protein profile matrix across disease categories
#    (This is what we cluster on.)
# ----------------------------
# Choose profile weights for clustering
DF_prof <- DFf %>%
  mutate(w_prof = if (CFG$profile_weight == "strength") strength else 1) %>%
  count(prot_term, dz_cat, wt = w_prof, name = "w") %>%
  group_by(prot_term) %>%
  mutate(w = w / sum(w)) %>%  # normalize each protein to a distribution over categories
  ungroup()

# wide matrix: rows = proteins, cols = dz_cat, values = normalized weights
W <- DF_prof %>%
  tidyr::pivot_wider(names_from = dz_cat, values_from = w, values_fill = 0)

prot_ids <- W$prot_term
X <- W %>% select(-prot_term) %>% as.matrix()

# optional: stabilize clustering by compressing outliers (rarely needed)
# X <- sqrt(X)

# ----------------------------
# 3) Cluster proteins into K modules (fast)
# ----------------------------
set.seed(CFG$seed)
km <- kmeans(X, centers = CFG$K_modules, nstart = 50, iter.max = 200)

prot2mod <- tibble(
  prot_term = prot_ids,
  module_id = km$cluster
) %>%
  mutate(module = sprintf("Module %02d", module_id))

# ----------------------------
# 4) Attach module to each link
# ----------------------------
DFm <- DFf %>%
  inner_join(prot2mod, by = "prot_term")

# ----------------------------
# 5) Make the ALLUVIAL table (aggregated; no spaghetti)
# ----------------------------
# Optionally keep only top disease categories (by total weight)
DFa0 <- DFm %>%
  mutate(w_all = if (CFG$alluvial_weight == "strength") strength else 1)

if (!is.null(CFG$keep_top_dz_cat)) {
  keep_cats <- DFa0 %>%
    group_by(dz_cat) %>%
    summarise(tw = sum(w_all), .groups = "drop") %>%
    arrange(desc(tw)) %>%
    slice_head(n = CFG$keep_top_dz_cat) %>%
    pull(dz_cat)
  
  DFa0 <- DFa0 %>%
    mutate(dz_cat = if_else(dz_cat %in% keep_cats, dz_cat, "Other"))
}

# aggregate flows: Disease category -> Module -> Driver
DFa <- DFa0 %>%
  count(dz_cat, module, driver, wt = w_all, name = "weight")

# order factors nicely
dz_levels <- DFa %>%
  group_by(dz_cat) %>%
  summarise(tw = sum(weight), .groups = "drop") %>%
  arrange(desc(tw)) %>%
  pull(dz_cat)

mod_levels <- DFa %>%
  group_by(module) %>%
  summarise(tw = sum(weight), .groups = "drop") %>%
  arrange(desc(tw)) %>%
  pull(module)

DFa <- DFa %>%
  mutate(
    dz_cat = factor(dz_cat, levels = dz_levels),
    module = factor(module, levels = mod_levels),
    driver = factor(driver, levels = c("E-dominant", "Gcis-dominant", "Gtrans-dominant"))
  )

# ----------------------------
# 6) Plot: clean module alluvial
# ----------------------------
p_alluvial_modules <- ggplot(DFa,
                             aes(axis1 = dz_cat, axis2 = module, axis3 = driver, y = weight)
) +
  geom_alluvium(aes(fill = driver), alpha = 0.80, width = 1/14) +
  geom_stratum(width = 1/10, fill = "white", color = "black") +
  geom_text(
    stat = "stratum",
    aes(label = after_stat(stratum)),
    size = 3
  ) +
  scale_x_discrete(
    limits = c("Disease category", "Protein mediation module", "Dominant driver"),
    expand = c(.03, .03)
  ) +
  labs(
    title = "Disease category → protein mediation modules → dominant driver of mediation",
    fill = "Driver"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold"),
    axis.title = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    legend.position = "right"
  )

print(p_alluvial_modules)

# ----------------------------
# 7) Module annotation table (for naming modules in caption/legend)
# ----------------------------
# (a) top proteins per module by total "support"
mod_top_prots <- DFm %>%
  mutate(w = if (CFG$alluvial_weight == "strength") strength else 1) %>%
  group_by(module, prot_term) %>%
  summarise(w = sum(w), .groups = "drop") %>%
  group_by(module) %>%
  arrange(desc(w), .by_group = TRUE) %>%
  slice_head(n = 8) %>%
  summarise(
    top_proteins = paste(prot_term, collapse = ", "),
    .groups = "drop"
  )

# (b) driver composition per module
mod_driver <- DFm %>%
  mutate(w = if (CFG$alluvial_weight == "strength") strength else 1) %>%
  group_by(module, driver) %>%
  summarise(w = sum(w), .groups = "drop") %>%
  group_by(module) %>%
  mutate(frac = w / sum(w)) %>%
  arrange(module, desc(frac))

# (c) disease-category composition per module
mod_dz <- DFm %>%
  mutate(w = if (CFG$alluvial_weight == "strength") strength else 1) %>%
  group_by(module, dz_cat) %>%
  summarise(w = sum(w), .groups = "drop") %>%
  group_by(module) %>%
  mutate(frac = w / sum(w)) %>%
  arrange(module, desc(frac)) %>%
  group_by(module) %>%
  slice_head(n = 3) %>%
  summarise(
    top_dz = paste0(dz_cat, " (", round(100*frac), "%)", collapse = "; "),
    .groups = "drop"
  )

module_summary <- tibble(module = levels(DFa$module)) %>%
  left_join(mod_top_prots, by = "module") %>%
  left_join(mod_dz, by = "module") %>%
  left_join(
    mod_driver %>%
      group_by(module) %>%
      summarise(
        driver_mix = paste0(driver, "=", round(100*frac), "%", collapse = "; "),
        .groups = "drop"
      ),
    by = "module"
  )

print(module_summary)

# OPTIONAL: write module annotations for later manual naming
# write.csv(module_summary, "protein_mediation_modules_summary.csv", row.names = FALSE)

# OPTIONAL: save plot
# ggsave("alluvial_modules.png", p_alluvial_modules, width = 12, height = 6.5, dpi = 300)


######## ALLUVIAL (BETTER VERSION 1) ####
