# ================================
# Weighted cross-platform correlation pipeline
# ================================

# Libraries
library(data.table)
library(tidyverse)
library(ggplot2)
library(ggpmisc)
library(pbapply)
library(dplyr)
library(purrr)
library(plotly)
library(htmlwidgets)
library(readxl)
library(qs)
library(psych)
library(weights)

# ----------------
# Load HEAP results
# ----------------
setwd("/n/groups/patel/shakson_ukb/UK_Biobank/")
HEAPassoc <- qread("./Output/HEAPres/HEAPassoc.qs")
# Use: HEAPassoc@HEAPlist

# ----------------
# Load Intervention Studies
# Exercise - JCI paper: https://pmc.ncbi.nlm.nih.gov/articles/PMC10132160/#sec13
# GLP1 agonist treatment: https://doi.org/10.1038/s41591-024-03355-2
# ----------------

# Read: first sheet w/ header start from line 3
file_path <- "/n/groups/patel/shakson_ukb/Motrpac/Related_Data/jciinsight_prot.xlsx"
data <- read_excel(file_path, sheet = excel_sheets(file_path)[1], skip = 2)
data$se <- data$`log(10) Fold Change`/data$`t-statistic`

# Read GLP1 sheets
file_path <- "/n/groups/patel/shakson_ukb/Motrpac/Related_Data/GLP1_proteomics.xlsx"
STEP1 <- read_excel(file_path, sheet = excel_sheets(file_path)[2], skip = 0)
STEP2 <- read_excel(file_path, sheet = excel_sheets(file_path)[3], skip = 0)

# ----------------
# Load / define cross-platform reliability per protein (Olink vs SomaScan)
# Expect a CSV with columns: EntrezGeneSymbol, r_cross  (values in [-1, 1])
# If you don't have it yet, this block will create a placeholder with NA's.
# ----------------

prot_rel <- fread("/n/groups/patel/IGLOO/UKB/OlinkSoma/OlinkSoma.csv", skip = 3)
prot_rel <- prot_rel %>% select(c("gene_name","olink_nonnorm_corr","olink_smpnorm_corr"))
colnames(prot_rel) <- c("EntrezGeneSymbol", "r_cross", "r_crossv2")
prot_rel <- na.omit(prot_rel)

# ----------------
# Function: process_data  (now adds reliability & weights)
# ----------------
process_data <- function(covar_spec_list, type) {
  data_v2 <- covar_spec_list[[type]]$test[[1]] %>% 
    filter(`Pr(>|t|)` < 0.05/n()) %>%
    select(ID, omicID, Estimate) %>%
    pivot_wider(names_from = ID, values_from = Estimate) %>%
    select(where(~sum(!is.na(.)) >= 3))
  
  # Harmonize key
  colnames(data_v2)[which(names(data_v2) == "omicID")] <- "EntrezGeneSymbol"
  
  # Intervention aggregates
  heritage <- data %>% 
    filter(`False Discovery Rate (q-value)` < 0.05) %>%
    rename(HERITAGE_effect = `log(10) Fold Change`) %>%
    select(c(EntrezGeneSymbol, HERITAGE_effect))
  
  glp1_step1 <- STEP1 %>% 
    filter(qvalue < 0.05) %>% 
    group_by(EntrezGeneSymbol) %>%
    mutate(GLP1_effect1 = mean(effect_size)) %>%
    select(c(EntrezGeneSymbol, GLP1_effect1)) %>%
    unique()
  
  glp1_step2 <- STEP2 %>% 
    filter(qvalue < 0.05) %>% 
    group_by(EntrezGeneSymbol) %>%
    mutate(GLP1_effect2 = mean(effect_size)) %>%
    select(c(EntrezGeneSymbol, GLP1_effect2)) %>%
    unique()
  
  # Merge all datasets
  merged_data <- list(data_v2, heritage, glp1_step1, glp1_step2) %>% reduce(full_join)
  
  # Add cross-platform reliability & compute reliability weight
  merged_data <- merged_data %>%
    left_join(prot_rel, by = "EntrezGeneSymbol") %>%
    mutate(
      # --- Imputation choice for missing reliability r_cross ---
      # Option A (default, conservative): treat missing as 0 (exclude)
      r_cross = ifelse(is.na(r_cross), mean(r_cross, na.rm = TRUE), r_cross)
      # Option B: shrink to global mean
      # r_cross = ifelse(is.na(r_cross), mean(r_cross, na.rm = TRUE), r_cross)
      # Option C: if you have groups, replace with group-wise shrinkage here
    ) %>%
    mutate(
      # Weight choice: emphasize strong positive agreement; ignore negatives
      w_rel = pmax(r_cross, 0) #^2
      # Alternative: keep negatives but down-weight by magnitude
      # w_rel = r_cross^2
      # Alternative: linear rescale [-1,1] -> [0,1]
      # w_rel = (r_cross + 1) / 2
    )
  
  #print(paste0("aver corr",mean(merged_data$r_cross)))
  
  return(merged_data)
}

# ----------------
# Function: weighted correlations + BH p-values
# ----------------
calculate_correlations <- function(merged_data) {
  # Identify columns
  prot_id_col <- "EntrezGeneSymbol"
  intervention_cols <- c("HERITAGE_effect", "GLP1_effect1", "GLP1_effect2")
  drop_cols <- c(prot_id_col, "w_rel", "r_cross")
  keep_cols <- setdiff(names(merged_data), drop_cols)
  exposure_cols <- setdiff(keep_cols, intervention_cols)  # UKB exposure columns (wide)
  
  # Helper: weighted Pearson r + p via effective N
  wtd_cor_and_p <- function(x, y, w) {
    ok <- is.finite(x) & is.finite(y) & is.finite(w) & (w > 0)
    x <- x[ok]; y <- y[ok]; w <- w[ok]
    if (length(x) < 3) return(c(r = NA_real_, p = NA_real_))
    r <- suppressWarnings(wtd.cor(x, y, weight = w)[1,1])
    # Effective sample size
    neff <- (sum(w)^2) / sum(w^2)
    if (!is.finite(r) || !is.finite(neff) || neff <= 2) return(c(r = NA_real_, p = NA_real_))
    tval <- r * sqrt((neff - 2) / pmax(1e-12, 1 - r^2))
    pval <- 2 * pt(-abs(tval), df = neff - 2)
    c(r = r, p = pval)
  }
  
  # Compute weighted correlations for each exposure vs each intervention
  out <- purrr::map_dfr(exposure_cols, function(exp_col){
    purrr::map_dfr(intervention_cols, function(int_col){
      res <- wtd_cor_and_p(
        merged_data[[exp_col]],
        merged_data[[int_col]],
        merged_data$w_rel
      )
      tibble::tibble(
        eID = exp_col,
        intervention = int_col,
        cor = as.numeric(res["r"]),
        pval = as.numeric(res["p"])
      )
    })
  })
  
  # BH adjust within intervention
  out <- out %>%
    group_by(intervention) %>%
    mutate(pval_adjust = p.adjust(pval, method = "BH")) %>%
    ungroup()
  
  # Keep significant rows (like your previous logic)
  selectIDs <- out %>% filter(pval_adjust < 0.05) %>% pull(eID) %>% unique()
  out_sig <- out %>% filter(eID %in% selectIDs)
  
  # Return in the same structure you use downstream
  UKBint_dfcor <- out_sig %>%
    select(eID, intervention, cor) %>%
    tidyr::pivot_wider(names_from = intervention, values_from = cor) %>%
    tibble::column_to_rownames(var = "eID")
  
  UKBint_dfpval <- out_sig %>%
    select(eID, intervention, pval_adjust) %>%
    tidyr::pivot_wider(names_from = intervention, values_from = pval_adjust) %>%
    tibble::column_to_rownames(var = "eID")
  
  return(list(cor_values = UKBint_dfcor, pval_adjusted = UKBint_dfpval))
}

# ----------------
# Function: createScatterDF (unchanged from your code)
# ----------------
createScatterDF <- function(CSpec){
  UKBspec <- HEAPassoc@HEAPlist[[CSpec]]$test[[1]] %>% 
    filter(`Pr(>|t|)` < 0.05/n()) %>%
    select(ID, omicID, Estimate, `Std. Error`)
  colnames(UKBspec)[which(names(UKBspec) == "omicID")] <- "EntrezGeneSymbol"
  
  HERITAGE <- data %>% 
    filter(`False Discovery Rate (q-value)` < 0.05) %>%
    mutate(HERITAGE_se = `log(10) Fold Change`/`t-statistic`) %>%
    rename(HERITAGE_effect = `log(10) Fold Change`) %>%
    select(c(EntrezGeneSymbol, HERITAGE_effect, HERITAGE_se))
  
  GLP1_STEP1 <- STEP1 %>% 
    filter(qvalue < 0.05) %>% 
    group_by(EntrezGeneSymbol) %>%
    mutate(GLP1_effect1 = mean(effect_size),
           GLP1_se1 = max(std_error)) %>%
    select(c(EntrezGeneSymbol, GLP1_effect1, GLP1_se1)) %>%
    unique()
  
  GLP1_STEP2 <- STEP2 %>% 
    filter(qvalue < 0.05) %>% 
    group_by(EntrezGeneSymbol) %>%
    mutate(GLP1_effect2 = mean(effect_size),
           GLP1_se2 = max(std_error)) %>%
    select(c(EntrezGeneSymbol, GLP1_effect2, GLP1_se2)) %>%
    unique()
  
  UKBscatter <- list(UKBspec, HERITAGE, GLP1_STEP1, GLP1_STEP2) %>% reduce(full_join)
  return(UKBscatter)
}

# ----------------
# Run through all specifications and get weighted correlations & p-values
# ----------------
corIntList <- lapply(names(HEAPassoc@HEAPlist), function(x){
  df <- process_data(HEAPassoc@HEAPlist, x)
  cor <- calculate_correlations(df)
  return(cor)
})
names(corIntList) <- names(HEAPassoc@HEAPlist)

# Combine across specifications
corList <- lapply(names(corIntList), function(x){
  df <- corIntList[[x]]$cor_values
  df$ID <- rownames(df)
  df$Type <- x
  return(df)
})
pList <- lapply(names(corIntList), function(x){
  df <- corIntList[[x]]$pval_adjusted
  df$ID <- rownames(df)
  df$Type <- x
  return(df)
})
names(corList) <- names(HEAPassoc@HEAPlist)
names(pList) <- names(HEAPassoc@HEAPlist)

# Create scatter data for plots
scatList <- lapply(names(HEAPassoc@HEAPlist), function(x){
  createScatterDF(x)
})
names(scatList) <- names(HEAPassoc@HEAPlist)

# ----------------
# Data structure
# ----------------
INTconstruct <- setClass(
  "INTconstruct",
  slots = c(
    sList = "list", # Individual protein comparisons with interventions with HEAP
    cList = "list", # Correlations of exposures and interventions
    pList = "list"  # P-values of correlations between exposures and interventions
  )
)

HEAPint <- INTconstruct(
  sList = scatList,
  cList = corList,
  pList = pList
)

# Save
gc()
class(HEAPint)
qsave(HEAPint, file = "./Output/HEAPres/HEAPintv2.qs")

# ----------------
# Notes:
# - To change weight behavior, edit the 'w_rel' line in process_data().
# - To keep proteins with negative cross-platform correlation but down-weight them,
#   use w_rel = r_cross^2 instead of pmax(r_cross, 0)^2.
# - To include precision (SE) later, multiply w_rel by a precision factor for each pair.
# ----------------


