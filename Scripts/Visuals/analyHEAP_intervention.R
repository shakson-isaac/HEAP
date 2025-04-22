#Libraries:
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

#'*Load Univariate UKB HEAP results*
setwd("/n/groups/patel/shakson_ukb/UK_Biobank/")
HEAPassoc <- qread("./Output/HEAPres/HEAPassoc.qs")
#Use this: HEAPassoc@HEAPlist

#'*Load Intervention Studies:*
#'Exercise - JCI paper: https://pmc.ncbi.nlm.nih.gov/articles/PMC10132160/#sec13
#'GLP1 agonist treatment: https://doi.org/10.1038/s41591-024-03355-2

# Read: first sheet w/ header start from line 3
file_path <- "/n/groups/patel/shakson_ukb/Motrpac/Related_Data/jciinsight_prot.xlsx"
data <- read_excel(file_path, sheet = excel_sheets(file_path)[1], skip = 2)
data$se <- data$`log(10) Fold Change`/data$`t-statistic`


# Read: first sheet w/ header start from line 3
file_path <- "/n/groups/patel/shakson_ukb/Motrpac/Related_Data/GLP1_proteomics.xlsx"
STEP1 <- read_excel(file_path, sheet = excel_sheets(file_path)[2], skip = 0)
STEP2 <- read_excel(file_path, sheet = excel_sheets(file_path)[3], skip = 0)


#### ^^Figure out how to aggregate info across all specifications: ####
# Function to process each data type (Type1 vs Type6)
process_data <- function(covar_spec_list, type) {
  data_v2 <- covar_spec_list[[type]]$test[[1]] %>% 
    filter(`Pr(>|t|)` < 0.05/n()) %>%
    select(ID, omicID, Estimate) %>%
    pivot_wider(names_from = ID, values_from = Estimate) %>%
    select(where(~sum(!is.na(.)) >= 3))
  
  # Merge with Intervention Studies
  colnames(data_v2)[which(names(data_v2) == "omicID")] <- "EntrezGeneSymbol"
  
  heritage <- data %>% 
    filter(`False Discovery Rate (q-value)` < 0.05) %>%
    rename(HERITAGE_effect = `log(10) Fold Change`) %>%
    select(c(EntrezGeneSymbol, HERITAGE_effect))
  
  glp1_step1 <- STEP1 %>% 
    filter(qvalue < 0.05) %>% 
    group_by(EntrezGeneSymbol) %>%
    mutate(GLP1_effect1 = mean(effect_size),
           GLP1_se1 = max(std_error),
           GLP1_qvalue1 = max(qvalue)) %>%
    select(c(EntrezGeneSymbol, GLP1_effect1)) %>%
    unique()
  
  glp1_step2 <- STEP2 %>% 
    filter(qvalue < 0.05) %>% 
    group_by(EntrezGeneSymbol) %>%
    mutate(GLP1_effect2 = mean(effect_size),
           GLP1_se2 = max(std_error),
           GLP1_qvalue2 = max(qvalue)) %>%
    select(c(EntrezGeneSymbol, GLP1_effect2)) %>%
    unique()
  
  # Merge all datasets
  merged_data <- list(data_v2, heritage, glp1_step1, glp1_step2) %>% reduce(full_join)
  
  return(merged_data)
}

# Function to calculate correlations and p-values
calculate_correlations <- function(merged_data) {
  
  UKBint <- merged_data
  UKBint_cor <- as.data.frame(cor(UKBint[,-1], use = "pairwise.complete.obs"))
  
  UKBint_cor <- UKBint_cor %>%
    select(c(HERITAGE_effect, GLP1_effect1,
             GLP1_effect2)) %>%
    na.omit()
  
  library(psych)
  UKBint_corv2 <- corr.test(UKBint[,-1], adjust = "none")
  UKBint_pval <- as.data.frame(UKBint_corv2$p) %>%
    select(c(HERITAGE_effect, GLP1_effect1,
             GLP1_effect2)) %>%
    na.omit()
  
  UKBint_pvalv2 <- sapply(UKBint_pval, function(x){
    p.adjust(x, method="BH")
  })
  UKBint_pvalv2 <- as.data.frame(UKBint_pvalv2)
  UKBint_pvalv2$eID <- rownames(UKBint_pval)
  
  
  UKBint_pvalv2 <- UKBint_pvalv2 %>%
    pivot_longer(cols = !eID,
                 names_to = "intervention", values_to = "pval_adjust") %>%
    filter(!(eID %in% c("HERITAGE_effect", 
                        "GLP1_effect1",
                        "GLP1_effect2")))
  head(UKBint_pvalv2)
  
  
  UKBint_corval <- as.data.frame(UKBint_corv2$r)
  UKBint_corval$eID <- rownames(UKBint_corval)
  
  UKBint_corval <- UKBint_corval %>%
    select(c(eID, HERITAGE_effect, GLP1_effect1,
             GLP1_effect2)) %>%
    pivot_longer(cols = !eID,
                 names_to = "intervention", values_to = "cor") %>%
    filter(!(eID %in% c("HERITAGE_effect", 
                        "GLP1_effect1",
                        "GLP1_effect2")))
  
  UKBint_df <- list(UKBint_pvalv2,UKBint_corval) %>% 
    reduce(full_join)
  
  selectIDs <- UKBint_df %>%
    filter(pval_adjust < 0.05) %>%
    pull(eID)
  
  UKBint_df <- UKBint_df %>% filter(eID %in% selectIDs)
  
  UKBint_dfcor <- UKBint_df %>%
    select(eID, intervention, cor) %>%
    pivot_wider(names_from = intervention, values_from = cor) %>%
    column_to_rownames(var = "eID")
  
  UKBint_dfpval <- UKBint_df %>%
    select(eID, intervention, pval_adjust) %>%
    pivot_wider(names_from = intervention, values_from = pval_adjust) %>%
    column_to_rownames(var = "eID")
  
  return(list(cor_values = UKBint_dfcor, pval_adjusted = UKBint_dfpval))
}

# Code to Runt hrough all specifications and get correlations and pvalues:
corIntList <- lapply(names(HEAPassoc@HEAPlist), function(x){
  df <- process_data(HEAPassoc@HEAPlist, x)
  cor <- calculate_correlations(df)
  return(cor)
})
names(corIntList) <- names(HEAPassoc@HEAPlist)

#Get correlation and pvalue df combining all specifications:
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

# Function to create dataframe for scatter plots:
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

# Create scatList (For correlation plots)
scatList <- lapply(names(HEAPassoc@HEAPlist), function(x){
  createScatterDF(x)
})
names(scatList) <- names(HEAPassoc@HEAPlist)

#'*Make Datastructure*
#'With: corList, pList, scatList


#'*Create DataStructure*
INTconstruct <- setClass(
  "INTconstruct",
  slots = c(
    sList = "list", #List: Individual protein comparisons with interventions with HEAP
    cList = "list", #List: Correlations of exposures and interventions
    pList = "list" #List: Pvalues of Correlations btw Exposures and Interventions
  )
)

#Create the datastructure to utilize:
HEAPint <- INTconstruct(
              sList = scatList,
              cList = corList,
              pList = pList
            )

#Save RDS file of PXSloader object
gc()
class(HEAPint)
qsave(HEAPint, file = "./Output/HEAPres/HEAPint.qs")

