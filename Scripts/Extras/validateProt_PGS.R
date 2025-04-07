library(data.table)
library(tidyverse)

#'*Load MetaboProteome Data:*
UKB_Omics <- readRDS("/n/scratch/users/s/shi872/UKB_intermediate/UKB_fc_metaboproteome.rds")
gc()


#Extract Proteomics:
UKBprot_df <- UKB_Omics@data %>% select(all_of(c("eid",UKB_Omics@protIDs)))
#colSums(!is.na(UKBprot_df))

#'*LOAD Covariates:*
#Load some covariates:
setwd("/n/groups/patel/shakson_ukb/UK_Biobank/RScripts/Extract_Raw/Finalized")
source("./dataloader_functions_parallel.R")
projID = 52887
load_project(projID)
UKBdict <- fread(file = paste0("/n/groups/patel/shakson_ukb/UK_Biobank/Data/Paths/",projID,"/allpaths.txt")) #'*UKBdict should be a global variable now!!*


#COVARIATES: ex. (Sex, BMI, Fasting Time, etc.)
covars = c(21003, 31 , 23104, 74)
covars_df <- fast_dataloader_viafield(UKBdict, UKBfieldIDs = covars, directoryInfo)
covars_df <- UKB_instances(covars_df, "_0_")
#Due to Medication Info + more:
covars_df <- UKB_multiarray_handle(covars_df)
covars_list <- colnames(covars_df)[-c(1)]




#'*Extract relevant PGS score from files:*
#'*NEED TO OUTPUT WARNING INSTEAD OF ERROR FOR ANY FILES MISSING!!!*
#How to extract files with output of .sscore
extract_protGS <- function(protID){
  omicpredIDs <- fread(file = "/n/groups/patel/shakson_ukb/UK_Biobank/Data/OMICSPRED/UKB_Olink_multi_ancestry_models_val_results_portal.csv")
  dir <- "/n/scratch/users/s/shi872/UKB_intermediate/Genetics/UKBProtGS/"
  opID <- omicpredIDs$OMICSPRED_ID[match(protID,omicpredIDs$Gene)]
  
  opFile <- paste0(opID, ".sscore")
  protGS <- fread(file = paste0(dir, opFile))
  colnames(protGS) <- c("eid", paste0(protID,"_GS"))
  
  #Convert dash to underscore:
  colnames(protGS) <- gsub("-", "_", colnames(protGS))
  return(protGS)
}

#Example: LDLR - To Validate PGS utility:
#'*I had to modify the GS_struct function + more to handle dashes and replace with underscores!*
GS_struct <- function(protID){
  omic_PGS <- extract_protGS(protID)
  
  #Convert dash to underscore: for HLA genes:
  protID <- gsub("-","_",protID)
  omic_orig <- UKBprot_df %>% select(all_of(c("eid", protID)))
  
  omicGS <- list()
  omicGS[["combo"]] <- na.omit(merge(omic_PGS, omic_orig, by = "eid"))
  omicGS[["solo"]] <- na.omit(omic_PGS)
  
  return(omicGS)
}

#Build Function to run association:
#Different Data Structure:
GxE_assoc <- function(df, E_list, res_var, G_var, covariates){
  full_stat <- list()
  stat1 <- list()
  stat2 <- list()
  
  #Quick scaling procedure: mean centered
  #allow for intercept to not be significant and around 0.
  df <- na.omit(df) # Remove rows with NA values
  df <- df %>% column_to_rownames("eid")
  
  numeric_cols <- sapply(df, is.numeric)
  numeric_col_names <- names(numeric_cols[numeric_cols])
  
  # Scale only numeric columns
  df[numeric_col_names] <- scale(df[numeric_col_names])
  
  for(i in E_list){
    pred_var <- c(i, G_var, paste(G_var, "*", i), covariates)
    formula <- as.formula(paste(c(res_var), paste(pred_var, collapse="+"), sep="~"))
    fit <- lm(formula, data = df)
    res <- summary(fit)
    
    #Coefficients for lm:
    CI_Edf <- res$coefficients
    CI_Edf <- as.data.frame(CI_Edf)
    
    #R2 & N for lm:
    CI_Edf$R2 <- summary(fit)$r.squared
    CI_Edf$adj.R2 <- summary(fit)$adj.r.squared
    CI_Edf$samplesize <- length(fit$fitted.values)
    
    stat_tbl <- CI_Edf %>%
      rownames_to_column(var = "id") %>%
      pivot_longer(cols = colnames(CI_Edf), 
                   names_to = "stats", 
                   values_to = "value")
    
    stat1[[i]] <- stat_tbl %>% filter(id == i)
    stat2[[i]] <- stat_tbl %>% filter(id == paste0(i,":",G_var))
  }
  
  #Print last res to double check:
  print(res)
  
  stat1 <- do.call(rbind, stat1)
  stat1 <- stat1 %>% pivot_wider(names_from = stats, values_from = value)
  
  stat2 <- do.call(rbind, stat2)
  stat2 <- stat2 %>% pivot_wider(names_from = stats, values_from = value)
  
  full_stat <- list(stat1, stat2)
  return(full_stat)
}


#Build Function to run association:
#Different Data Structure:
#Changed structure to handle differences in amount of missing Data:
#ONLY OMITS NA's for each Environmental Variable (NOT in totality:)
#'*NEED to TAKE INTO ACCOUNT CATEOGRICAL VARIABLES as FACTORS!!*
#'*NEED to TAKE INTO ACCOUNT ORDINAL VARIABLES as ORDERED:*
#'*^^factor( , levels = c(), ordered = T)*
GxE_assoc_ver2 <- function(df, E_list, res_var, G_var, covariates){
  full_stat <- list()
  stat1 <- list()
  stat2 <- list()
  
  #Quick scaling procedure: mean centered
  #allow for intercept to not be significant and around 0.
  #df <- na.omit(df) # Remove rows with NA values
  df <- df %>% column_to_rownames("eid")
  
  numeric_cols <- sapply(df, is.numeric)
  numeric_col_names <- names(numeric_cols[numeric_cols])
  
  # Scale only numeric columns
  df[numeric_col_names] <- scale(df[numeric_col_names])
  
  for(i in E_list){
    #Setup Formula:
    pred_var <- c(i, G_var, paste(G_var, "*", i), covariates)
    formula <- as.formula(paste(c(res_var), paste(pred_var, collapse="+"), sep="~"))
    
    #Subset df:
    df2 <- df %>% select(all_of(c(i,res_var,G_var,covariates)))
    
    #df <- na.omit(df)
    
    #Linear Regression:
    fit <- lm(formula, data = df2)
    res <- summary(fit)
    
    #Coefficients for lm:
    CI_Edf <- res$coefficients
    CI_Edf <- as.data.frame(CI_Edf)
    
    #R2 & N for lm:
    CI_Edf$R2 <- summary(fit)$r.squared
    CI_Edf$adj.R2 <- summary(fit)$adj.r.squared
    CI_Edf$samplesize <- length(fit$fitted.values)
    
    stat_tbl <- CI_Edf %>%
      rownames_to_column(var = "id") %>%
      pivot_longer(cols = colnames(CI_Edf), 
                   names_to = "stats", 
                   values_to = "value")
    
    stat1[[i]] <- stat_tbl %>% filter(id == i)
    stat2[[i]] <- stat_tbl %>% filter(id == paste0(i,":",G_var))
  }
  
  #Print last res to double check:
  print(res)
  
  stat1 <- do.call(rbind, stat1)
  stat1 <- stat1 %>% pivot_wider(names_from = stats, values_from = value)
  
  stat2 <- do.call(rbind, stat2)
  stat2 <- stat2 %>% pivot_wider(names_from = stats, values_from = value)
  
  full_stat <- list(stat1, stat2)
  return(full_stat)
}


#'*To Build Later: Partitioned R2 function:*


##### Nutrient Info and Diet (Axis 1) ####

#Check nutrient info:
nutr_metab <- readRDS(file = "/n/scratch/users/s/shi872/UKB_intermediate/UKB_nutr_metab_df.rds")
#'*Need to update this later:*
nutr_list = scan(file = "/n/groups/patel/shakson_ukb/UK_Biobank/RScripts/Env_Preprocess/nutr_names.txt",
                 what = as.character())
UKBnutr_df <- nutr_metab %>% select(all_of(c("eid",nutr_list)))

#70K individuals at instance 0: according to UKBshowcase - https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=26030
sum(!is.na(UKBnutr_df$total_weight_of_all_foods_and_beverages_f26000_0_0))



library(ggplot2)
library(ggpmisc)
check_GS_vis <- function(protID){
  omicDS <- GS_struct(protID) 
  
  protID <- gsub("-", "_", protID)
  protGS <- paste0(protID,"_GS")
  
  formula <- as.formula(paste(c(protID), paste(protGS, collapse="+"), sep="~"))
  
  #Function to Test if R2 is accurate for omicGS
  fit <- lm(formula, data = omicDS$combo)
  summary(fit)
  
  # Scatter Plot of Polygenic Score vs Original Omic:
  gg1 <- ggplot(data = omicDS$combo, aes(.data[[protGS]], .data[[protID]])) + geom_point() +
    stat_poly_line() +
    stat_poly_eq() + 
    theme_minimal()
  
  #In case people want to verify w/ scaled version:
  # gg3 <- ggplot(data = as.data.frame(scale(omicDS$combo)), aes(.data[[protGS]], .data[[protID]])) + geom_point() +
  #   stat_poly_line() +
  #   stat_poly_eq() + 
  #   theme_minimal()
  
  #QQ plot of the Residuals to check normality assumption
  # Create a data frame for plotting
  residuals_df <- data.frame(residuals = fit$residuals)
  
  # Generate the QQ plot
  gg2 <- ggplot(residuals_df, aes(sample = residuals)) +
    stat_qq() +
    stat_qq_line() +
    labs(title = "QQ Plot of Residuals", x = "Theoretical Quantiles", y = "Sample Quantiles") +
    theme_minimal()
  
  print(gg1)
  print(gg2)
  #print(gg3)
}
check_GS_stat <- function(protID){
  omicDS <- GS_struct(protID) 
  
  protID <- gsub("-", "_", protID)
  protGS <- paste0(protID,"_GS")
  
  formula <- as.formula(paste(c(protID), paste(protGS, collapse="+"), sep="~"))
  
  #Function to Test if R2 is accurate for omicGS
  fit <- lm(formula, data = omicDS$combo)
  summary(fit)
  
  #Kolmogrov-Smirnov Test
  stat <- ks.test(fit$residuals, "pnorm", 
                  mean = mean(fit$residuals), 
                  sd = sd(fit$residuals))
  #print(stat)
  
  #Return p-value from kolmogrov-smirnov test:
  return(stat)
}
GWIS_run <- function(protID){
  # LDLR Example:
  omicDS <- GS_struct(protID) #LDLR
  
  protID <- gsub("-", "_", protID)
  protGS <- paste0(protID,"_GS")
  
  #Merge relevant files together:
  exp1 <- merge(omicDS$combo, UKBnutr_df, by = "eid")
  exp2 <- merge(exp1, covars_df, by = "eid")
  
  # 7.1 K: 
  ans <- GxE_assoc(df = exp1, E_list = nutr_list, res_var = protID,
                   G_var = protGS, covariates = c())
  ans2 <- GxE_assoc(df = exp2, E_list = nutr_list, res_var = protID,
                    G_var = protGS, covariates = covars_list)
  
  res <- list(ans, ans2)
  return(res)
}
GWIS_run_ver2 <- function(protID){
  # Example:
  omicDS <- GS_struct(protID) 
  
  protID <- gsub("-", "_", protID)
  protGS <- paste0(protID,"_GS")
  
  #Merge relevant files together:
  exp1 <- merge(omicDS$combo, UKBnutr_df, by = "eid")
  exp2 <- merge(exp1, covars_df, by = "eid")
  
  # 7.1 K: 
  ans <- GxE_assoc_ver2(df = exp1, E_list = nutr_list, res_var = protID,
                        G_var = protGS, covariates = c())
  ans2 <- GxE_assoc_ver2(df = exp2, E_list = nutr_list, res_var = protID,
                         G_var = protGS, covariates = covars_list)
  
  res <- list(ans, ans2)
  return(res)
}



#'*For loop to check GS:*
omicpredIDs <- fread(file = "/n/groups/patel/shakson_ukb/UK_Biobank/Data/OMICSPRED/UKB_Olink_multi_ancestry_models_val_results_portal.csv")
#omicpredIDs$Gene[duplicated(omicpredIDs$Gene)]


GS_normpval <- c()
GS_Dstat <- c()
for(i in omicpredIDs$Gene){
  s <- check_GS_stat(protID = i)
  GS_normpval[i] <- s$p.value
  GS_Dstat[i] <- s$statistic
}
GS_normpval
GS_Dstat

GS_validation <- data.frame(
  Dstat = GS_Dstat,
  pval = GS_normpval 
)

GS_validation$ProtID <- rownames(GS_validation)
write.table(GS_validation, file = "/n/groups/patel/shakson_ukb/UK_Biobank/RScripts/Pure_StatGen/Prot_ExPGS/ProtGS_residnorm_KStest.txt")


#Example Test Stats of Normality Assumption: 
check_GS_stat(protID = "APOE")
check_GS_stat(protID = "GIP")
check_GS_stat(protID = "LDLR")



#Example of really bad PGS: APOE
check_GS_vis(protID = "LDLR")
check_GS_vis(protID = "APOE")
#'*CLEAR Rare-Variant Contribution Missing to the APOE score!!*
#'*APOE_GS might be a bad instrument*
check_GS_vis(protID = "GIP")
check_GS_vis(protID = "HLA-DRA")
#'*CLEAR Rare-Variant Contribution Missing to the HLA_DRA score!!*


