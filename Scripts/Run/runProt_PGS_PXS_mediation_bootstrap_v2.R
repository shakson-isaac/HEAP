library(data.table)
library(tidyverse)
library(ggplot2)
library(ggpmisc)
library(ggrepel)
library(gridExtra)
library(survival)
library(pbapply)

#'*TO-DO*
#'RUN THIS Bootstrap Procedure on potential findings of the mediation analysis.
#'Save the MODELS for each protein <--> disease pair for further visualizations.
#'Save the STATS (Mediation indirect/direct/etc.) of the bootstrap samples to use for visualizations.
#'SPECIFY the covarType in each.
#'CREATE a list to parallelize through for each covarType 
#'Write bash script for the mediation analysis
#'
#'
#'Make MDloader for each covarType!
#'Make sure the run occurs with specified covarType
#'Allow individual or batch runs
#'20 minutes per ~ covarType5
#'2 minutes per ~ covarType1/2

####'*NEW VERSION* 
PXSmediation <- setClass(
  "PXSmediation",
  slots = c(
    Protlist = "list", #E dataframes separated by category:
    protIDs = "character", #Names of the category:
    
    DZ_df = "data.frame", #First occurence age matrix + age of assessment, death, etc.
    DZ_ids = "character", #First occurence DZ age ids
    
    covars_df = "data.frame", #Covariates Dataframe
    covars_list = "character" #covars variables
  )
)

####'*Time to event info:* from "/RScripts/ICD10_Time2Event/time2event_runner.R"
####Setup:
# Load the Function and Libraries needed for UKB DataLoader
setwd("/n/groups/patel/shakson_ukb/UK_Biobank/RScripts/Extract_Raw/Finalized")
source("./dataloader_functions_parallel.R")
setwd("/n/groups/patel/shakson_ukb/UK_Biobank/RScripts/ICD10_Time2Event")
source("./time2event_functions.R")


##### Mediation Analysis ####
#'*Mediation Code*
# Function to load disease of interest
#'*Updated based on reconstructed DZ file!!*
mediation_analysis_DZ <- function(DZ_df, DZ_ID, Orig_df){
  
  #'Determine better organization
  #'survtime_name, inittime_name,status_name are fixed
  #'T2E_ID - "age_of_T2D" & T2D_df are variable - need as input to function.
  
  
  #Convert to survival data:
  T2E_ID = DZ_ID
  T2E_name = DZ_ID   #paste0(T2E_ID,"_0_0")
  status_name = "DZ_status" #paste0("status_",T2E_name)
  survtime_name = "DZ_survtime" #paste0("survtime_",T2E_name)
  survtime_altname = "survtime_standard" #paste0("survtime_yrs_",T2E_name)
  inittime_name = "recode_age_of_assessment_0_0"
  
  #Subset DZ_df to have only one disease:
  DZ_df <- DZ_df %>%
    select(all_of(c("eid",T2E_ID,
                    "recode_age_of_assessment_0_0",
                    "recode_age_of_death_0_0",
                    "age_of_removal_0_0",
                    "age_of_lastfollowup")))
  
  # Obtain Survival Time (Right Censoring):
  icd10_df <- survival_time(DZ_df, event_age = T2E_name, recode_status = status_name, recode_survtime = survtime_name)
  
  #'*Change from age to time after assessment 0:*
  icd10_df[[survtime_altname]] = icd10_df[[survtime_name]] - icd10_df[[inittime_name]]
  
  
  #Setup the regressions:
  df <- suppressMessages({list(Orig_df, icd10_df) %>% 
      reduce(inner_join)}) #inner join is correct here: want info from both proteins and disease
  
  return(df)
}

# Function to fit cox regression
coxph_fit <- function(pred_var, res_var, df){
  formula <- as.formula(paste(c(res_var), paste(pred_var, collapse="+"), sep="~"))
  res.cox <- coxph(formula, data = df)
  #cox_df <- summary(res.cox)
  return(res.cox)
}

# Function to fit mediator and outcome model
#'*BE transparent about variables like covarslist, etc.*
mediation_models <- function(protID, covars_list, CVdf){
  
  M_model <- list()
  
  #2 regressions:
  #1: PXS --> Mediator
  #2: Mediator --> Outcome
  #Be consistent: #Rule #1 :: Keep individuals if age_assessment < age_ICD10
  #433 is dropoff..
  protID <- gsub("-", "_", protID)
  protGS <- paste0(protID,"_GS")
  protE <- paste0(protID,"_E")
  protGxE <- paste0(protID,"_GxE")
  
  res_var <- protID
  pred_var <- c(protGS, protE, protGxE, covars_list)
  
  #Model 1:
  formula <- as.formula(paste(c(res_var), paste(pred_var, collapse="+"), sep="~"))
  
  CVdf1 <- CVdf %>% select(all_of(c(res_var, pred_var)))
  M_fit <- lm(formula, data = CVdf1)
  #print(summary(M_fit))
  
  
  #Model 2:
  ind_var <- c(protID, protGS, protE, protGxE, covars_list)
  
  #'*Determine how to include age_sex interactions when the timing is age in survival data*
  #'*Maybe switch to time after assessment 0 in years*
  CVdf2 <- CVdf %>% select(all_of(c(ind_var, "survtime_standard","DZ_status")))
  #print(paste0("mediation: ",ind_var))
  
  O_fit <- coxph_fit(pred_var = ind_var, res_var = "Surv(survtime_standard, DZ_status)", df = CVdf2)
  #print(summary(O_fit))
  
  M_model[["mediator"]] <- M_fit
  M_model[["outcome"]] <- O_fit
  
  return(M_model)
}

###G computation: Counterfactuals 
#General Structure: (G,E) ---> P; P_hat --> D

# Function to create newdata for g computation
newdata_gcomp <- function(statmodel, cf_list, cf_value){
  #Obtain dataset used for fit:
  df <- model.frame(statmodel)
  
  #Initialize Counterfactual List:
  names(cf_value) <- cf_list
  
  #Initialize List and add corresponding names in lm Model:
  newlist <- vector("list", length(colnames(df)) - 1)
  
  #Do not include response variable!!
  #print(setdiff(colnames(df), colnames(df)[-1]))
  names(newlist) <- colnames(df)[-1]
  
  for(i in names(newlist)){
    if(i %in% cf_list){
      
      newlist[[i]] = cf_value[names(cf_value) == i]
      
    } else if (is.numeric(df[[i]])) {
      newlist[[i]] = mean(df[[i]]) #median
      
    } else {
      newlist[[i]] = levels(df[[i]])[[1]] #Decide to just take the 1st factor arbitrarily.
      
    }
    
  }
  return(newlist)
}

# Function for G computation
g_computation <- function(protID, M_fit, O_fit, G_m, E_m, G_o, E_o){
  #G and E are counterfactuals of PGS and PXS respectively
  #_m is for the mediator model ('vary the mediator')
  #_o is for the outcome model ('vary the predictor')
  #M_fit is the mediator model
  #O_fit is the outcome model
  
  #Naming Conventions
  protID <- gsub("-", "_", protID)
  protGS <- paste0(protID,"_GS")
  protE <- paste0(protID,"_E")
  
  #Mediator Model
  newlist1 <- newdata_gcomp(M_fit, c(protGS,protE), c(G_m,E_m))
  M_df <- expand.grid(newlist1)
  #print(M_df)
  M_hat <- predict(M_fit, newdata = M_df)
  
  #Outcome Model
  newlist2 <- newdata_gcomp(O_fit, c(protID,protGS,protE), c(M_hat,G_o,E_o))
  O_df <- expand.grid(newlist2)
  hazard <- predict(O_fit, newdata = O_df, type = "risk")
  
  return(hazard)
}

# Function to get Direct and Indirect Effects
Mcounterfactuals <- function(protID, M_fit, O_fit){
  protID <- gsub("-", "_", protID)
  protGS <- paste0(protID,"_GS")
  protE <- paste0(protID,"_E")
  
  df <- model.frame(M_fit)
  
  #Counterfactual Table:
  G_p <- mean(df[[protGS]]) + sd(df[[protGS]])
  G_m <- mean(df[[protGS]]) - sd(df[[protGS]])
  E_p <- mean(df[[protE]]) + sd(df[[protE]])
  E_m <- mean(df[[protE]]) - sd(df[[protE]])
  
  #Counterfactuals Analyzed:
  #Realized:
  #Total Effect = Total NIE + Total NDE
  #Total NIE = NIE_frac1 + NIE_frac2
  #Total NDE = NDE_frac1 + NDE_frac2
  
  #Total Effect:
  TEc1 <- g_computation(protID, M_fit, O_fit, G_p, E_p, G_p, E_p)
  TEc2 <- g_computation(protID, M_fit, O_fit, G_m, E_m, G_m, E_m)
  TE <- log(TEc1/TEc2)
  
  #Total NIE: Allow Mediator to Vary with G and E
  TNIEc1 <- g_computation(protID, M_fit, O_fit, G_p, E_p, G_p, E_p)
  TNIEc2 <- g_computation(protID, M_fit, O_fit, G_m, E_m, G_p, E_p)
  TNIE <- log(TNIEc1/TNIEc2)
  
  #Total NDE: Fix Mediator, Vary G and E on outcome
  TNDEc1 <- g_computation(protID, M_fit, O_fit, G_m, E_m, G_p, E_p)
  TNDEc2 <- g_computation(protID, M_fit, O_fit, G_m, E_m, G_m, E_m)
  TNDE <- log(TNDEc1/TNDEc2)
  
  #For E only:
  ENIEc1 <- g_computation(protID, M_fit, O_fit, G_p, E_p, G_p, E_p)
  ENIEc2 <- g_computation(protID, M_fit, O_fit, G_p, E_m, G_p, E_p)
  ENIE <- log(ENIEc1/ENIEc2)
  
  ENDEc1 <- g_computation(protID, M_fit, O_fit, G_m, E_m, G_m, E_p)
  ENDEc2 <- g_computation(protID, M_fit, O_fit, G_m, E_m, G_m, E_m)
  ENDE <- log(ENDEc1/ENDEc2)
  
  
  # For G only:
  GNIEc1 <- g_computation(protID, M_fit, O_fit, G_p, E_p, G_p, E_p)
  GNIEc2 <- g_computation(protID, M_fit, O_fit, G_m, E_p, G_p, E_p)
  GNIE <- log(GNIEc1/GNIEc2)
  
  GNDEc1 <- g_computation(protID, M_fit, O_fit, G_m, E_m, G_p, E_m)
  GNDEc2 <- g_computation(protID, M_fit, O_fit, G_m, E_m, G_m, E_m)
  GNDE <- log(GNDEc1/GNDEc2)
  
  mediation_vec <- c(TE,TNIE,TNDE,GNIE, GNDE, ENIE, ENDE)
  names(mediation_vec) <- c("Total Effect",
                            "Total Indirect Effect","Total Direct Effect",
                            "Genetic Indirect Effect","Genetic Direct Effect",
                            "Exposure Indirect Effect","Exposure Direct Effect")
  return(mediation_vec)
  
}

# Run counterfactuals for each percentile
Mpercentiles <- function(protID, M_fit, O_fit){
  HRrisk <- list()
  
  
  protID <- gsub("-", "_", protID)
  protGS <- paste0(protID,"_GS")
  protE <- paste0(protID,"_E")
  
  df <- model.frame(M_fit)
  
  #Counterfactual Table:
  percentiles <- c(1:19/20) #(1:19/20)
  G_values <- quantile(df[[protGS]], probs = percentiles)
  E_values <- quantile(df[[protE]], probs = percentiles)
  
  n = length(percentiles)
  m = length(percentiles)
  
  risk <- matrix(0, nrow = n, ncol = m)
  
  for(i in 1:length(G_values)){
    for(j in 1:length(E_values)){
      #Counterfactuals Analyzed: every combo!
      G = G_values[i]
      E = E_values[j]
      
      G_o = 0
      E_o = 0
      
      risk[i,j] <- g_computation(protID, M_fit, O_fit, G, E, G_o, E_o) #as a hazard ratio
      #risk[i,j] <- g_computation(protID, M_fit, O_fit, G, E, G, E) #as a hazard ratio
      
      #FOR indirect effect OUTCOME model G_o, E_o must remain constant!
      
      #G and E are counterfactuals of PGS and PXS respectively
      #_m is for the mediator model ('vary the mediator')
      #_o is for the outcome model ('vary the predictor')
      #function(protID, M_fit, O_fit, G_m, E_m, G_o, E_o)
    }
  }
  
  
  HRrisk[["z"]] <- risk
  HRrisk[["x"]] <- G_values
  HRrisk[["y"]] <- E_values
  return(HRrisk)
}

#Get genetic scores:
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

# Bootstrapping Procedure: Obtain Indirect and Direct Effects.
bootstrap_MD <- function(protID, covars_list, df, nsamp){
  #Add genetic score for a protein:
  #omicGS <- extract_protGS(protID)
  #df <- merge(df, omicGS, by = "eid")
  protID <- gsub("-", "_", protID)
  
  #Number of samples
  n <- nrow(df)
  
  #Get original estimate:
  MDmodel_orig <- mediation_models(protID, covars_list, df)
  MDorig <- Mcounterfactuals(protID, MDmodel_orig$mediator, MDmodel_orig$outcome)
  MDorig <- as.data.frame(t(MDorig))
  #Stats about coefficients from outcome model: (Protein By Itself)
  MDorig$HR <- exp(coefficients(MDmodel_orig$outcome))[protID]
  MDorig$R2 <- summary(MDmodel_orig$mediator)$adj.r.squared
  MDorig$cindex <- MDmodel_orig$outcome$concordance["concordance"]
  
  

  # Use pblapply progress bar for bootstrapping procedure.
  MD <- pblapply(1:nsamp, function(x){
    df_sim <- df[sample(1:n, size = n, replace = TRUE), ]
    
    MDmodel_full <- mediation_models(protID, covars_list, df_sim)
    vec <- Mcounterfactuals(protID, MDmodel_full$mediator, MDmodel_full$outcome)
    vec <- as.data.frame(t(vec))
    
    #Stats about coefficients from outcome model: (Protein By Itself)
    vec$HR <- exp(coefficients(MDmodel_full$outcome))[protID]
    #vec$HR_lower95 <- exp(confint(MDmodel_full$outcome))[protID, "2.5 %"]
    #vec$HR_upper95 <- exp(confint(MDmodel_full$outcome))[protID, "97.5 %"]
    
    #Stats about models:
    vec$R2 <- summary(MDmodel_full$mediator)$adj.r.squared
    vec$cindex <- MDmodel_full$outcome$concordance["concordance"]
    #vec$cindex_std <- MDmodel_full$outcome$concordance["std"]
    #vec$ID <- protID
    return(vec)
  })
  MD <- do.call(rbind, MD)
  
  #Get percentile bootstrap 95% CIs
  MD <- MD %>%
    summarise(across(everything(), 
                     list(
                       q2.5 = ~ quantile(., 0.025),
                       q97.5 = ~ quantile(., 0.975)
                     ),
                     .names = "{.col} {.fn}")) #{.col}_{.fn}
  
  # Combine Bootstrap CIs with original point estimate
  MD <- cbind(MDorig,MD)
  MD$ID <- protID
  
  return(MD)
}

# Define a function to compute summary statistics (FOR normal stats CV)
get_CI <- function(column){
  stat <- data.frame()
  if(is.character(column)){
    c <- unique(column)
    return(c)
  } else{
    stat <- tibble(
      mean = mean(column),
      sd = sd(column),
      n = length(column),
      se = sd / sqrt(n),
      ci_margin = qt(0.975, df = n - 1) * se,
      #ci_low = mean - qt(0.975, df = n - 1) * se, #95% confidence interval:
      #ci_high = mean + qt(0.975, df = n - 1) * se
      ci = paste0(paste0(round(mean, 4), " ± ", round(ci_margin, 4)))
    )
    return(stat$ci)
  }
}

#Original Estimate:
protMediation_orig <- function(protID, DZ_df = T2D_df, DZ_ID = "age_of_T2D"){
  #Extract genetic scores
  omicGS <- extract_protGS(protID)
  
  #Convert protID to basename
  protID <- gsub("-", "_", protID)
  protEid <- paste0(protID,"_E")
  protGxEid <- paste0(protID,"_GxE")
  
  #Subset specific columns of interest: proteomics + covariates
  ProtScores <- ScoresCombo %>% select(all_of(c("eid",protEid,protGxEid)))
  UKBprot <- UKB_Omics@data %>% select(all_of(c("eid",protID)))
  MD_list <- list(omicGS, ProtScores, UKBprot, covars_df)
  MDdf <- suppressMessages({
    MD_list %>% reduce(inner_join)
  })
  MDdf <- as.data.frame(MDdf)
  
  #Z-score normalize protein, protein related scores, and covariates:
  skip_cols <- c("eid")
  numeric_cols <- setdiff(names(MDdf)[sapply(MDdf, is.numeric)], "eid")
  #MDdf <- as.data.frame(MDdf)
  MDdf[numeric_cols] <- lapply(MDdf[numeric_cols],function(x) as.numeric(scale(x)))
  
  #Add disease related column:
  #MDdf <- mediation_analysis_T2D(MDdf)
  MDdf <- mediation_analysis_DZ(DZ_df, DZ_ID, Orig_df = MDdf)
  
  MDmodel_full <- mediation_models(protID, MDdf)
  vec <- Mcounterfactuals(protID, MDmodel_full$mediator, MDmodel_full$outcome)
  vec <- as.data.frame(t(vec))
  vec$ID <- protID
  
  return(vec)
}

# Function to run mediation analysis: ORIGINAL ESTIMATE: 
protMediation_final <- function(protID, DZ_ID, MDobject){
  
  MD_list <- list(MDobject@Protlist[[protID]], 
                  MDobject@covars_df)
  MDdf <- suppressMessages({
    MD_list %>% reduce(inner_join)
  })
  MDdf <- as.data.frame(MDdf)
  
  
  #'*Z-score normalize protein, protein related scores, and covariates:*
  skip_cols <- c("eid")
  numeric_cols <- setdiff(names(MDdf)[sapply(MDdf, is.numeric)], "eid")
  #MDdf <- as.data.frame(MDdf)
  MDdf[numeric_cols] <- lapply(MDdf[numeric_cols],function(x) as.numeric(scale(x)))
  
  MDdf <- mediation_analysis_DZ(MDobject@DZ_df, DZ_ID, Orig_df = MDdf)
  
  #Mediation model:
  MDmodel_full <- mediation_models(protID, MDobject@covars_list, MDdf)
  
  vec <- Mcounterfactuals(protID, MDmodel_full$mediator, MDmodel_full$outcome)
  vec <- as.data.frame(t(vec))
  
  #'*Add protID to maintain consistency!*
  protID <- gsub("-", "_", protID)
  #Stats about coefficients from outcome model: (Protein By Itself)
  vec$HR <- exp(coefficients(MDmodel_full$outcome))[protID]
  vec$HR_lower95 <- exp(confint(MDmodel_full$outcome))[protID, "2.5 %"]
  vec$HR_upper95 <- exp(confint(MDmodel_full$outcome))[protID, "97.5 %"]
  
  #Stats about models:
  vec$R2 <- summary(MDmodel_full$mediator)$adj.r.squared
  vec$cindex <- MDmodel_full$outcome$concordance["concordance"]
  vec$cindex_std <- MDmodel_full$outcome$concordance["std"]
  
  #Add protein ID
  vec$ID <- protID
  
  return(vec)
}

# Function to run mediation analysis: W/ Bootstrapping
protMediation_Boot <- function(protID, DZ_ID, MDobject, nsim){
  
  # Setting up the datastructure for mediation analysis
  MD_list <- list(MDobject@Protlist[[protID]], 
                  MDobject@covars_df)
  MDdf <- suppressMessages({
    MD_list %>% reduce(inner_join)
  })
  MDdf <- as.data.frame(MDdf)
  
  
  #'*Z-score normalize protein, protein related scores, and covariates:*
  skip_cols <- c("eid")
  numeric_cols <- setdiff(names(MDdf)[sapply(MDdf, is.numeric)], "eid")
  #MDdf <- as.data.frame(MDdf)
  MDdf[numeric_cols] <- lapply(MDdf[numeric_cols],function(x) as.numeric(scale(x)))
  
  MDdf <- mediation_analysis_DZ(MDobject@DZ_df, DZ_ID, Orig_df = MDdf)
  
  
  # Run Analysis: Bootstrapping Procedure:
  MDres <- bootstrap_MD(protID, MDobject@covars_list, MDdf, nsamp = nsim)
  return(MDres)
}

# Function to run mediation analysis: Combos E and G:
protMediation_combos <- function(protID, DZ_ID, MDobject){
  
  MD_list <- list(MDobject@Protlist[[protID]], 
                  MDobject@covars_df)
  MDdf <- suppressMessages({
    MD_list %>% reduce(inner_join)
  })
  MDdf <- as.data.frame(MDdf)
  
  
  #'*Z-score normalize protein, protein related scores, and covariates:*
  skip_cols <- c("eid")
  numeric_cols <- setdiff(names(MDdf)[sapply(MDdf, is.numeric)], "eid")
  #MDdf <- as.data.frame(MDdf)
  MDdf[numeric_cols] <- lapply(MDdf[numeric_cols],function(x) as.numeric(scale(x)))
  
  MDdf <- mediation_analysis_DZ(MDobject@DZ_df, DZ_ID, Orig_df = MDdf)
  
  #Mediation model:
  MDmodel_full <- mediation_models(protID, MDobject@covars_list, MDdf)
  
  mat <- Mpercentiles(protID, MDmodel_full$mediator, MDmodel_full$outcome)
  return(mat)
}

# Function to run mediation analysis across Diseases:
MultiDisease_run <- function(protID, MDobject){
  # Use pblapply for progress bar, storing results in a list
  MDres_list <- pblapply(MDobject@DZ_ids, function(i) {
    res <- protMediation_final(protID, i, MDobject)
    res$DZid <- i
    return(res)
  })
  
  # Combine the list of data frames into a single data frame
  MDres <- do.call(rbind, MDres_list)
  
  return(MDres)
}

#'*Mediation Bootstrap RUNS*

#TODO: Arguments:
# Protein Names: ex.) protlist <- c("IGFBP1","IGFBP2","CKB","FABP4","FGF21","LEP")
# covarType: ex.) covarType <- "Type1"
# dzIDs: dz <- "age_e11_first_reported_non_insulin_dependent_diabetes_mellitus_f130708_0_0"


#LOAD:
covarType = "Type1"
#covarType = "Type2"
#covarType = "Type3"
#covarType = "Type4"
#covarType = "Type5"
#Load Mediation Data:
MDloader <- readRDS(file = paste0("/n/scratch/users/s/shi872/UKB_intermediate/UKB_MDstore_",
                                  covarType,".rds"))


#Example 1: Many Proteins - 1 Disease (For Real RUN do 1000 bootstraps)
#1000 bootstraps estimated time:
#1 protein / CovarSpec4  --> 32 minutes
#Original Estimate, 2.5%, 97.5% of bootstrap to get CI
#Percentile Bootstrap Method

# Run through selected proteins:
protlist <- c("IGFBP1","IGFBP2","CKB","FABP4","FGF21","LEP")
dzlist <- c("age_e11_first_reported_non_insulin_dependent_diabetes_mellitus_f130708_0_0",
            "age_n18_first_reported_chronic_renal_failure_f132032_0_0",
            "age_i25_first_reported_chronic_ischaemic_heart_disease_f131306_0_0")
dzNames <- c("T2D","CRF","CIHD")

MDloader@DZ_ids[grepl("heart",MDloader@DZ_ids)]

runBoot_oneMD <- function(protlist, dz, nsamp = 100){
  MDres_bootstrap <- pblapply(protlist, function(i) {
    res <- protMediation_Boot(protID = i, 
                              DZ_ID = dz,
                              MDobject = MDloader,
                              nsim = nsamp) 
    res$DZid <- dz
    return(res)
  })
  
  bootstrap_res <- do.call(rbind,MDres_bootstrap)
  
  return(bootstrap_res)
}
T2Dres <- runBoot_oneMD(protlist, dzlist[1], nsamp = 10)

runBoot_multiMD <- function(protlist, dzlist, nsamp = 100){
  
  multiMDres <- pblapply(dzlist, function(i){
    res <- runBoot_oneMD(protlist, i, nsamp)
    return(res)
  })

  multiBoot <- do.call(rbind,multiMDres)
  
  return(multiBoot)
}
DZres <- runBoot_multiMD(protlist, dzlist, nsamp = 10)



#Example PLOT: 1 Disease
bootstrap_res_upd <- T2Dres %>%
  mutate(E_HRi_orig = exp(`Exposure Indirect Effect`),
         E_HRi_2.5 = exp(`Exposure Indirect Effect q2.5`),
         E_HRi_97.5 = exp(`Exposure Indirect Effect q97.5`),
         G_HRi_orig = exp(`Genetic Indirect Effect`),
         G_HRi_2.5 = exp(`Genetic Indirect Effect q2.5`),
         G_HRi_97.5 = exp(`Genetic Indirect Effect q97.5`)
  ) %>% 
  select(all_of(c("ID","DZid",
                  "E_HRi_orig","E_HRi_2.5","E_HRi_97.5",
                  "G_HRi_orig","G_HRi_2.5","G_HRi_97.5"))) %>%
  pivot_longer(
    cols = -c(ID,DZid),
    names_to = c("Type",".value"),
    names_pattern = "(E|G)_HRi_(.*)"
  )

#Figures:
bootstrap_res_upd %>%
  filter(grepl("e11",DZid)) %>%
  ggplot(aes(x = ID, y = orig, color = Type)) +
  geom_point(position = position_dodge(width = 0.6), size = 3) +  # Point estimate
  geom_errorbar(aes(ymin = `2.5`, ymax = `97.5`), width = 0.2,    # Error bars
                position = position_dodge(width = 0.6))  +
  geom_hline(yintercept = 1, color = "black", size = 0.8, linetype = "dashed") +  # Reference line at 1
  #scale_color_manual(values = c("Group1" = "blue", "Group2" = "red")) +
  scale_color_discrete(labels = c("E", "G")) +
  labs(x = "Proteins", y = "Hazard Ratio: Mediation Effects") +
  coord_flip() +
  theme_minimal()


### Example Plot: Multiple Diseases
bootstrap_res_upd <- DZres %>%
  mutate(E_HRi_orig = exp(`Exposure Indirect Effect`),
         E_HRi_2.5 = exp(`Exposure Indirect Effect q2.5`),
         E_HRi_97.5 = exp(`Exposure Indirect Effect q97.5`),
         G_HRi_orig = exp(`Genetic Indirect Effect`),
         G_HRi_2.5 = exp(`Genetic Indirect Effect q2.5`),
         G_HRi_97.5 = exp(`Genetic Indirect Effect q97.5`)
  ) %>% 
  select(all_of(c("ID","DZid",
                  "E_HRi_orig","E_HRi_2.5","E_HRi_97.5",
                  "G_HRi_orig","G_HRi_2.5","G_HRi_97.5"))) %>%
  pivot_longer(
    cols = -c(ID,DZid),
    names_to = c("Type",".value"),
    names_pattern = "(E|G)_HRi_(.*)"
  )

# Replace DiseaseIDs with Abbreviations:
bootstrap_res_upd$DZid <- dzNames[match(bootstrap_res_upd$DZid, dzlist)]




#Figures:
bootstrap_res_upd %>%
  ggplot(aes(x = ID, y = orig, group = Type, color = DZid)) +
  geom_point(aes(shape = Type),
             position = position_dodge(width = 0.6), 
             size = 3) +  # Point estimate
  geom_errorbar(aes(ymin = `2.5`, ymax = `97.5`), width = 0.2,    # Error bars
                position = position_dodge(width = 0.6))  +
  geom_hline(yintercept = 1, color = "black", size = 0.8, linetype = "dashed") +  # Reference line at 1
  #scale_color_manual(values = c("Group1" = "blue", "Group2" = "red")) +
  #scale_color_discrete(labels = c("E", "G")) +
  labs(x = "Proteins", y = "Hazard Ratio: Mediation Effects") +
  coord_flip() +
  theme_minimal() +
  facet_wrap(~DZid)

bootstrap_res_upd %>%
  ggplot(aes(x = ID, y = orig, group = Type, color = DZid)) +
  geom_point(aes(shape = Type),
             position = position_dodge(width = 0.6), 
             size = 3) +  # Point estimate
  geom_errorbar(aes(ymin = `2.5`, ymax = `97.5`), width = 0.2,    # Error bars
                position = position_dodge(width = 0.6))  +
  geom_hline(yintercept = 1, color = "black", size = 0.8, linetype = "dashed") +  # Reference line at 1
  #scale_color_manual(values = c("Group1" = "blue", "Group2" = "red")) +
  #scale_color_discrete(labels = c("E", "G")) +
  labs(x = "Proteins", y = "Hazard Ratio: Mediation Effects") +
  theme_minimal()


bootstrap_res_upd %>%
  ggplot(aes(x = DZid, y = orig, group = Type, color = ID)) +
  geom_point(aes(shape = Type),
             position = position_dodge(width = 0.6), 
             size = 3) +  # Point estimate
  geom_errorbar(aes(ymin = `2.5`, ymax = `97.5`), width = 0.2,    # Error bars
                position = position_dodge(width = 0.6))  +
  geom_hline(yintercept = 1, color = "black", size = 0.8, linetype = "dashed") +  # Reference line at 1
  #scale_color_manual(values = c("Group1" = "blue", "Group2" = "red")) +
  #scale_color_discrete(labels = c("E", "G")) +
  labs(x = "Proteins", y = "Hazard Ratio: Mediation Effects") +
  coord_flip() +
  theme_minimal() 







##### ALTERNATIVE VISUALIZATION TO EXPLORE - SEPARATE SCRIPT #######
### New Visualization: Varying E and G 
# IGFBP2_diab_mat <- protMediation_combos(protID = "IGFBP2", 
#                                   DZ_ID = "age_e11_first_reported_non_insulin_dependent_diabetes_mellitus_f130708_0_0",
#                                   MDobject = MDloader)

#3D Plot and Visualization
# persp(x = 1:19/20, 
#       y = 1:19/20,
#       z = IGFBP2_diab_mat$z,
#       xlab = "PGS (%)", 
#       ylab = "PXS (%)",
#       zlab = "Mediated: HR",
#       theta = 120, phi = 0,
#       col = "cyan", shade = 0.1,
#       ltheta = 0,
#       ticktype = "detailed")

# Create a filled contour plot
# filled.contour(IGFBP2_diab_mat$x, 
#                IGFBP2_diab_mat$y, 
#                IGFBP2_diab_mat$z, 
#                xlab = "PGS", ylab = "PXS", 
#                color.palette = terrain.colors)