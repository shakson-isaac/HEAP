library(data.table)
library(tidyverse)
library(ggplot2)
library(ggpmisc)
library(ggrepel)
library(gridExtra)
library(survival)
library(pbapply)

#'*CHECKLIST*
#'Make MDloader for each covarType!
#'Make sure the run occurs with specified covarType
#'Allow individual or batch runs


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


# Bootstrap to determine mediation metrics: indirect and direct effects.
bootstrap_MD <- function(protID, df, nsamp){
  #Add genetic score for a protein:
  #omicGS <- extract_protGS(protID)
  #df <- merge(df, omicGS, by = "eid")
  #protID <- gsub("-", "_", protID)
  
  #Number of samples
  n <- nrow(df)
  
  #Store mediation
  MD <- data.frame()
  
  #Run mediation nsamp times
  for(i in 1:nsamp){
    df_sim <- df[sample(1:n, size = n, replace = TRUE), ]
    
    MDmodel_full <- mediation_models(protID, df_sim)
    vec <- Mcounterfactuals(protID, MDmodel_full$mediator, MDmodel_full$outcome)
    vec <- as.data.frame(t(vec))
    MD <- rbind(MD, vec)
  }
  
  MD$ID <- protID
  return(MD)
}

# Define a function to compute summary statistics
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

#With Bootstrapping
protMediation <- function(protID, DZ_df = T2D_df, DZ_ID = "age_of_T2D", nsim){
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
  MDdf <- MD_list %>% reduce(inner_join)
  MDdf <- as.data.frame(MDdf)
  
  #Z-score normalize protein, protein related scores, and covariates:
  skip_cols <- c("eid")
  numeric_cols <- setdiff(names(MDdf)[sapply(MDdf, is.numeric)], "eid")
  #MDdf <- as.data.frame(MDdf)
  MDdf[numeric_cols] <- lapply(MDdf[numeric_cols],function(x) as.numeric(scale(x)))
  
  #Add disease related column:
  #MDdf <- mediation_analysis_T2D(MDdf)
  MDdf <- mediation_analysis_DZ(DZ_df, DZ_ID, Orig_df = MDdf)
  
  
  #Run analysis:
  MDres <- bootstrap_MD(protID, MDdf, nsamp = nsim)
  return(MDres)
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

# Function to run mediation analysis: 
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


#'*Mediation Parallelization*
args = commandArgs(trailingOnly = T)
idx <- as.integer(args[1])
split_num <- as.integer(args[2])
covarType <- as.character(args[3])

#idx = 1
#split_num = 1000
#covarType = "Type1"

#Load Mediation Data:
MDloader <- readRDS(file = paste0("/n/scratch/users/s/shi872/UKB_intermediate/UKB_MDstore_",
                       covarType,".rds"))

#Get omiclist to parallelize through
omiclist <- MDloader@protIDs

# Create a grouping factor to split the vector 
groups <- cut(seq_along(omiclist), breaks = split_num, labels = FALSE)
split_vectors <- split(omiclist, groups)

# Function to run Mediation in batches and specific covariate type.
runMediation <- function(protlist, MDobject, covarType, idx){
  #Create folder if nonexistent
  dir.create(
    paste0("/n/groups/patel/shakson_ukb/UK_Biobank/Data/Parallel/Prot_PGSPXS_mediation/",covarType),
    showWarnings = FALSE
  )
  
  MDres_save <- pblapply(protlist, function(p) {
    res <- MultiDisease_run(protID = p, MDobject)
    return(res)
  })
  MDfinal <- do.call(rbind, MDres_save)
  
  write.table(MDfinal,
              file = paste0("/n/groups/patel/shakson_ukb/UK_Biobank/Data/Parallel/Prot_PGSPXS_mediation/",
                            covarType,"/MDres_",idx,".txt"),
              row.names = F)
}
runMediation(split_vectors[[idx]], MDloader, covarType, idx)

