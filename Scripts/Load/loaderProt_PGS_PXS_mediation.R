#Libraries:
library(data.table)
library(tidyverse)
library(ggplot2)
library(ggpmisc)

#TASKS:
#Remove proteins that do not have both a PXS and GIS score!

#Make a datastructure that is loaded once:
#Mediation analysis is occuring with proteins one-by-one
#Main Variables: PGS, PXS, GIS, protein measurement,  + covariates
#DataStructure:
# Prot1:
# eid, Prot1_expression, its PGS, PXS, GIS

# Common Structure:
# Covariates:
# Same version across proteins


# Common Structure:
# Diseases: 


##### PXS and GIS Score LOADING ####

load_PXSscores <- function(load, covarType){
  if(load == T){
    #'*Obtain PXS and GIS scores for all proteins*
    omiclist <- scan(file = "/n/groups/patel/shakson_ukb/UK_Biobank/BScripts/ProtPGS_PXS/OMICPREDproteins.txt",
                     as.character())
    
    PXS_list <- list()
    GIS_list <- list()
    
    for(p in omiclist){
      PXS <- read.csv(file=paste0("/n/groups/patel/shakson_ukb/UK_Biobank/Data/Parallel/Prot_PGSPXS_final/",covarType,"/","PXS_",p,".csv"))
      
      #Label variables and store
      E_name <- paste0(p,"_E")
      GxE_name <- paste0(p,"_GxE")
      
      PXS[[E_name]] <- PXS$E
      PXS[[GxE_name]] <- PXS$GxE
      
      PXS_list[[p]] <- PXS %>% select(all_of(c("eid",E_name)))
      GIS_list[[p]] <- PXS %>% select(all_of(c("eid",GxE_name)))
    }
    
    
    PXS_prot <- PXS_list %>% reduce(full_join)
    GIS_prot <- GIS_list %>% reduce(full_join)
    
    
    write.csv(PXS_prot, file = paste0("/n/groups/patel/shakson_ukb/UK_Biobank/Data/Parallel/Prot_PGSPXS_final/",covarType,"/","PXSprot.csv"), row.names = F)
    write.csv(GIS_prot, file = paste0("/n/groups/patel/shakson_ukb/UK_Biobank/Data/Parallel/Prot_PGSPXS_final/",covarType,"/","GISprot.csv"), row.names = F)
    
    return(PXS_prot)
  }
  
}
Type5 <- load_PXSscores(load = T, covarType = "Type5")
Type4 <- load_PXSscores(load = T, covarType = "Type4")
Type3 <- load_PXSscores(load = T, covarType = "Type3")
Type2 <- load_PXSscores(load = T, covarType = "Type2")
Type1 <- load_PXSscores(load = T, covarType = "Type1")


#### FIRST Occurence Disease Loading ####

# Setup:
# Load the Function and Libraries needed for UKB DataLoader
setwd("/n/groups/patel/shakson_ukb/UK_Biobank/RScripts/Extract_Raw/Finalized")
source("./dataloader_functions_parallel.R")
setwd("/n/groups/patel/shakson_ukb/UK_Biobank/RScripts/ICD10_Time2Event")
source("./time2event_functions.R")


#'*Source Project ID: *
projID = 52887
load_project(projID)
UKBdict <- fread(file = paste0("/n/groups/patel/shakson_ukb/UK_Biobank/Data/Paths/",projID,"/allpaths.txt"))
#^Use dictionary to guide selection of variables!!!

# LOAD:
fastload_baseline_info <- function(icd_df){
  #age_of_death:  f40007
  #date of death: f40000
  #ukbiobankassessment_centre: f54
  #ageassessment_center: f21003
  #dateassessment_center: f53
  #year_of_birth: f34
  #month_of_birth: f52
  #date_lost_to_followup: f191
  #reason_lost_to_followup: f190
  UKBdict <- fread(file = paste0("/n/groups/patel/shakson_ukb/UK_Biobank/Data/Paths/",projID,"/allpaths.txt"))
  baseline_metrics <- c(40007, 40000, 54, 21003, 53, 34, 52, 191, 190)
  baseline_df <- fast_dataloader_viafield(UKBdict, UKBfieldIDs = baseline_metrics, directoryInfo)
  baseline_df <- UKB_instances(baseline_df, "_0_")
  baseline_df <- baseline_df %>% mutate(month_of_birth_f52_0_0 = recode(month_of_birth_f52_0_0,
                                                                        January = 1, February = 2, March = 3,
                                                                        April = 4, May = 5, June = 6,
                                                                        July = 7, August = 8, September = 9,
                                                                        October = 10, November = 11, December = 12
  ))
  
  #Obtain Year-month and ("15" as day) Birth Date 
  baseline_df$birth_date <- as.Date(with(baseline_df, paste(year_of_birth_f34_0_0,
                                                            month_of_birth_f52_0_0, 
                                                            "15", sep="-")),"%Y-%m-%d")
  
  Time2Event_list <- list(icd_df, baseline_df)
  Time2Event_df <- Time2Event_list %>% reduce(full_join, by = "eid")
  return(Time2Event_df)
}

load_disease_T2E <- function(ICD10version, category, diseaseID, diseaseAGE){
  icd_df <- load_icd10_matrix(version = ICD10version, specific_category = category) #technically this can be done by field also (no need to load entire matrix)
  T2E_df <- fastload_baseline_info(icd_df)
  T2E_df <- time2event_ages(T2E_df)
  T2E_df <- ICD10_ages(version = ICD10version, T2E_df, disease_code = diseaseID,
                       disease_recode = diseaseAGE)
  T2E_df_subset <- T2E_df %>% select(all_of(c("eid","recode_age_of_assessment_0_0",
                                              "recode_age_of_death_0_0",
                                              "age_of_removal_0_0",
                                              "age_of_lastfollowup",
                                              paste0(diseaseAGE,"_0_0"))))
  return(T2E_df_subset)
}

# LOAD ENTIRE DISEASE FIRST OCCURENCES Table:
first_occur_categories <- c(2401,2403,2404,2405,
                            2406, 2407, 2408, 2409,
                            2410, 2411, 2412, 2413,
                            2414, 2415, 2416, 2417)

First_Occur_DZids <- UKBdict %>%
  filter(Category %in% first_occur_categories) %>%
  pull(FieldID)

#All First Occurence Diseases!
DZdf <- fast_dataloader_viafield(UKBdict, UKBfieldIDs = First_Occur_DZids, directoryInfo)
#^^Took about 10 minutes to complete!!^^

#Add baseline info!
T2E_df <- fastload_baseline_info(DZdf)

#Get Ages of Death, Last Followup, Removal, etc.
T2E_df <- time2event_ages(T2E_df)

#Get IDs for Date:
DZid_date <- UKBdict %>%
  filter(Category %in% first_occur_categories) %>%
  filter(grepl("date", descriptive_colnames)) %>%
  pull(descriptive_colnames)

#Get systematic ID for DZ age:
DZid_age <- gsub("date","age",DZid_date)


###Make sure to make this a function: Obtain Ages for 1st occurrence disease:
newdata <- T2E_df

# Updated Function to FIND Age of First Occurence Disease:
age_of_DZ <- function(field, recode, df){
  date_of_INST <- field
  age_of_INST <- recode
  df[[date_of_INST]] <- as.Date(df[[date_of_INST]])
  df[[age_of_INST]] <- round(as.numeric((df[[date_of_INST]] - df$birth_date) /365), digits = 1)
  return(df)
}

for(i in 1:length(DZid_date)){
  newdata <- age_of_DZ(DZid_date[i],
                       DZid_age[i], newdata)
}
#2-3 minutes

DZ_age <- newdata %>%
            select(all_of(c("eid",DZid_age,
                            "recode_age_of_assessment_0_0",
                            "recode_age_of_death_0_0",
                            "age_of_removal_0_0",
                            "age_of_lastfollowup")))

saveRDS(DZ_age, file = "/n/scratch/users/s/shi872/UKB_intermediate/UKB_DZage_load.rds")
#DZ_ageload <- readRDS(file = "/n/scratch/users/s/shi872/UKB_intermediate/UKB_DZage_load.rds")


#### Protein, Respective Scores, + Covariate Loading ####

#LOAD DATA used for PXS construction
PXSconstruct <- setClass(
  "PXSconstruct",
  slots = c(
    Elist = "list", #E dataframes separated by category:
    Elist_names = "character", #Names of the category:
    Eid_cat = "data.frame", #Dataframe with environmental IDs and corresponding category.
    ordinalIDs = "character", #list environmental ordinal names
    
    UKBprot_df = "data.frame", #Protein variables dataframe
    protIDs = "character", #protein variables
    
    covars_df = "data.frame", #Covariates Dataframe
    covars_list = "character" #covars variables
  )
)
PXSloader <- readRDS(file = "/n/scratch/users/s/shi872/UKB_intermediate/UKB_PGS_PXS_load.rds")


#PXS and GIS scores:
PXSscores <- fread("/n/groups/patel/shakson_ukb/UK_Biobank/Data/Parallel/Prot_PGSPXS_final/Type5/PXSprot.csv")
GISscores <- fread("/n/groups/patel/shakson_ukb/UK_Biobank/Data/Parallel/Prot_PGSPXS_final/Type5/GISprot.csv")

#E and GxE scores:
ScoresCombo <- merge(PXSscores, GISscores, by = "eid")
colnames(ScoresCombo) <- gsub("-","_",colnames(ScoresCombo)) #Make column names consistent for protIDs

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


omiclist <- scan(file = "/n/groups/patel/shakson_ukb/UK_Biobank/BScripts/ProtPGS_PXS/OMICPREDproteins.txt",
                 as.character())


#Remove specific proteins:
#Find proteins with no E or GxE score
NA_cols <- apply(ScoresCombo, 2, function(col) all(is.na(col)))
NA_prot <- colnames(ScoresCombo)[NA_cols]
NA_prot <- gsub("_GxE","",NA_prot)
NA_prot <- gsub("_E","",NA_prot)
omiclist <- omiclist[!(omiclist %in% NA_prot)]

library(progress)
# Set up progress bar
n <- length(omiclist)
# pb <- progress_bar$new(
#   format = "  [:bar] :percent :elapsed",
#   total = n, clear = FALSE, width = 60
# )

pb <- progress_bar$new(total = n)
UKBmediation_list <- list()
for(protID in omiclist){
  pb$tick()
  
  #Extract genetic scores
  omicGS <- extract_protGS(protID)
  
  #Convert protID to basename
  protID <- gsub("-", "_", protID)
  protEid <- paste0(protID,"_E")
  protGxEid <- paste0(protID,"_GxE")
  
  #Subset specific columns of interest: proteomics + covariates
  ProtScores <- ScoresCombo %>% select(all_of(c("eid",protEid,protGxEid)))
  UKBprot <- PXSloader@UKBprot_df %>% select(all_of(c("eid",protID)))
  MD_list <- list(omicGS, ProtScores, UKBprot, PXSloader@covars_df)
  
  MDdf <- suppressMessages({
    MD_list %>% reduce(inner_join)
  })
  MDdf <- as.data.frame(MDdf)
  
  #Z-score normalize protein, protein related scores, and covariates:
  skip_cols <- c("eid")
  numeric_cols <- setdiff(names(MDdf)[sapply(MDdf, is.numeric)], "eid")
  #MDdf <- as.data.frame(MDdf)
  MDdf[numeric_cols] <- lapply(MDdf[numeric_cols],function(x) as.numeric(scale(x)))
  
  #Store necessary features for each protein:
  UKBmediation_list[[protID]] <- MDdf
}

#tail(names(UKBmediation_list))
#tail(omiclist)
saveRDS(UKBmediation_list, file = "/n/scratch/users/s/shi872/UKB_intermediate/UKBmediation_load.rds")


















#Distributions time to event for diseases:
hist(MDdf %>% 
       filter(status_age_e11_first_reported_non_insulin_dependent_diabetes_mellitus_f130708_0_0 == 1) %>%
       pull(survtime_yrs_age_e11_first_reported_non_insulin_dependent_diabetes_mellitus_f130708_0_0)
)




#'*THE CODE can have a defined protID + everything else to improve runtime*
#'*HOWEVER, the combination of the protDF with disease has to be done separately since*
#'*Individuals will be removed when they have the specific disease that occurs at/before assessment 0*




####ALT ALT ALT
#ListProt
#Prot expr and scores
#Covarsdf
#Covarslist
#In the function need to normalize everything...


###### TRASH - EXTRA CODE #####



####'*TO SET UP THE MEDIATION!!!*
PXSmediation <- setClass(
  "PXSmediation",
  slots = c(
    
    # Elist = "list", #E dataframes separated by category:
    # Elist_names = "character", #Names of the category:
    # Eid_cat = "data.frame", #Dataframe with environmental IDs and corresponding category.
    # ordinalIDs = "character", #list environmental ordinal names
    UKBprot_Scores = "data.frame",
    
    UKBprot_df = "data.frame", #Protein variables dataframe
    protIDs = "character", #protein variables
    
    covars_df = "data.frame", #Covariates Dataframe
    covars_list = "character" #covars variables
  )
)


####TO DO:
##Final DF:
##Everything in terms of ages (one dataframe)
##Everything in terms of dates (second dataframe)
#Combine with PGS, PXS, and proteins to see final N of each disease for proteomics (for reporting)
#REMEMBER this N removes individuals who had the disease at or before assessment 0.


T2E_df$age_when_attended_assessment_centre_f21003_0_0
T2E_df$birth_date


#GET the age_firstreported:





# Function to FIND Age of Recruitment/Death Record/Censoring/etc.:
age_of_event <- function(assessment, field, recode, df){
  date_of_INST <- paste0(field, assessment)
  age_of_INST <- paste0(recode, assessment)
  df[[date_of_INST]] <- as.Date(df[[date_of_INST]])
  df[[age_of_INST]] <- round(as.numeric((df[[date_of_INST]] - df$birth_date) /365), digits = 1)
  return(df)
}


# Function to Obtain Ages to Events:
time2event_ages <- function(Time2Event_df, disease_code){
  #'*Obtain Ages of Events*
  #Assessment Center Instance 0: Age
  Time2Event_df <- age_of_event("_0_0", "date_of_attending_assessment_centre_f53", 
                                "recode_age_of_assessment", Time2Event_df)
  
  #Death: Age
  Time2Event_df <- age_of_event("_0_0", "date_of_death_f40000", 
                                "recode_age_of_death", Time2Event_df)
  
  #Censoring Left Study: Age
  Time2Event_df <- age_of_event("_0_0", "date_lost_to_follow_up_f191", 
                                "age_of_removal", Time2Event_df)
  
  #Censoring to End of Study Calculated via Modification Date and Modify UKB Category Date of ICD10 Codes
  end_of_study = as.Date("2022-10-03") #'*Double check if this is true for both dates*
  Time2Event_df[["age_of_lastfollowup"]] <- round(as.numeric((end_of_study - Time2Event_df$birth_date)/365), digits = 1)
  
  return(Time2Event_df)
}


###ONCE load consider date of disease, what type of report (source),
#whether any partipants got the disease after assessment 0.

baseline_metrics <- c(40007, 40000, 54, 21003, 53, 34, 52, 191, 190)








#Diseases of Interest: (Late-Onset)

#Cardiometabolic Traits:
#Cardiovascular Disease
#sr<-c('1075')  
#icd10<-c("I210 ", "I211", "I212", "I213", "I214", "I219", "I219",
#         "I220", "I221", "I228", "I219","I231","I232","I233","I236","I238",
#         "I241", "I252") 
#Type 2 Diabetes // Obesity
#Hypertension

#Lung Function:
#Asthma
#COPD


#Neurological:
#Alzheimers Disease & Dementia
#Schizophrenia
#Depression
#Parkinsons


#Autoimmune:
#MS (Multiple Sclerosis), 
# Rheumatoid Arthritis


#Cancer:
# Breast
# Colorectal
# Lung

#Use fieldID: 40008 for Age got Cancer:
#Use category 100092 for cancer registry:

#Other:
#Osteoporosis




#Type 2 Diabetes:
T2D_df <- load_disease_T2E(ICD10version = "ICD10_firstoccur", category = 2404,
                           diseaseID = "date_e11_first_reported_non_insulin_dependent_diabetes_mellitus_f130708",
                           diseaseAGE= "age_of_T2D")

#Alcoholic liver disease:
ALD_df <- load_disease_T2E(ICD10version = "ICD10_firstoccur", category = 2411,
                           diseaseID = "date_k70_first_reported_alcoholic_liver_disease_f131658",
                           diseaseAGE= "age_of_ALD")

#https://www.nature.com/articles/s41588-022-01199-5#Sec10
#Source for NAFLD in UKB.
NAFLD_df <- load_disease_T2E(ICD10version = "ICD10_firstoccur", category = 2411,
                             diseaseID = "date_k76_first_reported_other_diseases_of_liver_f131670",
                             diseaseAGE= "age_of_NAFLD")



