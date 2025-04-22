#Libraries:
library(data.table)
library(tidyverse)
library(ggplot2)
library(ggpmisc)

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
Type5 <- load_PXSscores(load = F, covarType = "Type5")
Type4 <- load_PXSscores(load = F, covarType = "Type4")
Type3 <- load_PXSscores(load = F, covarType = "Type3")
Type2 <- load_PXSscores(load = F, covarType = "Type2")
Type1 <- load_PXSscores(load = F, covarType = "Type1")
#'*Run above once with TRUE*


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


# Updated Function to FIND Age of First Occurence Disease:
age_of_DZ <- function(field, recode, df){
  date_of_INST <- field
  age_of_INST <- recode
  df[[date_of_INST]] <- as.Date(df[[date_of_INST]])
  df[[age_of_INST]] <- round(as.numeric((df[[date_of_INST]] - df$birth_date) /365), digits = 1)
  return(df)
}

#Obtain Ages for 1st occurrence disease:
for(i in 1:length(DZid_date)){
  T2E_df <- age_of_DZ(DZid_date[i],
                       DZid_age[i], T2E_df)
}

DZ_age <- T2E_df %>%
            select(all_of(c("eid",DZid_age,
                            "recode_age_of_assessment_0_0",
                            "recode_age_of_death_0_0",
                            "age_of_removal_0_0",
                            "age_of_lastfollowup")))


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


# Function to obtain PXS and GIS Scores for each protein:
PXSscores <- function(covarType){
  #PXS and GIS scores:
  PXSscores <- fread(paste0("/n/groups/patel/shakson_ukb/UK_Biobank/Data/Parallel/Prot_PGSPXS_final/",
                            covarType,"/PXSprot.csv"))
  GISscores <- fread(paste0("/n/groups/patel/shakson_ukb/UK_Biobank/Data/Parallel/Prot_PGSPXS_final/",
                            covarType,"/GISprot.csv"))
  
  #E and GxE scores:
  ScoresCombo <- merge(PXSscores, GISscores, by = "eid")
  colnames(ScoresCombo) <- gsub("-","_",colnames(ScoresCombo)) #Make column names consistent for protIDs
  return(ScoresCombo)
}

#'*Load Scores from Different Covariate Specifications*
ScoresCombo_Type1 <- PXSscores("Type1")
ScoresCombo_Type2 <- PXSscores("Type2")
ScoresCombo_Type3 <- PXSscores("Type3")
ScoresCombo_Type4 <- PXSscores("Type4")
ScoresCombo_Type5 <- PXSscores("Type5")


# Function to Obtain Genetic Scores for each protein:
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

library(progress)
library(pbapply)
# Function to Return List of Protein - Score Dataframes
ProtScoreList <- function(ScoresCombo){
  #Subset proteins with score info:
  omiclist <- scan(file = "/n/groups/patel/shakson_ukb/UK_Biobank/BScripts/ProtPGS_PXS/OMICPREDproteins.txt",
                   as.character())
  
  #Remove specific proteins:
  #Find proteins with no E or GxE score
  NA_cols <- apply(ScoresCombo, 2, function(col) all(is.na(col)))
  NA_prot <- colnames(ScoresCombo)[NA_cols]
  NA_prot <- gsub("_GxE","",NA_prot)
  NA_prot <- gsub("_E","",NA_prot)
  omiclist <- omiclist[!(omiclist %in% NA_prot)]
  
  UKBprotlist <- pblapply(omiclist, function(protID) {
    #Extract genetic scores
    omicGS <- extract_protGS(protID)
    
    #Convert protID to basename
    protID <- gsub("-", "_", protID)
    protEid <- paste0(protID,"_E")
    protGxEid <- paste0(protID,"_GxE")
    
    #Subset specific columns of interest: proteomics + covariates
    ProtScores <- ScoresCombo %>% select(all_of(c("eid",protEid,protGxEid)))
    UKBprot <- PXSloader@UKBprot_df %>% select(all_of(c("eid",protID)))
    Protlist <- list(omicGS, ProtScores, UKBprot)
    
    df <- suppressMessages({
      Protlist %>% reduce(inner_join)
    })
    df <- as.data.frame(df)
    
    return(df)
  })
  names(UKBprotlist) <- omiclist
  
  return(UKBprotlist)
}

UKBprotlist_Type1 <- ProtScoreList(ScoresCombo_Type1)
UKBprotlist_Type2 <- ProtScoreList(ScoresCombo_Type2)
UKBprotlist_Type3 <- ProtScoreList(ScoresCombo_Type3)
UKBprotlist_Type4 <- ProtScoreList(ScoresCombo_Type4)
UKBprotlist_Type5 <- ProtScoreList(ScoresCombo_Type5)




#Create R Object to store relevant info!
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

# Function to specify covariate specifications:
PXScovarSpec <- function(PXSdata, covars_subset){
  PXSdata@covars_list <- covars_subset
  PXSdata@covars_df <- PXSdata@covars_df %>% select(all_of(c("eid",covars_subset)))
  
  return(PXSdata)
}

CovarSpec <- list()
CovarSpec[["Type1"]] <- c("age_when_attended_assessment_centre_f21003_0_0",
                          "sex_f31_0_0")
CovarSpec[["Type2"]] <- c("age_when_attended_assessment_centre_f21003_0_0",
                          "sex_f31_0_0",
                          "body_mass_index_bmi_f23104_0_0",
                          "fasting_time_f74_0_0")
CovarSpec[["Type3"]] <- c("age_when_attended_assessment_centre_f21003_0_0",
                          "sex_f31_0_0",
                          "age2","age_sex", "age2_sex",
                          "body_mass_index_bmi_f23104_0_0",
                          "fasting_time_f74_0_0",
                          "uk_biobank_assessment_centre_f54_0_0",
                          paste0("genetic_principal_components_f22009_0_",1:20)
)
CovarSpec[["Type4"]] <- c("age_when_attended_assessment_centre_f21003_0_0",
                          "sex_f31_0_0",
                          "age2","age_sex", "age2_sex",
                          "body_mass_index_bmi_f23104_0_0",
                          "fasting_time_f74_0_0",
                          "uk_biobank_assessment_centre_f54_0_0",
                          paste0("genetic_principal_components_f22009_0_",1:20),
                          "combined_Blood_pressure_medication",
                          "combined_Hormone_replacement_therapy",
                          "combined_Oral_contraceptive_pill_or_minipill",
                          "combined_Insulin",
                          "combined_Cholesterol_lowering_medication",
                          "combined_Do_not_know",
                          "combined_None_of_the_above",
                          "combined_Prefer_not_to_answer"
)
CovarSpec[["Type5"]] <- PXSloader@covars_list


# Create Objects to use for Mediation Analysis
StoreMD <- function(UKBprotlist, ScoresCombo, covarType){
  MDstore <- PXSmediation(
    Protlist = UKBprotlist,
    protIDs = names(UKBprotlist)
  )
  
  #Find the count of cases for the proteomics:
  UKBsubset <- ScoresCombo$eid
  DZage_subset <- DZ_age[DZ_age$eid %in% UKBsubset,]
  
  #'*Able to Catalog amount of Cases out of approx 23000 proteins w/ G, E, and GxE measurements*
  DZcount <- lapply(DZid_age, function(i){
    info <- sum(DZage_subset[[i]] - DZage_subset$recode_age_of_assessment_0_0 > 0,
                na.rm = T)
    return(info)
  })
  names(DZcount) <- DZid_age
  
  count <- as.data.frame(do.call("rbind", DZcount))
  count$DZid <- rownames(count)
  colnames(count) <- c("count","DZid")
  
  #'*Determine optimal number of cases to have in cox regression using power analysis later*
  #'Use a heuristic of > 100 cases for now
  DZid_subset <- count[count$count >= 100,]$DZid
  
  MDstore@DZ_df <- DZ_age %>% select(all_of(c("eid",DZid_subset,
                                              "recode_age_of_assessment_0_0",
                                              "recode_age_of_death_0_0",
                                              "age_of_removal_0_0",
                                              "age_of_lastfollowup")))
  MDstore@DZ_df <- MDstore@DZ_df %>%
    filter(eid %in% UKBsubset)
  MDstore@DZ_ids <- DZid_subset
  
  #Use same covariates as PXSloader
  MDstore@covars_df <- PXSloader@covars_df %>%
    filter(eid %in% UKBsubset)
  MDstore@covars_list <- PXSloader@covars_list
  
  #Subset covariates based on CovarType choosen:
  MDstore <- PXScovarSpec(MDstore, CovarSpec[[covarType]])
  return(MDstore)
}

MD_Type1 <- StoreMD(UKBprotlist_Type1, ScoresCombo_Type1, "Type1")
MD_Type2 <- StoreMD(UKBprotlist_Type2, ScoresCombo_Type2, "Type2")
MD_Type3 <- StoreMD(UKBprotlist_Type3, ScoresCombo_Type3, "Type3")
MD_Type4 <- StoreMD(UKBprotlist_Type4, ScoresCombo_Type4, "Type4")
MD_Type5 <- StoreMD(UKBprotlist_Type5, ScoresCombo_Type5, "Type5")



saveRDS(MD_Type1, file = "/n/scratch/users/s/shi872/UKB_intermediate/UKB_MDstore_Type1.rds")
saveRDS(MD_Type2, file = "/n/scratch/users/s/shi872/UKB_intermediate/UKB_MDstore_Type2.rds")
saveRDS(MD_Type3, file = "/n/scratch/users/s/shi872/UKB_intermediate/UKB_MDstore_Type3.rds")
saveRDS(MD_Type4, file = "/n/scratch/users/s/shi872/UKB_intermediate/UKB_MDstore_Type4.rds")
saveRDS(MD_Type5, file = "/n/scratch/users/s/shi872/UKB_intermediate/UKB_MDstore_Type5.rds")


