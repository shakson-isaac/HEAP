#Libraries:
library(data.table)
library(tidyverse)
library(ggplot2)
library(ggpmisc)

##Organization: DataStructure -- PXSconstruct

#'*PXS Object:*
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

#'*Load Proteomics*
library(lubridate)

#'*Loading UKB Proteomics*
UKBprot <- fread(file = "/n/groups/patel/uk_biobank/olink_22881_52887/olink_data_52887.txt", header = T)

#Protein ID Conversion
protein_id <- read.table("/n/groups/patel/uk_biobank/project_22881_672185/protein_id_conv.txt", header = T, sep = "\t",
                         fill=TRUE, quote="", encoding="UTF-8")
#Separate by ";"
protein_id[c('prot_id', 'prot_name')] <- str_split_fixed(protein_id$meaning, ';', 2)

#Convert olink dataframe to wide format:
UKBprot <- dcast.data.table(UKBprot, eid + ins_index ~ protein_id, value.var = "result")  

#Switch protCoding to protID
UKBprot <- UKBprot %>%
  rename_at(vars(3:ncol(.)), ~ protein_id$prot_id[protein_id$coding %in% .])

#'Modify Protein Names to all have underscores*
prot_names <- colnames(UKBprot)[-c(1:2)]
prot_names <- gsub("-","_", prot_names)
colnames(UKBprot)[-c(1:2)] <- prot_names

#'ins_index for proteomics: subset only for instance 0*
UKBprot <- UKBprot[UKBprot$ins_index == 0, ]

PXSloader <- PXSconstruct(
  UKBprot_df = UKBprot,
  protIDs = prot_names
)

#'*Load Environmental Variables*
setwd("/n/groups/patel/shakson_ukb/UK_Biobank/RScripts/Extract_Raw/Finalized")
source("./dataloader_functions_parallel.R")
projID = 52887
load_project(projID)
UKBdict <- fread(file = paste0("/n/groups/patel/shakson_ukb/UK_Biobank/Data/Paths/",projID,"/allpaths.txt")) #'*UKBdict should be a global variable now!!*


#SES, Income, Deprivation
globalE = c(22189, 26411, 26428, 26418, 26410, 26427, 26426)
globalE_df <- fast_dataloader_viafield(UKBdict, UKBfieldIDs = globalE, directoryInfo)
globalE_df <- UKB_instances(globalE_df, "_0_")


#Lifestyle Factors:
# Load the Function and Libraries needed for UKB DataLoader
setwd("/n/groups/patel/shakson_ukb/UK_Biobank/RScripts/Extract_Raw/Finalized")
source("./dataloader_functions_upd.R")

#Load Relevant Info For a Given Project:
projID = 52887
load_project(projID)

#'*Lifestyle Dataframe*:
lifestyle_pathIDs <- pathIDs %>% filter(grepl("Lifestyle", Path)) %>% pull(Category)
lifestyle_pathIDs <- paste0("path_",lifestyle_pathIDs) #this got updates:

plan(multisession, workers = 1)

#Load processed Dataframes:
read_parquet_file <- function(path_id, features) {
  fileout <- paste0("/n/groups/patel/shakson_ukb/UK_Biobank/Data/UKB/",projID,"/",path_id,".parquet")
  arrow::read_parquet(fileout)
}

# Use furrr to parallelize the loop
df <- future_map(lifestyle_pathIDs, ~read_parquet_file(.x))
df <- setNames(df, lifestyle_pathIDs)
lifestyle_df <- df

#Vitamins and minerals: 6155, 6179
vit <- c(6155, 6179) 
vit_df <- fast_dataloader_viafield(UKBdict, UKBfieldIDs = vit, directoryInfo)
vit_df <- UKB_instances(vit_df, "_0_")
vit_df <- UKB_multiarray_handle(vit_df)
vit_df <- UKB_onehot_handle(vit_df)


#'*Organize Environmental Variables: Elist, Elist_names and Eids*
#Function to organize Category Names
orgElistnames <- function(Elist_names){
  #Make the naming structure consistent:
  Elist_names <- gsub(" ","_", Elist_names)
  Elist_names <- gsub(paste(c("[(]", "[)]"), collapse = "|"),"",Elist_names)
  
  return(Elist_names)
}

#Function to organize Eids Names
orgEids <- function(Elist){
  for(i in names(Elist)){
    Eids <- gsub(" ","_", colnames(Elist[[i]]))
    Eids <- make.names(Eids)
    colnames(Elist[[i]]) <- Eids
  }
  return(Elist)
}

#Function to link Category Names to Eids
orgEtoCat <- function(Elist, Elist_names){
  Ecategory <- lapply(Elist, function(x) setdiff(colnames(x), "eid"))
  #total_entries <- sum(sapply(Ecategory, length))
  
  names(Ecategory) <- Elist_names
  
  Eid_cat <- stack(Ecategory)
  colnames(Eid_cat) <- c("Eid", "Category")
  Eid_cat$Category <- as.character(Eid_cat$Category)
  
  return(Eid_cat)
}

#Function to quickly identify ordinal features:
ordinal_finder <- function(df){
  max_cols <- apply(df, 2, max, na.rm = T)
  ordinals <- names(max_cols[max_cols <= 5 & max_cols > 1])
  return(ordinals)
}


#Load base dataframes with E
PXSloader@Elist <- list(lifestyle_df$path_100051,lifestyle_df$path_100052,
                        lifestyle_df$path_100058,lifestyle_df$path_54,
                        lifestyle_df$path_100054,lifestyle_df$path_100053,
                        globalE_df,vit_df)
PXSloader@Elist_names <- c("Alcohol","Diet (Weekly)",
                           "Smoking","Exercise (MET)",
                           "Exercise (Freq)","Internet Usage",
                           "Deprivation Indices", "Vitamins")
PXSloader@Elist_names <- orgElistnames(PXSloader@Elist_names)
names(PXSloader@Elist) <- PXSloader@Elist_names

#Specify format of Environmental IDs and link to specific categories.
PXSloader@Elist <- orgEids(PXSloader@Elist)
PXSloader@Eid_cat <- orgEtoCat(PXSloader@Elist, PXSloader@Elist_names)


#Specify Ordinal Features
PXSloader@ordinalIDs <- unlist(lapply(PXSloader@Elist, function(x) ordinal_finder(x)))
PXSloader@ordinalIDs <- unname(PXSloader@ordinalIDs)


#'*Load Covariates*
###COVARIATES:

#BASE SET:
#Age: 21003
#SEX: 31
#BMI: 23104
#Fasting Time: 71

#Extended SET:
#Assessment Center: 54
#Genetic PCs: 22009
#ADD interaction and nonlinear terms: age^2, age * sex, age^2 * sex.

#Extended SET 2:
#Common Medications: 6177 & 6153

covars = c(21003, 31 , 23104, 74, 54, 22009)
covars_df <- fast_dataloader_viafield(UKBdict, UKBfieldIDs = covars, directoryInfo)
covars_df <- UKB_instances(covars_df, "_0_")

#AGE and SEX interaction TERMS
covars_df$age2 <- (covars_df$age_when_attended_assessment_centre_f21003_0_0)^2
covars_df$age_sex <- covars_df$age_when_attended_assessment_centre_f21003_0_0 * ifelse(covars_df$sex_f31_0_0 == "Male", 1, 0)
covars_df$age2_sex <- (covars_df$age_when_attended_assessment_centre_f21003_0_0)^2 * ifelse(covars_df$sex_f31_0_0 == "Male", 1, 0)

#Double check certain categorical variables encoded:
covars_df$uk_biobank_assessment_centre_f54_0_0 <- as.factor(covars_df$uk_biobank_assessment_centre_f54_0_0)


# COMMON Medications: 6177, 6153
medi <- c(6177, 6153)
medi_df <- fast_dataloader_viafield(UKBdict, UKBfieldIDs = medi, directoryInfo)
medi_df <- UKB_instances(medi_df, "_0_")
medi_df <- UKB_multiarray_handle(medi_df)

#Need to combine both columns since they include repeated responses in different questions:
#'*TO-DO change "combined" naming to somethign relevant*
medi_df_combo <- medi_df %>%
  pivot_longer(cols = -eid, names_to = "variable", values_to = "value") %>%
  group_by(eid) %>%
  summarise(combined = paste(unique(unlist(strsplit(value[!is.na(value)], ";"))), collapse = ";")) %>%
  ungroup()
medi_df <- UKB_onehot_handle(medi_df_combo)
colnames(medi_df) <- gsub(" ","_", colnames(medi_df))


covars_df <- merge(covars_df, medi_df, by = "eid")

#Define covars_list to control level of adjustment in the PXS construction later.
full_covars_list <- colnames(covars_df)[-c(1)]

PXSloader@covars_df <- covars_df 
PXSloader@covars_list <- full_covars_list



#Save RDS file of PXSloader object
gc()
class(PXSloader)
saveRDS(PXSloader, file = "/n/scratch/users/s/shi872/UKB_intermediate/UKB_PGS_PXS_load.rds")


