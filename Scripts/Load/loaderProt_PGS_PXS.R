#Libraries:
library(data.table)
library(tidyverse)
library(ggplot2)
library(ggpmisc)

#'*Loading:*
#Load MetaboProteome Data:
UKB_Omics <- readRDS("/n/scratch/users/s/shi872/UKB_intermediate/UKB_fc_metaboproteome.rds")
gc()


#Extract Proteomics:
UKBprot_df <- UKB_Omics@data %>% select(all_of(c("eid",UKB_Omics@protIDs)))
#colSums(!is.na(UKBprot_df))


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




#Check nutrient info: 24-hour survey
nutr_metab <- readRDS(file = "/n/scratch/users/s/shi872/UKB_intermediate/UKB_nutr_metab_df.rds")
#'*Need to update this later:*
nutr_list = scan(file = "/n/groups/patel/shakson_ukb/UK_Biobank/RScripts/Env_Preprocess/nutr_names.txt",
                 what = as.character())
UKBnutr_df <- nutr_metab %>% select(all_of(c("eid",nutr_list)))

#70K individuals at instance 0: according to UKBshowcase - https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=26030
sum(!is.na(UKBnutr_df$total_weight_of_all_foods_and_beverages_f26000_0_0))




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


# Medications: 6177, 6153
medi <- c(6177, 6153)
medi_df <- fast_dataloader_viafield(UKBdict, UKBfieldIDs = medi, directoryInfo)
medi_df <- UKB_instances(medi_df, "_0_")
medi_df <- UKB_multiarray_handle(medi_df)

#Old Version:
#medi_df2 <- UKB_onehot_handle(medi_df)

#Need to combine both columns since they include repeated responses in different questions:
#'*TO-DO change "combined" naming to somethign relevant*
medi_df_combo <- medi_df %>%
  pivot_longer(cols = -eid, names_to = "variable", values_to = "value") %>%
  group_by(eid) %>%
  summarise(combined = paste(unique(unlist(strsplit(value[!is.na(value)], ";"))), collapse = ";")) %>%
  ungroup()
medi_df <- UKB_onehot_handle(medi_df_combo)


# Medications: 6154 (pain relief) (TO-ADD Later)


#Vitamins and minerals: 6155, 6179
vit <- c(6155, 6179) 
vit_df <- fast_dataloader_viafield(UKBdict, UKBfieldIDs = vit, directoryInfo)
vit_df <- UKB_instances(vit_df, "_0_")
vit_df <- UKB_multiarray_handle(vit_df)
vit_df <- UKB_onehot_handle(vit_df)

# takes 7-10 minutes to load!!!

gc()
#Save this RDS file even if it is not super organized yet:
save.image(file = "/n/scratch/users/s/shi872/UKB_intermediate/UKB_PGS_PXS_data.rds")
##'*Change to .Rdata extension next time!!!*






#Additional Covariates:
#batch/plate: 30901 - ? need to also use this file: https://biobank.ndph.ox.ac.uk/showcase/ukb/auxdata/olink_batch_number.dat
#if want the correct batch number.
#well used: 30902 - ?
#consortium participants: 30903 - ?
#assessment centre: 54 - main confounder
#genetic PCs: 22009 - another confounder {populations tratification}
#ADD interaction and nonlinear terms: age^2, age * sex, age^2 * sex.

adjust = c(54, 22009)
adjust_df <- fast_dataloader_viafield(UKBdict, UKBfieldIDs = adjust, directoryInfo)
adjust_df <- UKB_instances(adjust_df, "_0_")

#Due to Assessment Centre Categorical Variable:
#adjust_df$plate_used_for_sample_run_f30901_0_0 <- as.factor(adjust_df$plate_used_for_sample_run_f30901_0_0)
#adjust_df$well_used_for_sample_run_f30902_0_0 <- as.factor(adjust_df$well_used_for_sample_run_f30902_0_0)
adjust_df$uk_biobank_assessment_centre_f54_0_0 <- as.factor(adjust_df$uk_biobank_assessment_centre_f54_0_0)

gPC_names <- paste0("genetic_principal_components_f22009_0_",1:20)
#"plate_used_for_sample_run_f30901_0_0","well_used_for_sample_run_f30902_0_0"
adjust_list <- c("uk_biobank_assessment_centre_f54_0_0", gPC_names)

covars_df <- merge(covars_df, adjust_df, by = "eid")
covars_list <- c(covars_list, adjust_list)

#Check # of samples:
sum(table(adjust_df$plate_used_for_sample_run_f30901_0_0)) > 50000
sum(table(adjust_df$well_used_for_sample_run_f30902_0_0)) > 50000

#not a good variable:
sum(table(adjust_df$ukb_ppp_consortium_selected_participant_f30903_0_0)) > 50000


#Additional Interactions:
covars_df$age2 <- (covars_df$age_when_attended_assessment_centre_f21003_0_0)^2
covars_df$age_sex <- covars_df$age_when_attended_assessment_centre_f21003_0_0 * ifelse(covars_df$sex_f31_0_0 == "Male", 1, 0)
covars_df$age2_sex <- (covars_df$age_when_attended_assessment_centre_f21003_0_0)^2 * ifelse(covars_df$sex_f31_0_0 == "Male", 1, 0)
covars_list <- c(covars_list, "age2", "age_sex","age2_sex")
colnames(covars_df)




gc()
#Save this RDS file even if it is not super organized yet:
save.image(file = "/n/scratch/users/s/shi872/UKB_intermediate/UKB_PGS_PXS_data_addncov.rds")




