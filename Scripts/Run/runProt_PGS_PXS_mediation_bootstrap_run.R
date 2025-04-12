# Script to Run Bootstrapping: Confidence Intervals for Mediation Analysis
#'*CHECK Bootstrapping Methods~~*
source("/n/groups/patel/shakson_ukb/UK_Biobank/RScripts/Pure_StatGen/Prot_ExPGS/runProt_PGS_PXS_mediation_bootstrap_source.R")

#'*Mediation Parallelization*
args = commandArgs(trailingOnly = T)
idx <- as.integer(args[1])
nSamp <- as.character(args[2])
cType <- as.character(args[3])


MDpairs <- fread(file = "/n/groups/patel/shakson_ukb/UK_Biobank/BScripts/ProtMediation/pairs.txt")
omicID <- as.character(MDpairs[idx,"omic"])
dzID <- as.character(MDpairs[idx,"dz"])



#Run through all protein/disease pairs:
#How long does 1000 bootstraps take for 1-protein, 1-disease pair: ~1 hour
MDbootdf <-runBoot_MDwrap(covarType = cType, 
                           protlist = omicID, 
                           dzlist = dzID, 
                           nSamp)

clean_results <- function(df){
  cols <- c("Exposure Indirect Effect", "Exposure Indirect Effect q2.5", "Exposure Indirect Effect q97.5",
            "Genetic Indirect Effect", "Genetic Indirect Effect q2.5", "Genetic Indirect Effect q97.5",
            "Exposure Direct Effect", "Exposure Direct Effect q2.5", "Exposure Direct Effect q97.5" ,
            "Genetic Direct Effect", "Genetic Direct Effect q2.5", "Genetic Direct Effect q97.5",
            "Total Indirect Effect", "Total Indirect Effect q2.5", "Total Indirect Effect q97.5",
            "Total Direct Effect", "Total Direct Effect q2.5", "Total Direct Effect q97.5",
            "Total Effect", "Total Effect q2.5", "Total Effect q97.5")
  new_names <- c("E_HRi_orig", "E_HRi_2.5", "E_HRi_97.5",
                 "G_HRi_orig", "G_HRi_2.5", "G_HRi_97.5",
                 "E_HRd_orig", "E_HRd_2.5", "E_HRd_97.5",
                 "G_HRd_orig", "G_HRd_2.5", "G_HRd_97.5",
                 "T_HRi_orig", "T_HRi_2.5", "T_HRi_97.5",
                 "T_HRd_orig", "T_HRd_2.5", "T_HRd_97.5",
                 "T_HR_orig", "T_HR_2.5", "T_HR_97.5")
  df <- df %>%
    mutate(across(
      all_of(cols),
      .fns = ~exp(.),
      .names = "tmp_{.col}"
    )) %>%
    rename_with(~new_names, starts_with("tmp_"))
  
  df <- df %>%
    select(c("ID","DZid",
             "E_HRi_orig", "E_HRi_2.5", "E_HRi_97.5",
             "G_HRi_orig", "G_HRi_2.5", "G_HRi_97.5",
             "E_HRd_orig", "E_HRd_2.5", "E_HRd_97.5",
             "G_HRd_orig", "G_HRd_2.5", "G_HRd_97.5",
             "T_HRi_orig", "T_HRi_2.5", "T_HRi_97.5",
             "T_HRd_orig", "T_HRd_2.5", "T_HRd_97.5",
             "T_HR_orig", "T_HR_2.5", "T_HR_97.5",
             "HR","HR q2.5","HR q97.5",
             "cindex","cindex q2.5","cindex q97.5",
             "R2","R2 q2.5","R2 q97.5"))
  return(df)
}
MDfinaldf <- clean_results(MDbootdf)
write.table(MDfinaldf,
            file = paste0("/n/groups/patel/shakson_ukb/UK_Biobank/Data/Parallel/Prot_PGSPXS_mediationboot/",
                          cType,"/","MDboot_",idx,".txt"),
            row.names = F
            )


##### Did below to populate txt files for protein and disease code IDs: ####
#Load Mediation Data:
# MDloader <- readRDS(file = paste0("/n/scratch/users/s/shi872/UKB_intermediate/UKB_MDstore_",
#                                   covarType,".rds"))

#Get omiclist to parallelize through
# omiclist <- MDloader@protIDs
# dzlist <- MDloader@DZ_ids
# 
# write.table(omiclist, "/n/groups/patel/shakson_ukb/UK_Biobank/BScripts/ProtMediation/protlist.txt",
#             row.names = F, col.names = F)
# write.table(dzlist, "/n/groups/patel/shakson_ukb/UK_Biobank/BScripts/ProtMediation/dzlist.txt",
#             row.names = F, col.names = F)

#Obtain a paired list to iterate: (Small EXAMPLE But can extend to all)
# omiclist_final <- c("MAMDC4","CEACAM16", "CXCL17","LAMP3","WFDC2","TNR","GDF15",
#                     "IGFBP1","IGFBP2","CKB","LPL","FABP4","FGF21","LEP","IL18R1","IGSF9")
# 
# dzlist_final <- c("age_f10_first_reported_mental_and_behavioural_disorders_due_to_use_of_alcohol_f130854_0_0",
#                   "age_j43_first_reported_emphysema_f131490_0_0", 
#                   "age_e11_first_reported_non_insulin_dependent_diabetes_mellitus_f130708_0_0")
# 
# # Generate all combinations
# pair_df <- expand.grid(omic = omiclist_final, dz = dzlist_final, stringsAsFactors = FALSE)
# write.table(pair_df, "/n/groups/patel/shakson_ukb/UK_Biobank/BScripts/ProtMediation/pairs.txt",
#             row.names = F, col.names = T)
