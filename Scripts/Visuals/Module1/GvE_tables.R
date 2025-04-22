library(data.table)
library(tidyverse)
library(ggplot2)
library(ggpmisc)
library(ggrepel)
library(gridExtra)
library(plotly)
library(htmlwidgets)

#Addn Libraries for Visualization:
library(patchwork)

#Output tables for specific components:
setwd("/n/groups/patel/shakson_ukb/UK_Biobank/")

#Load Summary Stats:
HEAPge <- readRDS("./Output/HEAPres/HEAPge.rds")


#'*Table: E, G, GvE R2 w/ mean and CIs*
GxE_R2table <- HEAPge@R2geCI %>%
  mutate(across(where(is.numeric), ~ signif(., digits = 2))) %>%
  mutate(G = G_stats_mean,
         G_CI = paste0("(",G_stats_ci_low,",",G_stats_ci_high,")"),
         E = E_stats_mean,
         E_CI = paste0("(",E_stats_ci_low,",", E_stats_ci_high,")"),
         GxE = GxE_stats_mean,
         GxE_CI = paste0("(",GxE_stats_ci_low,",", GxE_stats_ci_high,")")) %>%
  select(c("ID","G","G_CI","E","E_CI","GxE","GxE_CI"))
colnames(GxE_R2table) <- c("Protein","Genetics (G)","G: 95% CI",
                           "Exposures (E)","E: 95% CI",
                           "Interactions (GxE)",
                           "GxE: 95% CI")

# Public Facing App: CSV, Parquet
fwrite(GxE_R2table, file = "/n/groups/patel/shakson_ukb/UK_Biobank/Output/App/Tables/GxE_R2table.csv")
check <- fread(file = "/n/groups/patel/shakson_ukb/UK_Biobank/Output/App/Tables/GxE_R2table.csv")


#'*Table R2 of E categories: *
HEAPres <- HEAPge@HEAP

categories <- c("Alcohol","Diet_Weekly","Smoking","Exercise_MET",
                "Exercise_Freq","Internet_Usage","Deprivation_Indices","Vitamins")

# purrr::map2 iterate over the list created to generate one dataframe:
R2testCat <- map2_dfr(HEAPres, names(HEAPres), ~ .x$R2testCat %>%
                        select(all_of(c("ID", categories))) %>%
                        mutate(cID = .y))

R2testCat_CI <- R2testCat %>%
  group_by(ID, cID) %>%
  summarise(across(c(categories), mean, .names = "{col}")) %>%  #Mean of each Specification
  mutate(Exercise = Exercise_Freq + Exercise_MET) %>%
  mutate(across(c(categories,"Exercise"), mean, .names = "mean_{col}")) %>% #Mean across Specification
  mutate(across(c(categories,"Exercise"), 
                list(
                  stdErr = ~ sd(.) / sqrt(n()) #, 
                  #ci_low = ~ . - qt(0.975, df = n() - 1) * stdErr, 
                  #ci_high = ~ . + qt(0.975, df = n() - 1) * stdErr
                ), 
                .names = "{col}_{fn}")
  ) %>% #Obtain standard errors
  mutate(
    n = n()
  ) %>%
  mutate(across(starts_with("mean_"),
                list(
                  ci_low = ~ . - qt(0.975, df = n - 1) * get(paste0(sub("mean_", "", cur_column()), "_stdErr")), 
                  ci_high = ~ . + qt(0.975, df = n - 1) * get(paste0(sub("mean_", "", cur_column()), "_stdErr"))
                ),
                .names = "{col}_{fn}")
  ) #obtain CI margins for plotting


plot_testR2Cat_CI <- R2testCat_CI %>%
  select(-c("cID","Alcohol","Diet_Weekly","Smoking","Exercise_MET",
            "Exercise_Freq","Exercise","Internet_Usage",
            "Deprivation_Indices","Vitamins")) %>%
  unique()


GxE_Cat_R2table <- plot_testR2Cat_CI %>%
  mutate(across(where(is.numeric), ~ signif(., digits = 2))) %>%
  mutate(Alcohol = mean_Alcohol,
         Alcohol_CI = paste0("(",mean_Alcohol_ci_low,",",mean_Alcohol_ci_high,")"),
         Diet = mean_Diet_Weekly,
         Diet_CI = paste0("(",mean_Diet_Weekly_ci_low,",",mean_Diet_Weekly_ci_high,")"),
         Smoking = mean_Smoking,
         Smoking_CI = paste0("(",mean_Smoking_ci_low,",",mean_Smoking_ci_high,")"),
         Exercise = mean_Exercise,
         Exercise_CI = paste0("(",mean_Exercise_ci_low,",",mean_Exercise_ci_high,")"),
         InternetUse = mean_Internet_Usage,
         InternetUse_CI = paste0("(",mean_Internet_Usage_ci_low,",",mean_Internet_Usage_ci_high,")"),
         DeprivationIndices = mean_Deprivation_Indices,
         DeprivationIndices_CI = paste0("(",mean_Deprivation_Indices_ci_low,",",mean_Deprivation_Indices_ci_high,")"),
         Vitamins = mean_Vitamins,
         Vitamins_CI = paste0("(",mean_Vitamins_ci_low,",",mean_Vitamins_ci_high,")"),
  ) %>%
  select(ID,Alcohol,Alcohol_CI,
         Diet,Diet_CI,
         Smoking,Smoking_CI,
         Exercise,Exercise_CI,
         InternetUse,InternetUse_CI,
         DeprivationIndices,DeprivationIndices_CI,
         Vitamins,Vitamins_CI)

fwrite(GxE_Cat_R2table, file = "/n/groups/patel/shakson_ukb/UK_Biobank/Output/App/Tables/GxE_Cat_R2table.csv")
check2 <- fread(file = "/n/groups/patel/shakson_ukb/UK_Biobank/Output/App/Tables/GxE_Cat_R2table.csv")

#'*Table: Covariate R2 w/ last specification*
# Use the most adjusted model R2s
plot_R2test <- HEAPge@HEAP$Model5$R2test

plot_R2test_demo <- plot_R2test %>%
  mutate(sum_gPCs = rowSums(across(paste0("genetic_principal_components_f22009_0_",1:40))),
         medications = rowSums(across(c("combined_Blood_pressure_medication",
                                        "combined_Hormone_replacement_therapy",         
                                        "combined_Oral_contraceptive_pill_or_minipill",
                                        "combined_Insulin",
                                        "combined_Cholesterol_lowering_medication",
                                        "combined_Do_not_know",
                                        "combined_None_of_the_above",
                                        "combined_Prefer_not_to_answer" 
                                        
         ))))  %>%
  group_by(ID) %>%
  summarise(across(everything(), mean)) 

#%>% ##dplyr::select
#  mutate(totalR2 = rowSums(across(-c(ID, medications, sum_gPCs))))


# Select solumns for plotting:
CovR2_table <- plot_R2test_demo %>%
  select(c("ID","G","E","GxE","fasting_time_f74_0_0","body_mass_index_bmi_f23104_0_0",
           "age_when_attended_assessment_centre_f21003_0_0","sex_f31_0_0",
           "uk_biobank_assessment_centre_f54_0_0","sum_gPCs",
           "medications",
           "age2","age_sex","age2_sex")) %>% ##dplyr::select
  mutate(totalR2 = rowSums(across(-ID)))

colnames(CovR2_table) <- c("Protein","Genetics (G)","Exposures (E)",
                              "Interactions (GxE)", "Fasting Time",
                              "Body Mass Index", "Age", "Sex", "UKB Assessment Centre",
                              "Genetic PCs","Medications","Age Squared",
                              "Age-Sex Interaction", "Age Squared-Sex Interaction",
                              "Total R2")

fwrite(CovR2_table, file = "/n/groups/patel/shakson_ukb/UK_Biobank/Output/App/Tables/GxE_Covariates_R2table.csv")
check3 <- fread(file = "/n/groups/patel/shakson_ukb/UK_Biobank/Output/App/Tables/GxE_Covariates_R2table.csv")

