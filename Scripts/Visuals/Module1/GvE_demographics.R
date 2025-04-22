#Load Libraries and check for dependencies:
library(data.table)
library(tidyverse)
library(ggplot2)
library(ggpmisc)
library(ggrepel)
library(gridExtra)
library(patchwork)

#Specific Packages for this Script:
library(ggridges)


#Load Tables from Analysis:
setwd("/n/groups/patel/shakson_ukb/UK_Biobank/")
HEAPge <- readRDS("./Output/HEAPres/HEAPge.rds")
#source("./RScripts/Pure_StatGen/Prot_ExPGS/Visualizations/analyProt_PGS_PXS_final.R")

# =============================================
# R² Distributions of Demographics/Covariates
# =============================================

#'*Demographic/Covariate R2 Plots*
# Use the most adjusted model R2s
plot_R2test <- HEAPge@HEAP$Model5$R2test

# R2 of each feature:
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


# Select columns for plotting:
plot_R2test_ht <- plot_R2test_demo %>%
                        select(c("ID","G","E","GxE","fasting_time_f74_0_0","body_mass_index_bmi_f23104_0_0",
                                 "age_when_attended_assessment_centre_f21003_0_0","sex_f31_0_0",
                                 "uk_biobank_assessment_centre_f54_0_0","sum_gPCs",
                                 "medications",
                                 "age2","age_sex","age2_sex")) %>% ##dplyr::select
                        mutate(totalR2 = rowSums(across(-ID)))


#Goal: Show the top hits for high R2 proteins by topic:
plot_demoR2 <- function(df, label_id, label_name){
  gg1 <- df %>%
            mutate(label = ifelse(rank(-.data[[label_id]]) <= 20, ID, "")) %>%
            ggplot(aes(x = G, y = .data[[label_id]], label = label)) +
            geom_point() +
            geom_text_repel() +
            xlab(expression("G: PGS "~R^2)) +
            ylab(bquote(.(label_name) ~ R^2))
  
  return(gg1) 
}

d1 <- plot_demoR2(plot_R2test_ht, label_id = "age_when_attended_assessment_centre_f21003_0_0", label_name = "Age")
d2 <- plot_demoR2(plot_R2test_ht, label_id = "sex_f31_0_0", label_name = "Sex")
d3 <- plot_demoR2(plot_R2test_ht, label_id = "body_mass_index_bmi_f23104_0_0", label_name = "BMI")
d4 <- plot_demoR2(plot_R2test_ht, label_id = "fasting_time_f74_0_0", label_name = "Fasting Time")
d5 <- plot_demoR2(plot_R2test_ht, label_id = "uk_biobank_assessment_centre_f54_0_0", label_name = "UKB Assessment Center")
d6 <- plot_demoR2(plot_R2test_ht, label_id = "sum_gPCs", label_name = "Genetic PCs")

d7 <- plot_demoR2(plot_R2test_ht, label_id = "medications", label_name = "Medications")


demo1 <- (d1 | d2 | d3 ) / (d4 | d5 | d6)
demo2 <- (d1 | d2 | d3) / (d4 | d5 | d7)

lapply(c("png", "svg", "pdf"), 
       function(fmt) ggsave(paste0("./Figures/HEAP/SF2/", "DemographicR2_v1.", fmt), 
                            plot = demo1, width = 10, height = 8, dpi = 1000))

lapply(c("png", "svg", "pdf"), 
       function(fmt) ggsave(paste0("./Figures/HEAP/SF2/", "DemographicR2_v2.", fmt), 
                            plot = demo2, width = 10, height = 8, dpi = 1000))

# =============================================
# Ridgeplots of R² Distributions
# =============================================

#'*R2 Distribution Across All Features*
feat <- c(Genetics = "G", Exposures = "E", GxE = "GxE",
          `Fasting Time` = "fasting_time_f74_0_0", BMI = "body_mass_index_bmi_f23104_0_0",
          Age = "age_when_attended_assessment_centre_f21003_0_0",
          Sex = "sex_f31_0_0",
          Center = "uk_biobank_assessment_centre_f54_0_0", `Genetic PCs` = "sum_gPCs",
          Medications = "medications", Age2 = "age2", AgexSex = "age_sex",
          Age2xSex = "age2_sex", `Total R2` = "totalR2")

ridgeplotR2 <- plot_R2test_ht %>% 
                    rename(all_of(feat)) %>% 
                    pivot_longer(cols = c("Genetics","Exposures","GxE","Fasting Time","BMI",
                                          "Age","Sex",
                                          "Center","Genetic PCs",
                                          "Medications",
                                          "Age2","AgexSex","Age2xSex", "Total R2"))

highlight_categories <- c("Total R2","Genetics", "Exposures", "GxE")
all_categories = c("Genetics","Exposures","GxE","Fasting Time","BMI",
                   "Age","Sex",
                   "Center","Genetic PCs",
                   "Medications",
                   "Total R2")

# Version 1: Full with Total R2 Distribution
rp1 <- ridgeplotR2 %>%
            filter(!(name %in% c("Age2xSex", "Age2", "AgexSex"))) %>%
            ggplot(aes(x = value, y = fct_reorder(name, value, .fun = mean), 
                       fill = fct_reorder(name, value, .fun = mean))) +
            geom_density_ridges(
              jittered_points = TRUE,
              alpha = 0.5, scale = 0.8, point_size = 0.05,
              position = position_raincloud(width = 0.01, height = 0.15,
                                            ygap = 0.05)
            ) +
            scale_fill_manual(
              values = c(
                setNames(c("red","blue", "green", "magenta"), highlight_categories), # Highlighted categories
                setNames(rep("gray", length(setdiff(all_categories, highlight_categories))),
                         setdiff(all_categories, highlight_categories)) # Grayed-out categories
              )) +
            theme_minimal() +
            labs(x = bquote(R^2), y = NULL) +
            theme(legend.position="none")

lapply(c("png", "svg", "pdf"), 
       function(fmt) ggsave(paste0("./Figures/HEAP/F1/", "PartR2Dens_v1.", fmt), 
                            plot = rp1, width = 4, height = 4, dpi = 1000))



# Version 2: Version without Total R2
highlight_categories <- c("Genetics", "Exposures", "GxE")

rp2 <- ridgeplotR2 %>%
          filter(!(name %in% c("Age2xSex", "Age2", "AgexSex","Total R2"))) %>%
          ggplot(aes(x = value, y = fct_reorder(name, value, .fun = mean), 
                     fill = fct_reorder(name, value, .fun = mean))) +
          geom_density_ridges(
            panel_scaling = TRUE,
            jittered_points = TRUE,
            alpha = 0.7, scale = 0.9, point_size = 0.05,
            position = position_raincloud(width = 0.01, height = 0.15,
                                          ygap = 0.05)
          ) +
          scale_fill_manual(
            values = c(
              setNames(c("blue", "green", "magenta"), highlight_categories), # Highlighted categories
              setNames(rep("gray", length(setdiff(all_categories, highlight_categories))),
                       setdiff(all_categories, highlight_categories)) # Grayed-out categories
            )) +
          theme_minimal() +
          labs(x = bquote(R^2), y = NULL) +
          theme(legend.position="none")

lapply(c("png", "svg", "pdf"), 
       function(fmt) ggsave(paste0("./Figures/HEAP/F1/", "PartR2Dens_v2.", fmt), 
                            plot = rp2, width = 4, height = 4, dpi = 1000))
