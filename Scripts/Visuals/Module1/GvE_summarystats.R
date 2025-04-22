#Load Libraries
library(data.table)
library(tidyverse)
library(ggplot2)
library(ggpmisc)
library(ggrepel)
library(gridExtra)
library(patchwork)

#Load Tables from Analysis:
setwd("/n/groups/patel/shakson_ukb/UK_Biobank/")
HEAPge <- readRDS("./Output/HEAPres/HEAPge.rds")
#source("./RScripts/Pure_StatGen/Prot_ExPGS/Visualizations/analyProt_PGS_PXS_final.R")


#'*Simple Stats on R2 across Specifications*
# Summary Stats Table:
# Capture the summary output as strings
summary_G <- capture.output(summary(HEAPge@R2geCI$G_stats_mean))
summary_E <- capture.output(summary(HEAPge@R2geCI$E_stats_mean))
summary_GxE <- capture.output(summary(HEAPge@R2geCI$GxE_stats_mean))

# Combine all summaries into one text
summary_text <- c("Summary of G:\n", summary_G, "\n\n", 
                  "Summary of E:\n", summary_E, "\n\n", 
                  "Summary of GxE:\n", summary_GxE)
# Write the summaries to a text file
writeLines(summary_text, "./Figures/HEAP/F1/R2test_GvE_summary.txt")


# Pie Chart of Average R2 (G and E):
library(gridExtra)
plot_pie <- function(df, column, breaks, labels, title) {
  df <- df %>%
    mutate(binned = cut(.data[[column]], breaks = breaks, labels = FALSE)) %>%
    count(binned) %>%
    mutate(percentage = n / sum(n) * 100, 
           binned = labels[binned])
  
  ggplot(df, aes(x="", y=n, fill=binned)) +
    geom_bar(stat="identity", width=1, color="white") +
    coord_polar("y", start=0) +
    geom_text(aes(label = n), position = position_stack(vjust = 0.5), color = "black") +
    geom_col(aes(x=-1,y=0)) +
    theme_void() +
    ggtitle(title) +
    labs(fill = expression(R^2)) +
    theme(plot.title = element_text(hjust=0.5))
}
PieChart_GvEspec <- function(type = "AverageSpec"){
  # Genetics plot
  gg1 <- plot_pie(HEAPge@R2geCI, "G_stats_mean", c(0, 0.01, 0.1, 0.25, 0.5, 1), 
                  c("0-0.01", "0.01-0.1", "0.1-0.25", "0.25-0.5", "0.5-1"), "Genetics")
  
  # Environment plot
  gg2 <- plot_pie(HEAPge@R2geCI %>% filter(E_stats_mean > 0), "E_stats_mean", c(0, 0.01, 0.05, 0.2), 
                  c("0-0.01", "0.01-0.05", "0.05-0.2"), "Exposures")
  
  # Combine and save plots
  gg3 <- grid.arrange(gg1, gg2, nrow = 2)
  ggsave(paste0("./Figures/HEAP/F1/","EvsG_piechart_",type,".png"), 
         plot = gg3, width = 4, height = 6, dpi = 1000)
  ggsave(paste0("./Figures/HEAP/F1/","EvsG_piechart_",type,".svg"), 
         plot = gg3, width = 4, height = 6, dpi = 1000)
  ggsave(paste0("./Figures/HEAP/F1/","EvsG_piechart_",type,".pdf"), 
         plot = gg3, width = 4, height = 6, dpi = 1000)
}
PieChart_GvEspec()


#'*Potentially Explaining Low R2 Proteins:*
plot_R2test <- HEAPge@HEAP$Model5$R2test

# Create finalized matrix of R2s by each feature.
# Take average R2 across folds:
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
                                        
         )))) %>%
  group_by(ID) %>%
  summarise(across(everything(), mean))


# Convert wide to long format
plot_R2test_ht <- plot_R2test_demo %>%
  select(c("ID","G","E","GxE","fasting_time_f74_0_0","body_mass_index_bmi_f23104_0_0",
           "age_when_attended_assessment_centre_f21003_0_0","sex_f31_0_0",
           "uk_biobank_assessment_centre_f54_0_0","sum_gPCs",
           "medications",
           "age2","age_sex","age2_sex")) %>%
  mutate(totalR2 = rowSums(across(-ID)),
         GEtotR2 = rowSums(across(c("G","E","GxE"))))

# Summary of low R2 proteins and potential categories these proteins are involved in:
summary(plot_R2test_ht[c("totalR2", "GEtotR2")])
sum(plot_R2test_ht$GEtotR2 < 0.01)

lowR2prot <- plot_R2test_ht %>%
                filter(GEtotR2 < 0.01) %>%
                pull(ID) 

nlowR2prot <- length(lowR2prot)

covR2exp <- plot_R2test_ht %>%
                filter(ID %in% lowR2prot & 
                         (age_when_attended_assessment_centre_f21003_0_0 > 0.01 | 
                            sex_f31_0_0 > 0.01 | 
                            body_mass_index_bmi_f23104_0_0 > 0.01)) %>%
                nrow()

# Combine all summaries into one text
summary_text <- c("Number of Proteins: R2 < 0.01 Explained by Genetics or Environment:\n", 
                  nlowR2prot, "\n\n", 
                  "Of those xxx proteins are explained by Age, Sex, BMI with R2 > 0.01:\n", 
                  covR2exp)
# Write the summaries to a text file
writeLines(summary_text, "./Figures/HEAP/F1/R2test_GvE_lowR2summ.txt")


