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

#### COVARIATE SPECIFICATION VARIABILITY ANALYSIS ####
#'*Multiple covariate specification ANALYSIS*
library(dplyr)
library(purrr)

#'*Main GvE plot of HEAP manuscript (Encompasses Model Specification Variability)*
# Function to plot variability in GxE plot
CovarSpec_GvEplot <- function(plot_testR2_CI){

  #Highlight specific proteins with high E: R2
  subset_test <- plot_testR2_CI %>% 
    filter((E_stats_mean > G_stats_mean) & (E_stats_mean > 0.09))
  
  #Cutoff for y-axis:
  y_cutoff = 0.175
  
  #Without Error Bars:
  p1 <- plot_testR2_CI  %>% 
    ggplot(aes(x = G_stats_mean, y = E_stats_mean, label = ID)) + 
    geom_point() +
    geom_abline(slope=1, intercept = 0) +
    geom_ribbon(aes(ymin = 0, ymax = ifelse(G_stats_mean <= y_cutoff, G_stats_mean, y_cutoff), fill = "Genetics Driven"), alpha = 0.2) +
    geom_ribbon(aes(ymin = ifelse(G_stats_mean <= y_cutoff, G_stats_mean, y_cutoff), ymax = y_cutoff, fill = "Exposures Driven"), alpha = 0.2) +
    scale_fill_manual(values = c("Genetics Driven" = "blue", "Exposures Driven" = "green")) +
    ylim(c(0, y_cutoff)) +
    geom_text_repel(data=subset_test, aes(x = G_stats_mean,y = E_stats_mean,label= ID),
                    size = 4,
                    max.overlaps = 20,
                    min.segment.length = 0.025,
                    force = 10) +
    xlab(expression("G: PGS "~R^2)) +
    ylab(expression("E: PXS "~R^2)) +
    theme_minimal() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 14)) +
    guides(fill=guide_legend(title="Region")) +
    #ggtitle(expression("Test Set: Partitioned"~R^2)) +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "bottom") 
  
  ggsave(paste0("./Figures/HEAP/F1/","GvE_R2_CovarSpec.png"), 
         plot = p1, width = 6, height = 6, dpi = 1000)
  ggsave(paste0("./Figures/HEAP/F1/","GvE_R2_CovarSpec.svg"), 
         plot = p1, width = 6, height = 6, dpi = 1000)
  ggsave(paste0("./Figures/HEAP/F1/","GvE_R2_CovarSpec.pdf"), 
         plot = p1, width = 6, height = 6, dpi = 1000)
  
  
  #With Error Bars: All Covariate Specifications
  p2 <- plot_testR2_CI  %>% 
    ggplot(aes(x = G_stats_mean, y = E_stats_mean, label = ID)) + 
    geom_point() +
    geom_abline(slope=1, intercept = 0) +
    geom_ribbon(aes(ymin = 0, ymax = ifelse(G_stats_mean <= y_cutoff, G_stats_mean, y_cutoff), fill = "Genetics Driven"), alpha = 0.2) +
    geom_ribbon(aes(ymin = ifelse(G_stats_mean <= y_cutoff, G_stats_mean, y_cutoff), ymax = y_cutoff, fill = "Exposures Driven"), alpha = 0.2) +
    scale_fill_manual(values = c("Genetics Driven" = "blue", "Exposures Driven" = "green")) +
    ylim(c(0, y_cutoff)) +
    geom_text_repel(data=subset_test, aes(x = G_stats_mean,y = E_stats_mean,label= ID),
                    size = 3.5) +
    geom_pointrange(aes(ymin = E_stats_ci_low, ymax = E_stats_ci_high)) +
    geom_pointrange(aes(xmin = G_stats_ci_low, xmax = G_stats_ci_high)) +
    #geom_label() + 
    xlab(expression("G: PGS "~R^2)) +
    ylab(expression("E: PXS "~R^2)) +
    theme_minimal() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    guides(fill=guide_legend(title="Region")) +
    ggtitle(expression("Test Set: Partitioned"~R^2)) +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "bottom") 
  
  ggsave(paste0("./Figures/HEAP/A1/","GvE_R2_CovarSpec_wErrorBars.png"), 
         plot = p2, width = 6, height = 6, dpi = 1000)
  ggsave(paste0("./Figures/HEAP/A1/","GvE_R2_CovarSpec_wErrorBars.svg"), 
         plot = p2, width = 6, height = 6, dpi = 1000)
  ggsave(paste0("./Figures/HEAP/A1/","GvE_R2_CovarSpec_wErrorBars.pdf"), 
         plot = p2, width = 6, height = 6, dpi = 1000)
}
CovarSpec_GvEplot(HEAPge@R2geCI)


#'*Count number of proteins where E > G*
#Highlight specific proteins with high E: R2
Eproteins <- HEAPge@R2geCI %>% 
  filter(E_stats_mean > G_stats_mean)
write.csv(Eproteins, file = "./Figures/HEAP/T1/Eproteins.csv")

#Stats:
# Capture the summary output as strings
summary1 <- capture.output(summary(Eproteins$E_stats_mean))
summary2 <- capture.output(length(unique(Eproteins$ID)))

# Combine all summaries into one text
summary_text <- c("Summary of E Driven Proteins:\n", summary1, "\n\n", 
                  "Number of E Driven Proteins:\n", summary2, "\n\n")
# Write the summaries to a text file
writeLines(summary_text, "./Figures/HEAP/T1/EDrivenProteins_summary.txt")


#'*Plot variability of E, G, both separately*
#Renamed Standard Error to "Variability across Specification"
CovarSpec_varplot <- function(R2test){
  #'*Good Plot for overall E: R2 variability*
  E_varplot <- R2test %>%
    group_by(ID, cID) %>%
    summarise(
      mean_E = mean(E)
    ) %>%
    summarise(
      agg_E = mean(mean_E), 
      stdError_E = sd(mean_E)/sqrt(n()),
      CI_margin_E = qt(0.975, df = n() - 1) * stdError_E,
      n = n()
    )
  
  p1 <- E_varplot %>%
    ggplot(aes(x = agg_E, y = stdError_E, label = ID)) +
    geom_point() +
    xlab(expression("E: "~R^2)) +
    ylab("Std. Error") +
    geom_text_repel(aes(label = ifelse(agg_E > 0.07, ID, ""))) +
    theme_minimal()
  
  G_varplot <- R2test %>%
    group_by(ID, cID) %>%
    summarise(
      mean_G = mean(G)
    ) %>%
    summarise(
      agg_G = mean(mean_G), 
      stdError_G = sd(mean_G)/sqrt(n()),
      CI_margin_G = qt(0.975, df = n() - 1) * stdError_G,
      n = n()
    )
  
  p2 <- G_varplot %>%
    ggplot(aes(x = agg_G, y = stdError_G, label = ID)) +
    geom_point() +
    xlab(expression("G: "~R^2)) +
    ylab("Std. Error") +
    geom_text_repel(aes(label = ifelse(stdError_G > 5e-4, ID, ""))) +
    theme_minimal() 
  
  #'*Combine plots:*
  varplot <- merge(E_varplot, G_varplot, by = "ID")
  
  p3 <- varplot %>%
    ggplot(aes(x = stdError_G, y = stdError_E, color = agg_E)) +
    geom_point() +
    #scale_color_gradient(low = "", high = "green") +
    scale_color_viridis_c(option = "viridis") +
    #geom_text_repel(aes(label = ifelse(agg_E > 0.07, ID, ""))) +
    xlab("Std. Error: G") +
    ylab("Std. Error: E") +
    labs(color='E: R2')  +
    #ggtitle("Variability of G and E R2 across Covariate Specifications") +
    theme_minimal()
  
  #Save the plots:
  lapply(c("png", "svg", "pdf"), 
         function(fmt) ggsave(paste0("./Figures/HEAP/F1/","Evar.", fmt), 
                              plot = p1, width = 3, height = 3, dpi = 1000))
  
  lapply(c("png", "svg", "pdf"), 
         function(fmt) ggsave(paste0("./Figures/HEAP/F1/","Gvar.", fmt), 
                              plot = p2, width = 3, height = 3, dpi = 1000))
  
  lapply(c("png", "svg", "pdf"), 
         function(fmt) ggsave(paste0("./Figures/HEAP/F1/","EvsGvar.", fmt), 
                              plot = p3, width = 4, height = 3, dpi = 1000))
}
CovarSpec_varplot(HEAPge@R2ge)
