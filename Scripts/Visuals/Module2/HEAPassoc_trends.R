#Libraries:
library(data.table)
library(tidyverse)
library(ggplot2)
library(ggpmisc)
library(pbapply)
library(dplyr)
library(purrr)
library(plotly)
library(htmlwidgets)
library(qs) #For faster save/load of rds object.

#'*LOAD CovarSpecList*
setwd("/n/groups/patel/shakson_ukb/UK_Biobank/")
HEAPassoc <- qread("./Output/HEAPres/HEAPassoc.qs")

# FUNCTION: Replication Plot - Scatterplot
plot_replication <- function(assoc_df, title, fname){
  #Example Plot: Are Estimates similar between train & test?
  gg1 <- assoc_df %>%
    ggplot(aes(x = Estimate_train, y = Estimate_test)) +
    geom_point() +
    geom_smooth(method = "lm") +
    theme_minimal() +
    ggtitle(title)
  
  ggsave(paste0("./Figures/HEAP/SF3/",fname,".png"), 
         plot = gg1, width = 3, height = 3, dpi = 1000)
}

#Correlation of Associations in Train/Test (Only Replicated Ones)
RunRepPlots<- function(CovarSpecAssocList){
  lapply(names(CovarSpecAssocList), function(x){
   
     plot_replication(assoc_df = CovarSpecAssocList[[x]]$E$sigBOTH, 
                     title = paste0(x,": E Assoc Replication"),
                     fname = paste0(x,"_ReplicationCor"))
  })
}
RunRepPlots(HEAPassoc@HEAPsig)


#'*Understand Similarities and Differences across Covariate Specification*
# Function to Compare Association Concordance Across Specifications:
AssocCor_DF <- function(CovarSpecAssocList){
  repSTAT <- list()
  
  CovarSpecAssocCor <- lapply(names(CovarSpecAssocList), function(x){
    CovarSpecAssocList[[x]][["E"]][["sigBOTH"]]$CovarType <- x
    CovarSpecAssocList[[x]][["GxE"]][["sigBOTH"]]$CovarType <- x
    return(CovarSpecAssocList[[x]])
  })
  names(CovarSpecAssocCor) <- names(CovarSpecAssocList)
  
  # Aggregate Info across Specifications:
  Eassoc_AllType <- lapply(names(CovarSpecAssocCor), function(x){
    return(CovarSpecAssocCor[[x]][["E"]][["sigBOTH"]])
  }) %>% reduce(full_join)
  
  GxEassoc_AllType <- lapply(names(CovarSpecAssocCor), function(x){
    return(CovarSpecAssocCor[[x]][["GxE"]][["sigBOTH"]])
  }) %>% reduce(full_join)
  
  #Filter Associations Detected Across Specifications:
  
  #FIND THE margins for CIs
  # 95% confidence interval: use 1.96 as critical value since N > 100 for each association.
  Eassoc_AllType$CImargin_test <- Eassoc_AllType$Std.Error_test * 1.96
  GxEassoc_AllType$CImargin_test <- GxEassoc_AllType$Std.Error_test * 1.96
  #qt(0.975, df = samplesize - num_covars - 1)
  
  
  ###Dataframe for Association Tests:
  Eassoc_AllType_CI <- Eassoc_AllType %>%
    select(AssocID, CovarType, Estimate_test, Std.Error_test, CImargin_test, Category_test) %>% 
    pivot_wider(names_from = CovarType, values_from = c(Estimate_test, Std.Error_test, CImargin_test))
  
  
  GxEassoc_AllType_CI <- GxEassoc_AllType %>%
    select(AssocID, CovarType, Estimate_test, Std.Error_test, CImargin_test, Category_test) %>% 
    pivot_wider(names_from = CovarType, values_from = c(Estimate_test, Std.Error_test, CImargin_test))
  
  repSTAT[["E"]] <- Eassoc_AllType_CI
  repSTAT[["GxE"]] <- GxEassoc_AllType_CI
  
  return(repSTAT)
}

corDF <- AssocCor_DF(HEAPassoc@HEAPsig)

# Get STATS:
AllEassoc_cor <- corDF[["E"]]  %>%
  select(all_of(c("Estimate_test_Model1","Estimate_test_Model2",
                  "Estimate_test_Model3","Estimate_test_Model4",
                  "Estimate_test_Model5","Estimate_test_Model6")))

AllGxEassoc_cor <- corDF[["GxE"]] %>%
  select(all_of(c("Estimate_test_Model1","Estimate_test_Model2",
                  "Estimate_test_Model3","Estimate_test_Model4",
                  "Estimate_test_Model5","Estimate_test_Model6")))

# Function to check sign consistency across rows
corDF[["E"]]$Sign_Consistent <- apply(AllEassoc_cor, 1, function(row) {
  non_na_values <- row[!is.na(row)]     # Exclude NA values
  all(non_na_values > 0) || all(non_na_values < 0) # Check if all are positive or all are negative
})

corDF[["GxE"]]$Sign_Consistent <- apply(AllGxEassoc_cor , 1, function(row) {
  non_na_values <- row[!is.na(row)]     # Exclude NA values
  all(non_na_values > 0) || all(non_na_values < 0) # Check if all are positive or all are negative
})


# Proportion that have Consistent Sign across Specification
Esign_PropConsistent <- sum(corDF[["E"]]$Sign_Consistent)/length(corDF[["E"]]$Sign_Consistent)
GxEsign_PropConsistent <- sum(corDF[["GxE"]]$Sign_Consistent)/length(corDF[["GxE"]]$Sign_Consistent)

# Amount of Associations that are significant:
Ecount <- AllEassoc_cor %>%
  summarize(across(everything(), ~ sum(!is.na(.)), .names = "SigAssoc_{col}"))

GxEcount <- AllGxEassoc_cor %>%
  summarize(across(everything(), ~ sum(!is.na(.)), .names = "SigAssoc_{col}"))

HEAPcount <- rbind(Ecount,GxEcount)
HEAPcount$ID <- c("E","GxE")

fwrite(HEAPcount, file = "/n/groups/patel/shakson_ukb/UK_Biobank/Output/Univariate/ReplicationCount.txt")
fwrite(HEAPcount, file = "/n/groups/patel/shakson_ukb/UK_Biobank/Figures/HEAP/T1/AssocReplicationCount.txt")


# Function to get correlation of significant E estimates across Covariate Specifications:
corSpec_HT <- function(df, title_ht){
  cor_matrix <- cor(df,use='pairwise.complete.obs')
  
  # Shorten the row and column names
  rownames(cor_matrix) <- gsub("Estimate_test_", "", rownames(cor_matrix))
  colnames(cor_matrix) <- gsub("Estimate_test_", "", colnames(cor_matrix))
  
  # Convert the correlation matrix to a long format
  cor_long <- melt(cor_matrix) %>%
    filter(as.numeric(Var1) >= as.numeric(Var2))
  
  # Plot the heatmap
  g1 <- ggplot(cor_long, aes(x = Var1, y = Var2, fill = value)) +
    geom_tile(color = "white") + # Add grid lines for separation
    geom_text(aes(label = sprintf("%.3f", value)), size = 3) + # Add correlation values
    scale_fill_gradient2(low = "purple", high = "orange", mid = "white",
                         midpoint = 0.95, limit = c(0.90, 1),
                         name = "r") + # Narrow color range
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
      axis.title = element_blank()
    ) +
    labs(title = title_ht)
  return(g1)
}
corEplot <- corSpec_HT(AllEassoc_cor, "Correlation of Covariate Specifications: E")
corGplot <- corSpec_HT(AllGxEassoc_cor, "Correlation of Covariate Specifications: GxE")

ggsave(paste0("./Figures/HEAP/F2/ECorAcrossCSpec.png"), 
       plot = corEplot, width = 5, height = 3, dpi = 1000)
ggsave(paste0("./Figures/HEAP/F3/GxECorAcrossCSpec.png"), 
       plot = corGplot, width = 5, height = 3, dpi = 1000)

# Example: Correlation between Type 1 and Type 5 Specification.
library(ggrepel)
diffAssoc <- corDF[["E"]] %>%
  filter(Sign_Consistent == FALSE)

diffAssoc$AssocID <- gsub("_f\\d+_\\d_\\d{2}","",diffAssoc$AssocID)
diffAssoc$AssocID <- gsub(":"," ", diffAssoc$AssocID)
diffAssoc$AssocID <- gsub("_"," ", diffAssoc$AssocID)

# E Associations:
T1v5 <- corDF[["E"]] %>%
  ggplot(aes(x = Estimate_test_Model1, y = Estimate_test_Model5)) +
  geom_point(size = 1) +
  geom_linerange(aes(ymin = Estimate_test_Model5 - CImargin_test_Model5, 
                     ymax = Estimate_test_Model5 + CImargin_test_Model5)) +  # Vertical error bar
  geom_linerange(aes(xmin = Estimate_test_Model1 - CImargin_test_Model1, 
                     xmax = Estimate_test_Model1 + CImargin_test_Model1)) + # Horizontal error bar
  geom_text_repel(data = diffAssoc, # Only label the subset
                  aes(label = str_wrap(AssocID,
                                       width = 10)),
                  max.overlaps = 2,
                  min.segment.length = unit(0, 'lines'),
                  nudge_y = -0.5,
                  nudge_x = 0.3,
                  size = 2.5) +
  labs(x = "Beta - Model1", y = "Beta - Model5") +
  theme_minimal()

ggsave(paste0("./Figures/HEAP/F2/T1v5.png"), 
       plot = T1v5, width = 4, height = 4, dpi = 1000)


# Color Version
T1v5color <- corDF[["E"]]  %>%
  ggplot(aes(x = Estimate_test_Model1, y = Estimate_test_Model5, 
             color = Category_test)) +
  geom_point(size = 0.2, alpha = 1) +
  geom_linerange(aes(ymin = Estimate_test_Model5 - CImargin_test_Model5, 
                     ymax = Estimate_test_Model5 + CImargin_test_Model5),
                 alpha = 1,
                 size = 0.2) +  # Vertical error bar
  geom_linerange(aes(xmin = Estimate_test_Model1 - CImargin_test_Model1, 
                     xmax = Estimate_test_Model1 + CImargin_test_Model1),
                 alpha = 1,
                 size = 0.2) + # Horizontal error bar
  labs(x = "Beta - Model1", y = "Beta - Model5") +
  theme_minimal()

ggsave(paste0("./Figures/HEAP/F2/T1v5cat.png"), 
       plot = T1v5color, width = 6, height = 4, dpi = 1000)




T1v5colorsplit <- corDF[["E"]]  %>%
  ggplot(aes(x = Estimate_test_Model1, y = Estimate_test_Model5,
             color = Category_test)) +
  geom_point(size = 0.2) +
  geom_linerange(aes(ymin = Estimate_test_Model5 - CImargin_test_Model5,
                     ymax = Estimate_test_Model5 + CImargin_test_Model5),
                 size = 0.2) +  # Vertical error bar
  geom_linerange(aes(xmin = Estimate_test_Model1 - CImargin_test_Model1,
                     xmax = Estimate_test_Model1 + CImargin_test_Model1),
                 size = 0.2) + # Horizontal error bar
  geom_smooth(method = "lm", se = TRUE,
              linewidth = 0.5) +
  theme_minimal() +
  facet_wrap(~ Category_test)

ggsave(paste0("./Figures/HEAP/F2/T1v5cat_split.png"), 
       plot = T1v5colorsplit, width = 10, height = 6, dpi = 1000)


#### GxE Plots:
GxET1v5 <- corDF[["GxE"]]  %>%
  ggplot(aes(x = Estimate_test_Model1, y = Estimate_test_Model5)) +
  geom_point(size = 1) +
  geom_linerange(aes(ymin = Estimate_test_Model5 - CImargin_test_Model5, 
                     ymax = Estimate_test_Model5 + CImargin_test_Model5)) +  # Vertical error bar
  geom_linerange(aes(xmin = Estimate_test_Model1 - CImargin_test_Model1, 
                     xmax = Estimate_test_Model1 + CImargin_test_Model1)) + # Horizontal error bar
  labs(x = "Beta - Model1", y = "Beta - Model5") +
  xlim(-0.4, 0.4) +
  theme_minimal()
ggsave(paste0("./Figures/HEAP/F3/gxeT1v5.png"), 
       plot = GxET1v5, width = 3, height = 3, dpi = 1000)


GxET1v5_color <- corDF[["GxE"]]  %>%
  ggplot(aes(x = Estimate_test_Model1, y = Estimate_test_Model5,
             color = Category_test)) +
  geom_point(size = 1) +
  geom_linerange(aes(ymin = Estimate_test_Model5 - CImargin_test_Model5, 
                     ymax = Estimate_test_Model5 + CImargin_test_Model5),
                 width = 0.1) +  # Vertical error bar
  geom_linerange(aes(xmin = Estimate_test_Model1 - CImargin_test_Model1, 
                     xmax = Estimate_test_Model1 + CImargin_test_Model1), 
                 height = 0.1) + # Horizontal error bar
  labs(x = "Beta - Model1", y = "Beta - Model5") +
  xlim(-0.4, 0.4) +
  theme_minimal()
ggsave(paste0("./Figures/HEAP/F3/gxeT1v5_cat.png"), 
       plot = GxET1v5_color, width = 5, height = 3, dpi = 1000)


#'*DEALING WITH TRENDS - IDENTIFYING TRENDS OF SIGNFICANT and REPLICATED ASSOCIATIONS!*


#'# Function to Compare Association Concordance Across Specifications:
AssocTrends_DF <- function(CovarSpecAssocList){
  Trends <- list()
  
  CovarSpecAssocCor <- lapply(names(CovarSpecAssocList), function(x){
    CovarSpecAssocList[[x]][["E"]][["sigBOTH"]]$CovarType <- x
    CovarSpecAssocList[[x]][["GxE"]][["sigBOTH"]]$CovarType <- x
    return(CovarSpecAssocList[[x]])
  })
  names(CovarSpecAssocCor) <- names(CovarSpecAssocList)
  
  # Aggregate Info across Specifications:
  Eassoc_AllType <- lapply(names(CovarSpecAssocCor), function(x){
    return(CovarSpecAssocCor[[x]][["E"]][["sigBOTH"]])
  }) %>% reduce(full_join)
  
  GxEassoc_AllType <- lapply(names(CovarSpecAssocCor), function(x){
    return(CovarSpecAssocCor[[x]][["GxE"]][["sigBOTH"]])
  }) %>% reduce(full_join)
  
  Trends[["E"]] <- Eassoc_AllType
  Trends[["GxE"]] <- GxEassoc_AllType
  
  return(Trends)
}

HEAPtrends <- AssocTrends_DF(HEAPassoc@HEAPsig)

EassocCat_mean_abs <- HEAPtrends[["E"]] %>%
  group_by(Category_test, CovarType) %>%
  summarize(
    mean_effect_size = mean(abs(Estimate_test)),#weighted.mean(Estimate, w = 1 / `Std. Error`^2),
    se_mean_effect_size = sqrt(sum(`Std.Error_test`^2) / n()),
    n_observations = n(),
    .groups = 'drop'
  )

EassocCat_mean_abs %>%
  ggplot(aes(x=reorder(Category_test, -mean_effect_size), 
             y=mean_effect_size, 
             color = CovarType)) +
  geom_jitter(width = 0.3, height = 0.001) +
  ylim(0,0.4) +
  labs(
    x = NULL,
    y = "Absolute Mean Effect Size"
  ) +
  theme_minimal() +
  theme(axis.text.x=element_text(angle=60, hjust=1))



###GET Meta-Analysis Point Estimates of Effect Size for each Category!
library(metafor)

# Perform random-effects meta-analysis for each unique Category_test
TrendsRes <- EassocCat_mean_abs %>%
  group_by(Category_test) %>%
  summarise(
    Beta = as.numeric(rma(yi = mean_effect_size, 
                          sei = se_mean_effect_size, 
                          method = "REML", 
                          data = pick(everything()))$beta),
    SE = rma(yi = mean_effect_size, 
             sei = se_mean_effect_size, 
             method = "REML", 
             data = pick(everything()))$se,
    CI_lower =  rma(yi = mean_effect_size, 
                    sei = se_mean_effect_size, 
                    method = "REML", 
                    data = pick(everything()))$ci.lb,
    CI_upper = rma(yi = mean_effect_size, 
                   sei = se_mean_effect_size, 
                   method = "REML", 
                   data = pick(everything()))$ci.ub,
    .groups = 'drop'
  ) %>%
  rename(ID = Category_test)

# PLOT Averaged Main E Association Effects:
#qt(0.975, df = 5) - DO not do CIs here --> Not sure if VALID...

ECatTrends <- TrendsRes %>%
  ggplot(aes(x = reorder(ID, -Beta), y = Beta, color = ID)) +
  geom_point() +
  geom_linerange(aes(ymin = CI_lower, 
                     ymax = CI_upper)) +
  labs(x = NULL,
       y = "Absolute Mean Effect Size") +
  theme_minimal() +
  theme(axis.text.x=element_text(angle=60, hjust=1),
        legend.position="none")
ggsave(paste0("./Figures/HEAP/F2/EmeanabsTrends.png"), 
       plot = ECatTrends, width = 5, height = 4, dpi = 1000)



#'*TRENDS for certain E features*
# REMEMBER THIS IS FOR THE SET OF ALL SIGNIFICANT ASSOCIATIONS!!!
Eassoc_mean_abs <- HEAPtrends[["E"]]  %>%
  group_by(ID, CovarType) %>%
  summarize(
    mean_effect_size = mean(abs(Estimate_test)),#weighted.mean(Estimate, w = 1 / `Std. Error`^2),
    se_mean_effect_size = sqrt(sum(`Std.Error_test`^2) / n()),
    n_observations = n(),
    .groups = 'drop'
  )

E_TrendsRes <- Eassoc_mean_abs %>%
  filter(grepl("alcohol_intake_frequency",ID)) %>%
  group_by(ID) %>%
  summarise(
    Beta = as.numeric(rma(yi = mean_effect_size, 
                          sei = se_mean_effect_size, 
                          method = "REML", 
                          data = pick(everything()))$beta),
    SE = rma(yi = mean_effect_size, 
             sei = se_mean_effect_size, 
             method = "REML", 
             data = pick(everything()))$se,
    CI_lower =  rma(yi = mean_effect_size, 
                    sei = se_mean_effect_size, 
                    method = "REML", 
                    data = pick(everything()))$ci.lb,
    CI_upper = rma(yi = mean_effect_size, 
                   sei = se_mean_effect_size, 
                   method = "REML", 
                   data = pick(everything()))$ci.ub,
    .groups = 'drop'
  )

AlcTrend <- E_TrendsRes %>%
  ggplot(aes(x = reorder(ID, Beta), y = Beta)) +
  geom_point() +
  geom_linerange(aes(ymin = CI_lower, 
                     ymax = CI_upper)) +
  labs(x = "Alcohol Intake",
       y = "Absolute Mean Effect Size") +
  scale_x_discrete(labels=c("alcohol_intake_frequency_f1558_0_02" = "Special Occasions", 
                            "alcohol_intake_frequency_f1558_0_03" = "1-3 Times Month",
                            "alcohol_intake_frequency_f1558_0_04" = "1-2 Times Week",
                            "alcohol_intake_frequency_f1558_0_05" = "3-4 Times Week",
                            "alcohol_intake_frequency_f1558_0_06" = "Daily/Almost Daily")) +
  theme_minimal() +
  theme(axis.text.x=element_text(angle=60, hjust=1),
        legend.position="none")
ggsave(paste0("./Figures/HEAP/F2/AlcTrends.png"), 
       plot = AlcTrend, width = 5, height = 4, dpi = 1000)



#Other Trends: Exercise Intensity:
E_TrendsRes <- Eassoc_mean_abs %>%
  filter(grepl("types_of_physical_activity",ID)) %>%
  group_by(ID) %>%
  summarise(
    Beta = as.numeric(rma(yi = mean_effect_size, 
                          sei = se_mean_effect_size, 
                          method = "REML", 
                          data = pick(everything()))$beta),
    SE = rma(yi = mean_effect_size, 
             sei = se_mean_effect_size, 
             method = "REML", 
             data = pick(everything()))$se,
    CI_lower =  rma(yi = mean_effect_size, 
                    sei = se_mean_effect_size, 
                    method = "REML", 
                    data = pick(everything()))$ci.lb,
    CI_upper = rma(yi = mean_effect_size, 
                   sei = se_mean_effect_size, 
                   method = "REML", 
                   data = pick(everything()))$ci.ub,
    .groups = 'drop'
  )

ExerTrend <- E_TrendsRes %>%
                  ggplot(aes(x = reorder(ID, Beta), y = Beta)) +
                  geom_point() +
                  geom_linerange(aes(ymin = CI_lower, 
                                     ymax = CI_upper)) +
                  labs(x = "Exercise Type",
                       y = "Absolute Mean Effect Size") +
                  scale_x_discrete(labels=c("types_of_physical_activity_in_last_4_weeks_f6164_0_0.multi_Heavy_DIY_.eg._weeding._lawn_mowing._carpentry._digging." = "HeavyDIY: gardening, carpentry, digging",  
                                            "types_of_physical_activity_in_last_4_weeks_f6164_0_0.multi_Light_DIY_.eg._pruning._watering_the_lawn." = "LightDIY: prune, water lawn",              
                                            "types_of_physical_activity_in_last_4_weeks_f6164_0_0.multi_None_of_the_above" = "None",                                     
                                            "types_of_physical_activity_in_last_4_weeks_f6164_0_0.multi_Other_exercises_.eg._swimming._cycling._keep_fit._bowling." = "Swim/Cycle/KeepFit",
                                            "types_of_physical_activity_in_last_4_weeks_f6164_0_0.multi_Strenuous_sports" = "StrenuousSports",                                         
                                            "types_of_physical_activity_in_last_4_weeks_f6164_0_0.multi_Walking_for_pleasure_.not_as_a_means_of_transport." = "Walking")) +
                  theme_minimal() +
                  theme(axis.text.x=element_text(angle=60, hjust=1),
                        legend.position="none")

ggsave(paste0("./Figures/HEAP/F2/ExerciseTrends.png"), 
       plot = ExerTrend, width = 5, height = 6, dpi = 1000)


#'*Make function later*
# HEAPtrend <- function(label){
#   E_TrendsRes <- Eassoc_mean_abs %>%
#     filter(grepl(label,ID)) %>%
#     group_by(ID) %>%
#     summarise(
#       Beta = as.numeric(rma(yi = mean_effect_size, 
#                             sei = se_mean_effect_size, 
#                             method = "REML", 
#                             data = pick(everything()))$beta),
#       SE = rma(yi = mean_effect_size, 
#                sei = se_mean_effect_size, 
#                method = "REML", 
#                data = pick(everything()))$se,
#       CI_lower =  rma(yi = mean_effect_size, 
#                       sei = se_mean_effect_size, 
#                       method = "REML", 
#                       data = pick(everything()))$ci.lb,
#       CI_upper = rma(yi = mean_effect_size, 
#                      sei = se_mean_effect_size, 
#                      method = "REML", 
#                      data = pick(everything()))$ci.ub,
#       .groups = 'drop'
#     )
#   return(E_TrendsRes)
# }


