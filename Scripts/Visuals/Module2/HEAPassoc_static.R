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
library(qs)

#'*LOAD CovarSpecList*
setwd("/n/groups/patel/shakson_ukb/UK_Biobank/")
HEAPassoc <- qread("./Output/HEAPres/HEAPassoc.qs")

CovarSpecList <- HEAPassoc@HEAPlist
CovarSpecAssocList <- HEAPassoc@HEAPsig

#Also functions:
# FUNCTION: Miami Plot Single Omic
#'*Decide whether to adjust significance to both thresholds - per protein bonferonni or all proteins*
miamiplot_single_omic <- function(protID, assoc, title){
  total <- assoc %>% 
    filter(omicID == protID) %>%
    arrange(Category)
  
  #'*added to make sure the p-values are correctly shown in plot:*
  #total$`Pr(>|t|)` <- ifelse(total$`Pr(>|t|)` == 0, .Machine$double.xmin, total$`Pr(>|t|)`)
  
  
  order <- 1:nrow(total)
  p <- ggplot(total, aes(x = order, y = -log10(`Pr(>|t|)` + 1e-300) * sign(Estimate), color = Category)) +
    geom_point(aes(text = paste("ID:", ID, 
                                "<br>Beta:", signif(Estimate, digits = 2),
                                "<br>Pval:", signif(`Pr(>|t|)`, digits = 2),
                                "<br>N:", samplesize,
                                "<br>Category:", Category,
                                "<br>Total R2:", signif(adj.R2, digits = 2)))) +
    geom_abline(intercept = -log10((0.05/nrow(total)) + 1e-300), slope = 0, 
                color = "blue", linetype = "dashed") +
    geom_abline(intercept = log10((0.05/nrow(total)) + 1e-300), slope = 0, 
                color = "blue", linetype = "dashed") +
    theme_minimal() +
    labs(x = "Exposures", y = "-log10(P-value) * sign(Beta)") +
    ggtitle(title)
  
  print(p)
  
  ggplotly(p, tooltip = "text")
  return(p)
}

# FUNCTION: Save Interactive Plot:
plotsave <- function(gg, protID, type){
  gginteract <- ggplotly(gg, tooltip = "text")
  
  saveWidget(gginteract, 
             file = paste0("./Figures/HEAP/CaseExamples/Interactive/",
                           protID,"_",type,"assoc.html"))
  
  ggsave(paste0("./Figures/HEAP/CaseExamples/Individual/",
                protID,"_",type,"assoc.png"), 
         plot = gg, width = 10, height = 6, dpi = 1000)
}

#'*SINGLE Omic Signatures:*
saveOmicSignatures <- function(omicID, CovarType, plotname){
  gg1 <- miamiplot_single_omic(protID = omicID, 
                               assoc = CovarSpecList[[CovarType]]$train[[1]],
                               title = paste0(omicID,": ",plotname))
  plotsave(gg1, omicID, CovarType)
} 
saveOmicSignatures("GDF15","Model5","E Associations")
saveOmicSignatures("CXCL17","Model5","E Associations")
saveOmicSignatures("MAMDC4","Model5","E Associations")
saveOmicSignatures("LAMP3","Model5","E Associations")
saveOmicSignatures("WFDC2","Model5","E Associations")
saveOmicSignatures("CEACAM16","Model5","E Associations")
saveOmicSignatures("LEP","Model5","E Associations")
saveOmicSignatures("IGFBP2","Model5","E Associations")
saveOmicSignatures("IGFBP1","Model5","E Associations")
saveOmicSignatures("CKB","Model5","E Associations")
saveOmicSignatures("FABP4","Model5","E Associations")

#'*FIGURE out how to highlight replicated hits later.*

#'*Specific Associations E or GxE*
#Source univarInteract Script:
source("/n/groups/patel/shakson_ukb/UK_Biobank/RScripts/Pure_StatGen/Prot_ExPGS/runProt_PGS_univarInteract_v2.R")

#'*To Modify Later --> For # of days,etc. that got converted into z-score.*
#'*Discretize the above by days and find the respective quintiles.*

#Run specific associations:
saveGWISplot <- function(o,e,c,n,l=c(),p,f){
  gg <- runGWISplot(omic = o,
                    exposure = e,
                    CType = c,
                    ename = n,
                    levels = l,
                    plotname = p)
  ggsave(paste0("./Figures/HEAP/CaseExamples/GWIS/",paste0(o,"_",f),".png"), 
         plot = gg, width = 5, height = 4, dpi = 1000)
}

#Find top E hits that REPLICATE to show
Model5Eres <- CovarSpecAssocList$Model5$E$sigBOTH
sort(unique(Model5Eres$ID))

runGWISplot(omic = "IGFBP2",
            exposure = "pork_intake_f1389_0_0",
            CType = "Type6",
            ename = "Pork Intake",
            plotname = "Exposure Effect Conditioned on Genetics")


runGWISplot(omic = "CXCL17",
            exposure = "current_tobacco_smoking_f1239_0_0_Yes._on_most_or_all_days",
            CType = "Type6",
            ename = "Smoking Status",
            plotname = "Exposure Effect Conditioned on Genetics")

#Diet Related Proteins
#Hard to Find Solid Evidence for this based on trends
#Increase in exposure is not equivalent to increase in protein expr


### Alcohol Related Proteins
saveGWISplot(o = "MAMDC4",
             e = "alcohol_drinker_status_f20117_0_0_Current",
             c = "Type6",
             n = "Alcohol User",
             l = c("NO","YES"),
             p = "Exposure Effect Conditioned on Genetics",
             f = "AlcoholUse_E")
saveGWISplot(o = "MAMDC4",
             e = "alcohol_intake_frequency_f1558_0_0",
             c = "Type6",
             n = "Alcohol Freq.",
             l = c("Never","Special Occasions",
                   "1-3 times monthly","1-2 times weekly",
                   "3-4 times weekly","Daily"),
             p = "Exposure Effect Conditioned on Genetics",
             f = "AlcoholFreq_E")

saveGWISplot(o = "CEACAM16",
             e = "alcohol_drinker_status_f20117_0_0_Current",
             c = "Type6",
             n = "Alcohol User",
             l = c("NO","YES"),
             p = "Exposure Effect Conditioned on Genetics",
             f = "AlcoholUse_E")
saveGWISplot(o = "CEACAM16",
             e = "alcohol_intake_frequency_f1558_0_0",
             c = "Type6",
             n = "Alcohol Freq.",
             l = c("Never","Special Occasions",
                   "1-3 times monthly","1-2 times weekly",
                   "3-4 times weekly","Daily"),
             p = "Exposure Effect Conditioned on Genetics",
             f = "AlcoholFreq_E")





### Exercise Related proteins
saveGWISplot(o = "LEP",
             e = "types_of_physical_activity_in_last_4_weeks_f6164_0_0.multi_Strenuous_sports",
             c = "Type6",
             n = "Strenuous Sports",
             l = c("NO","YES"),
             p = "Exposure Effect Conditioned on Genetics",
             f = "StrenSports_E")
saveGWISplot(o = "LEP",
             e = "usual_walking_pace_f924_0_0",
             c = "Type6",
             n = "Walking Pace",
             #l = c("NO","YES"),
             p = "Exposure Effect Conditioned on Genetics",
             f = "WalkPace_E")
saveGWISplot(o = "LEP",
             e = "met_minutes_per_week_for_vigorous_activity_f22039_0_0" ,
             c = "Type6",
             n = "Vigorous Activity \n (Z-score)",
             #l = c("NO","YES"),
             p = "Exposure Effect Conditioned on Genetics",
             f = "Min_VigorousPhysAct_E")
saveGWISplot(o = "LEP",
             e = "number_of_days_week_of_vigorous_physical_activity_10_plus_minutes_f904_0_0",
             c = "Type6",
             n = "Vigorous Activity \n (Z-score)",
             #l = c("NO","YES"),
             p = "Exposure Effect Conditioned on Genetics",
             f = "Days_VigorousPhysAct_E")


saveGWISplot(o = "FABP4",
             e = "types_of_physical_activity_in_last_4_weeks_f6164_0_0.multi_Strenuous_sports",
             c = "Type6",
             n = "Strenuous Sports",
             l = c("NO","YES"),
             p = "Exposure Effect Conditioned on Genetics",
             f = "StrenSports_E")
saveGWISplot(o = "FABP4",
             e = "number_of_days_week_of_vigorous_physical_activity_10_plus_minutes_f904_0_0",
             c = "Type6",
             n = "Vigorous Activity \n (Z-score)",
             #l = c("NO","YES"),
             p = "Exposure Effect Conditioned on Genetics",
             f = "Days_VigorousPhysAct_E")


saveGWISplot(o = "IGFBP1",
             e = "number_of_days_week_of_vigorous_physical_activity_10_plus_minutes_f904_0_0",
             c = "Type6",
             n = "Vigorous Activity \n (Z-score)",
             #l = c("NO","YES"),
             p = "Exposure Effect Conditioned on Genetics",
             f = "Days_VigorousPhysAct_E")
saveGWISplot(o = "IGFBP1",
             e = "summed_days_activity_f22033_0_0",
             c = "Type6",
             n = "Days of Activity \n (Z-score)",
             #l = c("NO","YES"),
             p = "Exposure Effect Conditioned on Genetics",
             f = "Days_SummedDaysPhysAct_E")


### Smoking Related Proteins:
saveGWISplot(o = "CXCL17",
             e = "current_tobacco_smoking_f1239_0_0_Yes._on_most_or_all_days",
             c = "Type6",
             n = "Daily Smoker",
             l = c("NO","YES"),
             p = "Exposure Effect Conditioned on Genetics",
             f = "DailySmoker_E")

saveGWISplot(o = "LAMP3",
             e = "current_tobacco_smoking_f1239_0_0_Yes._on_most_or_all_days",
             c = "Type6",
             n = "Daily Smoker",
             l = c("NO","YES"),
             p = "Exposure Effect Conditioned on Genetics",
             f = "DailySmoker_E")

saveGWISplot(o = "TNR",
             e = "current_tobacco_smoking_f1239_0_0_Yes._on_most_or_all_days",
             c = "Type6",
             n = "Daily Smoker",
             l = c("NO","YES"),
             p = "Exposure Effect Conditioned on Genetics",
             f = "DailySmoker_E")

saveGWISplot(o = "WFDC2",
             e = "current_tobacco_smoking_f1239_0_0_Yes._on_most_or_all_days",
             c = "Type6",
             n = "Daily Smoker",
             l = c("NO","YES"),
             p = "Exposure Effect Conditioned on Genetics",
             f = "DailySmoker_E")

#Find top GxE hits that REPLICATE to show!!
Type6GxEres <- CovarSpecAssocList$Type6$GxE$sigBOTH

saveGWISplot(o = "APOF",
             e = "alcohol_intake_frequency_f1558_0_0",
             c = "Type6",
             n = "Alcohol Intake Freq.",
             l = c("Never","Special Occasions",
                   "1-3 times monthly","1-2 times weekly",
                   "3-4 times weekly","Daily"),
             p = "Polygenic GxE Interaction",
             f = "Alcohol_GxE")
saveGWISplot(o = "ALPP",
             e = "current_tobacco_smoking_f1239_0_0_Yes._on_most_or_all_days",
             c = "Type6",
             n = "Daily Smoker",
             l = c("NO","YES"),
             p = "Polygenic GxE Interaction",
             f = "Smoker_GxE")

