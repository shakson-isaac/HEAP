#Libraries:
library(data.table)
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(ggpmisc)
library(pbapply)
library(dplyr)
library(purrr)
setwd("/n/groups/patel/shakson_ukb/UK_Biobank/")

#TODO:
#Figures:
# General overall averaged result of mediation
# Variation of GEM statistic
# Highlight Specific Diseases (HR and c-index)
# Highlight E v G partitioning (Plot)

#'*TIPS VISUALS -- E vs G plots*
#'HAVE geom_point come after geom_ribbon to make sure the points POP OUT!!!

#'*To-ADD Algorithm Later*
#''*Need p-value for the HR of the protein itself (to get more accurate x-axis)*
#'*IDEA --> get variability across covariate specifications:*
#'*To Add Later -- > Get Variability of above via bootstrapping these associations many times*

#Save all summary stats as excel files:

# Load Files w/ Progress Bar:
loadMDAssoc <- function(covarType){
  error_files <- list()  # Track files w/ errors
  stat_all <- list()
  
  # Use pblapply for progress bar and iteration
  for(i in 1:1000){
    tryCatch({
      load <- fread(file=paste0("/n/groups/patel/shakson_ukb/UK_Biobank/Data/Parallel/Prot_PGSPXS_mediation/",
                    covarType,"/MDres_",i,".txt"), 
                    select = 2:16)
      
      stat_all[[i]] <- load
      
    }, error = function(e) {
      message(paste("Error reading file:", i))
      error_files <<- append(error_files, i)
    })
  }
  
  
  stat_all <- do.call("rbind", stat_all)
  stat_all <- stat_all %>%
                mutate(E_HRi = exp(`Exposure Indirect Effect`),
                       G_HRi = exp(`Genetic Indirect Effect`),
                       delta_HRi = exp(abs(`Exposure Indirect Effect` - `Genetic Indirect Effect`)),
                       total_HRi = exp(`Exposure Indirect Effect` + `Genetic Indirect Effect`),
                       Eprop_mediated = abs(`Exposure Indirect Effect`)/(abs(`Exposure Indirect Effect`) + abs(`Exposure Direct Effect`)),
                       Gprop_mediated = abs(`Genetic Indirect Effect`)/(abs(`Genetic Indirect Effect`) + abs(`Genetic Direct Effect`)))
  stat_all$CovarSpec <- covarType
  stat_all$HR_upper95 <- as.numeric(stat_all$HR_upper95)
  
  
  return(stat_all)
}

# Load results:
Type1 <- loadMDAssoc(covarType = "Type1")
Type2 <- loadMDAssoc(covarType = "Type2")
Type3 <- loadMDAssoc(covarType = "Type3")
Type4 <- loadMDAssoc(covarType = "Type4")
Type5 <- loadMDAssoc(covarType = "Type5")

AllCovarSpec <- list(Type1, Type2, Type3, Type4, Type5) %>% reduce(full_join)

# FINAL ANALYSIS ---------------------

# Obtain OLS slope across Specifications:
oSlope <- function(MDres){
  
  # Obtain GEM: Expected Indirect Effect Slope: log(E/G)
  protContext <- pblapply(unique(MDres$ID), function(p){
    df <- MDres %>% filter(ID == p) 
    
    fit <- lm(`Exposure Indirect Effect` ~ `Genetic Indirect Effect`,
              data = df)
    
    stats <- summary(fit)$coefficients

    maxE <- max(abs(df$`Exposure Indirect Effect`))
    maxG <- max(abs(df$`Genetic Indirect Effect`))

    df <- MDres %>%
      filter(ID == p) %>%
      filter(!((HR_lower95 < 1 & HR_upper95 > 1) | (HR_lower95 > 1 & HR_upper95 < 1)))
    
    num_diseases <- length(unique(df$DZid))
    
    protStats <- data.frame(
      ID = p,
      Estimate = stats[2,"Estimate"],
      Std.Error = stats[2,"Std. Error"],
      Pvalue = stats[2,"Pr(>|t|)"],
      max.E = maxE,
      max.G = maxG,
      Estimate_Norm = stats[2,"Estimate"] * (maxE/maxG),
      NumDiseases = num_diseases
    )
    
    #Below return the slope of the above regression
    return(protStats) 
  })
  protContext <- do.call(rbind, protContext)
  
  return(protContext)
}
Type1_o <- oSlope(Type1)
Type2_o <- oSlope(Type2)
Type3_o <- oSlope(Type3)
Type4_o <- oSlope(Type4)
Type5_o <- oSlope(Type5)

#Merge results together:
Type1_o$CovarType <- "Type1"
Type2_o$CovarType <- "Type2"
Type3_o$CovarType <- "Type3"
Type4_o$CovarType <- "Type4"
Type5_o$CovarType <- "Type5"


MDres <- list(Type1_o, Type2_o, Type3_o, Type4_o, Type5_o)
names(MDres) <- c("Type1", "Type2", "Type3", "Type4", "Type5")

MD_GvE_effects <- map2_dfr(MDres, names(MDres), ~ .x %>%
                            select(all_of(c("ID", "Estimate", "CovarType"))))

#GET GEM Stat:
Type5_o <- Type5_o %>%
  mutate(GEM = log(Estimate))



#Plot of std.error:
GEMres <- MD_GvE_effects %>% 
            group_by(ID) %>%
            summarise(meanGEM = mean(log(Estimate)),
                      stdErrorGEM = sd(log(Estimate))/sqrt(n()))

# GEMres <- MD_GvE_effects %>% 
#             filter(CovarType != "Type1") %>%
#             group_by(ID) %>%
#             summarise(meanGEM = mean(log(Estimate)),
#                       stdErrorGEM = sd(log(Estimate))/sqrt(n()))

subset_GEMres <- GEMres %>%
  filter(stdErrorGEM > 0.1)

ggVar <- GEMres %>%
  ggplot(aes(x = meanGEM, y = stdErrorGEM)) +
      geom_point() +
      labs(x = "GEM",
           y ="Std. Error",
           title = "Variation across Covariate Specification") +
      geom_text_repel(data = subset_GEMres,
                      aes(x = meanGEM, 
                      y = stdErrorGEM,
                      label = ID),
                      size = 2.5) +
      theme_minimal()
lapply(c("png", "svg", "pdf"), 
       function(fmt) ggsave(paste0("./Figures/HEAP/SF4/",
                                   "GEM_Variability.", fmt), 
                            plot = ggVar, width = 5, height = 4, dpi = 500))



library(purrr)

# Obtain correlations of Estimate with Num of Diseases:
datasets <- MDres

# Function for correlations
get_cor_results <- function(data) {
  cleaned_data <- na.omit(data)
  test_result <- cor.test(log(cleaned_data$Estimate), cleaned_data$NumDiseases)
  data.frame(
    correlation = test_result$estimate,
    p_value = test_result$p.value,
    ci_lower = test_result$conf.int[1],
    ci_upper = test_result$conf.int[2]
  )
}


results_df <- map_dfr(datasets, get_cor_results, .id = "specification")

#Plot: correlation
ggPL <- ggplot(results_df, aes(x = specification, y = correlation, ymin = ci_lower, ymax = ci_upper)) +
  geom_bar(stat = "identity", fill = "darkgray", color = "black") +
  geom_errorbar(width = 0.2) + 
  geom_text(aes(label = sprintf("p = %.e", p_value)), vjust = -7, size = 3) +  # Add p-value as text label
  labs(x = "Specifications", y = "Correlation", title = "Correlation: GEM vs Num. Diseases Across Specifications") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylim(0, 0.25)

lapply(c("png", "svg", "pdf"), 
       function(fmt) ggsave(paste0("./Figures/HEAP/SF4/",
                                   "GEM_Pleiotropic.", fmt), 
                            plot = ggPL, width = 6, height = 4, dpi = 500))



#Plot of GvE Mediation Effects: Only Type 5 for clarity in Manuscript:
subset_protContext <- Type5_o %>%
  filter(Estimate > 1 | Estimate < 0.8)


m1 <- Type5_o %>%
          ggplot(aes(x = NumDiseases, y = log(Estimate))) +
          #Color Environment and Genetics Driven Proteins:
          geom_ribbon(aes(ymin =-Inf, ymax = 0, fill = "Genetics Driven"), alpha = 0.2) +
          geom_ribbon(aes(ymin = 0, ymax = Inf, fill = "Environment Driven"), alpha = 0.2) +
          scale_fill_manual(values = c("Genetics Driven" = "blue", "Environment Driven" = "green")) +
          geom_point() +
          geom_text_repel(data=subset_protContext, 
                          aes(x = NumDiseases, y = log(Estimate), 
                              label= ID),
                          size = 2.5) +
          labs(x = "Number of Significant Disease Hits",
               y ="GEM",
               fill = "Region") +
          theme_minimal() +
          theme(legend.position = "bottom")

lapply(c("png", "svg", "pdf"), 
       function(fmt) ggsave(paste0("./Figures/HEAP/F4/",
                                   "GEM_MD.", fmt), 
                            plot = m1, width = 8, height = 5, dpi = 500))


# Function: Volcano Plot of Protein Associations to a Given Disease
protVolcanoPlot <- function(covarDF, dzID, R2cutoff, title){
  DZ_res <- covarDF %>%
    filter(DZid %in% dzID) %>%
    filter(R2 > R2cutoff)
  
  prot_labels <- DZ_res %>%
    filter(HR < 0.75 | HR > 1.5)
  
  # Find the maximum deviation from 1
  max_deviation <- max(abs(DZ_res$HR - 1))
  symmetric_limit <- 1 + max_deviation * c(-1, 1)
  
  
  gg1 <- DZ_res %>%
    ggplot(aes(x = HR, y = cindex)) +
    geom_point() +
    geom_text_repel(data = prot_labels, aes(label = ID)) +
    scale_x_continuous(limits = symmetric_limit)  +
    labs(x = "Hazard Ratio",
         y = "c-index") +
    ggtitle(title) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
  
  return(gg1)          
}


#For Figure 5:
library(patchwork)
m2.1 <- protVolcanoPlot(covarDF = Type5,
                dzID = "age_e11_first_reported_non_insulin_dependent_diabetes_mellitus_f130708_0_0",
                R2cutoff = 0.02,
                title = "Type 2 Diabetes")
m2.2 <- protVolcanoPlot(covarDF = Type5,
                dzID = "age_j43_first_reported_emphysema_f131490_0_0",
                R2cutoff = 0.02,
                title = "Emphysema")
m2.3 <- protVolcanoPlot(covarDF = Type5,
                        dzID = "age_f10_first_reported_mental_and_behavioural_disorders_due_to_use_of_alcohol_f130854_0_0",
                        R2cutoff = 0.02,
                        title = "Mental Disorders Due to Alcohol Usage")

MD1 <- (m2.1 | m2.2 | m2.3)


lapply(c("png", "svg", "pdf"), 
       function(fmt) ggsave(paste0("./Figures/HEAP/F4/", "ProtWideDiseaseAssoc.", fmt), 
                            plot = MD1, width = 12, height = 4, dpi = 500))

volcanoplots <- function(){

protVolcanoPlot(covarDF = Type5,
                dzID = "age_e11_first_reported_non_insulin_dependent_diabetes_mellitus_f130708_0_0",
                R2cutoff = 0.02,
                title = "Proteome Wide Analysis: Type 1 Diabetes")

protVolcanoPlot(covarDF = Type5,
                dzID = "age_f10_first_reported_mental_and_behavioural_disorders_due_to_use_of_alcohol_f130854_0_0",
                R2cutoff = 0.02,
                title = "Proteome Wide Analysis: Mental Disorders Due to Alcohol Usage")


protVolcanoPlot(covarDF = Type5,
                dzID = "age_n18_first_reported_chronic_renal_failure_f132032_0_0",
                R2cutoff = 0.02,
                title = "Proteome Wide Analysis: Chronic Renal Failure")

protVolcanoPlot(covarDF = Type1,
                dzID = "age_e66_first_reported_obesity_f130792_0_0",
                R2cutoff = 0.02,
                title = "Proteome Wide Analysis: Incident Obesity")


protVolcanoPlot(covarDF = Type5,
                dzID = "age_e66_first_reported_obesity_f130792_0_0",
                R2cutoff = 0.02,
                title = "Proteome Wide Analysis: Incident Obesity")

#age_j84_first_reported_other_interstitial_pulmonary_diseases_f131528_0_0
# PRSS8, CXCL17, ALPP

protVolcanoPlot(covarDF = Type5,
                dzID = "age_j84_first_reported_other_interstitial_pulmonary_diseases_f131528_0_0",
                R2cutoff = 0.02,
                title = "Proteome Wide Analysis: Interstitial Pulmonary Diseases")


protVolcanoPlot(covarDF = Type5,
                dzID = "age_j44_first_reported_other_chronic_obstructive_pulmonary_disease_f131492_0_0",
                R2cutoff = 0.02,
                title = "Proteome Wide Analysis: COPD")

protVolcanoPlot(covarDF = Type5,
                dzID = "age_j43_first_reported_emphysema_f131490_0_0",
                R2cutoff = 0.02,
                title = "Proteome Wide Analysis: Emphysema")

protVolcanoPlot(covarDF = Type5,
                dzID = "age_k76_first_reported_other_diseases_of_liver_f131670_0_0",
                R2cutoff = 0.02,
                title = "Proteome Wide Analysis: Liver Disease")


protVolcanoPlot(covarDF = Type5,
                dzID = "age_i50_first_reported_heart_failure_f131354_0_0",
                R2cutoff = 0.02,
                title = "Proteome Wide Analysis: Heart Failure")


protVolcanoPlot(covarDF = Type5,
                dzID = "age_i27_first_reported_other_pulmonary_heart_diseases_f131310_0_0",
                R2cutoff = 0.02,
                title = "Proteome Wide Analysis: Heart Failure")

protVolcanoPlot(covarDF = Type5,
                dzID = "age_i64_first_reported_stroke_not_specified_as_haemorrhage_or_infarction_f131368_0_0",
                R2cutoff = 0.02,
                title = "Proteome Wide Analysis: Stroke")

protVolcanoPlot(covarDF = Type5,
                dzID = "age_g30_first_reported_alzheimers_disease_f131036_0_0",
                R2cutoff = 0.02,
                title = "Proteome Wide Analysis: Alzheimers")

protVolcanoPlot(covarDF = Type5,
                dzID = "age_f03_first_reported_unspecified_dementia_f130842_0_0",
                R2cutoff = 0.02,
                title = "Proteome Wide Analysis: Dementia")

protVolcanoPlot(covarDF = Type5,
                dzID = "age_f00_first_reported_dementia_in_alzheimers_disease_f130836_0_0",
                R2cutoff = 0.02,
                title = "Proteome Wide Analysis: Dementia")

protVolcanoPlot(covarDF = Type5,
                dzID = "age_i70_first_reported_atherosclerosis_f131380_0_0",
                R2cutoff = 0.02,
                title = "Proteome Wide Analysis: Atherosclerosis")



#Suggestions: CKD, Heart Failure, Gout

protVolcanoPlot(covarDF = Type5,
                dzID = "age_n18_first_reported_chronic_renal_failure_f132032_0_0",
                R2cutoff = 0.02,
                title = "Proteome Wide Analysis: Chronic Renal Failure")

protVolcanoPlot(covarDF = Type5,
                dzID = "age_m10_first_reported_gout_f131858_0_0",
                R2cutoff = 0.02,
                title = "Proteome Wide Analysis: Gout")

protVolcanoPlot(covarDF = Type5,
                dzID = "age_i50_first_reported_heart_failure_f131354_0_0",
                R2cutoff = 0.02,
                title = "Proteome Wide Analysis: Heart Failure")

#"age_m10_first_reported_gout_f131858_0_0"
#"age_i70_first_reported_atherosclerosis_f131380_0_0"

}


### MAKE FUNCTION FOR BELOW!!!!
Type5 %>%
  filter(ID == "PRSS8") %>%
  slice_max(HR, n = 10)

Type5 %>%
  filter(ID == "CXCL17") %>%
  slice_max(HR, n = 10)

Type5 %>%
  filter(ID == "ALPP") %>%
  slice_max(HR, n = 10)
# PRSS8, CXCL17, ALPP, PIGR - Smoking
# PULMONARY diseases^^^

Type5 %>%
  filter(ID == "MAMDC4") %>%
  slice_max(HR, n = 10)

Type5 %>%
  filter(ID == "APOF") %>%
  slice_max(HR, n = 10)

Type5 %>%
  filter(ID == "SERPINA7") %>%
  slice_max(HR, n = 10)


#'*BOOTSTRAP VERSION HERE!!!!!*
###Obtain bootstrapped estimate for protein-disease connections:
#'*CHECK Bootstrapping Methods~~*
source("/n/groups/patel/shakson_ukb/UK_Biobank/RScripts/Pure_StatGen/Prot_ExPGS/runProt_PGS_PXS_mediation_bootstrap_source.R")

# Function to plot mediation results w/ Bootstrap.
mediation_plot <- function(MDres, title){
  bootstrap_res_upd <- MDres %>%
    mutate(E_HRi_orig = exp(`Exposure Indirect Effect`),
           E_HRi_2.5 = exp(`Exposure Indirect Effect q2.5`),
           E_HRi_97.5 = exp(`Exposure Indirect Effect q97.5`),
           G_HRi_orig = exp(`Genetic Indirect Effect`),
           G_HRi_2.5 = exp(`Genetic Indirect Effect q2.5`),
           G_HRi_97.5 = exp(`Genetic Indirect Effect q97.5`)
    ) %>% 
    select(all_of(c("ID","DZid",
                    "E_HRi_orig","E_HRi_2.5","E_HRi_97.5",
                    "G_HRi_orig","G_HRi_2.5","G_HRi_97.5"))) %>%
    pivot_longer(
      cols = -c(ID,DZid),
      names_to = c("Type",".value"),
      names_pattern = "(E|G)_HRi_(.*)"
    )
  
  gg1 <- bootstrap_res_upd  %>%
    ggplot(aes(x = ID, y = orig, color = Type)) +
    geom_point(position = position_dodge(width = 0.6), size = 3) +  # Point estimate
    geom_errorbar(aes(ymin = `2.5`, ymax = `97.5`), width = 0.2,    # Error bars
                  position = position_dodge(width = 0.6))  +
    geom_hline(yintercept = 1, color = "black", size = 0.8, linetype = "dashed") +  # Reference line at 1
    #scale_color_manual(values = c("Group1" = "blue", "Group2" = "red")) +
    scale_color_discrete(labels = c("E", "G")) +
    labs(x = "Proteins", y = "HR: Indirect Effects") +
    coord_flip() +
    theme_minimal() +
    ggtitle(title)
  
  return(gg1)
}
mediation_plot_save <- function(MDres, title){
  bootstrap_res_upd <- MDres %>%
    mutate(E_HRi_orig = exp(`Exposure Indirect Effect`),
           E_HRi_2.5 = exp(`Exposure Indirect Effect q2.5`),
           E_HRi_97.5 = exp(`Exposure Indirect Effect q97.5`),
           G_HRi_orig = exp(`Genetic Indirect Effect`),
           G_HRi_2.5 = exp(`Genetic Indirect Effect q2.5`),
           G_HRi_97.5 = exp(`Genetic Indirect Effect q97.5`)
    ) %>% 
    select(all_of(c("ID","DZid",
                    "E_HRi_orig","E_HRi_2.5","E_HRi_97.5",
                    "G_HRi_orig","G_HRi_2.5","G_HRi_97.5"))) %>%
    pivot_longer(
      cols = -c(ID,DZid),
      names_to = c("Type",".value"),
      names_pattern = "(E|G)_HRi_(.*)"
    )
  
  gg1 <- bootstrap_res_upd  %>%
    ggplot(aes(x = ID, y = orig, color = Type)) +
    geom_point(position = position_dodge(width = 0.6), size = 1.5) +  # Point estimate
    geom_errorbar(aes(ymin = `2.5`, ymax = `97.5`), width = 0.2,    # Error bars
                  position = position_dodge(width = 0.6))  +
    geom_hline(yintercept = 1, color = "black", size = 0.8, linetype = "dashed") +  # Reference line at 1
    #scale_color_manual(values = c("Group1" = "blue", "Group2" = "red")) +
    scale_color_discrete(labels = c("E", "G")) +
    labs(x = "Proteins", y = "HR: Indirect Effects") +
    coord_flip() +
    theme_minimal() +
    ggtitle(title)
  
  lapply(c("png", "svg", "pdf"), 
         function(fmt) ggsave(paste0("./Figures/HEAP/F4/", "Mediation",title,".", fmt), 
                              plot = gg1, width = 5, height = 4, dpi = 500))
  
}


#Find top hits for specific diseases:
# Mental Health, Neurological Disorders
# MAMDC4 - age_f10_first_reported_mental_and_behavioural_disorders_due_to_use_of_alcohol_f130854_0_0
#CEACAM16 -  age_f10_first_reported_mental_and_behavioural_disorders_due_to_use_of_alcohol_f130854_0_0
# PRSS8, CXCL17, ALPP, PIGR - Smoking
# PULMONARY diseases^^^
# FOLR1 - Vitamins
# age_l60_first_reported_nail_disorders_f131774_0_0
# ADM - Exercise == Cardiac Arrest

#IDEA:
AlcoholAbuse <- Type5 %>%
  filter(DZid %in% c("age_f10_first_reported_mental_and_behavioural_disorders_due_to_use_of_alcohol_f130854_0_0"))

#10 minutes runtime:

#Start: 10:28 am
#For 100 bootstraps estimated runtime = 5 minutes for each protein.
#End:  1 hour 30 minutes

# Mental Disorders - Alcohol Use: (10 minutes)
MentalDisorders_mediation <-  runBoot_MDwrap(covarType = "Type5", 
                                             protlist = c("MAMDC4","CEACAM16"), 
                                             dzlist = c("age_f10_first_reported_mental_and_behavioural_disorders_due_to_use_of_alcohol_f130854_0_0"), 
                                             100)
mediation_plot(MentalDisorders_mediation, title = "F10: Mental Disorders - Alcohol Use")
mediation_plot_save(MentalDisorders_mediation, title = "F10: Mental Disorders - Alcohol Use")

# Emphysema: (28 minutes)
Emphysema_mediation <-  runBoot_MDwrap(covarType = "Type5", 
                                       protlist = c("CXCL17","LAMP3","WFDC2","TNR","GDF15"), 
                                       dzlist = c("age_j43_first_reported_emphysema_f131490_0_0"), 
                                       100)
mediation_plot(Emphysema_mediation, title = "J43: Emphysema")
mediation_plot_save(Emphysema_mediation, title = "J43: Emphysema")

# Type 2 Diabetes: (43 minutes)
T2D_mediation <-  runBoot_MDwrap(covarType = "Type5", 
                                 protlist = c("IGFBP1","IGFBP2","CKB","LPL","FABP4","FGF21","LEP","IL18R1","IGSF9"), 
                                 dzlist = c("age_e11_first_reported_non_insulin_dependent_diabetes_mellitus_f130708_0_0"), 
                                 100)
mediation_plot(T2D_mediation, title = "E11: Type 2 Diabetes")
mediation_plot_save(T2D_mediation, title = "E11: Type 2 Diabetes")


### 
# Index: 10, 100, 1000, 100000
# Run it once (10000)
# 









######### TRASH BELOW -- SCRIPTS TO PULL FROM #######

# Obtain Deming Slope across Specifications:
#library(mcr) 
# dSlope <- function(MDres){
#   # Obtain Expected Indirect Effect Slope: log(E/G)
#   protContext <- pblapply(unique(MDres$ID), function(p){
#     df <- MDres %>% filter(ID == p) 
#     
#     df <- df %>% 
#               select(c(`Genetic Indirect Effect`,
#                        `Exposure Indirect Effect`)) %>%
#               drop_na()
#     # Perform Deming regression
#     deming_result <- mcreg(df$`Genetic Indirect Effect`, df$`Exposure Indirect Effect`, 
#                            method.reg = "Deming", error.ratio = 1)
# 
#     # Extract slope and intercept
#     #coef(deming_result)["Intercept","EST"]
#     #coef(deming_result)["Slope","EST"]
#     #coef(deming_result)["Slope","LCI"]
#     #coef(deming_result)["Slope","UCI"]
# 
#     maxE <- max(abs(df$`Exposure Indirect Effect`))
#     maxG <- max(abs(df$`Genetic Indirect Effect`))
#     
#     df <- MDres %>%
#       filter(ID == p) %>%
#       filter(!((HR_lower95 < 1 & HR_upper95 > 1) | (HR_lower95 > 1 & HR_upper95 < 1)))
#     
#     num_diseases <- length(unique(df$DZid))
#     
#     protStats <- data.frame(
#       ID = p,
#       Intercept = coef(deming_result)["Intercept","EST"],
#       Estimate = coef(deming_result)["Slope","EST"],
#       LCI = coef(deming_result)["Slope","LCI"],
#       UCI = coef(deming_result)["Slope","UCI"],
#       max.E = maxE,
#       max.G = maxG,
#       Estimate_Norm = coef(deming_result)["Slope","EST"] * (maxE/maxG),
#       NumDiseases = num_diseases
#     )
#     
#     #Below return the slope of the above regression
#     return(protStats) 
#   })
#   protContext <- do.call(rbind, protContext)
#   
#   return(protContext)
# }
#Type5_d <- dSlope(Type5)
#plot(Type5_d$Estimate, Type5_o$Estimate)
#Typeall <- merge(Type1_o, Type5_o, by = "ID")
#plot(Typeall$Estimate.x, Typeall$Estimate.y)


#'*Figure out how to switch to deming regression*
#'*REMOVE condition that significant associations have to be used for the slope...*
protContext <- pblapply(unique(Type5$ID), function(p){

  df <- Type5 %>%
    filter(ID == p) #%>%
  #filter(!((HR_lower95 < 1 & HR_upper95 > 1) | (HR_lower95 > 1 & HR_upper95 < 1)))
  #Second filt is only significant protein - disease connections!
  
  fit <- lm(`Exposure Indirect Effect` ~ `Genetic Indirect Effect`,
            data = df)
  
  stats <- summary(fit)$coefficients
  #as.data.frame(summary(fit))
  
  maxE <- max(abs(df$`Exposure Indirect Effect`))
  #diff(range(df$`Exposure Indirect Effect`))
  maxG <- max(abs(df$`Genetic Indirect Effect`))
  #diff(range(df$`Genetic Indirect Effect`))
  
  #df2 <- df %>%
  #  filter(!((HR_lower95 < 1 & HR_upper95 > 1) | (HR_lower95 > 1 & HR_upper95 < 1)))
  
  df <- Type5 %>%
    filter(ID == p) %>%
    filter(!((HR_lower95 < 1 & HR_upper95 > 1) | (HR_lower95 > 1 & HR_upper95 < 1)))
  
  num_diseases <- length(unique(df$DZid))
  
  protStats <- data.frame(
    ID = p,
    Estimate = stats[2,"Estimate"],
    Std.Error = stats[2,"Std. Error"],
    Pvalue = stats[2,"Pr(>|t|)"],
    max.E = maxE,
    max.G = maxG,
    Estimate_Norm = stats[2,"Estimate"] * (maxE/maxG),
    NumDiseases = num_diseases
  )
  
  #Below return the slope of the above regression
  return(protStats) 
})
#names(protContext) <- unique(Type5$ID)
protContext <- do.call(rbind, protContext)


library(ggrepel)
subset_protContext <- protContext %>%
  filter(Estimate > 1)

protContext %>%
  ggplot(aes(x = NumDiseases, y = log(Estimate))) +
  geom_point() +
  #Color Environment and Genetics Driven Proteins:
  geom_ribbon(aes(ymin =-Inf, ymax = 0, fill = "Genetics Driven"), alpha = 0.2) +
  geom_ribbon(aes(ymin = 0, ymax = Inf, fill = "Environment Driven"), alpha = 0.2) +
  scale_fill_manual(values = c("Genetics Driven" = "blue", "Environment Driven" = "green")) +
  geom_text_repel(data=subset_protContext, 
                  aes(x = NumDiseases, y = log(Estimate), label= ID)) +
  xlab("Number of Significant Disease Hits") +
  ylab("Expected Mixture of Indirect Effect log(E/G)") +
  theme_minimal()

#Plot with blue region:
subset_protContext <- protContext %>%
  filter(Estimate > 1.2 | Estimate < 0.8 |
           NumDiseases > 120)


protContext %>%
  ggplot(aes(x = NumDiseases, y = log(Estimate))) +
  #Color Environment and Genetics Driven Proteins:
  geom_ribbon(aes(ymin =-Inf, ymax = 0, fill = "Genetics Driven"), alpha = 0.2) +
  geom_ribbon(aes(ymin = 0, ymax = Inf, fill = "Environment Driven"), alpha = 0.2) +
  scale_fill_manual(values = c("Genetics Driven" = "blue", "Environment Driven" = "green")) +
  geom_point() +
  geom_text_repel(data=subset_protContext, 
                  aes(x = NumDiseases, y = log(Estimate), 
                      label= ID),
                  size = 3.5) +
  labs(x = "Number of Significant Disease Hits",
       y ="Expected Mixture of Indirect Effect log(E/G)",
       fill = "Region")+
  theme_minimal() +
  theme(legend.position = "bottom")


protContext %>%
  ggplot(aes(x = NumDiseases, y = log(Estimate))) +
  geom_point() +
  stat_smooth(method = "lm", 
              formula = y ~ x, 
              geom = "smooth")




# Obtain Deming Slope across Specifications:








# ANALYSIS and FIGURES --------

# FIND Robust Clusters 


### HELPFUL metrics to help me select:

# Disease selection:
# Number of Significant Associations with Proteins
# Min and Max HR for a given disease
# Min and Max E:HR for given disease
# Min and Max G:HR for given disease

# Protein selection:
# Number of Significant Associations with disease
# HR across diseases
# 




# Proteome-Wide Signal for Disease (from Outcome Model of Mediation Analysis)





# Function: Volcano Plot of Protein Associations to a Given Disease
protVolcanoPlot <- function(covarDF, dzID, R2cutoff, title){
  DZ_res <- covarDF %>%
                filter(DZid %in% dzID) %>%
                filter(R2 > R2cutoff)
  
  prot_labels <- DZ_res %>%
                    filter(HR < 0.75 | HR > 1.5)
  
  # Find the maximum deviation from 1
  max_deviation <- max(abs(DZ_res$HR - 1))
  symmetric_limit <- 1 + max_deviation * c(-1, 1)
  
  
  gg1 <- DZ_res %>%
    ggplot(aes(x = HR, y = cindex)) +
    geom_point() +
    geom_text_repel(data = prot_labels, aes(label = ID)) +
    scale_x_continuous(limits = symmetric_limit)  +
    labs(x = "Hazard Ratio",
         y = "c-index") +
    ggtitle(title) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
  
  return(gg1)          
}

protVolcanoPlot(covarDF = Type5,
                dzID = "age_e11_first_reported_non_insulin_dependent_diabetes_mellitus_f130708_0_0",
                R2cutoff = 0.02,
                title = "Proteome Wide Analysis: Type 2 Diabetes")

protVolcanoPlot(covarDF = Type5,
                dzID = "age_e11_first_reported_non_insulin_dependent_diabetes_mellitus_f130708_0_0",
                R2cutoff = 0.02,
                title = "Proteome Wide Analysis: Type 1 Diabetes")

protVolcanoPlot(covarDF = Type5,
                dzID = "age_f10_first_reported_mental_and_behavioural_disorders_due_to_use_of_alcohol_f130854_0_0",
                R2cutoff = 0.02,
                title = "Proteome Wide Analysis: Mental Disorders Due to Alcohol Usage")


protVolcanoPlot(covarDF = Type5,
                dzID = "age_n18_first_reported_chronic_renal_failure_f132032_0_0",
                R2cutoff = 0.02,
                title = "Proteome Wide Analysis: Chronic Renal Failure")

protVolcanoPlot(covarDF = Type1,
                dzID = "age_e66_first_reported_obesity_f130792_0_0",
                R2cutoff = 0.02,
                title = "Proteome Wide Analysis: Incident Obesity")


protVolcanoPlot(covarDF = Type5,
                dzID = "age_e66_first_reported_obesity_f130792_0_0",
                R2cutoff = 0.02,
                title = "Proteome Wide Analysis: Incident Obesity")

#age_j84_first_reported_other_interstitial_pulmonary_diseases_f131528_0_0
# PRSS8, CXCL17, ALPP

protVolcanoPlot(covarDF = Type5,
                dzID = "age_j84_first_reported_other_interstitial_pulmonary_diseases_f131528_0_0",
                R2cutoff = 0.02,
                title = "Proteome Wide Analysis: Interstitial Pulmonary Diseases")


protVolcanoPlot(covarDF = Type5,
                dzID = "age_j44_first_reported_other_chronic_obstructive_pulmonary_disease_f131492_0_0",
                R2cutoff = 0.02,
                title = "Proteome Wide Analysis: COPD")

protVolcanoPlot(covarDF = Type5,
                dzID = "age_j43_first_reported_emphysema_f131490_0_0",
                R2cutoff = 0.02,
                title = "Proteome Wide Analysis: Emphysema")

protVolcanoPlot(covarDF = Type5,
                dzID = "age_k76_first_reported_other_diseases_of_liver_f131670_0_0",
                R2cutoff = 0.02,
                title = "Proteome Wide Analysis: Liver Disease")


protVolcanoPlot(covarDF = Type5,
                dzID = "age_i50_first_reported_heart_failure_f131354_0_0",
                R2cutoff = 0.02,
                title = "Proteome Wide Analysis: Heart Failure")


protVolcanoPlot(covarDF = Type5,
                dzID = "age_i27_first_reported_other_pulmonary_heart_diseases_f131310_0_0",
                R2cutoff = 0.02,
                title = "Proteome Wide Analysis: Heart Failure")

protVolcanoPlot(covarDF = Type5,
                dzID = "age_i64_first_reported_stroke_not_specified_as_haemorrhage_or_infarction_f131368_0_0",
                R2cutoff = 0.02,
                title = "Proteome Wide Analysis: Stroke")

protVolcanoPlot(covarDF = Type1,
                dzID = "age_g30_first_reported_alzheimers_disease_f131036_0_0",
                R2cutoff = 0.02,
                title = "Proteome Wide Analysis: Alzheimers")

protVolcanoPlot(covarDF = Type1,
                dzID = "age_f03_first_reported_unspecified_dementia_f130842_0_0",
                R2cutoff = 0.02,
                title = "Proteome Wide Analysis: Dementia")

protVolcanoPlot(covarDF = Type1,
                dzID = "age_f00_first_reported_dementia_in_alzheimers_disease_f130836_0_0",
                R2cutoff = 0.02,
                title = "Proteome Wide Analysis: Dementia")

protVolcanoPlot(covarDF = Type5,
                dzID = "age_i70_first_reported_atherosclerosis_f131380_0_0",
                R2cutoff = 0.02,
                title = "Proteome Wide Analysis: Atherosclerosis")



#Suggestions: CKD, Heart Failure, Gout

protVolcanoPlot(covarDF = Type5,
                dzID = "age_n18_first_reported_chronic_renal_failure_f132032_0_0",
                R2cutoff = 0.02,
                title = "Proteome Wide Analysis: Chronic Renal Failure")

protVolcanoPlot(covarDF = Type5,
                dzID = "age_m10_first_reported_gout_f131858_0_0",
                R2cutoff = 0.02,
                title = "Proteome Wide Analysis: Gout")

protVolcanoPlot(covarDF = Type5,
                dzID = "age_i50_first_reported_heart_failure_f131354_0_0",
                R2cutoff = 0.02,
                title = "Proteome Wide Analysis: Heart Failure")

"age_m10_first_reported_gout_f131858_0_0"
"age_i70_first_reported_atherosclerosis_f131380_0_0"

### MAKE FUNCTION FOR BELOW!!!!
Type5 %>%
  filter(ID == "PRSS8") %>%
  slice_max(HR, n = 10)

Type5 %>%
  filter(ID == "CXCL17") %>%
  slice_max(HR, n = 10)

Type5 %>%
  filter(ID == "ALPP") %>%
  slice_max(HR, n = 10)
# PRSS8, CXCL17, ALPP, PIGR - Smoking
# PULMONARY diseases^^^

Type5 %>%
  filter(ID == "MAMDC4") %>%
  slice_max(HR, n = 10)

Type5 %>%
  filter(ID == "SERPINA7") %>%
  slice_max(HR, n = 10)



unique(Type5$DZid)[grepl("renal",unique(Type5$DZid))]
unique(Type5$DZid)[grepl("obesity",unique(Type5$DZid))]
unique(Type5$DZid)[grepl("diabetes",unique(Type5$DZid))]
unique(Type5$DZid)[grepl("liver",unique(Type5$DZid))]
unique(Type5$DZid)[grepl("alcohol",unique(Type5$DZid))]
unique(Type5$DZid)[grepl("heart",unique(Type5$DZid))]
unique(Type5$DZid)[grepl("stroke",unique(Type5$DZid))]
unique(Type5$DZid)[grepl("demen",unique(Type5$DZid))]
unique(Type5$DZid)[grepl("alz",unique(Type5$DZid))]
unique(Type5$DZid)[grepl("athe",unique(Type5$DZid))]
unique(Type5$DZid)[grepl("gout",unique(Type5$DZid))]


T2D <- Type5 %>%
  filter(DZid %in% c("age_e11_first_reported_non_insulin_dependent_diabetes_mellitus_f130708_0_0")) %>%
  filter(R2 > 0.01) %>%
  filter(R2 > 0.02)

T2D_label <- T2D %>%
  filter(HR < 0.75 | HR > 1.5 | cindex > 0.85)

T2D %>%
  ggplot(aes(x = HR, y = cindex)) +
  geom_point() +
  geom_text_repel(data = T2D_label, aes(label = ID)) +
  theme_minimal()










#'*BOOTSTRAP VERSION HERE!!!!!*
###Obtain bootstrapped estimate for protein-disease connections:
#'*CHECK Bootstrapping Methods~~*
source("/n/groups/patel/shakson_ukb/UK_Biobank/RScripts/Pure_StatGen/Prot_ExPGS/runProt_PGS_PXS_mediation_bootstrap_source.R")

# Function to plot mediation results w/ Bootstrap.
mediation_plot <- function(MDres, title){
  bootstrap_res_upd <- MDres %>%
    mutate(E_HRi_orig = exp(`Exposure Indirect Effect`),
           E_HRi_2.5 = exp(`Exposure Indirect Effect q2.5`),
           E_HRi_97.5 = exp(`Exposure Indirect Effect q97.5`),
           G_HRi_orig = exp(`Genetic Indirect Effect`),
           G_HRi_2.5 = exp(`Genetic Indirect Effect q2.5`),
           G_HRi_97.5 = exp(`Genetic Indirect Effect q97.5`)
    ) %>% 
    select(all_of(c("ID","DZid",
                    "E_HRi_orig","E_HRi_2.5","E_HRi_97.5",
                    "G_HRi_orig","G_HRi_2.5","G_HRi_97.5"))) %>%
    pivot_longer(
      cols = -c(ID,DZid),
      names_to = c("Type",".value"),
      names_pattern = "(E|G)_HRi_(.*)"
    )
  
  gg1 <- bootstrap_res_upd  %>%
    ggplot(aes(x = ID, y = orig, color = Type)) +
    geom_point(position = position_dodge(width = 0.6), size = 3) +  # Point estimate
    geom_errorbar(aes(ymin = `2.5`, ymax = `97.5`), width = 0.2,    # Error bars
                  position = position_dodge(width = 0.6))  +
    geom_hline(yintercept = 1, color = "black", size = 0.8, linetype = "dashed") +  # Reference line at 1
    #scale_color_manual(values = c("Group1" = "blue", "Group2" = "red")) +
    scale_color_discrete(labels = c("E", "G")) +
    labs(x = "Proteins", y = "HR: Indirect Effects") +
    coord_flip() +
    theme_minimal() +
    ggtitle(title)
  
  return(gg1)
}

#Find top hits for specific diseases:
# Mental Health, Neurological Disorders
# MAMDC4 - age_f10_first_reported_mental_and_behavioural_disorders_due_to_use_of_alcohol_f130854_0_0
#CEACAM16 -  age_f10_first_reported_mental_and_behavioural_disorders_due_to_use_of_alcohol_f130854_0_0
# PRSS8, CXCL17, ALPP, PIGR - Smoking
# PULMONARY diseases^^^
# FOLR1 - Vitamins
# age_l60_first_reported_nail_disorders_f131774_0_0
# ADM - Exercise == Cardiac Arrest

#IDEA:
AlcoholAbuse <- Type5 %>%
  filter(DZid %in% c("age_f10_first_reported_mental_and_behavioural_disorders_due_to_use_of_alcohol_f130854_0_0"))

#10 minutes runtime:

# Mental Disorders - Alcohol Use:
MentalDisorders_mediation <-  runBoot_MDwrap(covarType = "Type5", 
                                             protlist = c("MAMDC4","CEACAM16"), 
                                             dzlist = c("age_f10_first_reported_mental_and_behavioural_disorders_due_to_use_of_alcohol_f130854_0_0"), 
                                             10)
mediation_plot(MentalDisorders_mediation, title = "F10: Mental Disorders - Alcohol Use")


# Emphysema:
Emphysema_mediation <-  runBoot_MDwrap(covarType = "Type5", 
                                             protlist = c("CXCL17","LAMP3","WFDC2","TNR","GDF15"), 
                                             dzlist = c("age_j43_first_reported_emphysema_f131490_0_0"), 
                                             10)
mediation_plot(Emphysema_mediation, title = "J43: Emphysema")

# Type 2 Diabetes:
T2D_mediation <-  runBoot_MDwrap(covarType = "Type5", 
                                       protlist = c("IGFBP1","IGFBP2","CKB","LPL","FABP4","FGF21","LEP","IL18R1","IGSF9"), 
                                       dzlist = c("age_e11_first_reported_non_insulin_dependent_diabetes_mellitus_f130708_0_0"), 
                                       10)
mediation_plot(T2D_mediation, title = "J43: Type 2 Diabetes")




AlcoholAbuse %>%
  ggplot(aes(x = HR, y = cindex)) +
  geom_point() +
  geom_text_repel(aes(label = ID)) +
  theme_minimal()



T2D <- Type5 %>%
  filter(DZid %in% c("age_e11_first_reported_non_insulin_dependent_diabetes_mellitus_f130708_0_0")) %>%
  filter(R2 > 0.01) %>%
  filter(R2 > 0.02)

T2D_label <- T2D %>%
                filter(HR < 0.75 | HR > 1.5 | cindex > 0.85)

T2D %>%
  ggplot(aes(x = HR, y = cindex)) +
  geom_point() +
  geom_text_repel(data = T2D_label, aes(label = ID)) +
  theme_minimal()



AlcoholAbuse <- Type5 %>%
  filter(DZid %in% c("age_f10_first_reported_mental_and_behavioural_disorders_due_to_use_of_alcohol_f130854_0_0")) %>%
  filter(R2 > 0.02)

AlcoholAbuse_label <- AlcoholAbuse %>%
  filter(HR < 0.75 | HR > 1.5 | cindex > 0.85)

AlcoholAbuse %>%
  ggplot(aes(x = HR, y = cindex)) +
  geom_point() +
  geom_text_repel(data = AlcoholAbuse_label, aes(label = ID)) +
  theme_minimal()

#Example:
covarType = "Type1"
protlist <- c("IGFBP1","IGFBP2","CKB","FABP4","FGF21","LEP","GDF15")
dzlist <- c("age_e11_first_reported_non_insulin_dependent_diabetes_mellitus_f130708_0_0",
            "age_n18_first_reported_chronic_renal_failure_f132032_0_0",
            "age_i25_first_reported_chronic_ischaemic_heart_disease_f131306_0_0")
dzNames <- c("T2D","CRF","CIHD")
n <- 10

try2 <- runBoot_MDwrap(covarType, protlist, dzlist, 10)
View(try2)


T2D <- Type5 %>%
  filter(DZid == "age_e11_first_reported_non_insulin_dependent_diabetes_mellitus_f130708_0_0")

CKD <- Type5 %>%
  filter(DZid == "age_n18_first_reported_chronic_renal_failure_f132032_0_0")


T2D_e <- T2D %>%
  filter(E_HRi/G_HRi > 1)

T2D_g <- T2D %>%
  filter(E_HRi/G_HRi < 1)


slice_max(abs(E_HRi), n = 10)












#'*Figure out how to switch to deming regression*
#'*REMOVE condition that significant associations have to be used for the slope...*
protContext <- pblapply(unique(Type5$ID), function(p){
  #p = "LEP"
  #p = "SCGB2A2"
  
  df <- Type5 %>%
    filter(ID == p) #%>%
  #filter(!((HR_lower95 < 1 & HR_upper95 > 1) | (HR_lower95 > 1 & HR_upper95 < 1)))
  #Second filt is only significant protein - disease connections!
  
  fit <- lm(`Exposure Indirect Effect` ~ `Genetic Indirect Effect`,
            data = df)
  
  stats <- summary(fit)$coefficients
  #as.data.frame(summary(fit))
  
  maxE <- max(abs(df$`Exposure Indirect Effect`))
  #diff(range(df$`Exposure Indirect Effect`))
  maxG <- max(abs(df$`Genetic Indirect Effect`))
  #diff(range(df$`Genetic Indirect Effect`))
  
  #df2 <- df %>%
  #  filter(!((HR_lower95 < 1 & HR_upper95 > 1) | (HR_lower95 > 1 & HR_upper95 < 1)))
  
  df <- Type5 %>%
    filter(ID == p) %>%
    filter(!((HR_lower95 < 1 & HR_upper95 > 1) | (HR_lower95 > 1 & HR_upper95 < 1)))
  
  num_diseases <- length(unique(df$DZid))
  
  protStats <- data.frame(
    ID = p,
    Estimate = stats[2,"Estimate"],
    Std.Error = stats[2,"Std. Error"],
    Pvalue = stats[2,"Pr(>|t|)"],
    max.E = maxE,
    max.G = maxG,
    Estimate_Norm = stats[2,"Estimate"] * (maxE/maxG),
    NumDiseases = num_diseases
  )
  
  #Below return the slope of the above regression
  return(protStats) 
})
#names(protContext) <- unique(Type5$ID)
protContext <- do.call(rbind, protContext)


library(ggrepel)
subset_protContext <- protContext %>%
  filter(Estimate > 1)

protContext %>%
  ggplot(aes(x = NumDiseases, y = log(Estimate))) +
  geom_point() +
  #Color Environment and Genetics Driven Proteins:
  geom_ribbon(aes(ymin =-Inf, ymax = 0, fill = "Genetics Driven"), alpha = 0.2) +
  geom_ribbon(aes(ymin = 0, ymax = Inf, fill = "Environment Driven"), alpha = 0.2) +
  scale_fill_manual(values = c("Genetics Driven" = "blue", "Environment Driven" = "green")) +
  geom_text_repel(data=subset_protContext, 
                  aes(x = NumDiseases, y = log(Estimate), label= ID)) +
  xlab("Number of Significant Disease Hits") +
  ylab("Expected Mixture of Indirect Effect log(E/G)") +
  theme_minimal()

#Plot with blue region:
subset_protContext <- protContext %>%
  filter(Estimate > 1.2 | Estimate < 0.8 |
           NumDiseases > 120)


protContext %>%
  ggplot(aes(x = NumDiseases, y = log(Estimate))) +
  #Color Environment and Genetics Driven Proteins:
  geom_ribbon(aes(ymin =-Inf, ymax = 0, fill = "Genetics Driven"), alpha = 0.2) +
  geom_ribbon(aes(ymin = 0, ymax = Inf, fill = "Environment Driven"), alpha = 0.2) +
  scale_fill_manual(values = c("Genetics Driven" = "blue", "Environment Driven" = "green")) +
  geom_point() +
  geom_text_repel(data=subset_protContext, 
                  aes(x = NumDiseases, y = log(Estimate), 
                      label= ID),
                  size = 3.5) +
  labs(x = "Number of Significant Disease Hits",
      y ="Expected Mixture of Indirect Effect log(E/G)",
      fill = "Region")+
  theme_minimal() +
  theme(legend.position = "bottom")


protContext %>%
  ggplot(aes(x = NumDiseases, y = log(Estimate))) +
  geom_point() +
  stat_smooth(method = "lm", 
              formula = y ~ x, 
              geom = "smooth")
  

#Count number of E>G mediated effects
Edriven_mediation <- protContext %>%
  filter(Estimate > 1)

#Count number of diseases that were hit:
Disease_Hits_general <- Type5 %>%
  filter(!((HR_lower95 < 1 & HR_upper95 > 1) | (HR_lower95 > 1 & HR_upper95 < 1))) %>%
  count(DZid)

#Count number of diseases hit (with E>G proteins)
Disease_Hits_E <- Type5 %>%
  filter(!((HR_lower95 < 1 & HR_upper95 > 1) | (HR_lower95 > 1 & HR_upper95 < 1))) %>%
  filter(ID %in% Edriven_mediation$ID) %>%
  count(DZid)

subset_protContextv2 <- protContext %>%
  filter(NumDiseases > 100)

View(Type5)


#### GET Category Labels for Disease:
projID = 52887
UKBdict <- fread(file = paste0("/n/groups/patel/shakson_ukb/UK_Biobank/Data/Paths/",projID,"/allpaths.txt")) #'*UKBdict should be a global variable now!!*

# ids <- unique(sub("age","date",Type5$DZid)) #sub replaces first instance ONLY
# 
# ids_cat <- UKBdict[match(ids, UKBdict$descriptive_colnames), ]$Path
# ids_cat <- gsub(".*>\\s*", "", ids_cat)
# 
# DZmatch <- data.frame(
#   DZid = ids,
#   DZcat = ids_cat
# )

IDs <- sub("age","date",Type5$DZid)
Type5$DZid_cat <- UKBdict[match(IDs, UKBdict$descriptive_colnames), ]$Path
Type5$DZid_cat <- gsub(".*>\\s*", "", Type5$DZid_cat)

#### NEED Prot - Disease Mapping:

Type5_sig <- Type5 %>%
              filter(!((HR_lower95 < 1 & HR_upper95 > 1) | (HR_lower95 > 1 & HR_upper95 < 1))) %>%
              select(c(HR, ID, DZid, DZid_cat))

Type5_sigsel <- Type5_sig #%>%
                  #filter(HR < 0.5 | HR > 2)

Type5_context <- merge(Type5_sigsel, protContext, by = "ID")
#^^^^ USE THIS FOR SELECTION!!



Type5_context %>%
  ggplot(aes(x = DZid_cat, y = log(Estimate))) +
  #Color Environment and Genetics Driven Proteins:
  scale_fill_manual(values = c("Genetics Driven" = "blue", "Environment Driven" = "green")) +
  geom_point() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))


Type5_context %>%
  ggplot(aes(x = reorder(DZid,-log(Estimate)),
             y = log(Estimate))) +
  #Color Environment and Genetics Driven Proteins:
  scale_fill_manual(values = c("Genetics Driven" = "blue", "Environment Driven" = "green")) +
  geom_point() +
  theme_minimal() +
  theme(axis.text.x=element_blank(), 
        axis.ticks.x=element_blank())



#CONTEXT Proteins:
# MAMDC4, CEACAM16, SERPINA7 - Alcohol


# Mental Health, Neurological Disorders
# MAMDC4 - age_f10_first_reported_mental_and_behavioural_disorders_due_to_use_of_alcohol_f130854_0_0
#CEACAM16 -  age_f10_first_reported_mental_and_behavioural_disorders_due_to_use_of_alcohol_f130854_0_0

# PRSS8, CXCL17, ALPP, PIGR - Smoking
# PULMONARY diseases^^^

# FOLR1 - Vitamins
# age_l60_first_reported_nail_disorders_f131774_0_0

# ADM - Exercise == Cardiac Arrest
#








####### Potential FIGURES code to PULL from ######

#'*CHECK Bootstrapping Methods~~*
source("/n/groups/patel/shakson_ukb/UK_Biobank/RScripts/Pure_StatGen/Prot_ExPGS/runProt_PGS_PXS_mediation_bootstrap_source.R")

#Example:
covarType = "Type1"
protlist <- c("IGFBP1","IGFBP2","CKB","FABP4","FGF21","LEP","GDF15")
dzlist <- c("age_e11_first_reported_non_insulin_dependent_diabetes_mellitus_f130708_0_0",
            "age_n18_first_reported_chronic_renal_failure_f132032_0_0",
            "age_i25_first_reported_chronic_ischaemic_heart_disease_f131306_0_0")
dzNames <- c("T2D","CRF","CIHD")
n <- 10

try2 <- runBoot_MDwrap(covarType, protlist, dzlist, 10)
View(try2)


T2D <- Type5 %>%
          filter(DZid == "age_e11_first_reported_non_insulin_dependent_diabetes_mellitus_f130708_0_0")

CKD <- Type5 %>%
  filter(DZid == "age_n18_first_reported_chronic_renal_failure_f132032_0_0")


T2D_e <- T2D %>%
  filter(E_HRi/G_HRi > 1)

T2D_g <- T2D %>%
  filter(E_HRi/G_HRi < 1)

  
  slice_max(abs(E_HRi), n = 10)



### Get whether the proteins are E or G controlled across disease. 


### Check stability of assumptions: All associations:

# Is Exposure Indirect Effect independent of Genetic Indirect Effect? (Unless massive confounding)
# Obviously we do not see this and it is reduces as we adjust for potential confounding
# Other check: Direct Effects of G and E are independent so correlations should be close to 0.
fit <- lm(`Exposure Indirect Effect` ~ `Genetic Indirect Effect` , data = Type1)
summary(fit)

fit <- lm(`Exposure Indirect Effect` ~ `Genetic Indirect Effect`, data = Type2)
summary(fit)

fit <- lm(`Exposure Indirect Effect` ~ `Genetic Indirect Effect`, data = Type3)
summary(fit)

fit <- lm(`Exposure Indirect Effect` ~ `Genetic Indirect Effect`, data = Type4)
summary(fit)

fit <- lm(`Exposure Indirect Effect` ~ `Genetic Indirect Effect` , data = Type5)
summary(fit)

MD_cor <- function(Cspec){
  # Cspec is the association table from different covariate specifications:
  Cspec <- Cspec %>%
    filter_all(all_vars(!is.infinite(.))) 
  
  corDir <- cor.test(Cspec$`Genetic Direct Effect`, Cspec$`Exposure Direct Effect`)$estimate
  corIndir <- cor.test(Cspec$`Genetic Indirect Effect`, Cspec$`Exposure Indirect Effect`)$estimate
  
  print(corDir)
  print(corIndir)
}
MD_cor(Type1)
MD_cor(Type2)
MD_cor(Type3)
MD_cor(Type4)
MD_cor(Type5)

# Next questions: 
# Slope of E/G to understand these correlated effects across various proteins.
# Understand if certain disease have greater proportion of E/G proteins?

#Find Strong Assocations:
Type5_DZ <- Type5 %>%
  filter(!((HR_lower95 < 1 & HR_upper95 > 1) | (HR_lower95 > 1 & HR_upper95 < 1))) %>%
  filter(HR > 1.5 | HR < 0.75)

DZcount <- Type5_DZ %>%
              count(DZid)


#'*Figure out how to switch to deming regression*
#'*REMOVE condition that significant associations have to be used for the slope...*
protContext <- pblapply(unique(Type5$ID), function(p){
  #p = "LEP"
  #p = "SCGB2A2"
  
  df <- Type5 %>%
    filter(ID == p) #%>%
    #filter(!((HR_lower95 < 1 & HR_upper95 > 1) | (HR_lower95 > 1 & HR_upper95 < 1)))
  #Second filt is only significant protein - disease connections!
  
  fit <- lm(`Exposure Indirect Effect` ~ `Genetic Indirect Effect`,
            data = df)
  
  stats <- summary(fit)$coefficients
  #as.data.frame(summary(fit))
  
  maxE <- max(abs(df$`Exposure Indirect Effect`))
  #diff(range(df$`Exposure Indirect Effect`))
  maxG <- max(abs(df$`Genetic Indirect Effect`))
  #diff(range(df$`Genetic Indirect Effect`))
  
  #df2 <- df %>%
  #  filter(!((HR_lower95 < 1 & HR_upper95 > 1) | (HR_lower95 > 1 & HR_upper95 < 1)))
  
  df <- Type5 %>%
    filter(ID == p) %>%
    filter(!((HR_lower95 < 1 & HR_upper95 > 1) | (HR_lower95 > 1 & HR_upper95 < 1)))
  
  num_diseases <- length(unique(df$DZid))
  
  protStats <- data.frame(
    ID = p,
    Estimate = stats[2,"Estimate"],
    Std.Error = stats[2,"Std. Error"],
    Pvalue = stats[2,"Pr(>|t|)"],
    max.E = maxE,
    max.G = maxG,
    Estimate_Norm = stats[2,"Estimate"] * (maxE/maxG),
    NumDiseases = num_diseases
  )
  
  #Below return the slope of the above regression
  return(protStats) 
})
#names(protContext) <- unique(Type5$ID)
protContext <- do.call(rbind, protContext)


Eprot <- protContext %>%
  filter(Estimate > 1) %>%
  pull(ID)

Gprot <- protContext %>%
  filter(Estimate < 1) %>%
  pull(ID)

Type5_DZe <- Type5_DZ %>%
  filter(ID %in% Eprot) %>%
  count(DZid)

Type5_DZe_type <- Type5_DZ %>%
  filter(ID %in% Eprot) %>%
  count(ID)

Type5_DZg_type <- Type5_DZ %>%
  filter(ID %in% Gprot) %>%
  count(ID)

any(Type5$ID %in% Eprot)

Type5_DZg <- Type5_DZ %>%
  filter(ID %in% Gprot) %>%
  count(DZid)

Type5_counts <- merge(Type5_DZe, DZcount, by = "DZid")
Type5_counts$Erat <- Type5_counts$n.x/Type5_counts$n.y



#Examples to show:
Type5_DZe_type1 <- Type5_DZ %>%
  filter(ID %in% Eprot)
#CXCL17 protein involved with smoking. Has E>G. What diseases does this protein
#impact? 4 diseases emphysema, bronchietasis, interstitial pulmonary disease, alzheimers

#LEP protein: what diseases
#GDF15 protein:


#Weird proteins:
#SCGB2A2, CYB5A, NUDT15


p = "SCGB2A2"
p = "CYB5A"
p = "LEP"
p = "GDF15"
p = "IL1RN"
p = "TJP3"

#Only significant disease hits with proteins
df <- Type5 %>%
  filter(ID == p) %>%
  filter(!((HR_lower95 < 1 & HR_upper95 > 1) | (HR_lower95 > 1 & HR_upper95 < 1)))
#Second filt is only significant protein - disease connections!



#Call all diseases for a protein
#'*This one is more stable then only significant hits*
#'Should do All Diseases --> Calculate the corresponding statistics below
#'Then report number of significant disease hits for a protein using HR and CI
#'

p = "LEP"
df <- Type5 %>%
  filter(ID == p) 

fit <- lm(`Exposure Indirect Effect` ~ `Genetic Indirect Effect`,
          data = df)
summary(fit)

df %>%
  ggplot(aes(x = `Genetic Indirect Effect`, y = `Exposure Indirect Effect`)) +
  geom_point() +
  stat_smooth(method = "lm", 
              formula = y ~ x, 
              geom = "smooth") 

library(mcr)  
# Perform Deming regression
deming_result <- mcreg(df$`Genetic Indirect Effect`, df$`Exposure Indirect Effect`, 
                       method.reg = "Deming", error.ratio = 1)
printSummary(deming_result)

# Extract slope and intercept
intercept <- coef(deming_result)["Intercept","EST"]
slope <- coef(deming_result)["Slope","EST"]
coef(deming_result)["Slope","LCI"]
coef(deming_result)["Slope","UCI"]



# Create ggplot
df %>%
  ggplot(aes(x = `Genetic Indirect Effect`, y = `Exposure Indirect Effect`)) +
  geom_point() +  # Scatter plot
  geom_abline(intercept = intercept, slope = slope, color = "blue", linetype = "dashed") +
  labs(title = "Deming Regression in ggplot",
       x = "X-axis",
       y = "Y-axis") +
  theme_minimal()



####Direct Effects Differences:
df %>%
  ggplot(aes(x = `Genetic Direct Effect`, y = `Exposure Direct Effect`)) +
  geom_point() +
  stat_smooth(method = "lm", 
              formula = y ~ x, 
              geom = "smooth") 
cor.test(df$`Genetic Direct Effect`,df$`Exposure Direct Effect`)

df %>%
  ggplot(aes(x = `Genetic Direct Effect`, y = `Genetic Indirect Effect`)) +
  geom_point() +
  stat_smooth(method = "lm", 
              formula = y ~ x, 
              geom = "smooth") 


df %>%
  ggplot(aes(x = `Exposure Direct Effect`, y = `Exposure Indirect Effect`)) +
  geom_point() +
  stat_smooth(method = "lm", 
              formula = y ~ x, 
              geom = "smooth") 

df %>%
  ggplot(aes(x = `Genetic Direct Effect`, y = `Exposure Indirect Effect`)) +
  geom_point() +
  stat_smooth(method = "lm", 
              formula = y ~ x, 
              geom = "smooth") 

df %>%
  ggplot(aes(x = `Exposure Direct Effect`, y = `Genetic Indirect Effect`)) +
  geom_point() +
  stat_smooth(method = "lm", 
              formula = y ~ x, 
              geom = "smooth") 


#### Developing a statistic:
stat <- Type5$`Exposure Indirect Effect` - Type5$`Genetic Indirect Effect`

Type5 %>%
  ggplot(aes(sample=`Exposure Indirect Effect`)) +
  stat_qq() + 
  stat_qq_line()


### Variance Decomposition of Indirect Effects:
library(lme4)

Type5_stats <- Type5

# Fit the mixed-effects model
#model <- lmer(`Genetic Indirect Effect` ~ 1 + (1|ID) + (1|DZid), data = Type5_stats)
model <- lmer(`Exposure Indirect Effect` ~ 1 + (1|ID) + (1|DZid), data = Type5_stats)

var_comp <- VarCorr(model)
print(var_comp)

# Extract variances
#Mediator = ID
#Disease = DZid
var_mediator <- attr(var_comp$ID, "stddev")^2
var_disease <- attr(var_comp$DZid, "stddev")^2
var_residual <- attr(var_comp, "sc")^2

# Total variance
total_variance <- var_mediator + var_disease + var_residual

# Proportion of variance
prop_mediator <- var_mediator / total_variance
prop_disease <- var_disease / total_variance
prop_residual <- var_residual / total_variance

cat("Proportion of variance by Mediator:", prop_mediator, "\n")
cat("Proportion of variance by Disease:", prop_disease, "\n")
cat("Proportion of variance Residual:", prop_residual, "\n")


### TO DO Next





#'*Potetial 1st Mediation Analysis Figure:*

#Biologically: 
#Protein function should be similar across contexts (disease)
#But these proteins may be predictive of future disease risk
#'*Figure out how to switch to deming regression*
protContext <- pblapply(unique(Type5$ID), function(p){
  #p = "LEP"
  #p = "SCGB2A2"
  
  df <- Type5 %>%
    filter(ID == p) %>%
    filter(!((HR_lower95 < 1 & HR_upper95 > 1) | (HR_lower95 > 1 & HR_upper95 < 1)))
  #Second filt is only significant protein - disease connections!
  
  fit <- lm(`Exposure Indirect Effect` ~ `Genetic Indirect Effect`,
     data = df)
  
  stats <- summary(fit)$coefficients
    #as.data.frame(summary(fit))
  
  maxE <- max(abs(df$`Exposure Indirect Effect`))
    #diff(range(df$`Exposure Indirect Effect`))
  maxG <- max(abs(df$`Genetic Indirect Effect`))
    #diff(range(df$`Genetic Indirect Effect`))
  
  #df2 <- df %>%
  #  filter(!((HR_lower95 < 1 & HR_upper95 > 1) | (HR_lower95 > 1 & HR_upper95 < 1)))
  num_diseases <- length(unique(df$DZid))
  
  protStats <- data.frame(
    ID = p,
    Estimate = stats[2,"Estimate"],
    Std.Error = stats[2,"Std. Error"],
    Pvalue = stats[2,"Pr(>|t|)"],
    max.E = maxE,
    max.G = maxG,
    Estimate_Norm = stats[2,"Estimate"] * (maxE/maxG),
    NumDiseases = num_diseases
  )
  
  #Below return the slope of the above regression
  return(protStats) 
})
#names(protContext) <- unique(Type5$ID)
protContext <- do.call(rbind, protContext)

#Plot protContext
protContext %>%
  ggplot(aes(x = Estimate, y = max.G)) +
  xlab("Slope of Indirect Effects: Exposure/Genetics") +
  geom_point()

protContext %>%
  ggplot(aes(x = NumDiseases, y = Estimate)) +
  #xlab("Slope of Indirect Effects: Exposure/Genetics") +
  geom_point()

protContext %>%
  filter(Estimate > 1) %>%
  ggplot(aes(x = NumDiseases, y = Estimate)) +
  #xlab("Slope of Indirect Effects: Exposure/Genetics") +
  geom_point()


library(ggrepel)
subset_protContext <- protContext %>%
                        filter(Estimate > 1)

protContext %>%
  ggplot(aes(x = NumDiseases, y = log(Estimate))) +
  geom_point() +
  #Color Environment and Genetics Driven Proteins:
  geom_ribbon(aes(ymin =-Inf, ymax = 0, fill = "Genetics Driven"), alpha = 0.2) +
  geom_ribbon(aes(ymin = 0, ymax = Inf, fill = "Environment Driven"), alpha = 0.2) +
  scale_fill_manual(values = c("Genetics Driven" = "blue", "Environment Driven" = "green")) +
  geom_text_repel(data=subset_protContext, 
                  aes(x = NumDiseases, y = log(Estimate), label= ID)) +
  xlab("Number of Significant Disease Hits") +
  ylab("Expected Mixture of Indirect Effect log(E/G)") +
  theme_minimal()

#Plot with blue region:
subset_protContext <- protContext %>%
  filter(Estimate > 1.2 | Estimate < 0.8)

protContext %>%
  ggplot(aes(x = NumDiseases, y = log(Estimate))) +
  geom_point() +
  #Color Environment and Genetics Driven Proteins:
  geom_ribbon(aes(ymin =-Inf, ymax = 0, fill = "Genetics Driven"), alpha = 0.2) +
  geom_ribbon(aes(ymin = 0, ymax = Inf, fill = "Environment Driven"), alpha = 0.2) +
  scale_fill_manual(values = c("Genetics Driven" = "blue", "Environment Driven" = "green")) +
  geom_text_repel(data=subset_protContext, 
                  aes(x = NumDiseases, y = log(Estimate), label= ID)) +
  xlab("Number of Significant Disease Hits") +
  ylab("Expected Mixture of Indirect Effect log(E/G)") +
  theme_minimal()

#Count number of E>G mediated effects
Edriven_mediation <- protContext %>%
  filter(Estimate > 1)

#Count number of diseases that were hit:
Disease_Hits_general <- Type5 %>%
  filter(!((HR_lower95 < 1 & HR_upper95 > 1) | (HR_lower95 > 1 & HR_upper95 < 1))) %>%
  count(DZid)

#Count number of diseases hit (with E>G proteins)
Disease_Hits_E <- Type5 %>%
  filter(!((HR_lower95 < 1 & HR_upper95 > 1) | (HR_lower95 > 1 & HR_upper95 < 1))) %>%
  filter(ID %in% Edriven_mediation$ID) %>%
  count(DZid)

subset_protContextv2 <- protContext %>%
  filter(NumDiseases > 100)

##### IDEA #####
#Change x axis to different aspects about the protein
#Example: Protein more secreted? Protein more involved with Cardiometabolic Traits?
#etc.




protContext %>%
  ggplot(aes(x = NumDiseases, y = log(Estimate))) +
  geom_point() +
  #Color Environment and Genetics Driven Proteins:
  geom_ribbon(aes(ymin =-Inf, ymax = 0, fill = "Genetics Driven"), alpha = 0.2) +
  geom_ribbon(aes(ymin = 0, ymax = Inf, fill = "Environment Driven"), alpha = 0.2) +
  scale_fill_manual(values = c("Genetics Driven" = "blue", "Environment Driven" = "green")) +
  geom_text_repel(data=subset_protContextv2, 
                  aes(x = NumDiseases, y = log(Estimate), label= ID)) +
  xlab("Number of Significant Disease Hits") +
  ylab("Expected Mixture of Indirect Effect log(E/G)") +
  theme_minimal()

###Statistic to Record:
length(unique(Type5$DZid))
median(protContext$NumDiseases)
min(protContext$NumDiseases)



protContext %>%
  ggplot(aes(x = NumDiseases, y = max.E)) +
  geom_point() +
  #Color Environment and Genetics Driven Proteins:
  geom_ribbon(aes(ymin =-Inf, ymax = 0, fill = "Genetics Driven"), alpha = 0.2) +
  geom_ribbon(aes(ymin = 0, ymax = Inf, fill = "Environment Driven"), alpha = 0.2) +
  scale_fill_manual(values = c("Genetics Driven" = "blue", "Environment Driven" = "green")) +
  xlab("Number of Significant Disease Hits") +
  ylab("Expected Mixture of Indirect Effect log(E/G)") +
  theme_minimal()

protContext %>%
  ggplot(aes(x = NumDiseases, y = max.G)) +
  geom_point() +
  #Color Environment and Genetics Driven Proteins:
  geom_ribbon(aes(ymin =-Inf, ymax = 0, fill = "Genetics Driven"), alpha = 0.2) +
  geom_ribbon(aes(ymin = 0, ymax = Inf, fill = "Environment Driven"), alpha = 0.2) +
  scale_fill_manual(values = c("Genetics Driven" = "blue", "Environment Driven" = "green")) +
  xlab("Number of Significant Disease Hits") +
  ylab("Expected Mixture of Indirect Effect log(E/G)") +
  theme_minimal()




#Cartesian Plot: Specific Proteins
Type1 %>%
  filter(ID == "LEP") %>%
  ggplot(aes(x = G_HRi,
             y = E_HRi)) +
  geom_point() +
  ggtitle("LEP Indirect Effects of G and E on Multiple Diseases")

Type5 %>%
  filter(ID == "LEP") %>%
  ggplot(aes(x = `Genetic Indirect Effect`,
             y = `Exposure Indirect Effect`)) +
  geom_point() +
  ggtitle("LEP Indirect Effects of G and E on Multiple Diseases")

Type5 %>%
  filter(ID == "IGFBP2") %>%
  ggplot(aes(x = `Genetic Indirect Effect`,
             y = `Exposure Indirect Effect`)) +
  geom_point() +
  ggtitle("IGFBP2 Indirect Effects of G and E on Multiple Diseases")

Type5 %>%
  filter(ID == "GDF15") %>%
  ggplot(aes(x = `Genetic Indirect Effect`,
             y = `Exposure Indirect Effect`)) +
  geom_point() +
  ggtitle("GDF15 Indirect Effects of G and E on Multiple Diseases")

Type5 %>%
  filter(ID == "CA14") %>%
  ggplot(aes(x = `Genetic Indirect Effect`,
             y = `Exposure Indirect Effect`)) +
  geom_point() +
  ggtitle("CA14 Indirect Effects of G and E on Multiple Diseases")

Type5 %>%
  filter(ID == "APOF") %>%
  ggplot(aes(x = `Genetic Indirect Effect`,
             y = `Exposure Indirect Effect`)) +
  geom_point() +
  ggtitle("APOF Indirect Effects of G and E on Multiple Diseases")

Type5 %>%
  filter(ID == "REN") %>%
  ggplot(aes(x = `Genetic Indirect Effect`,
             y = `Exposure Indirect Effect`)) +
  geom_point() +
  ggtitle("REN Indirect Effects of G and E on Multiple Diseases")

Type5 %>%
  filter(ID == "TLR3") %>%
  ggplot(aes(x = `Genetic Indirect Effect`,
             y = `Exposure Indirect Effect`)) +
  geom_point() +
  ggtitle("TLR3 Indirect Effects of G and E on Multiple Diseases")


Type5 %>%
  filter(ID == "DKKL1") %>%
  ggplot(aes(x = `Genetic Indirect Effect`,
             y = `Exposure Indirect Effect`)) +
  geom_point() +
  ggtitle("DKKL1 Indirect Effects of G and E on Multiple Diseases")

#POLAR Coordinates:
# Vectorized version of the absolute angle function
absolute_angle_from_x <- function(x, y) {
  angle_rad <- atan2(y, x)
  angle_deg <- abs(angle_rad * (180 / pi))
  
  # Ensure angle is relative to x-axis (0 to 90 degrees): for any size vector
  angle_deg <- ifelse(angle_deg > 90, 180 - angle_deg, angle_deg)
  
  return(angle_deg)
}

# Test examples
absolute_angle_from_x(1, 0)   # Expect 0 degrees
absolute_angle_from_x(-1, 0)  # Expect 0 degrees
absolute_angle_from_x(0, 1)   # Expect 90 degrees
absolute_angle_from_x(0, -1)  # Expect 90 degrees
absolute_angle_from_x(1, 1)   # Expect 45 degrees
absolute_angle_from_x(-1, -1) # Expect 45 degrees


#Come up with more global plots: Cartesian Plots - multiple disease
Type1_angles <- Type1 %>%
  mutate(angle = absolute_angle_from_x(`Genetic Indirect Effect`,`Exposure Indirect Effect`),  # Convert radians to degrees
         magnitude = sqrt(`Exposure Indirect Effect`^2 + `Genetic Indirect Effect`^2)
         ) 


##Plotting: Cartesian and Radial Plot:
disease_types <- c("age_e11_first_reported_non_insulin_dependent_diabetes_mellitus_f130708_0_0",
                   "age_e66_first_reported_obesity_f130792_0_0",
                   "age_f00_first_reported_dementia_in_alzheimers_disease_f130836_0_0",
                   "age_i70_first_reported_atherosclerosis_f131380_0_0",
                   "age_n18_first_reported_chronic_renal_failure_f132032_0_0"
)
Type1_angles %>%
  filter(DZid %in% disease_types) %>%
  ggplot(aes(x = angle,
             y = magnitude)) +
  geom_point() +
  facet_wrap(~ DZid) +
  theme_minimal() +
  theme(legend.position = "none") 

Type1_angles <- Type1_angles %>%
  mutate(angle_rad = angle * (pi / 180))

# Create the radial plot
Type1_angles %>%
  filter(DZid %in% disease_types) %>%
  ggplot(aes(x = angle_rad, y = magnitude)) +
  geom_point(size = 4, color = "blue") +  # Plot the points
  geom_line(color = "blue") +  # Connect the points with a line
  coord_radial(start = 0, end = 0.5 * pi) +  # Use radial coordinates
  # scale_x_continuous(limits = c(0, pi/2),  # Limit the x-axis to 0 to 90 degrees (quarter circle)
  #                    labels = scales::number_format(scale = 1, accuracy = 1),
  #                    breaks = c(0, pi/4, pi/2)) +  # Define breaks (0°, 45°, 90°)
  facet_wrap(~ DZid) +
  theme_minimal() +  # Minimal theme
  labs(title = "Radial Plot with Angle and Magnitude", 
       x = "Angle (degrees)", y = "Magnitude") + 
  theme(axis.text.x = element_text(angle = 0))


Type1_angles %>%
  filter(DZid %in% disease_types) %>%
  ggplot(aes(x = `Genetic Indirect Effect`,
             y = `Exposure Indirect Effect`, color = DZid)) +
  geom_point() +
  facet_wrap(~ DZid) +
  theme_minimal() +
  theme(legend.position = "none") 






Type1_angles <- Type1 %>%
  mutate(angle = atan2(`Exposure Indirect Effect`, `Genetic Indirect Effect`) * 180 / pi,  # Convert radians to degrees
         magnitude = sqrt(`Exposure Indirect Effect`^2 + `Genetic Indirect Effect`^2)) 

length(unique(Type1_angles$DZid))

disease_types <- c("age_e11_first_reported_non_insulin_dependent_diabetes_mellitus_f130708_0_0",
                   "age_e66_first_reported_obesity_f130792_0_0",
                   "age_f00_first_reported_dementia_in_alzheimers_disease_f130836_0_0",
                   "age_h40_first_reported_glaucoma_f131186_0_0",
                   "age_i25_first_reported_chronic_ischaemic_heart_disease_f131306_0_0",
                   "age_i70_first_reported_atherosclerosis_f131380_0_0",
                   "age_n18_first_reported_chronic_renal_failure_f132032_0_0",
                   "age_m06_first_reported_other_rheumatoid_arthritis_f131850_0_0",
                   "age_j44_first_reported_other_chronic_obstructive_pulmonary_disease_f131492_0_0",
                   "age_j43_first_reported_emphysema_f131490_0_0",
                   "age_i26_first_reported_pulmonary_embolism_f131308_0_0"
                   )

Type1_angles %>%
  filter(DZid %in% disease_types) %>%
  ggplot(aes(x = `Genetic Indirect Effect`,
         y = `Exposure Indirect Effect`, color = DZid)) +
  geom_point() +
  facet_wrap(~ DZid) +
  theme_minimal() +
  theme(legend.position = "none") 

Type5 %>%
  filter(DZid %in% disease_types) %>%
  filter(`Genetic Indirect Effect` < 3) %>%
  ggplot(aes(x = `Genetic Indirect Effect`,
             y = `Exposure Indirect Effect`, color = DZid)) +
  geom_point() +
  facet_wrap(~ DZid) +
  theme_minimal() +
  theme(legend.position = "none") 


Type1_angles %>%
  filter(DZid %in% unique(Type1_angles$DZid)[10:20]) %>%
  ggplot(aes(x = DZid, y = angle, fill = DZid)) +
  geom_violin(trim = FALSE) +                      
  geom_jitter(width = 0.2, size = 2, alpha = 0.7) + 
  labs(
    title = "Distribution of Angles by Topic",
    x = "Topic",
    y = "Angle (Degrees)"
  ) +
  theme_minimal() +
  theme(legend.position = "none") 


Future_Obesity <- AllCovarSpec %>%
              filter(grepl("obesity",DZid))
CC <- AllCovarSpec %>%
  filter(grepl("cholecystitis",DZid))
unique(AllCovarSpec$DZid)

asthma <- AllCovarSpec %>%
  filter(grepl("asthma",DZid))

osteoporosis <- AllCovarSpec %>%
  filter(grepl("osteoporosis",DZid))
unique(AllCovarSpec$DZid)

# Clust_Disease <- AllCovarSpec %>%
#                   filter(CovarSpec == "Type5") %>%
#                   filter(!((HR_lower95 < 1 & HR_upper95 > 1) | (HR_lower95 > 1 & HR_upper95 < 1))) %>%
#                   select(all_of(c("Genetic Indirect Effect","Exposure Indirect Effect","ID","DZid")))
#                   
# # Figure out Clustering:
# 
# Clust_Disease_E <- Clust_Disease %>%
#                       select(all_of(c("Exposure Indirect Effect","ID","DZid"))) %>%
#                       pivot_wider(names_from = "DZid", 
#                                    values_from = "Exposure Indirect Effect") %>%
#                       column_to_rownames(var = "ID")
# Clust_Disease_E[is.na(Clust_Disease_E)] <- 0
# 
# Clust_Disease_G <- Clust_Disease %>%
#   select(all_of(c("Genetic Indirect Effect","ID","DZid"))) %>%
#   pivot_wider(names_from = "DZid", 
#               values_from = "Genetic Indirect Effect") %>%
#   column_to_rownames(var = "ID")
# Clust_Disease_G[is.na(Clust_Disease_G)] <- 0
# 
# 
# 
# library(ComplexHeatmap)
# ht1 <- Heatmap(as.matrix(Clust_Disease_E), 
#               name = "Mediation E",
#               show_row_names = FALSE, 
#               show_column_names = FALSE)
# 
# ht2 <- Heatmap(as.matrix(Clust_Disease_G), 
#               name = "Mediation G",
#               show_row_names = FALSE, 
#               show_column_names = FALSE)
# 
# draw(ht1)
# draw(ht2)
# 
# #PCA:
# pca <- prcomp(Clust_Disease_E, scale = T)
# pca_scores <- as.data.frame(pca$x)
# 
# var_explained <- pca$sdev^2 / sum(pca$sdev^2)
# 
# # Scree plot
# plot(var_explained, xlab = "Principal Component", ylab = "Proportion of Variance Explained",
#      type = "b", main = "Scree Plot")
# 
# # Plot PC1 vs. PC2
# ggplot(pca_scores, aes(x = PC1, y = PC2)) +
#   geom_point(color = "blue") +
#   labs(title = "PCA: PC1 vs PC2", x = "PC1", y = "PC2") +
#   theme_minimal()
# 
# 
# 
# library(ComplexHeatmap)
# ht <- Heatmap(as.matrix(Clust_Disease_E), 
#         name = "Mediation E",
#         show_row_names = FALSE, 
#         show_column_names = FALSE)
# 
# ht = draw(ht)
# # Extract column dendrogram (clustering)
# col_clusters <- column_dend(ht)
# 
# # Print or plot the column dendrogram to see the clusters
# print(col_clusters)
# 
# # Basic plot of the dendrogram
# plot(col_clusters, main = "Column Clusters Dendrogram")
# 
# library(dendextend)
# # Convert to dendrogram (if necessary)
# dend <- as.dendrogram(col_clusters)
# 
# # Now you can use cutree to split into k clusters (e.g., k = 3)
# k <- 5
# cluster_assignments <- cutree(dend, k = k)
# 
# # Print cluster assignments for each column
# print(cluster_assignments)
# 
# # Convert to a list of column names for each cluster
# cluster_names <- split(names(cluster_assignments), cluster_assignments)
# 
# # Print the clusters with corresponding column names
# print(cluster_names)



diabetes_rel <- AllCovarSpec %>%
  filter(grepl("diabetes", DZid))


diabetes_select <- diabetes_rel %>%
  filter(abs(`Exposure Indirect Effect`) > abs(`Genetic Indirect Effect`)) %>%
  filter(!((HR_lower95 < 1 & HR_upper95 > 1) | (HR_lower95 > 1 & HR_upper95 < 1)))

T2D <- diabetes_rel %>%
  filter(grepl("e11", DZid)) %>%
  filter(abs(`Exposure Indirect Effect`) > abs(`Genetic Indirect Effect`)) %>%
  filter(!((HR_lower95 < 1 & HR_upper95 > 1) | (HR_lower95 > 1 & HR_upper95 < 1)))

T2D_all <- AllCovarSpec %>%
  filter(grepl("e11", DZid))


Type1 <- loadMDAssoc(covarType = "Type1")
Type1 <- do.call("rbind",Type1)

Type1_selection <- Type1%>%
                  filter(abs(`Exposure Indirect Effect`) > abs(`Exposure Direct Effect`) &
                           abs(`Exposure Indirect Effect`) > abs(`Genetic Indirect Effect`)) %>%
                    mutate(E_HRi = exp(`Exposure Indirect Effect`),
                           G_HRi = exp(`Genetic Indirect Effect`),
                           delta_HRi = exp(abs(`Exposure Indirect Effect` - `Genetic Indirect Effect`)),
                           total_HRi = exp(`Exposure Indirect Effect` + `Genetic Indirect Effect`),
                           Eprop_mediated = abs(`Exposure Indirect Effect`)/(abs(`Exposure Indirect Effect`) + abs(`Exposure Direct Effect`)),
                           Gprop_mediated = abs(`Genetic Indirect Effect`)/(abs(`Genetic Indirect Effect`) + abs(`Genetic Direct Effect`)))

Type2 <- loadMDAssoc(covarType = "Type2")
Type2 <- do.call("rbind",Type2)

Type2_selection <- Type2 %>%
  filter(abs(`Exposure Indirect Effect`) > abs(`Exposure Direct Effect`) &
           abs(`Exposure Indirect Effect`) > abs(`Genetic Indirect Effect`)) %>%
  mutate(E_HRi = exp(`Exposure Indirect Effect`),
         G_HRi = exp(`Genetic Indirect Effect`),
         delta_HRi = exp(abs(`Exposure Indirect Effect` - `Genetic Indirect Effect`)),
         total_HRi = exp(`Exposure Indirect Effect` + `Genetic Indirect Effect`),
         Eprop_mediated = abs(`Exposure Indirect Effect`)/(abs(`Exposure Indirect Effect`) + abs(`Exposure Direct Effect`)),
         Gprop_mediated = abs(`Genetic Indirect Effect`)/(abs(`Genetic Indirect Effect`) + abs(`Genetic Direct Effect`)))


Type3 <- loadMDAssoc(covarType = "Type3")
Type3 <- do.call("rbind",Type3)

Type3_selection <- Type3%>%
  filter(abs(`Exposure Indirect Effect`) > abs(`Exposure Direct Effect`) &
           abs(`Exposure Indirect Effect`) > abs(`Genetic Indirect Effect`))  %>%
  mutate(E_HRi = exp(`Exposure Indirect Effect`),
         G_HRi = exp(`Genetic Indirect Effect`),
         delta_HRi = exp(abs(`Exposure Indirect Effect` - `Genetic Indirect Effect`)),
         total_HRi = exp(`Exposure Indirect Effect` + `Genetic Indirect Effect`),
         Eprop_mediated = abs(`Exposure Indirect Effect`)/(abs(`Exposure Indirect Effect`) + abs(`Exposure Direct Effect`)),
         Gprop_mediated = abs(`Genetic Indirect Effect`)/(abs(`Genetic Indirect Effect`) + abs(`Genetic Direct Effect`)))



Type5 <- loadMDAssoc(covarType = "Type5")
Type5 <- do.call("rbind",Type5)

Type5_selection <- Type5 %>%
  filter(abs(`Exposure Indirect Effect`) > abs(`Exposure Direct Effect`) &
           abs(`Exposure Indirect Effect`) > abs(`Genetic Indirect Effect`)) %>%
  filter(!((HR_lower95 < 1 & HR_upper95 > 1) | (HR_lower95 > 1 & HR_upper95 < 1)))

#Make sure hazard does not cross 1.
colnames(Type3_selection)

Type5_genetic_selection <- Type5 %>%
  filter(abs(`Genetic Indirect Effect`) > abs(`Genetic Direct Effect`) &
           abs(`Exposure Indirect Effect`) < abs(`Genetic Indirect Effect`)) 


#Calculate Proportion Mediated:

Type5_selection <- Type5_selection %>%
  mutate(E_HRi = exp(`Exposure Indirect Effect`),
        G_HRi = exp(`Genetic Indirect Effect`),
        delta_HRi = exp(abs(`Exposure Indirect Effect` - `Genetic Indirect Effect`)),
        total_HRi = exp(`Exposure Indirect Effect` + `Genetic Indirect Effect`),
        Eprop_mediated = abs(`Exposure Indirect Effect`)/(abs(`Exposure Indirect Effect`) + abs(`Exposure Direct Effect`)),
        Gprop_mediated = abs(`Genetic Indirect Effect`)/(abs(`Genetic Indirect Effect`) + abs(`Genetic Direct Effect`)))

Type5_upd <- Type5 %>%
  mutate(E_HRi = exp(`Exposure Indirect Effect`),
         G_HRi = exp(`Genetic Indirect Effect`),
         delta_HRi = abs(`Exposure Indirect Effect` - `Genetic Indirect Effect`), #exp(abs(`Exposure Indirect Effect` - `Genetic Indirect Effect`))
         total_HRi = exp(`Exposure Indirect Effect` + `Genetic Indirect Effect`),
         Eprop_mediated = abs(`Exposure Indirect Effect`)/(abs(`Exposure Indirect Effect`) + abs(`Exposure Direct Effect`)),
         Gprop_mediated = abs(`Genetic Indirect Effect`)/(abs(`Genetic Indirect Effect`) + abs(`Genetic Direct Effect`)))

Type5_genetic_selectionv2 <- Type5 %>%
  filter(abs(`Exposure Indirect Effect`) < abs(`Genetic Indirect Effect`))  %>%
  filter(!(cindex == 1)) %>%
  filter(R2 > 0.01) %>%
  mutate(E_HRi = exp(`Exposure Indirect Effect`),
         G_HRi = exp(`Genetic Indirect Effect`),
         delta_HRi = exp(abs(`Exposure Indirect Effect` - `Genetic Indirect Effect`)),
         total_HRi = exp(`Exposure Indirect Effect` + `Genetic Indirect Effect`),
         Eprop_mediated = abs(`Exposure Indirect Effect`)/(abs(`Exposure Indirect Effect`) + abs(`Exposure Direct Effect`)),
         Gprop_mediated = abs(`Genetic Indirect Effect`)/(abs(`Genetic Indirect Effect`) + abs(`Genetic Direct Effect`)))




#Selection Criteria for E influence:
Type5_updv1 <- Type5_upd %>%
                filter(abs(`Exposure Indirect Effect`) > abs(`Exposure Direct Effect`) &
                         abs(`Exposure Indirect Effect`) > abs(`Genetic Indirect Effect`)) %>%
                filter(!((HR_lower95 < 1 & HR_upper95 > 1) | (HR_lower95 > 1 & HR_upper95 < 1)))

Type5_updv2 <- Type5_upd %>%
  filter(abs(`Exposure Indirect Effect`) > abs(`Genetic Indirect Effect`)) %>%
  filter(!((HR_lower95 < 1 & HR_upper95 > 1) | (HR_lower95 > 1 & HR_upper95 < 1)))



#### PLOTS that make sense:

Type5_selection$HR_upper95 <- as.numeric(Type5_selection$HR_upper95)

unique(Type5_selection$ID)


Etophits_plot <- Type5_selection %>% 
  select(all_of(c("ID","DZid","E_HRi","G_HRi","HR"))) %>%
  pivot_longer(
    cols = c("E_HRi","G_HRi"),
    names_to = "Type",
    values_to = "Effects"
  )

#Without confidence intervals:
Etophits_plot %>%
  filter(grepl("e11",DZid)) %>%
  ggplot(aes(x = ID, y = Effects, fill = Type)) +
  geom_bar(position="dodge", stat="identity")




##'*Make below a function:*

Etophits_plot %>%
  filter(grepl("diabetes",DZid)) %>%
  ggplot(aes(x = ID, y = Effects, color = Type)) +
  geom_point(position = position_dodge(width = 0.6), size = 3) +  # Point estimate
  #geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2,  # Error bars
  #               position = position_dodge(width = 0.6)) +
  #scale_color_manual(values = c("Group1" = "blue", "Group2" = "red")) +
  scale_color_discrete(labels = c("E", "G")) +
  labs(x = "Proteins", y = "Hazard Ratio: Mediation Effects") +
  coord_flip() +
  theme_minimal() #+
  #theme(axis.title.x = element_blank())

Etophits_plot %>%
  filter(grepl("diabetes",DZid)) %>%
  ggplot(aes(x = ID, y = HR)) +
  geom_point(position = position_dodge(width = 0.6), size = 3) +  # Point estimate
  #geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2,  # Error bars
  #               position = position_dodge(width = 0.6)) +
  #scale_color_manual(values = c("Group1" = "blue", "Group2" = "red")) +
  labs(x = "Proteins", y = "Hazard Ratio: Protein Effects") +
  coord_flip() +
  theme_minimal()

Etophits_plot %>%
  filter(grepl("n18",DZid)) %>%
  ggplot(aes(x = ID, y = Effects, color = Type)) +
  geom_point(position = position_dodge(width = 0.6), size = 3) +  # Point estimate
  #geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2,  # Error bars
  #               position = position_dodge(width = 0.6)) +
  #scale_color_manual(values = c("Group1" = "blue", "Group2" = "red")) +
  scale_color_discrete(labels = c("E", "G")) +
  labs(x = "Proteins", y = "Hazard Ratio: Mediation Effects") +
  coord_flip() +
  theme_minimal()





Etophits_plot1 <- Etophits_plot %>%
  filter(Type %in% c("Genetic Indirect Effect","Exposure Indirect Effect")) 

Etophits_plot1$Type <- factor(Etophits_plot1$Type,
                              levels = c("Genetic Indirect Effect","Exposure Indirect Effect"))








#Slope plot
Etophits_plot1 %>%
  filter(Type %in% c("Genetic Indirect Effect","Exposure Indirect Effect")) %>%
  ggplot(aes(x = Type, y = Effects, group = ID)) +
  geom_line(aes(color = ID), size = 1) +  # Add lines connecting the effects
  geom_point(aes(color = ID), size = 3) + # Add points for each effect
  theme_minimal() +
  xlab("Mediation Effects") +
  ylab("Mediation Effect Size") +
  labs(color = "Protein") +
  theme(axis.title.x = element_blank())

















REN <- Type5 %>%
              filter(ID == "REN")

#GFAP:
GFAP <- Type5_upd %>%
            filter(ID == "GFAP")
LEP <- Type5_upd %>%
  filter(ID == "LEP")
CA14 <- Type5_upd %>%
  filter(ID == "CA14")
IGFBP2 <- Type5_upd %>%
  filter(ID == "IGFBP2")
CKB <- Type5_upd %>%
  filter(ID == "CKB")
IL1RN <- Type5_upd %>%
  filter(ID == "IL1RN")
FABP4 <- Type5_upd %>%
  filter(ID == "FABP4")
APOF <- Type5_upd %>%
  filter(ID == "APOF")

#Disease Specific:
e11 <- Type5_upd %>%
  filter(grepl("e11",DZid))
plot(e11$G_HRi, e11$E_HRi)

e78 <- Type5_upd %>%
  filter(grepl("e78",DZid))

n18 <- Type5_upd %>%
  filter(grepl("n18",DZid))

n18_E_effect <- Type5_selection %>%
                  filter(grepl("n18",DZid))
#FABP4

APOE <- Type5 %>%
          filter(ID == "APOE")


GFAP <- Type5 %>%
  filter(ID == "GFAP")

NEFL <- Type5 %>%
  filter(ID == "NEFL")


CKB <- Type5 %>%
  filter(ID == "CKB")





#'*Obtain Results from nested CV*
T2D_MD <- data.frame()
for(i in 1:50){ #200 is based on the split designated in the bash script.
  tryCatch({
    MD_idx <- fread(file=paste0("/n/groups/patel/shakson_ukb/UK_Biobank/Data/Parallel/Prot_PGSPXS_final/CovSpec2_v2/","T2D_Med_",i,".txt"))
    
    #Updates:
    T2D_MD  <- rbind(T2D_MD, MD_idx)
  }, error = function(e) {
    # Handle the error (e.g., print an error message)
    message(paste("Error reading file:",i))
  })
}


E_driven_T2D <- T2D_MD %>% 
  filter(abs(`Genetic Indirect Effect`) < abs(`Exposure Indirect Effect`))

Etophits <- E_driven_T2D %>% 
  filter(abs(`Total Direct Effect`) < abs(`Total Indirect Effect`)) %>%
  filter(abs(`Total Effect`) > 0.1)


G_driven_T2D <- T2D_MD %>% 
  filter(abs(`Genetic Indirect Effect`) > abs(`Exposure Indirect Effect`))

Gtophits <- G_driven_T2D %>% 
  filter(abs(`Total Direct Effect`) < abs(`Total Indirect Effect`)) %>%
  filter(abs(`Total Effect`) > 0.1)


Gtophits_v2 <- G_driven_T2D %>% 
  filter(abs(`Total Direct Effect`) < abs(`Total Indirect Effect`)) %>%
  filter(abs(`Total Effect`) > 0.1) %>%
  filter(abs(`Genetic Direct Effect`) < abs(`Genetic Indirect Effect`))



#'*Visualization of Top Hits:*
Etophits_plot <- Etophits %>% 
  pivot_longer(
    cols = -ID,
    names_to = "Type",
    values_to = "Effects"
  )
  
#Without confidence intervals:
Etophits_plot %>%
  ggplot(aes(x = ID, y = Effects, fill = Type)) +
  geom_bar(position="dodge", stat="identity")

Etophits_plot1 <- Etophits_plot %>%
  filter(Type %in% c("Genetic Indirect Effect","Exposure Indirect Effect")) 

Etophits_plot1$Type <- factor(Etophits_plot1$Type,
                              levels = c("Genetic Indirect Effect","Exposure Indirect Effect"))

#Slope plot
Etophits_plot1 %>%
  filter(Type %in% c("Genetic Indirect Effect","Exposure Indirect Effect")) %>%
  ggplot(aes(x = Type, y = Effects, group = ID)) +
  geom_line(aes(color = ID), size = 1) +  # Add lines connecting the effects
  geom_point(aes(color = ID), size = 3) + # Add points for each effect
  theme_minimal() +
  xlab("Mediation Effects") +
  ylab("Mediation Effect Size") +
  labs(color = "Protein") +
  theme(axis.title.x = element_blank())

#Grouped barchart
Etophits_plot1 %>%
  filter(Type %in% c("Genetic Indirect Effect","Exposure Indirect Effect")) %>%
  ggplot(aes(x = ID, y = Effects, fill = Type)) +
  geom_bar(position="dodge", stat="identity") +
  theme_minimal() +
  xlab("Mediation Effects") +
  ylab("Mediation Effect Size") +
  labs(color = "Protein") +
  theme(axis.title.x = element_blank())

#Grouped Forest plot
Etophits_plot1 %>%
  ggplot(aes(x = ID, y = Effects, color = Type)) +
  geom_point(position = position_dodge(width = 0.6), size = 3) +  # Point estimate
  #geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2,  # Error bars
  #               position = position_dodge(width = 0.6)) +
  #scale_color_manual(values = c("Group1" = "blue", "Group2" = "red")) +
  labs(x = "Proteins", y = "Mediation Effects") +
  theme_minimal() +
  theme(axis.title.x = element_blank())




Gtophits_plot <- Gtophits %>% 
  pivot_longer(
    cols = -ID,
    names_to = "Type",
    values_to = "Effects"
  )


Gtophits_plot1 <- Gtophits_plot %>%
  filter(Type %in% c("Genetic Indirect Effect","Exposure Indirect Effect")) 

Gtophits_plot1$Type <- factor(Gtophits_plot1$Type,
                              levels = c("Genetic Indirect Effect","Exposure Indirect Effect"))

Gtophits_plot1 %>%
  ggplot(aes(x = factor(Type, levels = c("Genetic Indirect Effect","Exposure Indirect Effect")), 
             y = Effects, group = ID)) +
  geom_line(aes(color = ID), size = 1) +  # Add lines connecting the effects
  geom_point(aes(color = ID), size = 3) + # Add points for each effect
  theme_minimal() +
  xlab("Mediation Effects") +
  ylab("Mediation Effect Size") +
  labs(color = "Protein")



##TO-DO: Maybe just do proteome wide analysis also...
##Visualization:


T2D_MD %>%
  ggplot(aes(x = `Genetic Indirect Effect`, y = `Exposure Indirect Effect`)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1)



T2D_MD %>%
  filter(abs(`Total Direct Effect`) < abs(`Total Indirect Effect`)) %>%
  ggplot(aes(x = `Genetic Indirect Effect`, y = `Exposure Indirect Effect`)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1)


T2D_MD %>%
  #filter(abs(`Total Direct Effect`) < abs(`Total Indirect Effect`)) %>%
  ggplot(aes(x = `Total Effect`, y = `Exposure Indirect Effect`)) +
  geom_point() 
  #+ geom_abline(intercept = 0, slope = 1)



T2D_MD %>%
  ggplot(aes(x = `Total Direct Effect`, y = `Total Indirect Effect`)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1)
  


fit <- lm(`Exposure Indirect Effect` ~ `Genetic Indirect Effect`, data = T2D_MD)
summary(fit)
#'*Maybe this gives intuition for other diseases?*
#'If T2D is 0.18 and ALZ is 0.01 then this suggest more of an E influence then G?

library(ggplot2)
library(ggpmisc)
library(ggpubr)
T2D_MD %>%
  filter(abs(`Total Direct Effect`) < abs(`Total Indirect Effect`)) %>%
  ggplot(aes(x = `Genetic Indirect Effect`, y = `Exposure Indirect Effect`)) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  stat_regline_equation(label.x = -0.5, label.y = 0.5, aes(label = ..eq.label..)) +
  stat_cor(label.x = -0.5, label.y = 0.4, aes(label = ..r.label..)) +
  #stat_poly_line() +
  #stat_poly_eq() +
  geom_point()

#geom_smooth(method = "lm", se = FALSE, color = "blue") +
#stat_cor(method = "spearman", label.x = -2, label.y = 3)
#stat_poly_line() +
#stat_poly_eq() +


T2D_MD %>%
  ggplot(aes(x = `Genetic Indirect Effect`, y = `Exposure Indirect Effect`)) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  stat_regline_equation(label.x = -0.5, label.y = 0.5, aes(label = ..eq.label..)) +
  stat_cor(label.x = -0.5, label.y = 0.4, aes(label = ..r.label..)) +
  geom_point()


T2D_MD %>%
  ggplot(aes(x = `Genetic Direct Effect`, y = `Exposure Direct Effect`)) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  stat_regline_equation(label.x = -0.5, label.y = 0.5, aes(label = ..eq.label..)) +
  stat_cor(label.x = -0.5, label.y = 0.4, aes(label = ..r.label..)) +
  geom_point()


#Notable hits: IGFBP2 & IGFBP1
#IDEA - FIND the best proteins here --> then run bootstrap to show confidence intervals.











#'*HEATMAPS find CLUSTERS of E --> Prot --> Disease vs G --> Prot --> Disease*
Clust_Disease <- AllCovarSpec %>%
                  filter(CovarSpec == "Type5") %>%
                  filter(!((HR_lower95 < 1 & HR_upper95 > 1) | (HR_lower95 > 1 & HR_upper95 < 1))) %>%
                  select(all_of(c("Genetic Indirect Effect","Exposure Indirect Effect","ID","DZid")))

# Figure out Clustering:
Clust_Disease_E <- Clust_Disease %>%
                      select(all_of(c("Exposure Indirect Effect","ID","DZid"))) %>%
                      pivot_wider(names_from = "DZid",
                                   values_from = "Exposure Indirect Effect") %>%
                      column_to_rownames(var = "ID")
Clust_Disease_E[is.na(Clust_Disease_E)] <- 0

Clust_Disease_G <- Clust_Disease %>%
  select(all_of(c("Genetic Indirect Effect","ID","DZid"))) %>%
  pivot_wider(names_from = "DZid",
              values_from = "Genetic Indirect Effect") %>%
  column_to_rownames(var = "ID")
Clust_Disease_G[is.na(Clust_Disease_G)] <- 0





library(ComplexHeatmap)
#Put in HAZARD scale:
Clust_Disease_Ehr <- exp(Clust_Disease_E)
Clust_Disease_Ghr <- exp(Clust_Disease_G)

ht1 <- Heatmap(as.matrix(Clust_Disease_Ehr),
              name = "Mediation E",
              show_row_names = FALSE,
              show_column_names = FALSE)
ht1



ht2 <- Heatmap(as.matrix(Clust_Disease_Ghr),
               name = "Mediation E",
               show_row_names = FALSE,
               show_column_names = FALSE)
ht2


#'*threshold by variability:*
#Put in HAZARD scale:
Clust_Disease_Ehr <- exp(Clust_Disease_E)
matrix <- Clust_Disease_Ehr

row_sd <- apply(matrix, 1, sd)  # Standard deviation for rows
col_sd <- apply(matrix, 2, sd)  # Standard deviation for columns

upper_bound_row <- quantile(row_sd, 0.99)
upper_bound_col <- quantile(col_sd, 0.99)

#Get thresholds
threshold_row <- median(row_sd) #median(row_sd)
threshold_col <- median(col_sd) #median(col_sd)

sum(row_sd < upper_bound_row)
sum(row_sd > threshold_row)

filtered_matrix <- matrix[row_sd > threshold_row & row_sd < upper_bound_row, ]  # Keep rows with SD > ?
filtered_matrix <- filtered_matrix[, col_sd > threshold_col & col_sd < upper_bound_col]  # Keep columns with SD > ?

# keep_rows <- row_sd > threshold
# keep_cols <- col_sd > threshold
# filtered_matrix <- matrix[keep_rows, keep_cols]

ht1 <- Heatmap(filtered_matrix,
               name = "Mediation E",
               show_row_names = FALSE,
               show_column_names = FALSE)
ht1

#### G version:
#Put in HAZARD scale:
Clust_Disease_Ghr <- exp(Clust_Disease_G)
matrix <- Clust_Disease_Ghr

row_sd <- apply(matrix, 1, sd)  # Standard deviation for rows
col_sd <- apply(matrix, 2, sd)  # Standard deviation for columns

upper_bound_row <- quantile(row_sd, 0.99)
upper_bound_col <- quantile(col_sd, 0.99)

#Get thresholds
threshold_row <- median(row_sd) #median(row_sd)
threshold_col <- median(col_sd) #median(col_sd)

sum(row_sd < upper_bound_row)
sum(row_sd > threshold_row)

filtered_matrix <- matrix[row_sd > threshold_row & row_sd < upper_bound_row, ]  # Keep rows with SD > ?
filtered_matrix <- filtered_matrix[, col_sd > threshold_col & col_sd < upper_bound_col]  # Keep columns with SD > ?

# keep_rows <- row_sd > threshold
# keep_cols <- col_sd > threshold
# filtered_matrix <- matrix[keep_rows, keep_cols]

ht1 <- Heatmap(filtered_matrix,
               name = "Mediation G",
               show_row_names = FALSE,
               show_column_names = FALSE)
ht1



Clust_Disease_Ehr <- exp(Clust_Disease_E)
Clust_Disease_Ghr <- exp(Clust_Disease_G)
Clust_Disease_Prothr <- Clust_Disease_Ehr * Clust_Disease_Ghr
matrix <- Clust_Disease_Prothr

row_sd <- apply(matrix, 1, sd)  # Standard deviation for rows
col_sd <- apply(matrix, 2, sd)  # Standard deviation for columns

upper_bound_row <- quantile(row_sd, 0.99)
upper_bound_col <- quantile(col_sd, 0.99)

#Get thresholds
threshold_row <- median(row_sd) #median(row_sd)
threshold_col <- median(col_sd) #median(col_sd)

sum(row_sd < upper_bound_row)
sum(row_sd > threshold_row)

filtered_matrix <- matrix[row_sd > threshold_row & row_sd < upper_bound_row, ]  # Keep rows with SD > ?
filtered_matrix <- filtered_matrix[, col_sd > threshold_col & col_sd < upper_bound_col]  # Keep columns with SD > ?

# keep_rows <- row_sd > threshold
# keep_cols <- col_sd > threshold
# filtered_matrix <- matrix[keep_rows, keep_cols]

ht1 <- Heatmap(filtered_matrix,
               name = "Prot HR",
               show_row_names = FALSE,
               show_column_names = FALSE)
ht1




#'*WHAT IS SHARED Structure of Assocations here from Mediation Analysis*

#PCA - E Indirect Effect

pca_stats <- function(df){
  PCstats <- list()
  
  pca <- prcomp(df, scale = T)
  
  pca_loadings <- as.data.frame(pca$rotation)
  pca_scores <- as.data.frame(pca$x)
  
  var_explained <- pca$sdev^2 / sum(pca$sdev^2)
  
  # Scree plot
  plot(var_explained, xlab = "Principal Component", ylab = "Proportion of Variance Explained",
       type = "b", main = "Scree Plot")
  
  # Plot PC1 vs. PC2
  p1 <- ggplot(pca_scores, aes(x = PC1, y = PC2)) +
    geom_point(color = "blue") +
    labs(title = "PCA: PC1 vs PC2", x = "PC1", y = "PC2") +
    theme_minimal()
  print(p1)
  
  PCstats[["loadings"]] <- pca_loadings
  PCstats[["scores"]] <- pca_scores
  
  return(PCstats)
}

Emediation <- pca_stats(t(Clust_Disease_E))
Gmediation <- pca_stats(t(Clust_Disease_G))

Emediation_prot <- pca_stats(Clust_Disease_E)
Gmediation_prot <- pca_stats(Clust_Disease_G)

#### DIFFERENCE MATRIX!!!!
Clust_GmE <- Clust_Disease_G - Clust_Disease_E
GmE_mediation <- pca_stats(t(Clust_GmE))
GmE_mediation_prot <- pca_stats(Clust_GmE)


Clust_EmG <- Clust_Disease_E -Clust_Disease_G
EmG_mediation <- pca_stats(t(Clust_EmG))
EmG_mediation_prot <- pca_stats(Clust_EmG)



