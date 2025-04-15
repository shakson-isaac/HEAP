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

#Figures:
# General overall averaged result of mediation
# Variation of GEM statistic
# Highlight Specific Diseases (HR and c-index)
# Highlight E v G partitioning (Plot) - Bootstrapping Procedure (Script #2: mediationboot.R)
# Highlight variability across covariate specifications

#Tables:
#Save significant summary stats as excel files for tables.

#Tips:
#'HAVE geom_point come after geom_ribbon to make sure the points POP OUT!!!
#'Next usage: Use p-value instead of HR for each protein - counting associations*


#'*Load Files of Mediation Results*
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


#'*Obtain OLS slope across Specifications:*
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
# Type5_o <- Type5_o %>%
#   mutate(GEM = log(Estimate))


#'*Plot Variation of GEM across Covariate Specifications*
GEMres <- MD_GvE_effects %>% 
            group_by(ID) %>%
            summarise(meanGEM = mean(log(Estimate)),
                      stdErrorGEM = sd(log(Estimate))/sqrt(n()))

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

#'*Plot correlation of GEM vs Number of Disease Associations*
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


#'*Plot of GEM vs number of disease hits: Maximal Specification (Type 5)*
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


#Function: Volcano Plot of Protein Associations to a Given Disease
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


#'*Plot Protein Volcano Plot For 3 Diseases: T2D, Emphysema, Mental Disorders/Alcohol Abuse*
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

# Additional Volcano Plots to Consider:
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
