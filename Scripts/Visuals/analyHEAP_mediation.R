#Libraries:
library(data.table)
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(ggpmisc)
library(pbapply)
library(dplyr)
library(purrr)
library(qs) #For faster save/load of rds object.
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

#AllCovarSpec <- list(Type1, Type2, Type3, Type4, Type5) %>% reduce(full_join)
MDSpec <- list(Type1, Type2, Type3, Type4, Type5)
names(MDSpec) <- c("Type1","Type2","Type3","Type4","Type5")

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

#'*Create DataStructure*
MDconstruct <- setClass(
  "MDconstruct",
  slots = c(
    MDlist = "list", #List: All association results
    GEMlist = "list" #List: Replicated results (Significant in both Train and Test)
  )
)

#Create the datastructure to utilize:
HEAPmd <- MDconstruct(
  MDlist = MDSpec,
  GEMlist = MDres
)

#Save RDS file of PXSloader object
gc()
class(HEAPmd)
qsave(HEAPmd, file = "./Output/HEAPres/HEAPmd.qs")