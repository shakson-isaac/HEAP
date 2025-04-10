#Libraries:
library(data.table)
library(tidyverse)
library(ggplot2)
library(polycor)

# Definition of PXSconstruct to make modifications to covars_list
PXSconstruct <- setClass(
  "PXSconstruct",
  slots = c(
    Elist = "list", #E dataframes separated by category:
    Elist_names = "character", #Names of the category:
    Eid_cat = "data.frame", #Dataframe with environmental IDs and corresponding category.
    ordinalIDs = "character", #list environmental ordinal names
    
    UKBprot_df = "data.frame", #Protein variables dataframe
    protIDs = "character", #protein variables
    
    covars_df = "data.frame", #Covariates Dataframe
    covars_list = "character" #covars variables
  )
)


#'*Loading:*
#1st ver:
#load(file = "/n/scratch/users/s/shi872/UKB_intermediate/UKB_PGS_PXS_data.rds")
PXSloader <- readRDS(file = "/n/scratch/users/s/shi872/UKB_intermediate/UKB_PGS_PXS_load.rds")

#By intra-category: Spearman Correlations (Takes 8 minutes)
ordinalIDs <- PXSloader@ordinalIDs
Ecor <- lapply(PXSloader@Elist, function(Edf){
  Ecor <- cor(Edf %>% select(!eid), 
      use = "pairwise.complete.obs",
      method = "spearman")
  return(Ecor)
})
names(Ecor) <- names(PXSloader@Elist)



# All E Categories: ~35 minutes
Eall <- PXSloader@Elist %>% reduce(full_join)
Ecorfin <- cor(Eall %>% select(!eid), 
            use = "pairwise.complete.obs",
            method = "spearman")

#Plot the correlations:
library(ggcorrplot)
Ecorfin[is.na(Ecorfin)] <- 0
gg1 <- ggcorrplot(Ecorfin, hc.order = TRUE, type = "lower",
           outline.col = "white") +  # for numbers inside tiles
            scale_fill_gradient2(
              low = "blue", mid = "white", high = "red",
              midpoint = 0, limit = c(-1, 1),
              name = expression(rho)  # for Greek ρ
              # Or use: name = "Spearman correlation"
            ) +
            theme(
              axis.text.x = element_text(angle = 90, size = 3),  # adjust as needed
              axis.text.y = element_text(size = 3)
            )
            

lapply(c("png", "svg", "pdf"), 
       function(fmt) ggsave(paste0("./Figures/HEAP/SF6/UKB_Ecor.", fmt), 
                            plot = gg1, width = 12, height = 8, dpi = 1000))



#Plot the correlations per category (Extra)
# Ecor$Alcohol[is.na(Ecor$Alcohol)] <- 0
# ggcorrplot(Ecor$Alcohol, hc.order = TRUE, type = "lower",
#            outline.col = "white") +  # for numbers inside tiles
#   theme(
#     axis.text.x = element_text(angle = 90, size = 3),  # adjust as needed
#     axis.text.y = element_text(size = 3)
#   )
# 
# Ecor$Smoking[is.na(Ecor$Smoking)] <- 0
# ggcorrplot(Ecor$Smoking, hc.order = TRUE, type = "lower",
#            outline.col = "white") +  # for numbers inside tiles
#   theme(
#     axis.text.x = element_text(angle = 90, size = 3),  # adjust as needed
#     axis.text.y = element_text(size = 3)
#   )
# 
# Ecor$Exercise_Freq[is.na(Ecor$Exercise_Freq)] <- 0
# ggcorrplot(Ecor$Exercise_Freq, hc.order = TRUE, type = "lower",
#            outline.col = "white") +  # for numbers inside tiles
#   theme(
#     axis.text.x = element_text(angle = 90, size = 3),  # adjust as needed
#     axis.text.y = element_text(size = 3)
#   )