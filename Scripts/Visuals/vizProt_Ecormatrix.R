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



#Plot the correlations:
library(ggcorrplot)
Ecor$Alcohol[is.na(Ecor$Alcohol)] <- 0
ggcorrplot(Ecor$Alcohol, hc.order = TRUE, type = "lower",
           outline.col = "white") +  # for numbers inside tiles
  theme(
    axis.text.x = element_text(angle = 90, size = 3),  # adjust as needed
    axis.text.y = element_text(size = 3)
  )

Ecor$Smoking[is.na(Ecor$Smoking)] <- 0
ggcorrplot(Ecor$Smoking, hc.order = TRUE, type = "lower",
           outline.col = "white") +  # for numbers inside tiles
  theme(
    axis.text.x = element_text(angle = 90, size = 3),  # adjust as needed
    axis.text.y = element_text(size = 3)
  )

Ecor$Exercise_Freq[is.na(Ecor$Exercise_Freq)] <- 0
ggcorrplot(Ecor$Exercise_Freq, hc.order = TRUE, type = "lower",
           outline.col = "white") +  # for numbers inside tiles
  theme(
    axis.text.x = element_text(angle = 90, size = 3),  # adjust as needed
    axis.text.y = element_text(size = 3)
  )


##### TRASH #####

#Check how correlated exposures are:
diet_example <- PXSloader@Elist$Diet_Weekly
sapply(diet_example, class)

ordinalIDs <- PXSloader@ordinalIDs
diet_examplev2 <- diet_example %>%
                      mutate(across(any_of(ordinalIDs),
                             ~ factor(.x, ordered = T)))
sapply(diet_examplev2, class)
 
#CHECKING:
#levels(diet_examplev2$beef_intake_f1369_0_0)

#Where everything is a number:
DietCorMatrix <- cor(diet_example %>% select(!eid), 
                     use = "pairwise.complete.obs",
                     method = "pearson")

#Takes around 6 minutes for spearman
DietCorMatrixSpear <- cor(diet_example %>% select(!eid), 
                     use = "pairwise.complete.obs",
                     method = "spearman")

DietCorMatrix[is.na(DietCorMatrix)] <- 0
DietCorMatrixSpear[is.na(DietCorMatrixSpear)] <- 0

library(ggcorrplot)
ggcorrplot(DietCorMatrix, hc.order = TRUE, type = "lower",
           outline.col = "white") +  # for numbers inside tiles
  theme(
    axis.text.x = element_text(angle = 90, size = 3),  # adjust as needed
    axis.text.y = element_text(size = 3)
  )


ggcorrplot(DietCorMatrixSpear, hc.order = TRUE, type = "lower",
           outline.col = "white") +  # for numbers inside tiles
  theme(
    axis.text.x = element_text(angle = 90, size = 3),  # adjust as needed
    axis.text.y = element_text(size = 3)
  )










##### HetCor doesnt work well: ####
diet_examplev3 <- diet_examplev2 %>% 
                      select(!eid) %>%
                      mutate(across(where(is.integer), as.numeric))
sapply(diet_examplev3, class)

diet_examplev3$beef_intake_f1369_0_0
diet_examplev3$milk_type_used_f1418_0_0_Do_not_know


diet_examplev4 <- na.omit(diet_examplev3)


try <- diet_examplev4[,2:3] 
sapply(try, class)

try <- as.data.frame(lapply(diet_examplev4[, 2:3], as.vector))

#'*Make sure all columns are base R atomic vectors*
#'
diet_example <- PXSloader@Elist$Diet_Weekly
diet_clean <- as.data.frame(lapply(diet_example, as.vector))
diet_clean <- diet_clean %>%
  #Remove eid id:
  select(!eid) %>%
  #Convert ordinal vectors appropriately
  mutate(across(any_of(ordinalIDs),
                ~ factor(.x, ordered = T))) %>%
  # Remove constant columns (zero variance)
  select(where(~ n_distinct(.) > 1)) %>%
  # Drop sparse levels from ordered factors
  mutate(across(where(is.ordered), ~ fct_lump_min(.x, min = 10000) %>% factor(ordered = TRUE))) %>%
  # Drop columns with "Prefer not to answer", etc.
  select(-matches("prefer_not_to_answer|do_not_know", ignore.case = TRUE)) %>%
  #Remove NAs 
  na.omit()

DietCorMatrixv2 <- polycor::hetcor(diet_clean,
                                   ML = FALSE,
                          use = "pairwise.complete.obs")




