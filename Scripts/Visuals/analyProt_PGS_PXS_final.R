library(data.table)
library(tidyverse)
library(ggplot2)
library(ggpmisc)
library(ggrepel)
library(gridExtra)
library(patchwork)


setwd("/n/groups/patel/shakson_ukb/UK_Biobank/")

#'*Create HEAPres Data Structure: Encompassing All Covariate Specifications*
loadPXSmetrics <- function(covarType){
  PXSmetrics <- list()
  
  Lasso <- data.frame()
  R2train <- data.frame()
  R2test <- data.frame()
  R2trainCat <- data.frame()
  R2testCat <- data.frame()
  error_files <- c()
  
  for(i in 1:400){
    tryCatch({
      base_path <- paste0("/n/groups/patel/shakson_ukb/UK_Biobank/Data/Parallel/Prot_PGSPXS_final/", covarType, "/")
      
      Lasso_idx <- fread(file=paste0(base_path, "lassofit_", i, ".txt"))
      R2train_idx <- fread(file=paste0(base_path, "R2train_", i, ".txt"))
      R2test_idx <- fread(file=paste0(base_path, "R2test_", i, ".txt"))
      R2trainCat_idx <- fread(file=paste0(base_path, "R2trainCat_", i, ".txt"))
      R2testCat_idx <- fread(file=paste0(base_path, "R2testCat_", i, ".txt"))
      
      Lasso <- rbind(Lasso, Lasso_idx)
      R2train <- rbind(R2train, R2train_idx)
      R2test <- rbind(R2test, R2test_idx)
      R2trainCat <- rbind(R2trainCat, R2trainCat_idx)
      R2testCat <- rbind(R2testCat, R2testCat_idx)
      
    }, error = function(e) {
      message(paste("Error reading file:", i))
      error_files <<- append(error_files, i)
    })
  }
  
  print(error_files)
  PXSmetrics[["Lasso"]] <- Lasso
  PXSmetrics[["R2train"]] <- R2train
  PXSmetrics[["R2test"]] <- R2test
  PXSmetrics[["R2trainCat"]] <- R2trainCat
  PXSmetrics[["R2testCat"]] <- R2testCat
  
  return(PXSmetrics)
}
covarTypes <- c("Type1", "Type2", "Type3", "Type4", "Type5") #Specify Types/Models utilized:
HEAPres <- lapply(covarTypes, loadPXSmetrics)
names(HEAPres) <- covarTypes

#Try this for naming:
names(HEAPres) <- c("Model1","Model2","Model3","Model4","Model5")

#'*Create DataStructure*
GvEconstruct <- setClass(
  "GvEconstruct",
  slots = c(
    HEAP = "list", #List: R2 results with Train v Test.
    R2ge = "data.frame", #Dataframe: R2 Test
    R2geCI = "data.frame", #Dataframe: R2 Test CI
    R2catCI = "data.frame", #Dataframe: R2 Test Category CI
    R2allCI = "data.frame", # Dataframe: R2 Test All CI
    Ecat = "character", #E Categories
    covarType = "character" #Covariate Specs (Models)
  )
)


#'*Combine Specifications into 1 Dataframe*
# Define a function to compute summary statistics
get_CI <- function(column){
  stat <- data.frame()
  if(is.character(column)){
    c <- unique(column)
    return(c)
  } else{
    stat <- tibble(
      mean = mean(column),
      sd = sd(column),
      n = length(column),
      se = sd / sqrt(n),
      ci_margin = qt(0.975, df = n - 1) * se,
      #ci_low = mean - qt(0.975, df = n - 1) * se, #95% confidence interval:
      #ci_high = mean + qt(0.975, df = n - 1) * se
      ci = paste0(paste0(round(mean, 4), " ± ", round(ci_margin, 4)))
    )
    return(stat$ci)
  }
}

# Function to split string and calculate numeric values
process_column <- function(column) {
  tibble(
    mean = as.numeric(str_extract(column, "^[^ ±]+")),
    ci_margin = as.numeric(str_extract(column, "(?<=± )[^ ]+$")),
    ci_low = mean - ci_margin,
    ci_high = mean + ci_margin
  )
}

# Create one structure for results from covariate specifications:
R2test <- map2_dfr(HEAPres, names(HEAPres), ~ .x$R2test %>%
                     select(all_of(c("ID", "G", "E", "GxE"))) %>%
                     mutate(cID = .y))


#'*Obtain Confidence Intervals*
R2test_CI <- R2test %>%
  group_by(ID, cID) %>% #Mean across test folds
  summarise(mean_G = mean(G),
            mean_E = mean(E), 
            mean_GxE = mean(GxE),
            .groups = 'drop') %>%
  group_by(ID) %>% #get CIs across covar spec
  summarize(across(c("mean_G","mean_E","mean_GxE"), get_CI))

plot_testR2_CI  <- R2test_CI %>%
  rowwise() %>%
  mutate(G_stats = list(process_column(mean_G)),
         E_stats = list(process_column(mean_E)),
         GxE_stats = list(process_column(mean_GxE))) %>%
  # Unnest the statistics into separate columns
  unnest(c(G_stats, E_stats, GxE_stats), names_sep = "_") 


#'*Category Specific Tables:*
categories <- c("Alcohol","Diet_Weekly","Smoking","Exercise_MET",
                "Exercise_Freq","Internet_Usage","Deprivation_Indices","Vitamins")

# purrr::map2 iterate over the list created to generate one dataframe:
R2testCat <- map2_dfr(HEAPres, names(HEAPres), ~ .x$R2testCat %>%
                        select(all_of(c("ID", categories))) %>%
                        mutate(cID = .y))

R2testCat_CI <- R2testCat %>%
  group_by(ID, cID) %>%
  summarise(across(c(categories), mean, .names = "{col}")) %>%  #Mean of each Specification
  mutate(across(c(categories), mean, .names = "mean_{col}")) %>% #Mean across Specification
  mutate(across(c(categories), 
                list(
                  stdErr = ~ sd(.) / sqrt(n()) #, 
                  #ci_low = ~ . - qt(0.975, df = n() - 1) * stdErr, 
                  #ci_high = ~ . + qt(0.975, df = n() - 1) * stdErr
                ), 
                .names = "{col}_{fn}")
  ) %>% #Obtain standard errors
  mutate(
    n = n()
  ) %>%
  mutate(across(starts_with("mean_"),
                list(
                  ci_low = ~ . - qt(0.975, df = n - 1) * get(paste0(sub("mean_", "", cur_column()), "_stdErr")), 
                  ci_high = ~ . + qt(0.975, df = n - 1) * get(paste0(sub("mean_", "", cur_column()), "_stdErr"))
                ),
                .names = "{col}_{fn}")
  ) #obtain CI margins for plotting

plot_testR2Cat_CI <- R2testCat_CI %>%
  select(-c("cID","Alcohol","Diet_Weekly","Smoking","Exercise_MET",
            "Exercise_Freq","Internet_Usage","Deprivation_Indices","Vitamins")) %>%
  unique()


plot_testR2_all <- merge(plot_testR2_CI, plot_testR2Cat_CI, by = "ID")


#Create the datastructure to utilize:
HEAPge <- GvEconstruct(
  HEAP = HEAPres,
  R2ge = R2test, #Dataframe: R2 Test
  R2geCI = plot_testR2_CI, #Dataframe: R2 Test CI
  R2catCI = plot_testR2Cat_CI, #Dataframe: R2 Test Category CI
  R2allCI = plot_testR2_all, # Dataframe: R2 Test All CI
  Ecat = categories, #E Categories
  covarType = names(HEAPres) #Covariate Specs (Models)
)

#Remove unnecessary dataframes:
rm(HEAPres, plot_testR2_all, plot_testR2_CI,
  plot_testR2Cat_CI, R2test, R2test_CI, 
  R2testCat, R2testCat_CI,
  categories,
  covarTypes)

#Save RDS file of PXSloader object
gc()
class(HEAPge)
saveRDS(HEAPge, file = "./Output/HEAPres/HEAPge.rds")


#To put into a datastructure:
#HEAPres --> Generic Data Structure.
#R2test --> Dataframe with all models.
#plot_testR2_CI --> For GvE plot, pathways etc.
#plot_testR2Cat_CI --> Category based GvE plot.
#plot_testR2_all --> Alternative Category based GvE plot --> put in the supplement.
