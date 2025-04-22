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

setwd("/n/groups/patel/shakson_ukb/UK_Biobank/")

# FUNCTION: Load Files w/ Progress Bar:
loadUniAssoc <- function(covarType){
  error_files <- list()  # Track files w/ errors
  stat_all <- list()
  train_data_list <- list()
  test_data_list <- list()
  
  # For loop through stored Univariate Associations in batches (400 batches here)
  for(i in 1:400){
    tryCatch({
      load <- readRDS(paste0("/n/groups/patel/shakson_ukb/UK_Biobank/Data/Parallel/Prot_PGSPXS_univar/",covarType,"/univar_assoc_",i,".rds"))
      
      train_data_list[[i]] <- load$train
      test_data_list[[i]]  <- load$test
      
    }, error = function(e) {
      message(paste("Error reading file:", i))
      error_files <<- append(error_files, i)
    })
  }
  
  # Store results
  stat_all[["train"]] <- train_data_list
  stat_all[["test"]] <- test_data_list
  
  return(stat_all)
}

# FUNCTION: List of list of dataframes to one list of dataframes:
combine_dataframes <- function(list_of_lists) {
  # Determine the number of elements in each inner list
  num_dfs <- length(list_of_lists[[1]])
  
  # For each position in the inner lists, bind the rows across all lists
  combined_dfs <- map(1:num_dfs, function(i) {
    bind_rows(map(list_of_lists, ~ .x[[i]]))
  })
  
  return(combined_dfs)
}

#'*LOADING Associations Across Specifications*
loadAssoc <- function(covarType){
  # Load specific Covar Spec:
  assoc <- loadUniAssoc(covarType)
  
  # Aggregating Train and Test Results:
  assoc$train <- combine_dataframes(assoc$train)
  assoc$test <- combine_dataframes(assoc$test)
  
  return(assoc)
}

covarTypes <- c("Type1", "Type2", "Type3", "Type4", "Type5", "Type6") #Specify Types/Models utilized:
CovarSpecList <- lapply(covarTypes, loadAssoc)
names(CovarSpecList) <- covarTypes
#Renaming:
names(CovarSpecList) <- c("Model1","Model2","Model3","Model4","Model5", "Model6")

#'*Only significant associations*
# FUNCTION: Obtain DF of replicated stats between Train and Test
replication_stats_final <- function(train, test){
  sigRES <- list()
  
  train_assoc <- train
  test_assoc <- test
  
  colnames(train_assoc)[!(colnames(train_assoc) %in% c("ID","omicID"))] <-  gsub(" ", "",
                                                                                 paste0(colnames(train_assoc)[!(colnames(train_assoc) %in% c("ID","omicID"))], "_train"))
  
  colnames(test_assoc)[!(colnames(test_assoc) %in% c("ID","omicID"))] <-  gsub(" ", "",
                                                                               paste0(colnames(test_assoc)[!(colnames(test_assoc) %in% c("ID","omicID"))], "_test"))
  
  
  all_assoc <- list(train_assoc, test_assoc) %>% reduce(full_join)
  
  #Add unique ID for each Exposure-Protein Pair:
  all_assoc$AssocID <- paste0(all_assoc$omicID,":", all_assoc$ID)
  
  
  
  pval_thresh <- 0.05/nrow(all_assoc) #Standard Bonferroni Correction for All ASSOCIATIONS. Ex. for 200 exposures, 3000 proteins: 0.05/(200*3000)
  
  allAssoc_significant <- all_assoc %>%
    filter(`Pr(>|t|)_test` < pval_thresh &
             `Pr(>|t|)_train` < pval_thresh)
  
  trainAssoc_significant <- all_assoc %>%
    filter(`Pr(>|t|)_train` < pval_thresh)
  
  testAssoc_significant <- all_assoc %>%
    filter(`Pr(>|t|)_test` < pval_thresh)
  
  sigRES[["sigTrain"]] <- trainAssoc_significant
  sigRES[["sigTest"]] <- testAssoc_significant
  sigRES[["sigBOTH"]] <- allAssoc_significant
  sigRES[["Bonf.correction"]] <- pval_thresh
  
  return(sigRES)
}

# CovarSpec: All Significant Assocations Datastructure
CovarSpecAssocList <- lapply(CovarSpecList, function(x){
  Assoc <- list()
  E <- replication_stats_final(x$train[[1]], x$test[[1]]) #E
  GxE <- replication_stats_final(x$train[[2]], x$test[[2]]) #GxE
  
  Assoc[["E"]] <- E
  Assoc[["GxE"]] <- GxE
  
  return(Assoc)
})

#'*Create DataStructure*
HEAPconstruct <- setClass(
  "HEAPconstruct",
  slots = c(
    HEAPlist = "list", #List: All association results
    HEAPsig = "list" #List: Replicated results (Significant in both Train and Test)
  )
)

#Create the datastructure to utilize:
HEAPassoc <- HEAPconstruct(
  HEAPlist = CovarSpecList,
  HEAPsig = CovarSpecAssocList
)

#Save RDS file of PXSloader object
gc()
class(HEAPassoc)
qsave(HEAPassoc, file = "./Output/HEAPres/HEAPassoc.qs")
