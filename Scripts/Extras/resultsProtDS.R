library(data.table)
library(tidyverse)

#'*Below is for Univariate E-Prot Links*
#### LOADING UNIVARIATE RESULTS #####
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

Type1 <- loadAssoc("Type1")
Type2 <- loadAssoc("Type2")
Type3 <- loadAssoc("Type3")
Type4 <- loadAssoc("Type4")
Type5 <- loadAssoc("Type5")
Type6 <- loadAssoc("Type6")
Type7 <- loadAssoc("Type7")

# CovarSpec: All ASSOCIATIONS DataStructure:
CovarSpecList <- list(Type1, Type2, Type3, Type4, Type5, Type6)
names(CovarSpecList) <- c("Type1","Type2","Type3","Type4","Type5","Type6")
