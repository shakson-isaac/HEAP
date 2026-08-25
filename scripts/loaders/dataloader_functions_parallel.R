
local({
  candidates <- c(
    Sys.getenv("HEAP_PATHS_FILE", unset = ""),
    file.path(getwd(), "workflow", "00_paths.R"),
    file.path(getwd(), "..", "workflow", "00_paths.R"),
    file.path(getwd(), "..", "..", "workflow", "00_paths.R"),
    file.path(getwd(), "..", "..", "..", "workflow", "00_paths.R")
  )
  candidates <- candidates[nzchar(candidates)]
  hit <- candidates[file.exists(candidates)][1]
  if (!is.na(hit)) source(hit)
})
#Overall List:
library(stringr)
library(tidyverse)
library(data.table)
library(ukbwranglr)
library(arrow)
library(furrr)
library(future)
library(purrr)
library(fastDummies)

#File Specific
library(future.apply)
library(future)
source("./dataloader_functions_upd.R")
#Globals
#projID
#directoryInfo
#CHECK dataloader_functions_upd.R for full list of globals.


# Function to load and process a single path:
load_and_process_path <- function(path_id, proj_id, directory_info) {
  path <- load_path(proj_id, path_id)
  result <- fast_dataloader(path, directory_info)
  return(result)
}


# Parallelized load paths:
parallel_UKB_load <- function(pathIDs){
  # Parallel load features: Combine the results into a list 
  features <- future_lapply(pathIDs, function(i) {
    load_and_process_path(i, projID, directoryInfo)
  }, future.seed = TRUE)
  
  # Rename list of dataframes based on pathID:
  name_pathIDs <- paste0("path_",pathIDs)
  features <- setNames(features, name_pathIDs)
  
  return(features)
}


# Function to preprocess a single path:
preprocess_UKB_df<- function(path_id, missingness, timepoint, df, feature_engineer){
  #Setup
  df <- UKB_instances(df, timepoint)
  df <- UKB_multiarray_handle(df)
  
  #'*CHECK the FILTERING STEP SOMETIMES DOESNT WORK WELL!!! especially metabolomics data*
  df <- UKB_filter_data(df, missingness) #I personally think this should be 70 or 80% complete!


  if(feature_engineer){
    codings(projID, path_id, df)
    #Feature Engineer the Specific Variable Types:
    df <- UKB_integer_handle(df)
    df <- UKB_ordinal_handle(df, ordinal_codings = ord_code)
    df <- scale_ordinal_columns(df, ord_code)
    df <- UKB_binary_handle(df)
    df <- UKB_onehot_handle(df)
  }
  
  return(df)
}


# Parallelized processing paths:
parallel_UKB_preprocess <- function(pathIDs, features, missingness, instance, feature_engineer){
  #missingness: perc complete  (ex. allow up to 0.3 missing)
  #instance: "_0_" is instance 0 for UKB
  
  features_processed <- future_lapply(pathIDs, function(i) {
    name_pathIDs <- paste0("path_",i) #this got updates:
    df <- features[[name_pathIDs]]
    
    preprocess_UKB_df(i, missingness, instance, df, feature_engineer)
  }, future.seed = TRUE)
  
  name_pathIDs <- paste0("path_",pathIDs)
  features_processed <- setNames(features_processed, name_pathIDs)
  
  return(features_processed)
}


# Function for writing Parquet files
# NOTE: This per-column parquet cache is NOT used by HEAP_loader.R (which reads
# whole-file parquet via fast_dataloader*). It is an optional preprocessed cache.
# Output retargeted off the legacy UK_Biobank tree to an IGLOO HEAP intermediate
# location so nothing is written under legacy paths.
write_parquet_file <- function(pathIDs, df) {
  outdir <- heap_project_intermediate("ukb_columns", projID)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  fileout <- file.path(outdir, paste0(pathIDs, ".parquet"))
  arrow::write_parquet(df, sink = fileout)
}


# Parallelized Write Parquet Files:
parallel_write_parquet <- function(pathIDs, features){
  # Use furrr to parallelize the loop
  name_pathIDs <- paste0("path_",pathIDs)
  future_map(name_pathIDs, ~write_parquet_file(.x, features[[.x]]))
}
