
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
library(stringr)
library(tidyverse)
library(data.table)
library(ukbwranglr)
library(arrow)
library(furrr)
library(future)
library(purrr)
library(fastDummies)

#'#*Develop a Try-Catch Block System*
#'DEFINE list of global variables (ex. projID) stays the same throughout

#### FUNCTIONS:
# Load Information for Given ProjID and PathID
load_project <- function(projID){
  # Directory Info from find_shared_variables.R (IGLOO RAW, legacy fallback)
  directoryInfo <<- fread(file = heap_raw_or_legacy(
    c("paths", "directoryInfo.txt"),
    legacy_ukb_path("Data", "Paths", "directoryInfo.txt")))

  # Dictionary of PathID to PathName (IGLOO RAW, legacy fallback)
  pathIDs <<- fread(file = heap_raw_or_legacy(
    c("paths", projID, "PathsID.txt"),
    legacy_ukb_path("Data", "Paths", projID, "PathsID.txt")))
  
  # All data related to projID
  data_files <<- directoryInfo[directoryInfo$projID == projID, ]$name
}

load_path <- function(projID, pathID){
  #'* Tip double arrow <<- means return to global env*
  # First path in projID (IGLOO RAW, legacy fallback)
  path <<- fread(file = heap_raw_or_legacy(
    c("paths", projID, paste0("paths_", pathID, ".txt")),
    legacy_ukb_path("Data", "Paths", projID, paste0("paths_", pathID, ".txt"))))
}

# Load Path Data
load_path_data <- function(projBinID, ncores){ #Assign different naming then projId to avoid confusion.
  # Get feather IDs
  featherIDs <<- directoryInfo[directoryInfo$projID == projBinID,]$name
  
  availableCores()
  plan(multisession, workers = ncores) 
  
  # Obtain dataframes for the specific UKB Category:
  UKB_df <- list()
  for(i in featherIDs){
    UKB_df[[i]] <- arrow::read_feather(
      heap_raw_or_legacy(paste0(i, ".feather"),
        legacy_ukb_path("Data", "Raw", paste0(i, ".feather"))),
      col_select = all_of(c("eid",path[grep(i, path$ProjectID),]$descriptive_colnames)
      ))
  }
  #'*Optimize loading of feather files...*
  
  # Merge the dataframes from separate files into 1:
  df <- UKB_df %>% reduce(full_join)
  
  return(df)
}

# Obtain info of most recent files:
order_pathInfo <- function(path, directoryInfo = directoryInfo){
  #Find the set of all IDs in path:
  setIDs <- unique(path$ProjectID)
  allIDs <- gsub("\\(", "", setIDs)
  allIDs <- gsub("\\)", "", allIDs)
  # Concatenate into a single string and split by comma
  all_ids <- unlist(strsplit(paste(allIDs, collapse = ","), ","))
  # Get unique set of IDs
  unique_ids <- unique(all_ids)
  
  #Order these IDs based on modification time:
  pathInfo <- directoryInfo[directoryInfo$name %in% unique_ids, ]
  pathInfo <- pathInfo[order(pathInfo$modifDate, decreasing = T), ]
  
  return(pathInfo)
}

# Load data in order of most recent files:
fast_dataloader <- function(path, directoryInfo = directoryInfo){
  #Concept: Make sure all variables are accounted for:
  #THIS WOULD BE A WHILE LOOP WHERE WE STOP WHEN NO MORE FIELDS ARE ACCOUNTED FOR!
  
  #Initialization:
  pathInfo <- order_pathInfo(path, directoryInfo)
  UKB_df <- list()
  fields <- path$descriptive_colnames
  
  #While loop
  c = 1
  while(length(fields) > 0){
    #Current directory:
    ukb_dir <- pathInfo$name[c]
    
    #Update path to existing fields left:
    path <- path[path$descriptive_colnames %in% fields, ]
    ukb_fields <- path[grepl(pathInfo$name[c], path$ProjectID), ]$descriptive_colnames
    #print(ukb_fields)
    
    UKB_df[[ukb_dir]] <- arrow::read_parquet(
      heap_raw_or_legacy(paste0(ukb_dir, ".parquet"),
        legacy_ukb_path("Data", "Raw", paste0(ukb_dir, ".parquet"))),
      col_select = all_of(c("eid",ukb_fields)
      ))
    fields <- setdiff(fields, ukb_fields)
    c = c + 1
  }
  
  #Merge the dataframes:
  df <- UKB_df %>% reduce(full_join)
  #'*Convert Empty Cells to NA*
  df[df == "",] <- NA
  
  # Remove rows that have missingness along all columns:
  df <- df[
    rowSums(is.na(df)) < (ncol(df) - 1), 
  ]
  
  return(df)
}

# Load data in order of most recent files: using field IDs
fast_dataloader_viafield <- function(allpath, UKBfieldIDs, directoryInfo = directoryInfo){
  #Concept: Make sure all variables are accounted for:
  #THIS WOULD BE A WHILE LOOP WHERE WE STOP WHEN NO MORE FIELDS ARE ACCOUNTED FOR!
  
  #Initialization:
  path <- allpath[allpath$FieldID %in% UKBfieldIDs, ] #'*Subset early on via fieldIDs listed (ONLY alteration DONE);  And provide full dictionary of paths*
  
  pathInfo <- order_pathInfo(path, directoryInfo)
  UKB_df <- list()
  fields <- path$descriptive_colnames
  
  #While loop
  c = 1
  while(length(fields) > 0){
    #Current directory:
    ukb_dir <- pathInfo$name[c]
    
    #Update path to existing fields left:
    path <- path[path$descriptive_colnames %in% fields, ]
    ukb_fields <- path[grepl(pathInfo$name[c], path$ProjectID), ]$descriptive_colnames
    #print(ukb_fields)
    
    UKB_df[[ukb_dir]] <- arrow::read_parquet(
      heap_raw_or_legacy(paste0(ukb_dir, ".parquet"),
        legacy_ukb_path("Data", "Raw", paste0(ukb_dir, ".parquet"))),
      col_select = all_of(c("eid",ukb_fields)
      ))
    fields <- setdiff(fields, ukb_fields)
    c = c + 1
  }
  
  #Merge the dataframes:
  df <- UKB_df %>% reduce(full_join)
  #'*Convert Empty Cells to NA*
  df[df == "",] <- NA
  
  # Remove rows that have missingness along all columns:
  df <- df[
    rowSums(is.na(df)) < (ncol(df) - 1), 
  ]
  
  return(df)
}


# Function to Obtain Data via Instances (Defined by UKB)
UKB_instances <- function(df, instance_expr){
  #"_0_", "_1_" "_2_" "_3_"
  #'*THIS Specifically looks at instance identifier in string*
  multiarray_instanceid <<- instance_expr #Global variable for multiarray handler
  instance_expr <- paste0("f\\d+",instance_expr,"\\d+$")
  
  instance_names <- colnames(df)[grepl(instance_expr,colnames(df))]
  df_inst <- subset(df, select = c("eid",instance_names))
  #print(colSums(is.na(df_inst)))
  
  final_df <- df_inst[
    rowSums(is.na(df_inst)) < (ncol(df_inst) - 1), 
  ]
  #print(colSums(is.na(final_df)))
  
  return(final_df)
}


# Function to Combine Variables that are multiarray:
UKB_multiarray_handle <- function(df){
  # Find variables with multiple arrays:
  UKB_fieldIDs <- gsub("_[^_]*_[^_]*$", "",colnames(df))
  fieldIDs_multiarray <- unique(UKB_fieldIDs[duplicated(UKB_fieldIDs)])
  
  for(i in fieldIDs_multiarray){
    # Combine columns with multiarray via semicolon
    cols <- colnames(df)[grepl(i, colnames(df))]
    
    # Rename column and combine the arrays together
    rename_col <- paste0(i, multiarray_instanceid ,"0.multi")
      #paste0(i, "_0_0.multi")
    df <- df %>% 
      unite("temp_col", all_of(cols), sep = ";", remove = TRUE, na.rm = TRUE) %>%
      rename_with(~rename_col, "temp_col")
    
    # Convert blanks to NA in the specified column
    df[[rename_col]] <- ifelse(df[[rename_col]] == "", NA, df[[rename_col]])
    
  }
  
  return(df)
}
# MAKE SURE TO NOT FILTER OUT COLUMNS IF HAVE MULTIPLE ARRAYS:
# EX. "2_0", "2_1' "2_2" "2_3"


# Function to Filter Out Variables by Missingness
UKB_filter_data <- function(df, missingness){
  # Make sure df is of class data.frameL
  df <- as.data.frame(df)
  
  # Threshold for Missingness (ex. 50%)
  missing_threshold <- missingness
  
  # Keep columns associated with low missingness
  keep_columns <- colSums(is.na(df))/nrow(df) < missing_threshold
  print(keep_columns)
  
  # Filter columns based on missigness
  filtered_data <- df[,keep_columns]
  
  print("DONE")
  return(filtered_data)
}


#'*Convert Categorical Variables and Process for Downstream Analysis*
"difficulty_not_smoking_for_1_day_f3476_0_0"
#vs
"time_from_waking_to_first_cigarette_f3466_1_0"

# Load dataframes that obtain descriptors for codings in a path
load_path_codings <- function(df){
  #Make dataframe with column names and associated category:
  fID <- gsub("\\..*$", "", colnames(df))
  valueType_df <- path[path$descriptive_colnames %in% fID,]
  valueType_df <<- subset(valueType_df, select = c("descriptive_colnames","ValueType","Coding","Notes"))
  
  # UKB Data Codings
  #ukb_datacodings <<- read_csv(url("https://biobank.ctsu.ox.ac.uk/~bbdatan/Codings.csv")) ## Data Codings file (Categorical Variables)
  #'*UPDATE: UKB no longer supports the file above*
  ukb_datacodings <<- read_csv(file = heap_raw_or_legacy(
    c("codings", "Codings.csv"),
    legacy_ukb_path("RScripts", "Extract_Raw", "Finalized", "Codings.csv")))
  
  # UKB Data Showcase Dictionary (NOT NEEDED below)
  #ukb_dictionary <<- read_csv(url("https://biobank.ctsu.ox.ac.uk/~bbdatan/Data_Dictionary_Showcase.csv")) ## dictionary file
  #ukb_dictionary <<- read_csv(file = "/n/groups/patel/shakson_ukb/UK_Biobank/RScripts/Extract_Raw/Finalized/Data_Dictionary_Showcase.csv")

}

# Define integer coding function
integer_coding <- function(x, col_name) {
  #Show Name of Columns edited:
  print(col_name)
  
  #Missingness for Coding: 100373
  x <- ifelse(x == -1,NA, x)
  x <- ifelse(x == -3,NA, x)
  
  #Lessthan1:
  x <- ifelse(x == -10, 0.5, x)
  
  return(x)
}

# Define ordinal coding function:
ordinal_coding <- function(x, col_name, varCode = valueType_df, ukb_coding = ukb_datacodings){
  #Show Name of Columns edited:
  print(col_name)
  
  # Obtain coding ID
  codingID <- varCode[varCode$descriptive_colnames == col_name, ]$Coding
  
  # Obtain coding Info to Convert
  codes <- ukb_coding[ukb_coding$Coding == codingID,]
  
  # Swap strings for integers for ordinal categorical variables
  x <- codes$Value[match(x, codes$Meaning)]
  
  #Turn Negative codings into NA: (Do NOT Know, Prefer Not the Answer)
  x <- ifelse(x < 0, NA, x)
  
  x <- as.numeric(x)
  
  return(x)
}

#'* Develop a system to split categorical variables into ordinal and unorder*
#'* Provide different ways of handling missingness*


##'*Integer Coding*
UKB_integer_handle <- function(df){
  integer_cols <- valueType_df[grepl("Integer",valueType_df$ValueType),]$descriptive_colnames
  
  #check the categories are in the dataframe
  transform_cols <- integer_cols[integer_cols %in% colnames(df)]
  
  if(length(transform_cols) > 0){
    # Use the functions across different types of columns.
    df <- df %>%
      mutate(across(all_of(transform_cols),
                    ~ integer_coding(., col_name = cur_column())))
  }
  
  print("DONE2")
  return(df)
}

##'*Ordinal Coding*
UKB_ordinal_handle <- function(df, ordinal_codings){
  #ordinal_codings <- c(100377, 100394, 100401)
  if(length(ordinal_codings) > 0){
    ordinal_cols <- valueType_df[valueType_df$Coding %in% ordinal_codings, ]$descriptive_colnames
    
    #check categories in dataframe.
    transform_cols <- ordinal_cols[ordinal_cols %in% colnames(df)]
    
    df <- df %>%
      mutate(across(all_of(transform_cols),
                    ~ ordinal_coding(., col_name = cur_column())))
  }
  
  print("DONE3")
  return(df)
}


##'*Categorical Encoding (Unordered) - OneHot*

# Function to OneHot Encode Rest of Categorical Variables
UKB_onehot_handle <- function(df){
  #rest of columns that are strings (characters) are categorical
  categorical_cols <- df %>% 
    select_if( function(col) { is.character(col) | is.factor(col)}) %>% 
    names()
  print(categorical_cols)
  
  if(length(categorical_cols) > 0){ #If categorical variables are remaining one-hot encode:
    # One-hot encode the Categories column
    #(Do NOT Know, Prefer Not the Answer) are also one_hot encoded columns
    df_onehot <- dummy_cols(df, select_columns = categorical_cols, 
                            split = ";", remove_selected_columns = TRUE,
                            ignore_na = TRUE)
  } else{
    df_onehot <- df
  }
  
  print("DONE6")
  return(df_onehot)
}

# Function to Factor Encode Rest of Categorical Variables
UKB_factor_handle <- function(df){
  #rest of columns that are strings (characters) are categorical
  categorical_cols <- df %>% select_if(is.character) %>% names()
  
  df[categorical_cols] <- lapply(df[categorical_cols], factor)
  
  return(df)
}


# Function to scale ordinal variables from 0 to max categories
scale_ordinal_columns <- function(df, ordinal_codings){
  ordinal_cols <- valueType_df[valueType_df$Coding %in% ordinal_codings, ]$descriptive_colnames
  
  if(length(ordinal_cols) > 0){
    df[ordinal_cols] <- lapply(df[ordinal_cols], function(x){
      unique_vals <- unique(x)
      sort_vals <- sort(unique_vals)
      rescale_vals <- match(x, sort_vals) - 1 #start the scale from 0
    })
  }
  
  print("DONE4")
  return(df)
}

# Function identify and encoding binary variables across all variables:
UKB_binary_handle <- function(df){
  #Check for binary columns
  binary_columns <- sapply(df, function(col){
    unique_vals <- unique(col)
    #'*Make sure to (DELETE Columns that are negative values and have Do NOT KNOW etc.)*
    #Make:
    #Prefer not to answer
    #Do not know
    #etc. into NA
    #OR SWITCH THE CODING NAME --> NUMERIC --> THEN COUNT only above and including 0.
    #here NA appears as a unique value. (how to remove NA)
    
    #1st Thing (Remove NAs):
    unique_vals <- unique_vals[!is.na(unique_vals)]
    #print(unique_vals)
    
    #length(unique_vals > 0)
    length(unique_vals) == 2
  })
  #print(binary_columns)
  
  #Apply transformation to binary columns
  cols_to_binarize <- names(binary_columns[binary_columns])
  print(cols_to_binarize)
  
  if(length(cols_to_binarize) > 0){
    df[cols_to_binarize] <- lapply(df[cols_to_binarize], function(col){
      binary <- as.integer(factor(col)) - 1 #, levels = unique_vals
    })
  }
  
  print("DONE5")
  return(df)
}

codings <- function(projID, pathID, df){
  load_path(projID, pathID)
  load_path_codings(df)
  
  # Annotated UKB codings:
  human_anno_codings <- read.table(file = heap_raw_or_legacy(
                                     c("codings", "ukb_codings_humananno_ver2.csv"),
                                     legacy_ukb_path("RScripts", "Extract_Raw", "ukb_codings_humananno_ver2.csv")),
                                   sep = ",", header = T, row.names = 1)
  #^Reads "T" and "F" and True and False
  
  # Obtain info on whether there are ordinal codes:
  codings_path <<- human_anno_codings[human_anno_codings$data_codes %in% valueType_df$Coding,]
  ord_code <<- codings_path$data_codes[codings_path$ordinal]
}

