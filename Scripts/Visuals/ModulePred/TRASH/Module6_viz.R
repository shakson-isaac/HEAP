#Libraries:
library(data.table)
library(tidyverse)
library(ggplot2)


# Read in specific files:
folder_name <- "/n/groups/patel/shakson_ukb/UK_Biobank/Data/Parallel/PES_Option1_fastv2/"
Viz <- list.files(folder_name)


cox_files <- Viz[grepl("Cox4_Type1",Viz)]
foldmetrics <- Viz[grepl("FoldMetrics.tsv",Viz)]
overallmetrics <- Viz[grepl("OverallMetrics.tsv",Viz)]


cox_res <- lapply(cox_files, function(x){
  x <- fread(paste0(folder_name,x))
})

fold_res <- lapply(foldmetrics, function(x){
  x <- fread(paste0(folder_name,x))
})

overall_res <- lapply(overallmetrics, function(x){
  x <- fread(paste0(folder_name,x))
})

#Combine cox_res
cox_df <- do.call("rbind",cox_res)









