#Libraries:
library(data.table)
library(tidyverse)
library(ggplot2)


# Read in specific files:
folder_name <- "/n/groups/patel/shakson_ukb/UK_Biobank/Data/Parallel/PES_Option1_fastv2/"

types <- c("Type1","Type2","Type3","Type4","Type5")


loadname <- paste0(folder_name,types[1],"/")
Viz <- list.files(loadname)


cox_files <- Viz[grepl("Cox4",Viz)]
foldmetrics <- Viz[grepl("FoldMetrics.tsv",Viz)]
overallmetrics <- Viz[grepl("OverallMetrics.tsv",Viz)]


cox_res <- lapply(cox_files, function(x){
  x <- fread(paste0(loadname,x))
})

fold_res <- lapply(foldmetrics, function(x){
  x <- fread(paste0(loadname,x))
})

overall_res <- lapply(overallmetrics, function(x){
  x <- fread(paste0(loadname,x))
})