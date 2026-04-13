library(data.table)
library(tidyverse)

x <- list()
for(i in 1:2000){
  x[["DE"]][[i]] <- fread(paste0("/n/groups/patel/shakson_ukb/UK_Biobank/Output/MR_edges/edges_DECODE/logs/DE/chunk_",i,".log"))
  x[["DP"]][[i]] <- fread(paste0("/n/groups/patel/shakson_ukb/UK_Biobank/Output/MR_edges/edges_DECODE/logs/DP/chunk_",i,".log"))
  x[["ED"]][[i]] <- fread(paste0("/n/groups/patel/shakson_ukb/UK_Biobank/Output/MR_edges/edges_DECODE/logs/ED/chunk_",i,".log"))
  x[["EP"]][[i]] <- fread(paste0("/n/groups/patel/shakson_ukb/UK_Biobank/Output/MR_edges/edges_DECODE/logs/EP/chunk_",i,".log"))
  x[["PD"]][[i]] <- fread(paste0("/n/groups/patel/shakson_ukb/UK_Biobank/Output/MR_edges/edges_DECODE/logs/PD/chunk_",i,".log"))
  x[["PE"]][[i]] <- fread(paste0("/n/groups/patel/shakson_ukb/UK_Biobank/Output/MR_edges/edges_DECODE/logs/PE/chunk_",i,".log"))
}


combined_data <- list.files(path = "/n/groups/patel/shakson_ukb/UK_Biobank/Output/MR_edges/edges_DECODE/logs/PD/", pattern = "*.log", full.names = TRUE) %>%
  map_dfr(fread)


library(data.table)
library(purrr)
library(dplyr)

files <- list.files(
  "/n/groups/patel/shakson_ukb/UK_Biobank/Output/MR_edges/edges_DECODE/logs/PD/",
  pattern = "\\.log$",
  full.names = TRUE
)

combined_data <- rbindlist(
  lapply(files, \(f) fread(f, header = FALSE, fill = TRUE)),
  fill = TRUE,
  use.names = FALSE
)

# optional: keep which file each row came from
combined_data[, source_file := rep(files, times = sapply(lapply(files, fread, header=FALSE, fill=TRUE), nrow))]



## Find the failed indices:

library(purrr)

edges <- c("EP","PE","PD","DP","ED","DE")
base  <- "/n/groups/patel/shakson_ukb/UK_Biobank/Output/MR_edges/edges_DECODE/logs"

expected <- 1:2000

missing_by_edge <- set_names(edges) %>%
  map(function(e) {
    files <- list.files(
      file.path(base, e),
      pattern = "\\.log$",
      full.names = TRUE
    )
    
    idx <- as.integer(sub(".*chunk_(\\d+)\\.log$", "\\1", files))
    idx <- idx[!is.na(idx)]
    
    setdiff(expected, idx)
  })

missing_by_edge







##### TRASH #####

edges <- c("EP","PE","PD","DP","ED","DE")


files_PD <- list.files(
  "/n/groups/patel/shakson_ukb/UK_Biobank/Output/MR_edges/edges_DECODE/logs/PD",
  pattern = "\\.log$",
  full.names = TRUE
)
files_PD


# extract indices like 1896 from ".../chunk_1896.log"
idx <- as.integer(sub(".*chunk_(\\d+)\\.log$", "\\1", files_PD))

# expected indices (change 2000 to whatever the true max is)
expected <- 1:2000

missing_idx <- setdiff(expected, idx)
missing_idx










# combined_data <- list.files(
#   path = "/n/groups/patel/shakson_ukb/UK_Biobank/Output/MR_edges/edges_DECODE/logs/PD/",
#   pattern = "\\.log$",
#   full.names = TRUE
# ) %>%
#   map_dfr(~fread(.x, header = FALSE))
# 
# combined_data2 <- list.files(
#   path = "/n/groups/patel/shakson_ukb/UK_Biobank/Output/MR_edges/edges_DECODE/logs/EP/",
#   pattern = "\\.log$",
#   full.names = TRUE
# ) %>%
#   map_dfr(~fread(.x, header = FALSE))
# 
# 
# combined_data$V2[grepl("FUR",combined_data$V2)]