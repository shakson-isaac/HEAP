library(data.table)
library(tidyverse)
# Extract PGS scores from scratch:
#
#
#

#How to extract files with output of .sscore
dir <- "/n/scratch/users/s/shi872/UKB_intermediate/Genetics/UKBProtGS/"
x <- list.files(path = dir,
           pattern = "\\.sscore$")
omicpredIDs <- fread(file = "/n/groups/patel/shakson_ukb/UK_Biobank/Data/OMICSPRED/UKB_Olink_multi_ancestry_models_val_results_portal.csv")


#Obtain matrix of Proteomics Genetic Scores from: https://www.omicspred.org/Scores/Olink/UKB_MULTI

UKBprotGS <- NULL
for(i in x){
  #Load PGS score for each omic file
  df <- fread(file = paste0(dir, i))
  
  #ID Conversions
  opID <- gsub("\\.sscore$", "", i)
  protID <- omicpredIDs$Gene[match(opID,omicpredIDs$OMICSPRED_ID)]
  
  #Change Colnames
  colnames(df) <- c("eid", protID)
  
  # Add to df
  if (is.null(UKBprotGS)) {
    # For the first data frame, just assign it to UKBprotGS
    UKBprotGS <- df
  } else {
    # Merge by 'eid'
    UKBprotGS <- merge(UKBprotGS, df, by = "eid", all = TRUE)
  }
}

#FIX THIS ERROR TOMORROW:
#Error in merge.data.table(UKBprotGS, df, by = "eid", all = TRUE) : 
#x has some duplicated column name(s): IDO1.x,IDO1.y. Please remove or rename the duplicate(s) and try again.
#In addition: Warning message:
#  In merge.data.table(UKBprotGS, df, by = "eid", all = TRUE) :
#  column names 'IDO1.x', 'IDO1.y' are duplicated in the result


