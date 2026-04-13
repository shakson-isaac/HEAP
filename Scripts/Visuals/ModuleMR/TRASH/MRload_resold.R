library(devtools)
library(TwoSampleMR)
library(ieugwasr)
library(data.table)
library(tidyverse)
library(readr)
library(arrow)
library(ggplot2)
library(genetics.binaRies)
library(pbapply)

#Wait Time: ~40 minutes

# Load config
MRconfig <- fread("/n/groups/patel/shakson_ukb/UK_Biobank/Output/App/Tables/MRpriority_new.csv")

# Base directory
outdirbase <- "/n/groups/patel/shakson_ukb/UK_Biobank/Output/MRcis"

# Function for one iteration
process_MR <- function(i) {
  x <- MRconfig[i, ]
  ExID      <- x[["Exposure"]]
  protID    <- x[["Protein"]]
  diseaseID <- x[["Disease"]]
  
  outdir <- file.path(outdirbase, ExID, protID, diseaseID)
  
  # Check if directory exists
  if (!dir.exists(outdir)) {
    warning(paste("Skipping:", outdir, " — directory does not exist"))
    return(NULL)
  }
  
  files <- list(
    X_to_M = file.path(outdir, "A_X_to_M_results.csv"),
    X_to_Y = file.path(outdir, "B_X_to_Y_results.csv"),
    M_to_Y = file.path(outdir, "C_M_to_Y_results.csv"),
    med    = file.path(outdir, "E_mediation_indirect_product.csv")
  )
  
  # Check if required files exist
  if (!all(file.exists(unlist(files)))) {
    warning(paste("Missing file(s) for:", outdir))
    return(NULL)
  }
  
  # Try reading files
  tryCatch({
    X_to_M <- fread(files$X_to_M)
    X_to_Y <- fread(files$X_to_Y)
    M_to_Y <- fread(files$M_to_Y)
    med <- fread(files$med)
    
    med$exposure <- ExID
    med$protein  <- protID
    med$disease  <- diseaseID
    
    list(XM = X_to_M, XY = X_to_Y, MY = M_to_Y, MED = med)
    
  }, error = function(e) {
    warning(paste("Error reading:", outdir, " —", e$message))
    return(NULL)
  })
}

# Run with progress bar (parallel if multicore available)
results <- pblapply(1:nrow(MRconfig), process_MR, cl = 1)

# Extract results safely
XM  <- lapply(results, function(r) if (!is.null(r)) r$XM)
XY  <- lapply(results, function(r) if (!is.null(r)) r$XY)
MY  <- lapply(results, function(r) if (!is.null(r)) r$MY)
MED <- lapply(results, function(r) if (!is.null(r)) r$MED)

# Combine all available results
XMdf  <- rbindlist(XM,  use.names = TRUE, fill = TRUE)
XYdf  <- rbindlist(XY,  use.names = TRUE, fill = TRUE)
MYdf  <- rbindlist(MY,  use.names = TRUE, fill = TRUE)
MEDdf <- rbindlist(MED, use.names = TRUE, fill = TRUE)


library(dplyr)

method_priority <- c(
  "Inverse variance weighted",
  "Wald ratio",
  "Weighted median",
  "MR Egger",
  "Weighted mode"
)

XMfin <- unique(XMdf)
# XMfin <- XMfin  %>%
#   mutate(method = as.character(method)) %>%
#   group_by(id.exposure, id.outcome) %>%
#   mutate(method_rank = match(method, method_priority)) %>%
#   filter(method_rank == min(method_rank, na.rm = TRUE)) %>%
#   slice_min(pval, with_ties = FALSE) %>%   # if multiple rows of same method, take smallest p
#   ungroup() %>%
#   mutate(
#     BHpval   = p.adjust(pval, method = "BH"),
#     Bonferroni = p.adjust(pval, method = "bonferroni")
#   )


MYfin <- unique(MYdf)
# MYfin <- MYfin  %>%
#   mutate(method = as.character(method)) %>%
#   group_by(id.exposure, id.outcome) %>%
#   mutate(method_rank = match(method, method_priority)) %>%
#   filter(method_rank == min(method_rank, na.rm = TRUE)) %>%
#   slice_min(pval, with_ties = FALSE) %>%   # if multiple rows of same method, take smallest p
#   ungroup() %>%
#   mutate(
#     BHpval   = p.adjust(pval, method = "BH"),
#     Bonferroni = p.adjust(pval, method = "bonferroni")
#   )

XYfin <- unique(XYdf)
# XYfin <- XYfin  %>%
#   mutate(method = as.character(method)) %>%
#   group_by(id.exposure, id.outcome) %>%
#   mutate(method_rank = match(method, method_priority)) %>%
#   filter(method_rank == min(method_rank, na.rm = TRUE)) %>%
#   slice_min(pval, with_ties = FALSE) %>%   # if multiple rows of same method, take smallest p
#   ungroup() %>%
#   mutate(
#     BHpval   = p.adjust(pval, method = "BH"),
#     Bonferroni = p.adjust(pval, method = "bonferroni")
#   )

#'*THE ABOVE IS INCORRECT NEED TO SPLIT THE bonferonni with TYPE of MR 'test'*

#'*TODO*
fwrite(XMfin, file ="/n/groups/patel/shakson_ukb/UK_Biobank/Output/MRres/XMres.csv")
fwrite(XYfin, file ="/n/groups/patel/shakson_ukb/UK_Biobank/Output/MRres/XYres.csv")
fwrite(MYfin, file ="/n/groups/patel/shakson_ukb/UK_Biobank/Output/MRres/cisMYres.csv")
fwrite(MEDdf, file ="/n/groups/patel/shakson_ukb/UK_Biobank/Output/MRres/MEDres.csv")


# THINGS IN TERMS OF EFFICIENCY:
#1.) I realized I ran the MY associations for each time there was a XM hit. 
# SAME occurs with X-->Y. I think the safest redo of MR is to split it into 3 groups
# M --> Y any significant indirect effect
# X --> Y any exposure-protein-disease triplet that exists. But limit to X-->Y to test
# X --> M any exposure-protein pairs that are found fromt he significant indirect effects

# This is overdoing it and slowing speeds
#2.) There is a mixture of wald ratio, IVW, etc. different MR methods happening which is interesting
#3.) Not sure how much to expand and minimize the actual MR screens.
#4.) Make sure the fwrite more useful tables than above^^
# Need bonferroni, FDR corrected p-values? I guess --> And unique entries
#5.) I see an issue with the MR
#'*MR for no cis-variant proteins -- completely doesn't consider the results for them*

# Calculate Percentage:








#Get the indirect product info:
MEDfin <- MEDdf %>% filter(grepl("indirect",component))
