library(data.table)
library(tidyverse)

#Get the ASGR1 variants (Used in MR)
ASGR1_decode <- fread("/n/groups/patel/shakson_ukb/UK_Biobank/Data/MR/ProteinInst_DECODE/protein_ASGR1_cis_clumped.tsv")

ASGR1_UKB <- fread("/n/groups/patel/shakson_ukb/UK_Biobank/Data/MR/ProteinInst/protein_ASGR1_cis_clumped.tsv")