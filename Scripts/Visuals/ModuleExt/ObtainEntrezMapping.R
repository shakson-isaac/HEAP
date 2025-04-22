# Rscript to convert UKB protein IDs to Entrez for Pathway Analysis:

#Detect Overlap: Convert ENTREZ to GeneSymbol
library(data.table)
library(tidyverse)
library(org.Hs.eg.db)
library(biomaRt)

#'*CHANGE useMart to useEnsembl*
# Function: Obtain Olink Protein Names to Entrez IDs map
entrezMap <- function(){
  # OMIC IDs
  omicpredIDs <- fread(file = "/n/groups/patel/shakson_ukb/UK_Biobank/Data/OMICSPRED/UKB_Olink_multi_ancestry_models_val_results_portal.csv")
  
  # Connect to the Ensembl database
  ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  
  #useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  #useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl", mirror = "useast")
  
  # Convert UniProt IDs to Entrez IDs
  map_prot_to_entrez <- function(uniprot_ids) {
    converted_ids <- getBM(attributes = c("hgnc_symbol", "entrezgene_id"),
                           filters = "hgnc_symbol",
                           values = uniprot_ids,
                           mart = ensembl)
    return(converted_ids)
  }
  
  # Split the protein names into separate genes - if needed:
  df_split <- omicpredIDs %>%
    mutate(genes = str_split(Gene, "_")) %>%
    unnest(genes)
  
  # Apply the mapping function to each gene
  entrez_results = map_prot_to_entrez(df_split$genes)
  
  # Combine results into a single data frame
  df_entrez <- df_split %>%
    left_join(entrez_results, by = c("genes" = "hgnc_symbol"))
  
  # Identify UniProt IDs that have not been converted
  OlinkIDs_noconvert <- df_entrez[is.na(df_entrez$entrezgene_id),]$Gene
  
  return(df_entrez)
}
UKBentrez <- entrezMap()

fwrite(UKBentrez, file = "/n/groups/patel/shakson_ukb/UK_Biobank/Data/GeneSets/OlinkEntrezConv.txt")

#Utility across other scripts is:
#df_entrez <- fread(file = "/n/groups/patel/shakson_ukb/UK_Biobank/Data/GeneSets/OlinkEntrezConv.txt")

