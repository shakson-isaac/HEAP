#Libraries:
library(data.table)
library(tidyverse)
library(ggplot2)
library(ggpmisc)
library(pbapply)
library(dplyr)
library(purrr)
library(plotly)
library(htmlwidgets)
library(qs)

#Library for Pathway Analysis:
library(org.Hs.eg.db)
library(biomaRt)
library(ggridges)
library(ReactomePA)
library(DOSE)
library(clusterProfiler)
library(enrichplot)

#'*LOAD E Associations*
setwd("/n/groups/patel/shakson_ukb/UK_Biobank/")
HEAPassoc <- qread("./Output/HEAPres/HEAPassoc.qs")


### Pathway Analysis: Univariate Associations
# Steps to Follow (For easy analysis)
# Create Genelists
# Import Tissue Specific Sets and Reactome (Most Informative Enrichments)
# Run GSEA

# Function to build genelist
CreateGenelists <- function(id,CovarType,Split){
  GeneLists <- list()
  
  covarSpec <- HEAPassoc@HEAPlist[[CovarType]][[Split]][[1]]
  bonferonni_thresh <- 0.05/nrow(covarSpec)
  
  Eassoc_df <- covarSpec %>% 
    filter(ID == id) #MUST not overlap id and Eid naming.
  
  # No directionality (Add directionality with 'Estimate > 0' or 
  # 'Estimate < 0')
  ORA_all <- Eassoc_df %>%
    filter(`Pr(>|t|)` < bonferonni_thresh) %>%
    pull(omicID) %>%
    unlist(strsplit(., "_")) #split complexes by individual prot.
  
  ORA_up <- Eassoc_df %>%
    filter(`Pr(>|t|)` < bonferonni_thresh & Estimate > 0) %>%
    pull(omicID) %>%
    unlist(strsplit(., "_")) #split complexes by individual prot.
  
  ORA_down <- Eassoc_df %>%
    filter(`Pr(>|t|)` < bonferonni_thresh & Estimate < 0) %>%
    pull(omicID) %>%
    unlist(strsplit(., "_")) #split complexes by individual prot.
  
  
  # USE the t value: Contains significance and direction. Instead of Beta Estimate.
  GSEAgenelist <- Eassoc_df$`t value`
  names(GSEAgenelist) <- Eassoc_df$omicID
  
  #ONLY keep the first name of set of complexes || to prevent ties in GSEA.
  names(GSEAgenelist) <- gsub("_.*", "", names(GSEAgenelist)) 
  GSEAgenelist = sort(GSEAgenelist, decreasing = TRUE)
  
  # Sorted Beta Estimates:
  Betas <- Eassoc_df$Estimate
  names(Betas) <- Eassoc_df$omicID
  
  #ONLY keep the first name of set of complexes || to prevent ties in GSEA.
  names(Betas) <- gsub("_.*", "", names(Betas)) 
  Betas = sort(Betas, decreasing = TRUE)
  
  
  GeneLists[["ORA_all"]] <- ORA_all
  GeneLists[["ORA_up"]] <- ORA_up
  GeneLists[["ORA_down"]] <- ORA_down
  GeneLists[["GSEA"]] <- GSEAgenelist
  GeneLists[["Betas"]] <- Betas
  
  return(GeneLists)
}

# Function to Run GSEA Analysis
PSEA_tissue <- function(GeneList){
  GTEX <- fread("/n/groups/patel/shakson_ukb/UK_Biobank/Data/GeneSets/GTEX_tissue.txt")
  omicpredIDs <- fread(file = "/n/groups/patel/shakson_ukb/UK_Biobank/Data/OMICSPRED/UKB_Olink_multi_ancestry_models_val_results_portal.csv")
  UKB_universe = unlist(strsplit(omicpredIDs$Gene, "_"))
  
  gsea <- GSEA(geneList= GeneList[["GSEA"]], 
               #nPerm=10000, 
               minGSSize = 10, 
               maxGSSize = 500, 
               pvalueCutoff = 0.05, 
               pAdjustMethod = "BH", 
               TERM2GENE=GTEX, 
               TERM2NAME=NA)
  return(gsea@result)
}

# Function: Obtain Olink Protein Names to Entrez IDs map
entrezMap <- function(){
  # OMIC IDs
  omicpredIDs <- fread(file = "/n/groups/patel/shakson_ukb/UK_Biobank/Data/OMICSPRED/UKB_Olink_multi_ancestry_models_val_results_portal.csv")
  
  # Connect to the Ensembl database
  ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  
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

# Function: Convert Genelists to include Entrez IDs
convertEntrez <- function(GeneList){
  
  df_entrez <- entrezMap()
  
  for(i in names(GeneList)){
    name <- paste0(i,"_entrez")
    if(is.null(names(GeneList[[i]]))){ #If not a named character (ORA)
      GeneList[[name]] <- df_entrez[match(GeneList[[i]], 
                                          df_entrez$Gene),]$entrezgene_id
    }else{
      GeneList[[name]] <- GeneList[[i]]
      names(GeneList[[name]]) <- df_entrez[match(names(GeneList[[i]]),
                                                 df_entrez$Gene),]$entrezgene_id
    }
  }
  
  return(GeneList)
}

# Function: GSEA of Pathways (Uses Entrez ID)
PSEA_paths <- function(GeneList){
  EnrichGSEA <- list()
  
  Reactome <- gsePathway(GeneList[["GSEA_entrez"]], 
                         minGSSize = 10, maxGSSize = 500,
                         pvalueCutoff = 0.05, pAdjustMethod = "BH", 
                         verbose = TRUE,
                         eps = 0)

  return(Reactome@result)
}

# Tissue + Reactome: GSEA:
TissueEnrichGSEA <- function(Eid, CType, Tsplit){
  tissue <- list()
  
  PA_SS <- CreateGenelists(id = Eid,
                           CovarType = CType,
                           Split = Tsplit)
  PSEA_res <- PSEA_tissue(PA_SS)
  
  return(PSEA_res)
}
PathwayEnrichGSEA <- function(Eid, CType, Tsplit){
  pathways <- list()
  
  PA_SS <- CreateGenelists(id = Eid,
                           CovarType = CType,
                           Split = Tsplit)
  PA_SS_entrez <- convertEntrez(PA_SS)

  GSEA_enrich <- PSEA_paths(PA_SS_entrez)

  return(GSEA_enrich)
}

df_entrez <- entrezMap()


# Run tissue and pathway enrichment across all exposures that had a signficant hit:
#Use maximal specification only:
Eid <- sort(unique(HEAPassoc@HEAPsig$Model6$E$sigBOTH$ID))

#Run Tissue Enrichments:
HEAPtissue <- pblapply(Eid, function(x) {
  Tgsea <- suppressMessages(suppressWarnings(
    TissueEnrichGSEA(Eid = x, CType = "Model6", Tsplit = "test")
  ))
  return(Tgsea)
})
names(HEAPtissue) <- Eid

#Run Pathway Enrichments:
HEAPpathway <- pblapply(Eid, function(x) {
  Pgsea <- suppressMessages(suppressWarnings(
    PathwayEnrichGSEA(Eid = x, CType = "Model6", Tsplit = "test")
  ))
  return(Pgsea)
})
names(HEAPpathway) <- Eid

#'*Store the Tissue and Pathway Enrichments in Object*
#'*Create DataStructure*
HEAPpathconstruct <- setClass(
  "HEAPpathconstruct",
  slots = c(
    HEAPtgsea = "list", #List: All association results
    HEAPpgsea = "list" #List: Replicated results (Significant in both Train and Test)
  )
)

#Create the datastructure to utilize:
HEAPgsea <- HEAPpathconstruct(
  HEAPtgsea = HEAPtissue,
  HEAPpgsea = HEAPpathway
)

#Save RDS file of PXSloader object
gc()
class(HEAPgsea)
qsave(HEAPgsea, file = "./Output/HEAPres/HEAPgsea.qs")
