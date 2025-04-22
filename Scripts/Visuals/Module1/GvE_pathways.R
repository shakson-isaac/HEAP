#Load Libraries
library(data.table)
library(tidyverse)
library(ggplot2)
library(ggpmisc)
library(ggrepel)
library(gridExtra)
library(patchwork)

#Load Tables from Analysis:
setwd("/n/groups/patel/shakson_ukb/UK_Biobank/")
HEAPge <- readRDS("./Output/HEAPres/HEAPge.rds")
#source("./RScripts/Pure_StatGen/Prot_ExPGS/Visualizations/analyProt_PGS_PXS_final.R")

#'*GSEA Analysis of G, E, and GvE Ranked Lists*
#'Using datastructure: HEAPge@R2geCI (GvE mean across specifications)
library(enrichplot)
library(clusterProfiler)
library(ReactomePA)
library(DOSE)
library(org.Hs.eg.db)

# EntrezID and ProteinID Conversions
df_entrez <- fread(file = "/n/groups/patel/shakson_ukb/UK_Biobank/Data/GeneSets/OlinkEntrezConv.txt")

# Create Unified Dataframe with HEAP results:
GSEAdf <- merge(df_entrez, HEAPge@R2geCI, by.x = "Gene", by.y = "ID")


# Create Ranked Lists:
rListCreate <- function(df){
  rankedList <- list()
  
  #Rank: Use ties.method - random to handle ties
  gene_df <- df %>% 
    mutate(rank = rank(E_stats_mean,  ties.method = "random")) %>%
    arrange(desc(rank))
  Eprotlist <- gene_df$rank
  names(Eprotlist) <- gene_df$entrezgene_id
  
  gene_df <- df %>% 
    mutate(rank = rank(G_stats_mean,  ties.method = "random")) %>%
    arrange(desc(rank))
  Gprotlist <- gene_df$rank
  names(Gprotlist) <- gene_df$entrezgene_id
  
  gene_df <- df %>% 
    mutate(rank = rank(GxE_stats_mean,  ties.method = "random")) %>%
    arrange(desc(rank))
  GxEprotlist <- gene_df$rank
  names(GxEprotlist) <- gene_df$entrezgene_id
  
  rankedList[["E"]] <- Eprotlist
  rankedList[["G"]] <- Gprotlist
  rankedList[["GxE"]] <- GxEprotlist
  
  return(rankedList)
}

RankedList <- rListCreate(GSEAdf)


# Function to run enrichments:
GSEA <- function(rlist){
  GSEAlist <- list()
  go <- gseGO(geneList     = rlist,
              OrgDb        = org.Hs.eg.db,
              ont          = "BP",
              minGSSize    = 10,
              maxGSSize    = 500,
              pvalueCutoff = 0.05,
              verbose      = FALSE,
              scoreType = "pos")
  
  
  go2 <- gseGO(geneList     = rlist,
               OrgDb        = org.Hs.eg.db,
               ont          = "MF",
               minGSSize    = 10,
               maxGSSize    = 500,
               pvalueCutoff = 0.05,
               verbose      = FALSE,
               scoreType = "pos")
  
  Reactome <- gsePathway(rlist, 
                  minGSSize = 10, maxGSSize = 500,
                  pvalueCutoff = 0.05, pAdjustMethod = "BH", 
                  verbose = TRUE, eps = 0,
                  scoreType = "pos")
  
  GSEAlist[["GO_BP"]] <- go
  GSEAlist[["GO_MF"]] <- go2
  GSEAlist[["Reactome"]] <- Reactome
  
  return(GSEAlist)
}

Egsea <- GSEA(RankedList$E)
Ggsea <- GSEA(RankedList$G)
GxEgsea <- GSEA(RankedList$GxE)


#Get tables:
Egobp <- Egsea$GO_BP@result %>% dplyr::select(c("ID","Description","setSize","NES","p.adjust"))
Egobp$ProteinSet <- "Exposures"
Egobp$GOType <- "Biological Process"

Egomf <- Egsea$GO_MF@result %>% dplyr::select(c("ID","Description","setSize","NES","p.adjust"))
Egomf$ProteinSet <- "Exposures"
Egomf$GOType <- "Molecular Function"


Ggobp <- Ggsea$GO_BP@result %>% dplyr::select(c("ID","Description","setSize","NES","p.adjust"))
Ggobp$ProteinSet <- "Genetics"
Ggobp$GOType <- "Biological Process"


Ggomf <- Ggsea$GO_MF@result %>% dplyr::select(c("ID","Description","setSize","NES","p.adjust"))
Ggomf$ProteinSet <- "Genetics"
Ggomf$GOType <- "Molecular Function"

#In one table:
GvEpathways <- list(Egobp, Egomf, Ggobp,Ggomf) %>% purrr::reduce(full_join)



## PUT pathway analysis results (table) in ./Figures/HEAP/T1/
fwrite(Egobp, file = "./Output/App/Tables/Extra/Eprot_GOBP.csv")
fwrite(Egomf, file = "./Output/App/Tables/Extra/Eprot_GOMF.csv")
fwrite(Ggobp, file = "./Output/App/Tables/Extra/Gprot_GOBP.csv")
fwrite(Ggomf, file = "./Output/App/Tables/Extra/Gprot_GOMF.csv")
fwrite(GvEpathways, file = "./Output/App/Tables/GvE_pathways.csv")




#Get plots:
ego <- gseaplot2(Egsea$GO_BP, 
                 geneSetID = c("GO:0010883","GO:0034369",
                               "GO:0055088","GO:0033344"),
                 pvalue_table = TRUE)
ego2 <- gseaplot2(Egsea$GO_MF, 
                  geneSetID = c("GO:0038024","GO:0034185","GO:0043178"),
                  pvalue_table = TRUE)
ggo <- gseaplot2(Ggsea$GO_BP, 
                 geneSetID = c("GO:0032609","GO:0046634",
                               "GO:0002250"),
                 pvalue_table = TRUE)
ggo2 <- gseaplot2(Ggsea$GO_MF, 
                  geneSetID = c("GO:0140375","GO:0032393",
                                "GO:0033691"),
                  pvalue_table = TRUE)


## PUT the GSEA plots in ./Figures/HEAP/A1/
lapply(c("png", "svg", "pdf"), 
       function(fmt) ggsave(paste0("./Figures/HEAP/A1/","Eprot_GOBP.", fmt), 
                            plot = ego, width = 10, height = 6, dpi = 1000))

lapply(c("png", "svg", "pdf"), 
       function(fmt) ggsave(paste0("./Figures/HEAP/A1/","Eprot_GOMF.", fmt), 
                            plot = ego2, width = 10, height = 6, dpi = 1000))

lapply(c("png", "svg", "pdf"), 
       function(fmt) ggsave(paste0("./Figures/HEAP/A1/","Gprot_GOBP.", fmt), 
                            plot = ggo, width = 10, height = 6, dpi = 1000))

lapply(c("png", "svg", "pdf"), 
       function(fmt) ggsave(paste0("./Figures/HEAP/A1/","Gprot_GOMF.", fmt), 
                            plot = ggo2, width = 10, height = 6, dpi = 1000))


#c("GO:0034368","GO:0010883","GO:0034369","GO:0055088","GO:0033344")




