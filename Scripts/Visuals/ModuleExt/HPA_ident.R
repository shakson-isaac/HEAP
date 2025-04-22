library(BiocManager)
#BiocManager::install("HPAanalyze")
library(tidyverse)
library(data.table)


#TODO:
#1.) Load HPA and GTEX
#2.) Link with Omics IDs
#3.) Create relevant datastructures for significance tests


#### OBTAIN RELEVANT TABLES for TISSUE SPECIFICITY #####

#'*Use the files below to generate appropriate genesets*
tissueinfo <- fread(file = "/n/groups/patel/shakson_ukb/UK_Biobank/Data/HPA/normal_ihc_data.tsv")
secreteinfo <- fread(file = "/n/groups/patel/shakson_ukb/UK_Biobank/Data/HPA/subcellular_location.tsv")


#Figure out heuristic to define tissue specific proteins.
#Keywords: "Medium", "High", "Low", "N/A", "Ascending","Descending","Not representative"

#Do 3 classifications:
#Hyper-Specific - "High" only
#Enriched - "High" + "Medium" 
#Expressed - "High" + "Medium" + "Low" + "Ascending" + "Descending"


#Figure out heuristic for defining secreted proteins.
#Keyword: "Predicted to be secreted"
#unique(tissueinfo$Tissue)
#unique(tissueinfo$`IHC tissue name`)
#unique(tissueinfo$`Cell type`)

SpecificSet <- tissueinfo %>%
                    filter(Level %in% c("High")) %>%
                    select(c("Gene","Gene name","Tissue"))
EnrichedSet <- tissueinfo %>%
                    filter(Level %in% c("High","Medium")) %>%
                    select(c("Gene","Gene name","Tissue"))
ExpressedSet <- tissueinfo %>%
                    filter(Level %in% c("High","Medium","Low",
                                        "Ascending","Descending")) %>%
                    select(c("Gene","Gene name","Tissue"))
Secretome <- secreteinfo %>%
                  filter(`Extracellular location` %in% 
                           c("Predicted to be secreted")) %>%
                  mutate(Pathway = "Secretome") %>%
                  select(c("Gene","Gene name","Pathway"))


#'*ALSO Load KEGG*
library(KEGGREST)

gene2pathway = keggLink("pathway", "hsa")

df3 = data.frame(
  gene_id    = gsub("hsa:", "", names(gene2pathway)),
  pathway_id = gsub("path:", "", gene2pathway)
)

pathway2name = keggList("pathway", "hsa")
df4 = data.frame(
  gene_id    = names(pathway2name),
  pathway_id = gsub(" - Homo sapiens \\(human\\)", "", pathway2name)
)


#WRITE FILES:
write.table(SpecificSet[,c("Tissue","Gene name")], 
            file = "/n/groups/patel/shakson_ukb/UK_Biobank/Data/GeneSets/HPA_specific.txt",
            sep = "\t", row.names = F)
write.table(EnrichedSet[,c("Tissue","Gene name")], 
            file = "/n/groups/patel/shakson_ukb/UK_Biobank/Data/GeneSets/HPA_enriched.txt",
            sep = "\t", row.names = F)
write.table(ExpressedSet[,c("Tissue","Gene name")], 
            file = "/n/groups/patel/shakson_ukb/UK_Biobank/Data/GeneSets/HPA_expressed.txt",
            sep = "\t", row.names = F)
write.table(Secretome[,c("Pathway","Gene name")], 
            file = "/n/groups/patel/shakson_ukb/UK_Biobank/Data/GeneSets/HPA_secretome.txt",
            sep = "\t", row.names = F)
write.table(df3[,c("pathway_id","gene_id")], 
            file = "/n/groups/patel/shakson_ukb/UK_Biobank/Data/GeneSets/KEGG_T2G.txt",
            sep = "\t", row.names = F)
write.table(df4[,c("pathway_id","gene_id")], 
            file = "/n/groups/patel/shakson_ukb/UK_Biobank/Data/GeneSets/KEGG_T2N.txt",
            sep = "\t", row.names = F)

