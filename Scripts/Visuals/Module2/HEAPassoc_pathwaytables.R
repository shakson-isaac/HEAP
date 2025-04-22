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

#For heatmaps
library(ComplexHeatmap)
library(circlize)

#'*LOAD GSEA results*
setwd("/n/groups/patel/shakson_ukb/UK_Biobank/")
HEAPht <- qread("./Output/HEAPres/HEAPgsea.qs")

#Tissue Effect Dataframe:
# purrr::map2 iterate over the list created to generate one dataframe:
TissueEffects <- map2_dfr(HEAPht@HEAPtgsea, names(HEAPht@HEAPtgsea), ~ .x %>%
                            select(all_of(c("ID","setSize","NES","p.adjust"))) %>%
                            mutate(cID = .y))
rownames(TissueEffects) <- NULL

#Pathway Effects Dataframe:
PathEffects <- map2_dfr(HEAPht@HEAPpgsea, names(HEAPht@HEAPpgsea), ~ .x %>%
                          select(all_of(c("Description","setSize","NES","p.adjust"))) %>%
                          mutate(cID = .y))
rownames(PathEffects) <- NULL

fwrite(TissueEffects, file = "./Output/App/Tables/EassocTissueEnrichment.csv")
fwrite(PathEffects, file = "./Output/App/Tables/EassocPathwayEnrichment.csv")

