#Libraries:
library(data.table)
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(ggpmisc)
library(pbapply)
library(dplyr)
library(purrr)
library(qs) #For faster save/load of rds object.

#Load MD structure:
setwd("/n/groups/patel/shakson_ukb/UK_Biobank/")
HEAPmd <- qread("./Output/HEAPres/HEAPmd.qs")

#'*Save tables for the website*
#'The GEM score can just be interactive and then create a file to download:

#Download the GEM score: MODEL 5
GEMweb <- HEAPmd@GEMlist$Type5 %>%
  mutate(Estimate = log(Estimate)) %>%
  na.omit()

GEMdownload <- GEMweb %>%
  select(c("ID","Estimate"))
colnames(GEMdownload) <- c("Protein","GEM")
fwrite(GEMdownload, 
       file = "/n/groups/patel/shakson_ukb/UK_Biobank/Output/App/Tables/GEMdownload.csv")


### Download Mediation Results: Model 5
hist(HEAPmd@MDlist$Type5$R2)
summary(HEAPmd@MDlist$Type5$R2)
quantile(HEAPmd@MDlist$Type5$R2, c(.01, .05, .95, 0.99)) 

#Only show proteins with valid GEM!
MDweb <- HEAPmd@MDlist$Type5 %>%
  filter(R2 > 0.045) %>%
  filter(if_all(where(is.numeric), is.finite)) %>%
  filter(ID %in% GEMweb$ID)
length(unique(MDweb$ID))

MDweb_fin <- MDweb %>%
  mutate(Indirect_G = exp(`Genetic Indirect Effect`),
         Indirect_E = exp(`Exposure Indirect Effect`),
         Direct_G = exp(`Genetic Direct Effect`),
         Direct_E = exp(`Exposure Direct Effect`) ,
         Total_Direct = exp(`Total Direct Effect`),
         Total_Indirect = exp(`Total Indirect Effect`)) %>%
  select(c("ID","DZid","Indirect_G","Indirect_E",
           "Direct_G","Direct_E","Total_Direct",
           "Total_Indirect","Eprop_mediated",
           "Gprop_mediated"))

colnames(MDweb_fin) <- c("Protein","Disease Code","G: Indirect",
                         "E: Indirect", "G: Direct",
                         "E: Direct","Total Direct",
                         "Total Indirect","E: Prop Mediate",
                         "G: Prop Mediate")

MDweb_fin <- MDweb_fin %>%
  mutate(across(where(is.numeric), ~ signif(., digits = 3)))

fwrite(MDweb_fin, 
       file = "/n/groups/patel/shakson_ukb/UK_Biobank/Output/App/Tables/MediationResults.csv")


### Extend to All Models:

MDtables <- function(Model){
  MDlist <- list()
  
  GEMweb <- HEAPmd@GEMlist[[Model]] %>%
    mutate(Estimate = log(Estimate)) %>%
    na.omit()
  
  GEMdownload <- GEMweb %>%
    select(c("ID","Estimate"))
  colnames(GEMdownload) <- c("Protein","GEM")
  
  MDweb <- HEAPmd@MDlist[[Model]] %>%
    filter(R2 > 0.045) %>%
    filter(if_all(where(is.numeric), is.finite)) %>%
    filter(ID %in% GEMweb$ID)
  length(unique(MDweb$ID))
  
  MDweb_fin <- MDweb %>%
    mutate(Indirect_G = exp(`Genetic Indirect Effect`),
           Indirect_E = exp(`Exposure Indirect Effect`),
           Direct_G = exp(`Genetic Direct Effect`),
           Direct_E = exp(`Exposure Direct Effect`) ,
           Total_Direct = exp(`Total Direct Effect`),
           Total_Indirect = exp(`Total Indirect Effect`)) %>%
    select(c("ID","DZid","Indirect_G","Indirect_E",
             "Direct_G","Direct_E","Total_Direct",
             "Total_Indirect","Eprop_mediated",
             "Gprop_mediated"))
  
  colnames(MDweb_fin) <- c("Protein","Disease Code","G: Indirect",
                           "E: Indirect", "G: Direct",
                           "E: Direct","Total Direct",
                           "Total Indirect","E: Prop Mediate",
                           "G: Prop Mediate")
  
  MDweb_fin <- MDweb_fin %>%
    mutate(across(where(is.numeric), ~ signif(., digits = 3)))
  
  MDlist[["GEMload"]] <-  GEMdownload
  MDlist[["MDweb"]] <- MDweb_fin
  return(MDlist)
}
models <- c("Type1","Type2","Type3","Type4","Type5")
MDfin <- lapply(models, MDtables)
names(MDfin) <- models

lapply(names(MDfin), function(x){
  gem <- MDfin[[x]]$GEMload
  md <- MDfin[[x]]$MDweb
  
  fwrite(gem, file = paste0("/n/groups/patel/shakson_ukb/UK_Biobank/Output/App/Tables/Models/",x,"GEMdownload.csv"))
  fwrite(md, file = paste0("/n/groups/patel/shakson_ukb/UK_Biobank/Output/App/Tables/Models/",x,"MediationResults.csv"))
})
