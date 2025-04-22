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

#'*LOAD CovarSpecList*
setwd("/n/groups/patel/shakson_ukb/UK_Biobank/")
HEAPassoc <- qread("./Output/HEAPres/HEAPassoc.qs")
CovarSpecAssocList <- HEAPassoc@HEAPsig


#'*FINISH BY:*
#'LOAD in
#'*CovarSpecAssocList*

##### Get Tables for Univariate Associations: E and GxE #####
# E Associations:
fwrite(CovarSpecAssocList$Model1$E$sigBOTH,
       file = "/n/groups/patel/shakson_ukb/UK_Biobank/Output/App/Tables/Models/Model1EassocReplicated.csv")
fwrite(CovarSpecAssocList$Model2$E$sigBOTH,
       file = "/n/groups/patel/shakson_ukb/UK_Biobank/Output/App/Tables/Models/Model2EassocReplicated.csv")
fwrite(CovarSpecAssocList$Model3$E$sigBOTH,
       file = "/n/groups/patel/shakson_ukb/UK_Biobank/Output/App/Tables/Models/Model3EassocReplicated.csv")
fwrite(CovarSpecAssocList$Model4$E$sigBOTH,
       file = "/n/groups/patel/shakson_ukb/UK_Biobank/Output/App/Tables/Models/Model4EassocReplicated.csv")
fwrite(CovarSpecAssocList$Model5$E$sigBOTH,
       file = "/n/groups/patel/shakson_ukb/UK_Biobank/Output/App/Tables/Models/Model5EassocReplicated.csv")
fwrite(CovarSpecAssocList$Model6$E$sigBOTH,
       file = "/n/groups/patel/shakson_ukb/UK_Biobank/Output/App/Tables/Models/Model6EassocReplicated.csv")

# GxE Associations:
fwrite(CovarSpecAssocList$Model1$GxE$sigBOTH,
       file = "/n/groups/patel/shakson_ukb/UK_Biobank/Output/App/Tables/Models/Model1GxEassocReplicated.csv")
fwrite(CovarSpecAssocList$Model2$GxE$sigBOTH,
       file = "/n/groups/patel/shakson_ukb/UK_Biobank/Output/App/Tables/Models/Model2GxEassocReplicated.csv")
fwrite(CovarSpecAssocList$Model3$GxE$sigBOTH,
       file = "/n/groups/patel/shakson_ukb/UK_Biobank/Output/App/Tables/Models/Model3GxEassocReplicated.csv")
fwrite(CovarSpecAssocList$Model4$GxE$sigBOTH,
       file = "/n/groups/patel/shakson_ukb/UK_Biobank/Output/App/Tables/Models/Model4GxEassocReplicated.csv")
fwrite(CovarSpecAssocList$Model5$GxE$sigBOTH,
       file = "/n/groups/patel/shakson_ukb/UK_Biobank/Output/App/Tables/Models/Model5GxEassocReplicated.csv")
fwrite(CovarSpecAssocList$Model6$GxE$sigBOTH,
       file = "/n/groups/patel/shakson_ukb/UK_Biobank/Output/App/Tables/Models/Model6GxEassocReplicated.csv")


fwrite(CovarSpecAssocList$Model6$E$sigBOTH,
       file = "/n/groups/patel/shakson_ukb/UK_Biobank/Output/App/Tables/ReplicatedEassoc.csv")
fwrite(CovarSpecAssocList$Model6$GxE$sigBOTH,
       file = "/n/groups/patel/shakson_ukb/UK_Biobank/Output/App/Tables/ReplicatedGxEassoc.csv")




