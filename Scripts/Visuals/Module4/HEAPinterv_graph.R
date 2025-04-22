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
library(readxl)
library(qs)

#Load INT structure:
setwd("/n/groups/patel/shakson_ukb/UK_Biobank/")

HEAPassoc <- qread("./Output/HEAPres/HEAPassoc.qs")
HEAPint <- qread("./Output/HEAPres/HEAPint.qs")
HEAPmd <- qread("./Output/HEAPres/HEAPmd.qs")

HEAPdf <- HEAPassoc@HEAPlist$Model6$test[[1]]
INTdf <- HEAPint@sList$Model6
MDdf <- HEAPmd@MDlist$Type5
GEMdf <- HEAPmd@GEMlist$Type5

#'*Exposure-Disease-Intervention Network:*
#'*Label HEAP associations in Network*
UKBnet <- HEAPdf %>% 
  filter(`Pr(>|t|)` < 0.05/n()) %>%
  select(ID, omicID, Estimate)

eID_net <- UKBnet %>%
  count(ID) %>%
  filter(n > 10) %>%
  pull(ID)

UKBnet <- UKBnet %>%
  filter(ID %in% eID_net)

eID_select <- c("alcohol_intake_frequency_f1558_0_06",
                "beef_intake_f1369_0_04",
                "fresh_fruit_intake_f1309_0_0",
                "bread_type_f1448_0_0_White",
                "current_tobacco_smoking_f1239_0_0_Yes._on_most_or_all_days",
                "number_of_days_week_of_vigorous_physical_activity_10_plus_minutes_f904_0_0",
                "types_of_physical_activity_in_last_4_weeks_f6164_0_0.multi_Strenuous_sports"
)

UKBnet_select <- UKBnet %>%
  filter(ID %in% eID_select) %>%
  mutate(Estimate = if_else(Estimate > 0, 1, -1))

colnames(UKBnet_select) <- c("Intervention","Protein","Effect")


#'*Significant effects in HERITAGE and GLP1*
colnames(INTdf)[which(names(INTdf) == "EntrezGeneSymbol")] <- "Protein"

Interventions <- INTdf %>% 
                    select(HERITAGE_effect, GLP1_effect1, GLP1_effect2, Protein) %>%
                    unique() %>%
                    pivot_longer(cols = !Protein,
                                 names_to = "Intervention", 
                                 values_to = "Effect") %>%
                    mutate(Effect = if_else(Effect > 0, 1, -1)) %>%
                    filter(Protein %in% unique(UKBnet_select$Protein)) %>%
                    na.omit()


UKB_IntNet <- list(UKBnet_select, Interventions) %>% reduce(full_join)




#'*Subset the above with mediation analysis results:*
T2D_assoc <- MDdf %>%
  filter(DZid == "age_e11_first_reported_non_insulin_dependent_diabetes_mellitus_f130708_0_0")

#Filter by nominal significance:
T2D_assoc1 <- T2D_assoc %>%
  select(ID, HR, HR_lower95, HR_upper95,
         cindex, E_HRi, G_HRi) %>%
  filter((HR_lower95 > 1 & HR_upper95 > 1) | 
           (HR_lower95 < 1 & HR_upper95 < 1))

#Get the GEM statistic Attribute
GEMstats <- GEMdf %>%
  mutate(GEM = log(Estimate)) %>%
  select(ID, GEM) %>%
  na.omit()

T2D_ProtAttr <- list(T2D_assoc1, GEMstats) %>% reduce(left_join)

T2D_ProtAttr <- T2D_ProtAttr %>%
  filter(HR > 1.2 | HR < 0.8)


#'*Subset UKB exposure effects for T2D specific effects:*
UKB_IntNet <- UKB_IntNet %>%
  filter(Protein %in% unique(T2D_ProtAttr$ID))

colnames(T2D_ProtAttr)[which(names(T2D_ProtAttr) == "ID")] <- "Protein"


UKB_IntNet_final <- list(UKB_IntNet, T2D_ProtAttr) %>% reduce(inner_join)
UKB_IntNet_finalv2 <- UKB_IntNet_final %>%
  filter(HR > 1.75 | HR < 0.75)

ProteinSet <- UKB_IntNet_finalv2 %>%
  filter(Intervention %in% c("HERITAGE_effect","GLP1_effect1","GLP1_effect2")) %>%
  pull(Protein) %>%
  unique()

UKB_IntNet_finalv2  <- UKB_IntNet_finalv2 %>%
  filter(Protein %in% ProteinSet)

#'*Add the Protein-T2D (links) rows:*
T2Dlinkage <- UKB_IntNet_finalv2 %>%
  select(Protein, HR) %>%
  unique() %>%
  mutate(Effect = if_else(HR > 1, 1, -1)) %>%
  select(Protein, Effect) %>%
  rename(Intervention = Protein)
T2Dlinkage$Protein <- "T2D"

UKB_IntNet_finalv3 <- list(T2Dlinkage, UKB_IntNet_finalv2) %>% reduce(full_join)
UKB_IntNet_finalv3[is.na(UKB_IntNet_finalv3)] <- 0

length(unique(UKB_IntNet_finalv2$Protein))

write.csv(UKBnet_select, file = "/n/groups/patel/shakson_ukb/UK_Biobank/Output/Network/UKBnet.csv",
          row.names =  F)
write.csv(T2D_ProtAttr, file = "/n/groups/patel/shakson_ukb/UK_Biobank/Output/Network/UKBnet_attributes.csv",
          row.names =  F)
write.csv(UKB_IntNet_final, file = "/n/groups/patel/shakson_ukb/UK_Biobank/Output/Network/UKBnet_final.csv",
          row.names =  F)
write.csv(UKB_IntNet_finalv3, file = "/n/groups/patel/shakson_ukb/UK_Biobank/Output/Network/UKBnet_finalv2.csv",
          row.names =  F)

### Network Construction Logic
# NODE1 - Intervention
# NODE2 - Protein
# NODE3 - Type 2 Diabetes

## Edges: 
# NODE1 - NODE2 (ONLY SIGNIFICANT Associations considered)
#     Sign direction (+1, -1) is it positive or negative correlation
#     UKB (Exercise, Diet, Smoking, etc.); HERITAGE; GLP RCTs
# Example:
# Having specific exposure links would be helpful.


# NODE2 - NODE3 (From mediation analysis)
#     Use the HR of the protein for positive vs negative links
#     Find if nominally significant (1, 0). From the HR confidence interval


## Attributes: (GET from mediation analysis)
# NODE1 - None
# NODE2 - 
#   GEM statistic - Use cytoscape to apply gradient coloring. Center on 0 if logFC
#   HR - protein for T2D
#   c-index - protein for T2D
# NODE3 - None

#Initial Steps:
#Get NODE1 - NODE2 in this script