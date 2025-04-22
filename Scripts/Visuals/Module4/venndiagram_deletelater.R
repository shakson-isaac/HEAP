##### VennDiagram Stuff: #####

### Preprocess:
library(ggVennDiagram)
library(ggplot2)
library(dplyr)

### Venn Diagram (Determine Overlaps)
unique(CovarSpecList[["Type6"]]$test[[1]]$ID)

# Get significant hits:
eID = "number_of_days_week_of_vigorous_physical_activity_10_plus_minutes_f904_0_0"
eID = "fresh_fruit_intake_f1309_0_0"
eID = "smoking_status_f20116_0_0_Current"
eID = "alcohol_intake_frequency_f1558_0_06"
UKB <- CovarSpecList[["Type6"]]$test[[1]] %>% 
  filter((ID == eID) & (`Pr(>|t|)` < 0.05/n()))
HERITAGE <- data %>% filter(`False Discovery Rate (q-value)` < 0.05)
GLP1 <- STEP1 %>% filter(qvalue < 0.05) %>% 
  group_by(EntrezGeneSymbol) %>%
  mutate(GLP1_effect = mean(effect_size),
         GLP1_se = max(std_error),
         GLP1_qvalue = max(qvalue)) %>%
  select(c(EntrezGeneSymbol, GLP1_effect, GLP1_se, GLP1_qvalue)) %>%
  unique()


# Overlapping Significant Proteins in each study
venn_data <- list(
  UKB = UKB$omicID,
  HERITAGE = HERITAGE$EntrezGeneSymbol,
  Study3 = GLP1$EntrezGeneSymbol)

# Create the 3-set Venn diagram
ggVennDiagram(venn_data) + scale_fill_gradient(low="grey90",high = "red")
#ggtitle("Concordant vs Non-Concordant Associations Between Studies")



### Need to merge the datasets:

interventions <- merge(HERITAGE, GLP1, by = "EntrezGeneSymbol")
ukb_validate <- merge(UKB, interventions, by.x = "omicID", by.y = "EntrezGeneSymbol")

sum(sign(ukb_validate$Estimate) == sign(ukb_validate$`log(10) Fold Change`))
sum(sign(ukb_validate$Estimate) == sign(ukb_validate$GLP1_effect))
sum(sign(ukb_validate$Estimate) == sign(ukb_validate$`log(10) Fold Change`) & 
      sign(ukb_validate$`log(10) Fold Change`) == sign(ukb_validate$GLP1_effect))

ukb_validate %>%
  filter(sign(ukb_validate$Estimate) == sign(ukb_validate$`log(10) Fold Change`) & 
           sign(ukb_validate$`log(10) Fold Change`) == sign(ukb_validate$GLP1_effect)) %>%
  pull(omicID)


UKB_HERITAGE <- merge(UKB, HERITAGE, by.x = "omicID", by.y = "EntrezGeneSymbol")
HERITAGE_GLP1 <- merge(HERITAGE, GLP1, by = "EntrezGeneSymbol")
UKB_GLP1 <- merge(UKB, GLP1, by.x = "omicID", by.y = "EntrezGeneSymbol")

a_b <- nrow(UKB_HERITAGE) 
ab <- sum(sign(UKB_HERITAGE$Estimate) == sign(UKB_HERITAGE$`log(10) Fold Change`))
b_c <- nrow(HERITAGE_GLP1)
bc <- sum(sign(HERITAGE_GLP1$`log(10) Fold Change`) == sign(HERITAGE_GLP1$GLP1_effect))
a_c <- nrow(UKB_GLP1)
ac <- sum(sign(UKB_GLP1$Estimate) == sign(UKB_GLP1$GLP1_effect))

#Understand Relationships: How many proteins concordant direction between studies:
UKB_con <- unique(c(UKB_HERITAGE %>% 
                      filter(sign(UKB_HERITAGE$Estimate) == sign(UKB_HERITAGE$`log(10) Fold Change`)) %>% 
                      pull(omicID),
                    UKB_GLP1 %>%
                      filter(sign(UKB_GLP1$Estimate) == sign(UKB_GLP1$GLP1_effect)) %>%
                      pull(omicID)))
HERITAGE_con <- unique(c(UKB_HERITAGE %>% 
                           filter(sign(UKB_HERITAGE$Estimate) == sign(UKB_HERITAGE$`log(10) Fold Change`)) %>% 
                           pull(omicID), 
                         HERITAGE_GLP1 %>%
                           filter(sign(HERITAGE_GLP1$`log(10) Fold Change`) == sign(HERITAGE_GLP1$GLP1_effect)) %>%
                           pull(EntrezGeneSymbol)))
GLP1_con <- unique(c(UKB_GLP1 %>%
                       filter(sign(UKB_GLP1$Estimate) == sign(UKB_GLP1$GLP1_effect)) %>%
                       pull(omicID), 
                     HERITAGE_GLP1 %>%
                       filter(sign(HERITAGE_GLP1$`log(10) Fold Change`) == sign(HERITAGE_GLP1$GLP1_effect)) %>%
                       pull(EntrezGeneSymbol)))
UKB_noncon <- unique(c(UKB_HERITAGE %>% 
                         filter(sign(UKB_HERITAGE$Estimate) != sign(UKB_HERITAGE$`log(10) Fold Change`)) %>% 
                         pull(omicID),
                       UKB_GLP1 %>%
                         filter(sign(UKB_GLP1$Estimate) != sign(UKB_GLP1$GLP1_effect)) %>%
                         pull(omicID)))
HERITAGE_noncon <- unique(c(UKB_HERITAGE %>% 
                              filter(sign(UKB_HERITAGE$Estimate) != sign(UKB_HERITAGE$`log(10) Fold Change`)) %>% 
                              pull(omicID), 
                            HERITAGE_GLP1 %>%
                              filter(sign(HERITAGE_GLP1$`log(10) Fold Change`) != sign(HERITAGE_GLP1$GLP1_effect)) %>%
                              pull(EntrezGeneSymbol)))
GLP1_noncon <- unique(c(UKB_GLP1 %>%
                          filter(sign(UKB_GLP1$Estimate) != sign(UKB_GLP1$GLP1_effect)) %>%
                          pull(omicID), 
                        HERITAGE_GLP1 %>%
                          filter(sign(HERITAGE_GLP1$`log(10) Fold Change`) != sign(HERITAGE_GLP1$GLP1_effect)) %>%
                          pull(EntrezGeneSymbol)))


venn_data <- list(
  UKB = c(UKB_noncon),
  HERITAGE = c(HERITAGE_noncon),
  GLP1 = c(GLP1_noncon))
ggVennDiagram(venn_data) + scale_fill_gradient(low="grey90",high = "red")

venn_data <- list(
  UKB = c(UKB_con),
  HERITAGE = c(HERITAGE_con),
  GLP1 = c(GLP1_con))
ggVennDiagram(venn_data) + scale_fill_gradient(low="grey90",high = "red")














# Example data for each study (protein identifiers, directionality, and significance)
study1 <- data.frame(protein = c("Protein1", "Protein2", "Protein3", "Protein4", "Protein5"),
                     direction_study1 = c("positive", "negative", "positive", "negative", "positive"),
                     significant = c(TRUE, TRUE, TRUE, FALSE, TRUE))

study2 <- data.frame(protein = c("Protein1", "Protein2", "Protein5", "Protein6", "Protein3"),
                     direction_study2 = c("positive", "positive", "negative", "positive", "positive"),
                     significant = c(TRUE, TRUE, TRUE, FALSE, TRUE))

study3 <- data.frame(protein = c("Protein1", "Protein2", "Protein3", "Protein6", "Protein5"),
                     direction_study3 = c("positive", "positive", "positive", "negative", "positive"),
                     significant = c(TRUE, TRUE, TRUE, FALSE, TRUE))

# Subset the data to include only significant proteins
study1_significant <- subset(study1, significant == TRUE)
study2_significant <- subset(study2, significant == TRUE)
study3_significant <- subset(study3, significant == TRUE)

# Merge the data into a combined dataframe
combined <- merge(study1_significant, study2_significant, by = "protein", suffixes = c("_study1", "_study2"))
combined <- merge(combined, study3_significant, by = "protein")

# Check for concordance (same direction across studies)
combined$concordant_12 <- with(combined, 
                               ifelse(direction_study1 == direction_study2, "Concordant", "Non-concordant"))

combined$concordant_13 <- with(combined, 
                               ifelse(direction_study1 == direction_study3, "Concordant", "Non-concordant"))

combined$concordant_23 <- with(combined, 
                               ifelse(direction_study2 == direction_study3, "Concordant", "Non-concordant"))

combined$concordant_123 <- with(combined, 
                                ifelse(direction_study1 == direction_study2 & direction_study2 == direction_study3, "Concordant", "Non-concordant"))

# Create subsets for each combination of concordant and non-concordant proteins
study1_2 <- unique(combined[combined$concordant_12 == "Concordant", "protein"])
study1_3 <- unique(combined[combined$concordant_13 == "Concordant", "protein"])
study2_3 <- unique(combined[combined$concordant_23 == "Concordant", "protein"])
study1_2_3 <- unique(combined[combined$concordant_123 == "Concordant", "protein"])

# Non-concordant sets (opposite directions)
study1_2_non <- unique(combined[combined$concordant_12 == "Non-concordant", "protein"])
study1_3_non <- unique(combined[combined$concordant_13 == "Non-concordant", "protein"])
study2_3_non <- unique(combined[combined$concordant_23 == "Non-concordant", "protein"])

# Combine everything into a list for the Venn diagram
venn_data <- list(
  `Study1 & Study2 (Concordant)` = study1_2,
  `Study1 & Study3 (Concordant)` = study1_3,
  `Study2 & Study3 (Concordant)` = study2_3,
  `Study1, Study2 & Study3 (Concordant)` = study1_2_3,
  `Study1 & Study2 (Non-Concordant)` = study1_2_non,
  `Study1 & Study3 (Non-Concordant)` = study1_3_non,
  `Study2 & Study3 (Non-Concordant)` = study2_3_non
)

# Create Venn diagram for the concordant and non-concordant associations
ggVennDiagram(venn_data)







