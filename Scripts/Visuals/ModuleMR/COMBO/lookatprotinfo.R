library(readr)
library(tidyverse)
Type5 <- read_csv("Output/Mediation/Results/Type5.csv")
colnames(Type5)

APOF <- Type5 %>% filter(prot_term == "APOF")
PRSS8 <- Type5 %>% filter(prot_term == "PRSS8")
FABP4 <- Type5 %>% filter(prot_term == "FABP4")
FURIN <- Type5 %>% filter(prot_term == "FURIN")
NCAN <- Type5 %>% filter(prot_term == "NCAN")

