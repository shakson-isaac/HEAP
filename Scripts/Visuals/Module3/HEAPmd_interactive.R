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



#'*GET interactive plot of GEM - for covarspec 5*

GEMweb <- HEAPmd@GEMlist$Type5 %>%
  mutate(Estimate = log(Estimate)) %>%
  na.omit()

p2 <- GEMweb %>%
  mutate(across(where(is.numeric), ~ signif(., digits = 3))) %>%
  ggplot(aes(x = NumDiseases, y = Estimate)) +
  geom_point(aes(text = paste("<br>Prot:", ID,
                              "<br>GEM:", paste0(Estimate)))) +
  geom_ribbon(aes(ymin =min(Estimate)*1.2, ymax = 0, fill = "Genetics Driven"), alpha = 0.2) +
  geom_ribbon(aes(ymin = 0, ymax = max(Estimate)*1.2, fill = "Environment Driven"), alpha = 0.2) +
  scale_fill_manual(values = c("Genetics Driven" = "blue", "Environment Driven" = "green")) +
  xlab("Number of Significant Disease Hits") +
  ylab("GEM") +
  theme_minimal() +
  guides(fill=guide_legend(title="Region")) +
  theme(legend.position = "bottom")
p2

gginteract <- ggplotly(p2, tooltip = "text")

saveWidget(gginteract, 
           file = paste0("./Output/App/Interactive/A1/GEMplot.html"))
