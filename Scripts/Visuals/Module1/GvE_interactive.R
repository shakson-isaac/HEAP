library(data.table)
library(tidyverse)
library(ggplot2)
library(ggpmisc)
library(ggrepel)
library(gridExtra)
library(plotly)
library(htmlwidgets)

#Addn Libraries for Visualization:
library(patchwork)

### SPECIFIC WARNINGS!!!
#'* We DO NOT save all html elements here to prevent massive files.*
#'*If selfcontained = TRUE the html files are 3.8 MB each!!!*

# Output tables for specific components:
setwd("/n/groups/patel/shakson_ukb/UK_Biobank/")

# Load Summary Stats:
HEAPge <- readRDS("./Output/HEAPres/HEAPge.rds")


# Interactive plot of GvE: 
library(plotly)
library(htmlwidgets)
CovarSpec_GvEplot_interact <- function(R2df){
  #Highlight specific proteins with high E: R2
  subset_test <- R2df %>% 
    filter((E_stats_mean > G_stats_mean) & (E_stats_mean > 0.09))
  
  #Cutoff for y-axis:
  y_cutoff = 0.175
  
  #Without Error Bars:
  p1 <- R2df  %>% 
    mutate(across(where(is.numeric), ~ signif(., digits = 2))) %>%
    ggplot(aes(x = G_stats_mean, y = E_stats_mean, label = ID)) + 
    geom_point(aes(text = paste("<br>Prot:", ID,
                                "<br>G:", paste0(G_stats_mean, ": (",
                                                 G_stats_ci_low,",",G_stats_ci_high,")"),
                                "<br>E:", paste0(E_stats_mean,": (",
                                                 E_stats_ci_low,",",E_stats_ci_high,")")
    )
    )
    ) +
    geom_abline(slope=1, intercept = 0) +
    geom_ribbon(aes(ymin = 0, ymax = ifelse(G_stats_mean <= y_cutoff, G_stats_mean, y_cutoff), fill = "Genetics Driven"), alpha = 0.2) +
    geom_ribbon(aes(ymin = ifelse(G_stats_mean <= y_cutoff, G_stats_mean, y_cutoff), ymax = y_cutoff, fill = "Environment Driven"), alpha = 0.2) +
    scale_fill_manual(values = c("Genetics Driven" = "blue", "Environment Driven" = "green")) +
    ylim(c(0, y_cutoff)) +
    #geom_text_repel(data=subset_test, aes(x = G_stats_mean,y = E_stats_mean,label= ID),
    #                size = 4) +
    xlab("G: PGS R<sup>2</sup>") +
    ylab("E: PXS R<sup>2</sup>") +
    theme_minimal() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    guides(fill=guide_legend(title="Region")) +
    ggtitle("Test Set: Partitioned R<sup>2</sup>") +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "bottom") 
  
  gginteract <- ggplotly(p1, tooltip = "text")
  
  saveWidget(gginteract, 
             file = paste0("./Output/App/Interactive/A1/GvEplot.html"))
  
}
CovarSpec_GvEplot_interact(HEAPge@R2geCI)


