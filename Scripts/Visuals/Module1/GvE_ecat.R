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

#### CATEGORY PLOTS ####
library(cowplot)
library(ggpubr)
category_specificR2plot <- function(plot_testR2_all, category, category_name, G, E, y_cutoff){
  subset_test <- plot_testR2_all  %>% 
    slice_max(.data[[category]], n = 10) 
  
  
  p1 <- plot_testR2_all %>% 
    ggplot(aes(x = .data[[G]], y = .data[[category]], label = ID)) + 
    geom_point(aes(color = .data[[E]])) + #aes(color = .data[[category]])) +
    geom_abline(slope=1, intercept = 0) +
    geom_ribbon(aes(ymin = 0, ymax = ifelse(.data[[G]] <= y_cutoff, .data[[G]], y_cutoff), fill = "Genetics Driven"), alpha = 0.2) +
    geom_ribbon(aes(ymin = ifelse(.data[[G]] <= y_cutoff, .data[[G]], y_cutoff), ymax = y_cutoff, fill = "Environment Driven"), alpha = 0.2) +
    scale_fill_manual(values = c("Genetics Driven" = "blue", "Environment Driven" = "green")) +
    ylim(c(0, y_cutoff)) +
    geom_text_repel(data=subset_test, aes(x = .data[[G]], y = .data[[category]], label= ID),
                    size = 4) +
    #geom_label() + 
    scale_color_gradient(low = "darkslategrey", high = "red") +
    xlab(expression("G: PGS"~R^2)) +
    ylab(bquote(~.(category_name)~":"~R^2)) +
    labs(color="E: PXS"~R^2) +
    theme_minimal() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) + #Adjust 90 to 0 for this plot:
    
    #Remove: the Environment vs Genetics Driven Legend for spacing:
    guides(fill = "none",  # Remove the fill legend
           color ="none"#guide_colorbar(title = bquote("PXS:" ~ R^2))  # Keep the color gradient legend
    ) +
    #ggtitle(bquote("Test Set:"~.(category_name))) +
    ggtitle(category_name) +
    theme(plot.title = element_text(hjust = 0.5, size = 14),
          legend.position = "bottom",
          legend.title = element_text(angle = 0),  # Optional: Set title orientation
          legend.text = element_text(angle = 90),  # Rotate legend text to vertical
          axis.text = element_text(size = 10),
          axis.title = element_text(size = 12)) #"right")
  p1
  return(p1)
}
get_legend <- function(plot_testR2_all, category, category_name, G, E, y_cutoff){
  p2 <- plot_testR2_all %>% 
    ggplot(aes(x = .data[[G]], y = .data[[category]], label = ID)) + 
    geom_point(aes(color = .data[[E]])) + #aes(color = .data[[category]])) +
    scale_color_gradient(low = "darkslategrey", high = "red", breaks=c(0,0.025,0.05,0.075,0.1,0.125,0.15)) +
    xlab(expression("G: PGS"~R^2)) +
    ylab(bquote(~.(category_name)~":"~R^2)) +
    labs(color="E: PXS"~R^2) +
    theme_minimal()
  
  library(ggpubr)
  leg <- get_legend(p2)
  p3 <- as_ggplot(leg)
  return(p3)
}

g1 <- category_specificR2plot(HEAPge@R2allCI, category = "mean_Exercise_Freq",
                              category_name = "Exercise_Freq",
                              G = "G_stats_mean",E = "E_stats_mean", y_cutoff = 0.1)

g2 <- category_specificR2plot(HEAPge@R2allCI, category = "mean_Exercise_MET",
                              category_name = "Exercise_MET",
                              G = "G_stats_mean",E = "E_stats_mean", y_cutoff = 0.06)

g3 <- category_specificR2plot(HEAPge@R2allCI, category = "mean_Diet_Weekly",
                              category_name = "Diet_Weekly",
                              G = "G_stats_mean",E = "E_stats_mean", y_cutoff = 0.06)

g4 <- category_specificR2plot(HEAPge@R2allCI, category = "mean_Alcohol",
                              category_name = "Alcohol",
                              G = "G_stats_mean",E = "E_stats_mean", y_cutoff = 0.06)

g5 <- category_specificR2plot(HEAPge@R2allCI, category = "mean_Smoking",
                              category_name = "Smoking",
                              G = "G_stats_mean",E = "E_stats_mean", y_cutoff = 0.045)

g6 <- category_specificR2plot(HEAPge@R2allCI, category = "mean_Vitamins",
                              category_name = "Vitamins",
                              G = "G_stats_mean",E = "E_stats_mean", y_cutoff = 0.025)

g7 <- category_specificR2plot(HEAPge@R2allCI, category = "mean_Internet_Usage",
                              category_name = "Internet_Usage",
                              G = "G_stats_mean",E = "E_stats_mean", y_cutoff = 0.02)

g8 <- category_specificR2plot(HEAPge@R2allCI, category = "mean_Deprivation_Indices",
                              category_name = "Deprivation_Indices",
                              G = "G_stats_mean",E = "E_stats_mean", y_cutoff = 0.02)


library(patchwork)
combined_plot <- (g1 | g2 | g3 | g4) /
  (g5 | g6 | g7 | g8)
#combined_plot
ggsave("./Figures/HEAP/A1/meanR2_Categories.png", 
       plot = combined_plot, width = 15, height = 8, units = "in")
ggsave("./Figures/HEAP/A1/meanR2_Categories.svg", 
       plot = combined_plot, width = 15, height = 8, units = "in")
ggsave("./Figures/HEAP/A1/meanR2_Categories.pdf", 
       plot = combined_plot, width = 15, height = 8, units = "in")


#'*HEATMAP of R2 by Category*
library(ComplexHeatmap)

catNames <- c(Alcohol = "mean_Alcohol",
              Diet_Weekly = "mean_Diet_Weekly",
              Smoking = "mean_Smoking",
              Exercise_MET = "mean_Exercise_MET",
              Exercise_Freq = "mean_Exercise_Freq",
              Internet_Usage = "mean_Internet_Usage",
              Deprivation_Indices = "mean_Deprivation_Indices",
              Vitamins = "mean_Vitamins")

testR2_heatmap <- HEAPge@R2catCI %>%
                    select(c("ID","mean_Alcohol","mean_Diet_Weekly",
                             "mean_Smoking","mean_Exercise_MET",
                             "mean_Exercise_Freq","mean_Internet_Usage",
                             "mean_Deprivation_Indices","mean_Vitamins")) %>%
                    rename(all_of(catNames))

#Top 70 proteins by R2:
topEprot <- HEAPge@R2geCI %>% 
              arrange(desc(E_stats_mean)) %>%
              head(70) %>%
              pull(ID)

#Select top 10 proteins for each category!
categories <- c("Alcohol","Diet_Weekly",
                "Smoking","Exercise_MET",
                "Exercise_Freq","Internet_Usage",
                "Deprivation_Indices","Vitamins")

protSelect <- lapply(categories, function(x){
  testR2_heatmap <- testR2_heatmap[testR2_heatmap[[x]] > 0.02,]
  as.character(testR2_heatmap[order(testR2_heatmap[[x]],
                                    decreasing = T),]$ID[1:30])
})
protSelect <- unique(unlist(protSelect))


testR2_heatmap <- testR2_heatmap %>%
  filter(ID %in% topEprot) %>%
  mutate(ID = factor(ID, levels = topEprot)) %>%
  arrange(ID) %>%
  column_to_rownames(var = "ID")

#Combine Exercise_MET and Exercise_Freq to just Exercise:
testR2_heatmap <- testR2_heatmap %>%
  mutate(Exercise = rowSums(across(c("Exercise_Freq","Exercise_MET"))))


library(circlize)
#Cluster Row and Columns
ht <- Heatmap(as.matrix(testR2_heatmap %>%
                          select(-c("Internet_Usage",
                                    "Deprivation_Indices",
                                    "Exercise_Freq",
                                    "Exercise_MET")) %>%
                          rename(Diet = "Diet_Weekly")),
              name = "R²",
              col = colorRamp2(c(0, 0.05, 0.1), c("blue", "white", "red")),
              #c("#440154", "#21908C", "#FDE725")
              width = ncol(testR2_heatmap)*unit(4, "mm"),
              height = nrow(testR2_heatmap)*unit(2.1, "mm"),
              show_row_names = TRUE,
              show_column_names = TRUE,
              cluster_rows = TRUE,
              cluster_columns = TRUE,
              show_column_dend = FALSE, 
              show_row_dend = FALSE,
              #row_split = 5,
              #row_km = 3,
              #column_gap = unit(8, "mm"),
              #row_gap = unit(20, "mm"),
              row_names_gp = gpar(fontsize = 7),
              column_names_gp = gpar(fontsize = 10),
              row_names_side = "left",
              heatmap_legend_param = list(
                title_gp = gpar(fontsize = 10), # Legend title font size
                labels_gp = gpar(fontsize = 10), # Legend labels font size
                legend_height = unit(2, "cm"), # Height of the legend
                legend_width = unit(1, "cm")   # Width of the legend
              )
)

draw(ht)

filename <- paste0("./Figures/HEAP/F1/Ecat_heatmap.png")
png(file=filename,
    height = 7, width = 4, units = "in", res = 1000)
draw(ht)
dev.off()

filename <- paste0("./Figures/HEAP/F1/Ecat_heatmap.svg")
svg(file=filename,
    height = 7, width = 4)
draw(ht)
dev.off()



###EXTRA::::
suppht <- as.matrix(testR2_heatmap %>%
            select(-c("Exercise_Freq",
                      "Exercise_MET")) %>%
            rename(Diet = "Diet_Weekly"))
ht <- Heatmap(suppht,
              name = "R2",
              col = colorRamp2(c(0, 0.04, 0.08), c("purple", "white", "orange")),
              width = ncol(testR2_heatmap)*unit(10, "mm"),
              height = nrow(testR2_heatmap)*unit(3, "mm"),
              show_row_names = TRUE,
              show_column_names = TRUE,
              cluster_rows = TRUE,
              cluster_columns = TRUE,
              show_column_dend = FALSE, 
              show_row_dend = FALSE,
              column_gap = unit(8, "mm"),
              row_names_gp = gpar(fontsize = 8),
              column_names_gp = gpar(fontsize = 8),
              row_names_side = "left",
              cell_fun = function(j, i, x, y, width, height, fill) {
                if (suppht[i, j] >= 0.02) {
                  grid.text(sprintf("%.2f", suppht[i, j]), x = x, y = y, 
                            gp = gpar(fontsize = 6, col = "black"))
                }
              },
              border = T,
              rect_gp = gpar(col = "black", lwd = 0.5)
)

filename <- paste0("./Figures/HEAP/A1/Ecat_heatmap.png")
png(file=filename,
    height = 10, width = 5, units = "in", res = 1000)
draw(ht)
dev.off()

filename <- paste0("./Figures/HEAP/A1/Ecat_heatmap.svg")
svg(file=filename,
    height = 10, width = 5)
draw(ht)
dev.off()
