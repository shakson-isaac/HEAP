#'*PLOT 5* Make a figure that shows these effects independent of the GvE plot



### Overall Average Plot across E: Categories
library(cowplot)
category_specificR2plot <- function(category, catR2, ycut){
  R2testCat_all <- merge(R2test, R2testCat, by = "ID")
  y_cutoff = ycut
  
  
  #Subset to label:
  subset_test <- R2testCat_all  %>% 
    filter(!ID %in% c(badfit_proteins)) %>%
    mutate_all(~replace_na(., 0)) %>%
    #filter(.data[[category]] > G) %>%
    filter(.data[[category]] > catR2) %>%
    top_n(10, .data[[category]]) 
  
  #Color by specific categories.
  p1 <- R2testCat_all %>% 
    filter(!ID %in% c(badfit_proteins)) %>%
    mutate_all(~replace_na(., 0)) %>% 
    ggplot(aes(x = G, y = .data[[category]], label = ID)) + 
    geom_point(aes(color = E)) + #aes(color = .data[[category]])) +
    geom_abline(slope=1, intercept = 0) +
    geom_ribbon(aes(ymin = 0, ymax = ifelse(G <= y_cutoff, G, y_cutoff), fill = "Genetics Driven"), alpha = 0.2) +
    geom_ribbon(aes(ymin = ifelse(G <= y_cutoff, G, y_cutoff), ymax = y_cutoff, fill = "Environment Driven"), alpha = 0.2) +
    scale_fill_manual(values = c("Genetics Driven" = "blue", "Environment Driven" = "green")) +
    ylim(c(0, y_cutoff)) +
    geom_text_repel(data=subset_test, aes(x = G,y = .data[[category]],label= ID),
                    size = 4) +
    #geom_label() + 
    scale_color_gradient(low = "darkslategrey", high = "red") +
    xlab(expression("G: PGS"~R^2)) +
    ylab(bquote(~.(category)~":"~R^2)) +
    theme_minimal() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) + #Adjust 90 to 0 for this plot:
    
    #Remove: the Environment vs Genetics Driven Legend for spacing:
    guides(fill = "none",  # Remove the fill legend
           color ="none"#guide_colorbar(title = bquote("PXS:" ~ R^2))  # Keep the color gradient legend
    ) +
    ggtitle(bquote("Test Set:"~.(category))) +
    theme(plot.title = element_text(hjust = 0.5, size = 14),
          legend.position = "bottom",
          legend.title = element_text(angle = 0),  # Optional: Set title orientation
          legend.text = element_text(angle = 90),  # Rotate legend text to vertical
          axis.text = element_text(size = 10),
          axis.title = element_text(size = 12)) #"right")
  
  return(p1)
}
get_legend <- function(){
  
}

g1 <- category_specificR2plot(category = "Exercise_Freq",catR2 = 0.02, ycut = 0.1)
g2 <- category_specificR2plot(category = "Exercise_MET",catR2 = 0.02, ycut = 0.05)
g3 <- category_specificR2plot(category = "Diet_Weekly",catR2 = 0.02, ycut = 0.075)
g4 <- category_specificR2plot(category = "Vitamins",catR2 = 0.01, ycut = 0.03)
g5 <- category_specificR2plot(category = "Alcohol",catR2 = 0.02, ycut = 0.05)
g6 <- category_specificR2plot(category = "Smoking",catR2 = 0.02, ycut = 0.05)
g7 <- category_specificR2plot(category = "Deprivation_Indices",catR2 = 0.01, ycut = 0.02)
g8 <- category_specificR2plot(category = "Internet_Usage",catR2 = 0.01, ycut = 0.02)
library(patchwork)
combined_plot <- (g1 | g2 | g3 | g4) /
  (g5 | g6 | g7 | g8)



#'*PLOT 4* Check the stability of R2 estimates across train and test:
partR2cat <- function(category){
  train <- R2trainCat %>% 
    rename(R2train = all_of(category)) %>%
    select(all_of(c("ID", "R2train")))
  test <- R2testCat %>% 
    rename(R2test = all_of(category)) %>%
    select(all_of(c("ID", "R2test")))
  R2cat <- merge(train, test, by = "ID")
  
  p1 <- R2cat %>%
    ggplot(aes(x = R2train, y = R2test)) +
    geom_point() +
    #geom_abline(slope=1, intercept = 0) +
    stat_poly_line() +
    stat_poly_eq() +
    theme_minimal() +
    xlab(bquote("Train" ~ R^2)) +
    ylab(bquote("Test" ~ R^2)) +
    ggtitle(bquote("Partitioned"~R^2~"-"~.(category)))
  
  return(p1)
}

g1 <- partR2cat(category = "Exercise_Freq")
g2 <- partR2cat(category = "Exercise_MET")
g3 <- partR2cat(category = "Diet_Weekly")
g4 <- partR2cat(category = "Vitamins")
g5 <- partR2cat(category = "Alcohol")
g6 <- partR2cat(category = "Smoking")
g7 <- partR2cat(category = "Deprivation_Indices")
g8 <- partR2cat(category = "Internet_Usage")

library(patchwork)
combined_plot <- (g1 | g2 | g3 | g4) /
  (g5 | g6 | g7 | g8)
# Save the plot to a file
ggsave("/n/groups/patel/shakson_ukb/UK_Biobank/Figures/protPGS_PXS_summ/R2_Categories.png", 
       plot = combined_plot, width = 15, height = 8, units = "in")

#'*^^REORDER the figures and categories later.*
#'*MAKE sure if you reorder the plot above to fix it for the plot here.*


#Remove Internet Usage + Deprivation Indices:
custom_colors <- colorRamp2(c(0, 0.04, 0.08), c("#440154", "#21908C", "#FDE725"))

library(ComplexHeatmap)
library(circlize)
library(grid)

# Assuming `testR2_heatmap` is your data matrix

ht <- Heatmap(
  as.matrix(testR2_heatmap),
  name = "R²",
  col = colorRamp2(c(0, 0.04, 0.08), c("blue", "white", "red")),  # Improved color gradient
  width = ncol(testR2_heatmap) * unit(10, "mm"),
  height = nrow(testR2_heatmap) * unit(5, "mm"),  # Adjusted height for better proportions
  show_row_names = TRUE,
  show_column_names = TRUE,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_column_dend = TRUE,  # Show dendrogram for context
  show_row_dend = TRUE,
  column_gap = unit(8, "mm"),
  row_names_gp = gpar(fontsize = 10, fontface = "bold"),  # Improved row label visibility
  column_names_gp = gpar(fontsize = 10, fontface = "bold"),  # Improved column label visibility
  row_names_side = "left",
  border = TRUE,
  rect_gp = gpar(col = "black", lwd = 0.5),
  cell_fun = function(j, i, x, y, width, height, fill) {
    value <- testR2_heatmap[i, j]
    if (value > 0.04) {
      grid.text(sprintf("%.2f", value), x = x, y = y, 
                gp = gpar(fontsize = 8, col = "black", fontface = "bold"))  # Text only for values > 0.04
    }
  },
  heatmap_legend_param = list(
    title = "R²",
    at = c(0, 0.04, 0.08),
    labels = c("Low", "Threshold", "High"),
    title_gp = gpar(fontsize = 12, fontface = "bold"),
    labels_gp = gpar(fontsize = 10)
  )
)

# Draw heatmap
draw(ht, heatmap_legend_side = "right")


# FUNCTION TO GO ACROSS ALL CATEGORIES AND PLOT E v G plot:
# FUNCTION TO GET LEGEND:
# & PROVIDE PLOTS FOR VARIABILITY of EACH CATEGORY.
