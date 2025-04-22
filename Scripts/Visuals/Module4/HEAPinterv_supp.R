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
library(ComplexHeatmap)
library(circlize)

#Load INT structure:
setwd("/n/groups/patel/shakson_ukb/UK_Biobank/")
HEAPint <- qread("./Output/HEAPres/HEAPint.qs")
#sList <- HEAPint@sList; corList <-  HEAPint@cList; pList <- HEAPint@pList

HeatmapCSpecInt <- function(df_cor, df_pval, title){
  ht <- Heatmap(as.matrix(df_cor),
                name = "cor",
                col = colorRamp2(c(-1, -1e-20, 0, 1e-20, 1), c("purple","white", "darkgrey","white", "orange")),
                width = ncol(df_cor)*unit(10, "mm"),
                height = nrow(df_cor)*unit(5, "mm"),
                show_row_names = TRUE,
                show_column_names = TRUE,
                cluster_rows = TRUE,
                show_row_dend = FALSE,
                cluster_columns = FALSE,
                column_gap = unit(8, "mm"),
                row_names_gp = gpar(fontsize = 12),
                column_names_gp = gpar(fontsize = 12),
                row_names_side = "left",
                cell_fun = function(j, i, x, y, width, height, fill) {
                  pval <- df_pval[i, j]
                  if (!is.na(pval)) {
                    if (pval < 0.001) {
                      grid.text("∗∗∗", x, y)
                    } else if (pval < 0.01) {
                      grid.text("∗∗", x, y)
                    } else if (pval < 0.05) {
                      grid.text("∗", x, y)
                    }
                  }
                },
                # cell_fun = function(j, i, x, y, width, height, fill) {
                #   grid.text(sprintf("%.2f", UKBint_df2[i, j]), x = x, y = y, 
                #             gp = gpar(fontsize = 10, col = "black"))
                # },
                border = T,
                rect_gp = gpar(col = "black", lwd = 0.5)
  )
  
  draw(ht, 
       column_title=title,
       column_title_gp=grid::gpar(fontsize=16))
}

IntTypeHT <- function(corList, pList, label, rowLabels, title){
  cor <- corList %>% 
          reduce(full_join) %>%
          select(c(label,"ID","Type")) %>%
          pivot_wider(names_from = Type, values_from = !!sym(label)) %>%
          filter(rowSums(!is.na(.)) > 3) %>%
          column_to_rownames(var = "ID") %>%
          replace(is.na(.), 0) #replace any NAs with 0 to get matrix.
  
  p <- pList %>% 
          reduce(full_join) %>%
          select(c(label,"ID","Type")) %>%
          pivot_wider(names_from = Type, values_from = !!sym(label)) %>%
          filter(rowSums(!is.na(.)) > 3) %>%
          column_to_rownames(var = "ID")
  
  rownames(cor) <- rowLabels
  
  HeatmapCSpecInt(cor, p, title)
}

UKBlabels = c("Bread Type: Wholegrain","Cereal Intake",
              "Exercise: Swimming, Cycling, Keep Fit, Bowling",
              "Exercise: Strenuous Sports","White Bread Intake",
              "# Days/WK of Vigorous Activity",
              "Bread Intake","Exercise: Cycling",
              "Incr. Alcohol Intake (vs. 10yrs ago)",
              "Beef Intake 1x/WK","Poulty Intake: 2-4x/WK",
              "High IPAQ Activity Group","Summed Activity #Days/MTH",
              "Moderate Vigorous Activity","MET min/WK Vigorous Activity",
              "Never Smoked","Alcohol Intake: Rarely",
              "Summed MET min/WK All Activity",
              "Ever Smoked","Former Daily Smoker",
              "1-3x/MTH Alcohol Intake","Fresh Fruit Intake")
IntTypeHT(HEAPint@cList, HEAPint@pList, label = "HERITAGE_effect", rowLabels = UKBlabels, title = "HERITAGE")
IntTypeHT(HEAPint@cList, HEAPint@pList, label = "GLP1_effect1", rowLabels = UKBlabels, title = "GLP1 STEP1")
IntTypeHT(HEAPint@cList, HEAPint@pList, label = "GLP1_effect2", rowLabels = UKBlabels, title = "GLP1 STEP2")

save_plot <- function(formats, path, label, cList, pList, rowLabels, title) {
  lapply(formats, function(fmt) {
    # Construct parameters for the graphics device
    params <- list(file = paste0(path, ".", fmt), height = 6, width = 7)
    
    # Add the `res` parameter only for 'png'
    if (fmt == "png") {
      params$res <- 500
      params$units <- "in"
    }
    
    # Dynamically call the correct graphics device (png, svg, pdf)
    do.call(fmt, params)
    
    # Generate the plot
    IntTypeHT(cList, pList, label = label, rowLabels = rowLabels, title)
    
    # Close the device
    dev.off()
  })
}


# Save Plots of HEAP & Intervention Study Correlations across Model Specifications
save_plot(c("png","svg","pdf"), "./Figures/HEAP/SF5/HERITAGE_CSpec", 
          "HERITAGE_effect", HEAPint@cList, HEAPint@pList, UKBlabels, title = "HERITAGE")
save_plot(c("png","svg","pdf"), "./Figures/HEAP/SF5/GLP1_STEP1_CSpec", 
          "GLP1_effect1", HEAPint@cList, HEAPint@pList, UKBlabels, title = "GLP1 STEP1")
save_plot(c("png","svg","pdf"), "./Figures/HEAP/SF5/GLP1_STEP2_CSpec", 
          "GLP1_effect2", HEAPint@cList, HEAPint@pList, UKBlabels, title = "GLP1 STEP2")


## Additional Scatterplots - FIND PROTEINS WITH CONCORDANCE
library(ggplot2)
library(ggpubr)
library(ggpmisc)
library(ggrepel)
# Function to create plot for each study
create_plot_supp <- function(df, x_col, y_col, eName, effect_name, cor_val, p_val, se_col) {
  p <- df %>%
    filter(!is.na(Estimate)) %>%
    ggplot(aes(x = !!sym(x_col), y = !!sym(y_col))) +
    geom_point() +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray") + 
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray") + 
    # Use dynamic `se_col` for error bars
    geom_linerange(aes(ymin = !!sym(y_col) - 1.96 * !!sym(se_col), ymax = !!sym(y_col) + 1.96 * !!sym(se_col)), 
                   color = "black", size = 0.5, alpha = 0.5) +
    geom_linerange(aes(xmin = Estimate - 1.96 * `Std. Error`, xmax = Estimate + 1.96 * `Std. Error`), 
                   color = "black", size = 0.5, alpha = 0.5) +
    theme_minimal() +
    labs(x = paste0("Beta (", eName, ") - UKB"), y = paste0("Effect Size - ", effect_name)) + 
    geom_text_repel(aes(label = EntrezGeneSymbol), color = "red", size = 3, max.overlaps = 10)
  
  # Get axis limits
  plot_build <- ggplot_build(p)
  x_range <- plot_build$layout$panel_scales_x[[1]]$range$range
  y_range <- plot_build$layout$panel_scales_y[[1]]$range$range
  
  # Dynamically position the correlation label in the top-left corner
  x_pos <- x_range[1] + 0.15 * (x_range[2] - x_range[1])  # 15% from the left
  y_pos <- y_range[2] - 0.05 * (y_range[2] - y_range[1])  # 5% from the top
  
  # Add the annotation
  p1 <- p + annotate("text",
                     x = x_pos, 
                     y = y_pos, 
                     label = paste0("R = ", signif(cor_val, 2), ", p = ", signif(p_val, 2)), 
                     color = "black")
  
  lapply(c("png", "svg", "pdf"), 
         function(fmt) ggsave(paste0("./Figures/HEAP/SF5/",
                                     make.names(eName),"_", make.names(effect_name), ".", fmt), 
                              plot = p1, width = 6, height = 4, dpi = 500))
  
  return(p1)
}

# Function to create and save scatterplots
IntCor_Plot <- function(scatList, cList, pvalList, eID, eName, CSpec){
  #Filter the UKBscatter dataframe by eID:
  UKBscatter <- scatList[[CSpec]] %>% filter(ID == eID)
  
  # Prepare the data for plots
  UKBint_cor <- cList[[CSpec]]
  UKBint_p <- pvalList[[CSpec]]
  
  corVal_HERITAGE <- UKBint_cor %>% filter(ID == eID) %>% pull(HERITAGE_effect)
  corVal_GLP1.1 <- UKBint_cor %>% filter(ID == eID) %>% pull(GLP1_effect1)
  corVal_GLP1.2 <- UKBint_cor %>% filter(ID == eID) %>% pull(GLP1_effect2)
  
  pVal_HERITAGE <- UKBint_p %>% filter(ID == eID) %>% pull(HERITAGE_effect)
  pVal_GLP1.1 <- UKBint_p %>% filter(ID == eID) %>% pull(GLP1_effect1)
  pVal_GLP1.2 <- UKBint_p %>% filter(ID == eID) %>% pull(GLP1_effect2)
  
  # Specify standard error column for each study (adapt this based on the dataset)
  se_col_HERITAGE <- "HERITAGE_se"
  se_col_GLP1_STEP1 <- "GLP1_se1"
  se_col_GLP1_STEP2 <- "GLP1_se2"
  
  # Create plots
  gg1 <- create_plot_supp(UKBscatter, "Estimate", "HERITAGE_effect", eName, "HERITAGE", 
                          corVal_HERITAGE, pVal_HERITAGE, se_col_HERITAGE)
  gg2 <- create_plot_supp(UKBscatter, "Estimate", "GLP1_effect1", eName, "GLP1 STEP1", 
                          corVal_GLP1.1, pVal_GLP1.1, se_col_GLP1_STEP1)
  gg3 <- create_plot_supp(UKBscatter, "Estimate", "GLP1_effect2", eName, "GLP1 STEP2", 
                          corVal_GLP1.2, pVal_GLP1.2, se_col_GLP1_STEP2)
  
}

IntCor_Plot(scatList = HEAPint@sList,
            cList = HEAPint@cList,
            pvalList = HEAPint@pList, 
            eID = "past_tobacco_smoking_f1249_0_04",
            eName = "Former Daily Smoking", CSpec = "Model6")
