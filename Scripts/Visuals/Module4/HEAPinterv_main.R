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


#'*Load Univariate UKB HEAP results*
UKBint_dfcor <- HEAPint@cList$Model6 %>%
  select(-c("ID","Type"))
UKBint_dfpval <-  HEAPint@pList$Model6 %>%
  select(-c("ID","Type"))

rownames(UKBint_dfcor) <- c("White Bread Intake",
                            "Exercise: (Swimming, Cycling, etc.)",
                            "Incr. Alcohol Intake (vs. 10 yrs ago)",
                            "Cereal Intake",
                            "Former Daily Smoking", #"Past Tobacco Smoking (Everyday)"
                            "Fresh Fruit Intake",
                            "Smoked in Lifetime (Yes)",
                            "Never Smoked",
                            "Beef Intake 1x/WK",
                            "MET min/WK of Vigorous Activity",
                            "# Days/WK of Vigorous Activity",
                            "Exercise: Strenuous Sports")
colnames(UKBint_dfcor) <- c("HERITAGE","GLP1_STEP1","GLP1_STEP2")

ht <- Heatmap(as.matrix(UKBint_dfcor),
              name = "cor",
              col = colorRamp2(c(-1, 0, 1), c("purple", "white", "darkorange")),
              width = ncol(UKBint_dfcor)*unit(10, "mm"),
              height = nrow(UKBint_dfcor)*unit(5, "mm"),
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
                pval <- UKBint_dfpval[i, j]
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
              border = T,
              rect_gp = gpar(col = "black", lwd = 0.5)
)

draw(ht)

# Save plot:
filename <- paste0("./Figures/HEAP/F5/Eintven_ht.png")
png(file=filename,
    height = 4, width = 6, units = "in", res = 500)
draw(ht)
dev.off()

filename <- paste0("./Figures/HEAP/F5/Eintven_ht.svg")
svg(file=filename,
    height = 4, width = 6)
draw(ht)
dev.off()

filename <- paste0("./Figures/HEAP/F5/Eintven_ht.pdf")
pdf(file=filename,
    height = 4, width = 6)
draw(ht)
dev.off()

#SAVE results:
write.table(UKBint_dfcor, file = "./Output/SupplementaryTables/UKB_INT_cor.txt")
write.table(UKBint_dfpval, file = "./Output/SupplementaryTables/UKB_INT_pval.txt")


###### ScatterPlot to Compare Effects across UKB and Interventions ######
library(ggplot2)
library(ggpubr)
library(ggpmisc)
library(ggrepel)
# Function to create plot for each study
create_plot <- function(df, x_col, y_col, eName, effect_name, cor_val, p_val, se_col) {
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
         function(fmt) ggsave(paste0("./Figures/HEAP/F5/",
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
  gg1 <- create_plot(UKBscatter, "Estimate", "HERITAGE_effect", eName, "HERITAGE", 
                          corVal_HERITAGE, pVal_HERITAGE, se_col_HERITAGE)
  gg2 <- create_plot(UKBscatter, "Estimate", "GLP1_effect1", eName, "GLP1 STEP1", 
                          corVal_GLP1.1, pVal_GLP1.1, se_col_GLP1_STEP1)
  gg3 <- create_plot(UKBscatter, "Estimate", "GLP1_effect2", eName, "GLP1 STEP2", 
                          corVal_GLP1.2, pVal_GLP1.2, se_col_GLP1_STEP2)
  
}

IntCor_Plot(scatList = HEAPint@sList,
            cList = HEAPint@cList,
            pvalList = HEAPint@pList,
            eID = "fresh_fruit_intake_f1309_0_0",
            eName = "Fresh Fruit Intake",
            CSpec = "Model6")

IntCor_Plot(scatList = HEAPint@sList,
            cList = HEAPint@cList,
            pvalList = HEAPint@pList,
            eID = "number_of_days_week_of_vigorous_physical_activity_10_plus_minutes_f904_0_0",
            eName = "# Days/WK of Vigorous Activity",
            CSpec = "Model6")

IntCor_Plot(scatList = HEAPint@sList,
            cList = HEAPint@cList,
            pvalList = HEAPint@pList,
            eID = "types_of_physical_activity_in_last_4_weeks_f6164_0_0.multi_Strenuous_sports",
            eName = "Strenuous Sports",
            CSpec = "Model6")

IntCor_Plot(scatList = HEAPint@sList,
            cList = HEAPint@cList,
            pvalList = HEAPint@pList,
            eID = "past_tobacco_smoking_f1249_0_04",
            eName = "Former Daily Smoking", CSpec = "Model6")
