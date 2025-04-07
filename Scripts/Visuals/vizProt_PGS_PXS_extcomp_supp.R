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

#TODO:
#'*Implement the correlation and pvalue calculations - simplified in the extcomp_final!*
#'*Copy and paste the code for correlation plots - it simplifies the process by using prior computations*
#'*Combine this and the final code into one script*


#'*Load Univariate UKB HEAP results*
setwd("/n/groups/patel/shakson_ukb/UK_Biobank/")
file.exists("RScripts/Pure_StatGen/Prot_ExPGS/resultsProtDS.R")
source("RScripts/Pure_StatGen/Prot_ExPGS/resultsProtDS.R")


#Intervention Studies:

#Compare results to:
#'*Exercise - JCI paper: https://pmc.ncbi.nlm.nih.gov/articles/PMC10132160/#sec13*
#'*GLP1 agonist treatment: https://doi.org/10.1038/s41591-024-03355-2*

# Read: first sheet w/ header start from line 3
file_path <- "/n/groups/patel/shakson_ukb/Motrpac/Related_Data/jciinsight_prot.xlsx"
data <- read_excel(file_path, sheet = excel_sheets(file_path)[1], skip = 2)
data$se <- data$`log(10) Fold Change`/data$`t-statistic`


# Read: first sheet w/ header start from line 3
file_path <- "/n/groups/patel/shakson_ukb/Motrpac/Related_Data/GLP1_proteomics.xlsx"
STEP1 <- read_excel(file_path, sheet = excel_sheets(file_path)[2], skip = 0)
STEP2 <- read_excel(file_path, sheet = excel_sheets(file_path)[3], skip = 0)

#### ^^Figure out how to aggregate info across all specifications: ####
# Function to process each data type (Type1 vs Type6)
process_data <- function(covar_spec_list, type) {
  data_v2 <- covar_spec_list[[type]]$test[[1]] %>% 
    filter(`Pr(>|t|)` < 0.05/n()) %>%
    select(ID, omicID, Estimate) %>%
    pivot_wider(names_from = ID, values_from = Estimate) %>%
    select(where(~sum(!is.na(.)) >= 3))
  
  # Merge with Intervention Studies
  colnames(data_v2)[which(names(data_v2) == "omicID")] <- "EntrezGeneSymbol"
  
  heritage <- data %>% 
    filter(`False Discovery Rate (q-value)` < 0.05) %>%
    rename(HERITAGE_effect = `log(10) Fold Change`) %>%
    select(c(EntrezGeneSymbol, HERITAGE_effect))
  
  glp1_step1 <- STEP1 %>% 
    filter(qvalue < 0.05) %>% 
    group_by(EntrezGeneSymbol) %>%
    mutate(GLP1_effect1 = mean(effect_size),
           GLP1_se1 = max(std_error),
           GLP1_qvalue1 = max(qvalue)) %>%
    select(c(EntrezGeneSymbol, GLP1_effect1)) %>%
    unique()
  
  glp1_step2 <- STEP2 %>% 
    filter(qvalue < 0.05) %>% 
    group_by(EntrezGeneSymbol) %>%
    mutate(GLP1_effect2 = mean(effect_size),
           GLP1_se2 = max(std_error),
           GLP1_qvalue2 = max(qvalue)) %>%
    select(c(EntrezGeneSymbol, GLP1_effect2)) %>%
    unique()
  
  # Merge all datasets
  merged_data <- list(data_v2, heritage, glp1_step1, glp1_step2) %>% reduce(full_join)
  
  return(merged_data)
}

# Function to calculate correlations and p-values
calculate_correlations <- function(merged_data) {
  
  UKBint <- merged_data
  UKBint_cor <- as.data.frame(cor(UKBint[,-1], use = "pairwise.complete.obs"))
  
  UKBint_cor <- UKBint_cor %>%
    select(c(HERITAGE_effect, GLP1_effect1,
             GLP1_effect2)) %>%
    na.omit()
  
  library(psych)
  UKBint_corv2 <- corr.test(UKBint[,-1], adjust = "none")
  UKBint_pval <- as.data.frame(UKBint_corv2$p) %>%
    select(c(HERITAGE_effect, GLP1_effect1,
             GLP1_effect2)) %>%
    na.omit()
  
  UKBint_pvalv2 <- sapply(UKBint_pval, function(x){
    p.adjust(x, method="BH")
  })
  UKBint_pvalv2 <- as.data.frame(UKBint_pvalv2)
  UKBint_pvalv2$eID <- rownames(UKBint_pval)
  
  
  UKBint_pvalv2 <- UKBint_pvalv2 %>%
    pivot_longer(cols = !eID,
                 names_to = "intervention", values_to = "pval_adjust") %>%
    filter(!(eID %in% c("HERITAGE_effect", 
                        "GLP1_effect1",
                        "GLP1_effect2")))
  head(UKBint_pvalv2)
  
  
  UKBint_corval <- as.data.frame(UKBint_corv2$r)
  UKBint_corval$eID <- rownames(UKBint_corval)
  
  UKBint_corval <- UKBint_corval %>%
    select(c(eID, HERITAGE_effect, GLP1_effect1,
             GLP1_effect2)) %>%
    pivot_longer(cols = !eID,
                 names_to = "intervention", values_to = "cor") %>%
    filter(!(eID %in% c("HERITAGE_effect", 
                        "GLP1_effect1",
                        "GLP1_effect2")))
  
  UKBint_df <- list(UKBint_pvalv2,UKBint_corval) %>% 
    reduce(full_join)
  
  selectIDs <- UKBint_df %>%
    filter(pval_adjust < 0.05) %>%
    pull(eID)
  
  UKBint_df <- UKBint_df %>% filter(eID %in% selectIDs)
  
  UKBint_dfcor <- UKBint_df %>%
    select(eID, intervention, cor) %>%
    pivot_wider(names_from = intervention, values_from = cor) %>%
    column_to_rownames(var = "eID")
  
  UKBint_dfpval <- UKBint_df %>%
    select(eID, intervention, pval_adjust) %>%
    pivot_wider(names_from = intervention, values_from = pval_adjust) %>%
    column_to_rownames(var = "eID")
  
  return(list(cor_values = UKBint_dfcor, pval_adjusted = UKBint_dfpval))
}

# Get the covariate specifications:
covar_spec_list <- CovarSpecList

# Code to Runt hrough all specifications and get correlations and pvalues:
corIntList <- lapply(names(CovarSpecList), function(x){
          df <- process_data(CovarSpecList, x)
          cor <- calculate_correlations(df)
          return(cor)
        })
names(corIntList) <- names(CovarSpecList)

#Get correlation and pvalue df combining all specifications:
corList <- lapply(names(corIntList), function(x){
  df <- corIntList[[x]]$cor_values
  df$ID <- rownames(df)
  df$Type <- x
  return(df)
})
pList <- lapply(names(corIntList), function(x){
  df <- corIntList[[x]]$pval_adjusted
  df$ID <- rownames(df)
  df$Type <- x
  return(df)
})
ans <-  corList %>% reduce(full_join)
ans2 <- pList %>% reduce(full_join)

library(ComplexHeatmap)
library(circlize)
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

IntTypeHT <- function(ans, ans2, label, rowLabels, title){
  cor <- ans %>%
    select(c(label,"ID","Type")) %>%
    pivot_wider(names_from = Type, values_from = !!sym(label)) %>%
    filter(rowSums(!is.na(.)) > 3) %>%
    column_to_rownames(var = "ID") %>%
    replace(is.na(.), 0) #replace any NAs with 0 to get matrix.
  print(rownames(cor))
  
  p <- ans2 %>%
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
IntTypeHT(ans, ans2, label = "HERITAGE_effect", rowLabels = UKBlabels, title = "HERITAGE")
IntTypeHT(ans, ans2, label = "GLP1_effect1", rowLabels = UKBlabels, title = "GLP1 STEP1")
IntTypeHT(ans, ans2, label = "GLP1_effect2", rowLabels = UKBlabels, title = "GLP1 STEP2")

save_plot <- function(formats, path, label, ans, ans2, rowLabels, title) {
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
    IntTypeHT(ans, ans2, label = label, rowLabels = rowLabels, title)
    
    # Close the device
    dev.off()
  })
}


# Example usage:
save_plot(c("png","svg","pdf"), "./Figures/HEAP/SF5/HERITAGE_CSpec", 
          "HERITAGE_effect", ans, ans2, UKBlabels, title = "HERITAGE")
save_plot(c("png","svg","pdf"), "./Figures/HEAP/SF5/GLP1_STEP1_CSpec", 
          "GLP1_effect1", ans, ans2, UKBlabels, title = "GLP1 STEP1")
save_plot(c("png","svg","pdf"), "./Figures/HEAP/SF5/GLP1_STEP2_CSpec", 
          "GLP1_effect2", ans, ans2, UKBlabels, title = "GLP1 STEP2")



#FIGURE OUT Labelling later:
## First, remove everything after the last '_0_' (including '_0_')
# x <- gsub("_f\\d+(_\\d+)*_0_", "_", ans$ID)
# 
# # Then, replace anything after the last '_0_' with just the digits (if non-zero)
# #y <- gsub("_0_00$", "", x)
# y <- gsub("_0$", "", x)
# z <- gsub(".*\\.multi", "", y)
# 
# 
# gsub("_f\\d+_\\d+_\\d+$", "", ans$ID)
# x <- gsub("_f\\d+(_\\d+)*$", "", ans$ID)
# x <- gsub(".*\\.multi", "",x)
# 
# gsub("_f\\d+(_\\d+)*", "", ans$ID)



#### FIND PROTEINS FOR INFLAMMATION ######
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
IntCor_Plot <- function(eID, eName, CSpec){
  UKBspec <- CovarSpecList[[CSpec]]$test[[1]] %>% 
    filter((ID == eID) & (`Pr(>|t|)` < 0.05/n())) %>%
    select(ID, omicID, Estimate, `Std. Error`)
  colnames(UKBspec)[which(names(UKBspec) == "omicID")] <- "EntrezGeneSymbol"
  
  HERITAGE <- data %>% 
    filter(`False Discovery Rate (q-value)` < 0.05) %>%
    mutate(HERITAGE_se = `log(10) Fold Change`/`t-statistic`) %>%
    rename(HERITAGE_effect = `log(10) Fold Change`) %>%
    select(c(EntrezGeneSymbol, HERITAGE_effect, HERITAGE_se))
  
  GLP1_STEP1 <- STEP1 %>% 
    filter(qvalue < 0.05) %>% 
    group_by(EntrezGeneSymbol) %>%
    mutate(GLP1_effect1 = mean(effect_size),
           GLP1_se1 = max(std_error)) %>%
    select(c(EntrezGeneSymbol, GLP1_effect1, GLP1_se1)) %>%
    unique()
  
  GLP1_STEP2 <- STEP2 %>% 
    filter(qvalue < 0.05) %>% 
    group_by(EntrezGeneSymbol) %>%
    mutate(GLP1_effect2 = mean(effect_size),
           GLP1_se2 = max(std_error)) %>%
    select(c(EntrezGeneSymbol, GLP1_effect2, GLP1_se2)) %>%
    unique()
  
  UKBscatter <- list(UKBspec, HERITAGE, GLP1_STEP1, GLP1_STEP2) %>% reduce(full_join)
  
  
  
  # Prepare the data for plots
  UKBint_cor <- corIntList[[CSpec]]$cor_values
  UKBint_p <- corIntList[[CSpec]]$pval_adjusted
  
  corVal_HERITAGE <- UKBint_cor[eID, "HERITAGE_effect"]
  corVal_GLP1.1 <- UKBint_cor[eID, "GLP1_effect1"]
  corVal_GLP1.2 <- UKBint_cor[eID, "GLP1_effect2"]
  
  pVal_HERITAGE <- UKBint_p[eID, "HERITAGE_effect"]
  pVal_GLP1.1 <- UKBint_p[eID, "GLP1_effect1"]
  pVal_GLP1.2 <- UKBint_p[eID, "GLP1_effect2"]
  
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

IntCor_Plot(eID = "past_tobacco_smoking_f1249_0_04",
            eName = "Former Daily Smoking", CSpec = "Type6")








###### PLOTTING ######

library(ComplexHeatmap)
library(circlize)
ht <- Heatmap(as.matrix(UKBint_dfcor),
              name = "cor",
              col = colorRamp2(c(-1, 0, 1), c("purple", "white", "orange")),
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
              # cell_fun = function(j, i, x, y, width, height, fill) {
              #   grid.text(sprintf("%.2f", UKBint_df2[i, j]), x = x, y = y, 
              #             gp = gpar(fontsize = 10, col = "black"))
              # },
              border = T,
              rect_gp = gpar(col = "black", lwd = 0.5)
)

draw(ht)

sort(unique(CovarSpecList[["Type6"]]$test[[1]]$ID))

rownames(UKBint_dfcor) <- c("White Bread Intake",
                            "Exercise: (Swimming, Cycling, etc.)",
                            "Incr. Alcohol Intake (vs. 10 yrs ago)",
                            "Cereal Intake",
                            "Past Tobacco Smoking (Everyday)",
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
              # cell_fun = function(j, i, x, y, width, height, fill) {
              #   grid.text(sprintf("%.2f", UKBint_df2[i, j]), x = x, y = y, 
              #             gp = gpar(fontsize = 10, col = "black"))
              # },
              border = T,
              rect_gp = gpar(col = "black", lwd = 0.5)
)


# Save plot:
filename <- paste0("./Figures/HEAP/SF5/Eintven_ht.png")
png(file=filename,
    height = 4, width = 6, units = "in", res = 500)
draw(ht)
dev.off()

filename <- paste0("./Figures/HEAP/SF5/Eintven_ht.svg")
svg(file=filename,
    height = 4, width = 6)
draw(ht)
dev.off()

filename <- paste0("./Figures/HEAP/SF5/Eintven_ht.pdf")
pdf(file=filename,
    height = 4, width = 6)
draw(ht)
dev.off()