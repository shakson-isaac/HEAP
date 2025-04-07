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

#### General Correlations ####
UKB <- CovarSpecList[["Type6"]]$test[[1]] %>% 
          filter(`Pr(>|t|)` < 0.05/n()) %>%
          select(ID, omicID, Estimate) %>%
          pivot_wider(names_from = ID, values_from = Estimate)

#Number of entries in a columns must be greater than 3 (for correlation estimates)
UKBv2 <- UKB %>%
          select(where(~sum(!is.na(.)) >= 3))

#Merge the above dataframe with Intervention studies:
colnames(UKBv2)[which(names(UKBv2) == "omicID")] <- "EntrezGeneSymbol"

HERITAGE <- data %>% 
            filter(`False Discovery Rate (q-value)` < 0.05) %>%
            rename(HERITAGE_effect = `log(10) Fold Change`) %>%
            select(c(EntrezGeneSymbol, HERITAGE_effect))

GLP1_STEP1 <- STEP1 %>% 
                filter(qvalue < 0.05) %>% 
                group_by(EntrezGeneSymbol) %>%
                mutate(GLP1_effect1 = mean(effect_size),
                       GLP1_se1 = max(std_error),
                       GLP1_qvalue1 = max(qvalue)) %>%
                select(c(EntrezGeneSymbol, GLP1_effect1)) %>%
                unique()

GLP1_STEP2 <- STEP2 %>% 
                filter(qvalue < 0.05) %>% 
                group_by(EntrezGeneSymbol) %>%
                mutate(GLP1_effect2 = mean(effect_size),
                       GLP1_se2 = max(std_error),
                       GLP1_qvalue2 = max(qvalue)) %>%
                select(c(EntrezGeneSymbol, GLP1_effect2)) %>%
                unique()


UKBint <- list(UKBv2, HERITAGE, GLP1_STEP1, GLP1_STEP2) %>% reduce(full_join)

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
              # cell_fun = function(j, i, x, y, width, height, fill) {
              #   grid.text(sprintf("%.2f", UKBint_df2[i, j]), x = x, y = y, 
              #             gp = gpar(fontsize = 10, col = "black"))
              # },
              border = T,
              rect_gp = gpar(col = "black", lwd = 0.5)
)


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
write.table(UKBint_dfcor, file = "/n/groups/patel/shakson_ukb/UK_Biobank/Output/SupplementaryTables/UKB_INT_cor.txt")
write.table(UKBint_dfpval, file = "/n/groups/patel/shakson_ukb/UK_Biobank/Output/SupplementaryTables/UKB_INT_pval.txt")

#/n/groups/patel/shakson_ukb/UK_Biobank/Output/SupplementaryTables/UKB_INT_cor.txt
#/n/groups/patel/shakson_ukb/UK_Biobank/Output/SupplementaryTables/UKB_INT_pval.txt

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

# Get correlation and p-value for eID
get_correlation_pvalue <- function(df, eID, study) {
  cor_pval <- df %>%
    filter(eFeat %in% eID & Study == study) %>%
    select(cor, pval_adjust) %>%
    slice(1)
  return(list(cor = cor_pval$cor, pval = cor_pval$pval_adjust))
}

# Function to create and save scatterplots
IntCor_Plot <- function(eID, eName){
  UKBspec <- CovarSpecList[["Type6"]]$test[[1]] %>% 
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
  colnames(UKBint_df) <- c("eFeat","Study","pval_adjust","cor")
  corVal_HERITAGE <- get_correlation_pvalue(UKBint_df, eID, "HERITAGE_effect")
  corVal_GLP1_STEP1 <- get_correlation_pvalue(UKBint_df, eID, "GLP1_effect1")
  corVal_GLP1_STEP2 <- get_correlation_pvalue(UKBint_df, eID, "GLP1_effect2")
  
  # Specify standard error column for each study (adapt this based on the dataset)
  se_col_HERITAGE <- "HERITAGE_se"
  se_col_GLP1_STEP1 <- "GLP1_se1"
  se_col_GLP1_STEP2 <- "GLP1_se2"
  
  # Create plots
  gg1 <- create_plot(UKBscatter, "Estimate", "HERITAGE_effect", eName, "HERITAGE", corVal_HERITAGE$cor, corVal_HERITAGE$pval, se_col_HERITAGE)
  gg2 <- create_plot(UKBscatter, "Estimate", "GLP1_effect1", eName, "GLP1 STEP1", corVal_GLP1_STEP1$cor, corVal_GLP1_STEP1$pval, se_col_GLP1_STEP1)
  gg3 <- create_plot(UKBscatter, "Estimate", "GLP1_effect2", eName, "GLP1 STEP2", corVal_GLP1_STEP2$cor, corVal_GLP1_STEP2$pval, se_col_GLP1_STEP2)
  
}
IntCor_Plot(eID = "fresh_fruit_intake_f1309_0_0",
            eName = "Fresh Fruit Intake")

IntCor_Plot(eID = "number_of_days_week_of_vigorous_physical_activity_10_plus_minutes_f904_0_0",
            eName = "# Days/WK of Vigorous Activity")

IntCor_Plot(eID = "types_of_physical_activity_in_last_4_weeks_f6164_0_0.multi_Strenuous_sports",
            eName = "Strenuous Sports")


##### Network Tables to Build for Cytoscape ######

## Network:
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


#Get NODE1 - NODE2 in this script

UKBnet <- CovarSpecList[["Type6"]]$test[[1]] %>% 
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

#Significant effects in HERITAGE and GLP1
Interventions <- list(HERITAGE, GLP1_STEP1, GLP1_STEP2) %>% reduce(full_join)
colnames(Interventions)[which(names(Interventions) == "EntrezGeneSymbol")] <- "Protein"

Interventions <- Interventions %>% 
  select(HERITAGE_effect, GLP1_effect1, GLP1_effect2, Protein) %>%
  pivot_longer(cols = !Protein,
               names_to = "Intervention", 
               values_to = "Effect") %>%
  mutate(Effect = if_else(Effect > 0, 1, -1)) %>%
  filter(Protein %in% unique(UKBnet_select$Protein)) %>%
  na.omit()

UKB_IntNet <- list(UKBnet_select, Interventions) %>% reduce(full_join)







#Subset the above with mediation analysis results:
# Load Files w/ Progress Bar:
loadMDAssoc <- function(covarType){
  error_files <- list()  # Track files w/ errors
  stat_all <- list()
  
  # Use pblapply for progress bar and iteration
  for(i in 1:1000){
    tryCatch({
      load <- fread(file=paste0("/n/groups/patel/shakson_ukb/UK_Biobank/Data/Parallel/Prot_PGSPXS_mediation/",
                                covarType,"/MDres_",i,".txt"), 
                    select = 2:16)
      
      stat_all[[i]] <- load
      
    }, error = function(e) {
      message(paste("Error reading file:", i))
      error_files <<- append(error_files, i)
    })
  }
  
  
  stat_all <- do.call("rbind", stat_all)
  stat_all <- stat_all %>%
    mutate(E_HRi = exp(`Exposure Indirect Effect`),
           G_HRi = exp(`Genetic Indirect Effect`),
           delta_HRi = exp(abs(`Exposure Indirect Effect` - `Genetic Indirect Effect`)),
           total_HRi = exp(`Exposure Indirect Effect` + `Genetic Indirect Effect`),
           Eprop_mediated = abs(`Exposure Indirect Effect`)/(abs(`Exposure Indirect Effect`) + abs(`Exposure Direct Effect`)),
           Gprop_mediated = abs(`Genetic Indirect Effect`)/(abs(`Genetic Indirect Effect`) + abs(`Genetic Direct Effect`)))
  stat_all$CovarSpec <- covarType
  stat_all$HR_upper95 <- as.numeric(stat_all$HR_upper95)
  
  
  return(stat_all)
}

# Load results:
Type5 <- loadMDAssoc(covarType = "Type5")

T2D_assoc <- Type5 %>%
                filter(DZid == "age_e11_first_reported_non_insulin_dependent_diabetes_mellitus_f130708_0_0")

#Filter by nominal significance:
T2D_assoc1 <- T2D_assoc %>%
                  select(ID, HR, HR_lower95, HR_upper95,
                         cindex, E_HRi, G_HRi) %>%
                  filter((HR_lower95 > 1 & HR_upper95 > 1) | 
                           (HR_lower95 < 1 & HR_upper95 < 1))
                  
# Obtain OLS slope across Specifications:
oSlope <- function(MDres){
  
  # Obtain GEM: Expected Indirect Effect Slope: log(E/G)
  protContext <- pblapply(unique(MDres$ID), function(p){
    df <- MDres %>% filter(ID == p) 
    
    fit <- lm(`Exposure Indirect Effect` ~ `Genetic Indirect Effect`,
              data = df)
    
    stats <- summary(fit)$coefficients
    
    maxE <- max(abs(df$`Exposure Indirect Effect`))
    maxG <- max(abs(df$`Genetic Indirect Effect`))
    
    df <- MDres %>%
      filter(ID == p) %>%
      filter(!((HR_lower95 < 1 & HR_upper95 > 1) | (HR_lower95 > 1 & HR_upper95 < 1)))
    
    num_diseases <- length(unique(df$DZid))
    
    protStats <- data.frame(
      ID = p,
      Estimate = stats[2,"Estimate"],
      Std.Error = stats[2,"Std. Error"],
      Pvalue = stats[2,"Pr(>|t|)"],
      max.E = maxE,
      max.G = maxG,
      Estimate_Norm = stats[2,"Estimate"] * (maxE/maxG),
      NumDiseases = num_diseases
    )
    
    #Below return the slope of the above regression
    return(protStats) 
  })
  protContext <- do.call(rbind, protContext)
  
  return(protContext)
}
Type5_o <- oSlope(Type5)
GEMstats <- Type5_o %>%
              mutate(GEM = log(Estimate)) %>%
              select(ID, GEM) %>%
              na.omit()

T2D_ProtAttr <- list(T2D_assoc1, GEMstats) %>% reduce(left_join)

T2D_ProtAttr <- T2D_ProtAttr %>%
                    filter(HR > 1.2 | HR < 0.8)



#Subset UKB exposure effects for T2D specific effects:
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

#Add the Protein-T2D (links) rows:
T2Dlinkage <- UKB_IntNet_finalv2 %>%
                select(Protein, HR) %>%
                unique() %>%
                mutate(Effect = if_else(HR > 1, 1, -1)) %>%
                select(Protein, Effect) %>%
                rename(Intervention = Protein)
T2Dlinkage$Protein <- "T2D"
# T2Dlinkage <- merge(T2Dlinkage, UKB_IntNet_finalv2 %>% select(Protein,
#                                                               HR, HR_lower95, 
#                                                               HR_upper95,
#                                                               cindex, E_HRi,
#                                                               G_HRi, GEM),
#                     by.x = "Intervention", by.y = "Protein")


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
 
  




 
##### OLDER VERSION: TRASH ####
#Obtain dataframe w/ correlations and std errors of each study:
eID = "number_of_days_week_of_vigorous_physical_activity_10_plus_minutes_f904_0_0"
eName = "# Days/WK of Vigorous Activity"

eID = "types_of_physical_activity_in_last_4_weeks_f6164_0_0.multi_Strenuous_sports"
eName = "Strenuous Sports"

eID = "fresh_fruit_intake_f1309_0_0"
eName = "Fresh Fruit Intake"

UKBspec <- CovarSpecList[["Type6"]]$test[[1]] %>% 
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

library(ggplot2)
library(ggpubr)
library(ggpmisc)
library(ggrepel)

### PLOT HERITAGE vs UKB:
gg1 <- UKBscatter %>%
  filter(!is.na(Estimate)) %>%
  ggplot(aes(x = Estimate, y = HERITAGE_effect)) +
  geom_point() +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +  # Vertical line for x = 0
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +  # Horizontal line for y = 0
  # Vertical error bars (Confidence intervals for Estimate)
  geom_linerange(aes(ymin = HERITAGE_effect - 1.96 * HERITAGE_se, 
                     ymax = HERITAGE_effect + 1.96 * HERITAGE_se),
                 color = "black", size = 0.5,
                 alpha = 0.5) +  
  # Horizontal error bars (Confidence intervals for FC)
  geom_linerange(aes(xmin = Estimate - 1.96 * `Std. Error`, 
                     xmax = Estimate + 1.96 * `Std. Error`),
                 color = "black", size = 0.5,
                 alpha = 0.5) +
  theme_minimal() +
  labs(x = paste0("Beta (", eName,") - UKB"),
       y = "log10(FC) Post Exercise 20wks - HERITAGE",) + 
  geom_text_repel(aes(label = EntrezGeneSymbol), color = "red", 
                  size = 3, max.overlaps = 10)
gg1

# Extract the x and y axis limits from the plot
plot_build <- ggplot_build(gg1)
x_range <- plot_build$layout$panel_scales_x[[1]]$range$range
y_range <- plot_build$layout$panel_scales_y[[1]]$range$range

# Dynamically position the correlation label in the top-left corner
x_pos <- x_range[1] + 0.15 * (x_range[2] - x_range[1])  # 15% from the left
y_pos <- y_range[2] - 0.05 * (y_range[2] - y_range[1])  # 5% from the top

# Final plot with dynamically positioned correlation label using stat_cor()
# Adjusting the position and size of the correlation coefficient text

#Get correlation and adjusted pvalue for eID.
colnames(UKBint_df) <- c("eFeat","Study","pval_adjust","cor")

corVal <- UKBint_df %>%
  filter(eFeat %in% eID & Study %in% "HERITAGE_effect") %>%
  pull(cor)
pVal <- UKBint_df %>%
  filter(eFeat %in% eID & Study %in% "HERITAGE_effect") %>%
  pull(pval_adjust)

gg1 + annotate("text",
               x=x_pos, 
               y=y_pos, 
               label=paste0("R = ",signif(corVal,2),
                            ", ",
                            "p = ",signif(pVal,2)),
               color = "black")

### GLP1 STEP1
gg2 <- UKBscatter %>%
  filter(!is.na(Estimate)) %>%
  ggplot(aes(x = Estimate, y = GLP1_effect1)) +
  geom_point() +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +  # Vertical line for x = 0
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +  # Horizontal line for y = 0
  # Vertical error bars (Confidence intervals for Estimate)
  geom_linerange(aes(ymin = GLP1_effect1 - 1.96 * GLP1_se1, 
                     ymax = GLP1_effect1 + 1.96 * GLP1_se1),
                 color = "black", size = 0.5,
                 alpha = 0.5) +  
  # Horizontal error bars (Confidence intervals for FC)
  geom_linerange(aes(xmin = Estimate - 1.96 * `Std. Error`, 
                     xmax = Estimate + 1.96 * `Std. Error`),
                 color = "black", size = 0.5,
                 alpha = 0.5) +
  theme_minimal() +
  labs(x = paste0("Beta (", eName,") - UKB"),
       y = "Effect Size - GLP1 RCT",) + 
  geom_text_repel(aes(label = EntrezGeneSymbol), color = "red", 
                  size = 3, max.overlaps = 10)
gg2

# Extract the x and y axis limits from the plot
plot_build <- ggplot_build(gg2)
x_range <- plot_build$layout$panel_scales_x[[1]]$range$range
y_range <- plot_build$layout$panel_scales_y[[1]]$range$range

# Dynamically position the correlation label in the top-left corner
x_pos <- x_range[1] + 0.15 * (x_range[2] - x_range[1])  # 15% from the left
y_pos <- y_range[2] - 0.05 * (y_range[2] - y_range[1])  # 5% from the top

# Final plot with dynamically positioned correlation label using stat_cor()
# Adjusting the position and size of the correlation coefficient text

#Get correlation and adjusted pvalue for eID.
colnames(UKBint_df) <- c("eFeat","Study","pval_adjust","cor")

corVal <- UKBint_df %>%
  filter(eFeat %in% eID & Study %in% "GLP1_effect1") %>%
  pull(cor)
pVal <- UKBint_df %>%
  filter(eFeat %in% eID & Study %in% "GLP1_effect1") %>%
  pull(pval_adjust)

gg2 + annotate("text",
               x=x_pos, 
               y=y_pos, 
               label=paste0("R = ",signif(corVal,2),
                            ", ",
                            "p = ",signif(pVal,2)),
               color = "black")


#### GLP1 STEP2
gg3 <- UKBscatter %>%
  filter(!is.na(Estimate)) %>%
  ggplot(aes(x = Estimate, y = GLP1_effect2)) +
  geom_point() +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +  # Vertical line for x = 0
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +  # Horizontal line for y = 0
  geom_linerange(aes(ymin = GLP1_effect2 - 1.96 * GLP1_se2, ymax = GLP1_effect2 + 1.96 * GLP1_se2),
                 color = "black", size = 0.5, alpha = 0.5) +  # Vertical error bars
  geom_linerange(aes(xmin = Estimate - 1.96 * `Std. Error`, xmax = Estimate + 1.96 * `Std. Error`),
                 color = "black", size = 0.5, alpha = 0.5) +  # Horizontal error bars
  theme_minimal() +
  labs(x = paste0("Beta (", eName,") - UKB"),
       y = "Effect Size - GLP1 RCT") + 
  geom_text_repel(aes(label = EntrezGeneSymbol), color = "red", size = 3, max.overlaps = 10)

# Extract correlation and p-value for eID
cor_pval <- UKBint_df %>%
  filter(eFeat %in% eID & Study == "GLP1_effect2") %>%
  select(cor, pval_adjust) %>%
  slice(1)

# Extract axis limits and calculate label position
plot_build <- ggplot_build(gg3)
x_range <- plot_build$layout$panel_scales_x[[1]]$range$range
y_range <- plot_build$layout$panel_scales_y[[1]]$range$range
x_pos <- x_range[1] + 0.15 * (x_range[2] - x_range[1])
y_pos <- y_range[2] - 0.05 * (y_range[2] - y_range[1])

# Add correlation label to plot
gg3 + annotate("text",
               x = x_pos, 
               y = y_pos, 
               label = paste0("R = ", signif(cor_pval$cor, 2), ", p = ", signif(cor_pval$pval_adjust, 2)),
               color = "black")