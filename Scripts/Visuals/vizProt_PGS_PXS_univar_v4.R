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

setwd("/n/groups/patel/shakson_ukb/UK_Biobank/")

# FUNCTIONS -------

# FUNCTION: Load Files w/ Progress Bar:
loadUniAssoc <- function(covarType){
  error_files <- list()  # Track files w/ errors
  stat_all <- list()
  train_data_list <- list()
  test_data_list <- list()
  
  # For loop through stored Univariate Associations in batches (400 batches here)
  for(i in 1:400){
    tryCatch({
      load <- readRDS(paste0("/n/groups/patel/shakson_ukb/UK_Biobank/Data/Parallel/Prot_PGSPXS_univar/",covarType,"/univar_assoc_",i,".rds"))
      
      train_data_list[[i]] <- load$train
      test_data_list[[i]]  <- load$test
    
    }, error = function(e) {
      message(paste("Error reading file:", i))
      error_files <<- append(error_files, i)
    })
  }
  
  # Store results
  stat_all[["train"]] <- train_data_list
  stat_all[["test"]] <- test_data_list
  
  return(stat_all)
}

# FUNCTION: List of list of dataframes to one list of dataframes:
combine_dataframes <- function(list_of_lists) {
  # Determine the number of elements in each inner list
  num_dfs <- length(list_of_lists[[1]])
  
  # For each position in the inner lists, bind the rows across all lists
  combined_dfs <- map(1:num_dfs, function(i) {
    bind_rows(map(list_of_lists, ~ .x[[i]]))
  })
  
  return(combined_dfs)
}

# FUNCTION: Check replication of associations:
replication_stats <- function(train, test, title){
  train_assoc <- train
  test_assoc <- test
  
  colnames(train_assoc)[!(colnames(train_assoc) %in% c("ID","omicID"))] <-  gsub(" ", "",
                                                                                 paste0(colnames(train_assoc)[!(colnames(train_assoc) %in% c("ID","omicID"))], "_train"))
  
  colnames(test_assoc)[!(colnames(test_assoc) %in% c("ID","omicID"))] <-  gsub(" ", "",
                                                                               paste0(colnames(test_assoc)[!(colnames(test_assoc) %in% c("ID","omicID"))], "_test"))
  
  
  all_assoc <- list(train_assoc, test_assoc) %>% reduce(full_join)
  
  pval_thresh <- 0.05/nrow(all_assoc)
  
  all_assoc_significant <- all_assoc %>%
    filter(`Pr(>|t|)_test` < pval_thresh &
             `Pr(>|t|)_train` < pval_thresh)
  
  #Example Plot: Are Estimates similar between train & test?
  gg1 <- all_assoc_significant %>%
    ggplot(aes(x = Estimate_train, y = Estimate_test)) +
    geom_point() +
    geom_smooth(method = "lm") +
    theme_minimal() +
    ggtitle(title)
  print(gg1)
  
  return(all_assoc_significant)
}

# FUNCTION: Miami Plot All Associations: (Train or Test Separately)
miamiplot_multi_omic <- function(assoc, title){
  #Arrange manhattan plot by E: Category
  assoc <- assoc %>%
    arrange(Category)
  
  #'*added to make sure the p-values are correctly shown in plot:*
  #assoc$`Pr(>|t|)` <- ifelse(assoc$`Pr(>|t|)` == 0, .Machine$double.xmin, assoc$`Pr(>|t|)`)
  
  
  ids <- gsub(":.*", "",assoc$ID) #accounts for ":" in the GxE ids
  
  num_feat <- length(unique(ids)) * length(unique(assoc$omicID))
  
  #get idx for a specific ordering
  assoc$idx <- as.numeric(factor(ids, levels = unique(ids)))
  
  #ggplot for exposures:
  p <- ggplot(assoc, aes(x = idx, y = -log10(`Pr(>|t|)` + 1e-300) * sign(Estimate), color = Category)) +
    geom_point() +
    geom_abline(intercept = -log10( (0.05/num_feat) + 1e-300), slope = 0, color = "blue", linetype = "dashed") +
    geom_abline(intercept = log10( (0.05/num_feat) + 1e-300), slope = 0, color = "blue", linetype = "dashed") +
    #facet_grid(rows = vars(CHR), scales = "free_x", space = "free_x") +
    #scale_y_continuous(trans = "reverse", limits = c(max(merged_df$logP1), 0)) +
    theme_minimal() +
    #theme(panel.spacing = unit(0.1, "lines")) +
    labs(x = "Exposures", y = bquote(-log[10]~"(P-value) * sign"~(beta)),
         color = "Category") + #expression(-log[10](P-value) * sign(Beta))) +
    ggtitle(title) +
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.y = element_text(size = 10),
          axis.title = element_text(size = 12),
          plot.title = element_text(hjust = 0.5, size = 14),
          legend.position = "right",
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 10)  
    )
  
  #print(p)
  ggsave(paste0("./Figures/HEAP/F2/",gsub(" ", "_",title),".png"), 
         plot = p, width = 12, height = 3, dpi = 400)
  #ggplotly(p, tooltip = "text")
}

# FUNCTION: Miami Plot Single Omic
#'*Decide whether to adjust significance to both thresholds - per protein bonferonni or all proteins*
miamiplot_single_omic <- function(protID, assoc, title){
  total <- assoc %>% 
              filter(omicID == protID) %>%
              arrange(Category)
  
  #'*added to make sure the p-values are correctly shown in plot:*
  #total$`Pr(>|t|)` <- ifelse(total$`Pr(>|t|)` == 0, .Machine$double.xmin, total$`Pr(>|t|)`)
  
  
  order <- 1:nrow(total)
  p <- ggplot(total, aes(x = order, y = -log10(`Pr(>|t|)` + 1e-300) * sign(Estimate), color = Category)) +
    geom_point(aes(text = paste("ID:", ID, 
                                "<br>Beta:", signif(Estimate, digits = 2),
                                "<br>Pval:", signif(`Pr(>|t|)`, digits = 2),
                                "<br>N:", samplesize,
                                "<br>Category:", Category,
                                "<br>Total R2:", signif(adj.R2, digits = 2)))) +
    geom_abline(intercept = -log10((0.05/nrow(total)) + 1e-300), slope = 0, 
                                        color = "blue", linetype = "dashed") +
    geom_abline(intercept = log10((0.05/nrow(total)) + 1e-300), slope = 0, 
                                        color = "blue", linetype = "dashed") +
    theme_minimal() +
    labs(x = "Exposures", y = "-log10(P-value) * sign(Beta)") +
    ggtitle(title)
  
  print(p)
  
  ggplotly(p, tooltip = "text")
  return(p)
}

# FUNCTION: Save Interactive Plot:
plotsave <- function(gg, protID, type){
  gginteract <- ggplotly(gg, tooltip = "text")
  
  saveWidget(gginteract, 
             file = paste0("./Figures/HEAP/CaseExamples/Interactive/",
                           protID,"_",type,"assoc.html"))
  
  ggsave(paste0("./Figures/HEAP/CaseExamples/Individual/",
                protID,"_",type,"assoc.png"), 
         plot = gg, width = 10, height = 6, dpi = 1000)
}


#  LOAD UNIVARIATE SUMMARY STATS -------
#'*LOADING Associations Across Specifications*
loadAssoc <- function(covarType){
  # Load specific Covar Spec:
  assoc <- loadUniAssoc(covarType)
  
  # Aggregating Train and Test Results:
  assoc$train <- combine_dataframes(assoc$train)
  assoc$test <- combine_dataframes(assoc$test)
  
  return(assoc)
}

Type1 <- loadAssoc("Type1")
Type2 <- loadAssoc("Type2")
Type3 <- loadAssoc("Type3")
Type4 <- loadAssoc("Type4")
Type5 <- loadAssoc("Type5")
Type6 <- loadAssoc("Type6")
Type7 <- loadAssoc("Type7")

# CovarSpec: All ASSOCIATIONS DataStructure:
CovarSpecList <- list(Type1, Type2, Type3, Type4, Type5, Type6)
names(CovarSpecList) <- c("Type1","Type2","Type3","Type4","Type5","Type6")

#Pull the R2 ranges of the univariate associations:
obtainR2_associations <- function(){
  Eprot <- CovarSpecList$Type6$test[[1]] %>%
    filter(`Pr(>|t|)` < 7e-8) %>%
    select(c("ID","omicID"))
  
  GxEprot <- CovarSpecList$Type6$test[[2]] %>%
    filter(`Pr(>|t|)` < 7e-8) %>%
    select(c("ID","omicID"))
  
  R2prot <- CovarSpecList$Type6$test[[3]] 
  
  
  # Find matching IDs and omicIDs across the three datasets
  merged_data1 <- Eprot %>%
    inner_join(R2prot, by = c("ID", "omicID"))
  
  merged_data2 <- GxEprot %>%
    inner_join(R2prot, by = c("ID", "omicID"))
  
  
  median(merged_data1$R2, na.rm = T)
  mean(merged_data1$R2, na.rm = T)
  max(merged_data1$R2, na.rm = T)
  summary(merged_data1$R2, na.rm = T)
  
  
  median(merged_data2$R2, na.rm = T)
  mean(merged_data2$R2, na.rm = T)
  max(merged_data2$R2, na.rm = T)
  summary(merged_data2$R2, na.rm = T)
  
}

# ANALYSIS and VISUALIZATIONS ------

#'*Miami Plot showing Replication Hits*

# FUNCTION: Obtain DF of replicated stats between Train and Test
replication_stats_final <- function(train, test){
  sigRES <- list()
  
  train_assoc <- train
  test_assoc <- test
  
  colnames(train_assoc)[!(colnames(train_assoc) %in% c("ID","omicID"))] <-  gsub(" ", "",
                                                                                 paste0(colnames(train_assoc)[!(colnames(train_assoc) %in% c("ID","omicID"))], "_train"))
  
  colnames(test_assoc)[!(colnames(test_assoc) %in% c("ID","omicID"))] <-  gsub(" ", "",
                                                                               paste0(colnames(test_assoc)[!(colnames(test_assoc) %in% c("ID","omicID"))], "_test"))
  
  
  all_assoc <- list(train_assoc, test_assoc) %>% reduce(full_join)
  
  #Add unique ID for each Exposure-Protein Pair:
  all_assoc$AssocID <- paste0(all_assoc$omicID,":", all_assoc$ID)
  
  
  
  pval_thresh <- 0.05/nrow(all_assoc) #Standard Bonferroni Correction for All ASSOCIATIONS. Ex. for 200 exposures, 3000 proteins: 0.05/(200*3000)
  
  allAssoc_significant <- all_assoc %>%
    filter(`Pr(>|t|)_test` < pval_thresh &
             `Pr(>|t|)_train` < pval_thresh)
  
  trainAssoc_significant <- all_assoc %>%
    filter(`Pr(>|t|)_train` < pval_thresh)
  
  testAssoc_significant <- all_assoc %>%
    filter(`Pr(>|t|)_test` < pval_thresh)
  
  sigRES[["sigTrain"]] <- trainAssoc_significant
  sigRES[["sigTest"]] <- testAssoc_significant
  sigRES[["sigBOTH"]] <- allAssoc_significant
  sigRES[["Bonf.correction"]] <- pval_thresh
  
  return(sigRES)
}

# FUNCTION: Replication Plot - Scatterplot (BASIC)
plot_replication <- function(assoc_df, title, fname){
  #Example Plot: Are Estimates similar between train & test?
  gg1 <- assoc_df %>%
    ggplot(aes(x = Estimate_train, y = Estimate_test)) +
    geom_point() +
    geom_smooth(method = "lm") +
    theme_minimal() +
    ggtitle(title)
  
  ggsave(paste0("./Figures/HEAP/SF3/",fname,".png"), 
         plot = gg1, width = 3, height = 3, dpi = 400)
}

# FUNCTION: DF for all stats between Train & Test:
replication_all <- function(train, test){
  train_assoc <- train
  test_assoc <- test
  
  colnames(train_assoc)[!(colnames(train_assoc) %in% c("ID","omicID"))] <-  gsub(" ", "",
                                                                                 paste0(colnames(train_assoc)[!(colnames(train_assoc) %in% c("ID","omicID"))], "_train"))
  
  colnames(test_assoc)[!(colnames(test_assoc) %in% c("ID","omicID"))] <-  gsub(" ", "",
                                                                               paste0(colnames(test_assoc)[!(colnames(test_assoc) %in% c("ID","omicID"))], "_test"))
  
  
  all_assoc <- list(train_assoc, test_assoc) %>% reduce(full_join)
  
  #Add unique ID for each Exposure-Protein Pair:
  all_assoc$AssocID <- paste0(all_assoc$omicID,":", all_assoc$ID)
  
  return(all_assoc)
}

# FUNCTION: Miami plot coloring associations replicated overlayed on Train Stats
miamiplot_multi_omic_rep <- function(assoc, title, filename){
  #Arrange manhattan plot by E: Category
  assoc <- assoc %>%
    arrange(Category_train)
  
  # Handle pvalue's that equal EXACTLY 0:
  #pval_min <- min(assoc$`Pr(>|t|)_train`[assoc$`Pr(>|t|)_train` > 0], na.rm = T)
  #assoc$`Pr(>|t|)_train` <- ifelse(assoc$`Pr(>|t|)_train` == 0, 
  #                                 pval_min, assoc$`Pr(>|t|)_train`)
  
  
  # Obtain significant results to plot first w/ COLOR!!
  pval_thresh <- 0.05/nrow(assoc) #Standard Bonferroni Correction for All ASSOCIATIONS. Ex. for 200 exposures, 3000 proteins: 0.05/(200*3000)
  
  
  # Get numeric ordering for x-axis
  ids <- gsub(":.*", "",assoc$ID) #accounts for ":" in the GxE ids
  
  num_feat <- length(unique(ids)) * length(unique(assoc$omicID))
  
  #get idx for a specific ordering
  assoc$idx <- as.numeric(factor(ids, levels = unique(ids)))
  
  # Define the factor in the combined dataframe
  assoc$Category_train <- factor(assoc$Category_train, 
                                    levels = unique(assoc$Category_train))
  
  # Subset Associations for Plotting:
  allAssoc_significant <- assoc %>%
    filter(`Pr(>|t|)_test` < pval_thresh &
             `Pr(>|t|)_train` < pval_thresh)
  
  allAssoc_nonsignificant <- assoc %>%
    filter(!(AssocID %in% allAssoc_significant$AssocID))
  
  
  #GGPLOT:
  #Added offset to -log10 to allow for inclusion of pvalues exactly = 0.
  p1 <- allAssoc_nonsignificant %>%
    ggplot(aes(x = idx, y = -log10(`Pr(>|t|)_train` + 1e-300) * sign(Estimate_train),
               color = Category_train)) +
    geom_point(size = 1, color = "gray") +
    scale_color_discrete(drop = FALSE)
  
  p2 <- p1 +
    geom_point(data = allAssoc_significant,
               size = 1) +
    geom_abline(intercept = -log10(pval_thresh + 1e-300), slope = 0, color = "blue", linetype = "dashed") +
    geom_abline(intercept = log10(pval_thresh + 1e-300), slope = 0, color = "blue", linetype = "dashed") +
    theme_minimal() +
    labs(x = "Exposures", y = bquote(-log[10]~"(P-value) * sign"~(beta)),
         color = "Category") + #expression(-log[10](P-value) * sign(Beta))) +
    ggtitle(title) +
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.y = element_text(size = 10),
          axis.title = element_text(size = 12),
          plot.title = element_text(hjust = 0.5, size = 14),
          legend.position = "right",
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 10)  
    )
  
  
  ggsave(paste0("./Figures/HEAP/F2/",gsub(" ", "_",filename),".png"), 
         plot = p2, width = 12, height = 3, dpi = 400)
}

# FUNCTION: Run through EVERY CovarSpec to Find Replicated Points
Replication_MiamiPlot <- function(CovarSpecAssocList){
  lapply(names(CovarSpecAssocList), function(x){
    # Obtain DF w/ both Train and Test Associations:
    Eall <- replication_all(train = CovarSpecList[[x]]$train[[1]],
                            test = CovarSpecList[[x]]$test[[1]])
    
    GxEall <- replication_all(train = CovarSpecList[[x]]$train[[2]],
                              test = CovarSpecList[[x]]$test[[2]])
    
    #Plot E and GxE Associations that are replicated:
    miamiplot_multi_omic_rep(Eall, "Entire Proteome: E Associations", 
                             filename = paste0("Eassoc_",x))
    
    miamiplot_multi_omic_rep(GxEall, "Entire Proteome: GxE Associations", 
                             filename = paste0("GxEassoc_",x))
  })
}


#'*MIAMI PLOTS of REPLICATION*
# CovarSpec: All Significant Assocations Datastructure
CovarSpecAssocList <- lapply(CovarSpecList, function(x){
  Assoc <- list()
  E <- replication_stats_final(x$train[[1]], x$test[[1]]) #E
  GxE <- replication_stats_final(x$train[[2]], x$test[[2]]) #GxE
  
  Assoc[["E"]] <- E
  Assoc[["GxE"]] <- GxE
  
  return(Assoc)
})

# Replication Miami Plots
#'*SWITCH CovarSpecAssocList --> CovarSpecList b/c it is MISLEADING HERE!!!!*
Replication_MiamiPlot(CovarSpecAssocList)

#Discovery Threshold: Defined as Threshold for significance in both Train and Test sets.
RunRepPlots<- function(){
  plot_replication(assoc_df = CovarSpecAssocList$Type1$E$sigBOTH, 
                   title = "Type 1: E Assoc Replication",
                   fname = "Type1_ReplicationCor")
  plot_replication(assoc_df = CovarSpecAssocList$Type2$E$sigBOTH, 
                   title = "Type 2: E Assoc Replication",
                   fname = "Type2_ReplicationCor")
  plot_replication(assoc_df = CovarSpecAssocList$Type3$E$sigBOTH, 
                   title = "Type 3: E Assoc Replication",
                   fname = "Type3_ReplicationCor")
  plot_replication(assoc_df = CovarSpecAssocList$Type4$E$sigBOTH, 
                   title = "Type 4: E Assoc Replication",
                   fname = "Type4_ReplicationCor")
  plot_replication(assoc_df = CovarSpecAssocList$Type5$E$sigBOTH, 
                   title = "Type 5: E Assoc Replication",
                   fname = "Type5_ReplicationCor")
  plot_replication(assoc_df = CovarSpecAssocList$Type6$E$sigBOTH, 
                   title = "Type 6: E Assoc Replication",
                   fname = "Type6_ReplicationCor")
}
RunRepPlots()
#'*Could make cor plot of above for next time*


#'*Understand Similarities and Differences across Covariate Specification*
# Function to Compare Association Concordance Across Specifications:
AssocCor_DF <- function(CovarSpecAssocList){
  repSTAT <- list()
  
  CovarSpecAssocCor <- lapply(names(CovarSpecAssocList), function(x){
    CovarSpecAssocList[[x]][["E"]][["sigBOTH"]]$CovarType <- x
    CovarSpecAssocList[[x]][["GxE"]][["sigBOTH"]]$CovarType <- x
    return(CovarSpecAssocList[[x]])
  })
  names(CovarSpecAssocCor) <- names(CovarSpecAssocList)
  
  # Aggregate Info across Specifications:
  Eassoc_AllType <- lapply(names(CovarSpecAssocCor), function(x){
    return(CovarSpecAssocCor[[x]][["E"]][["sigBOTH"]])
  }) %>% reduce(full_join)
  
  GxEassoc_AllType <- lapply(names(CovarSpecAssocCor), function(x){
    return(CovarSpecAssocCor[[x]][["GxE"]][["sigBOTH"]])
  }) %>% reduce(full_join)
  
  #Filter Associations Detected Across Specifications:
  
  #FIND THE margins for CIs
  # 95% confidence interval: use 1.96 as critical value since N > 100 for each association.
  Eassoc_AllType$CImargin_test <- Eassoc_AllType$Std.Error_test * 1.96
  GxEassoc_AllType$CImargin_test <- GxEassoc_AllType$Std.Error_test * 1.96
  #qt(0.975, df = samplesize - num_covars - 1)
  
  
  ###Dataframe for Association Tests:
  Eassoc_AllType_CI <- Eassoc_AllType %>%
    select(AssocID, CovarType, Estimate_test, Std.Error_test, CImargin_test, Category_test) %>% 
    pivot_wider(names_from = CovarType, values_from = c(Estimate_test, Std.Error_test, CImargin_test))
  
  
  GxEassoc_AllType_CI <- GxEassoc_AllType %>%
    select(AssocID, CovarType, Estimate_test, Std.Error_test, CImargin_test, Category_test) %>% 
    pivot_wider(names_from = CovarType, values_from = c(Estimate_test, Std.Error_test, CImargin_test))
  
  repSTAT[["E"]] <- Eassoc_AllType_CI
  repSTAT[["GxE"]] <- GxEassoc_AllType_CI
  
  return(repSTAT)
}

corDF <- AssocCor_DF(CovarSpecAssocList)

# Get STATS:
AllEassoc_cor <- corDF[["E"]]  %>%
                    select(all_of(c("Estimate_test_Type1","Estimate_test_Type2",
                                    "Estimate_test_Type3","Estimate_test_Type4",
                                    "Estimate_test_Type5","Estimate_test_Type6")))

AllGxEassoc_cor <- corDF[["GxE"]] %>%
                    select(all_of(c("Estimate_test_Type1","Estimate_test_Type2",
                                    "Estimate_test_Type3","Estimate_test_Type4",
                                    "Estimate_test_Type5","Estimate_test_Type6")))

# Function to check sign consistency across rows
corDF[["E"]]$Sign_Consistent <- apply(AllEassoc_cor, 1, function(row) {
  non_na_values <- row[!is.na(row)]     # Exclude NA values
  all(non_na_values > 0) || all(non_na_values < 0) # Check if all are positive or all are negative
})

corDF[["GxE"]]$Sign_Consistent <- apply(AllGxEassoc_cor , 1, function(row) {
  non_na_values <- row[!is.na(row)]     # Exclude NA values
  all(non_na_values > 0) || all(non_na_values < 0) # Check if all are positive or all are negative
})


# Proportion that have Consistent Sign across Specification
Esign_PropConsistent <- sum(corDF[["E"]]$Sign_Consistent)/length(corDF[["E"]]$Sign_Consistent)
GxEsign_PropConsistent <- sum(corDF[["GxE"]]$Sign_Consistent)/length(corDF[["GxE"]]$Sign_Consistent)

# Amount of Associations that are significant:
Ecount <- AllEassoc_cor %>%
  summarize(across(everything(), ~ sum(!is.na(.)), .names = "SigAssoc_{col}"))

GxEcount <- AllGxEassoc_cor %>%
  summarize(across(everything(), ~ sum(!is.na(.)), .names = "SigAssoc_{col}"))

HEAPcount <- rbind(Ecount,GxEcount)
HEAPcount$ID <- c("E","GxE")

fwrite(HEAPcount, file = "/n/groups/patel/shakson_ukb/UK_Biobank/Output/Univariate/ReplicationCount.txt")


# Function to get correlation of significant E estimates across Covariate Specifications:
corSpec_HT <- function(df, title_ht){
  cor_matrix <- cor(df,use='pairwise.complete.obs')
  
  # Shorten the row and column names
  rownames(cor_matrix) <- gsub("Estimate_test_", "", rownames(cor_matrix))
  colnames(cor_matrix) <- gsub("Estimate_test_", "", colnames(cor_matrix))
  
  # Convert the correlation matrix to a long format
  cor_long <- melt(cor_matrix) %>%
    filter(as.numeric(Var1) >= as.numeric(Var2))
  
  # Plot the heatmap
  g1 <- ggplot(cor_long, aes(x = Var1, y = Var2, fill = value)) +
    geom_tile(color = "white") + # Add grid lines for separation
    geom_text(aes(label = sprintf("%.3f", value)), size = 3) + # Add correlation values
    scale_fill_gradient2(low = "purple", high = "orange", mid = "white",
                         midpoint = 0.95, limit = c(0.90, 1),
                         name = "r") + # Narrow color range
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
      axis.title = element_blank()
    ) +
    labs(title = title_ht)
  return(g1)
}
corEplot <- corSpec_HT(AllEassoc_cor, "Correlation of Covariate Specifications: E")
corGplot <- corSpec_HT(AllGxEassoc_cor, "Correlation of Covariate Specifications: GxE")

ggsave(paste0("./Figures/HEAP/F2/ECorAcrossCSpec.png"), 
       plot = corEplot, width = 5, height = 3, dpi = 500)
ggsave(paste0("./Figures/HEAP/F3/GxECorAcrossCSpec.png"), 
       plot = corGplot, width = 5, height = 3, dpi = 500)

# Example: Correlation between Type 1 and Type 5 Specification.
library(ggrepel)
diffAssoc <- corDF[["E"]] %>%
  filter(Sign_Consistent == FALSE)

diffAssoc$AssocID <- gsub("_f\\d+_\\d_\\d{2}","",diffAssoc$AssocID)
diffAssoc$AssocID <- gsub(":"," ", diffAssoc$AssocID)
diffAssoc$AssocID <- gsub("_"," ", diffAssoc$AssocID)

# E Associations:
T1v5 <- corDF[["E"]] %>%
  ggplot(aes(x = Estimate_test_Type1, y = Estimate_test_Type5)) +
  geom_point(size = 1) +
  geom_linerange(aes(ymin = Estimate_test_Type5 - CImargin_test_Type5, 
                     ymax = Estimate_test_Type5 + CImargin_test_Type5)) +  # Vertical error bar
  geom_linerange(aes(xmin = Estimate_test_Type1 - CImargin_test_Type1, 
                     xmax = Estimate_test_Type1 + CImargin_test_Type1)) + # Horizontal error bar
  geom_text_repel(data = diffAssoc, # Only label the subset
                  aes(label = str_wrap(AssocID,
                      width = 10)),
                  max.overlaps = 2,
                  min.segment.length = unit(0, 'lines'),
                  nudge_y = -0.5,
                  nudge_x = 0.3,
                  size = 2.5) +
  labs(x = "Beta - Type1", y = "Beta - Type5") +
  theme_minimal()

ggsave(paste0("./Figures/HEAP/F2/T1v5.png"), 
       plot = T1v5, width = 4, height = 4, dpi = 500)



T1v5color <- corDF[["E"]]  %>%
  ggplot(aes(x = Estimate_test_Type1, y = Estimate_test_Type5, 
             color = Category_test)) +
  geom_point(size = 0.2, alpha = 1) +
  geom_linerange(aes(ymin = Estimate_test_Type5 - CImargin_test_Type5, 
                     ymax = Estimate_test_Type5 + CImargin_test_Type5),
                 alpha = 1,
                 size = 0.2) +  # Vertical error bar
  geom_linerange(aes(xmin = Estimate_test_Type1 - CImargin_test_Type1, 
                     xmax = Estimate_test_Type1 + CImargin_test_Type1),
                 alpha = 1,
                 size = 0.2) + # Horizontal error bar
  labs(x = "Beta - Type1", y = "Beta - Type5") +
  theme_minimal()

ggsave(paste0("./Figures/HEAP/F2/T1v5cat.png"), 
       plot = T1v5color, width = 6, height = 4, dpi = 500)




T1v5colorsplit <- corDF[["E"]]  %>%
  ggplot(aes(x = Estimate_test_Type1, y = Estimate_test_Type5,
             color = Category_test)) +
  geom_point(size = 0.2) +
  geom_linerange(aes(ymin = Estimate_test_Type5 - CImargin_test_Type5,
                     ymax = Estimate_test_Type5 + CImargin_test_Type5),
                 size = 0.2) +  # Vertical error bar
  geom_linerange(aes(xmin = Estimate_test_Type1 - CImargin_test_Type1,
                     xmax = Estimate_test_Type1 + CImargin_test_Type1),
                 size = 0.2) + # Horizontal error bar
  geom_smooth(method = "lm", se = TRUE,
              linewidth = 0.5) +
  theme_minimal() +
  facet_wrap(~ Category_test)

ggsave(paste0("./Figures/HEAP/F2/T1v5cat_split.png"), 
       plot = T1v5colorsplit, width = 10, height = 6, dpi = 500)


#### GxE Plots:
GxET1v5 <- corDF[["GxE"]]  %>%
  ggplot(aes(x = Estimate_test_Type1, y = Estimate_test_Type5)) +
  geom_point(size = 1) +
  geom_linerange(aes(ymin = Estimate_test_Type5 - CImargin_test_Type5, 
                     ymax = Estimate_test_Type5 + CImargin_test_Type5)) +  # Vertical error bar
  geom_linerange(aes(xmin = Estimate_test_Type1 - CImargin_test_Type1, 
                     xmax = Estimate_test_Type1 + CImargin_test_Type1)) + # Horizontal error bar
  labs(x = "Beta - Type1", y = "Beta - Type5") +
  xlim(-0.4, 0.4) +
  theme_minimal()
ggsave(paste0("./Figures/HEAP/F3/gxeT1v5.png"), 
       plot = GxET1v5, width = 3, height = 3, dpi = 500)


GxET1v5_color <- corDF[["GxE"]]  %>%
  ggplot(aes(x = Estimate_test_Type1, y = Estimate_test_Type5,
             color = Category_test)) +
  geom_point(size = 1) +
  geom_linerange(aes(ymin = Estimate_test_Type5 - CImargin_test_Type5, 
                     ymax = Estimate_test_Type5 + CImargin_test_Type5),
                 width = 0.1) +  # Vertical error bar
  geom_linerange(aes(xmin = Estimate_test_Type1 - CImargin_test_Type1, 
                     xmax = Estimate_test_Type1 + CImargin_test_Type1), 
                 height = 0.1) + # Horizontal error bar
  labs(x = "Beta - Type1", y = "Beta - Type5") +
  xlim(-0.4, 0.4) +
  theme_minimal()
ggsave(paste0("./Figures/HEAP/F3/gxeT1v5_cat.png"), 
       plot = GxET1v5_color, width = 5, height = 3, dpi = 500)


#'*DEALING WITH TRENDS - IDENTIFYING TRENDS OF SIGNFICANT and REPLICATED ASSOCIATIONS!*


#'# Function to Compare Association Concordance Across Specifications:
AssocTrends_DF <- function(CovarSpecAssocList){
  Trends <- list()
  
  CovarSpecAssocCor <- lapply(names(CovarSpecAssocList), function(x){
    CovarSpecAssocList[[x]][["E"]][["sigBOTH"]]$CovarType <- x
    CovarSpecAssocList[[x]][["GxE"]][["sigBOTH"]]$CovarType <- x
    return(CovarSpecAssocList[[x]])
  })
  names(CovarSpecAssocCor) <- names(CovarSpecAssocList)
  
  # Aggregate Info across Specifications:
  Eassoc_AllType <- lapply(names(CovarSpecAssocCor), function(x){
    return(CovarSpecAssocCor[[x]][["E"]][["sigBOTH"]])
  }) %>% reduce(full_join)
  
  GxEassoc_AllType <- lapply(names(CovarSpecAssocCor), function(x){
    return(CovarSpecAssocCor[[x]][["GxE"]][["sigBOTH"]])
  }) %>% reduce(full_join)
  
  Trends[["E"]] <- Eassoc_AllType
  Trends[["GxE"]] <- GxEassoc_AllType
  
  return(Trends)
}

HEAPtrends <- AssocTrends_DF(CovarSpecAssocList)

EassocCat_mean_abs <- HEAPtrends[["E"]] %>%
  group_by(Category_test, CovarType) %>%
  summarize(
    mean_effect_size = mean(abs(Estimate_test)),#weighted.mean(Estimate, w = 1 / `Std. Error`^2),
    se_mean_effect_size = sqrt(sum(`Std.Error_test`^2) / n()),
    n_observations = n(),
    .groups = 'drop'
  )

EassocCat_mean_abs %>%
  ggplot(aes(x=reorder(Category_test, -mean_effect_size), 
             y=mean_effect_size, 
             color = CovarType)) +
  geom_jitter(width = 0.3, height = 0.001) +
  ylim(0,0.4) +
  labs(
    x = NULL,
    y = "Absolute Mean Effect Size"
  ) +
  theme_minimal() +
  theme(axis.text.x=element_text(angle=60, hjust=1))



###GET Meta-Analysis Point Estimates of Effect Size for each Category!
library(metafor)

# Perform random-effects meta-analysis for each unique Category_test
TrendsRes <- EassocCat_mean_abs %>%
  group_by(Category_test) %>%
  summarise(
    Beta = as.numeric(rma(yi = mean_effect_size, 
                          sei = se_mean_effect_size, 
                          method = "REML", 
                          data = pick(everything()))$beta),
    SE = rma(yi = mean_effect_size, 
             sei = se_mean_effect_size, 
             method = "REML", 
             data = pick(everything()))$se,
    CI_lower =  rma(yi = mean_effect_size, 
                    sei = se_mean_effect_size, 
                    method = "REML", 
                    data = pick(everything()))$ci.lb,
    CI_upper = rma(yi = mean_effect_size, 
                   sei = se_mean_effect_size, 
                   method = "REML", 
                   data = pick(everything()))$ci.ub,
    .groups = 'drop'
  ) %>%
  rename(ID = Category_test)

# PLOT Averaged Main E Association Effects:
#qt(0.975, df = 5) - DO not do CIs here --> Not sure if VALID...

ECatTrends <- TrendsRes %>%
  ggplot(aes(x = reorder(ID, -Beta), y = Beta, color = ID)) +
  geom_point() +
  geom_linerange(aes(ymin = CI_lower, 
                     ymax = CI_upper)) +
  labs(x = NULL,
       y = "Absolute Mean Effect Size") +
  theme_minimal() +
  theme(axis.text.x=element_text(angle=60, hjust=1),
        legend.position="none")
ggsave(paste0("./Figures/HEAP/F2/EmeanabsTrends.png"), 
       plot = ECatTrends, width = 5, height = 4, dpi = 500)



#'*TRENDS for certain E features*
# REMEMBER THIS IS FOR THE SET OF ALL SIGNIFICANT ASSOCIATIONS!!!
Eassoc_mean_abs <- HEAPtrends[["E"]]  %>%
  group_by(ID, CovarType) %>%
  summarize(
    mean_effect_size = mean(abs(Estimate_test)),#weighted.mean(Estimate, w = 1 / `Std. Error`^2),
    se_mean_effect_size = sqrt(sum(`Std.Error_test`^2) / n()),
    n_observations = n(),
    .groups = 'drop'
  )

E_TrendsRes <- Eassoc_mean_abs %>%
  filter(grepl("alcohol_intake_frequency",ID)) %>%
  group_by(ID) %>%
  summarise(
    Beta = as.numeric(rma(yi = mean_effect_size, 
                          sei = se_mean_effect_size, 
                          method = "REML", 
                          data = pick(everything()))$beta),
    SE = rma(yi = mean_effect_size, 
             sei = se_mean_effect_size, 
             method = "REML", 
             data = pick(everything()))$se,
    CI_lower =  rma(yi = mean_effect_size, 
                    sei = se_mean_effect_size, 
                    method = "REML", 
                    data = pick(everything()))$ci.lb,
    CI_upper = rma(yi = mean_effect_size, 
                   sei = se_mean_effect_size, 
                   method = "REML", 
                   data = pick(everything()))$ci.ub,
    .groups = 'drop'
  )

AlcTrend <- E_TrendsRes %>%
  ggplot(aes(x = reorder(ID, Beta), y = Beta)) +
  geom_point() +
  geom_linerange(aes(ymin = CI_lower, 
                     ymax = CI_upper)) +
  labs(x = "Alcohol Intake",
       y = "Absolute Mean Effect Size") +
  scale_x_discrete(labels=c("alcohol_intake_frequency_f1558_0_02" = "Special Occasions", 
                            "alcohol_intake_frequency_f1558_0_03" = "1-3 Times Month",
                            "alcohol_intake_frequency_f1558_0_04" = "1-2 Times Week",
                            "alcohol_intake_frequency_f1558_0_05" = "3-4 Times Week",
                            "alcohol_intake_frequency_f1558_0_06" = "Daily/Almost Daily")) +
  theme_minimal() +
  theme(axis.text.x=element_text(angle=60, hjust=1),
        legend.position="none")
ggsave(paste0("./Figures/HEAP/F2/AlcTrends.png"), 
       plot = AlcTrend, width = 5, height = 4, dpi = 500)



#'*SINGLE Omic Signatures:*
saveOmicSignatures <- function(omicID, CovarType, plotname){
  gg1 <- miamiplot_single_omic(protID = omicID, 
                               assoc = CovarSpecList[[CovarType]]$train[[1]],
                               title = paste0(omicID,": ",plotname))
  plotsave(gg1, omicID, CovarType)
} 
saveOmicSignatures("GDF15","Type5","E Associations")
saveOmicSignatures("CXCL17","Type5","E Associations")
saveOmicSignatures("MAMDC4","Type5","E Associations")
saveOmicSignatures("LAMP3","Type5","E Associations")
saveOmicSignatures("WFDC2","Type5","E Associations")
saveOmicSignatures("CEACAM16","Type5","E Associations")
saveOmicSignatures("LEP","Type5","E Associations")
saveOmicSignatures("IGFBP2","Type5","E Associations")
saveOmicSignatures("IGFBP1","Type5","E Associations")
saveOmicSignatures("CKB","Type5","E Associations")
saveOmicSignatures("FABP4","Type5","E Associations")

#'*FIGURE out how to highlight replicated hits later.*

#'*Specific Associations E or GxE*
#Source univarInteract Script:
source("/n/groups/patel/shakson_ukb/UK_Biobank/RScripts/Pure_StatGen/Prot_ExPGS/runProt_PGS_univarInteract_v2.R")

#'*To Modify Later --> For # of days,etc. that got converted into z-score.*
#'*Discretize the above by days and find the respective quintiles.*

#Run specific associations:
saveGWISplot <- function(o,e,c,n,l=c(),p,f){
  gg <- runGWISplot(omic = o,
                    exposure = e,
                    CType = c,
                    ename = n,
                    levels = l,
                    plotname = p)
  ggsave(paste0("./Figures/HEAP/CaseExamples/GWIS/",paste0(o,"_",f),".png"), 
         plot = gg, width = 5, height = 4, dpi = 500)
}

#Find top E hits that REPLICATE to show
Type5Eres <- CovarSpecAssocList$Type5$E$sigBOTH
sort(unique(Type5Eres$ID))

runGWISplot(omic = "IGFBP2",
            exposure = "pork_intake_f1389_0_0",
            CType = "Type5",
            ename = "Pork Intake",
            plotname = "Exposure Effect Conditioned on Genetics")


runGWISplot(omic = "CXCL17",
            exposure = "current_tobacco_smoking_f1239_0_0_Yes._on_most_or_all_days",
            CType = "Type5",
            ename = "Smoking Status",
            plotname = "Exposure Effect Conditioned on Genetics")

#Diet Related Proteins
#Hard to Find Solid Evidence for this based on trends
#Increase in exposure is not equivalent to increase in protein expr


### Alcohol Related Proteins
saveGWISplot(o = "MAMDC4",
             e = "alcohol_drinker_status_f20117_0_0_Current",
             c = "Type5",
             n = "Alcohol User",
             l = c("NO","YES"),
             p = "Exposure Effect Conditioned on Genetics",
             f = "AlcoholUse_E")
saveGWISplot(o = "MAMDC4",
             e = "alcohol_intake_frequency_f1558_0_0",
             c = "Type5",
             n = "Alcohol Freq.",
             l = c("Never","Special Occasions",
                        "1-3 times monthly","1-2 times weekly",
                        "3-4 times weekly","Daily"),
             p = "Exposure Effect Conditioned on Genetics",
             f = "AlcoholFreq_E")

saveGWISplot(o = "CEACAM16",
             e = "alcohol_drinker_status_f20117_0_0_Current",
             c = "Type5",
             n = "Alcohol User",
             l = c("NO","YES"),
             p = "Exposure Effect Conditioned on Genetics",
             f = "AlcoholUse_E")
saveGWISplot(o = "CEACAM16",
             e = "alcohol_intake_frequency_f1558_0_0",
             c = "Type5",
             n = "Alcohol Freq.",
             l = c("Never","Special Occasions",
                   "1-3 times monthly","1-2 times weekly",
                   "3-4 times weekly","Daily"),
             p = "Exposure Effect Conditioned on Genetics",
             f = "AlcoholFreq_E")





### Exercise Related proteins
saveGWISplot(o = "LEP",
             e = "types_of_physical_activity_in_last_4_weeks_f6164_0_0.multi_Strenuous_sports",
             c = "Type5",
             n = "Strenuous Sports",
             l = c("NO","YES"),
             p = "Exposure Effect Conditioned on Genetics",
             f = "StrenSports_E")
saveGWISplot(o = "LEP",
             e = "usual_walking_pace_f924_0_0",
             c = "Type5",
             n = "Walking Pace",
             #l = c("NO","YES"),
             p = "Exposure Effect Conditioned on Genetics",
             f = "WalkPace_E")
saveGWISplot(o = "LEP",
             e = "met_minutes_per_week_for_vigorous_activity_f22039_0_0" ,
             c = "Type5",
             n = "Vigorous Activity \n (Z-score)",
             #l = c("NO","YES"),
             p = "Exposure Effect Conditioned on Genetics",
             f = "Min_VigorousPhysAct_E")
saveGWISplot(o = "LEP",
             e = "number_of_days_week_of_vigorous_physical_activity_10_plus_minutes_f904_0_0",
             c = "Type5",
             n = "Vigorous Activity \n (Z-score)",
             #l = c("NO","YES"),
             p = "Exposure Effect Conditioned on Genetics",
             f = "Days_VigorousPhysAct_E")


saveGWISplot(o = "FABP4",
             e = "types_of_physical_activity_in_last_4_weeks_f6164_0_0.multi_Strenuous_sports",
             c = "Type5",
             n = "Strenuous Sports",
             l = c("NO","YES"),
             p = "Exposure Effect Conditioned on Genetics",
             f = "StrenSports_E")
saveGWISplot(o = "FABP4",
             e = "number_of_days_week_of_vigorous_physical_activity_10_plus_minutes_f904_0_0",
             c = "Type5",
             n = "Vigorous Activity \n (Z-score)",
             #l = c("NO","YES"),
             p = "Exposure Effect Conditioned on Genetics",
             f = "Days_VigorousPhysAct_E")


saveGWISplot(o = "IGFBP1",
             e = "number_of_days_week_of_vigorous_physical_activity_10_plus_minutes_f904_0_0",
             c = "Type5",
             n = "Vigorous Activity \n (Z-score)",
             #l = c("NO","YES"),
             p = "Exposure Effect Conditioned on Genetics",
             f = "Days_VigorousPhysAct_E")
saveGWISplot(o = "IGFBP1",
             e = "summed_days_activity_f22033_0_0",
             c = "Type5",
             n = "Days of Activity \n (Z-score)",
             #l = c("NO","YES"),
             p = "Exposure Effect Conditioned on Genetics",
             f = "Days_SummedDaysPhysAct_E")


### Smoking Related Proteins:
saveGWISplot(o = "CXCL17",
             e = "current_tobacco_smoking_f1239_0_0_Yes._on_most_or_all_days",
             c = "Type5",
             n = "Daily Smoker",
             l = c("NO","YES"),
             p = "Exposure Effect Conditioned on Genetics",
             f = "DailySmoker_E")

saveGWISplot(o = "LAMP3",
             e = "current_tobacco_smoking_f1239_0_0_Yes._on_most_or_all_days",
             c = "Type5",
             n = "Daily Smoker",
             l = c("NO","YES"),
             p = "Exposure Effect Conditioned on Genetics",
             f = "DailySmoker_E")

saveGWISplot(o = "TNR",
             e = "current_tobacco_smoking_f1239_0_0_Yes._on_most_or_all_days",
             c = "Type5",
             n = "Daily Smoker",
             l = c("NO","YES"),
             p = "Exposure Effect Conditioned on Genetics",
             f = "DailySmoker_E")

saveGWISplot(o = "WFDC2",
             e = "current_tobacco_smoking_f1239_0_0_Yes._on_most_or_all_days",
             c = "Type5",
             n = "Daily Smoker",
             l = c("NO","YES"),
             p = "Exposure Effect Conditioned on Genetics",
             f = "DailySmoker_E")

#Find top GxE hits that REPLICATE to show!!
Type5GxEres <- CovarSpecAssocList$Type5$GxE$sigBOTH

saveGWISplot(o = "APOF",
             e = "alcohol_intake_frequency_f1558_0_0",
             c = "Type5",
             n = "Alcohol Intake Freq.",
             l = c("Never","Special Occasions",
                        "1-3 times monthly","1-2 times weekly",
                        "3-4 times weekly","Daily"),
             p = "Polygenic GxE Interaction",
             f = "Alcohol_GxE")
saveGWISplot(o = "ALPP",
             e = "current_tobacco_smoking_f1239_0_0_Yes._on_most_or_all_days",
             c = "Type5",
             n = "Daily Smoker",
             l = c("NO","YES"),
             p = "Polygenic GxE Interaction",
             f = "Smoker_GxE")



#'*Save CSV file of replicated associations in Output Folder!*

# E Associations:
fwrite(CovarSpecAssocList$Type1$E$sigBOTH,
       file = "/n/groups/patel/shakson_ukb/UK_Biobank/Output/SupplementaryTables/CovarSpec1EassocReplicated.txt")
fwrite(CovarSpecAssocList$Type2$E$sigBOTH,
       file = "/n/groups/patel/shakson_ukb/UK_Biobank/Output/SupplementaryTables/CovarSpec2EassocReplicated.txt")
fwrite(CovarSpecAssocList$Type3$E$sigBOTH,
       file = "/n/groups/patel/shakson_ukb/UK_Biobank/Output/SupplementaryTables/CovarSpec3EassocReplicated.txt")
fwrite(CovarSpecAssocList$Type4$E$sigBOTH,
       file = "/n/groups/patel/shakson_ukb/UK_Biobank/Output/SupplementaryTables/CovarSpec4EassocReplicated.txt")
fwrite(CovarSpecAssocList$Type5$E$sigBOTH,
       file = "/n/groups/patel/shakson_ukb/UK_Biobank/Output/SupplementaryTables/CovarSpec5EassocReplicated.txt")
fwrite(CovarSpecAssocList$Type6$E$sigBOTH,
       file = "/n/groups/patel/shakson_ukb/UK_Biobank/Output/SupplementaryTables/CovarSpec6EassocReplicated.txt")

# GxE Associations:
fwrite(CovarSpecAssocList$Type1$GxE$sigBOTH,
       file = "/n/groups/patel/shakson_ukb/UK_Biobank/Output/SupplementaryTables/CovarSpec1GxEassocReplicated.txt")
fwrite(CovarSpecAssocList$Type2$GxE$sigBOTH,
       file = "/n/groups/patel/shakson_ukb/UK_Biobank/Output/SupplementaryTables/CovarSpec2GxEassocReplicated.txt")
fwrite(CovarSpecAssocList$Type3$GxE$sigBOTH,
       file = "/n/groups/patel/shakson_ukb/UK_Biobank/Output/SupplementaryTables/CovarSpec3GxEassocReplicated.txt")
fwrite(CovarSpecAssocList$Type4$GxE$sigBOTH,
       file = "/n/groups/patel/shakson_ukb/UK_Biobank/Output/SupplementaryTables/CovarSpec4GxEassocReplicated.txt")
fwrite(CovarSpecAssocList$Type5$GxE$sigBOTH,
       file = "/n/groups/patel/shakson_ukb/UK_Biobank/Output/SupplementaryTables/CovarSpec5GxEassocReplicated.txt")
fwrite(CovarSpecAssocList$Type6$GxE$sigBOTH,
       file = "/n/groups/patel/shakson_ukb/UK_Biobank/Output/SupplementaryTables/CovarSpec6GxEassocReplicated.txt")

########## PULL CODE FROM BELOW ############################################################

##'*QQplot:*
CovarSpecList$Type1$train[[1]] %>%
  ggplot(aes(sample = -log10(`Pr(>|t|)`))) +
    stat_qq() + 
    stat_qq_line()

CovarSpecList$Type1$test[[1]] %>%
  ggplot(aes(sample = -log10(`Pr(>|t|)`))) +
  stat_qq() + 
  stat_qq_line()

CovarSpecList$Type5$train[[1]] %>%
  ggplot(aes(sample = -log10(`Pr(>|t|)`))) +
  stat_qq() + 
  stat_qq_line()

#Function from: https://slowkow.com/notes/ggplot2-qqplot/
gg_qqplot <- function(ps, ci = 0.95) {
  n  <- length(ps)
  df <- data.frame(
    observed = -log10(sort(ps)),
    expected = -log10(ppoints(n)),
    clower   = -log10(qbeta(p = (1 - ci) / 2, shape1 = 1:n, shape2 = n:1)),
    cupper   = -log10(qbeta(p = (1 + ci) / 2, shape1 = 1:n, shape2 = n:1))
  )
  log10Pe <- expression(paste("Expected -log"[10], plain(P)))
  log10Po <- expression(paste("Observed -log"[10], plain(P)))
  ggplot(df) +
    geom_ribbon(
      mapping = aes(x = expected, ymin = clower, ymax = cupper),
      alpha = 0.1
    ) +
    geom_point(aes(expected, observed), shape = 1, size = 3) +
    geom_abline(intercept = 0, slope = 1, alpha = 0.5) +
    # geom_line(aes(expected, cupper), linetype = 2, size = 0.5) +
    # geom_line(aes(expected, clower), linetype = 2, size = 0.5) +
    xlab(log10Pe) +
    ylab(log10Po)
}

gg_qqplot(CovarSpecList$Type1$test[[1]]$`Pr(>|t|)`) +
  theme_bw(base_size = 24) +
  theme(
    axis.ticks = element_line(size = 0.5),
    panel.grid = element_blank()
    # panel.grid = element_line(size = 0.5, color = "grey80")
  )

##qqplot:
inflation <- function(ps) {
  chisq <- qchisq(1 - ps, 1)
  lambda <- median(chisq) / qchisq(0.5, 1)
  lambda
}

CovarSpecList$Type1$test[[1]] %>%
  filter(omicID == "GDF15") %>%
  pull(`Pr(>|t|)`) %>%
  gg_qqplot() + 
  theme_bw(base_size = 24) +
  theme(
    axis.ticks = element_line(size = 0.5),
    panel.grid = element_blank()
    # panel.grid = element_line(size = 0.5, color = "grey80")
  )

CovarSpecList$Type1$test[[1]] %>%
  filter(omicID == "GDF15") %>%
  pull(`Pr(>|t|)`) %>%
  inflation()




gg_qqplot(CovarSpecList$Type6$test[[1]]$`Pr(>|t|)`) +
  theme_bw(base_size = 24) +
  theme(
    axis.ticks = element_line(size = 0.5),
    panel.grid = element_blank()
    # panel.grid = element_line(size = 0.5, color = "grey80")
  )


gg_qqplot(CovarSpecList$Type6$test[[2]]$`Pr(>|t|)`) +
  theme_bw(base_size = 24) +
  theme(
    axis.ticks = element_line(size = 0.5),
    panel.grid = element_blank()
    # panel.grid = element_line(size = 0.5, color = "grey80")
  )


inflation <- function(ps) {
  chisq <- qchisq(1 - ps, 1)
  lambda <- median(chisq) / qchisq(0.5, 1)
  lambda
}
set.seed(1234)
pvalue <- runif(1000, min=0, max=1)
inflation(pvalue)


Einflated <- c(inflation(CovarSpecList$Type1$test[[1]]$`Pr(>|t|)`),
               inflation(CovarSpecList$Type2$test[[1]]$`Pr(>|t|)`),
               inflation(CovarSpecList$Type3$test[[1]]$`Pr(>|t|)`),
               inflation(CovarSpecList$Type4$test[[1]]$`Pr(>|t|)`),
               inflation(CovarSpecList$Type5$test[[1]]$`Pr(>|t|)`),
               inflation(CovarSpecList$Type6$test[[1]]$`Pr(>|t|)`))
names(Einflated) <- names(CovarSpecList)

GxEinflated <- c(inflation(CovarSpecList$Type1$test[[2]]$`Pr(>|t|)`),
                  inflation(CovarSpecList$Type2$test[[2]]$`Pr(>|t|)`),
                  inflation(CovarSpecList$Type3$test[[2]]$`Pr(>|t|)`),
                  inflation(CovarSpecList$Type4$test[[2]]$`Pr(>|t|)`),
                  inflation(CovarSpecList$Type5$test[[2]]$`Pr(>|t|)`),
                  inflation(CovarSpecList$Type6$test[[2]]$`Pr(>|t|)`))
names(GxEinflated) <- names(CovarSpecList)


#GENOMIC CONTROL: of E associations:
#https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-1131-9#Sec9

CovarSpecList$Type1$test[[1]]$Tval_Adj <- CovarSpecList$Type1$test[[1]]$`t value`/2.809583
#Get pvalues;

# Recalculate corrected p-values
CovarSpecList$Type1$test[[1]]$p_corrected <- 2 * (1 - pt(abs(CovarSpecList$Type1$test[[1]]$Tval_Adj), 
                                                         df = CovarSpecList$Type1$test[[1]]$samplesize - 3 - 1))
###'*NEXT time get DF of model in summary stats!!*

# TODO: ffjfjfjfj 



#View(CovarSpecList$Type1$test[[1]])

#GOAL --> BY APPROPRIATE CONDITIONING ARE WE LOWERING THE INFLATION FACTOR of ASSOCIATIONS!!!


#Trends: from PLOT
# Test Set Threshold NOT well calibrated. 
# Some 'outlier' or behaviors occur considering only train set threshold.
# Should attach identifiers to these plots that behave differently.
# ^^FOR ABOVE DISCOVERY THRESHOLD MAKES THE MOST SENSE!!!

# OTHER IDEAS:
# Table of number of significance hits across specification.

replication_stats_metrics <- function(train, test){
  sigRES <- list()
  
  train_assoc <- train
  test_assoc <- test
  
  colnames(train_assoc)[!(colnames(train_assoc) %in% c("ID","omicID"))] <-  gsub(" ", "",
                                                                                 paste0(colnames(train_assoc)[!(colnames(train_assoc) %in% c("ID","omicID"))], "_train"))
  
  colnames(test_assoc)[!(colnames(test_assoc) %in% c("ID","omicID"))] <-  gsub(" ", "",
                                                                               paste0(colnames(test_assoc)[!(colnames(test_assoc) %in% c("ID","omicID"))], "_test"))
  
  
  all_assoc <- list(train_assoc, test_assoc) %>% reduce(full_join)
  
  pval_thresh <- 0.05/nrow(all_assoc) #Standard Bonferroni Correction for All ASSOCIATIONS. Ex. for 200 exposures, 3000 proteins: 0.05/(200*3000)
  
  allAssoc_significant <- all_assoc %>%
    filter(`Pr(>|t|)_test` < pval_thresh &
             `Pr(>|t|)_train` < pval_thresh)
  
  trainAssoc_significant <- all_assoc %>%
    filter(`Pr(>|t|)_train` < pval_thresh)
  
  testAssoc_significant <- all_assoc %>%
    filter(`Pr(>|t|)_test` < pval_thresh)
  
  sigRES[["sigTrain"]] <- trainAssoc_significant
  sigRES[["sigTest"]] <- testAssoc_significant
  sigRES[["sigBOTH"]] <- allAssoc_significant
  sigRES[["Bonf.correction"]] <- pval_thresh
  
  return(sigRES)
}


CovarSpec1_Euni <- replication_stats_metrics(Type1$train[[1]], Type1$test[[1]])
CovarSpec1_GxEuni <- replication_stats_metrics(Type1$train[[2]], Type1$test[[2]])

#### COUNTS: COUNT # of Significant Associations:

# Do difference in Density:
Type1_OmicRes <- CovarSpec1_Euni$sigBOTH %>% count(omicID)
Type5_OmicRes <- CovarSpec5_Euni$sigBOTH %>% count(omicID)

Type1_OmicRes %>%
    ggplot(aes(x = n)) +
  geom_density()

Type5_OmicRes %>%
  ggplot(aes(x = n)) +
  geom_density()


#Discovery: Significance in both Train and Test!
Type1_OmicRes %>%
  slice_max(n, n = 10) %>%
  ggplot(aes(x = reorder(omicID, -n), y = n)) +
    geom_bar(stat='identity') +
    #geom_point() +
    labs(x = "Proteins",
         y = "# of Significant Hits (Discovery)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  

Type5_OmicRes %>%
  slice_max(n, n = 10) %>%
  ggplot(aes(x = reorder(omicID, -n), y = n)) +
  geom_bar(stat='identity') +
  #geom_point() +
  labs(x = "Proteins",
       y = "# of Significant Hits (Discovery)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))



CovarSpec1_Euni$sigBOTH %>%
  count(Category_test)
CovarSpec5_Euni$sigBOTH %>%
  count(Category_test)


Type1_numAssoc <- CovarSpec1_Euni$sigBOTH %>%
  count(Category_test) %>%
  summarise(across(n, sum)) %>%
  pull(n)

Type5_numAssoc <- CovarSpec5_Euni$sigBOTH %>%
  count(Category_test) %>%
  summarise(across(n, sum)) %>%
  pull(n)

####Understand Difference Between Covariate Specifications: Correlation of Associations etc.

# CovarSpec: All ASSOCIATIONS DataStructure:
CovarSpecList <- list(Type1, Type2, Type3, Type4, Type5, Type6)
names(CovarSpecList) <- c("Type1","Type2","Type3","Type4","Type5","Type6")

# CovarSpec: All Significant Assocations Datastructure:
CovarSpecAssocList <- lapply(CovarSpecList, function(x){
  Assoc <- list()
  E <- replication_stats_final(x$train[[1]], x$test[[1]]) #E
  GxE <- replication_stats_final(x$train[[2]], x$test[[2]]) #GxE
  
  Assoc[["E"]] <- E
  Assoc[["GxE"]] <- GxE
  
  return(Assoc)
})

# Function to Compare Association Concordance Across Specifications:
CovarSpecAssocCor <- lapply(names(CovarSpecList), function(x){
  CovarSpecAssocList[[x]][["E"]][["sigBOTH"]]$CovarType <- x
  CovarSpecAssocList[[x]][["GxE"]][["sigBOTH"]]$CovarType <- x
  
  #Create Another Column to Highlight Unique Prot:Eid Associations
  Eassoc <- CovarSpecAssocList[[x]][["E"]][["sigBOTH"]]
  GxEassoc <- CovarSpecAssocList[[x]][["GxE"]][["sigBOTH"]]
  
  CovarSpecAssocList[[x]][["E"]][["sigBOTH"]]$AssocID <- paste0(Eassoc$omicID,":",Eassoc$ID)
  CovarSpecAssocList[[x]][["GxE"]][["sigBOTH"]]$AssocID <- paste0(GxEassoc$omicID,":",GxEassoc$ID)
  
  
  return(CovarSpecAssocList[[x]])
})
names(CovarSpecAssocCor) <- names(CovarSpecList)

#^^TO MITIGATE ABOVE SHOULD ADD THE COLUMN OF INTEREST BEFORE THIS STEP!!!


# Aggregate Info across Specifications:

CovarSpecAssocCor$Type1$E$sigBOTH

Eassoc_AllType <- lapply(names(CovarSpecAssocCor), function(x){
                          return(CovarSpecAssocCor[[x]][["E"]][["sigBOTH"]])
                        }) %>% reduce(full_join)

GxEassoc_AllType <- lapply(names(CovarSpecAssocCor), function(x){
                          return(CovarSpecAssocCor[[x]][["GxE"]][["sigBOTH"]])
                        }) %>% reduce(full_join)

#Filter Assocations Detected Across Specifications:

#FIND THE margins for CIs

# 95% confidence interval:
Eassoc_AllType$CImargin_test <- Eassoc_AllType$Std.Error_test * qt(0.975, 
                                                                   df = Eassoc_AllType$samplesize_test - 1)
GxEassoc_AllType$CImargin_test <- GxEassoc_AllType$Std.Error_test * qt(0.975, 
                                                                       df = GxEassoc_AllType$samplesize_test - 1)


###Dataframe for Association Tests:
Eassoc_AllType_CI <- Eassoc_AllType %>%
  select(AssocID, CovarType, Estimate_test, Std.Error_test, CImargin_test) %>% 
  pivot_wider(names_from = CovarType, values_from = c(Estimate_test, Std.Error_test, CImargin_test))


GxEassoc_AllType_CI <- GxEassoc_AllType %>%
  select(AssocID, CovarType, Estimate_test, Std.Error_test, CImargin_test) %>% 
  pivot_wider(names_from = CovarType, values_from = c(Estimate_test, Std.Error_test, CImargin_test))


Eassoc_AllType_CI %>%
  ggplot(aes(x = Estimate_test_Type1, y = Estimate_test_Type5)) +
  geom_point(size = 1) +
  geom_linerange(aes(ymin = Estimate_test_Type5 - CImargin_test_Type5, 
                     ymax = Estimate_test_Type5 + CImargin_test_Type5),
                 width = 0.05) +  # Vertical error bar
  geom_linerange(aes(xmin = Estimate_test_Type1 - CImargin_test_Type1, 
                     xmax = Estimate_test_Type1 + CImargin_test_Type1), 
                 height = 0.05) + # Horizontal error bar
  theme_minimal()

Eassoc_AllType_CI %>%
  ggplot(aes(x = Estimate_test_Type1, y = Estimate_test_Type2)) +
  geom_point(size = 1) +
  geom_linerange(aes(ymin = Estimate_test_Type2 - CImargin_test_Type2, 
                     ymax = Estimate_test_Type2 + CImargin_test_Type2),
                 width = 0.05) +  # Vertical error bar
  geom_linerange(aes(xmin = Estimate_test_Type1 - CImargin_test_Type1, 
                     xmax = Estimate_test_Type1 + CImargin_test_Type1), 
                 height = 0.05) + # Horizontal error bar
  theme_minimal()


#Examples: Thinking about if Categories of Exposures Hold Differing Patterns?!?!
Eassoc_AllType_CIv2 <- Eassoc_AllType %>%
  select(AssocID, CovarType, Estimate_test, Std.Error_test, CImargin_test,Category_test) %>% 
  pivot_wider(names_from = CovarType, values_from = c(Estimate_test, Std.Error_test, CImargin_test))

Eassoc_AllType_CIv2 %>%
  ggplot(aes(x = Estimate_test_Type1, y = Estimate_test_Type5, 
             color = Category_test)) +
  geom_point(size = 1) +
  geom_linerange(aes(ymin = Estimate_test_Type5 - CImargin_test_Type5, 
                     ymax = Estimate_test_Type5 + CImargin_test_Type5),
                 width = 0.05) +  # Vertical error bar
  geom_linerange(aes(xmin = Estimate_test_Type1 - CImargin_test_Type1, 
                     xmax = Estimate_test_Type1 + CImargin_test_Type1), 
                 height = 0.05) + # Horizontal error bar
  theme_minimal()


Eassoc_AllType_CIv2 %>%
  ggplot(aes(x = Estimate_test_Type1, y = Estimate_test_Type5, 
             color = Category_test)) +
  geom_point(size = 1) +
  geom_linerange(aes(ymin = Estimate_test_Type5 - CImargin_test_Type5, 
                     ymax = Estimate_test_Type5 + CImargin_test_Type5),
                 width = 0.05) +  # Vertical error bar
  geom_linerange(aes(xmin = Estimate_test_Type1 - CImargin_test_Type1, 
                     xmax = Estimate_test_Type1 + CImargin_test_Type1), 
                 height = 0.05) + # Horizontal error bar
  #geom_smooth(method = "lm", se = TRUE) + 
  theme_minimal() +
  facet_wrap(~ Category_test)

library(lme4)
Eassoc_AllType_CIv2 <- Eassoc_AllType_CIv2 %>%
                        filter(!(Category_test %in% c("Internet_Usage")))
                        
fit <- lmer(Estimate_test_Type5 ~ Estimate_test_Type1*Category_test + (1+Estimate_test_Type1|Category_test),
            data = Eassoc_AllType_CIv2)
summary(fit)

model <- lm(Estimate_test_Type5 ~ Estimate_test_Type1 * Category_test, data = Eassoc_AllType_CIv2)
summary(model)

model_no_interaction <- lm(Estimate_test_Type5 ~ Estimate_test_Type1 + Category_test, data = Eassoc_AllType_CIv2)
anova(model_no_interaction, model) #LRT of if Slopes (Correlation of Effects by CovarSpec is Dependent on Category of Exposure)


# Plot with regression lines per category
Eassoc_AllType_CIv2 %>%
  ggplot(aes(x = Estimate_test_Type1, y = Estimate_test_Type5, color = Category_test)) +
  geom_point(size = 1) +
  geom_line(aes(y = Predicted), size = 1) +  # Add predicted regression lines
  theme_minimal() +
  labs(title = "Differing Slopes by Category with Predicted Values",
       x = "Estimate Type 1",
       y = "Estimate Type 5")

Eassoc_AllType_CIv2 %>%
  ggplot(aes(x = Estimate_test_Type1, y = Estimate_test_Type5, 
             color = Category_test)) +
  geom_point(size = 0.05) +
  geom_linerange(aes(ymin = Estimate_test_Type5 - CImargin_test_Type5, 
                     ymax = Estimate_test_Type5 + CImargin_test_Type5),
                 width = 0.03) +  # Vertical error bar
  geom_linerange(aes(xmin = Estimate_test_Type1 - CImargin_test_Type1, 
                     xmax = Estimate_test_Type1 + CImargin_test_Type1), 
                 height = 0.03) + # Horizontal error bar
  #geom_smooth(method = "lm", se = TRUE) + 
  geom_abline(intercept = 0, slope = 1) +
  theme_minimal() +
  facet_wrap(~ Category_test)






#DO somesort of gaussian mixture model (or other form of clustering)
#Kmeans this:

AssocClusters <- na.omit(Eassoc_AllType_CI %>%
                          select(all_of(c("Estimate_test_Type1",
                                          "Estimate_test_Type5"))))
kmeans_result <- kmeans(AssocClusters, 
                        centers = 6)

# Add cluster information to the data
AssocClusters$Cluster <- as.factor(kmeans_result$cluster)

# Visualize the clusters using ggplot
AssocClusters %>%
  ggplot(aes(x = Estimate_test_Type1, 
             y = Estimate_test_Type5)) +
  #geom_point(size = 1) +
  geom_density_2d()

AssocClusters %>%
  ggplot(aes(x = Estimate_test_Type1, 
             y = Estimate_test_Type5, color = Cluster)) +
  geom_point(size = 1) +
  labs(
    title = "K-Means Clustering Visualization",
    x = "X Variable",
    y = "Y Variable",
    color = "Cluster"
  ) +
  theme_minimal()

AssocClusters %>%
  ggplot(aes(x = Estimate_test_Type1, 
             y = Estimate_test_Type5, color = Cluster)) +
  #geom_point(size = 1) +
  geom_density_2d() +
  labs(
    title = "K-Means Clustering Visualization",
    x = "X Variable",
    y = "Y Variable",
    color = "Cluster"
  ) +
  theme_minimal()





Eassoc_AllType_CI %>%
  ggplot(aes(x = Estimate_test_Type4, y = Estimate_test_Type5)) +
  geom_point(size = 1) +
  geom_linerange(aes(ymin = Estimate_test_Type5 - CImargin_test_Type5, 
                     ymax = Estimate_test_Type5 + CImargin_test_Type5),
                 width = 0.05) +  # Vertical error bar
  geom_linerange(aes(xmin = Estimate_test_Type2 - CImargin_test_Type2, 
                     xmax = Estimate_test_Type2 + CImargin_test_Type2), 
                 height = 0.05) + # Horizontal error bar
  theme_minimal()


GxEassoc_AllType_CI %>%
  ggplot(aes(x = Estimate_test_Type1, y = Estimate_test_Type5)) +
  geom_point(size = 1) +
  geom_linerange(aes(ymin = Estimate_test_Type5 - CImargin_test_Type5, 
                     ymax = Estimate_test_Type5 + CImargin_test_Type5),
                 width = 0.1) +  # Vertical error bar
  geom_linerange(aes(xmin = Estimate_test_Type1 - CImargin_test_Type1, 
                     xmax = Estimate_test_Type1 + CImargin_test_Type1), 
                 height = 0.1) + # Horizontal error bar
  theme_minimal() +
  geom_smooth(method = "lm")


cor(Eassoc_AllType_CI$Estimate_test_Type1, Eassoc_AllType_CI$Estimate_test_Type5, use='pairwise.complete.obs')
cor(Eassoc_AllType_CI$Estimate_test_Type4, Eassoc_AllType_CI$Estimate_test_Type5, use='pairwise.complete.obs')

cor(Eassoc_AllType_CI %>%
      select(all_of(c("Estimate_test_Type1",
                      "Estimate_test_Type2",
                      "Estimate_test_Type3",
                      "Estimate_test_Type4",
                      "Estimate_test_Type5",
                      "Estimate_test_Type6"))),use='pairwise.complete.obs')

cor(GxEassoc_AllType_CI %>%
      select(all_of(c("Estimate_test_Type1",
                      "Estimate_test_Type2",
                      "Estimate_test_Type3",
                      "Estimate_test_Type4",
                      "Estimate_test_Type5",
                      "Estimate_test_Type6"))),use='pairwise.complete.obs')


#Visualization w/ 5+ specifications outside of correlations:
Eassoc_all <- Eassoc_AllType_CI %>%
  select(all_of(c("AssocID", "Estimate_test_Type1",
                  "Estimate_test_Type2",
                  "Estimate_test_Type3",
                  "Estimate_test_Type4",
                  "Estimate_test_Type5",
                  "Estimate_test_Type6")))

# Reshape for plotting (with NAs intact)
long_data <- Eassoc_all %>%
  pivot_longer(cols = starts_with("Estimate"), names_to = "Type", values_to = "Value")


# Plot (lines skip missing values)
# Parallel coordinates plot
long_data2 <- na.omit(long_data)
long_data2 <- long_data2 %>%
                mutate(
                  Type = recode(Type,
                                "Estimate_test_Type1" = "Type1",
                                "Estimate_test_Type2" = "Type2",
                                "Estimate_test_Type3" = "Type3",
                                "Estimate_test_Type4" = "Type4",
                                "Estimate_test_Type5" = "Type5",
                                "Estimate_test_Type6" = "Type6"
                  )
                )


ggplot(long_data2, aes(x = Type, y = Value, group = AssocID)) +
  #geom_line(alpha = 0.1, color = "blue", na.rm = TRUE) +
  geom_line(alpha = 0.1, color = "blue") +
  geom_point() +
  theme_minimal() +
  labs(
    title = "Parallel Coordinates Plot: E Assocations",
    x = "Covariate Specifications",
    y = "Test Beta"
  )

long_data2 %>%
  filter(grepl("IGFBP2",AssocID)) %>%
  ggplot(aes(x = Type, y = Value, group = AssocID)) +
  #geom_line(alpha = 0.1, color = "blue", na.rm = TRUE) +
  geom_line(alpha = 0.1, color = "blue") +
  geom_point() +
  theme_minimal() +
  labs(
    title = "Parallel Coordinates Plot: E Assocations IGFBP2",
    x = "Covariate Specifications",
    y = "Test Beta"
  )

long_data2 %>%
  filter(grepl("LEP",AssocID)) %>%
  ggplot(aes(x = Type, y = Value, group = AssocID)) +
  #geom_line(alpha = 0.1, color = "blue", na.rm = TRUE) +
  geom_line(alpha = 0.1, color = "blue") +
  geom_point() +
  theme_minimal() +
  labs(
    title = "Parallel Coordinates Plot (Gaps for NAs)",
    x = "Data Types",
    y = "Values"
  )

long_data2 %>%
  filter(grepl("GDF15",AssocID)) %>%
  ggplot(aes(x = Type, y = Value, group = AssocID)) +
  #geom_line(alpha = 0.1, color = "blue", na.rm = TRUE) +
  geom_line(alpha = 0.1, color = "blue") +
  geom_point() +
  theme_minimal() +
  labs(
    title = "Parallel Coordinates Plot (Gaps for NAs)",
    x = "Data Types",
    y = "Values"
  )

long_data2 %>%
  filter(grepl("CKB",AssocID)) %>%
  ggplot(aes(x = Type, y = Value, group = AssocID)) +
  #geom_line(alpha = 0.1, color = "blue", na.rm = TRUE) +
  geom_line(alpha = 0.1, color = "blue") +
  geom_point() +
  theme_minimal() +
  labs(
    title = "Parallel Coordinates Plot (Gaps for NAs)",
    x = "Data Types",
    y = "Values"
  )


#### GxE Associations:
GxEassoc_all <- GxEassoc_AllType_CI %>%
  select(all_of(c("AssocID", "Estimate_test_Type1",
                  "Estimate_test_Type2",
                  "Estimate_test_Type3",
                  "Estimate_test_Type4",
                  "Estimate_test_Type5",
                  "Estimate_test_Type6")))

# Reshape for plotting (with NAs intact)
long_data <- GxEassoc_all %>%
  pivot_longer(cols = starts_with("Estimate"), names_to = "Type", values_to = "Value")


# Plot (lines skip missing values)
# Parallel coordinates plot
long_data2 <- na.omit(long_data)
long_data2 <- long_data2 %>%
  mutate(
    Type = recode(Type,
                  "Estimate_test_Type1" = "Type1",
                  "Estimate_test_Type2" = "Type2",
                  "Estimate_test_Type3" = "Type3",
                  "Estimate_test_Type4" = "Type4",
                  "Estimate_test_Type5" = "Type5",
                  "Estimate_test_Type6" = "Type6"
    )
  )


ggplot(long_data2, aes(x = Type, y = Value, group = AssocID)) +
  #geom_line(alpha = 0.1, color = "blue", na.rm = TRUE) +
  geom_line(alpha = 0.1, color = "blue") +
  geom_point() +
  theme_minimal() +
  labs(
    title = "Parallel Coordinates Plot: GxE Assocations",
    x = "Covariate Specifications",
    y = "Test Beta"
  )


#Visualization w/ 5+ specifications outside of correlations:
#'*Color by Cateogry:*
Eassoc_AllType_CIv2 <- Eassoc_AllType %>%
  select(AssocID, CovarType, Estimate_test, Std.Error_test, CImargin_test,Category_test) %>% 
  pivot_wider(names_from = CovarType, values_from = c(Estimate_test, Std.Error_test, CImargin_test))

Eassoc_all <- Eassoc_AllType_CIv2 %>%
  select(all_of(c("AssocID", "Category_test", "Estimate_test_Type1",
                  "Estimate_test_Type2",
                  "Estimate_test_Type3",
                  "Estimate_test_Type4",
                  "Estimate_test_Type5",
                  "Estimate_test_Type6")))

# Reshape for plotting (with NAs intact)
long_data <- Eassoc_all %>%
  pivot_longer(cols = starts_with("Estimate"), names_to = "Type", values_to = "Value")


# Plot (lines skip missing values)
# Parallel coordinates plot
long_data2 <- na.omit(long_data)
long_data2 <- long_data2 %>%
  mutate(
    Type = recode(Type,
                  "Estimate_test_Type1" = "Type1",
                  "Estimate_test_Type2" = "Type2",
                  "Estimate_test_Type3" = "Type3",
                  "Estimate_test_Type4" = "Type4",
                  "Estimate_test_Type5" = "Type5",
                  "Estimate_test_Type6" = "Type6"
    )
  )


ggplot(long_data2, aes(x = Type, y = Value, group = AssocID,
                       color = Category_test)) +
  #geom_line(alpha = 0.1, color = "blue", na.rm = TRUE) +
  geom_line(alpha = 0.5) +
  #geom_point() +
  theme_minimal() +
  labs(
    title = "Parallel Coordinates Plot: E Assocations",
    x = "Covariate Specifications",
    y = "Test Beta"
  )


## STATS to have
#Number of Replicated Findings across Specification
#Things Consistent Significance across Specification (Proportion)
#Things NOT Consistent across Specification (Proportion)
#Reminder - Coefficient of Interactions: POSITIVE IS SYNERGISTIC HERE! NEGATIVE IS COUNTERACTIVE!







#'*DEALING WITH TRENDS - IDENTIFYING TRENDS OF SIGNFICANT and REPLICATED ASSOCIATIONS!*
EassocCat_mean_abs <- Eassoc_AllType %>%
  group_by(Category_test, CovarType) %>%
  summarize(
    mean_effect_size = mean(abs(Estimate_test)),#weighted.mean(Estimate, w = 1 / `Std. Error`^2),
    se_mean_effect_size = sqrt(sum(`Std.Error_test`^2) / n()),
    n_observations = n(),
    .groups = 'drop'
  )

EassocCat_mean_abs %>%
  ggplot(aes(x=reorder(Category_test, -mean_effect_size), 
             y=mean_effect_size, 
             color = CovarType)) +
    geom_jitter(width = 0.3, height = 0.001) +
    ylim(0,0.4) +
    labs(
      x = NULL,
      y = "Absolute Mean Effect Size"
    ) +
    theme_minimal() +
    theme(axis.text.x=element_text(angle=60, hjust=1))



###GET Meta-Analysis Point Estimates of Effect Size for each Category!
library(metafor)

# Perform random-effects meta-analysis for each unique Category_test
TrendsRes <- EassocCat_mean_abs %>%
  group_by(Category_test) %>%
  summarise(
    Beta = as.numeric(rma(yi = mean_effect_size, 
                          sei = se_mean_effect_size, 
                          method = "REML", 
                          data = pick(everything()))$beta),
    SE = rma(yi = mean_effect_size, 
             sei = se_mean_effect_size, 
             method = "REML", 
             data = pick(everything()))$se,
    CI_lower =  rma(yi = mean_effect_size, 
                    sei = se_mean_effect_size, 
                    method = "REML", 
                    data = pick(everything()))$ci.lb,
    CI_upper = rma(yi = mean_effect_size, 
                   sei = se_mean_effect_size, 
                   method = "REML", 
                   data = pick(everything()))$ci.ub,
    .groups = 'drop'
  ) %>%
  rename(ID = Category_test)

# PLOT Averaged Main E Association Effects:
#qt(0.975, df = 5) - DO not do CIs here --> Not sure if VALID...

TrendsRes %>%
  ggplot(aes(x = reorder(ID, -Beta), y = Beta, color = ID)) +
  geom_point() +
  geom_linerange(aes(ymin = CI_lower, 
                ymax = CI_upper)) +
  labs(x = NULL,
    y = "Absolute Mean Effect Size") +
  theme_minimal() +
  theme(axis.text.x=element_text(angle=60, hjust=1),
        legend.position="none")


#'*TRENDS for certain E features*
# REMEMBER THIS IS FOR THE SET OF ALL SIGNIFICANT ASSOCIATIONS!!!
Eassoc_mean_abs <- Eassoc_AllType  %>%
  group_by(ID, CovarType) %>%
  summarize(
    mean_effect_size = mean(abs(Estimate_test)),#weighted.mean(Estimate, w = 1 / `Std. Error`^2),
    se_mean_effect_size = sqrt(sum(`Std.Error_test`^2) / n()),
    n_observations = n(),
    .groups = 'drop'
  )

E_TrendsRes <- Eassoc_mean_abs %>%
  filter(grepl("alcohol_intake_frequency",ID)) %>%
  group_by(ID) %>%
  summarise(
    Beta = as.numeric(rma(yi = mean_effect_size, 
                          sei = se_mean_effect_size, 
                          method = "REML", 
                          data = pick(everything()))$beta),
    SE = rma(yi = mean_effect_size, 
             sei = se_mean_effect_size, 
             method = "REML", 
             data = pick(everything()))$se,
    CI_lower =  rma(yi = mean_effect_size, 
                    sei = se_mean_effect_size, 
                    method = "REML", 
                    data = pick(everything()))$ci.lb,
    CI_upper = rma(yi = mean_effect_size, 
                   sei = se_mean_effect_size, 
                   method = "REML", 
                   data = pick(everything()))$ci.ub,
    .groups = 'drop'
  )

E_TrendsRes %>%
  ggplot(aes(x = reorder(ID, Beta), y = Beta)) +
  geom_point() +
  geom_linerange(aes(ymin = CI_lower, 
                     ymax = CI_upper)) +
  labs(x = "Alcohol Intake",
       y = "Absolute Mean Effect Size") +
  scale_x_discrete(labels=c("alcohol_intake_frequency_f1558_0_02" = "Special Occasions", 
                            "alcohol_intake_frequency_f1558_0_03" = "1-3 Times Month",
                            "alcohol_intake_frequency_f1558_0_04" = "1-2 Times Week",
                            "alcohol_intake_frequency_f1558_0_05" = "3-4 Times Week",
                            "alcohol_intake_frequency_f1558_0_06" = "Daily/Almost Daily")) +
  theme_minimal() +
  theme(axis.text.x=element_text(angle=60, hjust=1),
        legend.position="none")


#ALL other Variables:



# TrendsRes <- data.frame(
#   ID = c(),
#   Beta = c(),
#   SE = c()
# )
# c = 1
# for(i in unique(EassocCat_mean_abs$Category_test)){
#   data <- EassocCat_mean_abs %>%
#             filter(Category_test == i)
#   # Run a random-effects meta-analysis
#   meta_result <- rma(yi = mean_effect_size, sei = se_mean_effect_size, method = "REML", data = data)
# 
#   TrendsRes[c, "ID"] <- i
#   TrendsRes[c, "Beta"] <- as.numeric(meta_result$beta)
#   TrendsRes[c, "SE"] <- meta_result$se
# 
#   c = c + 1
# }






#NOT USEFUL SINCE ANY applies to any combo:
# Eassoc_Types1_6 <- Eassoc_AllType %>%
#                       group_by(AssocID) %>%
#                       filter(any(c("Type1","Type2","Type3",
#                                    "Type4","Type5","Type6") %in% CovarType)) %>% #any is used since we want to find connections across multiple specifications
#                       ungroup()
# 
# GxEassoc_Types1_6 <- GxEassoc_AllType %>%
#                         group_by(AssocID) %>%
#                         filter(any(c("Type1","Type2","Type3",
#                                      "Type4","Type5","Type6") %in% CovarType)) %>% #any is used since we want to find connections across multiple specifications
#                         ungroup()


# First filter associations detected in both Specifications:
Eassoc_Type1_5 <- Eassoc_AllType %>%
  group_by(AssocID) %>%
  filter(any(c("Type1", "Type5") %in% CovarType)) %>%
  #filter(all(c("Type1", "Type5") %in% CovarType)) %>%
  ungroup()





###PULL from discovery set:

#NEED Structure w/ all Associations and correspondign covariates:
CovarSpec1_Euni <- replication_stats_final(Type1$train[[1]], Type1$test[[1]])
CovarSpec1_GxEuni <- replication_stats_final(Type1$train[[2]], Type1$test[[2]])

CovarSpec5_Euni <- replication_stats_final(Type5$train[[1]], Type5$test[[1]])
CovarSpec5_GxEuni <- replication_stats_final(Type5$train[[2]], Type5$test[[2]])

#Question in multiverse analysis do you correct for # of specifications?
CovarSpec1_Edisc <- CovarSpec1_Euni$sigBOTH
CovarSpec5_Edisc <- CovarSpec5_Euni$sigBOTH

CovarSpec1_Edisc$CovarType <- "Type1"
CovarSpec5_Edisc$CovarType <- "Type5"

Eassoc_AllType <- list(CovarSpec1_Edisc, CovarSpec5_Edisc) %>% reduce(full_join)


#Correlations:
Eassoc_AllType$AssocID <- paste0(Eassoc_AllType$omicID,":",Eassoc_AllType$ID)

# First filter associations detected in both Specifications:
Eassoc_Type1_5 <- Eassoc_AllType %>%
  group_by(AssocID) %>%
  filter(any(c("Type1", "Type5") %in% CovarType)) %>%
  #filter(all(c("Type1", "Type5") %in% CovarType)) %>%
  ungroup()

# Pivot the data to wide format: Both Train and Test Set Estimates!!
Type1_5_train <- Eassoc_Type1_5 %>%
  select(AssocID, CovarType, Estimate_train) %>% 
  pivot_wider(names_from = CovarType, values_from = Estimate_train)


Type1_5_test <- Eassoc_Type1_5 %>%
  select(AssocID, CovarType, Estimate_test) %>% 
  pivot_wider(names_from = CovarType, values_from = Estimate_test)


#Get standard errors of each estimate:
Type1_5_test_se <- Eassoc_Type1_5 %>%
                      select(AssocID, CovarType, Estimate_test, Std.Error_test) %>% 
                      pivot_wider(names_from = CovarType, values_from = c(Estimate_test, Std.Error_test))

#'*Make confidence interval in Essoc_Type1_5*
#'Reminder N is differing for each association:

#qt function is vectorized for differing size n!
critical_values <- qt(0.975, df = Eassoc_Type1_5$samplesize_test - 1)


#95% confidence interval:
Eassoc_Type1_5$CImargin_test <- Eassoc_Type1_5$Std.Error_test * qt(0.975, df = Eassoc_Type1_5$samplesize_test - 1)



Type1_5_test_CI <- Eassoc_Type1_5 %>%
  select(AssocID, CovarType, Estimate_test, Std.Error_test, CImargin_test) %>% 
  pivot_wider(names_from = CovarType, values_from = c(Estimate_test, Std.Error_test, CImargin_test))



###PLOTS:
Type1_5_train %>%
  ggplot(aes(x = Type1, y = Type5)) +
    geom_point() +
    geom_smooth(method = "lm") +
    theme_minimal()


Type1_5_test %>%
  ggplot(aes(x = Type1, y = Type5)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_minimal()

#With Confidence Intervals:
Type1_5_test_CI %>%
  ggplot(aes(x = Estimate_test_Type1, y = Estimate_test_Type5)) +
  geom_point(size = 1) +
  geom_linerange(aes(ymin = Estimate_test_Type5 - CImargin_test_Type5, 
                     ymax = Estimate_test_Type5 + CImargin_test_Type5),
                     width = 0.2) +  # Vertical error bar
  geom_linerange(aes(xmin = Estimate_test_Type1 - CImargin_test_Type1, 
                     xmax = Estimate_test_Type1 + CImargin_test_Type1), 
                     height = 0.2) + # Horizontal error bar
  theme_minimal()


#CHANGE in sign (Type 1 vs Type 5):
Type1_5_test_CI %>%
  filter(!(sign(Estimate_test_Type1) == sign(Estimate_test_Type5)))
sign(-1)
sign(-1.2)


##MAKE Datastructure FOR ABOVE FOR EASIER PLOTTING!!!
  


#With Standard Errors:





test_assoc1 <- Type1$test[[1]]
test_assoc2 <- Type5$test[[1]]

colnames(test_assoc1)[!(colnames(test_assoc1) %in% c("ID","omicID"))] <-  gsub(" ", "",
                                                                               paste0(colnames(test_assoc1)[!(colnames(test_assoc1) %in% c("ID","omicID"))], "_assoc1"))

colnames(test_assoc2)[!(colnames(test_assoc2) %in% c("ID","omicID"))] <-  gsub(" ", "",
                                                                               paste0(colnames(test_assoc2)[!(colnames(test_assoc2) %in% c("ID","omicID"))], "_assoc2"))

all_assoc <- list(test_assoc1, test_assoc2) %>% reduce(full_join)

pval_thresh <- 0.05/(188 * 2686)

all_assoc_significant <- all_assoc %>%
  filter(`Pr(>|t|)_assoc1` < pval_thresh &
           `Pr(>|t|)_assoc2` < pval_thresh)

#Example Plot: Are Estimates similar between train & test?
all_assoc_significant %>%
  ggplot(aes(x = Estimate_assoc1, y = Estimate_assoc2)) +
  geom_point() +
  geom_smooth(method = "lm")






#Sign Flip Example:
miamiplot_single_omic("FABP4", Type1$train[[1]], "FABP4 E Associations")
miamiplot_single_omic("FABP4", Type5$train[[1]], "FABP4 E Associations")

miamiplot_single_omic("IGFBP2", Type1$train[[1]], "IGFBP2 E Associations")

#Train vs Test Pvalues and Betas.
miamiplot_single_omic("IGFBP2", Type5$train[[1]], "IGFBP2 E Associations")
miamiplot_single_omic("IGFBP2", Type5$test[[1]], "IGFBP2 E Associations")


miamiplot_single_omic("GDF15", Type5$train[[1]], "GDF15 E Associations Train Set")
miamiplot_single_omic("GDF15", Type5$test[[1]], "GDF15 E Associations Test Set")

miamiplot_single_omic("FABP4", Type5$train[[1]], "FABP4 E Associations Train Set")
miamiplot_single_omic("FABP4", Type5$test[[1]], "FABP4 E Associations Test Set")




miamiplot_single_omic("LEP", Type5$train[[1]], "LEP E Associations")
miamiplot_single_omic("APOE", Type5$train[[1]], "APOE E Associations")

miamiplot_single_omic("GFAP", Type1$train[[1]], "GFAP E Associations")
miamiplot_single_omic("GFAP", Type1$test[[1]], "GFAP E Associations")

miamiplot_single_omic("GFAP", Type5$train[[1]], "GFAP E Associations")
miamiplot_single_omic("GFAP", Type5$test[[1]], "GFAP E Associations")

miamiplot_single_omic("NEFL", Type5$train[[1]], "NEFL E Associations")



# Example usage
train <- Type7$train
test <- Type7$test
covarName <- "Covarspec7"

train <- Type6$train
test <- Type6$test
covarName <- "Covarspec6"


# Replication Plot
Eassoc_Type5 <- replication_stats(train[[1]], 
                                  test[[1]],
                                  title = paste0(covarName," E Associations"))
GxEassoc_Type5 <- replication_stats(train[[2]], 
                                  test[[2]],
                                  title = paste0(covarName," GxE Associations"))

#'*Adjust this*

# Miami Plot
miamiplot_multi_omic(train[[1]], paste0(covarName," (Train) - Entire Proteome: E Associations"))
miamiplot_multi_omic(train[[2]], paste0(covarName," (Train) - Entire Proteome: GxE Associations"))


miamiplot_multi_omic(test[[1]], paste0(covarName," (Test) - Entire Proteome: E Associations"))
miamiplot_multi_omic(test[[2]], paste0(covarName," (Test) - Entire Proteome: GxE Associations"))
#What to plot? Train vs Test.
#IDEA plot the x-axis by how correlated the proteins are (1.2.3,...) for each category
#Or make sure the order is consist for proteins across category.



#'*To-DO:*
#'Visualizations: Each Category should have a specific color code!!! Helps to make everything match between covarSpec.


### Function to Find Significant Hits shared across Specifications.

# Check replication:
replication_stats_tbl <- function(train, test){
  train_assoc <- train
  test_assoc <- test
  
  colnames(train_assoc)[!(colnames(train_assoc) %in% c("ID","omicID"))] <-  gsub(" ", "",
                                                                                 paste0(colnames(train_assoc)[!(colnames(train_assoc) %in% c("ID","omicID"))], "_train"))
  
  colnames(test_assoc)[!(colnames(test_assoc) %in% c("ID","omicID"))] <-  gsub(" ", "",
                                                                               paste0(colnames(test_assoc)[!(colnames(test_assoc) %in% c("ID","omicID"))], "_test"))
  
  
  all_assoc <- list(train_assoc, test_assoc) %>% reduce(full_join)
  
  pval_thresh <- 0.05/nrow(all_assoc)
  
  all_assoc_significant <- all_assoc %>%
    filter(`Pr(>|t|)_test` < pval_thresh &
             `Pr(>|t|)_train` < pval_thresh)
  
  return(all_assoc_significant)
}

# Replication Plot
Eassoc_Type1 <- replication_stats_tbl(Type1$train[[1]], 
                                      Type1$test[[1]])
Eassoc_Type2 <- replication_stats_tbl(Type2$train[[1]], 
                                      Type2$test[[1]])
Eassoc_Type3 <- replication_stats_tbl(Type3$train[[1]], 
                                      Type3$test[[1]])
Eassoc_Type4 <- replication_stats_tbl(Type4$train[[1]], 
                                      Type4$test[[1]])
Eassoc_Type5 <- replication_stats_tbl(Type5$train[[1]], 
                                      Type5$test[[1]])
Eassoc_Type6 <- replication_stats_tbl(Type6$train[[1]], 
                                      Type6$test[[1]])
Eassoc_Type7 <- replication_stats_tbl(Type7$train[[1]], 
                                  Type7$test[[1]])


#Combine effects across Type 1 - Type 5 first:
Eassoc_Type1$CovarType <- "Type1"
Eassoc_Type2$CovarType <- "Type2"
Eassoc_Type3$CovarType <- "Type3"
Eassoc_Type4$CovarType <- "Type4"
Eassoc_Type5$CovarType <- "Type5"
Eassoc_Type6$CovarType <- "Type6"
Eassoc_Type7$CovarType <- "Type7"

Eassoc_AllType <- list(Eassoc_Type1, Eassoc_Type2, Eassoc_Type3, Eassoc_Type4,
     Eassoc_Type5, Eassoc_Type6, Eassoc_Type7) %>% reduce(full_join)

#Look at
CovarSpecCompare <- Eassoc_AllType %>%
                        filter(CovarType %in% c("Type1","Type2","Type3","Type4","Type5"))



#Make new dataframe for CovarSpec Analysis:
CovarSpecCompare$ExID <- paste0(CovarSpecCompare$ID,"_",CovarSpecCompare$omicID)

CovarSpecSign <- CovarSpecCompare %>%
  select(all_of(c("ExID","CovarType","Estimate_train"))) %>%
  pivot_wider(
    names_from = CovarType,
    values_from = Estimate_train
  )

CovarSpecSign$SignAgreement <- apply(CovarSpecSign[,-1],  1, function(row) {
  if (all(row > 0, na.rm = TRUE)) {
    return(1)  # All positive
  } else if (all(row < 0, na.rm = TRUE)) {
    return(1)  # All negative
  } else {
    return(0)  # Mixed signs or NAs
  }
})


v1 <- CovarSpecCompare %>% filter(CovarType == "Type1")
v2 <- CovarSpecCompare %>% filter(CovarType == "Type2")


CovarSpecCompare_1 

#Simpsons Paradox: https://en.wikipedia.org/wiki/Simpson%27s_paradox

#Count number of times sign flip occurs across the associations Type1 --> Type5.

###With Conditional on Income/Deprivation:
#Look at
CovarSpecCompare <- Eassoc_AllType %>%
  filter(CovarType %in% c("Type1","Type2","Type3","Type4","Type5","Type6"))



#Make new dataframe for CovarSpec Analysis:
CovarSpecCompare$ExID <- paste0(CovarSpecCompare$ID,"_",CovarSpecCompare$omicID)

CovarSpecSign <- CovarSpecCompare %>%
  select(all_of(c("ExID","CovarType","Estimate_train"))) %>%
  pivot_wider(
    names_from = CovarType,
    values_from = Estimate_train
  )

CovarSpecSign$SignAgreement <- apply(CovarSpecSign[,-1],  1, function(row) {
  if (all(row > 0, na.rm = TRUE)) {
    return(1)  # All positive
  } else if (all(row < 0, na.rm = TRUE)) {
    return(1)  # All negative
  } else {
    return(0)  # Mixed signs or NAs
  }
})


### Make distance matrix:
CovarSpecDist <- as.matrix(dist(t(CovarSpecSign %>%
                                    select(-c("ExID","SignAgreement")))))

mds_result <- cmdscale(CovarSpecDist, k = 2, eig = T)

variance_explained <- mds_result$eig / sum(mds_result$eig)

lookup <- c(CovarSpec = "ID", MDS1 = "V1", MDS2 = "V2")

p1 <- mds_result$points %>%
  as.data.frame() %>%
  rownames_to_column(var = "ID") %>%
  rename(all_of(lookup)) %>%
  ggplot(aes(x = MDS1, y = MDS2, color = CovarSpec)) +
  geom_point() +
  xlab(paste0("PC1 ","(", round(variance_explained[1], 2) * 100,"%)")) +
  ylab(paste0("PC2 ","(", round(variance_explained[2], 2) * 100,"%)")) +
  labs(color = "Covar Spec") +
  #geom_text_repel() +
  theme_minimal()
p1


###Function to visualize R2s across proteome for individual exposures:
train_assoc1 <- train[[3]]
test_assoc1 <- test[[3]]

colnames(train_assoc1)[!(colnames(train_assoc1) %in% c("ID","omicID"))] <-  gsub(" ", "",
                                                                                 paste0(colnames(train_assoc1)[!(colnames(train_assoc1) %in% c("ID","omicID"))], "_train"))

colnames(test_assoc1)[!(colnames(test_assoc1) %in% c("ID","omicID"))] <-  gsub(" ", "",
                                                                               paste0(colnames(test_assoc1)[!(colnames(test_assoc1) %in% c("ID","omicID"))], "_test"))


all_assoc <- list(train_assoc1, test_assoc1) %>% reduce(full_join)


#Are R2 similar across train/test
all_assoc %>%
  ggplot(aes(x = R2_train, y = R2_test)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_minimal()


GxE <- all_assoc %>%
  filter(grepl("_GS", ID))
GxE %>%
  ggplot(aes(x = R2_train, y = R2_test)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_minimal()



E <- all_assoc %>%
  filter(!grepl("_GS", ID))
E %>%
  ggplot(aes(x = R2_train, y = R2_test)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_minimal()


#'*TO-DO: Make plot showing correlations among differing covariate specifications*
#'*PLOT Correlations*

#Covariate Specification: SIGNIFICANT ESTIMATES CORRELATIONS!!!
##IDEA ADD ERROR BARS TO THE POINT ESTIMATES
## USE DEMING REGRESSION OR LM TO GET SLOPE, CORRELATION etc.

test_assoc1 <- Type1$test[[1]]
test_assoc2 <- Type5$test[[1]]

colnames(test_assoc1)[!(colnames(test_assoc1) %in% c("ID","omicID"))] <-  gsub(" ", "",
                                                                                 paste0(colnames(test_assoc1)[!(colnames(test_assoc1) %in% c("ID","omicID"))], "_assoc1"))

colnames(test_assoc2)[!(colnames(test_assoc2) %in% c("ID","omicID"))] <-  gsub(" ", "",
                                                                               paste0(colnames(test_assoc2)[!(colnames(test_assoc2) %in% c("ID","omicID"))], "_assoc2"))

all_assoc <- list(test_assoc1, test_assoc2) %>% reduce(full_join)

pval_thresh <- 0.05/(188 * 2686)

all_assoc_significant <- all_assoc %>%
  filter(`Pr(>|t|)_assoc1` < pval_thresh &
           `Pr(>|t|)_assoc2` < pval_thresh)

#Example Plot: Are Estimates similar between train & test?
all_assoc_significant %>%
  ggplot(aes(x = Estimate_assoc1, y = Estimate_assoc2)) +
  geom_point() +
  geom_smooth(method = "lm")



#Train/Test Correlations FOR EACH SPECIFICATION
train_assoc1 <- Type5$train[[1]]
test_assoc1 <- Type5$test[[1]]

colnames(train_assoc1)[!(colnames(train_assoc1) %in% c("ID","omicID"))] <-  gsub(" ", "",
                                                                                 paste0(colnames(train_assoc1)[!(colnames(train_assoc1) %in% c("ID","omicID"))], "_train"))

colnames(test_assoc1)[!(colnames(test_assoc1) %in% c("ID","omicID"))] <-  gsub(" ", "",
                                                                               paste0(colnames(test_assoc1)[!(colnames(test_assoc1) %in% c("ID","omicID"))], "_test"))


all_assoc <- list(train_assoc1, test_assoc1) %>% reduce(full_join)

pval_thresh <- 0.05/(188 * 2686)

##'*TO TRY - Significant Associations in Discovery but not Replication Set and Vice-Verse:*
all_assoc_significant <- all_assoc %>%
  filter(`Pr(>|t|)_test` < pval_thresh &
           `Pr(>|t|)_train` < pval_thresh)

#Example Plot: Are Estimates similar between train & test?
all_assoc_significant %>%
  ggplot(aes(x = Estimate_train, y = Estimate_test)) +
  geom_point() +
  geom_smooth(method = "lm")



#States to Have:
#Number replicated in train/test split
#Regression of Betas in Train and test.

# Type1_assoc <- loadUniAssoc(covarType = "Type1")
# Type2_assoc <- loadUniAssoc(covarType = "Type2")
# Type3_assoc <- loadUniAssoc(covarType = "Type3")
# Type4_assoc <- loadUniAssoc(covarType = "Type4")
# Type6_assoc <- loadUniAssoc(covarType = "Type6")
# Type7_assoc <- loadUniAssoc(covarType = "Type7")



# Join train and test statistics with `purrr::map2`
# stat_train <<- suppressMessages(
#   map2(stat_train %||% vector("list", length(load$train)), 
#        load$train, 
#        ~ .x %||% .y %>% full_join(.y, by = NULL))
# )
# 
# stat_test <<- suppressMessages(
#   map2(stat_test %||% vector("list", length(load$test)), 
#        load$test, 
#        ~ .x %||% .y %>% full_join(.y, by = NULL))
# )

#First Thing: Miami Plot of E and GxE for a Specification








#Second Thing: Replication across train and test set.


# loadUniAssoc <- function(covarType){
#   stat_all <- list()
#   stat_train <- NULL
#   stat_test <- NULL
#   
#   for(i in 1:400){
#     tryCatch({
#       load <- readRDS(paste0("/n/groups/patel/shakson_ukb/UK_Biobank/Data/Parallel/Prot_PGSPXS_univar/",covarType,"/univar_assoc_",i,".rds"))
#       
#       stat_train <- suppressMessages(map2(stat_train %||% vector("list", length(load$train)), 
#                                           load$train, ~ {
#                                             .x %||% .y %>% full_join(.y, by = NULL)
#                                           }))
#       stat_test <- suppressMessages(map2(stat_test %||% vector("list", length(load$test)), 
#                                          load$test, ~ {
#                                            .x %||% .y %>% full_join(.y, by = NULL)
#                                          }))
#       
#       
#     }, error = function(e) {
#       # Handle the error (e.g., print an error message)
#       message(paste("Error reading file:",i))
#       error_files <<- append(error_files, i)
#     })
#   }
#   
#   stat_all[["train"]] <- stat_batch_train
#   stat_all[["test"]] <- stat_batch_test
#   
#   return(stat_all)
# }
# Type1_assoc <- loadUniAssoc(covarType = "Type1")



#### Function 6: Establish if Replication Occurs between Train and Test

#Using example batch see if 
replication_run <- example_run

train_assoc1 <- replication_run$train[[1]]
test_assoc1 <- replication_run$test[[1]]


colnames(train_assoc1)[!(colnames(train_assoc1) %in% c("ID","omicID"))] <-  gsub(" ", "",
                                                                                 paste0(colnames(train_assoc1)[!(colnames(train_assoc1) %in% c("ID","omicID"))], "_train"))

colnames(test_assoc1)[!(colnames(test_assoc1) %in% c("ID","omicID"))] <-  gsub(" ", "",
                                                                               paste0(colnames(test_assoc1)[!(colnames(test_assoc1) %in% c("ID","omicID"))], "_test"))


all_assoc <- list(train_assoc1, test_assoc1) %>% reduce(full_join)

pval_thresh <- 0.05/(188 * 2686)

all_assoc_significant <- all_assoc %>%
  filter(`Pr(>|t|)_test` < pval_thresh &
           `Pr(>|t|)_train` < pval_thresh)

#Example Plot: Are Estimates similar between train & test?
all_assoc_significant %>%
  ggplot(aes(x = Estimate_train, y = Estimate_test)) +
  geom_point() +
  geom_smooth(method = "lm")


####Ver2
#Using example batch see if 
replication_run <- example_run

train_assoc1 <- replication_run$train[[2]]
test_assoc1 <- replication_run$test[[2]]


colnames(train_assoc1)[!(colnames(train_assoc1) %in% c("ID","omicID"))] <-  gsub(" ", "",
                                                                                 paste0(colnames(train_assoc1)[!(colnames(train_assoc1) %in% c("ID","omicID"))], "_train"))

colnames(test_assoc1)[!(colnames(test_assoc1) %in% c("ID","omicID"))] <-  gsub(" ", "",
                                                                               paste0(colnames(test_assoc1)[!(colnames(test_assoc1) %in% c("ID","omicID"))], "_test"))


all_assoc <- list(train_assoc1, test_assoc1) %>% reduce(full_join)

pval_thresh <- 0.05/(188 * 2686)

all_assoc_significant <- all_assoc %>%
  filter(`Pr(>|t|)_test` < pval_thresh &
           `Pr(>|t|)_train` < pval_thresh)

#Example Plot: Are Estimates similar between train & test?
all_assoc_significant %>%
  ggplot(aes(x = Estimate_train, y = Estimate_test)) +
  geom_point() +
  geom_smooth(method = "lm")


####Ver3
#Using example batch see if 
replication_run <- example_run

train_assoc1 <- replication_run$train[[3]]
test_assoc1 <- replication_run$test[[3]]


colnames(train_assoc1)[!(colnames(train_assoc1) %in% c("ID","omicID"))] <-  gsub(" ", "",
                                                                                 paste0(colnames(train_assoc1)[!(colnames(train_assoc1) %in% c("ID","omicID"))], "_train"))

colnames(test_assoc1)[!(colnames(test_assoc1) %in% c("ID","omicID"))] <-  gsub(" ", "",
                                                                               paste0(colnames(test_assoc1)[!(colnames(test_assoc1) %in% c("ID","omicID"))], "_test"))


all_assoc <- list(train_assoc1, test_assoc1) %>% reduce(full_join)



pval_thresh <- 0.05/(188 * 2686)

# all_assoc_significant <- all_assoc %>%
#   filter(`Pr(>|t|)_test` < pval_thresh &
#            `Pr(>|t|)_train` < pval_thresh)

#Example Plot: Are Estimates similar between train & test?
all_assoc %>%
  ggplot(aes(x = R2_train, y = R2_test)) +
  geom_point() +
  geom_smooth(method = "lm")

#plot(all_assoc_significant$Estimate_train, all_assoc_significant$Estimate_test)



#Set up function for varying p-values
#Try different methods:

train_assoc1 <- trial_stat_train[[1]]
test_assoc1 <- trial_stat_test[[1]]


view(trial_stat_test[[3]])
view(trial_stat_test[[4]])

train_assoc1 <- trial_stat_train[[2]]
test_assoc1 <- trial_stat_test[[2]]

##Look for significant hits from both train and test
##See if effects are correlated well

colnames(train_assoc1)[-1] <- gsub(" ", "",
                                   paste0(colnames(train_assoc1)[-1], "_train"))
colnames(test_assoc1)[-1] <- gsub(" ", "",
                                  paste0(colnames(test_assoc1)[-1], "_test"))

all_assoc <- merge(train_assoc1, test_assoc1, by = "ID")

pval_thresh <- 0.05/(188 * 2686)

all_assoc_significant <- all_assoc %>%
  filter(`Pr(>|t|)_test` < pval_thresh &
           `Pr(>|t|)_train` < pval_thresh)

plot(all_assoc_significant$Estimate_train, 
     all_assoc_significant$Estimate_test)


###### Visualization and Handling Outputs to Record #######


# mutate(label = ifelse(is.na(label), 
#                       label, 
#                       str_wrap(label, width = 10))) 

library(ggrepel)
R2univar_plot <- function(df,omic){
  R2omic <- df %>% filter(omicID == omic)
  R2omic$idx <- 1:nrow(R2omic)
  
  
  p1 <- R2omic %>%
    mutate(label = ifelse(rank(-R2) <= 5, ID, NA)) %>%
    mutate(label = gsub("f[0-9]+_0_0", "", label)) %>%
    mutate(label = gsub("combined", "", label)) %>%
    mutate(label = gsub("_"," ", label)) %>%
    mutate(label = str_wrap(label, width = 40)) %>%
    ggplot(aes(x = reorder(idx, -R2), y = R2)) +
    geom_point() +
    geom_text_repel(aes(label = label), 
                    na.rm = TRUE,
                    box.padding = 0.3,
                    size = 3) +
    theme_classic() +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) +
    ggtitle(paste0(omic," univariate R2"))
  return(p1)
}

p1 <- R2univar_plot(topOmicHits$R2, "GDF15")
p2 <- R2univar_plot(topOmicHits$R2, "LEP")
p3 <- R2univar_plot(topOmicHits$R2, "CA14")
p4 <- R2univar_plot(topOmicHits$R2, "FABP4")
p5 <- R2univar_plot(topOmicHits$R2, "IL1RN")
p6 <- R2univar_plot(topOmicHits$R2, "HGF")

library(patchwork)
(p1 | p2 | p3) / (p4 | p5 | p6)

e1 <- R2univar_plot(topOmicHits$R2, "CXCL17")
e2 <- R2univar_plot(topOmicHits$R2, "ALPP")
e3 <- R2univar_plot(topOmicHits$R2, "MAMDC4")
e4 <- R2univar_plot(topOmicHits$R2, "SERPINA7")
(e1 | e2) / (e3 | e4)

#Extra for later:
R2univar_plot(topOmicHits$R2, "APOF")
R2univar_plot(topOmicHits$R2, "REN")
R2univar_plot(topOmicHits$R2, "FGF21")
R2univar_plot(topOmicHits$R2, "ADM")
R2univar_plot(topOmicHits$R2, "ELN")



#Simple plot: Function from "viz_PGS_PXS_univar.R"
library(plotly)
library(htmlwidgets)
miamiplot_single_omic <- function(protID, assoc, title){
  total <- assoc %>% filter(omicID == protID)
  
  order <- 1:nrow(total)
  p <- ggplot(total, aes(x = order, y = -log10(`Pr(>|t|)`) * sign(Estimate), color = category)) +
    geom_point(aes(text = paste("ID:", id, 
                                "<br>Beta:", signif(Estimate, digits = 2),
                                "<br>Pval:", signif(`Pr(>|t|)`, digits = 2),
                                "<br>N:", samplesize,
                                "<br>Category:", category,
                                "<br>Total R2:", signif(adj.R2, digits = 2)))) +
    geom_abline(intercept = -log10(0.05/nrow(total)), slope = 0, color = "blue", linetype = "dashed") +
    geom_abline(intercept = log10(0.05/nrow(total)), slope = 0, color = "blue", linetype = "dashed") +
    #facet_grid(rows = vars(CHR), scales = "free_x", space = "free_x") +
    #scale_y_continuous(trans = "reverse", limits = c(max(merged_df$logP1), 0)) +
    theme_minimal() +
    #theme(panel.spacing = unit(0.1, "lines")) +
    labs(x = "Exposures", y = "-log10(P-value) * sign(Beta)") +
    ggtitle(title)
  
  print(p)
  
  ggplotly(p, tooltip = "text")
}

#APOF example with GxE associations:
miamiplot_single_omic(protID = "GDF15", assoc = try1$E, title = "GDF15: E Associations w/o medication adj")
miamiplot_single_omic(protID = "GDF15", assoc = try2$E, title = "GDF15: E Associations w/ medication adj")
miamiplot_single_omic(protID = "GDF15", assoc = try2$GxE, title = "GDF15: GxE Associations w/ medication adj")


miamiplot_single_omic(protID = "LEP", assoc = try1$E, title = "LEP: E Associations w/o medication adj")
miamiplot_single_omic(protID = "LEP", assoc = try2$E, title = "LEP: E Associations w/ medication adj")



miamiplot_single_omic(protID = "FABP4", assoc = try1$E, title = "FABP4: E Associations w/o medication adj")
miamiplot_single_omic(protID = "FABP4", assoc = try2$E, title = "FABP4: E Associations w/ medication adj")

