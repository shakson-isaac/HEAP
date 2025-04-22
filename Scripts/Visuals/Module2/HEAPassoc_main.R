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
library(qs)

#'*LOAD CovarSpecList*
setwd("/n/groups/patel/shakson_ukb/UK_Biobank/")
HEAPassoc <- qread("./Output/HEAPres/HEAPassoc.qs")

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
         plot = p2, width = 12, height = 3, dpi = 1000)
}

# FUNCTION: Run through EVERY CovarSpec to Find Replicated Points
Replication_MiamiPlot <- function(CovarSpecList){
  lapply(names(CovarSpecList), function(x){
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

# Replication Miami Plots
#'*SWITCH CovarSpecAssocList --> CovarSpecList b/c it is MISLEADING HERE!!!!*
Replication_MiamiPlot(HEAPassoc@HEAPlist)
