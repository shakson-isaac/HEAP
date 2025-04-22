library(data.table)
library(tidyverse)
library(ggplot2)
library(ggpmisc)
library(ggrepel)
library(gridExtra)
library(plotly)
library(htmlwidgets)
library(pbapply)
library(dplyr)
library(purrr)

#'*LOAD CovarSpecList*
setwd("/n/groups/patel/shakson_ukb/UK_Biobank/")
HEAPassoc <- qread("./Output/HEAPres/HEAPassoc.qs")

CovarSpecList <- HEAPassoc@HEAPlist
CovarSpecAssocList <- HEAPassoc@HEAPsig

#Requires names to be "Type" not "Model": for Website
names(CovarSpecList) <- c("Type1", "Type2", "Type3", "Type4", "Type5", "Type6")
names(CovarSpecAssocList) <- c("Type1", "Type2", "Type3", "Type4", "Type5", "Type6")


#'*HTMLs of Univariate Associations*

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

# FUNCTION: Obtain interactive plot - where we highlight associations that replicate:
miamiplot_single_omic_rep <- function(protID, Eall, title){
  # GxEall <- replication_all(train = CovarSpecList[[CovarType]]$train[[2]],
  #                           test = CovarSpecList[[CovarType]]$test[[2]])
  
  #Arrange manhattan plot by E: Category
  assoc <- Eall %>%
    filter(omicID == protID) %>%
    arrange(Category_train)
  
  # Obtain significant results to plot first w/ COLOR!!
  pval_thresh <- 0.05/nrow(Eall) #Standard Bonferroni Correction for All ASSOCIATIONS. Ex. for 200 exposures, 3000 proteins: 0.05/(200*3000)
  
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
  
  #Specify colors for each category: Regardless of Specification!!
  category_colors <- c("Alcohol" = "#F8766D", 
                       "Deprivation_Indices" = "#CD9600", 
                       "Diet_Weekly" = "#7CAE00",
                       "Exercise_Freq" = "#00BE67",
                       "Exercise_MET"  = "#00BFC4",
                       "Internet_Usage" = "#00A9FF",
                       "Smoking" = "#C77CFF",
                       "Vitamins" = "#FF61CC")
  
  #GGPLOT:
  #Added offset to -log10 to allow for inclusion of pvalues exactly = 0.
  p1 <- allAssoc_nonsignificant %>%
    ggplot(aes(x = idx, y = -log10(`Pr(>|t|)_train` + 1e-300) * sign(Estimate_train),
               color = Category_train)) +
    geom_point(size = 2, color = "gray",
               aes(text = paste("ID:", ID, 
                                "<br>Beta:", signif(Estimate_train, digits = 2),
                                "<br>Pval:", signif(`Pr(>|t|)_train`, digits = 2),
                                "<br>N:", samplesize_train,
                                "<br>Category:", Category_train,
                                "<br>Total R2:", signif(adj.R2_train, digits = 2)))) +
    scale_color_discrete(drop = FALSE)
  
  p2 <- p1 +
    geom_point(data = allAssoc_significant,
               size = 2,
               aes(text = paste("ID:", ID, 
                                "<br>Beta:", signif(Estimate_train, digits = 2),
                                "<br>Pval:", signif(`Pr(>|t|)_train`, digits = 2),
                                "<br>N:", samplesize_train,
                                "<br>Category:", Category_train,
                                "<br>Total R2:", signif(adj.R2_train, digits = 2)))) +
    geom_abline(intercept = -log10(pval_thresh + 1e-300), slope = 0, color = "blue", linetype = "dashed") +
    geom_abline(intercept = log10(pval_thresh + 1e-300), slope = 0, color = "blue", linetype = "dashed") +
    theme_minimal() +
    scale_color_manual(values = category_colors) +
    labs(x = "Exposures", y = "-log10(P-value) * sign(Beta)",
         color = "Category") +
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
  
  return(p2)
  
}

# FUNCTION: Save Interactive Plot:
plotsave_app <- function(gg, protID, type){
  gginteract <- ggplotly(gg, tooltip = "text") 
  html_file <- paste0("./Output/App/Interactive/A2/",type,"/",
                      protID,"_",type,"assoc.html")
  
  # Save the Plotly widget without embedding resources (Plotly will be fetched from the CDN)
  saveWidget(gginteract, 
             file = html_file,
             selfcontained = FALSE)
  
  # Make sure all the html files refer to static package folder:
  html_content <- readLines(html_file)
  
  fname <- basename(html_file)
  fname <- gsub(".html","_files", fname)
  html_content <- gsub(fname, 'static', html_content)  # Removes fill.css
  
  # Write the modified HTML content back to the file
  writeLines(html_content, html_file)
  
  # Remove the extra folder created by saveWidget
  #unlink(paste0("./Output/App/Interactive/A2/", type, 
  #              "/", protID, "_", type, "assoc_files"), recursive = TRUE)
  
}


# FUNCTION: Single Omic Signatures
saveOmicSignaturesApp <- function(omicID, EType, Type, plotname){
  suppressWarnings({
    suppressMessages({
      gg1 <- miamiplot_single_omic_rep(protID = omicID, 
                                       Eall = EType,
                                       title = paste0(omicID,": ",plotname))
      plotsave_app(gg1, omicID, Type)
    })
  })
} 
#Example usage:
#saveOmicSignaturesApp("GDF15",Eall, "Type6","E Associations")
#saveOmicSignaturesApp("GDF15", Eall, "Type5","E Associations")

#Save interactive plots for all specifications:
runAllInteractiveAssoc <- function(CovarType){
  # Obtain DF w/ both Train and Test Associations:
  Eall <- replication_all(train = CovarSpecList[[CovarType]]$train[[1]],
                          test = CovarSpecList[[CovarType]]$test[[1]])
  protList <- unique(Eall$omicID)
  
  pblapply(protList, function(x){
    saveOmicSignaturesApp(x, Eall, CovarType,"E Associations")
  })
}

runAllInteractiveAssoc("Type1")
runAllInteractiveAssoc("Type2")
runAllInteractiveAssoc("Type3")
runAllInteractiveAssoc("Type4")
runAllInteractiveAssoc("Type5")
runAllInteractiveAssoc("Type6") #Expected Runtime is 1 Hour...: took 1 hour 10 minutes total

#Extra Info:
#Figured out color scheme for categories:
#library(scales)
#hex <- hue_pal()(8)
