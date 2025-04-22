#'*Miami Plot showing Replication Hits*

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