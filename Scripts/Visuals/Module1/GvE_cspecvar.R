#Load Libraries and check for dependencies:
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

HEAPres <- HEAPge@HEAP
R2test <- HEAPge@R2ge

#'*FIGURE OUT BETTER VISUAL FOR RESULTS w/ and w/o condition on BMI*
# Function to show CovarSpec variability in top proteins:
CovarSpec_topProt <- function(R2test, nProt = 10){
  top_protein_feat <- R2test %>%
    group_by(ID, cID) %>%
    summarise(mean_E = mean(E), .groups = 'drop') %>%
    group_by(ID) %>%
    summarise(
      agg_E = mean(mean_E), 
      stdError = sd(mean_E) / sqrt(n()),
      CI_margin = qt(0.975, df = n() - 1) * stdError,
      n = n(),
      .groups = 'drop'
    ) %>%
    top_n(nProt, agg_E) %>%  # Select top 10 IDs mean R2
    pull(ID)
  
  p1 <- R2test %>%
    group_by(ID, cID) %>%
    summarise(E = mean(E))  %>% #get mean E: R2 for each specification
    mutate(mean_E = mean(E),
           n = n(),
           stdErr_E = sd(E)/sqrt(n),
           ci_low = mean_E - qt(0.975, df = n-1) * stdErr_E,
           ci_high = mean_E + qt(0.975, df = n-1) * stdErr_E,
    ) %>% #get CI of the E: R2 across specifications
    filter(ID %in% top_protein_feat) %>%
    ggplot(aes(x = reorder(ID, -mean_E), y = mean_E)) +
    #geom_pointrange(aes(y = mean_E, ymin = ci_low, ymax = ci_high), color = "black") +
    geom_linerange(aes(ymin = ci_low, ymax = ci_high), color = "black",
                   size = 0.5) + 
    geom_linerange(aes(xmin = as.numeric(reorder(ID, -mean_E)) - 0.2, 
                       xmax = as.numeric(reorder(ID, -mean_E)) + 0.2),
                   color = "black",
                   size = 1) +
    geom_jitter(aes(x = reorder(ID, -E),
                    y = E,
                    color = cID), width = 0.2, height = 0)  +
    ylim(0, NA) +
    labs(x = "ID", 
         y = "E: R2", 
         color = "Covar Spec") +
    #ggtitle("Top 10 proteins across 10 folds and CovarSpec - E: R2") +
    theme_minimal() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
  #+
  #theme_minimal()
  
  #Save the plots:
  lapply(c("png", "svg", "pdf"), 
         function(fmt) ggsave(paste0("./Figures/HEAP/SF2/","CovarSpec_topProt_",nProt,".", fmt), 
                              plot = p1, width = 6, height = 4, dpi = 1000))
}
CovarSpec_topProt(R2test)
CovarSpec_topProt(R2test, nProt = 20)
CovarSpec_topProt(R2test, nProt = 30)
CovarSpec_topProt(R2test, nProt = 50)


# Function to highlight Similarity and Differences between CovarSpec Variability for E: R2
CovarSpec_PCoA <- function(R2test){
  CovarSpec <- R2test %>%
    group_by(ID, cID) %>%
    summarise(mean_E = mean(E), .groups = 'drop') %>%
    pivot_wider(names_from = cID, values_from = mean_E) %>% #Make matrix of meanE and cID type.
    column_to_rownames(var = "ID")
  
  
  # MDS/PCoA Scaling using Distance Matrix:
  m_dist <- as.matrix(dist(t(CovarSpec), method = "euclidean"))
  #m_dist <- as.matrix(dist(t(CovarSpec), method = "manhattan"))
  mds_result <- cmdscale(m_dist, k = 2, eig = T)
  
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
  
  
  #Save the plots:
  lapply(c("png", "svg", "pdf"), 
         function(fmt) ggsave(paste0("./Figures/HEAP/SF2/",
                                     "CovarSpec_PCoA_E.", fmt), 
                              plot = p1, width = 4, height = 3, dpi = 1000))
  
}
CovarSpec_PCoA(R2test)


#'*Plotting Functions to Understand Covariate Specification Variability of E: R2 across Categories*
categories <- c("Alcohol","Diet_Weekly","Smoking","Exercise_MET",
                "Exercise_Freq","Internet_Usage","Deprivation_Indices","Vitamins")

# purrr::map2 iterate over the list created to generate one dataframe:
R2testCat <- map2_dfr(HEAPres, names(HEAPres), ~ .x$R2testCat %>%
                        select(all_of(c("ID", categories))) %>%
                        mutate(cID = .y))

R2testCat_CI <- R2testCat %>%
  group_by(ID, cID) %>%
  summarise(across(c(categories), mean, .names = "{col}")) %>%  #Mean of each Specification
  mutate(across(c(categories), mean, .names = "mean_{col}")) %>% #Mean across Specification
  mutate(across(c(categories), 
                list(
                  stdErr = ~ sd(.) / sqrt(n()) #, 
                  #ci_low = ~ . - qt(0.975, df = n() - 1) * stdErr, 
                  #ci_high = ~ . + qt(0.975, df = n() - 1) * stdErr
                ), 
                .names = "{col}_{fn}")
  ) %>% #Obtain standard errors
  mutate(
    n = n()
  ) %>%
  mutate(across(starts_with("mean_"),
                list(
                  ci_low = ~ . - qt(0.975, df = n - 1) * get(paste0(sub("mean_", "", cur_column()), "_stdErr")), 
                  ci_high = ~ . + qt(0.975, df = n - 1) * get(paste0(sub("mean_", "", cur_column()), "_stdErr"))
                ),
                .names = "{col}_{fn}")
  ) #obtain CI margins for plotting


library(cowplot)
plot_r2_column <- function(data, column_name, nProt = 10) {
  #R2filt = 0.03
  
  # Dynamically build aesthetic mappings for mean, ci_low, ci_high
  mean_col <- paste0("mean_", column_name)
  ci_low_col <- paste0(mean_col, "_ci_low")
  ci_high_col <- paste0(mean_col, "_ci_high")
  
  #Define amount of proteins to show
  filt_prot <- unique(data[order(data[[mean_col]], decreasing = T),]$ID)[1:nProt]
  
  
  # Create the plot
  p1 <- data %>%
    #filter(!!sym(mean_col) > R2filt) %>%
    filter(ID %in% filt_prot) %>%
    ggplot(aes(x = reorder(ID, -!!sym(mean_col)), y = !!sym(mean_col))) +
    geom_linerange(aes(ymin = !!sym(ci_low_col), ymax = !!sym(ci_high_col)), 
                   color = "black", size = 0.5) +
    geom_linerange(aes(xmin = as.numeric(reorder(ID, -!!sym(mean_col))) - 0.2, 
                       xmax = as.numeric(reorder(ID, -!!sym(mean_col))) + 0.2),
                   color = "black", size = 1) +
    geom_jitter(aes(x = reorder(ID, -!!sym(column_name)),
                    y = !!sym(column_name),
                    color = cID), width = 0.2, height = 0) +
    ylim(0, 0.1) +
    labs(x = "ID", 
         y = paste(column_name, ": R2"), 
         color = "Covar Spec") +
    theme_minimal() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          legend.position = "none")
  
  legend <- get_legend(
    p1 + theme(legend.position = "right")  # Set legend position
  )
  
  
  lapply(c("png", "svg", "pdf"), 
         function(fmt) ggsave(paste0("./Figures/HEAP/SF2/",
                                     "CovarSpec_legend.", fmt), 
                              plot = p1, width = 1, height = 4, dpi = 1000))
  
  p2 <- p1 + 
    theme(axis.title.x = element_blank())
  return(p2)
}
#"Alcohol","Diet_Weekly","Smoking","Exercise_MET",
#"Exercise_Freq","Internet_Usage","Deprivation_Indices","Vitamins"
library(patchwork)
p1 <- plot_r2_column(R2testCat_CI , "Exercise_Freq", nProt = 5)
p2 <- plot_r2_column(R2testCat_CI , "Exercise_MET", nProt = 5) 
p3 <- plot_r2_column(R2testCat_CI , "Diet_Weekly", nProt = 5)


p4 <- plot_r2_column(R2testCat_CI , "Alcohol", nProt = 5)
p5 <- plot_r2_column(R2testCat_CI , "Smoking", nProt = 5)
p6 <- plot_r2_column(R2testCat_CI , "Vitamins", nProt = 5) 
#p7 <- plot_r2_column(R2testCat_CI , "Deprivation_Indices", nProt = 1)
#p8 <- plot_r2_column(R2testCat_CI , "Internet_Usage", nProt = 5)

cow1 <- (p1 | p2 | p3) / (p4 | p5 | p6)

lapply(c("png", "svg", "pdf"), 
       function(fmt) ggsave(paste0("./Figures/HEAP/SF2/",
                                   "CovarSpec_Catvar.", fmt), 
                            plot = cow1, width = 8, height = 6, dpi = 1000))
