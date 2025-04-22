#Load Libraries
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
#source("./RScripts/Pure_StatGen/Prot_ExPGS/Visualizations/analyProt_PGS_PXS_final.R")

# Define HEAPres:
HEAPres <- HEAPge@HEAP

###### Generalization Metrics #####
#'*Supplementary Figure: GENERALIZATION METRIC FIGURES FOR HEAP*
#'Make separate figures but also have a dataframe with (metric + CI) for use later

# Supplementary Figure 1: Generalization of PXS across train and test splits.

# Function to plot slopes and R2 across train/test splits
plotPXS_generalization <- function(Lasso, R2train, R2test, filename){
  
  # Obtain mean R2 for overall lasso.
  plot_lassoR2 <- Lasso %>%
    group_by(omic) %>%
    summarize(
      mean_trainR2 = mean(train_lasso),
      mean_testR2 = mean(test_lasso)
    )
  
  # Categorize any bad fit proteins:
  outlier_threshold <- 2 * sd(plot_lassoR2$mean_testR2)
  outlier_indices <- which(abs(plot_lassoR2$mean_testR2 - plot_lassoR2$mean_trainR2) > outlier_threshold)
  badfit_proteins <- plot_lassoR2[outlier_indices, ]$omic
  print(badfit_proteins)
  
  #Obtain mean R2 for train and test:
  plot_R2train<- R2train %>%
    group_by(ID) %>%
    summarise(across(everything(), mean, na.rm = TRUE))
  
  plot_R2test <- R2test %>%
    group_by(ID) %>%
    summarise(across(everything(), mean, na.rm = TRUE))
  

  #Combine train and test R2's
  plot_partR2 <- merge(plot_R2train, plot_R2test, by = "ID")
  
  p1 <- plot_partR2 %>% 
    filter(!ID %in% c(badfit_proteins)) %>%
    ggplot(aes(x = G.x, y = G.y)) +
    geom_point(color = "black", size = 3)  +
    stat_poly_line() +
    stat_poly_eq(aes(label = paste(after_stat(eq.label), after_stat(adj.rr.label), sep = "~~~")),
                 formula = y ~ x,  # Formula for regression
                 parse = TRUE,  # Parse the labels
                 color = "black",
                 rr.digits = 2,
                 coef.digits = 2,
                 size = 3) +
    theme_minimal() +
    xlab(expression("Train R"^2)) +
    ylab(expression("Test R"^2)) +
    ggtitle(bquote("G Partitioned"~R^2))
  
  
  p2 <- plot_partR2 %>% 
    filter(!ID %in% c(badfit_proteins)) %>%
    ggplot(aes(x = E.x, y = E.y)) +
    geom_point(color = "black", size = 3)  +
    stat_poly_line() +
    stat_poly_eq(aes(label = paste(after_stat(eq.label), after_stat(adj.rr.label), sep = "~~~")),
                 formula = y ~ x,
                 parse = TRUE,
                 color = "black",
                 rr.digits = 2,
                 coef.digits = 2,
                 size = 3) +
    theme_minimal() +
    xlab(expression("Train R"^2)) +
    ylab(expression("Test R"^2)) +
    ggtitle(bquote("E Partitioned"~R^2))
  
  
  p3 <- plot_partR2 %>% 
    filter(!ID %in% c(badfit_proteins)) %>%
    ggplot(aes(x = GxE.x, y = GxE.y)) +
    geom_point(color = "black", size = 3)  +
    stat_poly_line() +
    stat_poly_eq(aes(label = paste(after_stat(eq.label), 
                                   after_stat(adj.rr.label), 
                                   sep = "~~~")),
                 formula = y ~ x,
                 parse = TRUE,
                 color = "black",
                 rr.digits = 2,
                 coef.digits = 2,
                 size = 3) +
    theme_minimal() +
    xlab(expression("Train R"^2)) +
    ylab(expression("Test R"^2)) +
    ggtitle(bquote("GxE Partitioned"~R^2))
  
  p4 <- (p1 | p2 | p3)
  ggsave(paste0("./Figures/HEAP/SF1/","GeneralizationR2_CovarSpec_",filename,".png"), 
         plot = p4, width = 10, height = 3, dpi = 1000)
  ggsave(paste0("./Figures/HEAP/SF1/","GeneralizationR2_CovarSpec_",filename,".svg"), 
         plot = p4, width = 10, height = 3, dpi = 1000)
  ggsave(paste0("./Figures/HEAP/SF1/","GeneralizationR2_CovarSpec_",filename,".pdf"), 
         plot = p4, width = 10, height = 3, dpi = 1000)
  
  return(p4)
}

lapply(names(HEAPres), function(x) {
  data <- HEAPres[[x]]
  plotPXS_generalization(data$Lasso, data$R2train, data$R2test, x)
})


# Supplementary Figure 1: Generalization Metric: Slope (Created Metric)

# Function to get Slopes of Train/Test split as generalization:
slope_R2 <- function(R2train, R2test, typeID){
  #ID to mark all entries.
  R2train$rID <- 1:nrow(R2train)
  R2test$rID <- 1:nrow(R2test)
  
  #ID to mark data from specific folds.
  R2train$gID <- rep(1:10, nrow(R2train)/10)
  R2test$gID <- rep(1:10, nrow(R2test)/10)
  
  partR2 <- merge(R2train, R2test, by = c("rID", "gID","ID"))
  
  Ecoefs <- sapply(1:10, function(i) {
    fit <- lm(E.y ~ E.x, data = partR2 %>% filter(gID == i))
    unname(fit$coefficients[2])
  })
  
  Gcoefs <- sapply(1:10, function(i) {
    fit <- lm(G.y ~ G.x, data = partR2 %>% filter(gID == i))
    unname(fit$coefficients[2])
  })
  
  GxEcoefs <- sapply(1:10, function(i) {
    fit <- lm(GxE.y ~ GxE.x, data = partR2 %>% filter(gID == i))
    unname(fit$coefficients[2])
  })
  
  storeCoefs <- data.frame(
    G = Gcoefs,
    E = Ecoefs,
    GxE = GxEcoefs
  )
  storeCoefs$ID <- typeID
  return(storeCoefs)
}

# Function to get Confidence Intervals across Folds
get_CI <- function(column){
  stat <- data.frame()
  if(is.character(column)){
    c <- unique(column)
    return(c)
  } else{
    stat <- tibble(
      mean = mean(column),
      sd = sd(column),
      n = length(column),
      se = sd / sqrt(n),
      ci_margin = qt(0.975, df = n - 1) * se,
      #ci_low = mean - qt(0.975, df = n - 1) * se, #95% confidence interval:
      #ci_high = mean + qt(0.975, df = n - 1) * se
      ci = paste0(paste0(round(mean, 4), " ± ", round(ci_margin, 4)))
    )
    return(stat$ci)
  }
}

# Function to process dataframes into ggplot format
process_column <- function(column) {
  tibble(
    mean = as.numeric(str_extract(column, "^[^ ±]+")),
    ci_margin = as.numeric(str_extract(column, "(?<=± )[^ ]+$")),
    ci_low = mean - ci_margin,
    ci_high = mean + ci_margin
  )
}

# Across folds: 1 - 10: Calculate slope of train R2 vs test R2.
plotPXS_slope <- function(HEAPres){
  
  # Obtain slope metrics
  sR2 <- lapply(names(HEAPres), function(x) {
    data <- HEAPres[[x]] 
    slope_R2(data$R2train, data$R2test, typeID = x)
  }) %>% reduce(full_join)
  
  covarSpecR2 <- sR2 %>%
    group_by(ID) %>%
    summarize(across(everything(), get_CI))
  
  
  data_long <- covarSpecR2 %>%
    pivot_longer(cols = c(G, E, GxE), names_to = "Type", values_to = "Value") %>%
    separate(Value, into = c("Mean", "Error"), sep = " ± ") %>%
    mutate(
      Mean = as.numeric(Mean),
      Error = as.numeric(Error)
    )
  
  # Create the ggplot
  slope1 <- ggplot(data_long, aes(x = ID, y = Mean, fill = Type)) +
    geom_bar(stat = "identity", position = position_dodge(), width = 0.7) +
    geom_errorbar(aes(ymin = Mean - Error, ymax = Mean + Error), 
                  position = position_dodge(0.7), width = 0.2) +
    labs(title = "Generalization of Train vs Test across Proteome",
         x = "Covariate Specifications",
         y = "Slope (Coefficient)",
         fill = "Scores") +
    facet_wrap(~ Type, scales = "fixed") + 
    theme_minimal()
  
  ggsave(paste0("./Figures/HEAP/SF1/",
                "GeneralizationR2_Slope.png"), 
         plot = slope1, width = 10, height = 3, dpi = 1000)
  ggsave(paste0("./Figures/HEAP/SF1/",
                "GeneralizationR2_Slope.svg"), 
         plot = slope1, width = 10, height = 3, dpi = 1000)
  ggsave(paste0("./Figures/HEAP/SF1/",
                "GeneralizationR2_Slope.pdf"), 
         plot = slope1, width = 10, height = 3, dpi = 1000)
}
plotPXS_slope(HEAPres)


# Supplementary Figure 1: Additional Generalization Metrics: NMAE, MSE, etc.
# Normalization for MSE, MAE
# Scaled Absolute differences... ^normalization constant.
HEAPmetrics <- function(train, test){
  MSE <- mean((test - train)^2)
  MAE <- mean(abs(test - train))
  NMSE <- MSE / var(test)
  NMAE <- MAE / mean(test)
  
  Hmetrics <- c(MSE,MAE, NMSE, NMAE)
  names(Hmetrics) <- c("MSE","MAE","NMSE","NMAE")
  return(Hmetrics)
}

addn_metrics <- function(R2train, R2test, typeID){
  
  #ID to mark all entries.
  R2train$rID <- 1:nrow(R2train)
  R2test$rID <- 1:nrow(R2test)
  
  #ID to mark data from specific folds.
  R2train$gID <- rep(1:10, nrow(R2train)/10)
  R2test$gID <- rep(1:10, nrow(R2test)/10)
  
  partR2 <- merge(R2train, R2test, by = c("rID", "gID","ID"))
  
  G <- lapply(1:10, function(i) {
    df <- partR2 %>% filter(gID == i)
    HEAPmetrics(df$G.x, df$G.y)
  })
  G <- as.data.frame(do.call("rbind",G))
  colnames(G) <- paste0("G_",colnames(G))
  
  E <- lapply(1:10, function(i) {
    df <- partR2 %>% filter(gID == i)
    HEAPmetrics(df$E.x, df$E.y)
  })
  E <- as.data.frame(do.call("rbind",E))
  colnames(E) <- paste0("E_",colnames(E))
  
  
  GxE <- lapply(1:10, function(i) {
    df <- partR2 %>% filter(gID == i)
    HEAPmetrics(df$GxE.x, df$GxE.y)
  })
  GxE <- as.data.frame(do.call("rbind",GxE))
  colnames(GxE) <- paste0("GxE_",colnames(GxE))
  
  metrics <- as.data.frame(do.call("cbind",list(G,E,GxE)))
  metrics$ID <- typeID
  return(metrics)
}

plotPXS_metrics <- function(HEAPres){
  # Obtain slope metrics
  sR2 <- lapply(names(HEAPres), function(x) {
    data <- HEAPres[[x]] 
    addn_metrics(data$R2train, data$R2test, typeID = x)
  }) %>% reduce(full_join)
  
  covarSpecR2 <- sR2 %>%
    group_by(ID) %>%
    summarize(across(everything(), get_CI))
  
  data_long <- covarSpecR2 %>%
    pivot_longer(cols = c(G_NMAE, E_NMAE, GxE_NMAE), names_to = "Type", values_to = "Value") %>%
    separate(Value, into = c("Mean", "Error"), sep = " ± ") %>%
    mutate(
      Mean = as.numeric(Mean),
      Error = as.numeric(Error)
    )
  
  # Create the ggplot
  NMAE <- ggplot(data_long, aes(x = ID, y = Mean, fill = Type)) +
    geom_bar(stat = "identity", position = position_dodge(), width = 0.7) +
    geom_errorbar(aes(ymin = Mean - Error, ymax = Mean + Error), 
                  position = position_dodge(0.7), width = 0.2) +
    labs(title = "Generalization of Train vs Test NMAE across Proteome",
         x = "Covariate Specifications",
         y = "Normalized Mean Absolute Error") +
    facet_wrap(~ Type, scales = "fixed") + #free_y for change in y axis, fixed for equivalent scales.
    theme_minimal() 
  
  ggsave(paste0("./Figures/HEAP/SF1/",
                "GeneralizationR2_NMAE.png"), 
         plot = NMAE, width = 10, height = 3, dpi = 1000)
  ggsave(paste0("./Figures/HEAP/SF1/",
                "GeneralizationR2_NMAE.svg"), 
         plot = NMAE, width = 10, height = 3, dpi = 1000)
  ggsave(paste0("./Figures/HEAP/SF1/",
                "GeneralizationR2_NMAE.pdf"), 
         plot = NMAE, width = 10, height = 3, dpi = 1000)
  
  return(NMAE)
}
plotPXS_metrics(HEAPres)


# Supplementary Figure 1: Assess Metrics for Total R2 across CovarSpec
#'*Total R2 across CovarSpec*

# Remember: 
# Lasso R2 can be negative - Comparison across CovarSpec for R2 >= 0 is fine. <= 0 is poor fit.
# The partitioned R2 are all >= 0 since we do not penalize for amount of covariates here. (standard R2)

#Comparing R2 train vs test lasso (entire)
addn_metrics_lasso <- function(Lasso, typeID){
  # Obtain mean R2 for overall lasso.
  plot_lassoR2 <- Lasso %>%
    group_by(omic) %>%
    summarize(
      mean_trainR2 = mean(train_lasso),
      mean_testR2 = mean(test_lasso)
    )
  
  # Categorize any bad fit proteins:
  outlier_threshold <- 2 * sd(plot_lassoR2$mean_testR2)
  outlier_indices <- which(abs(plot_lassoR2$mean_testR2 - plot_lassoR2$mean_trainR2) > outlier_threshold)
  badfit_proteins <- plot_lassoR2[outlier_indices, ]$omic
  
  Lasso <- Lasso %>%
    filter(!omic %in% c(badfit_proteins))
  
  #ID to mark all entries.
  Lasso$rID <- 1:nrow(Lasso)
  
  #ID to mark data from specific folds.
  Lasso$gID <- rep(1:10, nrow(Lasso)/10)
  
  All <- lapply(1:10, function(i) {
    df <- Lasso %>% filter(gID == i)
    HEAPmetrics(df$train_lasso, df$test_lasso)
  })
  All <- as.data.frame(do.call("rbind",All))
  colnames(All) <- paste0("All_",colnames(All))
  
  metrics <- All
  metrics$ID <- typeID
  return(metrics)
}

plotPXS_lassofit <- function(HEAPres){
  # Obtain lasso metrics
  sR2 <- lapply(names(HEAPres), function(x) {
    data <- HEAPres[[x]] 
    addn_metrics_lasso(data$Lasso, typeID = x)
  }) %>% reduce(full_join)
  
  covarSpecR2 <- sR2 %>%
    group_by(ID) %>%
    summarize(across(everything(), get_CI))
  
  data_long <- covarSpecR2 %>%
    pivot_longer(cols = c(All_NMAE), names_to = "Type", values_to = "Value") %>%
    separate(Value, into = c("Mean", "Error"), sep = " ± ") %>%
    mutate(
      Mean = as.numeric(Mean),
      Error = as.numeric(Error)
    )
  
  # Create the ggplot
  NMAE_all <- ggplot(data_long, aes(x = ID, y = Mean, fill = Type)) +
    geom_bar(stat = "identity", position = position_dodge(), width = 0.7) +
    geom_errorbar(aes(ymin = Mean - Error, ymax = Mean + Error), 
                  position = position_dodge(0.7), width = 0.2) +
    labs(title = "Generalization: Train/Test R2 across Proteome",
         x = "Covariate Specifications",
         y = "Normalized Mean Absolute Error") +
    facet_wrap(~ Type, scales = "fixed") + #free_y for change in y axis, fixed for equivalent scales.
    theme_minimal() +
    theme(legend.position="none")
  
  ggsave(paste0("./Figures/HEAP/SF1/",
                "GeneralizationR2_NMAE_all.png"), 
         plot = NMAE_all, width = 5, height = 3, dpi = 1000)
  ggsave(paste0("./Figures/HEAP/SF1/",
                "GeneralizationR2_NMAE_all.svg"), 
         plot = NMAE_all, width = 5, height = 3, dpi = 1000)
  ggsave(paste0("./Figures/HEAP/SF1/",
                "GeneralizationR2_NMAE_all.pdf"), 
         plot = NMAE_all, width = 5, height = 3, dpi = 1000)
  
  return(NMAE_all)
}
plotPXS_lassofit(HEAPres)


# Supplementary Figure 1: Reporting Density Shifts of R2:
## Interpretation: As increase addition of covariates --> total test R2 increases.
## Add we condition on more variables the attributed variance of E: R2 decreases.

library(ggridges)
plotPXS_densities <- function(HEAPres){
  
  # purrr::map2 iterate over the list created to generate one dataframe:
  Lasso <- map2_dfr(HEAPres, names(HEAPres), ~ .x$Lasso %>%
                      select(all_of(c("omic","train_lasso","test_lasso"))) %>%
                      mutate(cID = .y))
  
  #Total R2 distribution: ALL
  Lasso_CovarSpec <- Lasso %>%
    group_by(omic, cID) %>%
    summarize(
      mean_trainR2 = mean(train_lasso),
      mean_testR2 = mean(test_lasso)
    )
  
  Lasso_ds1 <- Lasso_CovarSpec %>%
    ggplot(aes(x = mean_testR2, color = cID)) +
    geom_density(show.legend = F) +
    stat_density(aes(x = mean_testR2, color = cID),
                 geom="line",position="identity") +
    xlim(0, 1) +
    theme_minimal() +
    labs(x = "Test: Total R2",
         y = "Density",
         color = "Covar Spec") 
  
  ggsave(paste0("./Figures/HEAP/SF1/",
                "DensityR2_all.png"), 
         plot = Lasso_ds1, width = 5, height = 3, dpi = 1000)
  ggsave(paste0("./Figures/HEAP/SF1/",
                "DensityR2_all.svg"), 
         plot = Lasso_ds1, width = 5, height = 3, dpi = 1000)
  ggsave(paste0("./Figures/HEAP/SF1/",
                "DensityR2_all.pdf"), 
         plot = Lasso_ds1, width = 5, height = 3, dpi = 1000)
  
  
  R2test <- map2_dfr(HEAPres, names(HEAPres), ~ .x$R2test %>%
                       select(all_of(c("ID", "G", "E", "GxE"))) %>%
                       mutate(cID = .y))
  
  R2test_CovarSpec <- R2test %>%
    group_by(ID, cID) %>%
    summarize(
      mean_E = mean(E),
      mean_G = mean(G),
      mean_GxE = mean(GxE)
    )
  
  E_ds1 <- R2test_CovarSpec %>%
    ggplot(aes(x = mean_E, color = cID)) +
    geom_density(show.legend = F) +
    stat_density(aes(x = mean_E, color = cID),
                 geom="line",position="identity") +
    xlim(0, 0.17) +
    theme_minimal() +
    labs(x = "Test: E R2",
         y = "Density",
         color = "Covar Spec")
  
  ggsave(paste0("./Figures/HEAP/SF1/",
                "DensityR2_E.png"), 
         plot = E_ds1, width = 5, height = 3, dpi = 1000)
  ggsave(paste0("./Figures/HEAP/SF1/",
                "DensityR2_E.svg"), 
         plot = E_ds1, width = 5, height = 3, dpi = 1000)
  ggsave(paste0("./Figures/HEAP/SF1/",
                "DensityR2_E.pdf"), 
         plot = E_ds1, width = 5, height = 3, dpi = 1000)
}
plotPXS_densities(HEAPres)

