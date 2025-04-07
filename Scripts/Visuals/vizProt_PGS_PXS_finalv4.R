library(data.table)
library(tidyverse)
library(ggplot2)
library(ggpmisc)
library(ggrepel)
library(gridExtra)

#Addn Libraries for Visualization:
library(patchwork)

setwd("/n/groups/patel/shakson_ukb/UK_Biobank/")

#'*Option*
#SAVE AS VECTOR GRAPHICS
#Examples:
#library(svglite)
#ggsave(paste0(".svg"),plot = E_ds1, width = 5, height = 3, dpi = 500)
#ggsave(paste0(".pdf"), plot = E_ds1, width = 5, height = 3, dpi = 500)
##Visualization THINGS: TITLES DO I NEED TO SPECIFY THE 'TEST SET' ALL THE TIME

#'*GOALS*
#Things to look at: Distribution of Total R2 from lasso
#Train vs Test R2 across covariate specifications
#Central tendency of R2 for G, E, and GxE terms
#Function types: summarise vs summarize
# purrr::map2 iterate over the list created to generate one dataframe:

#'*Additional Notes*'
#State in paper polygenic GxE is not well powered in this study
#- given analysis of univariate tests.



#'*Runtime Comments on Results from HEAP: using nested CV*
#'General Comments about Runtime of Covariate Specifications:
#'Type 1-3 - works with 1 hour, 10 minutes
#'Type 4-5 - needs more time to run: set to maybe 1 hour, 30 minutes.
loadPXSmetrics <- function(covarType){
  PXSmetrics <- list()
  
  Lasso <- data.frame()
  R2train <- data.frame()
  R2test <- data.frame()
  R2trainCat <- data.frame()
  R2testCat <- data.frame()
  error_files <- c()
  
  for(i in 1:400){ #200 is based on the split designated in the bash script.
    tryCatch({
      Lasso_idx <- fread(file=paste0("/n/groups/patel/shakson_ukb/UK_Biobank/Data/Parallel/Prot_PGSPXS_final/", covarType,"/","lassofit_",i,".txt"))
      R2train_idx <- fread(file=paste0("/n/groups/patel/shakson_ukb/UK_Biobank/Data/Parallel/Prot_PGSPXS_final/", covarType,"/","R2train_",i,".txt"))
      R2test_idx <- fread(file=paste0("/n/groups/patel/shakson_ukb/UK_Biobank/Data/Parallel/Prot_PGSPXS_final/", covarType,"/","R2test_",i,".txt"))
      R2trainCat_idx <- fread(file=paste0("/n/groups/patel/shakson_ukb/UK_Biobank/Data/Parallel/Prot_PGSPXS_final/", covarType,"/","R2trainCat_",i,".txt"))
      R2testCat_idx <- fread(file=paste0("/n/groups/patel/shakson_ukb/UK_Biobank/Data/Parallel/Prot_PGSPXS_final/", covarType,"/","R2testCat_",i,".txt"))
      
      
      #Updates:
      Lasso <- rbind(Lasso, Lasso_idx)
      R2train <- rbind(R2train, R2train_idx)
      R2test <- rbind(R2test, R2test_idx)
      R2trainCat <- rbind(R2trainCat, R2trainCat_idx)
      R2testCat <- rbind(R2testCat, R2testCat_idx)
      
    }, error = function(e) {
      # Handle the error (e.g., print an error message)
      message(paste("Error reading file:",i))
      error_files <<- append(error_files, i)
    })
  }
  
  print(error_files)
  PXSmetrics[["Lasso"]] <- Lasso
  PXSmetrics[["R2train"]] <- R2train
  PXSmetrics[["R2test"]] <- R2test
  PXSmetrics[["R2trainCat"]] <- R2trainCat
  PXSmetrics[["R2testCat"]] <- R2testCat
  
  return(PXSmetrics)
}
PXS_type1 <- loadPXSmetrics(covarType = "Type1")
PXS_type2 <- loadPXSmetrics(covarType = "Type2")
PXS_type3 <- loadPXSmetrics(covarType = "Type3")
PXS_type4 <- loadPXSmetrics(covarType = "Type4")
PXS_type5 <- loadPXSmetrics(covarType = "Type5")


#'*HEAPres Data Structure*
HEAPres <- list(PXS_type1, PXS_type2, PXS_type3, PXS_type4, PXS_type5)
names(HEAPres) <- c("Type1", "Type2", "Type3", "Type4", "Type5")


#### Simple Statistics & Plots on Average R2 across Specifications  ######

## Get mean of R2 across folds and specifications:
R2test <- map2_dfr(HEAPres, names(HEAPres), ~ .x$R2test %>%
                     select(all_of(c("ID", "G", "E", "GxE"))) %>%
                     mutate(cID = .y))

R2test_GvE <- R2test %>%
                group_by(ID, cID) %>%
                summarise(across(c(G, E, GxE), mean), .groups = 'drop') %>%
                group_by(ID) %>%
                summarise(across(c(G, E, GxE), mean))

# Get summary stats:
# Capture the summary output as strings
summary_G <- capture.output(summary(R2test_GvE$G))
summary_E <- capture.output(summary(R2test_GvE$E))
summary_GxE <- capture.output(summary(R2test_GvE$GxE))

# Combine all summaries into one text
summary_text <- c("Summary of G:\n", summary_G, "\n\n", 
                  "Summary of E:\n", summary_E, "\n\n", 
                  "Summary of GxE:\n", summary_GxE)
# Write the summaries to a text file
writeLines(summary_text, "./Figures/HEAP/F1/R2test_GvE_summary.txt")



# Pie Chart of Average R2 (G and E):
library(gridExtra)
plot_pie <- function(df, column, breaks, labels, title) {
  df <- df %>%
    mutate(binned = cut(.data[[column]], breaks = breaks, labels = FALSE)) %>%
    count(binned) %>%
    mutate(percentage = n / sum(n) * 100, 
           binned = labels[binned])
  
  ggplot(df, aes(x="", y=n, fill=binned)) +
    geom_bar(stat="identity", width=1, color="white") +
    coord_polar("y", start=0) +
    geom_text(aes(label = n), position = position_stack(vjust = 0.5), color = "black") +
    geom_col(aes(x=-1,y=0)) +
    theme_void() +
    ggtitle(title) +
    labs(fill = expression(R^2)) +
    theme(plot.title = element_text(hjust=0.5))
}
PieChart_GvEspec <- function(type = "AverageSpec"){
  # Genetics plot
  gg1 <- plot_pie(R2test_GvE, "G", c(0, 0.01, 0.1, 0.25, 0.5, 1), 
                  c("0-0.01", "0.01-0.1", "0.1-0.25", "0.25-0.5", "0.5-1"), "Genetics")
  
  # Environment plot
  gg2 <- plot_pie(R2test_GvE %>% filter(E > 0), "E", c(0, 0.01, 0.05, 0.2), 
                  c("0-0.01", "0.01-0.05", "0.05-0.2"), "Environment")
  
  # Combine and save plots
  gg3 <- grid.arrange(gg1, gg2, nrow = 2)
  ggsave(paste0("./Figures/HEAP/F1/","EvsG_piechart_",type,".png"), 
         plot = gg3, width = 4, height = 6, dpi = 1000)
  ggsave(paste0("./Figures/HEAP/F1/","EvsG_piechart_",type,".svg"), 
         plot = gg3, width = 4, height = 6, dpi = 1000)
  ggsave(paste0("./Figures/HEAP/F1/","EvsG_piechart_",type,".pdf"), 
         plot = gg3, width = 4, height = 6, dpi = 1000)
}
PieChart_GvEspec()


#'*Get Covariates Specific Metrics:*
plot_R2test <- PXS_type5$R2test

# Create finalized matrix of R2s by each feature.
# Take average R2 across folds:
plot_R2test_demo <- plot_R2test %>%
    mutate(sum_gPCs = rowSums(across(paste0("genetic_principal_components_f22009_0_",1:40))),
           medications = rowSums(across(c("combined_Blood_pressure_medication",
                                          "combined_Hormone_replacement_therapy",         
                                          "combined_Oral_contraceptive_pill_or_minipill",
                                          "combined_Insulin",
                                          "combined_Cholesterol_lowering_medication",
                                          "combined_Do_not_know",
                                          "combined_None_of_the_above",
                                          "combined_Prefer_not_to_answer" 
                                          
           )))) %>%
    group_by(ID) %>%
    summarise(across(everything(), mean))


# Convert wide to long format
plot_R2test_ht <- plot_R2test_demo %>%
  select(c("ID","G","E","GxE","fasting_time_f74_0_0","body_mass_index_bmi_f23104_0_0",
           "age_when_attended_assessment_centre_f21003_0_0","sex_f31_0_0",
           "uk_biobank_assessment_centre_f54_0_0","sum_gPCs",
           "medications",
           "age2","age_sex","age2_sex")) %>%
  mutate(totalR2 = rowSums(across(-ID)),
         GEtotR2 = rowSums(across(c("G","E","GxE"))))

# Summary of low R2 proteins and potential categories these proteins are involved in:
summary(plot_R2test_ht[c("totalR2", "GEtotR2")])
sum(plot_R2test_ht$GEtotR2 < 0.01)

lowR2prot <- plot_R2test_ht %>%
                filter(GEtotR2 < 0.01) %>%
                pull(ID)

plot_R2test_ht %>%
  filter(ID %in% lowR2prot & 
           (age_when_attended_assessment_centre_f21003_0_0 > 0.01 | 
              sex_f31_0_0 > 0.01 | 
              body_mass_index_bmi_f23104_0_0 > 0.01)) %>%
  nrow()



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
    #geom_abline(slope=1, intercept = 0) +
    stat_poly_line() +
    stat_poly_eq(aes(label = paste(after_stat(eq.label), after_stat(adj.rr.label), sep = "~~~")),
                 formula = y ~ x,  # Specify the formula for regression
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
    #geom_abline(slope=1, intercept = 0) +
    stat_poly_line() +
    stat_poly_eq(aes(label = paste(after_stat(eq.label), after_stat(adj.rr.label), sep = "~~~")),
                 formula = y ~ x,  # Specify the formula for regression
                 parse = TRUE,  # Parse the labels
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
    #geom_abline(slope=1, intercept = 0) +
    stat_poly_line() +
    stat_poly_eq(aes(label = paste(after_stat(eq.label), 
                                   after_stat(adj.rr.label), 
                                   sep = "~~~")),
                 #paste(after_stat(eq.label), "~~~", after_stat(adj.rr.label), sep = "\n")
                 #paste(after_stat(eq.label), after_stat(adj.rr.label), sep = "~~~")
                 formula = y ~ x,  # Specify the formula for regression
                 parse = TRUE,  # Parse the labels
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
    plot = p4, width = 10, height = 3, dpi = 500)
  ggsave(paste0("./Figures/HEAP/SF1/","GeneralizationR2_CovarSpec_",filename,".svg"), 
         plot = p4, width = 10, height = 3, dpi = 500)
  ggsave(paste0("./Figures/HEAP/SF1/","GeneralizationR2_CovarSpec_",filename,".pdf"), 
         plot = p4, width = 10, height = 3, dpi = 500)
  
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
         plot = slope1, width = 10, height = 3, dpi = 500)
  ggsave(paste0("./Figures/HEAP/SF1/",
                "GeneralizationR2_Slope.svg"), 
         plot = slope1, width = 10, height = 3, dpi = 500)
  ggsave(paste0("./Figures/HEAP/SF1/",
                "GeneralizationR2_Slope.pdf"), 
         plot = slope1, width = 10, height = 3, dpi = 500)
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
         plot = NMAE, width = 10, height = 3, dpi = 500)
  ggsave(paste0("./Figures/HEAP/SF1/",
                "GeneralizationR2_NMAE.svg"), 
         plot = NMAE, width = 10, height = 3, dpi = 500)
  ggsave(paste0("./Figures/HEAP/SF1/",
                "GeneralizationR2_NMAE.pdf"), 
         plot = NMAE, width = 10, height = 3, dpi = 500)
  
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
         plot = NMAE_all, width = 5, height = 3, dpi = 500)
  ggsave(paste0("./Figures/HEAP/SF1/",
                "GeneralizationR2_NMAE_all.svg"), 
         plot = NMAE_all, width = 5, height = 3, dpi = 500)
  ggsave(paste0("./Figures/HEAP/SF1/",
                "GeneralizationR2_NMAE_all.pdf"), 
         plot = NMAE_all, width = 5, height = 3, dpi = 500)
  
  return(NMAE_all)
}
plotPXS_lassofit(HEAPres)


# Supplementary Figure 1: Reporting Density Shifts of R2:
## Interpretation: As increase addition of covariates --> total test R2 increases.
## Add we condition on more variables the attributed variance of E: R2 decreases.

##TO-DO later think about statistical tests (when the R2 are nested estimates and correlated)
#Alternatives below rely on independence
#ks.test(list1, list2)
#kruskal.test(mean_testR2 ~ cID, data = Lasso_CovarSpec)
#wilcox.test() 

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
         plot = Lasso_ds1, width = 5, height = 3, dpi = 500)
  ggsave(paste0("./Figures/HEAP/SF1/",
                "DensityR2_all.svg"), 
         plot = Lasso_ds1, width = 5, height = 3, dpi = 500)
  ggsave(paste0("./Figures/HEAP/SF1/",
                "DensityR2_all.pdf"), 
         plot = Lasso_ds1, width = 5, height = 3, dpi = 500)
  
  
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
         plot = E_ds1, width = 5, height = 3, dpi = 500)
  ggsave(paste0("./Figures/HEAP/SF1/",
                "DensityR2_E.svg"), 
         plot = E_ds1, width = 5, height = 3, dpi = 500)
  ggsave(paste0("./Figures/HEAP/SF1/",
                "DensityR2_E.pdf"), 
         plot = E_ds1, width = 5, height = 3, dpi = 500)
}
plotPXS_densities(HEAPres)


#'*TO-DO: Fix Warning Messages*
#####Warning to fix later:
# Previously
#across(a:b, mean, na.rm = TRUE)
# Now
#across(a:b, \(x) mean(x, na.rm = TRUE))

#### Plots on Demographic/Covariate Info ####

#'*Figures on R2s of Demographic Info*
# plot_R2test1 <- PXS_type5$R2test
# 
# #Take average of R2 across folds:
# plot_R2test1 <- plot_R2test1 %>%
#                   group_by(ID) %>%
#                   summarise(across(everything(), mean))
# 
# cormat <- cor(plot_R2test1[,-1])
# cor.test(plot_R2test1$E, plot_R2test1$body_mass_index_bmi_f23104_0_0)
# cor.test(plot_R2test1$E, plot_R2test1$G)
# cor.test(plot_R2test1$E, plot_R2test1$age_when_attended_assessment_centre_f21003_0_0)

#'*Understand Distributions of Covariate R2s:*
#'Use the most adjusted model:
plot_R2test <- PXS_type5$R2test

# Create finalized matrix of R2s by each feature.
plot_R2test_demo <- plot_R2test %>%
  mutate(sum_gPCs = rowSums(across(paste0("genetic_principal_components_f22009_0_",1:40))),
         medications = rowSums(across(c("combined_Blood_pressure_medication",
                                        "combined_Hormone_replacement_therapy",         
                                        "combined_Oral_contraceptive_pill_or_minipill",
                                        "combined_Insulin",
                                        "combined_Cholesterol_lowering_medication",
                                        "combined_Do_not_know",
                                        "combined_None_of_the_above",
                                        "combined_Prefer_not_to_answer" 
           
         ))))  %>%
        group_by(ID) %>%
        summarise(across(everything(), mean))


# Convert wide to long format
plot_R2test_ht <- plot_R2test_demo %>%
  select(c("ID","G","E","GxE","fasting_time_f74_0_0","body_mass_index_bmi_f23104_0_0",
           "age_when_attended_assessment_centre_f21003_0_0","sex_f31_0_0",
           "uk_biobank_assessment_centre_f54_0_0","sum_gPCs",
           "medications",
           "age2","age_sex","age2_sex")) %>%
  mutate(totalR2 = rowSums(across(-ID)))


#Goal: Show the top hits for high R2 proteins by topic:
plot_demoR2 <- function(df, label_id, label_name){
  gg1 <- df %>%
          mutate(label = ifelse(rank(-.data[[label_id]]) <= 20, ID, "")) %>%
          ggplot(aes(x = G, y = .data[[label_id]], label = label)) +
          geom_point() +
          geom_text_repel() +
          #xlab(expression("Polygenic Score "~R^2)) +
          xlab(expression("G: PGS "~R^2)) +
          ylab(bquote(.(label_name) ~ R^2))
 
  return(gg1) 
}

d1 <- plot_demoR2(plot_R2test_ht, label_id = "age_when_attended_assessment_centre_f21003_0_0", label_name = "Age")
d2 <- plot_demoR2(plot_R2test_ht, label_id = "sex_f31_0_0", label_name = "Sex")
d3 <- plot_demoR2(plot_R2test_ht, label_id = "body_mass_index_bmi_f23104_0_0", label_name = "BMI")
d4 <- plot_demoR2(plot_R2test_ht, label_id = "fasting_time_f74_0_0", label_name = "Fasting Time")
d5 <- plot_demoR2(plot_R2test_ht, label_id = "uk_biobank_assessment_centre_f54_0_0", label_name = "UKB Assessment Center")
d6 <- plot_demoR2(plot_R2test_ht, label_id = "sum_gPCs", label_name = "Genetic PCs")

d7 <- plot_demoR2(plot_R2test_ht, label_id = "medications", label_name = "Medications")


demo1 <- (d1 | d2 | d3 ) /
  (d4 | d5 | d6)


lapply(c("png", "svg", "pdf"), 
       function(fmt) ggsave(paste0("./Figures/HEAP/SF2/", "DemographicR2_v1.", fmt), 
                            plot = demo1, width = 10, height = 8, dpi = 500))

demo2 <- (d1 | d2 | d3) /
  (d4 | d5 | d7)

lapply(c("png", "svg", "pdf"), 
       function(fmt) ggsave(paste0("./Figures/HEAP/SF2/", "DemographicR2_v2.", fmt), 
                            plot = demo2, width = 10, height = 8, dpi = 500))

#Distribution of total R2 across all components:
##ggridges:
library(ggridges)

feat <- c(Genetics = "G", Exposures = "E", GxE = "GxE",
          `Fasting Time` = "fasting_time_f74_0_0", BMI = "body_mass_index_bmi_f23104_0_0",
          Age = "age_when_attended_assessment_centre_f21003_0_0",
          Sex = "sex_f31_0_0",
          Center = "uk_biobank_assessment_centre_f54_0_0", `Genetic PCs` = "sum_gPCs",
          Medications = "medications", Age2 = "age2", AgexSex = "age_sex",
          Age2xSex = "age2_sex", `Total R2` = "totalR2")

ridgeplotR2 <- plot_R2test_ht %>% 
  rename(all_of(feat)) %>%
  pivot_longer(cols = c("Genetics","Exposures","GxE","Fasting Time","BMI",
                        "Age","Sex",
                        "Center","Genetic PCs",
                        "Medications",
                        "Age2","AgexSex","Age2xSex", "Total R2"))

highlight_categories <- c("Total R2","Genetics", "Exposures", "GxE")
all_categories = c("Genetics","Exposures","GxE","Fasting Time","BMI",
                   "Age","Sex",
                   "Center","Genetic PCs",
                   "Medications",
                   "Total R2")

rp1 <- ridgeplotR2 %>%
  filter(!(name %in% c("Age2xSex", "Age2", "AgexSex"))) %>%
  ggplot(aes(x = value, y = fct_reorder(name, value, .fun = mean), 
             fill = fct_reorder(name, value, .fun = mean))) +
  geom_density_ridges(
    jittered_points = TRUE,
    alpha = 0.5, scale = 0.8, point_size = 0.05,
    position = position_raincloud(width = 0.01, height = 0.15,
                                  ygap = 0.05)
  ) +
  scale_fill_manual(
    values = c(
      setNames(c("red","blue", "green", "magenta"), highlight_categories), # Highlighted categories
      setNames(rep("gray", length(setdiff(all_categories, highlight_categories))),
               setdiff(all_categories, highlight_categories)) # Grayed-out categories
    )) +
  theme_minimal() +
  labs(x = bquote(R^2), y = NULL) +
  theme(legend.position="none")

lapply(c("png", "svg", "pdf"), 
       function(fmt) ggsave(paste0("./Figures/HEAP/F1/", "PartR2Dens_v1.", fmt), 
                            plot = rp1, width = 4, height = 4, dpi = 500))



#'*In Paper*
highlight_categories <- c("Genetics", "Exposures", "GxE")

rp2 <- ridgeplotR2 %>%
  filter(!(name %in% c("Age2xSex", "Age2", "AgexSex","Total R2"))) %>%
  ggplot(aes(x = value, y = fct_reorder(name, value, .fun = mean), 
             fill = fct_reorder(name, value, .fun = mean))) +
  geom_density_ridges(
    panel_scaling = TRUE,
    jittered_points = TRUE,
    alpha = 0.7, scale = 0.9, point_size = 0.05,
    position = position_raincloud(width = 0.01, height = 0.15,
                                  ygap = 0.05)
  ) +
  scale_fill_manual(
    values = c(
      setNames(c("blue", "green", "magenta"), highlight_categories), # Highlighted categories
      setNames(rep("gray", length(setdiff(all_categories, highlight_categories))),
               setdiff(all_categories, highlight_categories)) # Grayed-out categories
    )) +
  theme_minimal() +
  labs(x = bquote(R^2), y = NULL) +
  theme(legend.position="none")

lapply(c("png", "svg", "pdf"), 
       function(fmt) ggsave(paste0("./Figures/HEAP/F1/", "PartR2Dens_v2.", fmt), 
                            plot = rp2, width = 4, height = 4, dpi = 500))


#### COVARIATE SPECIFICATION VARIABILITY ANALYSIS ####
#'*Multiple covariate specification ANALYSIS*
library(dplyr)
library(purrr)

# Define a function to compute summary statistics
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

# Function to split string and calculate numeric values
process_column <- function(column) {
  tibble(
    mean = as.numeric(str_extract(column, "^[^ ±]+")),
    ci_margin = as.numeric(str_extract(column, "(?<=± )[^ ]+$")),
    ci_low = mean - ci_margin,
    ci_high = mean + ci_margin
  )
}

# Create one structure for results from covariate specifications:
R2test <- map2_dfr(HEAPres, names(HEAPres), ~ .x$R2test %>%
                     select(all_of(c("ID", "G", "E", "GxE"))) %>%
                     mutate(cID = .y))


#'*Plotting Functions to Understand Covariate Specification Variability of E: R2*
# Function to plot variability in GxE plot
CovarSpec_GvEplot <- function(R2test){
  #R2test is dataframe of E and G R2 across proteins and covarSpec type:
  R2test_GvE <- R2test %>%
    group_by(ID, cID) %>%
    summarise(mean_G = mean(G),
              mean_E = mean(E), 
              mean_GxE = mean(GxE),
              .groups = 'drop')
  
  R2test_CI <- R2test_GvE %>%
    group_by(ID) %>%
    summarize(across(c("mean_G","mean_E","mean_GxE"), get_CI))
  
  plot_testR2_CI  <- R2test_CI %>%
    rowwise() %>%
    mutate(G_stats = list(process_column(mean_G)),
           E_stats = list(process_column(mean_E)),
           GxE_stats = list(process_column(mean_GxE))) %>%
    # Unnest the statistics into separate columns
    unnest(c(G_stats, E_stats, GxE_stats), names_sep = "_") 
  
  #Highlight specific proteins with high E: R2
  subset_test <- plot_testR2_CI %>% 
    filter((E_stats_mean > G_stats_mean) & (E_stats_mean > 0.09))
  
  #Cutoff for y-axis:
  y_cutoff = 0.175
  
  #Without Error Bars:
  p1 <- plot_testR2_CI  %>% 
    ggplot(aes(x = G_stats_mean, y = E_stats_mean, label = ID)) + 
    geom_point() +
    geom_abline(slope=1, intercept = 0) +
    geom_ribbon(aes(ymin = 0, ymax = ifelse(G_stats_mean <= y_cutoff, G_stats_mean, y_cutoff), fill = "Genetics Driven"), alpha = 0.2) +
    geom_ribbon(aes(ymin = ifelse(G_stats_mean <= y_cutoff, G_stats_mean, y_cutoff), ymax = y_cutoff, fill = "Environment Driven"), alpha = 0.2) +
    scale_fill_manual(values = c("Genetics Driven" = "blue", "Environment Driven" = "green")) +
    ylim(c(0, y_cutoff)) +
    geom_text_repel(data=subset_test, aes(x = G_stats_mean,y = E_stats_mean,label= ID),
                    size = 4,
                    max.overlaps = 20,
                    min.segment.length = 0.025,
                    force = 10) +
    xlab(expression("G: PGS "~R^2)) +
    ylab(expression("E: PXS "~R^2)) +
    theme_minimal() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    guides(fill=guide_legend(title="Region")) +
    ggtitle(expression("Test Set: Partitioned"~R^2)) +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "bottom") 

  ggsave(paste0("./Figures/HEAP/F1/","GvE_R2_CovarSpec.png"), 
         plot = p1, width = 6, height = 6, dpi = 500)
  ggsave(paste0("./Figures/HEAP/F1/","GvE_R2_CovarSpec.svg"), 
         plot = p1, width = 6, height = 6, dpi = 500)
  ggsave(paste0("./Figures/HEAP/F1/","GvE_R2_CovarSpec.pdf"), 
         plot = p1, width = 6, height = 6, dpi = 500)
  
  
  #With Error Bars: All Covariate Specifications
  p2 <- plot_testR2_CI  %>% 
    ggplot(aes(x = G_stats_mean, y = E_stats_mean, label = ID)) + 
    geom_point() +
    geom_abline(slope=1, intercept = 0) +
    geom_ribbon(aes(ymin = 0, ymax = ifelse(G_stats_mean <= y_cutoff, G_stats_mean, y_cutoff), fill = "Genetics Driven"), alpha = 0.2) +
    geom_ribbon(aes(ymin = ifelse(G_stats_mean <= y_cutoff, G_stats_mean, y_cutoff), ymax = y_cutoff, fill = "Environment Driven"), alpha = 0.2) +
    scale_fill_manual(values = c("Genetics Driven" = "blue", "Environment Driven" = "green")) +
    ylim(c(0, y_cutoff)) +
    geom_text_repel(data=subset_test, aes(x = G_stats_mean,y = E_stats_mean,label= ID),
                    size = 3.5) +
    geom_pointrange(aes(ymin = E_stats_ci_low, ymax = E_stats_ci_high)) +
    geom_pointrange(aes(xmin = G_stats_ci_low, xmax = G_stats_ci_high)) +
    #geom_label() + 
    xlab(expression("G: PGS "~R^2)) +
    ylab(expression("E: PXS "~R^2)) +
    theme_minimal() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    guides(fill=guide_legend(title="Region")) +
    ggtitle(expression("Test Set: Partitioned"~R^2)) +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "bottom") 
  
  ggsave(paste0("./Figures/HEAP/A1/","GvE_R2_CovarSpec_wErrorBars.png"), 
         plot = p2, width = 6, height = 6, dpi = 500)
  ggsave(paste0("./Figures/HEAP/A1/","GvE_R2_CovarSpec_wErrorBars.svg"), 
         plot = p2, width = 6, height = 6, dpi = 500)
  ggsave(paste0("./Figures/HEAP/A1/","GvE_R2_CovarSpec_wErrorBars.pdf"), 
         plot = p2, width = 6, height = 6, dpi = 500)
}
CovarSpec_GvEplot(R2test)


#'*Function to count number of proteins where E > G*
#R2test is dataframe of E and G R2 across proteins and covarSpec type:
R2test_CI <- R2test %>%
                group_by(ID, cID) %>% #Mean across test folds
                summarise(mean_G = mean(G),
                          mean_E = mean(E), 
                          mean_GxE = mean(GxE),
                          .groups = 'drop') %>%
                group_by(ID) %>% #get CIs across covar spec
                summarize(across(c("mean_G","mean_E","mean_GxE"), get_CI))

plot_testR2_CI  <- R2test_CI %>%
  rowwise() %>%
  mutate(G_stats = list(process_column(mean_G)),
         E_stats = list(process_column(mean_E)),
         GxE_stats = list(process_column(mean_GxE))) %>%
  # Unnest the statistics into separate columns
  unnest(c(G_stats, E_stats, GxE_stats), names_sep = "_") 

#Highlight specific proteins with high E: R2
Eproteins <- plot_testR2_CI %>% 
                filter(E_stats_mean > G_stats_mean)
write.csv(Eproteins, file = "./Figures/HEAP/T1/Eproteins.csv")

#Stats:
# Capture the summary output as strings
summary1 <- capture.output(summary(Eproteins$E_stats_mean))
summary2 <- capture.output(length(unique(Eproteins$ID)))

# Combine all summaries into one text
summary_text <- c("Summary of E Driven Proteins:\n", summary1, "\n\n", 
                  "Number of E Driven Proteins:\n", summary2, "\n\n")
# Write the summaries to a text file
writeLines(summary_text, "./Figures/HEAP/T1/EDrivenProteins_summary.txt")



#Obtain proteins:

# Function to plot variability of E, G, both separately
#Renamed Standard Error to "Variability across Specification"
CovarSpec_varplot <- function(R2test){
  #'*Good Plot for overall E: R2 variability*
  E_varplot <- R2test %>%
    group_by(ID, cID) %>%
    summarise(
      mean_E = mean(E)
    ) %>%
    summarise(
      agg_E = mean(mean_E), 
      stdError_E = sd(mean_E)/sqrt(n()),
      CI_margin_E = qt(0.975, df = n() - 1) * stdError_E,
      n = n()
    )
  
  p1 <- E_varplot %>%
    ggplot(aes(x = agg_E, y = stdError_E, label = ID)) +
    geom_point() +
    xlab(expression("E: "~R^2)) +
    ylab("Std. Error") +
    geom_text_repel(aes(label = ifelse(agg_E > 0.07, ID, ""))) +
    theme_minimal()
  
  G_varplot <- R2test %>%
    group_by(ID, cID) %>%
    summarise(
      mean_G = mean(G)
    ) %>%
    summarise(
      agg_G = mean(mean_G), 
      stdError_G = sd(mean_G)/sqrt(n()),
      CI_margin_G = qt(0.975, df = n() - 1) * stdError_G,
      n = n()
    )
  
  p2 <- G_varplot %>%
    ggplot(aes(x = agg_G, y = stdError_G, label = ID)) +
    geom_point() +
    xlab(expression("G: "~R^2)) +
    ylab("Std. Error") +
    geom_text_repel(aes(label = ifelse(stdError_G > 5e-4, ID, ""))) +
    theme_minimal() 
  
  #'*Combine plots:*
  varplot <- merge(E_varplot, G_varplot, by = "ID")
  
  p3 <- varplot %>%
    ggplot(aes(x = stdError_G, y = stdError_E, color = agg_E)) +
    geom_point() +
    #scale_color_gradient(low = "", high = "green") +
    scale_color_viridis_c(option = "viridis") +
    #geom_text_repel(aes(label = ifelse(agg_E > 0.07, ID, ""))) +
    xlab("Std. Error: G") +
    ylab("Std. Error: E") +
    labs(color='E: R2')  +
    #ggtitle("Variability of G and E R2 across Covariate Specifications") +
    theme_minimal()
  
  #Save the plots:
  lapply(c("png", "svg", "pdf"), 
         function(fmt) ggsave(paste0("./Figures/HEAP/F1/","Evar.", fmt), 
                              plot = p1, width = 3, height = 3, dpi = 500))
  
  lapply(c("png", "svg", "pdf"), 
         function(fmt) ggsave(paste0("./Figures/HEAP/F1/","Gvar.", fmt), 
                              plot = p2, width = 3, height = 3, dpi = 500))
  
  lapply(c("png", "svg", "pdf"), 
         function(fmt) ggsave(paste0("./Figures/HEAP/F1/","EvsGvar.", fmt), 
                              plot = p3, width = 4, height = 3, dpi = 500))
}
CovarSpec_varplot(R2test)



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
                              plot = p1, width = 6, height = 4, dpi = 500))
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
                              plot = p1, width = 4, height = 3, dpi = 500))
  
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
                              plot = p1, width = 1, height = 4, dpi = 500))
  
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
                            plot = cow1, width = 8, height = 6, dpi = 500))


#'*FIGURE OUT BETTER VISUAL FOR RESULTS w/ and w/o condition on BMI*


#### CATEGORY PLOTS ####
plot_testR2Cat_CI <- R2testCat_CI %>%
                 select(-c("cID","Alcohol","Diet_Weekly","Smoking","Exercise_MET",
                           "Exercise_Freq","Internet_Usage","Deprivation_Indices","Vitamins")) %>%
                  unique()
  

plot_testR2_all <- merge(plot_testR2_CI, plot_testR2Cat_CI, by = "ID")


library(cowplot)
library(ggpubr)
category_specificR2plot <- function(plot_testR2_all, category, category_name, G, E, y_cutoff){
  subset_test <- plot_testR2_all  %>% 
    slice_max(.data[[category]], n = 10) 
  
  
  p1 <- plot_testR2_all %>% 
    ggplot(aes(x = .data[[G]], y = .data[[category]], label = ID)) + 
    geom_point(aes(color = .data[[E]])) + #aes(color = .data[[category]])) +
    geom_abline(slope=1, intercept = 0) +
    geom_ribbon(aes(ymin = 0, ymax = ifelse(.data[[G]] <= y_cutoff, .data[[G]], y_cutoff), fill = "Genetics Driven"), alpha = 0.2) +
    geom_ribbon(aes(ymin = ifelse(.data[[G]] <= y_cutoff, .data[[G]], y_cutoff), ymax = y_cutoff, fill = "Environment Driven"), alpha = 0.2) +
    scale_fill_manual(values = c("Genetics Driven" = "blue", "Environment Driven" = "green")) +
    ylim(c(0, y_cutoff)) +
    geom_text_repel(data=subset_test, aes(x = .data[[G]], y = .data[[category]], label= ID),
                    size = 4) +
    #geom_label() + 
    scale_color_gradient(low = "darkslategrey", high = "red") +
    xlab(expression("G: PGS"~R^2)) +
    ylab(bquote(~.(category_name)~":"~R^2)) +
    labs(color="E: PXS"~R^2) +
    theme_minimal() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) + #Adjust 90 to 0 for this plot:
    
    #Remove: the Environment vs Genetics Driven Legend for spacing:
    guides(fill = "none",  # Remove the fill legend
           color ="none"#guide_colorbar(title = bquote("PXS:" ~ R^2))  # Keep the color gradient legend
    ) +
    #ggtitle(bquote("Test Set:"~.(category_name))) +
    ggtitle(category_name) +
    theme(plot.title = element_text(hjust = 0.5, size = 14),
          legend.position = "bottom",
          legend.title = element_text(angle = 0),  # Optional: Set title orientation
          legend.text = element_text(angle = 90),  # Rotate legend text to vertical
          axis.text = element_text(size = 10),
          axis.title = element_text(size = 12)) #"right")
  p1
  return(p1)
}
get_legend <- function(plot_testR2_all, category, category_name, G, E, y_cutoff){
  p2 <- plot_testR2_all %>% 
    ggplot(aes(x = .data[[G]], y = .data[[category]], label = ID)) + 
    geom_point(aes(color = .data[[E]])) + #aes(color = .data[[category]])) +
    scale_color_gradient(low = "darkslategrey", high = "red", breaks=c(0,0.025,0.05,0.075,0.1,0.125,0.15)) +
    xlab(expression("G: PGS"~R^2)) +
    ylab(bquote(~.(category_name)~":"~R^2)) +
    labs(color="E: PXS"~R^2) +
    theme_minimal()
  
  library(ggpubr)
  leg <- get_legend(p2)
  p3 <- as_ggplot(leg)
  return(p3)
}

g1 <- category_specificR2plot(plot_testR2_all, category = "mean_Exercise_Freq",
                              category_name = "Exercise_Freq",
                              G = "G_stats_mean",E = "E_stats_mean", y_cutoff = 0.1)

g2 <- category_specificR2plot(plot_testR2_all, category = "mean_Exercise_MET",
                              category_name = "Exercise_MET",
                              G = "G_stats_mean",E = "E_stats_mean", y_cutoff = 0.06)

g3 <- category_specificR2plot(plot_testR2_all, category = "mean_Diet_Weekly",
                              category_name = "Diet_Weekly",
                              G = "G_stats_mean",E = "E_stats_mean", y_cutoff = 0.06)

g4 <- category_specificR2plot(plot_testR2_all, category = "mean_Alcohol",
                              category_name = "Alcohol",
                              G = "G_stats_mean",E = "E_stats_mean", y_cutoff = 0.06)

g5 <- category_specificR2plot(plot_testR2_all, category = "mean_Smoking",
                              category_name = "Smoking",
                              G = "G_stats_mean",E = "E_stats_mean", y_cutoff = 0.045)

g6 <- category_specificR2plot(plot_testR2_all, category = "mean_Vitamins",
                              category_name = "Vitamins",
                              G = "G_stats_mean",E = "E_stats_mean", y_cutoff = 0.025)

g7 <- category_specificR2plot(plot_testR2_all, category = "mean_Internet_Usage",
                              category_name = "Internet_Usage",
                              G = "G_stats_mean",E = "E_stats_mean", y_cutoff = 0.02)

g8 <- category_specificR2plot(plot_testR2_all, category = "mean_Deprivation_Indices",
                              category_name = "Deprivation_Indices",
                              G = "G_stats_mean",E = "E_stats_mean", y_cutoff = 0.02)



library(patchwork)
combined_plot <- (g1 | g2 | g3 | g4) /
  (g5 | g6 | g7 | g8)
#combined_plot
ggsave("./Figures/HEAP/A1/meanR2_Categories.png", 
       plot = combined_plot, width = 15, height = 8, units = "in")
ggsave("./Figures/HEAP/A1/meanR2_Categories.svg", 
       plot = combined_plot, width = 15, height = 8, units = "in")
ggsave("./Figures/HEAP/A1/meanR2_Categories.pdf", 
       plot = combined_plot, width = 15, height = 8, units = "in")


####Alternative HEATMAP:
library(ComplexHeatmap)


catNames <- c(Alcohol = "mean_Alcohol",
              Diet_Weekly = "mean_Diet_Weekly",
              Smoking = "mean_Smoking",
              Exercise_MET = "mean_Exercise_MET",
              Exercise_Freq = "mean_Exercise_Freq",
              Internet_Usage = "mean_Internet_Usage",
              Deprivation_Indices = "mean_Deprivation_Indices",
              Vitamins = "mean_Vitamins")

testR2_heatmap <- plot_testR2Cat_CI %>%
                    select(c("ID","mean_Alcohol","mean_Diet_Weekly",
                             "mean_Smoking","mean_Exercise_MET",
                             "mean_Exercise_Freq","mean_Internet_Usage",
                             "mean_Deprivation_Indices","mean_Vitamins")) %>%
                    rename(all_of(catNames))

#Top 70 proteins by R2:
topEprot <- plot_testR2_CI %>% 
                arrange(desc(E_stats_mean)) %>%
                head(70) %>%
                pull(ID)

# filter(E_stats_mean > 0.055) %>%
# pull(ID)

#Select top 10 proteins for each category!
categories <- c("Alcohol","Diet_Weekly",
                "Smoking","Exercise_MET",
                "Exercise_Freq","Internet_Usage",
                "Deprivation_Indices","Vitamins")

protSelect <- lapply(categories, function(x){
  testR2_heatmap <- testR2_heatmap[testR2_heatmap[[x]] > 0.02,]
  as.character(testR2_heatmap[order(testR2_heatmap[[x]],
                                    decreasing = T),]$ID[1:30])
})
protSelect <- unique(unlist(protSelect))


testR2_heatmap <- testR2_heatmap %>%
                      filter(ID %in% topEprot) %>%
                      mutate(ID = factor(ID, levels = topEprot)) %>%
                      arrange(ID) %>%
                      column_to_rownames(var = "ID")

#Combine Exercise_MET and Exercise_Freq to just Exercise:
testR2_heatmap <- testR2_heatmap %>%
                    mutate(Exercise = rowSums(across(c("Exercise_Freq","Exercise_MET"))))


library(circlize)


#Cluster Row and Columns
ht <- Heatmap(as.matrix(testR2_heatmap %>%
                          select(-c("Internet_Usage",
                                    "Deprivation_Indices",
                                    "Exercise_Freq",
                                    "Exercise_MET")) %>%
                rename(Diet = "Diet_Weekly")),
              name = "R²",
              col = colorRamp2(c(0, 0.05, 0.1), c("blue", "white", "red")),
              #c("#440154", "#21908C", "#FDE725")
              width = ncol(testR2_heatmap)*unit(4, "mm"),
              height = nrow(testR2_heatmap)*unit(2.1, "mm"),
              show_row_names = TRUE,
              show_column_names = TRUE,
              cluster_rows = TRUE,
              cluster_columns = TRUE,
              show_column_dend = FALSE, 
              show_row_dend = FALSE,
              #row_split = 5,
              #row_km = 3,
              #column_gap = unit(8, "mm"),
              #row_gap = unit(20, "mm"),
              row_names_gp = gpar(fontsize = 7),
              column_names_gp = gpar(fontsize = 10),
              row_names_side = "left",
              heatmap_legend_param = list(
                title_gp = gpar(fontsize = 10), # Legend title font size
                labels_gp = gpar(fontsize = 10), # Legend labels font size
                legend_height = unit(2, "cm"), # Height of the legend
                legend_width = unit(1, "cm")   # Width of the legend
              )
)

draw(ht)

filename <- paste0("./Figures/HEAP/F1/Ecat_heatmap.png")
png(file=filename,
    height = 7, width = 4, units = "in", res = 500)
draw(ht)
dev.off()

filename <- paste0("./Figures/HEAP/F1/Ecat_heatmap.svg")
svg(file=filename,
    height = 7, width = 4)
draw(ht)
dev.off()



#Cluster only Columns
ht <- Heatmap(as.matrix(testR2_heatmap %>%
                          select(-c("Internet_Usage",
                                    "Deprivation_Indices",
                                    "Exercise_Freq",
                                    "Exercise_MET")) %>%
                          rename(Diet = "Diet_Weekly")),
              name = "R²",
              col = colorRamp2(c(0, 0.05, 0.1), c("blue", "white", "red")),
              #c("#440154", "#21908C", "#FDE725")
              width = ncol(testR2_heatmap)*unit(4, "mm"),
              height = nrow(testR2_heatmap)*unit(2.1, "mm"),
              show_row_names = TRUE,
              show_column_names = TRUE,
              cluster_rows = FALSE,
              cluster_columns = TRUE,
              show_column_dend = FALSE, 
              show_row_dend = FALSE,
              #row_split = 5,
              #row_km = 3,
              #column_gap = unit(8, "mm"),
              #row_gap = unit(20, "mm"),
              row_names_gp = gpar(fontsize = 7),
              column_names_gp = gpar(fontsize = 10),
              row_names_side = "left",
              heatmap_legend_param = list(
                title_gp = gpar(fontsize = 10), # Legend title font size
                labels_gp = gpar(fontsize = 10), # Legend labels font size
                legend_height = unit(2, "cm"), # Height of the legend
                legend_width = unit(1, "cm")   # Width of the legend
              )
)

draw(ht)





###EXTRA::::
ht <- Heatmap(as.matrix(testR2_heatmap),
              name = "R2",
              col = colorRamp2(c(0, 0.04, 0.08), c("purple", "white", "orange")),
              width = ncol(testR2_heatmap)*unit(10, "mm"),
              height = nrow(testR2_heatmap)*unit(3, "mm"),
              show_row_names = TRUE,
              show_column_names = TRUE,
              cluster_rows = TRUE,
              cluster_columns = TRUE,
              show_column_dend = FALSE, 
              show_row_dend = FALSE,
              column_gap = unit(8, "mm"),
              row_names_gp = gpar(fontsize = 8),
              column_names_gp = gpar(fontsize = 8),
              row_names_side = "left",
              cell_fun = function(j, i, x, y, width, height, fill) {
                if (testR2_heatmap[i, j] >= 0.02) {
                  grid.text(sprintf("%.2f", testR2_heatmap[i, j]), x = x, y = y, 
                            gp = gpar(fontsize = 6, col = "black"))
                }
              },
              border = T,
              rect_gp = gpar(col = "black", lwd = 0.5)
)

draw(ht)

#Remove Internet Usage + Deprivation Indices:
custom_colors <- colorRamp2(c(0, 0.04, 0.08), c("#440154", "#21908C", "#FDE725"))

library(ComplexHeatmap)
library(circlize)
library(grid)

# Assuming `testR2_heatmap` is your data matrix

ht <- Heatmap(
  as.matrix(testR2_heatmap),
  name = "R²",
  col = colorRamp2(c(0, 0.04, 0.08), c("blue", "white", "red")),  # Improved color gradient
  width = ncol(testR2_heatmap) * unit(10, "mm"),
  height = nrow(testR2_heatmap) * unit(5, "mm"),  # Adjusted height for better proportions
  show_row_names = TRUE,
  show_column_names = TRUE,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_column_dend = TRUE,  # Show dendrogram for context
  show_row_dend = TRUE,
  column_gap = unit(8, "mm"),
  row_names_gp = gpar(fontsize = 10, fontface = "bold"),  # Improved row label visibility
  column_names_gp = gpar(fontsize = 10, fontface = "bold"),  # Improved column label visibility
  row_names_side = "left",
  border = TRUE,
  rect_gp = gpar(col = "black", lwd = 0.5),
  cell_fun = function(j, i, x, y, width, height, fill) {
    value <- testR2_heatmap[i, j]
    if (value > 0.04) {
      grid.text(sprintf("%.2f", value), x = x, y = y, 
                gp = gpar(fontsize = 8, col = "black", fontface = "bold"))  # Text only for values > 0.04
    }
  },
  heatmap_legend_param = list(
    title = "R²",
    at = c(0, 0.04, 0.08),
    labels = c("Low", "Threshold", "High"),
    title_gp = gpar(fontsize = 12, fontface = "bold"),
    labels_gp = gpar(fontsize = 10)
  )
)

# Draw heatmap
draw(ht, heatmap_legend_side = "right")


# FUNCTION TO GO ACROSS ALL CATEGORIES AND PLOT E v G plot:
# FUNCTION TO GET LEGEND:
# & PROVIDE PLOTS FOR VARIABILITY of EACH CATEGORY.




#'*PLOT 5* Make a figure that shows these effects independent of the GvE plot



### Overall Average Plot across E: Categories
library(cowplot)
category_specificR2plot <- function(category, catR2, ycut){
  R2testCat_all <- merge(R2test, R2testCat, by = "ID")
  y_cutoff = ycut
  
  
  #Subset to label:
  subset_test <- R2testCat_all  %>% 
    filter(!ID %in% c(badfit_proteins)) %>%
    mutate_all(~replace_na(., 0)) %>%
    #filter(.data[[category]] > G) %>%
    filter(.data[[category]] > catR2) %>%
    top_n(10, .data[[category]]) 
  
  #Color by specific categories.
  p1 <- R2testCat_all %>% 
    filter(!ID %in% c(badfit_proteins)) %>%
    mutate_all(~replace_na(., 0)) %>% 
    ggplot(aes(x = G, y = .data[[category]], label = ID)) + 
    geom_point(aes(color = E)) + #aes(color = .data[[category]])) +
    geom_abline(slope=1, intercept = 0) +
    geom_ribbon(aes(ymin = 0, ymax = ifelse(G <= y_cutoff, G, y_cutoff), fill = "Genetics Driven"), alpha = 0.2) +
    geom_ribbon(aes(ymin = ifelse(G <= y_cutoff, G, y_cutoff), ymax = y_cutoff, fill = "Environment Driven"), alpha = 0.2) +
    scale_fill_manual(values = c("Genetics Driven" = "blue", "Environment Driven" = "green")) +
    ylim(c(0, y_cutoff)) +
    geom_text_repel(data=subset_test, aes(x = G,y = .data[[category]],label= ID),
                    size = 4) +
    #geom_label() + 
    scale_color_gradient(low = "darkslategrey", high = "red") +
    xlab(expression("G: PGS"~R^2)) +
    ylab(bquote(~.(category)~":"~R^2)) +
    theme_minimal() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) + #Adjust 90 to 0 for this plot:
    
    #Remove: the Environment vs Genetics Driven Legend for spacing:
    guides(fill = "none",  # Remove the fill legend
           color ="none"#guide_colorbar(title = bquote("PXS:" ~ R^2))  # Keep the color gradient legend
    ) +
    ggtitle(bquote("Test Set:"~.(category))) +
    theme(plot.title = element_text(hjust = 0.5, size = 14),
          legend.position = "bottom",
          legend.title = element_text(angle = 0),  # Optional: Set title orientation
          legend.text = element_text(angle = 90),  # Rotate legend text to vertical
          axis.text = element_text(size = 10),
          axis.title = element_text(size = 12)) #"right")
  
  return(p1)
}
get_legend <- function(){
  
}

g1 <- category_specificR2plot(category = "Exercise_Freq",catR2 = 0.02, ycut = 0.1)
g2 <- category_specificR2plot(category = "Exercise_MET",catR2 = 0.02, ycut = 0.05)
g3 <- category_specificR2plot(category = "Diet_Weekly",catR2 = 0.02, ycut = 0.075)
g4 <- category_specificR2plot(category = "Vitamins",catR2 = 0.01, ycut = 0.03)
g5 <- category_specificR2plot(category = "Alcohol",catR2 = 0.02, ycut = 0.05)
g6 <- category_specificR2plot(category = "Smoking",catR2 = 0.02, ycut = 0.05)
g7 <- category_specificR2plot(category = "Deprivation_Indices",catR2 = 0.01, ycut = 0.02)
g8 <- category_specificR2plot(category = "Internet_Usage",catR2 = 0.01, ycut = 0.02)
library(patchwork)
combined_plot <- (g1 | g2 | g3 | g4) /
  (g5 | g6 | g7 | g8)



#'*PLOT 4* Check the stability of R2 estimates across train and test:
partR2cat <- function(category){
  train <- R2trainCat %>% 
    rename(R2train = all_of(category)) %>%
    select(all_of(c("ID", "R2train")))
  test <- R2testCat %>% 
    rename(R2test = all_of(category)) %>%
    select(all_of(c("ID", "R2test")))
  R2cat <- merge(train, test, by = "ID")
  
  p1 <- R2cat %>%
    ggplot(aes(x = R2train, y = R2test)) +
    geom_point() +
    #geom_abline(slope=1, intercept = 0) +
    stat_poly_line() +
    stat_poly_eq() +
    theme_minimal() +
    xlab(bquote("Train" ~ R^2)) +
    ylab(bquote("Test" ~ R^2)) +
    ggtitle(bquote("Partitioned"~R^2~"-"~.(category)))
  
  return(p1)
}

g1 <- partR2cat(category = "Exercise_Freq")
g2 <- partR2cat(category = "Exercise_MET")
g3 <- partR2cat(category = "Diet_Weekly")
g4 <- partR2cat(category = "Vitamins")
g5 <- partR2cat(category = "Alcohol")
g6 <- partR2cat(category = "Smoking")
g7 <- partR2cat(category = "Deprivation_Indices")
g8 <- partR2cat(category = "Internet_Usage")

library(patchwork)
combined_plot <- (g1 | g2 | g3 | g4) /
  (g5 | g6 | g7 | g8)
# Save the plot to a file
ggsave("/n/groups/patel/shakson_ukb/UK_Biobank/Figures/protPGS_PXS_summ/R2_Categories.png", 
       plot = combined_plot, width = 15, height = 8, units = "in")

#'*^^REORDER the figures and categories later.*
#'*MAKE sure if you reorder the plot above to fix it for the plot here.*




#### BELOW IS TRASH BIN #####
#### HELPFUL Code for Plot Generation ####







#Get variability of E and G...
library(tidyr)
library(broom)

Lasso <- PXS_type5$Lasso
R2train <- PXS_type5$R2train
R2test <- PXS_type5$R2test
R2trainCat <- PXS_type5$R2trainCat
R2testCat <- PXS_type5$R2testCat


# Define a function to compute summary statistics
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

Lasso_CI <- Lasso %>%
  group_by(omic) %>%
  summarize(across(everything(), get_CI)) #used summarize previously

R2train_CI <- R2train %>%
  group_by(ID) %>%
  summarize(across(everything(), get_CI))

R2test_CI <- R2test %>%
  group_by(ID) %>%
  summarize(across(everything(), get_CI))

R2trainCat_CI <- R2trainCat %>%
  group_by(ID) %>%
  summarize(across(everything(), get_CI))

R2testCat_CI <- R2testCat %>%
  group_by(ID) %>%
  summarize(across(everything(), get_CI))


# Function to split string and calculate numeric values
process_column <- function(column) {
  tibble(
    mean = as.numeric(str_extract(column, "^[^ ±]+")),
    ci_margin = as.numeric(str_extract(column, "(?<=± )[^ ]+$")),
    ci_low = mean - ci_margin,
    ci_high = mean + ci_margin
  )
}

df <- R2test_CI %>%
  # Apply the processing function to each column and store the results
  rowwise() %>%
  mutate(
    G_stats = list(process_column(G)),
    E_stats = list(process_column(E)),
    GxE_stats = list(process_column(GxE))
  ) %>%
  # Unnest the statistics into separate columns
  unnest(c(G_stats, E_stats, GxE_stats), names_sep = "_") 
#%>%
# Remove NA rows (if needed)
#drop_na() %>%
# Group by ID and summarize statistics
#group_by(ID) %>%





#'*Plotting example:*
library(ggrepel)
df %>%
  ggplot(aes(x = G_stats_mean, y = E_stats_mean, label = ID)) +
  geom_point() +
  geom_text_repel() +
  geom_pointrange(aes(ymin = E_stats_ci_low, ymax = E_stats_ci_high)) +
  geom_pointrange(aes(xmin = G_stats_ci_low, xmax = G_stats_ci_high)) +
  geom_abline(slope=1, intercept = 0) +
  labs(
    x = "G",
    y = "E",
    title = "E-driven Proteins Error bars"
  ) +
  theme_minimal()

plot_testR2_CI <- df

#Plot for entire E vs G:
y_cutoff = 0.2
subset_test <- plot_testR2_CI %>% 
  filter(!ID %in% c(badfit_proteins)) %>%
  mutate_all(~replace_na(., 0)) %>%
  filter((E_stats_mean > G_stats_mean) & (E_stats_mean > 0.07))


p1 <- plot_testR2_CI %>% 
  filter(!ID %in% c(badfit_proteins)) %>%
  mutate_all(~replace_na(., 0)) %>% 
  ggplot(aes(x = G_stats_mean, y = E_stats_mean, label = ID)) + 
  geom_point() +
  geom_abline(slope=1, intercept = 0) +
  geom_ribbon(aes(ymin = 0, ymax = ifelse(G_stats_mean <= y_cutoff, G_stats_mean, y_cutoff), fill = "Genetics Driven"), alpha = 0.2) +
  geom_ribbon(aes(ymin = ifelse(G_stats_mean <= y_cutoff, G_stats_mean, y_cutoff), ymax = y_cutoff, fill = "Environment Driven"), alpha = 0.2) +
  scale_fill_manual(values = c("Genetics Driven" = "blue", "Environment Driven" = "green")) +
  ylim(c(0, y_cutoff)) +
  geom_text_repel(data=subset_test, aes(x = G_stats_mean,y = E_stats_mean,label= ID),
                  size = 3.5) +
  geom_pointrange(aes(ymin = E_stats_ci_low, ymax = E_stats_ci_high)) +
  geom_pointrange(aes(xmin = G_stats_ci_low, xmax = G_stats_ci_high)) +
  #geom_label() + 
  xlab(expression("G: PGS "~R^2)) +
  ylab(expression("E: PXS "~R^2)) +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  guides(fill=guide_legend(title="Region")) +
  ggtitle(expression("Test Set: Partitioned"~R^2)) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "bottom") 
p1
















#'*Specific Categories*
####Exercise Category:
df2 <- R2testCat_CI %>%
  # Apply the processing function to each column and store the results
  rowwise() %>%
  mutate(
    Alcohol_stats = list(process_column(Alcohol)),
    Diet_Weekly_stats = list(process_column(Diet_Weekly)),
    Smoking_stats = list(process_column(Smoking)),
    Exercise_MET_stats = list(process_column(Exercise_MET)),
    Exercise_Freq_stats = list(process_column(Exercise_Freq)),
    Internet_Usage_stats = list(process_column(Internet_Usage)),
    Deprivation_Indices_stats = list(process_column(Deprivation_Indices)),
    Vitamins_stats = list(process_column(Vitamins))
  ) %>%
  # Unnest the statistics into separate columns
  unnest(c(Alcohol_stats, Diet_Weekly_stats,
           Smoking_stats, Exercise_MET_stats,
           Exercise_Freq_stats, Internet_Usage_stats,
           Deprivation_Indices_stats, Vitamins_stats), names_sep = "_") 

df3 <- merge(df,df2,by="ID")
df3 %>%
  ggplot(aes(x = G_stats_mean, y = Exercise_Freq_stats_mean, label = ID)) +
  geom_point() +
  geom_text_repel() +
  geom_pointrange(aes(ymin = Exercise_Freq_stats_ci_low, 
                      ymax = Exercise_Freq_stats_ci_high)) +
  geom_pointrange(aes(xmin = G_stats_ci_low, xmax = G_stats_ci_high)) +
  geom_abline(slope=1, intercept = 0) +
  labs(
    x = "G",
    y = "E",
    title = "E-driven Proteins Error bars"
  ) +
  theme_minimal()


df3 %>%
  ggplot(aes(x = G_stats_mean, y = Diet_Weekly_stats_mean, label = ID)) +
  geom_point() +
  geom_text_repel() +
  geom_pointrange(aes(ymin = Diet_Weekly_stats_ci_low, 
                      ymax = Diet_Weekly_stats_ci_high)) +
  geom_pointrange(aes(xmin = G_stats_ci_low, xmax = G_stats_ci_high)) +
  geom_abline(slope=1, intercept = 0) +
  labs(
    x = "G",
    y = "E",
    title = "E-driven Proteins Error bars"
  ) +
  theme_minimal()


df3 %>%
  ggplot(aes(x = G_stats_mean, y = Smoking_stats_mean, label = ID)) +
  geom_point() +
  geom_text_repel() +
  geom_pointrange(aes(ymin = Smoking_stats_ci_low, 
                      ymax = Smoking_stats_ci_high)) +
  geom_pointrange(aes(xmin = G_stats_ci_low, xmax = G_stats_ci_high)) +
  geom_abline(slope=1, intercept = 0) +
  labs(
    x = "G",
    y = "E",
    title = "E-driven Proteins Error bars"
  ) +
  theme_minimal()


df3 %>%
  ggplot(aes(x = G_stats_mean, y = Alcohol_stats_mean, label = ID)) +
  geom_point() +
  geom_text_repel() +
  geom_pointrange(aes(ymin = Alcohol_stats_ci_low, 
                      ymax = Alcohol_stats_ci_high)) +
  geom_pointrange(aes(xmin = G_stats_ci_low, xmax = G_stats_ci_high)) +
  geom_abline(slope=1, intercept = 0) +
  labs(
    x = "G",
    y = "E",
    title = "E-driven Proteins Error bars"
  ) +
  theme_minimal()


#%>%
# Remove NA rows (if needed)
#drop_na() %>%
# Group by ID and summarize statistics
#group_by(ID) %>%




#'*Plotting example:*
library(ggrepel)
df %>%
  ggplot(aes(x = G_stats_mean, y = E_stats_mean, label = ID)) +
  geom_point() +
  geom_text_repel() +
  geom_pointrange(aes(ymin = E_stats_ci_low, ymax = E_stats_ci_high)) +
  geom_pointrange(aes(xmin = G_stats_ci_low, xmax = G_stats_ci_high)) +
  geom_abline(slope=1, intercept = 0) +
  labs(
    x = "G",
    y = "E",
    title = "E-driven Proteins Error bars"
  ) +
  theme_minimal()



#### NOT USED: HELPFUL CODE INDIVIDUAL CovarSpec Visualization w/ Nested CV Variability (Organize Later) ####

#Just look at R2 test: G, E, GxE

R2test <- PXS_type5$R2test

# Define a function to compute summary statistics
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
R2test_CI <- R2test %>%
  select(all_of(c("ID","G","E","GxE"))) %>%
  group_by(ID) %>%
  summarize(across(everything(), get_CI))

# Function to split string and calculate numeric values
process_column <- function(column) {
  tibble(
    mean = as.numeric(str_extract(column, "^[^ ±]+")),
    ci_margin = as.numeric(str_extract(column, "(?<=± )[^ ]+$")),
    ci_low = mean - ci_margin,
    ci_high = mean + ci_margin
  )
}

df <- R2test_CI %>%
  # Apply the processing function to each column and store the results
  rowwise() %>%
  mutate(
    G_stats = list(process_column(G)),
    E_stats = list(process_column(E)),
    GxE_stats = list(process_column(GxE))
  ) %>%
  # Unnest the statistics into separate columns
  unnest(c(G_stats, E_stats, GxE_stats), names_sep = "_") 




df %>%
  ggplot(aes(x = G_stats_mean, y = E_stats_mean, label = ID)) +
  geom_point() +
  geom_text_repel() +
  geom_pointrange(aes(ymin = E_stats_ci_low, ymax = E_stats_ci_high)) +
  geom_pointrange(aes(xmin = G_stats_ci_low, xmax = G_stats_ci_high)) +
  geom_abline(slope=1, intercept = 0) +
  labs(
    x = "G",
    y = "E",
    title = "E-driven Proteins Error bars"
  ) +
  theme_minimal()

#PLOT of R2 and variance...
#'*Plot of R2 and variance*
plot(df$E_stats_mean, df$E_stats_ci_margin)
plot(df$GxE_stats_mean, df$GxE_stats_ci_margin)


#Aggregate R2's
plot(df$E_stats_mean - df$G_stats_mean, df$E_stats_ci_margin)
plot(df$E_stats_mean - df$G_stats_mean, df$G_stats_ci_margin)


#'*Finalized GvE plot across multiple covariate specification*
#NEED to also do these GvE plots for each covariate specification separately.
