#Libraries:
library(data.table)
library(tidyverse)
library(ggplot2)
library(ggpmisc)

#For Lasso:
library(glmnet)
library(car)
set.seed(123)

# Definition of PXSconstruct to make modifications to covars_list
PXSconstruct <- setClass(
  "PXSconstruct",
  slots = c(
    Elist = "list", #E dataframes separated by category:
    Elist_names = "character", #Names of the category:
    Eid_cat = "data.frame", #Dataframe with environmental IDs and corresponding category.
    ordinalIDs = "character", #list environmental ordinal names
    
    UKBprot_df = "data.frame", #Protein variables dataframe
    protIDs = "character", #protein variables
    
    covars_df = "data.frame", #Covariates Dataframe
    covars_list = "character" #covars variables
  )
)


#'*Loading:*
PXSloader <- readRDS(file = "/n/scratch/users/s/shi872/UKB_intermediate/UKB_PGS_PXS_load.rds")


#'*Functions: Loadings Genetic Score*
#'Extract relevant PGS score from files:
#'NEED TO OUTPUT WARNING INSTEAD OF ERROR FOR ANY FILES MISSING!!!
#How to extract files with output of .sscore
extract_protGS <- function(protID){
  omicpredIDs <- fread(file = "/n/groups/patel/shakson_ukb/UK_Biobank/Data/OMICSPRED/UKB_Olink_multi_ancestry_models_val_results_portal.csv")
  dir <- "/n/scratch/users/s/shi872/UKB_intermediate/Genetics/UKBProtGS/"
  opID <- omicpredIDs$OMICSPRED_ID[match(protID,omicpredIDs$Gene)]
  
  opFile <- paste0(opID, ".sscore")
  protGS <- fread(file = paste0(dir, opFile))
  colnames(protGS) <- c("eid", paste0(protID,"_GS"))
  
  #Convert dash to underscore:
  colnames(protGS) <- gsub("-", "_", colnames(protGS))
  return(protGS)
}


#'I had to modify the GS_struct function + more to handle dashes and replace with underscores!
GS_struct <- function(protID, UKBprot_df){
  omic_PGS <- extract_protGS(protID)
  
  #Convert dash to underscore: for HLA genes:
  protID <- gsub("-","_",protID)
  omic_orig <- UKBprot_df %>% select(all_of(c("eid", protID)))
  
  omicGS <- list()
  omicGS[["combo"]] <- na.omit(merge(omic_PGS, omic_orig, by = "eid"))
  omicGS[["solo"]] <- na.omit(omic_PGS)
  
  return(omicGS)
}


#'*Functions: Preprocessing*
continuous_finder <- function(df){
  max_cols <- apply(df, 2, max, na.rm = TRUE)
  min_cols <- apply(df, 2, min, na.rm = TRUE) 
  unique_vals <- sapply(df, function(x) length(unique(x[!is.na(x)])))  # Number of unique values
  #^unique_vals handles cases of continuous features that range from 0 to 1 but are not one-hot encoded.
  # Identify continuous columns: not one-hot encoded (0 and 1) or values greater than > 5 as max
  continuous <- names(max_cols[max_cols > 5 | unique_vals > 2])
  
  return(continuous)
}

# Remove columns with too much missing data
remove_missing_cols <- function(df, miss_rate = 0.2){
  #NAcols <- colSums(is.na(df))
  NAcols <- colMeans(is.na(df))
  
  removecols <- names(NAcols[NAcols > miss_rate]) #look here
  df2 <- df %>% select(!all_of(c(removecols)))
  return(df2)
}

remove_cols_name <- function(df, miss_rate = 0.2){
  #NAcols <- colSums(is.na(df))
  NAcols <- colMeans(is.na(df))
  
  removecols <- names(NAcols[NAcols > miss_rate]) #look here
  return(removecols)
}

# Function to Deal with Categorical Variables
categorical_handler <- function(df, ordinal_names, ordinal_contrast = "treatment"){
  
  # :: Categorical Variables ::
  #Make some of the data categorical to help with building code:
  #'*Make ordinal factors and MAKE sure these levels remain consistent in train and test*
  for(i in ordinal_names){
    df[[i]] <- factor(df[[i]], ordered = T)
    
    # Set treatment contrasts: instead of polynomial contrasts:
    #if (length(combined_levels) >= 2) {
    
    if(ordinal_contrast == "treatment"){
      contrasts(df[[i]]) <- contr.treatment(length(levels(df[[i]])))
    } else if (ordinal_contrast == "sum"){
      contrasts(df[[i]]) <- contr.sum(length(levels(df[[i]])))
    } 
  }
  
  return(df)
}

# Function to Build Formula
build_formula <- function(protID, covars_list, E_ids){
  protID <- gsub("-", "_", protID)
  G_var <- paste0(protID,"_GS")
  
  covariates <- covars_list
  
  #For loop to build formula: E, GxE components
  pred_E <- c()
  for(i in E_ids){
    #Setup Formula:
    pred_Ei <- c(i, paste(G_var, "*", i))
    pred_E <- c(pred_E, pred_Ei)
  }
  
  #Final formula with G, E, and GxE components with covariates:
  pred_var <- c(G_var, pred_E, covariates,"0")
  formula <- as.formula(paste(c(protID), paste(pred_var, collapse="+"), sep="~"))
  
  return(formula)
}

# Function to Preprocess
preprocess_data <- function(protID, covariates,
                            trainData, testData, ordinal_contrast = "treatment",
                            ordinal_names){
  # Setup protID and Genetic Score ID:
  protID <- gsub("-", "_", protID)
  G_var <- paste0(protID,"_GS")
  
  
  #na.omit: make complete data
  trainData <- na.omit(trainData)
  testData <- na.omit(testData)
  
  trainData <- as.data.frame(trainData)
  testData <- as.data.frame(testData)
  
  # :: Continuous Variables ::
  
  #Scale only the omics, genetic score, and covariates:
  rel_columns <- c(protID, G_var, covariates)
  numeric_cols <- sapply(trainData[, rel_columns], 
                         function(col) is.numeric(col) && !(all(col %in% c(0, 1)))) #Only continuous feature no binaries
  DG_cont_names <- names(numeric_cols[numeric_cols])
  
  #'*KEY TO ADD in runProt_PGS_ver2*
  #'Scale other environments with continuous features:
  #Only obtain Env variables (NOT ordinal, NOT binary, NOT covariates, NOT genetic score or Omic)
  df <- trainData %>% select(!all_of(c("eid", rel_columns, ordinal_names)))
  E_cont_names <- continuous_finder(df)
  
  #Combine names to Scale:
  numeric_col_names <- c(DG_cont_names, E_cont_names)
  
  
  #Scale numeric columns of continuous features::
  mean_train <- lapply(trainData[numeric_col_names], function(x) mean(x,na.rm = T))
  sd_train <- lapply(trainData[numeric_col_names], function(x) sd(x,na.rm = T))
  
  for(i in numeric_col_names){
    #Z-score train Dataset
    trainData[[i]] <- (trainData[[i]] - mean_train[[i]]) / sd_train[[i]]
    
    #Z-score test Dataset using same parameters as train:
    testData[[i]] <- (testData[[i]] - mean_train[[i]]) / sd_train[[i]]
  }
  
  Data <- list()
  Data[["train"]] <- trainData
  Data[["test"]] <- testData
  
  return(Data)
}


#'*Functions for Lasso Fit, Score Generation and Metrics*

# Function to extract contrasts for all factor columns in a data frame
extract_contrasts <- function(data) {
  contrasts_list <- list()
  for (colname in names(data)) {
    if (is.factor(data[[colname]])) {
      contrasts_list[[colname]] <- contrasts(data[[colname]], contrasts = TRUE)
    }
  }
  return(contrasts_list)
}

# Function to fit lasso and choose the best hyperparameter lamba
# Default 10 fold cross validation on training dataset
lasso_fit <- function(protID, trainData, formula){
  #Account for protID different naming schemes:
  protID <- gsub("-", "_", protID)
  
  # Obtain Independent and Dependent variables
  contrast_list <- extract_contrasts(trainData) #'*added to prevent Errors*
  
  x <- model.matrix(formula, data = trainData, contrasts.arg = contrast_list)
  y <- trainData[[protID]]
  
  # Fit Lasso model with cross-validation
  cv_lasso <- cv.glmnet(x, y, alpha = 1) #alpha = 1 is lasso
  best_lambda <- cv_lasso$lambda.min #'**
  
  # Fit the final Lasso model with the best lambda
  lasso_model <- glmnet(x, y, alpha = 1, lambda = best_lambda)
  
  # Extract coefficients from the Lasso model
  lasso_coef <- coef(lasso_model)
  non_zero_coef <- which(lasso_coef != 0)[-1]  # Exclude the intercept
  
  # Names of features with non-zero coefficients
  selected_features <- rownames(lasso_coef)[non_zero_coef]
  
  # Coefficients of the selected features
  selected_coefficients <- lasso_coef[non_zero_coef]
  
  fit <- list()
  fit[["lasso"]] <- lasso_model
  fit[["lambda.min"]] <- best_lambda
  fit[["select.features"]] <- selected_features
  fit[["select.coefficients"]] <- selected_coefficients
  
  return(fit)
}

# Function to calculate R2:
lassoR2 <- function(y, ypred){
  R2 <- 1 - sum((y - ypred)^2) / sum((y - mean(y))^2)
  return(R2)
}

# Function to obtain R2 of train and test set:
R2_lassofit <- function(protID, lasso_model, trainData, testData, formula){
  #Deal with different protID symbols:
  protID <- gsub("-", "_", protID)
  
  #Build Model Matrices for Prediction and Score Generation!
  contrast_list <- extract_contrasts(trainData)
  X_train <- model.matrix(formula, data = trainData, contrasts.arg = contrast_list)
  X_test <- model.matrix(formula, data = testData, contrasts.arg = contrast_list)
  
  #Predict on Train and Test Set: Worked after MAKING sure the LEVELS are consistent between Train and Test:
  ytrain_pred <- predict(lasso_model$lasso, newx = X_train, s = lasso_model$lambda.min)
  ytest_pred <- predict(lasso_model$lasso, newx = X_test, s = lasso_model$lambda.min)
  
  R2_train <- lassoR2(trainData[[protID]], ytrain_pred)
  R2_test <- lassoR2(testData[[protID]], ytest_pred)
  
  R2 <- data.frame(
    omic = protID,
    train_lasso = R2_train,
    test_lasso = R2_test
  )
  return(R2)
}

# Function to build formula to extract coefficients (especially for categorical variables with contrast)
get_featurenames <- function(protID, covars_list, E_ids){
  protID <- gsub("-", "_", protID)
  G_var <- paste0(protID,"_GS")
  
  covariates <- covars_list
  
  #For loop to build formula: E
  pred_E <- c()
  for(i in E_ids){
    #Setup Formula:
    pred_Ei <- i
    pred_E <- c(pred_E, pred_Ei)
  }
  
  #For loop to build formula: GxE components
  pred_GxE <- c()
  for(i in E_ids){
    #Setup Formula:
    pred_Ei <- paste0(G_var, ":", i)
    pred_GxE <- c(pred_GxE, pred_Ei)
  }
  
  #Final formula with G, E, and GxE components with covariates:
  formulaD <- as.formula(paste(c(protID), paste(covariates, collapse="+"), sep="~"))
  formulaE <- as.formula(paste(c(protID), paste(pred_E, collapse="+"), sep="~"))
  formulaGxE <- as.formula(paste(c(protID), paste(pred_GxE, collapse="+"), sep="~"))
  
  formula <- list()
  formula[["D"]] <- formulaD
  formula[["E"]] <- formulaE
  formula[["GxE"]] <- formulaGxE
  
  return(formula)
}

#'*Can obtain list of features extracted and coefficients for each protein (if necessary here)*
subset_featurecoefs <- function(coef_forms, trainData, selected_features, selected_coefficients){
  #Adding contrast_list and making the model.matrix consistent!!!
  #contrast_list <- extract_contrasts(trainData)
  design_matrix <- model.matrix(coef_forms, data = trainData) #, contrasts.arg = contrast_list)
  column_names <- colnames(design_matrix)
  
  features <- selected_features[selected_features %in% column_names]
  coefs <- selected_coefficients[selected_features %in% features]
  
  block <- list()
  block[["features"]] <- features
  block[["coefs"]] <- coefs
  
  return(block)
}

# Create score by multiplying coefficients and summing (Additive)
create_score <- function(df, features, coefs){
  #Build Model Matrices for Prediction and Score Generation!
  subset <-  df[,features,drop = FALSE]
  score <- subset %*% as.vector(coefs)
  return(score)
}

# Obtain Lasso Score for Train and Test Set
lasso_score_gen <- function(formula, trainData, testData,
                            features, coefs){
  contrast_list <- extract_contrasts(trainData)
  X_train <- model.matrix(formula, data = trainData, contrasts.arg = contrast_list)
  X_test <- model.matrix(formula, data = testData, contrasts.arg = contrast_list)
  
  train_score <- create_score(X_train, features, coefs)
  test_score <- create_score(X_test, features, coefs)
  
  Data <- list()
  Data[["score.train"]] <- train_score
  Data[["score.test"]] <- test_score
  
  return(Data)
}


# Function to update train and test Datasets with corresponding Scores
metrics_df <- function(trainData, testData, 
                       Dscore, Escore, GxEscore){
  trainScores <- data.frame(eid = trainData[["eid"]],
                            demo = Dscore$score.train,
                            E = Escore$score.train,
                            GxE = GxEscore$score.train)
  
  testScores <- data.frame(eid = testData[["eid"]],
                           demo = Dscore$score.test,
                           E = Escore$score.test,
                           GxE = GxEscore$score.test)
  
  #Scale train and test scores appropriately!
  #Scale numeric columns of continuous features::
  col_names <- c("demo","E","GxE")
  mean_train <- lapply(trainScores[col_names], function(x) mean(x,na.rm = T))
  sd_train <- lapply(trainScores[col_names], function(x) sd(x,na.rm = T))
  
  for(i in col_names){
    #Z-score train Scores:
    trainScores[[i]] <- (trainScores[[i]] - mean_train[[i]]) / sd_train[[i]]
    
    #Z-score test Scores using same parameters as train:
    testScores[[i]] <- (testScores[[i]] - mean_train[[i]]) / sd_train[[i]]
  }
  
  
  Data <- list()
  Data[["train"]] <- merge(trainData, trainScores, by = "eid")
  Data[["test"]] <- merge(testData, testScores, by = "eid")
  
  return(Data)
}

# Function to calculate partitioned R2 of joint model
calculate_partitioned_r2 <- function(model, anovaType) {
  # Get the ANOVA table
  if(anovaType > 1){
    anova_table <- Anova(model, type = anovaType)
  } else {
    anova_table <- anova(model)
  }
  
  # Extract the sums of squares
  ss_total <- sum(anova_table$`Sum Sq`)
  ss_residual <- anova_table["Residuals", "Sum Sq"]
  ss_factors <- anova_table[rownames(anova_table) != "Residuals", "Sum Sq"]
  
  # Calculate the total variance explained (R-squared)
  r2_total <- 1 - ss_residual / ss_total
  
  # Calculate the partitioned R-squared for each factor
  r2_factors <- ss_factors / ss_total
  
  # Combine the results into a data frame
  r2_results <- data.frame(
    Factor = rownames(anova_table)[rownames(anova_table) != "Residuals"],
    `Sum Sq` = ss_factors,
    `Partitioned R2` = r2_factors
  )
  
  # Add total R-squared to the results
  r2_results <- rbind(r2_results, data.frame(
    Factor = "Total",
    `Sum Sq` = ss_total,
    `Partitioned R2` = r2_total
  ))
  
  r2_results <- r2_results %>% column_to_rownames("Factor") %>% select("Partitioned.R2")
  
  return(r2_results)
}
joint_model_lm <- function(protID, data, covars){
  protID <- gsub("-", "_", protID)
  protGS <- paste0(protID,"_GS")
  
  res_var <- protID
  pred_var <- c(protGS, "E", "GxE", covars)
  
  #Errors in E and GxE are derived from not having the E or GxE score present:
  #Pseudocode:
  #If E or GxE is not an available score - remove this feature from the analysis:
  if(all(is.na(data[["GxE"]])) & all(is.na(data[["E"]]))){
    pred_var <- c(protGS, covars)
  } else if(all(is.na(data[["E"]]))){
    pred_var <- c(protGS,"GxE", covars)
  } else if(all(is.na(data[["GxE"]]))){
    pred_var <- c(protGS, "E", covars)
  }
  #print(pred_var)
  
  #Was used to determine that there might be errors in E and GxE
  # print(dim(data))
  # for(i in pred_var){
  #   print(summary(data[[i]]))
  # }
  formula <- as.formula(paste(c(res_var), paste(pred_var, collapse="+"), sep="~"))
  fit <- lm(formula, data = data)
  return(fit)
}

# Function to calculate paritioned R2 marginals (Very Biased)
#'*May not be needed since GxE term is collinear to G term and will cause issues!! - inflate R2*
marginal_model_lm <- function(protID, data){
  protID <- gsub("-", "_", protID)
  protGS <- paste0(protID,"_GS")
  
  res_var <- protID
  pred_var <- protGS
  formula1 <- as.formula(paste(c(res_var), paste(pred_var, collapse="+"), sep="~"))
  fit <- lm(formula1, data = data)
  print(summary(fit)$r.squared)
  
  pred_var <- "E"
  formula2 <- as.formula(paste(c(res_var), paste(pred_var, collapse="+"), sep="~"))
  fit <- lm(formula2, data = data)
  print(summary(fit)$r.squared)
  
  
  pred_var <- "GxE"
  formula3 <- as.formula(paste(c(res_var), paste(pred_var, collapse="+"), sep="~"))
  fit <- lm(formula3, data = data)
  print(summary(fit)$r.squared)
  
  return(fit)
  
}

# Function to organize R2 values to plot from joint model:
get_R2_terms <- function(protID, df){
  protID <- gsub("-", "_", protID)
  G_var <- paste0(protID,"_GS")
  
  #Now figure out what to do if G_var, E or GxE is missing
  #Current autohandling to NA
  vec <- c(protID, df[G_var,], df["E",], df["GxE",])
  names(vec) <- c("ID","G","E","GxE")
  vec
  
  R <- as.data.frame(t(vec))
  R[c("G","E","GxE")] <- sapply(R[c("G","E","GxE")],as.numeric)
  
  return(R)
}

# Function to organize R2 values to plot from joint model: Generalized
get_R2_terms_upd <- function(protID, df, covars){
  # protID <- gsub("-", "_", protID)
  # G_var <- paste0(protID,"_GS")
  # 
  # #Show the R2 for all features
  # R2 <- as.data.frame(t(df))
  # 
  # #Simple renaming
  # R2$ID <- protID
  # names(R2)[names(R2) == G_var] <- "G"
  # rownames(R2) <- NULL
  # 
  # return(R2)
  protID <- gsub("-", "_", protID)
  G_var <- paste0(protID, "_GS")
  
  features <- c("G","E","GxE",covars)
    
  # Create a vector with the protID and the feature values
  vec <- c(protID)
  
  # Add values for each feature, handle missing features by setting them to NA
  for (f in features) {
    feature_var <- if (f == "G") G_var else f
    vec[f] <- if (feature_var %in% rownames(df)) df[feature_var, ] else 0
  }
  
  # Convert to a data frame
  R <- as.data.frame(t(vec))
  colnames(R) <- c("ID", features)
  #R <- as.data.frame(lapply(R, function(x) as.numeric(as.character(x))))
  R[features] <- sapply(R[features],as.numeric)
  
  return(R)
}

# Function to update train and test Datasets with corresponding Scores for categories
metrics_df_cat <- function(trainData, testData, scores_list) {
  # Create train and test scores data frames
  trainScores <- data.frame(eid = trainData[["eid"]])
  testScores <- data.frame(eid = testData[["eid"]])
  
  # Add each score to the trainScores and testScores data frames
  for (score_name in names(scores_list)) {
    trainScores[[score_name]] <- scores_list[[score_name]]$score.train
    testScores[[score_name]] <- scores_list[[score_name]]$score.test
  }
  
  # Scale train and test scores appropriately
  col_names <- names(scores_list)
  mean_train <- lapply(trainScores[col_names], function(x) mean(x, na.rm = TRUE))
  sd_train <- lapply(trainScores[col_names], function(x) sd(x, na.rm = TRUE))
  
  for (i in col_names) {
    # Z-score train Scores
    trainScores[[i]] <- (trainScores[[i]] - mean_train[[i]]) / sd_train[[i]]
    
    # Z-score test Scores using same parameters as train
    testScores[[i]] <- (testScores[[i]] - mean_train[[i]]) / sd_train[[i]]
  }
  
  Data <- list()
  Data[["train"]] <- merge(trainData, trainScores, by = "eid")
  Data[["test"]] <- merge(testData, testScores, by = "eid")
  
  return(Data)
}

#Function just for categories:
joint_model_lm_cat <- function(protID, data, covars, scoreIDs = c("E", "GxE")) {
  # Replace "-" with "_" in protID
  protID <- gsub("-", "_", protID)
  protGS <- paste0(protID, "_GS")
  
  res_var <- protID
  pred_var <- c(protGS, scoreIDs, covars)
  
  # Check availability of scores in the data
  available_scores <- scoreIDs[!sapply(scoreIDs, function(feat) all(is.na(data[[feat]])))]
  
  # Update predictors based on available scores
  if (length(available_scores) == 0) {
    pred_var <- c(protGS, covars)
  } else {
    pred_var <- c(protGS, available_scores, covars)
  }
  
  # Create formula and fit the model
  formula <- as.formula(paste(res_var, paste(pred_var, collapse = "+"), sep = "~"))
  fit <- lm(formula, data = data)
  
  return(fit)
}

#Function to get R2 terms for categories:
get_R2_terms_cat <- function(protID, df, features) {
  protID <- gsub("-", "_", protID)
  G_var <- paste0(protID, "_GS")
  
  # Create a vector with the protID and the feature values
  vec <- c(protID)
  
  # Add values for each feature, handle missing features by setting them to NA
  for (feature in features) {
    feature_var <- if (feature == "G") G_var else feature
    vec[feature] <- if (feature_var %in% rownames(df)) df[feature_var, ] else 0
  }
  
  # Convert to a data frame
  R <- as.data.frame(t(vec))
  colnames(R) <- c("ID", features)
  #R <- as.data.frame(lapply(R, function(x) as.numeric(as.character(x))))
  R[features] <- sapply(R[features],as.numeric)
  
  return(R)
}


##'*Assess Variability: Nested CV approach*

#Out-of-Fold (Test-set) Predictions of PXS:
OOF_scores <- function(protID, scores){
  protID <- gsub("-", "_", protID)
  
  test_predictions <- scores$test
  
  # Make sure we are storing E and GxE PXS
  required_cols <- c("eid", "E", "GxE")
  
  # Setup to make sure combining PXS across folds works properly
  omic_PXS <- test_predictions %>% 
    select(any_of(c(required_cols))) %>% #Select available E and GxE scores.
    add_column(!!!setNames(rep(list(NA), length(setdiff(required_cols, names(.)))), 
                           setdiff(required_cols, names(.)))) %>% #NA missing variables
    select(all_of(required_cols)) #Reorder columns to ensure correct rbind
  
  # Add protID info for easy access later
  omic_PXS$ID <- protID
  
  return(omic_PXS)
}


# K fold CV:
library(caret) #for k-fold split
CV_split <-  function(protID, PXSdata, kfold){
  #References in PXSdata
  # Elist = PXSdata@Elist
  # Elist_names = PXSdata@Elist_names
  # Eid_cat = PXSdata@Eid_cat
  # ordinalIDs = PXSdata@ordinalIDs
  # covars_df = PXSdata@covars_df
  # covars_list = PXSdata@covars_list
  # UKBprot_df = PXSdata@UKBprot_df
  
  
  #Load proteomics data:
  omicDS <- GS_struct(protID, PXSdata@UKBprot_df)
  protID <- gsub("-", "_", protID)
  protGS <- paste0(protID,"_GS")
  
  E_df <- PXSdata@Elist %>% reduce(full_join)

  df <- merge(omicDS$combo, E_df, by = "eid")
  df <- merge(df, PXSdata@covars_df, by = "eid")
  
  # Remove columns with too much missingness:
  removecols <- remove_cols_name(df)
  df <- remove_missing_cols(df)
  print(paste0("Too many missing values",removecols))
  
  # CHECK 1: Split occurs with the maximal amount of proteins available:
  # A.) CHECKED YES!!
  
  #Environmental ids & Scaling || Dealing with Ordinals / Categoricals Appropriately:
  E_ids <- colnames(E_df)[-1]
  E_ids <- E_ids[!(E_ids %in% removecols)]
  
  # Categorical Handler BEFORE Splitting:
  df <- categorical_handler(df, PXSdata@ordinalIDs, ordinal_contrast = "treatment")
  
  # Create k-folds
  folds <- createFolds(seq_len(nrow(df)), k = kfold, list = TRUE)
  
  CVdata <- list()
  for(i in seq_along(folds)){
    valIndex <- folds[[i]]
    
    trainData <- df[-valIndex, ]
    testData <- df[valIndex, ]
    
    CVdata[["train"]][[i]] <- trainData
    CVdata[["test"]][[i]] <- testData
  }
  
  # Store important values:
  CVdata$Eids <- E_ids
  CVdata$ordinalVar <- PXSdata@ordinalIDs #Assumes ordinalIDs do not contain a lot of missing info.
  
  return(CVdata)
}


# Function to obtain stats using lasso derived PXS score:
CV_protPXSrun <- function(protID, folds, PXSdata, idx){
  #Specific Features that apply to all folds:
  Elist_names = PXSdata@Elist_names
  Eid_cat = PXSdata@Eid_cat
  covars_list = PXSdata@covars_list
  
  #Formula and Dataset:
  protData <- list()
  protData[["train"]] <- folds[["train"]][[idx]]
  protData[["test"]] <- folds[["test"]][[idx]]
  protData[["Eids"]] <- folds[["Eids"]]
  protData[["ordinalVar"]] <- folds[["ordinalVar"]]
  
  formula <- build_formula(protID, covars_list, protData$Eids)
  #print(formula)
  
  pData <- preprocess_data(protID, covariates = covars_list,
                           trainData = protData$train,
                           testData = protData$test,
                           ordinal_contrast = "treatment",
                           ordinal_names = protData$ordinalVar) 
  
  
  # Fit Lasso
  lasso <- lasso_fit(protID, trainData = pData$train, formula)
  #print(lasso)
  
  # Build G and GxE scores:
  featnames <- get_featurenames(protID, covars_list, E_ids = protData$Eids)
  
  
  
  #D, E & GxE Features:
  Dblock <- subset_featurecoefs(coef_forms = featnames$D, trainData = pData$train, 
                                selected_features = lasso$select.features,
                                selected_coefficients = lasso$select.coefficients)
  Eblock <- subset_featurecoefs(coef_forms = featnames$E, trainData = pData$train,
                                selected_features = lasso$select.features, 
                                selected_coefficients = lasso$select.coefficients)
  GxEblock <- subset_featurecoefs(coef_forms = featnames$GxE, trainData = pData$train,
                                  selected_features = lasso$select.features, 
                                  selected_coefficients = lasso$select.coefficients)
  
  
  Dscore <- lasso_score_gen(formula, trainData = pData$train,
                            testData = pData$test,
                            features = Dblock$features,
                            coefs = Dblock$coefs)
  Escore <- lasso_score_gen(formula, trainData = pData$train,
                            testData = pData$test,
                            features = Eblock$features,
                            coefs = Eblock$coefs)
  GxEscore <- lasso_score_gen(formula, trainData = pData$train,
                              testData = pData$test,
                              features = GxEblock$features,
                              coefs = GxEblock$coefs)
  
  #Updated Dataframe w/ Scores:
  scores <- metrics_df(trainData = pData$train,
                       testData = pData$test,
                       Dscore, Escore, GxEscore)
  
  
  #'*Out of Fold Predictions*
  test_PXS <- OOF_scores(protID, scores)
  
  
  
  # :: Stats for Metrics using the Scores ::
  
  joint_train <- joint_model_lm(protID, scores$train, covars = covars_list)
  joint_test <- joint_model_lm(protID, scores$test, covars = covars_list)
  
  train_partR2_v1 <- calculate_partitioned_r2(joint_train, anovaType = 1)
  test_partR2_v1 <- calculate_partitioned_r2(joint_test, anovaType = 1)
  
  
  #'*Split by Category: Not including GxE*
  Ecategories <- list()
  for(i in Elist_names){
    pred_E <- Eid_cat[Eid_cat$Category == i, ]$Eid
    
    #If not appropriate still use:
    #pred_E <- gsub(" ","_", pred_E)
    #pred_E <- make.names(pred_E)
    pred_E <- intersect(pred_E, protData$Eids) #make sure the calls contain only filtered Eids...
    
    #R2 for every E category:
    protID <- gsub("-", "_", protID)
    formulaE_cat <- as.formula(paste(c(protID), paste(pred_E, collapse="+"), sep="~"))
    
    catEblock <- subset_featurecoefs(coef_forms = formulaE_cat, trainData = pData$train, 
                                     selected_features = lasso$select.features,
                                     selected_coefficients = lasso$select.coefficients)
    catEscore <- lasso_score_gen(formula, trainData = pData$train,
                                 testData = pData$test,
                                 features = catEblock$features,
                                 coefs = catEblock$coefs)
    Ecategories[[i]] <- catEscore
  }
  #PXS_bycat <<- Ecategories
  #PXS <<- Escore
  
  scores_cat <- metrics_df_cat(trainData = pData$train,
                               testData = pData$test,
                               Ecategories)
  
  fit.train <- joint_model_lm_cat(protID, data = scores_cat$train, covars = covars_list, scoreIDs = Elist_names)
  fit.test <- joint_model_lm_cat(protID, data = scores_cat$test, covars = covars_list, scoreIDs = Elist_names)
  
  train_R2cat <- calculate_partitioned_r2(fit.train, anovaType = 1)
  test_R2cat <- calculate_partitioned_r2(fit.test, anovaType = 1)
  
  
  
  #'*STATS:: R2 on train and test set (lasso,PXS,etc.)*
  lasso_fit_metrics <- R2_lassofit(protID,
                                   lasso_model = lasso,
                                   trainData = pData$train,
                                   testData = pData$test,
                                   formula)
  
  metrics <- list()
  metrics[["R2.lasso.fit"]] <- lasso_fit_metrics
  
  #Report Train and Test R2 (G, E, GxE, etc.)
  metrics[["partR2.train"]] <- get_R2_terms_upd(protID, train_partR2_v1, covars_list)
  metrics[["partR2.test"]] <- get_R2_terms_upd(protID, test_partR2_v1, covars_list)
  
  
  #Report Train and Test R2 for E Categories
  metrics[["partR2.train.cat"]] <- get_R2_terms_cat(protID, train_R2cat, features = Elist_names)
  metrics[["partR2.test.cat"]] <- get_R2_terms_cat(protID, test_R2cat, features = Elist_names)
  
  #Report test set PXS:
  metrics[["test.PXS"]] <- test_PXS
  
  return(metrics)
}


# Run everything together:
prot_PXS_CVrun <- function(protID, PXSdata, k){
  #Assumes Elist, Elist_names, and Eid_cat are properly encoded from dataloader:
  # Elist = PXSdata@Elist
  # Elist_names = PXSdata@Elist_names
  # Eid_cat = PXSdata@Eid_cat
  # ordinalIDs = PXSdata@ordinalIDs
  # covars_df = PXSdata@covars_df
  # covars_list = PXSdata@covars_list
  # UKBprot_df = PXSdata@UKBprot_df
  
  #'*Make folds observable in global env.*
  folds <- CV_split(protID, PXSdata, kfold = k) 
  
  Lasso <- data.frame()
  R2train <- data.frame()
  R2test <- data.frame()
  R2trainCat <- data.frame()
  R2testCat <- data.frame()
  OOF_PXS <- data.frame()
  for(i in 1:k){
    model <- CV_protPXSrun(protID, folds, PXSdata, idx = i)
    
    #Updates
    Lasso <- rbind(Lasso, model$R2.lasso.fit)
    R2train <- rbind(R2train, model$partR2.train)
    R2test <- rbind(R2test, model$partR2.test)
    R2trainCat <- rbind(R2trainCat, model$partR2.train.cat)
    R2testCat <- rbind(R2testCat, model$partR2.test.cat)
    OOF_PXS <- rbind(OOF_PXS, model$test.PXS)
    
    
    print(paste0("Fold: ", i))
  }
  
  metrics <- list()
  metrics[["Lasso"]] <- Lasso
  metrics[["R2train"]] <- R2train
  metrics[["R2test"]] <- R2test
  metrics[["R2trainCat"]] <- R2trainCat
  metrics[["R2testCat"]] <- R2testCat
  metrics[["OOF_PXS"]] <- OOF_PXS 
  
  return(metrics)
}


#'*Running Nested CV across Proteome*
## PGS/PXS run with partitioned R2 for each category.
omiclist <- scan(file = "/n/groups/patel/shakson_ukb/UK_Biobank/BScripts/ProtPGS_PXS/OMICPREDproteins.txt",
                 as.character())

args = commandArgs(trailingOnly = T)
idx <- as.integer(args[1])
split_num <- as.integer(args[2])
covarType <- as.character(args[3])

# Create a grouping factor to split the vector 
groups <- cut(seq_along(omiclist), breaks = split_num, labels = FALSE)
split_vectors <- split(omiclist, groups)

# Function to create PXS scores and metrics:
runProt_PGS_PXS_multi <- function(protlist, PXSdata, folder_id, idx){
  #Create folder if nonexistent
  dir.create(
    paste0("/n/groups/patel/shakson_ukb/UK_Biobank/Data/Parallel/Prot_PGSPXS_final/",folder_id),
    showWarnings = FALSE
    )
  
  
  Lasso <- data.frame()
  R2train <- data.frame()
  R2test <- data.frame()
  R2trainCat <- data.frame()
  R2testCat <- data.frame()
  
  for(i in protlist){
    model <- prot_PXS_CVrun(protID = i, 
                            PXSdata,
                            k = 10)
    Lasso <- rbind(Lasso, model$Lasso)
    R2train <- rbind(R2train, model$R2train)
    R2test <- rbind(R2test, model$R2test)
    R2trainCat <- rbind(R2trainCat, model$R2trainCat)
    R2testCat <- rbind(R2testCat, model$R2testCat)

    #Output: PXS for each protein.
    write.csv(model$OOF_PXS,
                file = paste0("/n/groups/patel/shakson_ukb/UK_Biobank/Data/Parallel/Prot_PGSPXS_final/",folder_id,"/","PXS_",i,".csv"),
              row.names = F)
  }
  
  
  
  
  #Write metrics for grouped proteins.
  write.table(Lasso,
              file = paste0("/n/groups/patel/shakson_ukb/UK_Biobank/Data/Parallel/Prot_PGSPXS_final/",folder_id,"/","lassofit_",idx,".txt"),
              row.names = F)
  write.table(R2train,
              file = paste0("/n/groups/patel/shakson_ukb/UK_Biobank/Data/Parallel/Prot_PGSPXS_final/",folder_id,"/","R2train_",idx,".txt"),
              row.names = F)
  write.table(R2test,
              file = paste0("/n/groups/patel/shakson_ukb/UK_Biobank/Data/Parallel/Prot_PGSPXS_final/",folder_id,"/","R2test_",idx,".txt"),
              row.names = F)
  write.table(R2trainCat,
              file = paste0("/n/groups/patel/shakson_ukb/UK_Biobank/Data/Parallel/Prot_PGSPXS_final/",folder_id,"/","R2trainCat_",idx,".txt"),
              row.names = F)
  write.table(R2testCat,
              file = paste0("/n/groups/patel/shakson_ukb/UK_Biobank/Data/Parallel/Prot_PGSPXS_final/",folder_id,"/","R2testCat_",idx,".txt"),
              row.names = F)
} 

# Function to update PXSloader to subset of covariates selected.
PXScovarSpec <- function(PXSdata, covars_subset){
  PXSdata@covars_list <- covars_subset
  PXSdata@covars_df <- PXSdata@covars_df %>% select(all_of(c("eid",covars_subset)))
  
  return(PXSdata)
}

CovarSpec <- list()
CovarSpec[["Type1"]] <- c("age_when_attended_assessment_centre_f21003_0_0",
                           "sex_f31_0_0")
CovarSpec[["Type2"]] <- c("age_when_attended_assessment_centre_f21003_0_0",
                           "sex_f31_0_0",
                           "body_mass_index_bmi_f23104_0_0",
                           "fasting_time_f74_0_0")
CovarSpec[["Type3"]] <- c("age_when_attended_assessment_centre_f21003_0_0",
                           "sex_f31_0_0",
                           "age2","age_sex", "age2_sex",
                           "body_mass_index_bmi_f23104_0_0",
                           "fasting_time_f74_0_0",
                           "uk_biobank_assessment_centre_f54_0_0",
                           paste0("genetic_principal_components_f22009_0_",1:20)
                           )
CovarSpec[["Type4"]] <- c("age_when_attended_assessment_centre_f21003_0_0",
                           "sex_f31_0_0",
                           "age2","age_sex", "age2_sex",
                           "body_mass_index_bmi_f23104_0_0",
                           "fasting_time_f74_0_0",
                           "uk_biobank_assessment_centre_f54_0_0",
                           paste0("genetic_principal_components_f22009_0_",1:20),
                           "combined_Blood_pressure_medication",
                           "combined_Hormone_replacement_therapy",
                           "combined_Oral_contraceptive_pill_or_minipill",
                           "combined_Insulin",
                           "combined_Cholesterol_lowering_medication",
                           "combined_Do_not_know",
                           "combined_None_of_the_above",
                           "combined_Prefer_not_to_answer"
                           )
CovarSpec[["Type5"]] <- PXSloader@covars_list



runProt_PGS_PXS_multi(split_vectors[[idx]], 
                      PXSdata = PXScovarSpec(PXSloader, CovarSpec[[covarType]]), 
                      folder_id = covarType,
                      idx)