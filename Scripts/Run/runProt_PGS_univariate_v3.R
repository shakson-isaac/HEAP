#Libraries:
library(data.table)
library(tidyverse)
library(ggplot2)
library(ggpmisc)
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


#'*Functions: Preprocessing*
# Function to find continuous features
continuous_finder <- function(df){
  max_cols <- apply(df, 2, max, na.rm = TRUE)
  min_cols <- apply(df, 2, min, na.rm = TRUE) 
  unique_vals <- sapply(df, function(x) length(unique(x[!is.na(x)])))  # Number of unique values
  #^unique_vals handles cases of continuous features that range from 0 to 1 but are not one-hot encoded.
  # Identify continuous columns: not one-hot encoded (0 and 1) or values greater than > 5 as max
  continuous <- names(max_cols[max_cols > 5 | unique_vals > 2])
  
  return(continuous)
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

# To-DO: OUTPUT WARNING INSTEAD OF ERROR FOR ANY FILES MISSING
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

# Function to Obtain PGS of Protein
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


######FUNCTION 1: LOADING and Simple Preprocess:
#Idea: For train only put: trainIndex amount as 100%

# Function to create Train/Test Split:
create_data <- function(protID, PXSdata, split_ratio){
  #Load proteomics data: + more
  omicDS <- GS_struct(protID, PXSdata@UKBprot_df)
  protID <- gsub("-", "_", protID)
  protGS <- paste0(protID,"_GS")
  
  E_df <- PXSdata@Elist %>% reduce(full_join)
  
  df <- merge(omicDS$combo, E_df, by = "eid")
  df <- merge(df, PXSdata@covars_df, by = "eid")
  
  #Environmental ids & Scaling || Dealing with Ordinals / Categoricals Appropriately:
  E_ids <- colnames(E_df)[-1]
  
  # Categorical Handler BEFORE Splitting:
  df <- categorical_handler(df, PXSdata@ordinalIDs, ordinal_contrast = "treatment")
  
  # Train/Test split:
  Edata <- list()
  #May need to try 70/30 and 60/40 split based on power analysis (apriori)
  trainIndex <- sample(seq_len(nrow(df)), size = split_ratio * nrow(df))
  Edata[["train"]] <- df[trainIndex, ]
  Edata[["test"]] <- df[-trainIndex, ]
  
  # Store important values:
  Edata$Eids <- E_ids
  Edata$ordinalVar <- PXSdata@ordinalIDs #Assumes ordinalIDs do not contain a lot of missing info.
  
  return(Edata)
}
#TEST CASE to have for above function:
#nrow(Edata$train) + nrow(Edata$test) == length(unique(c(Edata$train$eid, Edata$test$eid)))


######FUNCTION 2: MORE PREPROCESSING!!

# Function to preprocess data under train/test split:
preprocess_data <- function(protID, covariates,
                            trainData, testData, ordinal_contrast = "treatment",
                            ordinal_names){
  # Setup protID and Genetic Score ID:
  protID <- gsub("-", "_", protID)
  G_var <- paste0(protID,"_GS")
  
  
  #na.omit: make complete data (NO need since this is univariate!!)
  #trainData <- na.omit(trainData)
  #testData <- na.omit(testData)
  
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



###### FUNCTION 3: lm Univariate Analysis

# Function to run E and GxE univariate analysis:
GxE_assoc_fin <- function(df, protID, covariates, Eids){
  full_stat <- list()
  stat1 <- list() # G terms
  stat2 <- list() # GxE terms
  stat3 <- list() # Partitioned R2
  stat4 <- list() # F tests
  stat5 <- list() # F tests
  
  # Setup protID and Genetic Score ID:
  protID <- gsub("-", "_", protID)
  res_var <- protID
  G_var <- paste0(protID,"_GS")
  
  
  for(i in Eids){
    #Setup Formula:
    pred_var <- c(i, G_var, paste(G_var, "*", i), covariates)
    formula <- as.formula(paste(c(res_var), paste(pred_var, collapse="+"), sep="~"))
    
    #Subset df:
    df2 <- df %>% select(all_of(c(i,res_var,G_var,covariates)))
    
    #df <- na.omit(df)
    
    #Linear Regression:
    fit <- lm(formula, data = df2)
    res <- summary(fit)
    
    #Coefficients for lm:
    CI_Edf <- res$coefficients
    CI_Edf <- as.data.frame(CI_Edf)
    
    #R2 & N for lm:
    CI_Edf$R2 <- summary(fit)$r.squared
    CI_Edf$adj.R2 <- summary(fit)$adj.r.squared
    CI_Edf$samplesize <- length(fit$fitted.values)
    
    stat_tbl <- CI_Edf %>%
      rownames_to_column(var = "ID") %>%
      pivot_longer(cols = colnames(CI_Edf), 
                   names_to = "stats", 
                   values_to = "value")
    
    stat1[[i]] <- stat_tbl %>% filter(grepl(i,ID) & !grepl(G_var,ID)) #Get E terms
    stat2[[i]] <- stat_tbl %>% filter(grepl(paste0(":",G_var), ID)) #Get GxE Terms
    
    #Partitioned R2:
    partR2 <- calculate_partitioned_r2(fit, anovaType = 1)
    
    partR2_df <- data.frame(
      ID = c(i, paste0(i,":",G_var)),
      R2 = c(partR2[i,], partR2[paste0(i,":",G_var),])
    )
    
    stat3[[i]] <- partR2_df
    
    # F test: Type 1 Anova
    anova_res <- anova(fit)
    
    E_ftest <- i
    GxE_ftest <- paste0(i,":",G_var)

    aov_type1 <- anova_res[rownames(anova_res) %in% c(E_ftest, GxE_ftest), ]
    aov_type1$ID <- rownames(aov_type1)
    stat4[[i]] <- aov_type1
    
    # F test: Type 2 Anova
    anova_res2 <- Anova(fit, type = 2)
    
    aov_type2 <- anova_res2[rownames(anova_res2) %in% c(E_ftest, GxE_ftest), ]
    aov_type2$ID <- rownames(aov_type2)
    
    stat5[[i]] <- aov_type2
    
    
    #Old Version:
    #stat1[[i]] <- stat_tbl %>% filter(id == i)
    #stat2[[i]] <- stat_tbl %>% filter(id == paste0(i,":",G_var))
  }
  #print(stat3)
  
  #Aggregate the stats to make a dataframe
  stat1 <- do.call(rbind, stat1)
  stat1 <- stat1 %>% pivot_wider(names_from = stats, values_from = value)
  
  stat2 <- do.call(rbind, stat2)
  stat2 <- stat2 %>% pivot_wider(names_from = stats, values_from = value)
  
  stat3 <- do.call(rbind,stat3)
  rownames(stat3) <- NULL
  
  stat4 <- do.call(rbind,stat4)
  
  stat5 <- do.call(rbind,stat5)
  
  full_stat <- list(stat1, stat2, stat3, stat4, stat5)
  return(full_stat)
}



##### FUNCTION 4: Aggregate results to store:

# Create a function to find all matches
find_all_matches <- function(partial_id, unique_ids) {
  matches <- unique_ids[str_detect(unique_ids, partial_id)]
  
  #str_detect(unique, pattern)
  if (length(matches) > 0) {
    return(matches)  # Return all matches
  } else {
    return(NA  # Return NA if no match is found
    )
  }
}

# Create ID matching for stat columns:
add_cat_info <- function(Eid_cat, stat_df){
  df <- Eid_cat
  
  # Get list of matches:
  df$ID <- lapply(df$Eid, find_all_matches, unique_ids = stat_df$ID)
  
  df_long <- df %>%
    unnest(ID)  # Expands the list into multiple rows
  
  # Merge stat_df to get Category IDs
  stat_df <- merge(df_long, stat_df, by = "ID")
  
  return(stat_df)
}



#### Function 5: Combine all univariate R2s together.

# Function to RUN Univariate Associations in batches:
runProt_univar <- function(protlist, PXSdata){
  stat_all <- list()
  stat_batch_train <- NULL
  stat_batch_test <- NULL
  
  for(i in protlist){
    protData <- create_data(protID = i, PXSdata, split_ratio = 0.8)
    
    protPP <- preprocess_data(protID = i, covariates = PXSdata@covars_list,
                              trainData = protData$train,
                              testData = protData$test,
                              ordinal_contrast = "treatment",
                              ordinal_names = protData$ordinalVar)
    
    stat_train <- GxE_assoc_fin(df = protPP$train, protID = i, 
                                covariates = PXSdata@covars_list, 
                                Eids = protData$Eids)
    
    stat_test <- GxE_assoc_fin(df = protPP$test, protID = i, 
                               covariates = PXSdata@covars_list, 
                               Eids = protData$Eids)
    
    # Attach Category Info into Stats
    stat_train <- lapply(stat_train, 
                         function(x) add_cat_info(PXSdata@Eid_cat, x))
    
    stat_test <- lapply(stat_test, 
                        function(x) add_cat_info(PXSdata@Eid_cat, x))
    
    # attach omic ID info:
    stat_train <- map(stat_train, ~ mutate(.x, omicID = i))
    stat_test <- map(stat_test, ~ mutate(.x, omicID = i))
    
    
    #Combine stats of each protein together:
    stat_batch_train <- map2(stat_batch_train %||% vector("list", length(stat_train)), 
                            stat_train, ~ {
                              .x %||% .y %>% full_join(.y, by = NULL)
                            })
    
    stat_batch_test <- map2(stat_batch_test %||% vector("list", length(stat_test)), 
                            stat_test, ~ {
                              .x %||% .y %>% full_join(.y, by = NULL)
                              })
    
  }
  
  stat_all[["train"]] <- stat_batch_train
  stat_all[["test"]] <- stat_batch_test
  
  return(stat_all)
}


#### Function 6: Run Different Covariate Specifications and Simple SWAPs (of E list vs Covars list)

# Function to organize Category Names:
# Function to label Eids, Categories, Etc.*
orgElistnames <- function(Elist_names){
  #Make the naming structure consistent:
  Elist_names <- gsub(" ","_", Elist_names)
  Elist_names <- gsub(paste(c("[(]", "[)]"), collapse = "|"),"",Elist_names)
  
  return(Elist_names)
}

#Function to organize Eids Names
orgEids <- function(Elist){
  for(i in names(Elist)){
    Eids <- gsub(" ","_", colnames(Elist[[i]]))
    Eids <- make.names(Eids)
    colnames(Elist[[i]]) <- Eids
  }
  return(Elist)
}

#Function to link Category Names to Eids
orgEtoCat <- function(Elist, Elist_names){
  Ecategory <- lapply(Elist, function(x) setdiff(colnames(x), "eid"))
  #total_entries <- sum(sapply(Ecategory, length))
  
  names(Ecategory) <- Elist_names
  
  Eid_cat <- stack(Ecategory)
  colnames(Eid_cat) <- c("Eid", "Category")
  Eid_cat$Category <- as.character(Eid_cat$Category)
  
  return(Eid_cat)
}


# Function to move: Covariate Feats --> E Feats
covariate_to_E <- function(PXSnew, FeatIds, FeatCat){
  #FeatIds - list of ids corresponding to new category
  #FeatCat - category name
  
  #Obtain df to full join later.
  df <- PXSnew@covars_df %>% select(all_of(c("eid",FeatIds)))
  
  # Attach new dataframe to list of E variables:
  PXSnew@Elist[[FeatCat]] <- df
  
  # Update list of category names
  PXSnew@Elist_names <- orgElistnames(names(PXSnew@Elist))
  
  # Update Eids to keep proper formatting:
  PXSnew@Elist <- orgEids(PXSnew@Elist)
  
  # Update Eid_to_Category
  PXSnew@Eid_cat <- orgEtoCat(PXSnew@Elist, PXSnew@Elist_names)
  
  # Remove FeatIds from Covariate Df and Covariate list names
  PXSnew@covars_df <- PXSnew@covars_df %>% 
    select(-all_of(c(FeatIds)))
  
  PXSnew@covars_list <- PXSnew@covars_list[!(PXSnew@covars_list %in% FeatIds)]
  
  #To-ADD Later: Update of Ordinal IDs - Will not be useful for this updating scheme
  return(PXSnew)
  
}


# Function to move: E Feats --> Covariate Feats
E_to_Covariate <- function(PXSnew, Eids, Ecat){
  
  # Obtain df of Efeatures to transfer to covars
  Edf <- PXSnew@Elist[[Ecat]] %>%
    select(all_of(c("eid",Eids)))
  
  # Add these Efeatures to covars dataframe:
  PXSnew@covars_df <- merge(PXSnew@covars_df, Edf, by = "eid")
  
  # Add Efeature names to covariates: Append
  PXSnew@covars_list <- c(PXSnew@covars_list, Eids)
  
  
  # Remove Category from Elist to Assess:
  PXSnew@Elist[[Ecat]] <- NULL
  
  # Remove this Category from list of potential names:
  PXSnew@Elist_names <- PXSnew@Elist_names[!(PXSnew@Elist_names %in% Ecat)]
  
  #Update Eid_to_Category dataframe after removal:
  PXSnew@Eid_cat <- orgEtoCat(PXSnew@Elist, PXSnew@Elist_names)
  
  #To-ADD Later: Update of Ordinal IDs - Will not be useful for this updating scheme
  return(PXSnew)
}


# Function to update PXSloader to subset of covariates selected.
PXScovarSpec <- function(PXSdata, covars_subset){
  PXSdata@covars_list <- covars_subset
  PXSdata@covars_df <- PXSdata@covars_df %>% select(all_of(c("eid",covars_subset)))
  
  return(PXSdata)
}


#'*List of CovarSpecs*
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



#'*RUN Univariate Associations w/ Corresponding Specifications:*

# Obtain OmicList:
omiclist <- scan(file = "/n/groups/patel/shakson_ukb/UK_Biobank/BScripts/ProtPGS_PXS/OMICPREDproteins.txt",
                 as.character())

# Obtain arguments from bash script:
args = commandArgs(trailingOnly = T)
idx <- as.integer(args[1])
split_num <- as.integer(args[2])
covarType <- as.character(args[3])

#Example arguments:
#idx = 1
#split_num = 400
#covarType = "Type6"

# Create a grouping factor to split the vector 
groups <- cut(seq_along(omiclist), breaks = split_num, labels = FALSE)
split_vectors <- split(omiclist, groups)

# Function to run and store Univarate Associations:
runProt_univar_assoc <- function(protlist, PXSdata, folder_id, idx){
  #Create folder if nonexistent
  dir.create(
    paste0("/n/groups/patel/shakson_ukb/UK_Biobank/Data/Parallel/Prot_PGSPXS_univar/",folder_id),
    showWarnings = FALSE
  )
  
  batch_run <- runProt_univar(protlist, PXSdata)
  
  saveRDS(batch_run, file = paste0("/n/groups/patel/shakson_ukb/UK_Biobank/Data/Parallel/Prot_PGSPXS_univar/",
                                     folder_id,"/univar_assoc_",idx,".rds"))
}


### CovarSpec IDs:
# Type 1 -- Type 5: Covariate Specifications (General E list)
# Type 6: Deprivation Indices (England)
# Type 7: Variant of Type5 where medications is E not a confounder.

# Set condition on the PXSloader: Current is ALL covariates:
# Later can introduce: PXScond = PXScovarSpec(PXSloader, CovarSpec[[covarType]])
# PXScond <- PXSloader


if(covarType %in% c(paste0("Type",1:5))){
  runProt_univar_assoc(split_vectors[[idx]], 
                       PXSdata = PXScovarSpec(PXSloader, CovarSpec[[covarType]]), 
                       folder_id = covarType,
                       idx)
  } else if(covarType == "Type6") {
  
    #E to Covars: Deprivation Indices
    PXS_SES <- E_to_Covariate(PXSnew = PXSloader, 
                              c("income_score_england_f26411_0_0",
                                "index_of_multiple_deprivation_england_f26410_0_0"), 
                              "Deprivation_Indices")
    
    runProt_univar_assoc(split_vectors[[idx]], 
                         PXSdata = PXS_SES, 
                         folder_id = covarType,
                         idx)
  } else if(covarType == "Type7"){
    
    #Covars Medications to E:
    medication_ids <- PXSloader@covars_list[grepl("combined",PXSloader@covars_list)]
    medication_cat <- "Medications"
    
    PXSmedn <- covariate_to_E(PXSloader, medication_ids, medication_cat)
    
    runProt_univar_assoc(split_vectors[[idx]], 
                         PXSdata = PXSmedn, 
                         folder_id = covarType,
                         idx)
}
