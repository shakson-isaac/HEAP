#Build Function to run association:
GxE_assoc_ver2 <- function(df, E_list, res_var, G_var, covariates, 
                           feature_eng = "ordinal", ordinal_vars = c(),
                           ordinal_contrast = c()){
  full_stat <- list()
  stat1 <- list()
  stat2 <- list()
  
  #Quick scaling procedure: mean centered
  #allow for intercept to not be significant and around 0.
  #df <- na.omit(df) # Remove rows with NA values
  df <- df %>% column_to_rownames("eid")
  
  
  #Current Version: Scale Everything
  if(feature_eng == "scale"){
    numeric_cols <- sapply(df, is.numeric)
    numeric_col_names <- names(numeric_cols[numeric_cols])
    
    # Scale only numeric columns
    df[numeric_col_names] <- scale(df[numeric_col_names])
  } else if(feature_eng == "ordinal"){
    #Scale only the omics, genetic score, and covariates:
    rel_columns <- c(res_var, G_var, covariates)
    numeric_cols <- sapply(df[,rel_columns], is.numeric)
    numeric_col_names <- names(numeric_cols[numeric_cols])
    
    #Scale numeric columns of continuous features::
    df[numeric_col_names] <- scale(df[numeric_col_names])
    
    
    #Declare Ordinal Variables:
    for(i in ordinal_vars){
      df[[i]] <- factor(df[[i]], ordered = T)
      
      # Set treatment contrasts: instead of polynomial contrasts:
      if(ordinal_contrast == "treatment"){
        contrasts(df[[i]]) <- contr.treatment(length(levels(df[[i]])))
      } else if (ordinal_contrast == "sum"){
        contrasts(df[[i]]) <- contr.sum(length(levels(df[[i]])))
      } else {
        #df[[i]] <- df[[i]]
      }
    }
    
  } else {
    #Scale only the omics, genetic score, and covariates:
    rel_columns <- c(res_var, G_var, covariates)
    numeric_cols <- sapply(df[,rel_columns], is.numeric)
    numeric_col_names <- names(numeric_cols[numeric_cols])
    
    #Scale numeric columns of continuous features::
    df[numeric_col_names] <- scale(df[numeric_col_names])
  }
  
  for(i in E_list){
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
      rownames_to_column(var = "id") %>%
      pivot_longer(cols = colnames(CI_Edf), 
                   names_to = "stats", 
                   values_to = "value")
    
    stat1[[i]] <- stat_tbl %>% filter(grepl(i,id) & !grepl(G_var,id)) #Get E terms
    stat2[[i]] <- stat_tbl %>% filter(grepl(paste0(":",G_var), id)) #Get GxE Terms
    
    #Old Version:
    #stat1[[i]] <- stat_tbl %>% filter(id == i)
    #stat2[[i]] <- stat_tbl %>% filter(id == paste0(i,":",G_var))
  }
  
  #Print last res to double check:
  print(res)
  
  stat1 <- do.call(rbind, stat1)
  stat1 <- stat1 %>% pivot_wider(names_from = stats, values_from = value)
  
  stat2 <- do.call(rbind, stat2)
  stat2 <- stat2 %>% pivot_wider(names_from = stats, values_from = value)
  
  full_stat <- list(stat1, stat2)
  return(full_stat)
}

GxE_assoc_model <- function(df, E_list, res_var, G_var, covariates, 
                            feature_eng = "ordinal", ordinal_vars = c(),
                            ordinal_contrast = c()){
  models <- list()
  
  #Quick scaling procedure: mean centered
  #allow for intercept to not be significant and around 0.
  #df <- na.omit(df) # Remove rows with NA values
  df <- df %>% column_to_rownames("eid")
  
  #Current Version: Scale Everything
  if(feature_eng == "scale"){
    numeric_cols <- sapply(df, is.numeric)
    numeric_col_names <- names(numeric_cols[numeric_cols])
    
    # Scale only numeric columns
    df[numeric_col_names] <- scale(df[numeric_col_names])
  } else if(feature_eng == "ordinal"){
    #Scale only the omics, genetic score, and covariates:
    rel_columns <- c(res_var, G_var, covariates)
    numeric_cols <- sapply(df[,rel_columns], is.numeric)
    numeric_col_names <- names(numeric_cols[numeric_cols])
    
    #Scale numeric columns of continuous features::
    df[numeric_col_names] <- scale(df[numeric_col_names])
    
    #Declare Ordinal Variables:
    for(i in ordinal_vars){
      df[[i]] <- factor(df[[i]], ordered = T)
      
      # Set treatment contrasts: instead of polynomial contrasts:
      if(ordinal_contrast == "treatment"){
        contrasts(df[[i]]) <- contr.treatment(length(levels(df[[i]])))
      } else if (ordinal_contrast == "sum"){
        contrasts(df[[i]]) <- contr.sum(length(levels(df[[i]])))
      } else {
        #df[[i]] <- df[[i]]
      }
    }
  } else {
    #Scale only the omics, genetic score, and covariates:
    rel_columns <- c(res_var, G_var, covariates)
    numeric_cols <- sapply(df[,rel_columns], is.numeric)
    numeric_col_names <- names(numeric_cols[numeric_cols])
    
    #Scale numeric columns of continuous features::
    df[numeric_col_names] <- scale(df[numeric_col_names])
  }
  
  for(i in E_list){
    #Setup Formula:
    pred_var <- c(i, G_var, paste(G_var, "*", i), covariates)
    formula <- as.formula(paste(c(res_var), paste(pred_var, collapse="+"), sep="~"))
    
    #Subset df:
    df2 <- df %>% select(all_of(c(i,res_var,G_var,covariates)))
    
    #df <- na.omit(df)
    
    #Linear Regression:
    fit <- lm(formula, data = df2)
    models[[i]] <- fit
  }
  
  return(models)
}
#'*To Build Later: Partitioned R2 function:*


#'*Generalized Version of this function*
GWIS_run_generalized <- function(protID, E_df, covars_df, covars_list, feature_spec, 
                                 ordinal_list, contrast_type){
  omicDS <- GS_struct(protID) #LDLR
  
  protID <- gsub("-", "_", protID)
  protGS <- paste0(protID,"_GS")
  
  
  #Merge relevant files together:
  #'*Important distinction with colnames to check:*
  #'*Need to build function to check names first to see if appropriate then modify if needed:*
  colnames(E_df) <- gsub(" ","_", colnames(E_df))
  #If not appropriate still use:
  colnames(E_df) <- make.names(colnames(E_df))
  
  exp1 <- merge(omicDS$combo, E_df, by = "eid")
  exp2 <- merge(exp1, covars_df, by = "eid")
  
  
  E_ids <- colnames(E_df)[-1]
  
  ans <- GxE_assoc_ver2(df = exp1, E_list = E_ids, res_var = protID,
                        G_var = protGS, covariates = c(), 
                        feature_eng = feature_spec, ordinal_vars = ordinal_list, 
                        ordinal_contrast = contrast_type)
  ans2 <- GxE_assoc_ver2(df = exp2, E_list = E_ids, res_var = protID,
                         G_var = protGS, covariates = covars_list,
                         feature_eng = feature_spec,  ordinal_vars = ordinal_list, 
                         ordinal_contrast = contrast_type)
  
  res <- list(ans, ans2)
  return(res)
}



#Visualize a specific Environmental Factor and Association w/ Genetic Score:
#'*build func for this visualization later:*
#'*Get Model Output*
GWIS_run_modelOut <- function(protID, E_df, covars_df, covars_list){
  omicDS <- GS_struct(protID) #LDLR
  
  protID <- gsub("-", "_", protID)
  protGS <- paste0(protID,"_GS")
  
  
  #Merge relevant files together:
  #'*Important distinction with colnames to check:*
  #'*Need to build function to check names first to see if appropriate then modify if needed:*
  colnames(E_df) <- gsub(" ","_", colnames(E_df))
  #If not appropriate still use:
  colnames(E_df) <- make.names(colnames(E_df))
  
  exp1 <- merge(omicDS$combo, E_df, by = "eid")
  exp2 <- merge(exp1, covars_df, by = "eid")
  
  
  E_ids <- colnames(E_df)[-1]
  
  ans <- GxE_assoc_model(df = exp1, E_list = E_ids, res_var = protID,
                         G_var = protGS, covariates = c())
  ans2 <- GxE_assoc_model(df = exp2, E_list = E_ids, res_var = protID,
                          G_var = protGS, covariates = covars_list)
  
  res <- list(ans, ans2)
  return(res)
}

#Obtain newdata || for counterfactual plot:
newdata_GWIS <- function(statmodel, protID, protGS, Evar){
  df <- model.frame(statmodel)
  
  #Initialize List and add corresponding names in lm Model:
  newlist <- vector("list", length(colnames(df)) - 1)
  
  names(newlist) <- colnames(df)[!colnames(df) == protID]
  names(newlist)
  
  for(i in names(newlist)){
    if(i == protGS){
      newlist[[i]] = seq(min(df[[i]]), 
                         max(df[[i]]),
                         length.out = 10)
      
    } else if (i == Evar){
      newlist[[i]] = quantile(df[[i]], probs = c(0.1, 0.9))
      
    } else if (is.numeric(df[[i]])) {
      newlist[[i]] = median(df[[i]])
      
    } else {
      newlist[[i]] = levels(df[[i]])
      
    }
    
  }
  return(newlist)
}

#ONLY works for idx == 2 currently:
GWIS_GxEplot <- function(statmodel, idx, Evar, protID){
  protID <- gsub("-", "_", protID)
  protGS <- paste0(protID,"_GS")
  model <- statmodel[[idx]][[Evar]]
  GxEdf <- model$model
  
  #Change from quantile to line plot:
  #Oh I guess the other way to determine GxE is nonlinear interactions..??
  newlist <- newdata_GWIS(model, protID, protGS, Evar)
  newdata <- expand.grid(newlist)
  
  newdata$y_pred <- predict(model, newdata, type = "response")
  newdata$se <- predict(model, newdata, se.fit = TRUE)$se.fit  
  
  #Plot scatterplot:
  #Curves with vs without T2D:
  gg1 <- ggplot() + 
    geom_point(data = newdata[newdata$sex_f31_0_0 == "Male",], aes(x = .data[[protGS]], 
                                                                   y = y_pred, color = .data[[Evar]]))
  print(gg1)
  
  gg2 <- ggplot() + 
    geom_point(data = newdata, aes(x = .data[[protGS]], 
                                   y = y_pred, color = .data[[Evar]])) +
    facet_wrap(~ .data[["sex_f31_0_0"]])
  print(gg2)
}
