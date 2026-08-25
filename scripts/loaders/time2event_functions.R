
local({
  candidates <- c(
    Sys.getenv("HEAP_PATHS_FILE", unset = ""),
    file.path(getwd(), "workflow", "00_paths.R"),
    file.path(getwd(), "..", "workflow", "00_paths.R"),
    file.path(getwd(), "..", "..", "workflow", "00_paths.R"),
    file.path(getwd(), "..", "..", "..", "workflow", "00_paths.R")
  )
  candidates <- candidates[nzchar(candidates)]
  hit <- candidates[file.exists(candidates)][1]
  if (!is.na(hit)) source(hit)
})
#'*Functions for Obtaining Age of Event Occuring*
# Function to FIND Age of Recruitment/Death Record/Censoring/etc.:
age_of_event <- function(assessment, field, recode, df){
  date_of_INST <- paste0(field, assessment)
  age_of_INST <- paste0(recode, assessment)
  df[[date_of_INST]] <- as.Date(df[[date_of_INST]])
  df[[age_of_INST]] <- round(as.numeric((df[[date_of_INST]] - df$birth_date) /365), digits = 1)
  return(df)
}

# Function to FIND people with a specific ICD10 code
UKB_patients_ICDcode <- function(df, ICD_code){
  #Expects a dataframe that IS ALREADY multiarray handled!!
  #NEED A check for multiarray handlings later:
  
  diagICD <- list()
  
  # Output patients with ICDcode 1.) ever 2.) primary 3.) secondary
  diagICD[["All"]] <- UKB_instance0[grepl(ICD_code, UKB_instance0$diagnoses_icd10_f41270_0_0.multi), ]$eid
  diagICD[["Main"]] <- UKB_instance0[grepl(ICD_code, UKB_instance0$diagnoses_main_icd10_f41202_0_0.multi), ]$eid
  diagICD[["Secondary"]] <- UKB_instance0[grepl(ICD_code, UKB_instance0$diagnoses_secondary_icd10_f41204_0_0.multi), ]$eid
  return(diagICD)
}

# Function to FIND Date of Diagnosis of specific ICD10 code *full matching*
find_date <- function(letter, col1, col2) {
  letters <- strsplit(col1, ';')
  dates <- strsplit(col2, ';')
  
  match_index <- sapply(letters, function(x) letter %in% x)
  corresponding_date <- sapply(seq_along(letters), function(i) ifelse(match_index[i], dates[[i]][which(letters[[i]] == letter)], NA))
  
  return(corresponding_date)
}

# Function to FIND Date of Diagnosis of specific ICD10 code *partial matching*
find_date2 <- function(partial_letter, col1, col2) {
  letters <- strsplit(col1, ';')
  dates <- strsplit(col2, ';')
  
  corresponding_date <- sapply(seq_along(letters), function(i) {
    index <- grep(partial_letter, letters[[i]])
    if (length(index) > 0) dates[[i]][index] else NA
  })
  
  return(corresponding_date)
}

# Function to FIND Age of Diagnosis on merged df *df should contain at least info from Categories [2002 and 100094]*
age_of_diagnosis <- function(icd10_code, df){
  date_of_ICD <- paste0("date_of_",icd10_code) 
  df[[date_of_ICD]] <- find_date2(icd10_code, df$diagnoses_icd10_f41270_0_0.multi, 
                                  df$date_of_first_in_patient_diagnosis_icd10_f41280_0_0.multi)
  df[[date_of_ICD]] <- as.Date(df[[date_of_ICD]])
  
  age_of_ICD <- paste0("age_of_",icd10_code) 
  df[[age_of_ICD]] <- round(as.numeric((df[[date_of_ICD]] - df$birth_date)/365), digits = 1)
  
  return(df)
}





#'*FUNCTIONS in ICD10 Loader*
# Function to Load ICD10 Info:
load_icd10_matrix <- function(version, specific_category){
  #Helpful tips:
  #obtain specific_category via pathIDs
  #ex: version = "ICD10_firstoccur"; specific_category = "path_2404"
  
  if(version == "ICD10_summarydiag"){
    pathID = 2002
    load_path(projID, pathID)
    icd_df <- fast_dataloader(path, directoryInfo) 
    icd_df <- UKB_instances(icd_df, "_0_")
    icd_df <- UKB_multiarray_handle(icd_df) #'*Need to Speed up this process*
  } else if(version == "ICD10_firstoccur"){
    pathID = specific_category
    load_path(projID, pathID)
    icd_df <- fast_dataloader(path, directoryInfo)
    
    # #Setup Parallelization:
    # plan(multisession, workers = availableCores() - 1)
    # #Increase the RAM for each GLOBAL OBJECT:
    # maxSize = 2000*1024^2
    # options(future.globals.maxSize= maxSize)
    # 
    # #Load via Parallelization: *Why load all pathways if not necessary???*
    # icd10_pathIDs <- pathIDs %>% filter(grepl("First occurrences", Path)) %>% pull(Category)
    # icd10_features <- parallel_UKB_load(icd10_pathIDs)
    # icd_df <- icd10_features[[specific_category]]
  } else{
    icd_df <- NULL
  }
  return(icd_df)
}

# Function to Load Baseline Measurement Info:
load_baseline_info <- function(icd_df){
  
  #'###*BIRTH DATE of INDIVIDUALS*###
  #Category 100094; #Population Char » Baseline Char;
  pathID = 100094
  load_path(projID, pathID)
  baseline_df <- fast_dataloader(path, directoryInfo)
  baseline_df <- UKB_instances(baseline_df, "_0_") #gives age of recruitment for Instances
  baseline_df <- baseline_df %>% mutate(month_of_birth_f52_0_0 = recode(month_of_birth_f52_0_0,
                                                                        January = 1, February = 2, March = 3,
                                                                        April = 4, May = 5, June = 6,
                                                                        July = 7, August = 8, September = 9,
                                                                        October = 10, November = 11, December = 12
  ))
  
  #Obtain Year-month and ("15" as day) Birth Date 
  baseline_df$birth_date <- as.Date(with(baseline_df, paste(year_of_birth_f34_0_0,
                                                            month_of_birth_f52_0_0, 
                                                            "15", sep="-")),"%Y-%m-%d")
  
  
  #'###*Get Data For CENSORING*###
  #'*Age at Each Instance of Recruitment*
  pathID = 100024
  load_project(projID)
  load_path(projID, pathID)
  reception_df <- fast_dataloader(path, directoryInfo) 
  reception_df <- UKB_instances(reception_df, "_0_")
  
  #'*Death Registry*
  pathID = 100093
  load_project(projID)
  load_path(projID, pathID)
  death_df <- fast_dataloader(path, directoryInfo)
  death_df <- UKB_instances(death_df, "_0_") #To obtain registry for Instance 0 participants 
  death_df <- UKB_multiarray_handle(death_df)
  
  #'*Lost to Followup*
  pathID = 2
  load_project(projID)
  load_path(projID, pathID)
  censor_df <- fast_dataloader(path, directoryInfo)
  
  Time2Event_list <- list(icd_df, baseline_df, reception_df, death_df, censor_df)
  Time2Event_df <- Time2Event_list %>% reduce(full_join, by = "eid")
  return(Time2Event_df)
}

# Function to Obtain Ages to Events:
time2event_ages <- function(Time2Event_df, disease_code){
  #'*Obtain Ages of Events*
  #Assessment Center Instance 0: Age
  Time2Event_df <- age_of_event("_0_0", "date_of_attending_assessment_centre_f53", 
                                "recode_age_of_assessment", Time2Event_df)
  
  #Death: Age
  Time2Event_df <- age_of_event("_0_0", "date_of_death_f40000", 
                                "recode_age_of_death", Time2Event_df)
  
  #Censoring Left Study: Age
  Time2Event_df <- age_of_event("_0_0", "date_lost_to_follow_up_f191", 
                                "age_of_removal", Time2Event_df)
  
  #Censoring to End of Study Calculated via Modification Date and Modify UKB Category Date of ICD10 Codes
  end_of_study = as.Date("2022-10-03") #'*Double check if this is true for both dates*
  Time2Event_df[["age_of_lastfollowup"]] <- round(as.numeric((end_of_study - Time2Event_df$birth_date)/365), digits = 1)
  
  return(Time2Event_df)
}

# Function to Obtain Age of Diagnosis
ICD10_ages <- function(version, Time2Event_df, disease_code, disease_recode){
  #disease code can be ICD10 (If from summarydiag category) or
  #the name of the column in ICD (firstoccurence category)
  if(version == "ICD10_summarydiag"){
    Time2Event_df<- age_of_diagnosis(disease_code, df = Time2Event_df)
  }
  else if(version == "ICD10_firstoccur"){
    #'*Picking ICD10 code to do survival/Time-to-Event Analysis*
    # Recode into survival: Ex. E11 is Non-Insulin Dependent T2D
    Time2Event_df <- age_of_event("_0_0", disease_code,
                                  disease_recode, Time2Event_df)
  }
  return(Time2Event_df)
}


#'*FUNCTIONS for CoxPH Model*
library(survival)
coxph_obtain_stats<- function(df, metab){
# Obtain conf.int, hazard ratio, etc.
v1 <- df$conf.int[metab,c("lower .95", "upper .95")]
v2 <- df$coefficients[metab,c("exp(coef)", "se(coef)", "Pr(>|z|)")]
v3 <- df$concordance

# Output coxph results
names(metab) <- "predictor"
coxph_stats <- c(metab,v1,v2,v3)
return(coxph_stats)
}

coxph_model <- function(pred_var, res_var, df){
  #res_var <- "Surv(surv_time, T2D_status)"
  #pred_var <- c("recode_age_of_assessment_0_0", "sex_f31_0_0", metab)
  formula <- as.formula(paste(c(res_var), paste(pred_var, collapse="+"), sep="~"))
  res.cox <- coxph(formula, data = df)
  cox_df <- summary(res.cox)
  return(cox_df)
}


#'*FUNCTIONS in Model Specification RScript*
# Function to obtain survival status and time of an event:
# Reminder: age_of_lastfollowup is predicted via date of files (reason we take minimum of all values)
survival_time <- function(Time2Event_df, event_age, recode_status, recode_survtime){
  ##Survival Function
  ##Rules:
  #1.) Keep individuals if age_assessment < age_ICD10
  #2.) All individuals with ICD10 code is listed as 1; 0 otherwise
  #3.) For the individuals with censoring: (Find the appropriate date)
  T2E_df <- as.data.frame(Time2Event_df)
  
  #Rule #1 :: Keep individuals if age_assessment < age_ICD10
  T2E_df <- T2E_df[T2E_df$recode_age_of_assessment_0_0 < T2E_df[[event_age]] | is.na(T2E_df[[event_age]]), ]
  
  #Rule #2 ::
  T2E_df[[recode_status]] <- ifelse(is.na(T2E_df[[event_age]]), 0, 1)
  
  #Rule #3 ::
  censored_indv <- T2E_df[T2E_df[[recode_status]] == 0, ]$eid
  
  #remove individuals without birth date (Equiv to Indivudals who dont have a recode_age_assesment)
  T2E_df <- T2E_df[!(is.na(T2E_df$recode_age_of_assessment_0_0)),]
  
  #GET THE CORRECT AGES FOR TIME: Solution take the minimum of ages (death, diagnosis, and last followup)
  # Taking the minimum age of censored data
  T2E_df[[recode_survtime]] <- apply(T2E_df[,c(event_age,"recode_age_of_death_0_0",
                                               "age_of_removal_0_0","age_of_lastfollowup")], 1, min, na.rm = TRUE)
  
  return(T2E_df)
}

# Function to load covariate matrix based on specified UKB field ids
load_covariate_matrix <- function(fields, instance = "_0_"){
  #Dictionary with all fields and their location to specific files/folders
  UKBdict <- fread(file = heap_raw_or_legacy(paste0("allpaths_", projID, ".txt"),
    legacy_ukb_path("Data", "Paths", projID, "allpaths.txt")))
  df <- fast_dataloader_viafield(UKBdict, UKBfieldIDs = fields, directoryInfo)
  df <- UKB_instances(df, instance)
  return(df)
}

# Function to get the covariate names:
covariate_names <- function(fields, instance_time = 0){
  #Dictionary with all fields and their location to specific files/folders
  UKBdict <- fread(file = heap_raw_or_legacy(paste0("allpaths_", projID, ".txt"),
    legacy_ukb_path("Data", "Paths", projID, "allpaths.txt")))
  names <- UKBdict[UKBdict$FieldID %in% fields & UKBdict$instance == instance_time,]$descriptive_colnames
  
  return(names)
}

# Function to scale and factor specific variables:
scale_features <- function(df, scaled_fields, factor_fields, transformed_fields, transformation_function = NULL){
  df <- df %>%
    mutate(across(all_of(names(factor_fields)), as.factor))
  
  if (!is.null(transformation_function)) {
    df <- df %>%
      mutate(across(all_of(transformed_fields), transformation_function))
  }
  
  scaled_data <- df %>%
    mutate(across(all_of(scaled_fields), scale))
  
  return(scaled_data)
}

# Function for Inverse Normal Transformation:
inv_normal <- function(x){
  x <- qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))
  return(x)
}

# Functions to Run Different Versions of COXPH: 
run_coxph <- function(data, predictors, covariates, disease_status, survtime){
  model <- data.frame()
  for(p in predictors){
    #Subset dataframe with proper information:
    coxph_df <- data %>% select(all_of(c(survtime, disease_status, p, covariates)))
    coxph_df <- na.omit(coxph_df)
    
    #Run coxph model:
    full_predictor <- c(p, covariates)
    response <- paste0("Surv(",survtime,",",disease_status,")")
    
    cox_df <- coxph_model(full_predictor, response, coxph_df)
    
    coxph_res <- coxph_obtain_stats(cox_df, p)
    
    model <- rbind(model, coxph_res)
    colnames(model) <- names(coxph_res)
  }
  
  return(model)
}
run_coxph_ver2 <- function(data, predictors, covariates, disease_status, inittime, survtime){
  #Function assumes:
  #Predictors and Covariates are scaled/transformed such that the hazard ratios reflect change in 1 SD unit.
  model <- data.frame()
  for(p in predictors){
    #Subset dataframe with proper information:
    coxph_df <- data %>% select(all_of(c(survtime, inittime, disease_status, p, covariates)))
    coxph_df <- na.omit(coxph_df) #'*remove participants with missing values*
  
    #Run coxph model:
    full_predictor <- c(p, covariates)
    response <- paste0("Surv(time =", inittime, ", time2 =", survtime, ", event =", disease_status,")")
    
    cox_df <- coxph_model(full_predictor, response, coxph_df)
    
    coxph_res <- coxph_obtain_stats(cox_df, p)
    
    model <- rbind(model, coxph_res)
    colnames(model) <- names(coxph_res)
  }
  
  return(model)
}
run_coxph_parallel <- function(data, predictors, covariates, disease_status, survtime) {
  #set number of cores
  plan(multisession, workers = availableCores() - 1)  # Set the number of workers
  #Increase the RAM for each GLOBAL OBJECT:
  maxSize = 2000*1024^2
  options(future.globals.maxSize= maxSize)
  
  future_results <- future_map_dfr(predictors, ~ {
    p <- .
    # Subset dataframe with proper information:
    coxph_df <- data %>% select(all_of(c(survtime, disease_status, p, covariates)))
    
    #Scale/Standardize the predictor variable:
    #1.) log(1+x) transformation to get approx. normal
    #2.) z-score 
    #3.) hazard ratio reflect a change in 1SD unit
    coxph_df[[p]] <- scale(log1p(coxph_df[[p]]))
    
    
    # Run coxph model:
    full_predictor <- c(p, covariates)
    response <- paste0("Surv(",survtime,",",disease_status,")")
    
    cox_df <- coxph_model(full_predictor, response, coxph_df)
    
    coxph_res <- coxph_obtain_stats(cox_df, p)
    
    return(coxph_res)
  })
  
  return(future_results)
}

#'*Functions from time2event_runner.R*
#'*CHECK FOR INSTANCE ISSUES LATER!!!*
fastload_baseline_info <- function(icd_df){
  #age_of_death:  f40007
  #date of death: f40000
  #ukbiobankassessment_centre: f54
  #ageassessment_center: f21003
  #dateassessment_center: f53
  #year_of_birth: f34
  #month_of_birth: f52
  #date_lost_to_followup: f191
  #reason_lost_to_followup: f190
  UKBdict <- fread(file = heap_raw_or_legacy(paste0("allpaths_", projID, ".txt"),
    legacy_ukb_path("Data", "Paths", projID, "allpaths.txt")))
  baseline_metrics <- c(40007, 40000, 54, 21003, 53, 34, 52, 191, 190)
  baseline_df <- fast_dataloader_viafield(UKBdict, UKBfieldIDs = baseline_metrics, directoryInfo)
  baseline_df <- UKB_instances(baseline_df, "_0_")
  baseline_df <- baseline_df %>% mutate(month_of_birth_f52_0_0 = recode(month_of_birth_f52_0_0,
                                                                        January = 1, February = 2, March = 3,
                                                                        April = 4, May = 5, June = 6,
                                                                        July = 7, August = 8, September = 9,
                                                                        October = 10, November = 11, December = 12
  ))
  
  #Obtain Year-month and ("15" as day) Birth Date 
  baseline_df$birth_date <- as.Date(with(baseline_df, paste(year_of_birth_f34_0_0,
                                                            month_of_birth_f52_0_0, 
                                                            "15", sep="-")),"%Y-%m-%d")
  
  Time2Event_list <- list(icd_df, baseline_df)
  Time2Event_df <- Time2Event_list %>% reduce(full_join, by = "eid")
  return(Time2Event_df)
}
load_disease_T2E <- function(ICD10version, category, diseaseID, diseaseAGE){
  icd_df <- load_icd10_matrix(version = ICD10version, specific_category = category) #technically this can be done by field also (no need to load entire matrix)
  T2E_df <- fastload_baseline_info(icd_df)
  T2E_df <- time2event_ages(T2E_df)
  T2E_df <- ICD10_ages(version = ICD10version, T2E_df, disease_code = diseaseID,
                       disease_recode = diseaseAGE)
  T2E_df_subset <- T2E_df %>% select(all_of(c("eid","recode_age_of_assessment_0_0",
                                              "recode_age_of_death_0_0",
                                              "age_of_removal_0_0",
                                              "age_of_lastfollowup",
                                              paste0(diseaseAGE,"_0_0"))))
  return(T2E_df_subset)
}



