library(data.table)
library(tidyverse)
library(ggplot2)
library(ggpmisc)
library(ggrepel)
library(gridExtra)
library(plotly)
library(htmlwidgets)

#Addn Libraries for Visualization:
library(patchwork)

### SPECIFIC WARNINGS!!!
#'* We DO NOT save all html elements here to prevent massive files.*
#'*If selfcontained = TRUE the html files are 3.8 MB each!!!*


#Output tables for specific components:
setwd("/n/groups/patel/shakson_ukb/UK_Biobank/")

###### G v E - Lasso Version #######
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


### GET the G vs E component across specifications:
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



#'*Organize the table: to have mean and confidence intervals*
GxE_R2table <- plot_testR2_CI %>%
                    mutate(across(where(is.numeric), ~ signif(., digits = 2))) %>%
                    mutate(G = G_stats_mean,
                           G_CI = paste0("(",G_stats_ci_low,",",G_stats_ci_high,")"),
                           E = E_stats_mean,
                           E_CI = paste0("(",E_stats_ci_low,",", E_stats_ci_high,")"),
                           GxE = GxE_stats_mean,
                           GxE_CI = paste0("(",GxE_stats_ci_low,",", GxE_stats_ci_high,")")) %>%
                    select(c("ID","G","G_CI","E","E_CI","GxE","GxE_CI"))
colnames(GxE_R2table) <- c("Protein","Genetics (G)","G: 95% CI",
                      "Exposures (E)","E: 95% CI",
                      "Interactions (GxE)",
                      "GxE: 95% CI")

# Public Facing App: CSV, Parquet
fwrite(GxE_R2table, file = "/n/groups/patel/shakson_ukb/UK_Biobank/Output/App/Tables/GxE_R2table.csv")
check <- fread(file = "/n/groups/patel/shakson_ukb/UK_Biobank/Output/App/Tables/GxE_R2table.csv")


### Get Interactive plot of GvE: 
library(plotly)
library(htmlwidgets)
CovarSpec_GvEplot_interact <- function(plot_testR2_CI){
  #Highlight specific proteins with high E: R2
  subset_test <- plot_testR2_CI %>% 
    filter((E_stats_mean > G_stats_mean) & (E_stats_mean > 0.09))
  
  #Cutoff for y-axis:
  y_cutoff = 0.175
  
  #Without Error Bars:
  p1 <- plot_testR2_CI  %>% 
    mutate(across(where(is.numeric), ~ signif(., digits = 2))) %>%
    ggplot(aes(x = G_stats_mean, y = E_stats_mean, label = ID)) + 
    geom_point(aes(text = paste("<br>Prot:", ID,
                                "<br>G:", paste0(G_stats_mean, ": (",
                                                 G_stats_ci_low,",",G_stats_ci_high,")"),
                                "<br>E:", paste0(E_stats_mean,": (",
                                                 E_stats_ci_low,",",E_stats_ci_high,")")
                                )
                   )
               ) +
    geom_abline(slope=1, intercept = 0) +
    geom_ribbon(aes(ymin = 0, ymax = ifelse(G_stats_mean <= y_cutoff, G_stats_mean, y_cutoff), fill = "Genetics Driven"), alpha = 0.2) +
    geom_ribbon(aes(ymin = ifelse(G_stats_mean <= y_cutoff, G_stats_mean, y_cutoff), ymax = y_cutoff, fill = "Environment Driven"), alpha = 0.2) +
    scale_fill_manual(values = c("Genetics Driven" = "blue", "Environment Driven" = "green")) +
    ylim(c(0, y_cutoff)) +
    #geom_text_repel(data=subset_test, aes(x = G_stats_mean,y = E_stats_mean,label= ID),
    #                size = 4) +
    xlab("G: PGS R<sup>2</sup>") +
    ylab("E: PXS R<sup>2</sup>") +
    theme_minimal() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    guides(fill=guide_legend(title="Region")) +
    ggtitle("Test Set: Partitioned R<sup>2</sup>") +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "bottom") 

  gginteract <- ggplotly(p1, tooltip = "text")
  
  saveWidget(gginteract, 
             file = paste0("./Output/App/Interactive/A1/GvEplot.html"))
  
}
CovarSpec_GvEplot_interact(plot_testR2_CI)



#### For specific E categories: #####
categories <- c("Alcohol","Diet_Weekly","Smoking","Exercise_MET",
                "Exercise_Freq","Internet_Usage","Deprivation_Indices","Vitamins")

# purrr::map2 iterate over the list created to generate one dataframe:
R2testCat <- map2_dfr(HEAPres, names(HEAPres), ~ .x$R2testCat %>%
                        select(all_of(c("ID", categories))) %>%
                        mutate(cID = .y))

R2testCat_CI <- R2testCat %>%
  group_by(ID, cID) %>%
  summarise(across(c(categories), mean, .names = "{col}")) %>%  #Mean of each Specification
  mutate(Exercise = Exercise_Freq + Exercise_MET) %>%
  mutate(across(c(categories,"Exercise"), mean, .names = "mean_{col}")) %>% #Mean across Specification
  mutate(across(c(categories,"Exercise"), 
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


plot_testR2Cat_CI <- R2testCat_CI %>%
                        select(-c("cID","Alcohol","Diet_Weekly","Smoking","Exercise_MET",
                                  "Exercise_Freq","Exercise","Internet_Usage",
                                  "Deprivation_Indices","Vitamins")) %>%
                        unique()


GxE_Cat_R2table <- plot_testR2Cat_CI %>%
                        mutate(across(where(is.numeric), ~ signif(., digits = 2))) %>%
                        mutate(Alcohol = mean_Alcohol,
                               Alcohol_CI = paste0("(",mean_Alcohol_ci_low,",",mean_Alcohol_ci_high,")"),
                               Diet = mean_Diet_Weekly,
                               Diet_CI = paste0("(",mean_Diet_Weekly_ci_low,",",mean_Diet_Weekly_ci_high,")"),
                               Smoking = mean_Smoking,
                               Smoking_CI = paste0("(",mean_Smoking_ci_low,",",mean_Smoking_ci_high,")"),
                               Exercise = mean_Exercise,
                               Exercise_CI = paste0("(",mean_Exercise_ci_low,",",mean_Exercise_ci_high,")"),
                               InternetUse = mean_Internet_Usage,
                               InternetUse_CI = paste0("(",mean_Internet_Usage_ci_low,",",mean_Internet_Usage_ci_high,")"),
                               DeprivationIndices = mean_Deprivation_Indices,
                               DeprivationIndices_CI = paste0("(",mean_Deprivation_Indices_ci_low,",",mean_Deprivation_Indices_ci_high,")"),
                               Vitamins = mean_Vitamins,
                               Vitamins_CI = paste0("(",mean_Vitamins_ci_low,",",mean_Vitamins_ci_high,")"),
                               ) %>%
                        select(ID,Alcohol,Alcohol_CI,
                               Diet,Diet_CI,
                               Smoking,Smoking_CI,
                               Exercise,Exercise_CI,
                               InternetUse,InternetUse_CI,
                               DeprivationIndices,DeprivationIndices_CI,
                               Vitamins,Vitamins_CI)

fwrite(GxE_Cat_R2table, file = "/n/groups/patel/shakson_ukb/UK_Biobank/Output/App/Tables/GxE_Cat_R2table.csv")
check2 <- fread(file = "/n/groups/patel/shakson_ukb/UK_Biobank/Output/App/Tables/GxE_Cat_R2table.csv")



#### For covariates:


##### MEDIATION RESULTS #####
library(pbapply)
library(dplyr)
library(purrr)
#'*Get the Mediation Results*
# Load Files w/ Progress Bar:
loadMDAssoc <- function(covarType){
  error_files <- list()  # Track files w/ errors
  stat_all <- list()
  
  # Use pblapply for progress bar and iteration
  for(i in 1:1000){
    tryCatch({
      load <- fread(file=paste0("/n/groups/patel/shakson_ukb/UK_Biobank/Data/Parallel/Prot_PGSPXS_mediation/",
                                covarType,"/MDres_",i,".txt"), 
                    select = 2:16)
      
      stat_all[[i]] <- load
      
    }, error = function(e) {
      message(paste("Error reading file:", i))
      error_files <<- append(error_files, i)
    })
  }
  
  
  stat_all <- do.call("rbind", stat_all)
  stat_all <- stat_all %>%
    mutate(E_HRi = exp(`Exposure Indirect Effect`),
           G_HRi = exp(`Genetic Indirect Effect`),
           delta_HRi = exp(abs(`Exposure Indirect Effect` - `Genetic Indirect Effect`)),
           total_HRi = exp(`Exposure Indirect Effect` + `Genetic Indirect Effect`),
           Eprop_mediated = abs(`Exposure Indirect Effect`)/(abs(`Exposure Indirect Effect`) + abs(`Exposure Direct Effect`)),
           Gprop_mediated = abs(`Genetic Indirect Effect`)/(abs(`Genetic Indirect Effect`) + abs(`Genetic Direct Effect`)))
  stat_all$CovarSpec <- covarType
  stat_all$HR_upper95 <- as.numeric(stat_all$HR_upper95)
  
  
  return(stat_all)
}

# Load results:
Type1 <- loadMDAssoc(covarType = "Type1")
Type2 <- loadMDAssoc(covarType = "Type2")
Type3 <- loadMDAssoc(covarType = "Type3")
Type4 <- loadMDAssoc(covarType = "Type4")
Type5 <- loadMDAssoc(covarType = "Type5")

AllCovarSpec <- list(Type1, Type2, Type3, Type4, Type5) %>% reduce(full_join)

# FINAL ANALYSIS 

# Obtain OLS slope across Specifications:
oSlope <- function(MDres){
  
  # Obtain GEM: Expected Indirect Effect Slope: log(E/G)
  protContext <- pblapply(unique(MDres$ID), function(p){
    df <- MDres %>% filter(ID == p) 
    
    fit <- lm(`Exposure Indirect Effect` ~ `Genetic Indirect Effect`,
              data = df)
    
    stats <- summary(fit)$coefficients
    
    maxE <- max(abs(df$`Exposure Indirect Effect`))
    maxG <- max(abs(df$`Genetic Indirect Effect`))
    
    df <- MDres %>%
      filter(ID == p) %>%
      filter(!((HR_lower95 < 1 & HR_upper95 > 1) | (HR_lower95 > 1 & HR_upper95 < 1)))
    
    num_diseases <- length(unique(df$DZid))
    
    protStats <- data.frame(
      ID = p,
      Estimate = stats[2,"Estimate"],
      Std.Error = stats[2,"Std. Error"],
      Pvalue = stats[2,"Pr(>|t|)"],
      max.E = maxE,
      max.G = maxG,
      Estimate_Norm = stats[2,"Estimate"] * (maxE/maxG),
      NumDiseases = num_diseases
    )
    
    #Below return the slope of the above regression
    return(protStats) 
  })
  protContext <- do.call(rbind, protContext)
  
  return(protContext)
}
Type1_o <- oSlope(Type1)
Type2_o <- oSlope(Type2)
Type3_o <- oSlope(Type3)
Type4_o <- oSlope(Type4)
Type5_o <- oSlope(Type5)

#Merge results together:
Type1_o$CovarType <- "Type1"
Type2_o$CovarType <- "Type2"
Type3_o$CovarType <- "Type3"
Type4_o$CovarType <- "Type4"
Type5_o$CovarType <- "Type5"


MDres <- list(Type1_o, Type2_o, Type3_o, Type4_o, Type5_o)
names(MDres) <- c("Type1", "Type2", "Type3", "Type4", "Type5")



#'*GET interactive plot of GEM - for covarspec 5*

GEMweb <- Type5_o %>%
  mutate(Estimate = log(Estimate)) %>%
  na.omit()

p2 <- GEMweb %>%
  mutate(across(where(is.numeric), ~ signif(., digits = 3))) %>%
  ggplot(aes(x = NumDiseases, y = Estimate)) +
  geom_point(aes(text = paste("<br>Prot:", ID,
                              "<br>GEM:", paste0(Estimate)))) +
  geom_ribbon(aes(ymin =min(Estimate)*1.2, ymax = 0, fill = "Genetics Driven"), alpha = 0.2) +
  geom_ribbon(aes(ymin = 0, ymax = max(Estimate)*1.2, fill = "Environment Driven"), alpha = 0.2) +
  scale_fill_manual(values = c("Genetics Driven" = "blue", "Environment Driven" = "green")) +
  xlab("Number of Significant Disease Hits") +
  ylab("GEM") +
  theme_minimal() +
  guides(fill=guide_legend(title="Region")) +
  theme(legend.position = "bottom")
p2

gginteract <- ggplotly(p2, tooltip = "text")

saveWidget(gginteract, 
           file = paste0("./Output/App/Interactive/A1/GEMplot.html"))


#'*Save tables for the website*
#'The GEM score can just be interactive and then create a file to download:

#Download the GEM score:
GEMdownload <- GEMweb %>%
                  select(c("ID","Estimate"))
colnames(GEMdownload) <- c("Protein","GEM")
fwrite(GEMdownload, 
       file = "/n/groups/patel/shakson_ukb/UK_Biobank/Output/App/Tables/GEMdownload.csv")


hist(AllCovarSpec$R2)
summary(AllCovarSpec$R2)
quantile(AllCovarSpec$R2, c(.01, .05, .95, 0.99)) 



#Download the Type 5 specification disease mediation results:
#Only show proteins with valid GEM!
MDweb <- AllCovarSpec %>%
              filter(CovarSpec == "Type5") %>%
              filter(R2 > 0.045) %>%
              filter(if_all(where(is.numeric), is.finite)) %>%
              filter(ID %in% GEMweb$ID)
length(unique(MDweb$ID))

MDweb_fin <- MDweb %>%
                mutate(Indirect_G = exp(`Genetic Indirect Effect`),
                       Indirect_E = exp(`Exposure Indirect Effect`),
                       Direct_G = exp(`Genetic Direct Effect`),
                       Direct_E = exp(`Exposure Direct Effect`) ,
                       Total_Direct = exp(`Total Direct Effect`),
                       Total_Indirect = exp(`Total Indirect Effect`)) %>%
                select(c("ID","DZid","Indirect_G","Indirect_E",
                         "Direct_G","Direct_E","Total_Direct",
                         "Total_Indirect","Eprop_mediated",
                         "Gprop_mediated"))

colnames(MDweb_fin) <- c("Protein","Disease Code","G: Indirect",
                         "E: Indirect", "G: Direct",
                         "E: Direct","Total Direct",
                         "Total Indirect","E: Prop Mediate",
                         "G: Prop Mediate")

MDweb_fin <- MDweb_fin %>%
                mutate(across(where(is.numeric), ~ signif(., digits = 3)))

fwrite(MDweb_fin, 
       file = "/n/groups/patel/shakson_ukb/UK_Biobank/Output/App/Tables/MediationResults.csv")


##### Intervention results #####
#'*Get the Interventions*

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
library(readxl)


#'*Load Univariate UKB HEAP results*
setwd("/n/groups/patel/shakson_ukb/UK_Biobank/")
file.exists("RScripts/Pure_StatGen/Prot_ExPGS/resultsProtDS.R")
source("RScripts/Pure_StatGen/Prot_ExPGS/resultsProtDS.R")


#Intervention Studies:

#Compare results to:
#'*Exercise - JCI paper: https://pmc.ncbi.nlm.nih.gov/articles/PMC10132160/#sec13*
#'*GLP1 agonist treatment: https://doi.org/10.1038/s41591-024-03355-2*

# Read: first sheet w/ header start from line 3
file_path <- "/n/groups/patel/shakson_ukb/Motrpac/Related_Data/jciinsight_prot.xlsx"
data <- read_excel(file_path, sheet = excel_sheets(file_path)[1], skip = 2)
data$se <- data$`log(10) Fold Change`/data$`t-statistic`


# Read: first sheet w/ header start from line 3
file_path <- "/n/groups/patel/shakson_ukb/Motrpac/Related_Data/GLP1_proteomics.xlsx"
STEP1 <- read_excel(file_path, sheet = excel_sheets(file_path)[2], skip = 0)
STEP2 <- read_excel(file_path, sheet = excel_sheets(file_path)[3], skip = 0)

#### General Correlations 
UKB <- CovarSpecList[["Type6"]]$test[[1]] %>% 
  filter(`Pr(>|t|)` < 0.05/n()) %>%
  select(ID, omicID, Estimate) %>%
  pivot_wider(names_from = ID, values_from = Estimate)

#Number of entries in a columns must be greater than 3 (for correlation estimates)
UKBv2 <- UKB %>%
  select(where(~sum(!is.na(.)) >= 3))

#Merge the above dataframe with Intervention studies:
colnames(UKBv2)[which(names(UKBv2) == "omicID")] <- "EntrezGeneSymbol"

HERITAGE <- data %>% 
  filter(`False Discovery Rate (q-value)` < 0.05) %>%
  rename(HERITAGE_effect = `log(10) Fold Change`) %>%
  select(c(EntrezGeneSymbol, HERITAGE_effect))

GLP1_STEP1 <- STEP1 %>% 
  filter(qvalue < 0.05) %>% 
  group_by(EntrezGeneSymbol) %>%
  mutate(GLP1_effect1 = mean(effect_size),
         GLP1_se1 = max(std_error),
         GLP1_qvalue1 = max(qvalue)) %>%
  select(c(EntrezGeneSymbol, GLP1_effect1)) %>%
  unique()

GLP1_STEP2 <- STEP2 %>% 
  filter(qvalue < 0.05) %>% 
  group_by(EntrezGeneSymbol) %>%
  mutate(GLP1_effect2 = mean(effect_size),
         GLP1_se2 = max(std_error),
         GLP1_qvalue2 = max(qvalue)) %>%
  select(c(EntrezGeneSymbol, GLP1_effect2)) %>%
  unique()


UKBint <- list(UKBv2, HERITAGE, GLP1_STEP1, GLP1_STEP2) %>% reduce(full_join)

UKBint_cor <- as.data.frame(cor(UKBint[,-1], use = "pairwise.complete.obs"))

UKBint_cor <- UKBint_cor %>%
  select(c(HERITAGE_effect, GLP1_effect1,
           GLP1_effect2)) %>%
  na.omit()

library(psych)
UKBint_corv2 <- corr.test(UKBint[,-1], adjust = "none")
UKBint_pval <- as.data.frame(UKBint_corv2$p) %>%
  select(c(HERITAGE_effect, GLP1_effect1,
           GLP1_effect2)) %>%
  na.omit()

UKBint_pvalv2 <- sapply(UKBint_pval, function(x){
  p.adjust(x, method="BH")
})
UKBint_pvalv2 <- as.data.frame(UKBint_pvalv2)
UKBint_pvalv2$eID <- rownames(UKBint_pval)


UKBint_pvalv2 <- UKBint_pvalv2 %>%
  pivot_longer(cols = !eID,
               names_to = "intervention", values_to = "pval_adjust") %>%
  filter(!(eID %in% c("HERITAGE_effect", 
                      "GLP1_effect1",
                      "GLP1_effect2")))



UKBint_corval <- as.data.frame(UKBint_corv2$r)
UKBint_corval$eID <- rownames(UKBint_corval)

UKBint_corval <- UKBint_corval %>%
  select(c(eID, HERITAGE_effect, GLP1_effect1,
           GLP1_effect2)) %>%
  pivot_longer(cols = !eID,
               names_to = "intervention", values_to = "cor") %>%
  filter(!(eID %in% c("HERITAGE_effect", 
                      "GLP1_effect1",
                      "GLP1_effect2")))

UKBint_df <- list(UKBint_pvalv2,UKBint_corval) %>% 
  reduce(full_join)

selectIDs <- UKBint_df %>%
  filter(pval_adjust < 0.05) %>%
  pull(eID)

UKBint_df <- UKBint_df %>% filter(eID %in% selectIDs)

UKBint_dfcor <- UKBint_df %>%
  select(eID, intervention, cor) %>%
  pivot_wider(names_from = intervention, values_from = cor) %>%
  column_to_rownames(var = "eID")

UKBint_dfpval <- UKBint_df %>%
  select(eID, intervention, pval_adjust) %>%
  pivot_wider(names_from = intervention, values_from = pval_adjust) %>%
  column_to_rownames(var = "eID")

#SAVE results:
#write.table(UKBint_dfcor, file = "/n/groups/patel/shakson_ukb/UK_Biobank/Output/SupplementaryTables/UKB_INT_cor.txt")
#write.table(UKBint_dfpval, file = "/n/groups/patel/shakson_ukb/UK_Biobank/Output/SupplementaryTables/UKB_INT_pval.txt")

#/n/groups/patel/shakson_ukb/UK_Biobank/Output/SupplementaryTables/UKB_INT_cor.txt
#/n/groups/patel/shakson_ukb/UK_Biobank/Output/SupplementaryTables/UKB_INT_pval.txt

###### ScatterPlot to Compare Effects across UKB and Interventions 
library(ggplot2)
library(ggpubr)
library(ggpmisc)
library(ggrepel)
library(plotly)
library(htmlwidgets)

# Function to create plot for each study
create_plot_interact <- function(df, x_col, y_col, eName, effect_name, cor_val, p_val, se_col) {
  
  p <- df %>%
    filter(!is.na(Estimate)) %>%
    mutate(across(where(is.numeric), ~ signif(., digits = 2))) %>%
    ggplot(aes(x = !!sym(x_col), y = !!sym(y_col))) +
    geom_point(aes(text = paste("<br>Prot:", EntrezGeneSymbol,
                                "<br>Exposure:", paste0(!!sym(x_col),"(",signif(!!sym(x_col) - 1.96 * `Std. Error`, digits = 2),
                                                        ",",
                                                        signif(!!sym(x_col) + 1.96 * `Std. Error`,digits = 2),")"),
                                "<br>Intervention:", paste0(!!sym(y_col),"(",signif(!!sym(y_col) - 1.96 * !!sym(se_col), digits = 2),
                                                            ",",
                                                            signif(!!sym(y_col) + 1.96 * !!sym(se_col), digits = 2),")")
                                ))) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray") + 
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray") + 
    # Use dynamic `se_col` for error bars
    geom_linerange(aes(ymin = !!sym(y_col) - 1.96 * !!sym(se_col), ymax = !!sym(y_col) + 1.96 * !!sym(se_col)), 
                   color = "black", size = 0.5, alpha = 0.5) +
    geom_linerange(aes(xmin = Estimate - 1.96 * `Std. Error`, xmax = Estimate + 1.96 * `Std. Error`), 
                   color = "black", size = 0.5, alpha = 0.5) +
    theme_minimal() +
    labs(x = paste0("Beta (", eName, ") - UKB"), y = paste0("Effect Size - ", effect_name)) + 
    geom_text_repel(aes(label = EntrezGeneSymbol), color = "red", size = 3, max.overlaps = 10)
  
  # Get axis limits
  plot_build <- ggplot_build(p)
  x_range <- plot_build$layout$panel_scales_x[[1]]$range$range
  y_range <- plot_build$layout$panel_scales_y[[1]]$range$range
  
  # Dynamically position the correlation label in the top-left corner
  x_pos <- x_range[1] + 0.15 * (x_range[2] - x_range[1])  # 15% from the left
  y_pos <- y_range[2] - 0.05 * (y_range[2] - y_range[1])  # 5% from the top
  
  # Add the annotation
  p1 <- p + annotate("text",
                     x = x_pos, 
                     y = y_pos, 
                     label = paste0("R = ", signif(cor_val, 2), ", p = ", signif(p_val, 2)), 
                     color = "black")
  
  # lapply(c("png", "svg", "pdf"), 
  #        function(fmt) ggsave(paste0("./Figures/HEAP/F5/",
  #                                    make.names(eName),"_", make.names(effect_name), ".", fmt), 
  #                             plot = p1, width = 6, height = 4, dpi = 500))
  
  return(p1)
}

# Get correlation and p-value for eID
get_correlation_pvalue <- function(df, eID, study) {
  cor_pval <- df %>%
    filter(eFeat %in% eID & Study == study) %>%
    select(cor, pval_adjust) %>%
    slice(1)
  return(list(cor = cor_pval$cor, pval = cor_pval$pval_adjust))
}

# Save HTML efficiently:
savePlotly <- function(gg, html_file){
  gginteract <- ggplotly(gg, tooltip = "text") 
  #%>% partial_bundle(local = FALSE)
  
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
}

# Function to create and save scatterplots
IntCor_Plot <- function(eID, eName, sID){
  UKBspec <- CovarSpecList[["Type6"]]$test[[1]] %>% 
    filter((ID == eID) & (`Pr(>|t|)` < 0.05/n())) %>%
    select(ID, omicID, Estimate, `Std. Error`)
  colnames(UKBspec)[which(names(UKBspec) == "omicID")] <- "EntrezGeneSymbol"
  
  HERITAGE <- data %>% 
    filter(`False Discovery Rate (q-value)` < 0.05) %>%
    mutate(HERITAGE_se = `log(10) Fold Change`/`t-statistic`) %>%
    rename(HERITAGE_effect = `log(10) Fold Change`) %>%
    select(c(EntrezGeneSymbol, HERITAGE_effect, HERITAGE_se))
  
  GLP1_STEP1 <- STEP1 %>% 
    filter(qvalue < 0.05) %>% 
    group_by(EntrezGeneSymbol) %>%
    mutate(GLP1_effect1 = mean(effect_size),
           GLP1_se1 = max(std_error)) %>%
    select(c(EntrezGeneSymbol, GLP1_effect1, GLP1_se1)) %>%
    unique()
  
  GLP1_STEP2 <- STEP2 %>% 
    filter(qvalue < 0.05) %>% 
    group_by(EntrezGeneSymbol) %>%
    mutate(GLP1_effect2 = mean(effect_size),
           GLP1_se2 = max(std_error)) %>%
    select(c(EntrezGeneSymbol, GLP1_effect2, GLP1_se2)) %>%
    unique()
  
  UKBscatter <- list(UKBspec, HERITAGE, GLP1_STEP1, GLP1_STEP2) %>% reduce(full_join)
  
  
  
  # Prepare the data for plots
  colnames(UKBint_df) <- c("eFeat","Study","pval_adjust","cor")
  corVal_HERITAGE <- get_correlation_pvalue(UKBint_df, eID, "HERITAGE_effect")
  corVal_GLP1_STEP1 <- get_correlation_pvalue(UKBint_df, eID, "GLP1_effect1")
  corVal_GLP1_STEP2 <- get_correlation_pvalue(UKBint_df, eID, "GLP1_effect2")
  
  # Specify standard error column for each study (adapt this based on the dataset)
  se_col_HERITAGE <- "HERITAGE_se"
  se_col_GLP1_STEP1 <- "GLP1_se1"
  se_col_GLP1_STEP2 <- "GLP1_se2"
  
  # Create plots
  gg1 <- create_plot_interact(UKBscatter, "Estimate", "HERITAGE_effect", eName, "HERITAGE", corVal_HERITAGE$cor, corVal_HERITAGE$pval, se_col_HERITAGE)
  gg2 <- create_plot_interact(UKBscatter, "Estimate", "GLP1_effect1", eName, "GLP1 STEP1", corVal_GLP1_STEP1$cor, corVal_GLP1_STEP1$pval, se_col_GLP1_STEP1)
  gg3 <- create_plot_interact(UKBscatter, "Estimate", "GLP1_effect2", eName, "GLP1 STEP2", corVal_GLP1_STEP2$cor, corVal_GLP1_STEP2$pval, se_col_GLP1_STEP2)
  
  
  #Individual Plots:

  #Combine plots into 1 interactive plot:
  #gginteract1 <- ggplotly(gg1, tooltip = "text")
  savePlotly(gg1, paste0("./Output/App/Interactive/A3/",
                                 sID,"_","HERITAGE",".html"))
  
  
  # saveWidget(gginteract1,
  #            file = paste0("./Output/App/Interactive/A3/",
  #            sID,"_","HERITAGE",".html"),
  #            selfcontained = FALSE)
  #gginteract2 <- ggplotly(gg2, tooltip = "text")
  
  savePlotly(gg2, paste0("./Output/App/Interactive/A3/",
                           sID,"_","GLP1_STEP1",".html"))
  
  # saveWidget(gginteract2,
  #            file = paste0("./Output/App/Interactive/A3/",
  #                          sID,"_","GLP1_STEP1",".html"),
  #            selfcontained = FALSE)
  #gginteract3 <- ggplotly(gg3, tooltip = "text")
  
  savePlotly(gg3, paste0("./Output/App/Interactive/A3/",
                           sID,"_","GLP1_STEP2",".html"))
  
  # saveWidget(gginteract3,
  #            file = paste0("./Output/App/Interactive/A3/",
  #                          sID,"_","GLP1_STEP2",".html"),
  #            selfcontained = FALSE)
}

#For loop through these exposures:
ExposureDict <- data.frame(
  origID = rownames(UKBint_dfcor),
  sID =  paste0("E",1:length(rownames(UKBint_dfcor))),
  name = c("White Bread Intake",
           "Exercise: Swimming, Cycling, etc.",
           "Incr. Alcohol Intake (vs. 10 yrs ago)",
           "Cereal Intake",
           "Former Daily Smoking",
           "Fresh Fruit Intake",
           "Smoked in Lifetime (Yes)",
           "Never Smoked",
           "Beef Intake 1x/WK",
           "MET min/WK of Vigorous Activity",
           "# Days/WK of Vigorous Activity",
           "Exercise: Strenuous Sports")
)

fwrite(ExposureDict, file = "/n/groups/patel/shakson_ukb/UK_Biobank/Output/App/Interactive/A3/exposure_interv.csv")

lapply(1:nrow(ExposureDict), function(x){
  IntCor_Plot(eID = ExposureDict[x,"origID"],
              eName = ExposureDict[x,"name"],
              sID = ExposureDict[x,"sID"])
})


### Create plots across relevant significant exposure and intervention study pairs:






#### Obtain HTMLs of Univariate Associations #####
#Use most specified model:

#  LOAD UNIVARIATE SUMMARY STATS -
# FUNCTION: Load Files w/ Progress Bar:
loadUniAssoc <- function(covarType){
  error_files <- list()  # Track files w/ errors
  stat_all <- list()
  train_data_list <- list()
  test_data_list <- list()
  
  # For loop through stored Univariate Associations in batches (400 batches here)
  for(i in 1:400){
    tryCatch({
      load <- readRDS(paste0("/n/groups/patel/shakson_ukb/UK_Biobank/Data/Parallel/Prot_PGSPXS_univar/",covarType,"/univar_assoc_",i,".rds"))
      
      train_data_list[[i]] <- load$train
      test_data_list[[i]]  <- load$test
      
    }, error = function(e) {
      message(paste("Error reading file:", i))
      error_files <<- append(error_files, i)
    })
  }
  
  # Store results
  stat_all[["train"]] <- train_data_list
  stat_all[["test"]] <- test_data_list
  
  return(stat_all)
}

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

#'*LOADING Associations Across Specifications*
loadAssoc <- function(covarType){
  # Load specific Covar Spec:
  assoc <- loadUniAssoc(covarType)
  
  # Aggregating Train and Test Results:
  assoc$train <- combine_dataframes(assoc$train)
  assoc$test <- combine_dataframes(assoc$test)
  
  return(assoc)
}

Type1 <- loadAssoc("Type1")
Type2 <- loadAssoc("Type2")
Type3 <- loadAssoc("Type3")
Type4 <- loadAssoc("Type4")
Type5 <- loadAssoc("Type5")
Type6 <- loadAssoc("Type6")
Type7 <- loadAssoc("Type7")

# CovarSpec: All ASSOCIATIONS DataStructure:
CovarSpecList <- list(Type1, Type2, Type3, Type4, Type5, Type6)
names(CovarSpecList) <- c("Type1","Type2","Type3","Type4","Type5","Type6")

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


#'*SINGLE Omic Signatures:*
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
#saveOmicSignaturesApp("IGFBP2",Eall, "Type6","E Associations")


library(pbapply)
library(dplyr)
library(purrr)

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

#Figured out color scheme for categories:
#library(scales)
#hex <- hue_pal()(8)

runAllInteractiveAssoc("Type1")
runAllInteractiveAssoc("Type2")
runAllInteractiveAssoc("Type3")
runAllInteractiveAssoc("Type4")
runAllInteractiveAssoc("Type5")
runAllInteractiveAssoc("Type6") #Expected Runtime is 1 Hour...: took 1 hour 10 minutes total



##### Get Tables for Univariate Associations: E and GxE #####

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

# CovarSpec: All Significant Assocations Datastructure
CovarSpecAssocList <- lapply(CovarSpecList, function(x){
  Assoc <- list()
  E <- replication_stats_final(x$train[[1]], x$test[[1]]) #E
  GxE <- replication_stats_final(x$train[[2]], x$test[[2]]) #GxE
  
  Assoc[["E"]] <- E
  Assoc[["GxE"]] <- GxE
  
  return(Assoc)
})


# E Associations:
fwrite(CovarSpecAssocList$Type1$E$sigBOTH,
       file = "/n/groups/patel/shakson_ukb/UK_Biobank/Output/SupplementaryTables/CovarSpec1EassocReplicated.txt")
fwrite(CovarSpecAssocList$Type2$E$sigBOTH,
       file = "/n/groups/patel/shakson_ukb/UK_Biobank/Output/SupplementaryTables/CovarSpec2EassocReplicated.txt")
fwrite(CovarSpecAssocList$Type3$E$sigBOTH,
       file = "/n/groups/patel/shakson_ukb/UK_Biobank/Output/SupplementaryTables/CovarSpec3EassocReplicated.txt")
fwrite(CovarSpecAssocList$Type4$E$sigBOTH,
       file = "/n/groups/patel/shakson_ukb/UK_Biobank/Output/SupplementaryTables/CovarSpec4EassocReplicated.txt")
fwrite(CovarSpecAssocList$Type5$E$sigBOTH,
       file = "/n/groups/patel/shakson_ukb/UK_Biobank/Output/SupplementaryTables/CovarSpec5EassocReplicated.txt")
fwrite(CovarSpecAssocList$Type6$E$sigBOTH,
       file = "/n/groups/patel/shakson_ukb/UK_Biobank/Output/SupplementaryTables/CovarSpec6EassocReplicated.txt")

# GxE Associations:
fwrite(CovarSpecAssocList$Type1$GxE$sigBOTH,
       file = "/n/groups/patel/shakson_ukb/UK_Biobank/Output/SupplementaryTables/CovarSpec1GxEassocReplicated.txt")
fwrite(CovarSpecAssocList$Type2$GxE$sigBOTH,
       file = "/n/groups/patel/shakson_ukb/UK_Biobank/Output/SupplementaryTables/CovarSpec2GxEassocReplicated.txt")
fwrite(CovarSpecAssocList$Type3$GxE$sigBOTH,
       file = "/n/groups/patel/shakson_ukb/UK_Biobank/Output/SupplementaryTables/CovarSpec3GxEassocReplicated.txt")
fwrite(CovarSpecAssocList$Type4$GxE$sigBOTH,
       file = "/n/groups/patel/shakson_ukb/UK_Biobank/Output/SupplementaryTables/CovarSpec4GxEassocReplicated.txt")
fwrite(CovarSpecAssocList$Type5$GxE$sigBOTH,
       file = "/n/groups/patel/shakson_ukb/UK_Biobank/Output/SupplementaryTables/CovarSpec5GxEassocReplicated.txt")
fwrite(CovarSpecAssocList$Type6$GxE$sigBOTH,
       file = "/n/groups/patel/shakson_ukb/UK_Biobank/Output/SupplementaryTables/CovarSpec6GxEassocReplicated.txt")


fwrite(CovarSpecAssocList$Type6$E$sigBOTH,
       file = "/n/groups/patel/shakson_ukb/UK_Biobank/Output/App/Tables/Signif_Eassoc.txt")
fwrite(CovarSpecAssocList$Type6$GxE$sigBOTH,
       file = "/n/groups/patel/shakson_ukb/UK_Biobank/Output/App/Tables/Signif_GxEassoc.txt")




#### Get Plots for Significant GxE Associations ####
#Source univarInteract Script:
source("/n/groups/patel/shakson_ukb/UK_Biobank/RScripts/Pure_StatGen/Prot_ExPGS/runProt_PGS_univarInteract_v2.R")

#'*To Modify Later --> For # of days,etc. that got converted into z-score.*
#'*Discretize the above by days and find the respective quintiles.*

#Run specific associations:
saveGWISplot_app <- function(o,e,c,n,l=c(),p,f){
  gg <- runGWISplot(omic = o,
                    exposure = e,
                    CType = c,
                    ename = n,
                    levels = l,
                    plotname = p)
  ggsave(paste0("./Output/App/Interactive/P1/",paste0(o,"_",f),".png"), 
         plot = gg, width = 5, height = 4, dpi = 500)
}

#Find top E hits that REPLICATE to show:
GxEappType5 <- fread("/n/groups/patel/shakson_ukb/UK_Biobank/Output/SupplementaryTables/CovarSpec5GxEassocReplicated.txt")
GxEappType6 <- fread("/n/groups/patel/shakson_ukb/UK_Biobank/Output/SupplementaryTables/CovarSpec6GxEassocReplicated.txt")


GxEids <- GxEappType6 %>%
              select(c("Eid_train","omicID"))

IDmatching <- data.frame(
   ID = c("alcohol_intake_frequency_f1558_0_0",
          "plays_computer_games_f2237_0_0",                                                                                     
          "types_of_transport_used_excluding_work_f6162_0_0.multi_Public_transport",                                            
          "alcohol_drinker_status_f20117_0_0_Current",                                                                          
          "alcohol_drinker_status_f20117_0_0_Never",                                                                            
          "bread_intake_f1438_0_0",                                                                                             
          "coffee_intake_f1498_0_0",                                                                                            
          "pork_intake_f1389_0_0",                                                                                              
          "time_spent_driving_f1090_0_0",                                                                                       
          "current_tobacco_smoking_f1239_0_0_No",                                                                               
          "current_tobacco_smoking_f1239_0_0_Yes._on_most_or_all_days",                                                         
          "smoking_status_f20116_0_0_Current",                                                                                  
          "types_of_physical_activity_in_last_4_weeks_f6164_0_0.multi_Heavy_DIY_.eg._weeding._lawn_mowing._carpentry._digging.",
          "types_of_transport_used_excluding_work_f6162_0_0.multi_Car.motor_vehicle",                                           
          "water_intake_f1528_0_0",                                                                                             
          "cheese_intake_f1408_0_0"),
   Name = c("Alcohol Intake Freq.",
            "Plays Computer Games",
            "Uses Public Transport",
            "Current Alcohol Drinker",
            "Never Drinks Alcohol",
            "Bread Intake",
            "Coffee Intake",
            "Pork Intake",
            "Time Spent Driving",
            "Non Tobacco Smoker",
            "Current Tobacco Smoker (Most Days)",
            "Current Smoker",
            "Physical Activity: weeding, lawn mower, carpentry",
            "Motor Vehicle Transport",
            "Water Intake",
            "Cheese Intake")
)

GxEmatch <- merge(GxEids, IDmatching, by.x = "Eid_train", by.y = "ID")


### Find the scoring for specific IDs.
IDscores <- data.frame(
  ID = c("alcohol_drinker_status_f20117_0_0_Current",
         "alcohol_drinker_status_f20117_0_0_Never",
         "alcohol_intake_frequency_f1558_0_0",
         "bread_intake_f1438_0_0",
         "cheese_intake_f1408_0_0",
         "coffee_intake_f1498_0_0",
         "current_tobacco_smoking_f1239_0_0_No",
         "current_tobacco_smoking_f1239_0_0_Yes._on_most_or_all_days",
         "pork_intake_f1389_0_0",
         "plays_computer_games_f2237_0_0",
         "smoking_status_f20116_0_0_Current",
         "time_spent_driving_f1090_0_0",
         "types_of_physical_activity_in_last_4_weeks_f6164_0_0.multi_Heavy_DIY_.eg._weeding._lawn_mowing._carpentry._digging.",
         "types_of_transport_used_excluding_work_f6162_0_0.multi_Car.motor_vehicle",
         "types_of_transport_used_excluding_work_f6162_0_0.multi_Public_transport",
         "water_intake_f1528_0_0"),
  Scores = I(list(c("NO","YES"),
             c("NO","YES"),
             c("Never","Special Occasions",
               "1-3 times monthly","1-2 times weekly",
               "3-4 times weekly","Daily"),
             c(),
             c(),
             c(),
             c("NO","YES"),
             c("NO","YES"),
             c(),
             c("Never","Sometimes","Often"),
             c("NO","YES"),
             c(),
             c("NO","YES"),
             c("NO","YES"),
             c("NO","YES"),
             c()))
)

GxEmatchv2 <- merge(GxEmatch, IDscores, by.x = "Eid_train", by.y="ID")

lapply(1:nrow(GxEmatch), function(x){
  oID <- GxEmatchv2[x, ][["omicID"]]
  eID <- GxEmatchv2[x, ][["Eid_train"]]
  nID <- GxEmatchv2[x, ][["Name"]]
  lID <- GxEmatchv2[x,][["Scores"]][[1]]
  
  saveGWISplot_app(o = oID,
                   e = eID,
                   c = "Type6",
                   n = nID,
                   l = lID,
                   p = "Polygenic GxE Interaction",
                   f = paste0(nID,"_GxE"))
}) 



c("Never","Special Occasions",
  "1-3 times monthly","1-2 times weekly",
  "3-4 times weekly","Daily")


saveGWISplot_app(o = "APOF",
                 e = "alcohol_intake_frequency_f1558_0_0",
                 c = "Type6",
                 n = "Alcohol Intake Freq.",
                 l = c("Never","Special Occasions",
                       "1-3 times monthly","1-2 times weekly",
                       "3-4 times weekly","Daily"),
                 p = "Polygenic GxE Interaction",
                 f = "Alcohol_GxE")








