library(data.table)
library(tidyverse)
library(ggplot2)
library(ggpmisc)
library(ggrepel)
library(gridExtra)

#Addn Libraries for Visualization:
library(patchwork)

setwd("/n/groups/patel/shakson_ukb/UK_Biobank/")

#Load All Results:
#'*Obtain Results from nested CV*
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


# purrr::map2 iterate over the list created to generate one dataframe:
R2test <- map2_dfr(HEAPres, names(HEAPres), ~ .x$R2test %>%
                     select(all_of(c("ID", "G", "E", "GxE"))) %>%
                     mutate(cID = .y))

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


median(plot_testR2_CI$G_stats_mean)
median(plot_testR2_CI$E_stats_mean)

summary(plot_testR2_CI$G_stats_mean)
summary(plot_testR2_CI$E_stats_mean)
percenti

quantile(plot_testR2_CI$E_stats_mean, probs = seq(.1, .95, by = .05))
quantile(plot_testR2_CI$G_stats_mean, probs = seq(.1, .95, by = .05))
#ntile(plot_testR2_CI$G_stats_mean, 10)


#Define genesets to test:
Eprot <- plot_testR2_CI %>%
            filter(E_stats_mean > 0.05) %>%
            pull(ID)

Gprot <- plot_testR2_CI %>%
            filter(G_stats_mean > 0.40) %>% #0.25 orig --> 0.3 led to less significant hits.
            pull(ID)

#Convert to EntrezID:
df_entrez <- fread(file = "/n/groups/patel/shakson_ukb/UK_Biobank/Data/GeneSets/OlinkEntrezConv.txt")

Eprot <- df_entrez[match(Eprot, df_entrez$Gene),]$entrezgene_id
Gprot <- df_entrez[match(Gprot, df_entrez$Gene),]$entrezgene_id


# Function: ORA of Pathways (Uses Entrez ID)
library(clusterProfiler)
library(ReactomePA)
library(DOSE)
ORA_PGSPXS_paths <- function(GeneList) {
  EnrichORA <- list()
  
  # Load KEGG Pathways:
  KEGG_T2G <- fread("/n/groups/patel/shakson_ukb/UK_Biobank/Data/GeneSets/KEGG_T2G.txt")
  KEGG_T2N <- fread("/n/groups/patel/shakson_ukb/UK_Biobank/Data/GeneSets/KEGG_T2N.txt")
  KEGG_T2N <- KEGG_T2N[,c(2,1)]
  
  # Setup Universe w/ Entrez IDs
  omicpredIDs <- fread(file = "/n/groups/patel/shakson_ukb/UK_Biobank/Data/OMICSPRED/UKB_Olink_multi_ancestry_models_val_results_portal.csv")
  df_entrez <- fread(file = "/n/groups/patel/shakson_ukb/UK_Biobank/Data/GeneSets/OlinkEntrezConv.txt")
  
  UKB_universe = unlist(strsplit(omicpredIDs$Gene, "_"))
  UKB_universe_entrez <- df_entrez[match(UKB_universe,
                                         df_entrez$Gene),]$entrezgene_id
  UKB_universe_entrez <- as.character(UKB_universe_entrez)
  
  # KEGG enrichment
  KEGG_ora <- tryCatch({
    enricher(gene = GeneList,
             pvalueCutoff = 0.05,
             pAdjustMethod = "BH",
             universe = UKB_universe_entrez,
             qvalueCutoff = 0.2,
             minGSSize = 10, maxGSSize = 500,
             TERM2GENE = KEGG_T2G, TERM2NAME = KEGG_T2N)
  }, error = function(e) NULL)  # Catch error if KEGG_ora fails
  
  # GO enrichment
  GO_BP <- tryCatch({
    enrichGO(gene = GeneList,
             OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
             ont = "BP", pvalueCutoff = 0.05,
             pAdjustMethod = "BH",
             universe = UKB_universe_entrez,
             qvalueCutoff = 0.2,
             minGSSize = 10, maxGSSize = 500,
             readable = TRUE)
  }, error = function(e) NULL)  # Catch error if GO_BP fails
  
  GO_CC <- tryCatch({
    enrichGO(gene = GeneList,
             OrgDb = org.Hs.eg.db,
             keyType = "ENTREZID",
             ont = "CC", pvalueCutoff = 0.05,
             pAdjustMethod = "BH",
             universe = UKB_universe_entrez,
             qvalueCutoff = 0.2,
             minGSSize = 10, maxGSSize = 500,
             readable = TRUE)
  }, error = function(e) NULL)  # Catch error if GO_CC fails
  
  GO_MF <- tryCatch({
    enrichGO(gene = GeneList,
             OrgDb = org.Hs.eg.db,
             keyType = "ENTREZID",
             ont = "MF", pvalueCutoff = 0.05,
             pAdjustMethod = "BH",
             universe = UKB_universe_entrez,
             qvalueCutoff = 0.2,
             minGSSize = 10, maxGSSize = 500,
             readable = TRUE)
  }, error = function(e) NULL)  # Catch error if GO_MF fails
  
  # Reactome enrichment
  Reactome <- tryCatch({
    enrichPathway(gene = GeneList, 
                  pvalueCutoff = 0.05, readable = TRUE,
                  universe = UKB_universe_entrez,
                  qvalueCutoff = 0.2,
                  minGSSize = 10, maxGSSize = 500)
  }, error = function(e) NULL)  # Catch error if Reactome fails
  
  # DO enrichment
  DO <- tryCatch({
    enrichDO(gene = GeneList,
             ont = "DO", pvalueCutoff = 0.05,
             pAdjustMethod = "BH",
             universe = UKB_universe_entrez,
             qvalueCutoff = 0.2,
             minGSSize = 10, maxGSSize = 500,
             readable = TRUE)
  }, error = function(e) NULL)  # Catch error if DO fails
  
  # Assign results to EnrichORA, check for NULL
  EnrichORA[["KEGG"]] <- if (!is.null(KEGG_ora)) KEGG_ora@result else NULL
  EnrichORA[["GO_BP"]] <- if (!is.null(GO_BP)) GO_BP@result else NULL
  EnrichORA[["GO_CC"]] <- if (!is.null(GO_CC)) GO_CC@result else NULL
  EnrichORA[["GO_MF"]] <- if (!is.null(GO_MF)) GO_MF@result else NULL
  EnrichORA[["Reactome"]] <- if (!is.null(Reactome)) Reactome@result else NULL
  EnrichORA[["DO"]] <- if (!is.null(DO)) DO@result else NULL
  
  return(EnrichORA)
}

Epaths <- ORA_PGSPXS_paths(Eprot)
Gpaths <- ORA_PGSPXS_paths(Gprot)

# Function to Run Over-Representation Analysis
ORA_tissue <- function(GeneList){
  GTEX <- fread("/n/groups/patel/shakson_ukb/UK_Biobank/Data/GeneSets/GTEX_tissue.txt")
  omicpredIDs <- fread(file = "/n/groups/patel/shakson_ukb/UK_Biobank/Data/OMICSPRED/UKB_Olink_multi_ancestry_models_val_results_portal.csv")
  UKB_universe = unlist(strsplit(omicpredIDs$Gene, "_"))
  
  ora <- enricher(gene = GeneList,
                  pvalueCutoff = 0.05,
                  pAdjustMethod = "BH",
                  universe = UKB_universe, #universe is set of proteins surveyed in biobank
                  minGSSize = 10,
                  maxGSSize = 500,
                  qvalueCutoff = 0.2,
                  TERM2GENE = GTEX,
                  TERM2NAME = NA)
  return(ora@result)
}

Eprotid <- plot_testR2_CI %>%
  filter(E_stats_mean > 0.05) %>%
  pull(ID)
Etissue <- ORA_tissue(Eprotid)

Gprotid <- plot_testR2_CI %>%
  filter(G_stats_mean > 0.25) %>%
  pull(ID)
Gtissue <- ORA_tissue(Gprotid)


#'*Function to Run GSEA:*
GSEAdf <- merge(df_entrez, plot_testR2_CI, by.x = "Gene", by.y = "ID")

# #Econtrolled:
# Eprotlist <- GSEAdf$E_stats_mean
# names(Eprotlist) <- GSEAdf$entrezgene_id
# #Gcontrolled:
# Gprotlist <- GSEAdf$G_stats_mean
# names(Gprotlist) <- GSEAdf$entrezgene_id
# #Both should be TRUE:
# sum(is.na(Gprotlist)) == 0
# sum(is.na(Eprotlist)) == 0
# #Sort in Decreasing Order:
# Eprotlist = sort(Eprotlist, decreasing = TRUE)
# Gprotlist = sort(Gprotlist, decreasing = TRUE)

#Rank (handle ties like:)
gene_df <- GSEAdf %>% 
                mutate(rank = rank(E_stats_mean,  ties.method = "random")) %>%
                arrange(desc(rank))
Eprotlist <- gene_df$rank
names(Eprotlist) <- gene_df$entrezgene_id


#Exposome Version:
library(enrichplot)
go <- gseGO(geneList     = Eprotlist,
            OrgDb        = org.Hs.eg.db,
            ont          = "BP",
            minGSSize    = 10,
            maxGSSize    = 500,
            pvalueCutoff = 0.05,
            verbose      = FALSE,
            scoreType = "pos")
head(go)
Ego <- go@result
ego <- gseaplot2(go, geneSetID = c("GO:0034368","GO:0010883","GO:0055088"))


go2 <- gseGO(geneList     = Eprotlist,
            OrgDb        = org.Hs.eg.db,
            ont          = "MF",
            minGSSize    = 10,
            maxGSSize    = 500,
            pvalueCutoff = 0.05,
            verbose      = FALSE,
            scoreType = "pos")
head(go2)
Ego2 <- go2@result
ego2 <- gseaplot2(go2, geneSetID = c("GO:0038024","GO:0034185","GO:0043178"))


#### GENETICS ####
#Rank (handle ties like:)
gene_df <- GSEAdf %>% 
  mutate(rank = rank(G_stats_mean,  ties.method = "random")) %>%
  arrange(desc(rank))
Gprotlist <- gene_df$rank
names(Gprotlist) <- gene_df$entrezgene_id


#Genetics Version:
library(enrichplot)
go1.1 <- gseGO(geneList     = Gprotlist,
            OrgDb        = org.Hs.eg.db,
            ont          = "BP",
            minGSSize    = 10,
            maxGSSize    = 500,
            pvalueCutoff = 0.05,
            verbose      = FALSE,
            scoreType = "pos")
head(go1.1)
Ggo <- go1.1@result
ggo <- gseaplot2(go1.1, geneSetID = c("GO:0032609","GO:0046634"))


go2.1 <- gseGO(geneList     = Gprotlist,
             OrgDb        = org.Hs.eg.db,
             ont          = "MF",
             minGSSize    = 10,
             maxGSSize    = 500,
             pvalueCutoff = 0.05,
             verbose      = FALSE,
             scoreType = "pos")
head(go2.1)
Ggo2 <- go2.1@result
ggo2 <- gseaplot2(go2.1, geneSetID = c("GO:0140375","GO:0032393"))




### Summarize and output:

## PUT pathway analysis results (table) in ./Figures/HEAP/T1/
write.csv(Ego, file = "./Figures/HEAP/T1/Eprot_GOBP.txt")
write.csv(Ego2, file = "./Figures/HEAP/T1/Eprot_GOMF.txt")
write.csv(Ggo, file = "./Figures/HEAP/T1/Gprot_GOBP.txt")
write.csv(Ggo2, file = "./Figures/HEAP/T1/Gprot_GOMF.txt")



## PUT the GSEA plots in ./Figures/HEAP/A1/
lapply(c("png", "svg", "pdf"), 
       function(fmt) ggsave(paste0("./Figures/HEAP/A1/","Eprot_GOBP.", fmt), 
                            plot = ego, width = 6, height = 4, dpi = 500))

lapply(c("png", "svg", "pdf"), 
       function(fmt) ggsave(paste0("./Figures/HEAP/A1/","Eprot_GOMF.", fmt), 
                            plot = ego2, width = 6, height = 4, dpi = 500))

lapply(c("png", "svg", "pdf"), 
       function(fmt) ggsave(paste0("./Figures/HEAP/A1/","Gprot_GOBP.", fmt), 
                            plot = ggo, width = 6, height = 4, dpi = 500))

lapply(c("png", "svg", "pdf"), 
       function(fmt) ggsave(paste0("./Figures/HEAP/A1/","Gprot_GOMF.", fmt), 
                            plot = ggo2, width = 6, height = 4, dpi = 500))






#### Category Pathways: NOT SUPER RELEVANT HERE ####
### If below doesnt work as well -- move on



categories <- c("Alcohol","Diet_Weekly","Smoking","Exercise_MET",
                "Exercise_Freq","Internet_Usage","Deprivation_Indices","Vitamins")

# purrr::map2 iterate over the list created to generate one dataframe:
R2testCat <- map2_dfr(HEAPres, names(HEAPres), ~ .x$R2testCat %>%
                        select(all_of(c("ID", categories))) %>%
                        mutate(cID = .y))


#In this scenario we just want to look at the mean R2 
#associated with each category across specifications
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

plot_testR2Cat_CI <- R2testCat_CI %>%
  select(-c("cID","Alcohol","Diet_Weekly","Smoking","Exercise_MET",
            "Exercise_Freq","Internet_Usage","Deprivation_Indices","Vitamins")) %>%
  unique()


catNames <- c(Alcohol = "mean_Alcohol",
              Diet_Weekly = "mean_Diet_Weekly",
              Smoking = "mean_Smoking",
              Exercise_MET = "mean_Exercise_MET",
              Exercise_Freq = "mean_Exercise_Freq",
              Internet_Usage = "mean_Internet_Usage",
              Deprivation_Indices = "mean_Deprivation_Indices",
              Vitamins = "mean_Vitamins")


testR2_CatEnrich <- plot_testR2Cat_CI %>%
  select(c("ID","mean_Alcohol","mean_Diet_Weekly",
           "mean_Smoking","mean_Exercise_MET",
           "mean_Exercise_Freq","mean_Internet_Usage",
           "mean_Deprivation_Indices","mean_Vitamins")) %>%
  rename(all_of(catNames))


# Function: ORA of Pathways (Uses Entrez ID)
library(clusterProfiler)
library(ReactomePA)
library(DOSE)

# Function to Run Over-Representation Analysis for Pathways and Tissues:
ORA_PGSPXS_paths <- function(GeneList) {
  EnrichORA <- list()
  
  # Load KEGG Pathways:
  KEGG_T2G <- fread("/n/groups/patel/shakson_ukb/UK_Biobank/Data/GeneSets/KEGG_T2G.txt")
  KEGG_T2N <- fread("/n/groups/patel/shakson_ukb/UK_Biobank/Data/GeneSets/KEGG_T2N.txt")
  KEGG_T2N <- KEGG_T2N[,c(2,1)]
  
  # Setup Universe w/ Entrez IDs
  omicpredIDs <- fread(file = "/n/groups/patel/shakson_ukb/UK_Biobank/Data/OMICSPRED/UKB_Olink_multi_ancestry_models_val_results_portal.csv")
  df_entrez <- fread(file = "/n/groups/patel/shakson_ukb/UK_Biobank/Data/GeneSets/OlinkEntrezConv.txt")
  
  UKB_universe = unlist(strsplit(omicpredIDs$Gene, "_"))
  UKB_universe_entrez <- df_entrez[match(UKB_universe,
                                         df_entrez$Gene),]$entrezgene_id
  UKB_universe_entrez <- as.character(UKB_universe_entrez)
  
  # KEGG enrichment
  KEGG_ora <- tryCatch({
    enricher(gene = GeneList,
             pvalueCutoff = 0.05,
             pAdjustMethod = "BH",
             universe = UKB_universe_entrez,
             qvalueCutoff = 0.2,
             minGSSize = 10, maxGSSize = 500,
             TERM2GENE = KEGG_T2G, TERM2NAME = KEGG_T2N)
  }, error = function(e) NULL)  # Catch error if KEGG_ora fails
  
  # GO enrichment
  GO_BP <- tryCatch({
    enrichGO(gene = GeneList,
             OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
             ont = "BP", pvalueCutoff = 0.05,
             pAdjustMethod = "BH",
             universe = UKB_universe_entrez,
             qvalueCutoff = 0.2,
             minGSSize = 10, maxGSSize = 500,
             readable = TRUE)
  }, error = function(e) NULL)  # Catch error if GO_BP fails
  
  GO_CC <- tryCatch({
    enrichGO(gene = GeneList,
             OrgDb = org.Hs.eg.db,
             keyType = "ENTREZID",
             ont = "CC", pvalueCutoff = 0.05,
             pAdjustMethod = "BH",
             universe = UKB_universe_entrez,
             qvalueCutoff = 0.2,
             minGSSize = 10, maxGSSize = 500,
             readable = TRUE)
  }, error = function(e) NULL)  # Catch error if GO_CC fails
  
  GO_MF <- tryCatch({
    enrichGO(gene = GeneList,
             OrgDb = org.Hs.eg.db,
             keyType = "ENTREZID",
             ont = "MF", pvalueCutoff = 0.05,
             pAdjustMethod = "BH",
             universe = UKB_universe_entrez,
             qvalueCutoff = 0.2,
             minGSSize = 10, maxGSSize = 500,
             readable = TRUE)
  }, error = function(e) NULL)  # Catch error if GO_MF fails
  
  # Reactome enrichment
  Reactome <- tryCatch({
    enrichPathway(gene = GeneList, 
                  pvalueCutoff = 0.05, readable = TRUE,
                  universe = UKB_universe_entrez,
                  qvalueCutoff = 0.2,
                  minGSSize = 10, maxGSSize = 500)
  }, error = function(e) NULL)  # Catch error if Reactome fails
  
  # DO enrichment
  DO <- tryCatch({
    enrichDO(gene = GeneList,
             ont = "DO", pvalueCutoff = 0.05,
             pAdjustMethod = "BH",
             universe = UKB_universe_entrez,
             qvalueCutoff = 0.2,
             minGSSize = 10, maxGSSize = 500,
             readable = TRUE)
  }, error = function(e) NULL)  # Catch error if DO fails
  
  # Assign results to EnrichORA, check for NULL
  EnrichORA[["KEGG"]] <- if (!is.null(KEGG_ora)) KEGG_ora@result else NULL
  EnrichORA[["GO_BP"]] <- if (!is.null(GO_BP)) GO_BP@result else NULL
  EnrichORA[["GO_CC"]] <- if (!is.null(GO_CC)) GO_CC@result else NULL
  EnrichORA[["GO_MF"]] <- if (!is.null(GO_MF)) GO_MF@result else NULL
  EnrichORA[["Reactome"]] <- if (!is.null(Reactome)) Reactome@result else NULL
  EnrichORA[["DO"]] <- if (!is.null(DO)) DO@result else NULL
  
  return(EnrichORA)
}
ORA_tissue <- function(GeneList){
  GTEX <- fread("/n/groups/patel/shakson_ukb/UK_Biobank/Data/GeneSets/GTEX_tissue.txt")
  omicpredIDs <- fread(file = "/n/groups/patel/shakson_ukb/UK_Biobank/Data/OMICSPRED/UKB_Olink_multi_ancestry_models_val_results_portal.csv")
  UKB_universe = unlist(strsplit(omicpredIDs$Gene, "_"))
  
  ora <- enricher(gene = GeneList,
                  pvalueCutoff = 0.05,
                  pAdjustMethod = "BH",
                  universe = UKB_universe, #universe is set of proteins surveyed in biobank
                  minGSSize = 10,
                  maxGSSize = 500,
                  qvalueCutoff = 0.2,
                  TERM2GENE = GTEX,
                  TERM2NAME = NA)
  return(ora@result)
}

Eprotid <- testR2_CatEnrich %>%
              filter(Alcohol > 0.02) %>%
              pull(ID)

#Convert to EntrezID:
df_entrez <- fread(file = "/n/groups/patel/shakson_ukb/UK_Biobank/Data/GeneSets/OlinkEntrezConv.txt")
Eprot_ent <- df_entrez[match(Eprotid, df_entrez$Gene),]$entrezgene_id

ExerPaths <- ORA_PGSPXS_paths(Eprot_ent)
ExerTissue <- ORA_tissue(Eprotid)



