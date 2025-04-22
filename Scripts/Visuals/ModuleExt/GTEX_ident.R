library(tidyverse)
library(data.table)

#'*MAKE 2 TERM2GENE files one with GeneSymbol other with ENSEMBL ID*

#TODO:
#Goal is to aggregate all proteins in UK Biobank
#Find respective Gene Name // Transcript Info
#Subset everything
#Split by tissue
#Categorize genes based on where they are expressed
#Heuristic Based^^^

## Downlaod GTEx ------
# Get a list of all files in the folder
file_list <- list.files(path = "/n/groups/patel/shakson_ukb/UK_Biobank/Data/GTEX/RNA/", full.names = TRUE)

# Read all files (assuming they are CSV files in this case)
data_list <- lapply(file_list, fread)

names(data_list) <- list.files(path = "/n/groups/patel/shakson_ukb/UK_Biobank/Data/GTEX/RNA/", full.names = F)
names(data_list) <- gsub("gene_tpm_v10_", "", names(data_list))
names(data_list) <- gsub(".gct", "", names(data_list))

#General Transcriptome:
median_expr <- function(df, tissueID){
  #Get median expression:
  df2 <- df %>%
    pivot_longer(-c(Name, Description)) %>%
    group_by(Name) %>%
    summarize(median_value = median(value)) %>%
    left_join(df, by = "Name") %>%
    select(Name, Description, median_value)
  
  #Rename:
  colnames(df2) <- c("Name","Description",tissueID)
  
  return(df2)
}
data_list <- lapply(names(data_list), function(id){
  median_expr(data_list[[id]],id)
})

#Create main GTEx DataFrame
GTEX_df <- data_list %>% reduce(full_join)
GTEX_spec <- GTEX_df %>%
  pivot_longer(-c(Name, Description)) %>%
  group_by(Name) %>%
  mutate(value = log1p(value)) %>% #to make sure all expression values are positive.
  summarize(tau_score = sum(1 - (value / max(value))) / (length(value) - 1),
            max_expr = max(value)) %>%
  left_join(GTEX_df, by = "Name") %>%
  select(Name, Description, tau_score, max_expr)

#Get plot of tau score overlap between all transcriptome and UKB proteome.
#Tau score files to compare: is on HPA atlas website 
library(ggplot2)

omicpredIDs <- fread(file = "/n/groups/patel/shakson_ukb/UK_Biobank/Data/OMICSPRED/UKB_Olink_multi_ancestry_models_val_results_portal.csv")
GTEX_spec %>%
  ggplot(aes(x = tau_score)) +
  geom_density(aes(fill = "GTEx"), size = 1, alpha = 0.4) +
  geom_density(data = subset(GTEX_spec, 
                             Description %in% unlist(strsplit(omicpredIDs$Gene, 
                                                              "_"))), 
               aes(fill = "UKB Proteomics"), 
               size = 1, alpha = 0.4) +
  geom_vline(xintercept = 0.8, linetype="dashed", 
             color = "darkgoldenrod", size=1) +
  labs(title = "Quantifying Tissue Specificity of Plasma Proteomics",
       x = "Tissue Specificity (Tau score)",
       y = "Density",
       fill = "Data Type") +
  theme_minimal()



#'*Calculate log fold change:*

# Get expressiond ata
expression_data <- GTEX_df %>%
  select(-Name, -Description)  # Select only the tissue columns

# Calculate fold change: x_i / mean(x_-i)
# Optimized fold change calculation
fold_change_results <- apply(expression_data, 1, function(x) {
  # Calculate the total sum of all tissue expressions
  total_sum <- sum(x)
  
  # Calculate fold change for each tissue
  fold_changes <- x / ((total_sum - x) / (length(x) - 1))
  #B/c
  # x(i) / sum(x) - x(i)/(n - 1) ==> x(i) * (n-1) / (sum(x) - x)
  
  # Return the fold changes
  return(fold_changes)
})


# Convert the list of fold changes into a data frame for easier manipulation
fold_change_df <- as.data.frame(t(fold_change_results))
fold_change_df$Name <- GTEX_df$Name
fold_change_df$Description <- GTEX_df$Description


GTEx_GeneSet <- lapply(colnames(fold_change_df)[1:54], function(x){
  # List genes with >4 FC for each tissues
  fold_change_df %>%
    filter(.[[x]] > 4) %>% 
    pull(Description)
})
names(GTEx_GeneSet) <- colnames(fold_change_df)[1:54]


# Initialize an empty data frame
Term2Gene <- data.frame(Term = character(0), Gene = character(0))

# Loop through the list and convert it into the required format
for (term in names(GTEx_GeneSet)) {
  genes <- GTEx_GeneSet[[term]]  # Extract the vector of genes
  term_data <- data.frame(Term = rep(term, length(genes)), Gene = genes)
  Term2Gene <- rbind(Term2Gene, term_data)
}

# Now Term2Gene is in the correct format for ClusterProfiler
head(Term2Gene)

write.table(Term2Gene, 
            file = "/n/groups/patel/shakson_ukb/UK_Biobank/Data/GeneSets/GTEX_tissue.txt",
            sep = "\t", row.names = F)

