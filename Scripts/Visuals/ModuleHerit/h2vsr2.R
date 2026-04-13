library(readxl)

## TODO:
## LOAD in H2 stats from Ben Sun Paper
## PREPARE the R2 genetics column from the HEAP runs
## Compare the R2 genetics to H2 in scatterplot with the regression coefficient


#Load in excel file and go to page 19.
UKBprotH2 <- read_excel(
  "/n/groups/patel/IGLOO/UKB/pQTLmetadata/BenSunNature2023.xlsx",
  sheet = "ST19",
  skip = 3            # skips first 4 rows, so reading starts at row 5
)
colnames(UKBprotH2) <- c("UKBprotID", "cispQTLs", "transpQTLs",
                         "allpQTLs", "polygen_comp", "total_herit",   
                         "THcis", "THtrans" )


library(data.table)
library(dplyr)
library(stringr)
library(tidyr)
library(purrr)

root <- "/n/groups/patel/shakson_ukb/UK_Biobank/Data/Parallel/Module1"
covarType <- "Type5"
indir <- file.path(root, covarType)

# -------------------------
# Load all chunk outputs
# -------------------------
r2_files <- list.files(indir, pattern = "^R2groups_\\d+\\.txt$", full.names = TRUE)
lasso_files <- list.files(indir, pattern = "^lassofit_\\d+\\.txt$", full.names = TRUE)

R2 <- r2_files %>%
  map_dfr(~ fread(.x)) %>%
  mutate(omic = as.character(omic))

LASSO <- lasso_files %>%
  map_dfr(~ fread(.x)) %>%
  mutate(omic = as.character(omic))

# -------------------------
# Collapse Shapley groups -> components
# -------------------------
library(dplyr)
library(tidyr)
library(stringr)

# helper: add missing columns (if absent) with a default value
add_missing_cols <- function(df, cols, fill = 0) {
  missing <- setdiff(cols, names(df))
  if (length(missing) > 0) {
    for (m in missing) df[[m]] <- fill
  }
  df
}

R2_components_by_fold <- R2 %>%
  mutate(
    component = case_when(
      group %in% c("Gcis","Gtrans") ~ "G",
      group == "Covars" ~ "Covars",
      str_starts(group, "E_") ~ "E",
      str_starts(group, "GxEcis_") | str_starts(group, "GxEtrans_") ~ "GxE",
      TRUE ~ "Other"
    )
  ) %>%
  group_by(omic, fold, component) %>%
  summarise(r2 = sum(r2_main, na.rm = TRUE), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = component, values_from = r2, values_fill = 0) %>%
  { add_missing_cols(., c("G","E","GxE","Covars"), fill = 0) }

R2_components <- R2_components_by_fold %>%
  group_by(omic) %>%
  summarise(
    R2_G      = mean(G, na.rm = TRUE),
    R2_E      = mean(E, na.rm = TRUE),
    R2_GxE    = mean(GxE, na.rm = TRUE),
    R2_Covars = mean(Covars, na.rm = TRUE),
    .groups = "drop"
  )
# -------------------------
# Add total predictive R2 (from test_lasso)
# -------------------------
LASSO_by_prot <- LASSO %>%
  group_by(omic) %>%
  summarise(
    R2_test_lasso = mean(test_lasso, na.rm = TRUE),
    R2_train_lasso = mean(train_lasso, na.rm = TRUE),
    .groups = "drop"
  )

final_tbl <- R2_components %>%
  left_join(LASSO_by_prot, by = "omic") %>%
  arrange(desc(R2_E))

### Get the regression/correlation:

UKBprotH2$UKBprotID

library(stringr)
UKBprotH2$protID <- str_extract(UKBprotH2$UKBprotID, "^[^:]+")

comparetbl <- merge(UKBprotH2,final_tbl, by.x = "protID", by.y = "omic")

cor.test(comparetbl$R2_G, comparetbl$allpQTLs)
cor.test(comparetbl$R2_G, comparetbl$total_herit)

# Get the cis- and trans- effects:

R2_cis_trans <- R2 %>%
  filter(group %in% c("Gcis","Gtrans")) %>%
  group_by(omic, group) %>%
  tidyr::pivot_wider(names_from = group, values_from = r2_main, values_fill = 0) %>%
  summarise(
    R2_Gcis      = mean(Gcis, na.rm = TRUE),
    R2_Gtrans      = mean(Gtrans, na.rm = TRUE),
    .groups = "drop"
  )

comparetblv2 <- merge(UKBprotH2,R2_cis_trans, by.x = "protID", by.y = "omic")
cor.test(comparetblv2$R2_Gcis, comparetblv2$cispQTLs)
cor.test(comparetblv2$R2_Gcis, comparetblv2$THcis)

cor.test(comparetblv2$R2_Gtrans, comparetblv2$transpQTLs)
cor.test(comparetblv2$R2_Gtrans, comparetblv2$THtrans)



