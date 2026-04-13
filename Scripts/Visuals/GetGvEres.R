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

fwrite(final_tbl, file.path(indir, "R2_GE_summary_all_proteins.txt"), sep = "\t")

# quick sanity summaries
cat("n proteins:", nrow(final_tbl), "\n")
cat("n E>G:", sum(final_tbl$GE_dominant == "E > G", na.rm = TRUE), "\n")
summary(final_tbl$R2_E)
summary(final_tbl$R2_G)