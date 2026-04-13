library(readxl)
library(data.table)
library(dplyr)
library(stringr)
library(tidyr)
library(purrr)
library(ggplot2)
library(scales)
library(broom)

options(stringsAsFactors = FALSE)

# =========================================================
# 1. Load Sun et al. heritability / pQTL summary table
# =========================================================

UKBprotH2 <- read_excel(
  "/n/groups/patel/IGLOO/UKB/pQTLmetadata/BenSunNature2023.xlsx",
  sheet = "ST19",
  skip = 3
)

colnames(UKBprotH2) <- c(
  "UKBprotID", "cispQTLs", "transpQTLs",
  "allpQTLs", "polygen_comp", "total_herit",
  "THcis", "THtrans"
)

# Extract the protein ID before the colon
UKBprotH2 <- UKBprotH2 %>%
  mutate(
    protID = str_extract(UKBprotID, "^[^:]+")
  )

# Coerce relevant columns to numeric just in case
num_cols_h2 <- c("cispQTLs", "transpQTLs", "allpQTLs",
                 "polygen_comp", "total_herit", "THcis", "THtrans")

UKBprotH2 <- UKBprotH2 %>%
  mutate(across(all_of(num_cols_h2), as.numeric))

# =========================================================
# 2. Load HEAP R2 and lasso outputs
# =========================================================

root <- "/n/groups/patel/shakson_ukb/UK_Biobank/Data/Parallel/Module1"
covarType <- "Type5"
indir <- file.path(root, covarType)
outdir = "/n/groups/patel/shakson_ukb/UK_Biobank/RScripts/Pure_StatGen/Prot_ExPGS/Visualizations/ModuleHerit/"

r2_files <- list.files(indir, pattern = "^R2groups_[0-9]+\\.txt$", full.names = TRUE)
lasso_files <- list.files(indir, pattern = "^lassofit_[0-9]+\\.txt$", full.names = TRUE)

if (length(r2_files) == 0) stop("No R2groups_*.txt files found.")
if (length(lasso_files) == 0) stop("No lassofit_*.txt files found.")

R2 <- r2_files %>%
  map_dfr(~ fread(.x)) %>%
  mutate(omic = as.character(omic))

LASSO <- lasso_files %>%
  map_dfr(~ fread(.x)) %>%
  mutate(omic = as.character(omic))

# =========================================================
# 3. Helper function
# =========================================================

add_missing_cols <- function(df, cols, fill = 0) {
  missing <- setdiff(cols, names(df))
  if (length(missing) > 0) {
    for (m in missing) df[[m]] <- fill
  }
  df
}

# =========================================================
# 4. Collapse HEAP Shapley groups into G / E / GxE / Covars
# =========================================================

R2_components_by_fold <- R2 %>%
  mutate(
    component = case_when(
      group %in% c("Gcis", "Gtrans") ~ "G",
      group == "Covars" ~ "Covars",
      str_starts(group, "E_") ~ "E",
      str_starts(group, "GxEcis_") | str_starts(group, "GxEtrans_") ~ "GxE",
      TRUE ~ "Other"
    )
  ) %>%
  group_by(omic, fold, component) %>%
  summarise(r2 = sum(r2_main, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = component, values_from = r2, values_fill = 0) %>%
  { add_missing_cols(., c("G", "E", "GxE", "Covars"), fill = 0) }

R2_components <- R2_components_by_fold %>%
  group_by(omic) %>%
  summarise(
    R2_G = mean(G, na.rm = TRUE),
    R2_E = mean(E, na.rm = TRUE),
    R2_GxE = mean(GxE, na.rm = TRUE),
    R2_Covars = mean(Covars, na.rm = TRUE),
    .groups = "drop"
  )

# =========================================================
# 5. Add predictive lasso R2
# =========================================================

LASSO_by_prot <- LASSO %>%
  group_by(omic) %>%
  summarise(
    R2_test_lasso = mean(test_lasso, na.rm = TRUE),
    R2_train_lasso = mean(train_lasso, na.rm = TRUE),
    .groups = "drop"
  )

final_tbl <- R2_components %>%
  left_join(LASSO_by_prot, by = "omic") %>%
  arrange(desc(R2_G))

# =========================================================
# 6. Compute cis and trans genetic components correctly
# =========================================================

R2_cis_trans <- R2 %>%
  filter(group %in% c("Gcis", "Gtrans")) %>%
  group_by(omic, fold, group) %>%
  summarise(r2 = sum(r2_main, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = group, values_from = r2, values_fill = 0) %>%
  { add_missing_cols(., c("Gcis", "Gtrans"), fill = 0) } %>%
  group_by(omic) %>%
  summarise(
    R2_Gcis = mean(Gcis, na.rm = TRUE),
    R2_Gtrans = mean(Gtrans, na.rm = TRUE),
    .groups = "drop"
  )

# =========================================================
# 7. Merge Sun et al. with HEAP summaries
# =========================================================

comparetbl <- UKBprotH2 %>%
  inner_join(final_tbl, by = c("protID" = "omic"))

comparetblv2 <- UKBprotH2 %>%
  inner_join(R2_cis_trans, by = c("protID" = "omic"))

cat("Number of proteins in total comparison:", nrow(comparetbl), "\n")
cat("Number of proteins in cis/trans comparison:", nrow(comparetblv2), "\n")

# =========================================================
# 8. Correlations
# =========================================================

safe_cor_test <- function(df, x, y) {
  sub <- df %>% select(all_of(c(x, y))) %>% filter(is.finite(.data[[x]]), is.finite(.data[[y]]))
  if (nrow(sub) < 3) return(NULL)
  cor.test(sub[[x]], sub[[y]], method = "pearson")
}

cor_R2G_allpQTLs   <- safe_cor_test(comparetbl, "R2_G", "allpQTLs")
cor_R2G_totalherit <- safe_cor_test(comparetbl, "R2_G", "total_herit")

cor_Gcis_cispQTLs <- safe_cor_test(comparetblv2, "R2_Gcis", "cispQTLs")
cor_Gcis_THcis    <- safe_cor_test(comparetblv2, "R2_Gcis", "THcis")

cor_Gtrans_transpQTLs <- safe_cor_test(comparetblv2, "R2_Gtrans", "transpQTLs")
cor_Gtrans_THtrans    <- safe_cor_test(comparetblv2, "R2_Gtrans", "THtrans")

print(cor_R2G_allpQTLs)
print(cor_R2G_totalherit)
print(cor_Gcis_cispQTLs)
print(cor_Gcis_THcis)
print(cor_Gtrans_transpQTLs)
print(cor_Gtrans_THtrans)

# =========================================================
# 9. Linear regressions
# =========================================================

safe_lm <- function(df, x, y) {
  sub <- df %>% select(all_of(c(x, y))) %>% filter(is.finite(.data[[x]]), is.finite(.data[[y]]))
  if (nrow(sub) < 3) return(NULL)
  lm(reformulate(x, response = y), data = sub)
}

lm_R2G_allpQTLs   <- safe_lm(comparetbl, "R2_G", "allpQTLs")
lm_R2G_totalherit <- safe_lm(comparetbl, "R2_G", "total_herit")

lm_Gcis_cispQTLs <- safe_lm(comparetblv2, "R2_Gcis", "cispQTLs")
lm_Gcis_THcis    <- safe_lm(comparetblv2, "R2_Gcis", "THcis")

lm_Gtrans_transpQTLs <- safe_lm(comparetblv2, "R2_Gtrans", "transpQTLs")
lm_Gtrans_THtrans    <- safe_lm(comparetblv2, "R2_Gtrans", "THtrans")

summary(lm_R2G_allpQTLs)
summary(lm_R2G_totalherit)
summary(lm_Gcis_cispQTLs)
summary(lm_Gcis_THcis)
summary(lm_Gtrans_transpQTLs)
summary(lm_Gtrans_THtrans)

# Tidy regression summaries
reg_results <- bind_rows(
  tidy(lm_R2G_allpQTLs) %>% mutate(model = "allpQTLs ~ R2_G"),
  tidy(lm_R2G_totalherit) %>% mutate(model = "total_herit ~ R2_G"),
  tidy(lm_Gcis_cispQTLs) %>% mutate(model = "cispQTLs ~ R2_Gcis"),
  tidy(lm_Gcis_THcis) %>% mutate(model = "THcis ~ R2_Gcis"),
  tidy(lm_Gtrans_transpQTLs) %>% mutate(model = "transpQTLs ~ R2_Gtrans"),
  tidy(lm_Gtrans_THtrans) %>% mutate(model = "THtrans ~ R2_Gtrans")
)

write.table(
  reg_results,
  file = file.path(outdir, paste0("heritability_regression_results_", covarType, ".tsv")),
  sep = "\t", row.names = FALSE, quote = FALSE
)

# =========================================================
# 10. Summary statistics to show attenuation of HEAP R2_G
# =========================================================

summary_stats_total <- comparetbl %>%
  summarise(
    n = n(),
    mean_R2_G = mean(R2_G, na.rm = TRUE),
    median_R2_G = median(R2_G, na.rm = TRUE),
    mean_allpQTLs = mean(allpQTLs, na.rm = TRUE),
    median_allpQTLs = median(allpQTLs, na.rm = TRUE),
    mean_total_herit = mean(total_herit, na.rm = TRUE),
    median_total_herit = median(total_herit, na.rm = TRUE)
  )

summary_stats_cis_trans <- comparetblv2 %>%
  summarise(
    n = n(),
    mean_R2_Gcis = mean(R2_Gcis, na.rm = TRUE),
    median_R2_Gcis = median(R2_Gcis, na.rm = TRUE),
    mean_cispQTLs = mean(cispQTLs, na.rm = TRUE),
    median_cispQTLs = median(cispQTLs, na.rm = TRUE),
    mean_THcis = mean(THcis, na.rm = TRUE),
    median_THcis = median(THcis, na.rm = TRUE),
    mean_R2_Gtrans = mean(R2_Gtrans, na.rm = TRUE),
    median_R2_Gtrans = median(R2_Gtrans, na.rm = TRUE),
    mean_transpQTLs = mean(transpQTLs, na.rm = TRUE),
    median_transpQTLs = median(transpQTLs, na.rm = TRUE),
    mean_THtrans = mean(THtrans, na.rm = TRUE),
    median_THtrans = median(THtrans, na.rm = TRUE)
  )

print(summary_stats_total)
print(summary_stats_cis_trans)

write.table(
  summary_stats_total,
  file = file.path(outdir, paste0("heritability_summary_total_", covarType, ".tsv")),
  sep = "\t", row.names = FALSE, quote = FALSE
)

write.table(
  summary_stats_cis_trans,
  file = file.path(outdir, paste0("heritability_summary_cistrans_", covarType, ".tsv")),
  sep = "\t", row.names = FALSE, quote = FALSE
)

# =========================================================
# 11. Plotting helper
# =========================================================

make_scatter <- function(df, x, y, xlab, ylab, title_txt, outfile) {
  sub <- df %>%
    select(all_of(c(x, y))) %>%
    filter(is.finite(.data[[x]]), is.finite(.data[[y]]))
  
  ct <- cor.test(sub[[x]], sub[[y]], method = "pearson")
  fit <- lm(reformulate(x, response = y), data = sub)
  
  xmax <- max(sub[[x]], na.rm = TRUE)
  ymax <- max(sub[[y]], na.rm = TRUE)
  limmax <- max(xmax, ymax)
  
  p <- ggplot(sub, aes(x = .data[[x]], y = .data[[y]])) +
    geom_point(alpha = 0.5, size = 1.5) +
    geom_smooth(method = "lm", se = TRUE) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    annotate(
      "text",
      x = limmax * 0.05,
      y = limmax * 0.95,
      hjust = 0,
      vjust = 1,
      label = paste0(
        "Pearson r = ", round(unname(ct$estimate), 3),
        "\nP = ", signif(ct$p.value, 3),
        "\nSlope = ", round(coef(fit)[2], 3),
        "\nIntercept = ", round(coef(fit)[1], 3),
        "\nN = ", nrow(sub)
      ),
      size = 4
    ) +
    labs(
      x = xlab,
      y = ylab,
      title = title_txt
    ) +
    theme_bw(base_size = 13)
  
  ggsave(outfile, p, width = 6.5, height = 5.5, dpi = 300)
  return(p)
}

# =========================================================
# 12. Make plots
# =========================================================

p1 <- make_scatter(
  comparetbl,
  x = "R2_G",
  y = "total_herit",
  xlab = "HEAP genetic R2 (PGS-captured)",
  ylab = "Sun et al. total heritability",
  title_txt = "HEAP genetic R2 vs total protein heritability",
  outfile = file.path(outdir, paste0("scatter_R2G_vs_totalherit_", covarType, ".png"))
)

p2 <- make_scatter(
  comparetbl,
  x = "R2_G",
  y = "allpQTLs",
  xlab = "HEAP genetic R2 (PGS-captured)",
  ylab = "Sun et al. variance explained by all pQTLs",
  title_txt = "HEAP genetic R2 vs all pQTL variance explained",
  outfile = file.path(outdir, paste0("scatter_R2G_vs_allpQTLs_", covarType, ".png"))
)

p3 <- make_scatter(
  comparetblv2,
  x = "R2_Gcis",
  y = "THcis",
  xlab = "HEAP cis genetic R2",
  ylab = "Sun et al. cis heritability",
  title_txt = "HEAP cis R2 vs cis heritability",
  outfile = file.path(outdir, paste0("scatter_R2Gcis_vs_THcis_", covarType, ".png"))
)

p4 <- make_scatter(
  comparetblv2,
  x = "R2_Gcis",
  y = "cispQTLs",
  xlab = "HEAP cis genetic R2",
  ylab = "Sun et al. variance explained by cis pQTLs",
  title_txt = "HEAP cis R2 vs cis pQTL variance explained",
  outfile = file.path(outdir, paste0("scatter_R2Gcis_vs_cispQTLs_", covarType, ".png"))
)

p5 <- make_scatter(
  comparetblv2,
  x = "R2_Gtrans",
  y = "THtrans",
  xlab = "HEAP trans genetic R2",
  ylab = "Sun et al. trans heritability",
  title_txt = "HEAP trans R2 vs trans heritability",
  outfile = file.path(outdir, paste0("scatter_R2Gtrans_vs_THtrans_", covarType, ".png"))
)

p6 <- make_scatter(
  comparetblv2,
  x = "R2_Gtrans",
  y = "transpQTLs",
  xlab = "HEAP trans genetic R2",
  ylab = "Sun et al. variance explained by trans pQTLs",
  title_txt = "HEAP trans R2 vs trans pQTL variance explained",
  outfile = file.path(outdir, paste0("scatter_R2Gtrans_vs_transpQTLs_", covarType, ".png"))
)

print(p1)
print(p2)
print(p3)
print(p4)
print(p5)
print(p6)

# =========================================================
# 13. Save merged tables for inspection
# =========================================================

write.table(
  comparetbl,
  file = file.path(outdir, paste0("HEAP_vs_Sun_total_comparison_", covarType, ".tsv")),
  sep = "\t", row.names = FALSE, quote = FALSE
)

write.table(
  comparetblv2,
  file = file.path(outdir, paste0("HEAP_vs_Sun_cistrans_comparison_", covarType, ".tsv")),
  sep = "\t", row.names = FALSE, quote = FALSE
)

# =========================================================
# 14. Optional ratio summaries
# =========================================================
# These are useful if you want a direct sentence saying that
# HEAP PGS-captured R2 is smaller than published heritability.

ratio_tbl <- comparetbl %>%
  mutate(
    ratio_R2G_to_totalherit = R2_G / total_herit,
    ratio_R2G_to_allpQTLs = R2_G / allpQTLs
  ) %>%
  summarise(
    median_ratio_R2G_to_totalherit = median(ratio_R2G_to_totalherit[is.finite(ratio_R2G_to_totalherit)], na.rm = TRUE),
    mean_ratio_R2G_to_totalherit = mean(ratio_R2G_to_totalherit[is.finite(ratio_R2G_to_totalherit)], na.rm = TRUE),
    median_ratio_R2G_to_allpQTLs = median(ratio_R2G_to_allpQTLs[is.finite(ratio_R2G_to_allpQTLs)], na.rm = TRUE),
    mean_ratio_R2G_to_allpQTLs = mean(ratio_R2G_to_allpQTLs[is.finite(ratio_R2G_to_allpQTLs)], na.rm = TRUE)
  )

print(ratio_tbl)

write.table(
  ratio_tbl,
  file = file.path(outdir, paste0("heritability_ratio_summary_", covarType, ".tsv")),
  sep = "\t", row.names = FALSE, quote = FALSE
)

cat("Done.\n")