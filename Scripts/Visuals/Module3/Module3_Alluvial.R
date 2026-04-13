# ============================================================
# HEAP Mediation — Protein modules based on (PXS/Gcis/Gtrans) NIE patterns
# + Proportion mediated summaries
# + Clean Alluvial: Disease category -> Mediation module -> Dominant driver
#
# What this does:
#  1) Filter significant links (your Bonferroni-any-component rule)
#  2) Compute per-link totals: NIE_total, NDE_total, TE_total (lp/logHR scale)
#  3) Compute "mediation strength":
#       - PM_abs = |NIE| / (|NIE| + |NDE|)   (stable, in [0,1])
#       - PM_signed = NIE / TE              (only if |TE| > eps)
#  4) Build per-protein feature vectors from the significant set:
#       - mean abs NIE components (E/cis/trans)
#       - mean signed NIE components (E/cis/trans)
#       - fractions of dominant driver across links
#       - mean PM_abs and fraction PM_abs>0.5
#       - n_links
#  5) Cluster proteins into K modules (kmeans; fast)
#  6) Alluvial using modules (aggregated; no spaghetti)
#  7) Output a module annotation table (top proteins, driver mix, PM mix)
#
# Inputs:
#  - Type5 data.frame/tibble/data.table with required columns:
#    prot_term, DZ_ID,
#    pxs_NIE_logHR, gcis_NIE_logHR, gtrans_NIE_logHR,
#    pxs_NDE_logHR, gcis_NDE_logHR, gtrans_NDE_logHR,
#    total_TE_logHR (optional but preferred; else TE = NDE_total + NIE_total),
#    pxs_NIE_logHR_delta_p, gcis_NIE_logHR_delta_p, gtrans_NIE_logHR_delta_p
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(ggalluvial)
})

# ----------------------------
# 0) Config
# ----------------------------
CFG <- list(
  alpha_global = 0.05,
  denom_min = 0.01,          # drops tiny |NIE| sums (plot stability)
  te_eps = 0.01,             # TE threshold for PM_signed stability (logHR units)
  K_modules = 12,            # try 10–15 as defaults
  seed = 1,
  alluvial_weight = c("count", "strength")[1],  # "count" or "strength"
  keep_top_dz_cat = 10       # NULL to keep all
)

# ----------------------------
# 1) Filter significant links (your rule)
# ----------------------------
DF <- Type5 %>% as_tibble()
alpha <- CFG$alpha_global / nrow(DF)

DFf <- DF %>%
  filter(
    pxs_NIE_logHR_delta_p    < alpha |
      gcis_NIE_logHR_delta_p   < alpha |
      gtrans_NIE_logHR_delta_p < alpha
  ) %>%
  mutate(
    # disease category from DZ_ID (ICD10 chapter letter embedded in your string)
    icd_letter = str_to_upper(str_match(DZ_ID, "^age_([a-z])")[,2]),
    dz_cat = case_when(
      icd_letter %in% c("A","B") ~ "Infectious",
      icd_letter == "D" ~ "Blood/Immune",
      icd_letter == "E" ~ "Endocrine/Metab",
      icd_letter == "F" ~ "Mental",
      icd_letter == "G" ~ "Nervous",
      icd_letter == "H" ~ "Eye/Ear",
      icd_letter == "I" ~ "Circulatory",
      icd_letter == "J" ~ "Respiratory",
      icd_letter == "K" ~ "Digestive",
      icd_letter == "L" ~ "Skin",
      icd_letter == "M" ~ "MSK",
      icd_letter == "N" ~ "Genitourinary",
      TRUE ~ "Other"
    ),
    
    # per-component NIE & NDE (logHR scale)
    NIE_E    = pxs_NIE_logHR,
    NIE_cis  = gcis_NIE_logHR,
    NIE_tr   = gtrans_NIE_logHR,
    
    NDE_E    = pxs_NDE_logHR,
    NDE_cis  = gcis_NDE_logHR,
    NDE_tr   = gtrans_NDE_logHR,
    
    # totals
    NIE_total = NIE_E + NIE_cis + NIE_tr,
    NDE_total = NDE_E + NDE_cis + NDE_tr,
    
    # TE: prefer provided total_TE_logHR; otherwise reconstruct
    TE_total  = if_else(!is.na(total_TE_logHR), total_TE_logHR, NDE_total + NIE_total),
    
    # stability: drop numerical dust for plotting & clustering
    denom_abs = abs(NIE_E) + abs(NIE_cis) + abs(NIE_tr)
  ) %>%
  filter(denom_abs > CFG$denom_min) %>%
  mutate(
    # dominant driver for THIS link (based on |NIE|)
    driver_link = case_when(
      abs(NIE_E)   >= abs(NIE_cis) & abs(NIE_E)   >= abs(NIE_tr) ~ "E-dominant",
      abs(NIE_cis) >= abs(NIE_E)   & abs(NIE_cis) >= abs(NIE_tr) ~ "Gcis-dominant",
      TRUE ~ "Gtrans-dominant"
    ),
    
    # evidence strength (for optional weighting)
    p_min = pmin(pxs_NIE_logHR_delta_p, gcis_NIE_logHR_delta_p, gtrans_NIE_logHR_delta_p, na.rm = TRUE),
    strength = -log10(pmax(p_min, 1e-300)),
    
    # mediated direction on total NIE
    nie_dir = if_else(NIE_total >= 0, "Mediated risk ↑", "Mediated risk ↓"),
    
    # proportion mediated (stable version)
    PM_abs = abs(NIE_total) / (abs(NIE_total) + abs(NDE_total) + 1e-12),
    
    # signed proportion mediated (use only when TE not tiny)
    PM_signed = if_else(abs(TE_total) >= CFG$te_eps, NIE_total / TE_total, NA_real_),
    
    # flag inconsistent mediation (often yields PM outside [0,1] if you used NIE/TE)
    inconsistent = (sign(NIE_total) != sign(TE_total)) & (abs(TE_total) >= CFG$te_eps)
  )

message("Filtered links (after denom_min): n = ", nrow(DFf),
        " | proteins: ", n_distinct(DFf$prot_term),
        " | diseases: ", n_distinct(DFf$DZ_ID))

# ----------------------------
# 2) Build per-protein features for clustering
#    (computed only from the significant set DFf)
# ----------------------------
prot_feat <- DFf %>%
  group_by(prot_term) %>%
  summarise(
    n_links = n(),
    
    # magnitude summaries per component
    mean_abs_NIE_E   = mean(abs(NIE_E),   na.rm = TRUE),
    mean_abs_NIE_cis = mean(abs(NIE_cis), na.rm = TRUE),
    mean_abs_NIE_tr  = mean(abs(NIE_tr),  na.rm = TRUE),
    
    # signed summaries per component (directional tendency)
    mean_NIE_E   = mean(NIE_E,   na.rm = TRUE),
    mean_NIE_cis = mean(NIE_cis, na.rm = TRUE),
    mean_NIE_tr  = mean(NIE_tr,  na.rm = TRUE),
    
    # driver fractions (how often each dominates)
    frac_E_dom   = mean(driver_link == "E-dominant"),
    frac_cis_dom = mean(driver_link == "Gcis-dominant"),
    frac_tr_dom  = mean(driver_link == "Gtrans-dominant"),
    
    # mediation strength summaries
    mean_PM_abs = mean(PM_abs, na.rm = TRUE),
    frac_PM_abs_gt50 = mean(PM_abs >= 0.5, na.rm = TRUE),
    
    # overall NIE direction tendency
    frac_NIE_pos = mean(NIE_total >= 0, na.rm = TRUE),
    
    # robustness / evidence
    mean_strength = mean(strength, na.rm = TRUE),
    
    .groups = "drop"
  )

# Convert to matrix for clustering
X <- prot_feat %>%
  select(-prot_term) %>%
  as.matrix()

# Standardize features so scales don't dominate (important!)
Xz <- scale(X)

# ----------------------------
# 3) Cluster proteins into modules
# ----------------------------
set.seed(CFG$seed)
km <- kmeans(Xz, centers = CFG$K_modules, nstart = 50, iter.max = 200)

prot2mod <- prot_feat %>%
  mutate(
    module_id = km$cluster,
    module = sprintf("Module %02d", module_id)
  ) %>%
  select(prot_term, module_id, module)

# Attach modules back to links
DFm <- DFf %>%
  inner_join(prot2mod, by = "prot_term")

# ----------------------------
# 4) Build alluvial table (aggregated)
#    Disease category -> Module -> Dominant driver (link-level)
# ----------------------------
DFa0 <- DFm %>%
  mutate(w_all = if (CFG$alluvial_weight == "strength") strength else 1)

# Optionally collapse small disease categories to "Other"
if (!is.null(CFG$keep_top_dz_cat)) {
  keep_cats <- DFa0 %>%
    group_by(dz_cat) %>%
    summarise(tw = sum(w_all), .groups = "drop") %>%
    arrange(desc(tw)) %>%
    slice_head(n = CFG$keep_top_dz_cat) %>%
    pull(dz_cat)
  
  DFa0 <- DFa0 %>%
    mutate(dz_cat = if_else(dz_cat %in% keep_cats, dz_cat, "Other"))
}

DFa <- DFa0 %>%
  count(dz_cat, module, driver_link, wt = w_all, name = "weight") %>%
  rename(driver = driver_link)

# Order factors for aesthetics
dz_levels <- DFa %>%
  group_by(dz_cat) %>%
  summarise(tw = sum(weight), .groups = "drop") %>%
  arrange(desc(tw)) %>%
  pull(dz_cat)

mod_levels <- DFa %>%
  group_by(module) %>%
  summarise(tw = sum(weight), .groups = "drop") %>%
  arrange(desc(tw)) %>%
  pull(module)

DFa <- DFa %>%
  mutate(
    dz_cat = factor(dz_cat, levels = dz_levels),
    module = factor(module, levels = mod_levels),
    driver = factor(driver, levels = c("E-dominant", "Gcis-dominant", "Gtrans-dominant"))
  )

# ----------------------------
# 5) Plot: clean alluvial with modules
# ----------------------------
p_alluvial_modules <- ggplot(DFa,
                             aes(axis1 = dz_cat, axis2 = module, axis3 = driver, y = weight)
) +
  geom_alluvium(aes(fill = driver), alpha = 0.85, width = 1/14) +
  geom_stratum(width = 1/10, fill = "white", color = "black") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
  scale_x_discrete(
    limits = c("Disease category", "Mediation module", "Dominant driver"),
    expand = c(.03, .03)
  ) +
  labs(
    title = "Disease category → mediation modules → dominant driver of mediation",
    fill = "Driver"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold"),
    axis.title = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    legend.position = "right"
  )

print(p_alluvial_modules)

# ----------------------------
# 6) Optional: Focus the figure on "strongly mediated" links only (PM_abs >= 0.5)
#    (This often looks MUCH cleaner and matches your narrative.)
# ----------------------------
DFa_strong <- DFa0 %>%
  filter(PM_abs >= 0.5) %>%
  count(dz_cat, module, driver_link, wt = w_all, name = "weight") %>%
  rename(driver = driver_link) %>%
  mutate(
    dz_cat = factor(dz_cat, levels = dz_levels),
    module = factor(module, levels = mod_levels),
    driver = factor(driver, levels = c("E-dominant", "Gcis-dominant", "Gtrans-dominant"))
  )

p_alluvial_modules_strong <- ggplot(DFa_strong,
                                    aes(axis1 = dz_cat, axis2 = module, axis3 = driver, y = weight)
) +
  geom_alluvium(aes(fill = driver), alpha = 0.85, width = 1/14) +
  geom_stratum(width = 1/10, fill = "white", color = "black") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
  scale_x_discrete(
    limits = c("Disease category", "Mediation module", "Dominant driver"),
    expand = c(.03, .03)
  ) +
  labs(
    title = "Disease category → mediation modules → dominant driver (PM_abs ≥ 0.5)",
    fill = "Driver"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold"),
    axis.title = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    legend.position = "right"
  )

print(p_alluvial_modules_strong)

# ----------------------------
# 7) Module annotation tables (for naming modules later)
# ----------------------------
# (a) Driver mix per module
mod_driver <- DFm %>%
  mutate(w = if (CFG$alluvial_weight == "strength") strength else 1) %>%
  group_by(module, driver_link) %>%
  summarise(w = sum(w), .groups = "drop") %>%
  group_by(module) %>%
  mutate(frac = w / sum(w)) %>%
  arrange(module, desc(frac)) %>%
  summarise(driver_mix = paste0(driver_link, "=", round(100*frac), "%", collapse = "; "), .groups="drop")

# (b) Mediation strength per module
mod_pm <- DFm %>%
  group_by(module) %>%
  summarise(
    n_links = n(),
    mean_PM_abs = mean(PM_abs, na.rm = TRUE),
    frac_PM_abs_gt50 = mean(PM_abs >= 0.5, na.rm = TRUE),
    frac_inconsistent = mean(inconsistent, na.rm = TRUE),
    .groups = "drop"
  )

# (c) Top proteins per module (by #links; you can swap to sum(strength))
mod_top_prots <- DFm %>%
  count(module, prot_term, name = "n") %>%
  group_by(module) %>%
  arrange(desc(n), .by_group = TRUE) %>%
  slice_head(n = 8) %>%
  summarise(top_proteins = paste(prot_term, collapse = ", "), .groups = "drop")

# (d) Top disease categories per module (by #links)
mod_top_dz <- DFm %>%
  count(module, dz_cat, name = "n") %>%
  group_by(module) %>%
  arrange(desc(n), .by_group = TRUE) %>%
  slice_head(n = 3) %>%
  summarise(top_dz = paste0(dz_cat, " (", n, ")", collapse = "; "), .groups = "drop")

module_summary <- tibble(module = sort(unique(DFm$module))) %>%
  left_join(mod_pm, by = "module") %>%
  left_join(mod_driver, by = "module") %>%
  left_join(mod_top_dz, by = "module") %>%
  left_join(mod_top_prots, by = "module")

print(module_summary)

# OPTIONAL: Save for manual renaming of modules
# write.csv(module_summary, "mediation_modules_summary.csv", row.names = FALSE)

# OPTIONAL: Save plots
# ggsave("alluvial_modules.png", p_alluvial_modules, width = 12, height = 6.5, dpi = 300)
# ggsave("alluvial_modules_PMabs_ge_0.5.png", p_alluvial_modules_strong, width = 12, height = 6.5, dpi = 300)
