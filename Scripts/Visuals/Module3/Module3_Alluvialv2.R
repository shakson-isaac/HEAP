# ============================================================
# HEAP Mediation Visualization (FINAL)
# - Build significant link set
# - Compute mediation strength (PM_abs, PM_signed)
# - Cluster proteins into mediation modules using NIE patterns
# - Flip alluvial: Driver -> Module -> Disease category
# - Choose K via silhouette (subsample)
# - Protein-level and disease-level mediation summary plots
#
# INPUT: Type5 (data.frame/tibble/data.table)
# Required columns:
#   prot_term, DZ_ID,
#   pxs_NIE_logHR, gcis_NIE_logHR, gtrans_NIE_logHR,
#   pxs_NDE_logHR, gcis_NDE_logHR, gtrans_NDE_logHR,
#   pxs_NIE_logHR_delta_p, gcis_NIE_logHR_delta_p, gtrans_NIE_logHR_delta_p
# Optional (preferred):
#   total_TE_logHR
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(ggalluvial)
  library(cluster)
})

# ----------------------------
# 0) Config
# ----------------------------
CFG <- list(
  alpha_global = 0.05,
  denom_min = 0.01,          # drop tiny |NIE| sums (stability)
  te_eps = 0.01,             # TE threshold for PM_signed stability (logHR)
  seed = 1,
  
  # K selection
  K_grid = 6:20,
  K_subsample_n = 800,       # silhouette computed on this many proteins (max)
  K_nstart = 25,
  K_itermax = 200,
  
  # Final K (if NULL, pick best by silhouette)
  K_modules = NULL,
  
  # Alluvial options
  alluvial_weight = c("count", "strength")[1], # "count" or "strength"
  keep_top_dz_cat = 10,        # NULL to keep all dz_cat
  driver_levels = c("E-dominant","Gcis-dominant","Gtrans-dominant")
)

# ----------------------------
# 1) Build filtered link set DFf
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
    # ICD10 chapter letter embedded in your DZ_ID string (age_<letter>...)
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
    
    # per-component effects on logHR scale
    NIE_E   = pxs_NIE_logHR,
    NIE_cis = gcis_NIE_logHR,
    NIE_tr  = gtrans_NIE_logHR,
    
    NDE_E   = pxs_NDE_logHR,
    NDE_cis = gcis_NDE_logHR,
    NDE_tr  = gtrans_NDE_logHR,
    
    # totals
    NIE_total = NIE_E + NIE_cis + NIE_tr,
    NDE_total = NDE_E + NDE_cis + NDE_tr,
    
    TE_total = if_else(!is.na(total_TE_logHR), total_TE_logHR, NIE_total + NDE_total),
    
    # dust filter
    denom_abs = abs(NIE_E) + abs(NIE_cis) + abs(NIE_tr)
  ) %>%
  filter(denom_abs > CFG$denom_min) %>%
  mutate(
    # dominant driver for this link (based on |NIE|)
    driver = case_when(
      abs(NIE_E)   >= abs(NIE_cis) & abs(NIE_E)   >= abs(NIE_tr) ~ "E-dominant",
      abs(NIE_cis) >= abs(NIE_E)   & abs(NIE_cis) >= abs(NIE_tr) ~ "Gcis-dominant",
      TRUE ~ "Gtrans-dominant"
    ),
    
    # evidence strength
    p_min = pmin(pxs_NIE_logHR_delta_p, gcis_NIE_logHR_delta_p, gtrans_NIE_logHR_delta_p, na.rm = TRUE),
    strength = -log10(pmax(p_min, 1e-300)),
    
    # direction of mediated effect
    nie_dir = if_else(NIE_total >= 0, "Mediated risk ↑", "Mediated risk ↓"),
    
    # Proportion mediated:
    # stable, bounded measure: how mediated vs direct in magnitude space
    PM_abs = abs(NIE_total) / (abs(NIE_total) + abs(NDE_total) + 1e-12),
    
    # signed PM (classical NIE/TE), only if TE not tiny
    PM_signed = if_else(abs(TE_total) >= CFG$te_eps, NIE_total / TE_total, NA_real_),
    
    # inconsistent mediation flag (NIE opposite sign from TE)
    inconsistent = (abs(TE_total) >= CFG$te_eps) & (sign(NIE_total) != sign(TE_total))
  )

message("Filtered links: n=", nrow(DFf),
        " | proteins=", n_distinct(DFf$prot_term),
        " | diseases=", n_distinct(DFf$DZ_ID),
        " | dz_cat=", n_distinct(DFf$dz_cat))

# ----------------------------
# 2) Protein feature table for clustering (NIE patterns + PM)
# ----------------------------
prot_feat <- DFf %>%
  group_by(prot_term) %>%
  summarise(
    n_links = n(),
    
    # NIE magnitude per component
    mean_abs_NIE_E   = mean(abs(NIE_E),   na.rm = TRUE),
    mean_abs_NIE_cis = mean(abs(NIE_cis), na.rm = TRUE),
    mean_abs_NIE_tr  = mean(abs(NIE_tr),  na.rm = TRUE),
    
    # NIE signed tendency per component
    mean_NIE_E   = mean(NIE_E,   na.rm = TRUE),
    mean_NIE_cis = mean(NIE_cis, na.rm = TRUE),
    mean_NIE_tr  = mean(NIE_tr,  na.rm = TRUE),
    
    # dominant-driver frequency
    frac_E_dom   = mean(driver == "E-dominant"),
    frac_cis_dom = mean(driver == "Gcis-dominant"),
    frac_tr_dom  = mean(driver == "Gtrans-dominant"),
    
    # mediation strength
    mean_PM_abs = mean(PM_abs, na.rm = TRUE),
    median_PM_abs = median(PM_abs, na.rm = TRUE),
    frac_PM_abs_gt50 = mean(PM_abs >= 0.5, na.rm = TRUE),
    
    # mediated direction tendency
    frac_NIE_pos = mean(NIE_total >= 0, na.rm = TRUE),
    
    # evidence
    mean_strength = mean(strength, na.rm = TRUE),
    
    .groups = "drop"
  )

# Feature matrix (standardized)
X <- prot_feat %>% select(-prot_term) %>% as.matrix()
Xz <- scale(X)

# ----------------------------
# 3) Choose K via silhouette (subsample proteins for speed)
# ----------------------------
set.seed(CFG$seed)
nP <- nrow(Xz)
idx <- sample(seq_len(nP), size = min(CFG$K_subsample_n, nP))
Xsub <- Xz[idx, , drop = FALSE]

Ks <- CFG$K_grid
sil <- numeric(length(Ks))

set.seed(CFG$seed)
for (i in seq_along(Ks)) {
  k <- Ks[i]
  km <- kmeans(Xsub, centers = k, nstart = CFG$K_nstart, iter.max = CFG$K_itermax)
  sil[i] <- mean(silhouette(km$cluster, dist(Xsub))[, 3])
}

k_tbl <- tibble(K = Ks, silhouette = sil)
print(k_tbl)

p_k <- ggplot(k_tbl, aes(K, silhouette)) +
  geom_line() + geom_point() +
  theme_bw() +
  labs(title = "Silhouette vs K (kmeans on protein NIE/PM features)",
       y = "Mean silhouette (subsample)")
print(p_k)

bestK <- k_tbl %>% arrange(desc(silhouette)) %>% slice(1) %>% pull(K)
K_final <- if (is.null(CFG$K_modules)) bestK else CFG$K_modules
message("Using K = ", K_final, " (best by silhouette = ", bestK, ")")

# ----------------------------
# 4) Cluster all proteins into K modules (kmeans)
# ----------------------------
set.seed(CFG$seed)
km_full <- kmeans(Xz, centers = K_final, nstart = 50, iter.max = CFG$K_itermax)

prot2mod <- prot_feat %>%
  mutate(
    module_id = km_full$cluster,
    module = sprintf("Module %02d", module_id)
  ) %>%
  select(prot_term, module_id, module)

DFm <- DFf %>% inner_join(prot2mod, by = "prot_term")

# ----------------------------
# 5) Build aggregated alluvial table and FLIP it:
#    Driver -> Module -> Disease category
# ----------------------------
DFa0 <- DFm %>%
  mutate(w_all = if (CFG$alluvial_weight == "strength") strength else 1)

# Optional: keep only top disease categories by flow weight
if (!is.null(CFG$keep_top_dz_cat)) {
  keep_cats <- DFa0 %>%
    group_by(dz_cat) %>%
    summarise(tw = sum(w_all), .groups = "drop") %>%
    arrange(desc(tw)) %>%
    slice_head(n = CFG$keep_top_dz_cat) %>%
    pull(dz_cat)
  
  DFa0 <- DFa0 %>% mutate(dz_cat = if_else(dz_cat %in% keep_cats, dz_cat, "Other"))
}

DFa <- DFa0 %>%
  count(driver, module, dz_cat, wt = w_all, name = "weight")

# Order drivers
DFa <- DFa %>% mutate(driver = factor(driver, levels = CFG$driver_levels))

# Order modules to reduce crossings:
# group modules by their dominant driver, then by strength within that group
mod_order <- DFa %>%
  group_by(module, driver) %>% summarise(w = sum(weight), .groups="drop") %>%
  group_by(module) %>%
  summarise(
    wE = sum(w[driver=="E-dominant"]),
    wC = sum(w[driver=="Gcis-dominant"]),
    wT = sum(w[driver=="Gtrans-dominant"]),
    dom = CFG$driver_levels[which.max(c(wE,wC,wT))],
    dom_w = max(c(wE,wC,wT)),
    .groups="drop"
  ) %>%
  arrange(factor(dom, levels = CFG$driver_levels), desc(dom_w)) %>%
  pull(module)

# Order disease categories
dz_order <- DFa %>%
  group_by(dz_cat) %>% summarise(tw = sum(weight), .groups="drop") %>%
  arrange(desc(tw)) %>%
  pull(dz_cat)

DFa <- DFa %>%
  mutate(
    module = factor(module, levels = mod_order),
    dz_cat = factor(dz_cat, levels = dz_order)
  )

# Plot flipped alluvial
p_alluvial_flip <- ggplot(DFa,
                          aes(axis1 = driver, axis2 = module, axis3 = dz_cat, y = weight)
) +
  geom_alluvium(aes(fill = driver), alpha = 0.85, width = 1/14) +
  geom_stratum(width = 1/10, fill = "white", color = "black") +
  geom_text(stat="stratum", aes(label = after_stat(stratum)), size = 3) +
  scale_x_discrete(
    limits = c("Dominant driver", "Mediation module", "Disease category"),
    expand = c(.03, .03)
  ) +
  labs(
    title = "Dominant driver → mediation modules → disease categories",
    fill = "Driver"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face="bold"),
    axis.title = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    legend.position = "right"
  )
print(p_alluvial_flip)

# ----------------------------
# 6) Protein-level mediation amount plots
# ----------------------------
prot_pm <- DFm %>%
  group_by(prot_term) %>%
  summarise(
    n_links = n(),
    PM_median = median(PM_abs, na.rm = TRUE),
    PM_mean = mean(PM_abs, na.rm = TRUE),
    frac_strong = mean(PM_abs >= 0.5, na.rm = TRUE),
    module = first(module),
    .groups="drop"
  )

# Distribution across proteins
p_prot_hist <- ggplot(prot_pm, aes(PM_median)) +
  geom_histogram(bins = 40) +
  theme_bw() +
  labs(
    title = "Protein-level mediation strength",
    subtitle = "Median PM_abs across significant protein–disease links",
    x = "Median PM_abs per protein",
    y = "Number of proteins"
  )
print(p_prot_hist)

# Top proteins by mediation strength (requires min n_links for stability)
topP <- prot_pm %>%
  filter(n_links >= 10) %>%
  arrange(desc(PM_median)) %>%
  slice_head(n = 30)

p_prot_top <- ggplot(topP, aes(x = reorder(prot_term, PM_median), y = PM_median)) +
  geom_point() +
  coord_flip() +
  theme_bw() +
  labs(
    title = "Top proteins by mediation strength",
    subtitle = "Median PM_abs (only proteins with ≥10 significant links)",
    x = NULL, y = "Median PM_abs"
  )
print(p_prot_top)

# Optional: module summary of mediation strength
mod_pm <- DFm %>%
  group_by(module) %>%
  summarise(
    n_links = n(),
    PM_median = median(PM_abs, na.rm = TRUE),
    frac_strong = mean(PM_abs >= 0.5, na.rm = TRUE),
    .groups="drop"
  ) %>%
  arrange(desc(PM_median))

p_mod_pm <- ggplot(mod_pm, aes(x = reorder(module, PM_median), y = PM_median)) +
  geom_col() +
  coord_flip() +
  theme_bw() +
  labs(
    title = "Mediation strength by module",
    x = NULL, y = "Median PM_abs"
  )
print(p_mod_pm)

# ----------------------------
# 7) Disease-specific mediation strength plots
# ----------------------------
# Per disease (DZ_ID)
dz_pm <- DFm %>%
  group_by(DZ_ID, dz_cat) %>%
  summarise(
    n_links = n(),
    PM_median = median(PM_abs, na.rm = TRUE),
    PM_mean = mean(PM_abs, na.rm = TRUE),
    frac_strong = mean(PM_abs >= 0.5, na.rm = TRUE),
    .groups="drop"
  )

# Top diseases by mediation (requires enough links)
dz_top <- dz_pm %>%
  filter(n_links >= 20) %>%
  arrange(desc(PM_median)) %>%
  slice_head(n = 30)

p_dz_top <- ggplot(dz_top, aes(x = reorder(DZ_ID, PM_median), y = PM_median)) +
  geom_point() +
  coord_flip() +
  theme_bw() +
  labs(
    title = "Diseases with strongest mediation signal",
    subtitle = "Median PM_abs across significant links (only diseases with ≥20 links)",
    x = NULL, y = "Median PM_abs"
  )
print(p_dz_top)

# Per disease category (dz_cat)
dzcat_pm <- DFm %>%
  group_by(dz_cat) %>%
  summarise(
    n_links = n(),
    PM_median = median(PM_abs, na.rm = TRUE),
    PM_mean = mean(PM_abs, na.rm = TRUE),
    frac_strong = mean(PM_abs >= 0.5, na.rm = TRUE),
    .groups="drop"
  ) %>%
  arrange(desc(PM_median))

print(dzcat_pm)

p_dzcat_pm <- ggplot(dzcat_pm, aes(x = reorder(dz_cat, PM_median), y = PM_median)) +
  geom_col() +
  coord_flip() +
  theme_bw() +
  labs(
    title = "Mediation strength differs by disease category",
    subtitle = "Median PM_abs across significant protein–disease links",
    x = NULL, y = "Median PM_abs"
  )
print(p_dzcat_pm)

# ----------------------------
# 8) Module annotation table (for naming modules)
# ----------------------------
mod_driver <- DFm %>%
  mutate(w = if (CFG$alluvial_weight == "strength") strength else 1) %>%
  group_by(module, driver) %>%
  summarise(w = sum(w), .groups="drop") %>%
  group_by(module) %>%
  mutate(frac = w/sum(w)) %>%
  summarise(driver_mix = paste0(driver, "=", round(100*frac), "%", collapse="; "), .groups="drop")

mod_top_dz <- DFm %>%
  count(module, dz_cat, name="n") %>%
  group_by(module) %>%
  arrange(desc(n), .by_group = TRUE) %>%
  slice_head(n = 3) %>%
  summarise(top_dz = paste0(dz_cat, " (", n, ")", collapse="; "), .groups="drop")

mod_top_prot <- DFm %>%
  count(module, prot_term, name="n") %>%
  group_by(module) %>%
  arrange(desc(n), .by_group = TRUE) %>%
  slice_head(n = 8) %>%
  summarise(top_proteins = paste(prot_term, collapse=", "), .groups="drop")

module_summary <- tibble(module = levels(DFa$module)) %>%
  left_join(mod_pm, by="module") %>%
  left_join(mod_driver, by="module") %>%
  left_join(mod_top_dz, by="module") %>%
  left_join(mod_top_prot, by="module")

print(module_summary)

# OPTIONAL: save outputs
# ggsave("alluvial_flip_modules.png", p_alluvial_flip, width = 12, height = 6.5, dpi = 300)
# ggsave("K_silhouette.png", p_k, width = 6, height = 4, dpi = 300)
# ggsave("protein_PM_hist.png", p_prot_hist, width = 6, height = 4, dpi = 300)
# ggsave("protein_PM_top.png", p_prot_top, width = 7, height = 6, dpi = 300)
# ggsave("module_PM.png", p_mod_pm, width = 6, height = 4, dpi = 300)
# ggsave("disease_PM_top.png", p_dz_top, width = 10, height = 7, dpi = 300)
# ggsave("dzcat_PM.png", p_dzcat_pm, width = 6, height = 4, dpi = 300)
# write.csv(module_summary, "mediation_module_summary.csv", row.names = FALSE)
