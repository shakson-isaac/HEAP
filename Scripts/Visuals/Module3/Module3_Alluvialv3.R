# ============================================================
# HEAP Mediation — Publication-ready Sankey/Alluvial (STATIC)
# Driver -> Module -> Disease Category
#
# Fixes vs your last version:
#  - FIXED: module_label not found (no eval_tidy; use rlang::sym + explicit module2)
#  - FIXED: stratum labels no longer sink to bottom
#      * remove manual y=n/2 label layers
#      * use geom_text(stat="stratum") so labels are centered in each stratum
#      * still shows "N (%)" in each stratum by relabeling factor levels
#  - Keeps y-axis visible (your theme does not blank y text/ticks)
#  - Guard: if edge filtering removes everything, returns NULL cleanly
# ============================================================

setwd("/n/groups/patel/shakson_ukb/UK_Biobank/Output/Mediation/")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(ggalluvial)
  library(cluster)
  library(scales)
  library(rlang)
})

# ----------------------------
# 0) Config
# ----------------------------
CFG <- list(
  # significance filter
  alpha_global = 0.05,
  
  # stability filters
  denom_min = 0.01,    # drop numerical dust in |NIE| sum
  te_eps = 0.01,       # TE threshold for PM_signed
  
  # clustering
  do_choose_K = TRUE,
  K_grid = 6:20,
  K_subsample_n = 800,
  K_modules = NULL,    # if NULL uses best by silhouette when do_choose_K=TRUE; else fixed K
  seed = 1,
  
  # alluvial
  alluvial_weight = c("count", "strength")[1], # main = "count"
  keep_top_dz_cat = 10,        # collapse others into "Other"; set NULL to keep all
  min_edge_weight = 5,         # drop edges with fewer than this many links
  rename_modules = TRUE,       # Module 04 (Gcis 71%) style labels
  
  # optional panels
  make_strong_panel = TRUE,    # PM_abs >= 0.5
  strong_cut = 0.5,
  
  # optional extra plots
  make_protein_pm_plots = TRUE,
  make_dzcat_pm_plot = TRUE
)

driver_levels <- c("E-dominant","Gcis-dominant","Gtrans-dominant")

# ----------------------------
# Theme helper (fix clipping)
# ----------------------------
theme_heap <- function(base_size = 14) {
  theme_bw(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0),
      plot.subtitle = element_text(hjust = 0),
      plot.margin = margin(t = 10, r = 30, b = 10, l = 30),
      panel.grid = element_blank(),
      axis.title = element_blank(),
      legend.position = "right",
      legend.title = element_text(face = "bold")
    )
}

# ============================================================
# 1) Build significant link set + PM_abs
# ============================================================
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
    
    NIE_total = NIE_E + NIE_cis + NIE_tr,
    NDE_total = NDE_E + NDE_cis + NDE_tr,
    
    TE_total = if_else(!is.na(total_TE_logHR), total_TE_logHR, NIE_total + NDE_total),
    
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
    
    # evidence strength (optional weighting)
    p_min = pmin(pxs_NIE_logHR_delta_p, gcis_NIE_logHR_delta_p, gtrans_NIE_logHR_delta_p, na.rm = TRUE),
    strength = -log10(pmax(p_min, 1e-300)),
    
    # Proportion mediated (stable, bounded)
    PM_abs = abs(NIE_total) / (abs(NIE_total) + abs(NDE_total) + 1e-12),
    
    # Signed proportion mediated (classical NIE/TE) only if TE not tiny
    PM_signed = if_else(abs(TE_total) >= CFG$te_eps, NIE_total / TE_total, NA_real_)
  )

message("Filtered links: n=", nrow(DFf),
        " | proteins=", n_distinct(DFf$prot_term),
        " | diseases=", n_distinct(DFf$DZ_ID),
        " | dz_cat=", n_distinct(DFf$dz_cat))

# ============================================================
# 2) Protein feature table -> kmeans modules
# ============================================================
prot_feat <- DFf %>%
  group_by(prot_term) %>%
  summarise(
    n_links = n(),
    
    # NIE magnitude per component
    mean_abs_NIE_E   = mean(abs(NIE_E),   na.rm = TRUE),
    mean_abs_NIE_cis = mean(abs(NIE_cis), na.rm = TRUE),
    mean_abs_NIE_tr  = mean(abs(NIE_tr),  na.rm = TRUE),
    
    # signed tendency
    mean_NIE_E   = mean(NIE_E,   na.rm = TRUE),
    mean_NIE_cis = mean(NIE_cis, na.rm = TRUE),
    mean_NIE_tr  = mean(NIE_tr,  na.rm = TRUE),
    
    # driver fractions
    frac_E_dom   = mean(driver == "E-dominant"),
    frac_cis_dom = mean(driver == "Gcis-dominant"),
    frac_tr_dom  = mean(driver == "Gtrans-dominant"),
    
    # mediation strength
    mean_PM_abs = mean(PM_abs, na.rm = TRUE),
    median_PM_abs = median(PM_abs, na.rm = TRUE),
    frac_PM_abs_gt50 = mean(PM_abs >= 0.5, na.rm = TRUE),
    
    # evidence
    mean_strength = mean(strength, na.rm = TRUE),
    
    .groups = "drop"
  )

X <- prot_feat %>% select(-prot_term) %>% as.matrix()
Xz <- scale(X)

# ---- Choose K (optional) ----
if (CFG$do_choose_K) {
  set.seed(CFG$seed)
  nP <- nrow(Xz)
  idx <- sample(seq_len(nP), size = min(CFG$K_subsample_n, nP))
  Xsub <- Xz[idx, , drop = FALSE]
  
  Ks <- CFG$K_grid
  sil <- numeric(length(Ks))
  
  set.seed(CFG$seed)
  for (i in seq_along(Ks)) {
    k <- Ks[i]
    km <- kmeans(Xsub, centers = k, nstart = 25, iter.max = 200)
    sil[i] <- mean(silhouette(km$cluster, dist(Xsub))[, 3])
  }
  
  k_tbl <- tibble(K = Ks, silhouette = sil)
  print(k_tbl)
  
  p_k <- ggplot(k_tbl, aes(K, silhouette)) +
    geom_line() + geom_point() +
    theme_heap(14) +
    coord_cartesian(clip = "off") +
    labs(title = "Silhouette vs K (kmeans on protein NIE/PM features)",
         y = "Mean silhouette (subsample)")
  print(p_k)
  ggsave("K_silhouette.png", p_k, width = 12, height = 4, dpi = 300, bg = "white")
  
  bestK <- k_tbl %>% arrange(desc(silhouette)) %>% slice(1) %>% pull(K)
  if (is.null(CFG$K_modules)) CFG$K_modules <- bestK
  message("Using K = ", CFG$K_modules, " (best by silhouette = ", bestK, ")")
} else {
  if (is.null(CFG$K_modules)) CFG$K_modules <- 12
}

# ---- Cluster all proteins ----
set.seed(CFG$seed)
km_full <- kmeans(Xz, centers = CFG$K_modules, nstart = 50, iter.max = 200)

prot2mod <- prot_feat %>%
  mutate(
    module_id = km_full$cluster,
    module = sprintf("Module %02d", module_id)
  ) %>%
  select(prot_term, module_id, module)

DFm <- DFf %>% inner_join(prot2mod, by = "prot_term")

# ============================================================
# 3) Helper: build refined alluvial from a link-level DF
# ============================================================
build_alluvial_refined <- function(DF_link, file_out, title, subtitle_extra = "") {
  
  # Weight definition
  DFa0 <- DF_link %>%
    mutate(w_all = if (CFG$alluvial_weight == "strength") strength else 1)
  
  # Collapse disease categories
  if (!is.null(CFG$keep_top_dz_cat)) {
    keep_cats <- DFa0 %>%
      count(dz_cat, wt = w_all, name = "tw") %>%
      arrange(desc(tw)) %>%
      slice_head(n = CFG$keep_top_dz_cat) %>%
      pull(dz_cat)
    DFa0 <- DFa0 %>% mutate(dz_cat = if_else(dz_cat %in% keep_cats, dz_cat, "Other"))
  }
  
  # Edge aggregation + drop tiny flows
  DFa_edge <- DFa0 %>%
    count(driver, module, dz_cat, wt = w_all, name = "weight") %>%
    filter(weight >= CFG$min_edge_weight)
  
  if (nrow(DFa_edge) == 0) {
    warning("No edges left after filtering. Try lowering CFG$min_edge_weight.")
    return(invisible(NULL))
  }
  
  # Factor ordering
  DFa_edge <- DFa_edge %>%
    mutate(driver = factor(driver, levels = driver_levels))
  
  # Module order to reduce crossings (group by dominant driver)
  mod_order <- DFa_edge %>%
    group_by(module, driver) %>% summarise(w = sum(weight), .groups="drop") %>%
    group_by(module) %>%
    summarise(
      wE = sum(w[driver=="E-dominant"]),
      wC = sum(w[driver=="Gcis-dominant"]),
      wT = sum(w[driver=="Gtrans-dominant"]),
      dom = driver_levels[which.max(c(wE,wC,wT))],
      dom_w = max(c(wE,wC,wT)),
      .groups="drop"
    ) %>%
    arrange(factor(dom, levels = driver_levels), desc(dom_w)) %>%
    pull(module)
  
  dz_order <- DFa_edge %>%
    group_by(dz_cat) %>% summarise(tw = sum(weight), .groups="drop") %>%
    arrange(desc(tw)) %>%
    pull(dz_cat)
  
  DFa_edge <- DFa_edge %>%
    mutate(
      module = factor(module, levels = mod_order),
      dz_cat = factor(dz_cat, levels = dz_order)
    )
  
  # Optional module renaming: "Module 04 (Gcis 71%)"
  if (CFG$rename_modules) {
    mod_labels <- DFa_edge %>%
      group_by(module, driver) %>% summarise(w = sum(weight), .groups="drop") %>%
      group_by(module) %>%
      mutate(frac = w/sum(w)) %>%
      summarise(
        dom = as.character(driver[which.max(frac)]),
        dom_frac = max(frac),
        .groups="drop"
      ) %>%
      mutate(module_label = paste0(as.character(module), " (", dom, " ", percent(dom_frac, accuracy=1), ")"))
    
    DFa_edge <- DFa_edge %>%
      left_join(mod_labels %>% select(module, module_label), by="module")
    
    # preserve ordering
    label_order <- mod_labels$module_label[match(mod_order, mod_labels$module)]
    DFa_edge <- DFa_edge %>%
      mutate(module_label = factor(module_label, levels = label_order))
  }
  
  # ---- Define axis2 variable explicitly (NO eval_tidy) ----
  DFa_edge <- DFa_edge %>%
    mutate(module2 = if (CFG$rename_modules) module_label else module)
  
  # ---- Build "N (%)" stratum labels by relabeling factor levels ----
  total_w <- sum(DFa_edge$weight)
  
  # driver labels
  drv_tbl <- DFa_edge %>% count(driver, wt = weight, name = "n")
  drv_lv  <- driver_levels
  drv_n   <- drv_tbl$n[match(drv_lv, as.character(drv_tbl$driver))]
  drv_n[is.na(drv_n)] <- 0
  drv_lab <- setNames(
    paste0(drv_lv, "\n", drv_n, " (", percent(drv_n / total_w, accuracy = 1), ")"),
    drv_lv
  )
  
  # module labels
  mod_tbl <- DFa_edge %>% count(module2, wt = weight, name = "n")
  mod_lv  <- levels(factor(DFa_edge$module2))
  mod_n   <- mod_tbl$n[match(mod_lv, as.character(mod_tbl$module2))]
  mod_n[is.na(mod_n)] <- 0
  mod_lab <- setNames(
    paste0(mod_lv, "\n", mod_n, " (", percent(mod_n / total_w, accuracy = 1), ")"),
    mod_lv
  )
  
  # disease category labels
  dz_tbl <- DFa_edge %>% count(dz_cat, wt = weight, name = "n")
  dz_lv  <- levels(DFa_edge$dz_cat)
  dz_n   <- dz_tbl$n[match(dz_lv, as.character(dz_tbl$dz_cat))]
  dz_n[is.na(dz_n)] <- 0
  dz_lab <- setNames(
    paste0(dz_lv, "\n", dz_n, " (", percent(dz_n / total_w, accuracy = 1), ")"),
    dz_lv
  )
  
  DFa_edge <- DFa_edge %>%
    mutate(
      driver  = factor(as.character(driver),  levels = drv_lv, labels = drv_lab[drv_lv]),
      module2 = factor(as.character(module2), levels = mod_lv, labels = mod_lab[mod_lv]),
      dz_cat  = factor(as.character(dz_cat),  levels = dz_lv,  labels = dz_lab[dz_lv])
    )
  
  # Plot
  p <- ggplot(
    DFa_edge,
    aes(axis1 = driver, axis2 = module2, axis3 = dz_cat, y = weight)
  ) +
    geom_alluvium(aes(fill = driver), alpha = 0.85, width = 1/14) +
    geom_stratum(width = 1/10, fill = "white", color = "black") +
    # Centered labels per stratum (fixes "all text goes to bottom")
    geom_text(stat = "stratum", aes(label = after_stat(stratum)),
              size = 3, lineheight = 0.9) +
    scale_x_discrete(limits = c("Dominant driver", "Mediation module", "Disease category"),
                     expand = c(.04, .04)) +
    labs(
      title = title,
      subtitle = paste0(
        "Flow width = ",
        if (CFG$alluvial_weight == "strength") "evidence strength (−log10 p)" else "# significant protein–disease links",
        " | edges < ", CFG$min_edge_weight, " removed",
        ifelse(subtitle_extra == "", "", paste0(" | ", subtitle_extra))
      ),
      fill = "Driver"
    ) +
    coord_cartesian(clip = "off") +
    theme_heap(14)
  
  print(p)
  ggsave(file_out, p, width = 16, height = 8, dpi = 300, bg = "white")
  invisible(p)
}

# ============================================================
# 4) Main refined alluvial
# ============================================================
p_main <- build_alluvial_refined(
  DF_link = DFm,
  file_out = "alluvial_driver_module_dz_refined.png",
  title = "Dominant driver → mediation modules → disease categories"
)

# ============================================================
# 5) Optional: strong mediation-only alluvial (PM_abs >= cut)
# ============================================================
if (CFG$make_strong_panel) {
  DFm_strong <- DFm %>% filter(PM_abs >= CFG$strong_cut)
  
  p_strong <- build_alluvial_refined(
    DF_link = DFm_strong,
    file_out = "alluvial_driver_module_dz_refined_strong.png",
    title = "Dominant driver → mediation modules → disease categories (strongly mediated links)",
    subtitle_extra = paste0("PM_abs ≥ ", CFG$strong_cut)
  )
}

# ============================================================
# 6) Optional: protein-level mediation amount plots
# ============================================================
if (CFG$make_protein_pm_plots) {
  
  prot_pm <- DFm %>%
    group_by(prot_term) %>%
    summarise(
      n_links = n(),
      PM_median = median(PM_abs, na.rm = TRUE),
      PM_mean = mean(PM_abs, na.rm = TRUE),
      frac_strong = mean(PM_abs >= 0.5, na.rm = TRUE),
      .groups="drop"
    )
  
  p_prot_hist <- ggplot(prot_pm, aes(PM_median)) +
    geom_histogram(bins = 40) +
    coord_cartesian(clip = "off") +
    theme_heap(14) +
    labs(
      title = "Protein-level mediation strength",
      subtitle = "Median PM_abs across significant protein–disease links",
      x = "Median PM_abs per protein",
      y = "Number of proteins"
    )
  print(p_prot_hist)
  ggsave("protein_PM_hist.png", p_prot_hist, width = 12, height = 5, dpi = 300, bg = "white")
  
  topP <- prot_pm %>%
    filter(n_links >= 10) %>%
    arrange(desc(PM_median)) %>%
    slice_head(n = 30)
  
  p_prot_top <- ggplot(topP, aes(x = reorder(prot_term, PM_median), y = PM_median)) +
    geom_point() +
    coord_flip(clip = "off") +
    theme_heap(14) +
    labs(
      title = "Top proteins by mediation strength",
      subtitle = "Median PM_abs (proteins with ≥10 significant links)",
      x = NULL, y = "Median PM_abs"
    )
  print(p_prot_top)
  ggsave("protein_PM_top30.png", p_prot_top, width = 10, height = 9, dpi = 300, bg = "white")
}

# ============================================================
# 7) Optional: disease-category mediation strength plot
# ============================================================
if (CFG$make_dzcat_pm_plot) {
  
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
  
  p_dzcat <- ggplot(dzcat_pm, aes(x = reorder(dz_cat, PM_median), y = PM_median)) +
    geom_col() +
    coord_flip(clip = "off") +
    theme_heap(14) +
    labs(
      title = "Mediation strength differs by disease category",
      subtitle = "Median PM_abs across significant protein–disease links",
      x = NULL, y = "Median PM_abs"
    )
  print(p_dzcat)
  ggsave("dzcat_PM.png", p_dzcat, width = 10, height = 6, dpi = 300, bg = "white")
}

message("Done. Files saved in working directory:
- alluvial_driver_module_dz_refined.png
- alluvial_driver_module_dz_refined_strong.png (if enabled)
- protein_PM_hist.png, protein_PM_top30.png (if enabled)
- dzcat_PM.png (if enabled)
- K_silhouette.png (if enabled)")

