# ============================================================
# HEAP Mediation — Publication-ready Alluvial + Protein-centric views
#
# MAIN FIGURE (unweighted, link-mass):
#   Dominant driver (per link) -> Protein module (kmeans, biology-ish) -> Disease category
#
# PROTEIN-CENTRIC VIEWS:
#   (B) Heatmap: Top hub proteins (by # links) x disease categories (counts)
#   (C) Composition plot: per-protein mediated composition (E/cis/trans) with hub size
#
# Notes:
# - The alluvial is intentionally UNWEIGHTED so proteins with many disease links dominate ribbon mass.
# - Modules are defined WITHOUT explicit driver-fraction features (avoids tautological driver buckets).
# - Driver on axis1 is link-level (per protein–disease link), based on max(|NIE_E|,|NIE_cis|,|NIE_tr|).
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
  denom_min = 0.01,   # drop numerical dust in |NIE| sum
  te_eps = 0.01,      # TE threshold for PM_signed
  
  # clustering
  do_choose_K = TRUE,
  K_grid = 6:20,
  K_subsample_n = 800,
  K_modules = NULL,   # if NULL uses best by silhouette when do_choose_K=TRUE; else fixed K
  seed = 1,
  
  # alluvial
  alluvial_weight = c("count", "strength")[1], # "count" recommended for main
  keep_top_dz_cat = 10,     # collapse others into "Other"; set NULL to keep all
  min_edge_weight = 5,      # drop edges with fewer than this many links (or strength mass)
  rename_modules = FALSE,    # label modules with dominant composition fraction
  
  # optional panels
  make_strong_panel = TRUE,
  strong_cut = 0.5,
  
  # protein-centric plots
  make_protein_hub_table = TRUE,
  make_heatmap = TRUE,
  heatmap_topN = 60,        # top proteins by n_links
  heatmap_log1p = TRUE,     # log1p(count) for fill
  
  make_composition_plot = TRUE,
  comp_topN = 500           # top proteins by n_links for composition scatter
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
# 1) Build significant link set + PM_abs + link-level driver
# ============================================================
DF <- Type5 %>% as_tibble()
DF <- Type5 %>% as_tibble()

# Force key columns to be plain numerics (handles list-cols safely)
num_cols <- c(
  "pxs_NIE_logHR_delta_p","gcis_NIE_logHR_delta_p","gtrans_NIE_logHR_delta_p",
  "pxs_NIE_logHR","gcis_NIE_logHR","gtrans_NIE_logHR",
  "pxs_NDE_logHR","gcis_NDE_logHR","gtrans_NDE_logHR",
  "total_TE_logHR"
)

DF <- DF %>%
  mutate(
    across(any_of(num_cols), ~ suppressWarnings(as.numeric(unlist(.x)))),
    DZ_ID = as.character(DZ_ID),
    prot_term = as.character(prot_term)
  )



alpha <- CFG$alpha_global / nrow(DF)

DFf <- DF %>%
  filter(
    pxs_NIE_logHR_delta_p      < alpha |
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
    # dominant driver per LINK (based on |NIE|)
    driver = case_when(
      abs(NIE_E)   >= abs(NIE_cis) & abs(NIE_E)   >= abs(NIE_tr) ~ "E-dominant",
      abs(NIE_cis) >= abs(NIE_E)   & abs(NIE_cis) >= abs(NIE_tr) ~ "Gcis-dominant",
      TRUE ~ "Gtrans-dominant"
    ),
    
    # evidence strength (optional weighting)
    p_min = pmin(pxs_NIE_logHR_delta_p, gcis_NIE_logHR_delta_p, gtrans_NIE_logHR_delta_p, na.rm = TRUE),
    #strength = -log10(pmax(p_min, 1e-300))
    strength = as.numeric(-log10(pmax(p_min, 1e-300))),
    
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
#    (No explicit driver-fraction variables used for clustering)
# ============================================================
prot_feat <- DFf %>%
  group_by(prot_term) %>%
  summarise(
    n_links = n(),
    
    # mediated composition (pattern)
    frac_abs_E   = mean(abs(NIE_E)   / (abs(NIE_E)+abs(NIE_cis)+abs(NIE_tr) + 1e-12), na.rm = TRUE),
    frac_abs_cis = mean(abs(NIE_cis) / (abs(NIE_E)+abs(NIE_cis)+abs(NIE_tr) + 1e-12), na.rm = TRUE),
    frac_abs_tr  = mean(abs(NIE_tr)  / (abs(NIE_E)+abs(NIE_cis)+abs(NIE_tr) + 1e-12), na.rm = TRUE),
    
    # mediation strength + evidence
    mean_PM_abs   = mean(PM_abs, na.rm = TRUE),
    median_PM_abs = median(PM_abs, na.rm = TRUE),
    mean_strength = mean(strength, na.rm = TRUE),
    
    # signed tendencies (optional, helps separate + vs -)
    mean_NIE_E   = mean(NIE_E,   na.rm = TRUE),
    mean_NIE_cis = mean(NIE_cis, na.rm = TRUE),
    mean_NIE_tr  = mean(NIE_tr,  na.rm = TRUE),
    
    # magnitude (helps split "big" vs "small" proteins)
    mean_abs_NIE_total = mean(abs(NIE_total), na.rm = TRUE),
    mean_abs_NDE_total = mean(abs(NDE_total), na.rm = TRUE),
    
    .groups = "drop"
  )

# clustering matrix
X <- prot_feat %>%
  select(
    frac_abs_E, frac_abs_cis, frac_abs_tr,
    mean_PM_abs, median_PM_abs,
    mean_strength,
    mean_NIE_E, mean_NIE_cis, mean_NIE_tr,
    mean_abs_NIE_total, mean_abs_NDE_total
  ) %>% as.matrix()

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
    labs(title = "Silhouette vs K (protein modules from composition/PM/evidence)",
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
  select(prot_term, module_id, module,
         n_links,
         frac_abs_E, frac_abs_cis, frac_abs_tr,
         mean_PM_abs, median_PM_abs, mean_strength)

# add link-level module assignments
DFm <- DFf %>% inner_join(prot2mod %>% select(prot_term, module_id, module, n_links), by = "prot_term")

# ============================================================
# 3) Helper: refined alluvial (UNWEIGHTED / link-mass)
# ============================================================
build_alluvial_refined <- function(DF_link, file_out, title, subtitle_extra = "") {
  
  # ---- force atomic columns (this kills list-cols) ----
  DFa0 <- DF_link %>%
    transmute(
      driver = as.character(unlist(driver)),
      module = as.character(unlist(module)),
      dz_cat = as.character(unlist(dz_cat)),
      strength = suppressWarnings(as.numeric(unlist(strength)))
    )
  
  DFa0$strength[is.na(DFa0$strength)] <- 0
  
  # weight: either 1 per link (count) or -log10(p) strength
  if (CFG$alluvial_weight == "strength") {
    DFa0$w_all <- DFa0$strength
  } else {
    DFa0$w_all <- 1
  }
  DFa0$w_all <- as.numeric(DFa0$w_all)
  DFa0$w_all[is.na(DFa0$w_all)] <- 0
  
  # ---- collapse disease categories (top N by total weight) ----
  if (!is.null(CFG$keep_top_dz_cat)) {
    dz_tot <- aggregate(w_all ~ dz_cat, data = DFa0, sum)
    dz_tot <- dz_tot[order(dz_tot$w_all, decreasing = TRUE), , drop = FALSE]
    keep_cats <- head(dz_tot$dz_cat, CFG$keep_top_dz_cat)
    DFa0$dz_cat <- ifelse(DFa0$dz_cat %in% keep_cats, DFa0$dz_cat, "Other")
  }
  
  # ---- aggregate edges with base R (no dplyr count) ----
  edge <- aggregate(w_all ~ driver + module + dz_cat, data = DFa0, sum)
  names(edge)[names(edge) == "w_all"] <- "weight"
  
  # drop tiny flows
  edge <- edge[edge$weight >= CFG$min_edge_weight, , drop = FALSE]
  if (nrow(edge) == 0) {
    warning("No edges left after filtering. Lower CFG$min_edge_weight.")
    return(invisible(NULL))
  }
  
  # ---- factor ordering ----
  edge$driver <- factor(edge$driver, levels = driver_levels)
  
  # module order: by dominant driver mass
  mod_tab <- aggregate(weight ~ module + driver, data = edge, sum)
  # wide-ish manual dominance
  getw <- function(mod, drv) {
    ii <- mod_tab$module == mod & as.character(mod_tab$driver) == drv
    if (any(ii)) mod_tab$weight[ii][1] else 0
  }
  mods <- unique(edge$module)
  dom_drv <- sapply(mods, function(m) {
    wE <- getw(m, "E-dominant")
    wC <- getw(m, "Gcis-dominant")
    wT <- getw(m, "Gtrans-dominant")
    driver_levels[which.max(c(wE, wC, wT))]
  })
  dom_w <- sapply(mods, function(m) {
    wE <- getw(m, "E-dominant")
    wC <- getw(m, "Gcis-dominant")
    wT <- getw(m, "Gtrans-dominant")
    max(c(wE, wC, wT))
  })
  mod_order <- mods[order(match(dom_drv, driver_levels), -dom_w)]
  
  # disease order by total
  dz_tot2 <- aggregate(weight ~ dz_cat, data = edge, sum)
  dz_order <- dz_tot2$dz_cat[order(dz_tot2$weight, decreasing = TRUE)]
  
  edge$module <- factor(edge$module, levels = mod_order)
  edge$dz_cat <- factor(edge$dz_cat, levels = dz_order)
  
  # ---- optional module rename (still safe; uses prot2mod, which is clean) ----
  if (CFG$rename_modules) {
    mod_comp <- prot2mod %>%
      group_by(module) %>%
      summarise(
        fracE = mean(frac_abs_E, na.rm = TRUE),
        fracC = mean(frac_abs_cis, na.rm = TRUE),
        fracT = mean(frac_abs_tr, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(
        dom = ifelse(fracE >= fracC & fracE >= fracT, "E",
                     ifelse(fracC >= fracE & fracC >= fracT, "cis", "trans")),
        dom_frac = pmax(fracE, fracC, fracT),
        module_label = paste0(module, " (", dom, " ", scales::percent(dom_frac, accuracy = 1), ")")
      )
    
    edge <- merge(edge, mod_comp[, c("module", "module_label")], by = "module", all.x = TRUE)
    edge$module2 <- factor(edge$module_label, levels = mod_comp$module_label[match(mod_order, mod_comp$module)])
  } else {
    edge$module2 <- edge$module
  }
  
  # ---- plot ----
  p <- ggplot(edge, aes(axis1 = driver, axis2 = module2, axis3 = dz_cat, y = weight)) +
    geom_alluvium(aes(fill = driver), alpha = 0.85, width = 1/14) +
    geom_stratum(width = 1/10, fill = "white", color = "black") +
    geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3, lineheight = 0.9) +
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
# 4) Main refined alluvial (UNWEIGHTED)
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
# 6) Protein-centric summary table (hub proteins)
# ============================================================
prot_summary <- DFm %>%
  group_by(prot_term) %>%
  summarise(
    n_links = n(),
    n_dz = n_distinct(DZ_ID),
    n_dz_cat = n_distinct(dz_cat),
    
    # driver composition across links
    frac_E_dom   = mean(driver == "E-dominant"),
    frac_cis_dom = mean(driver == "Gcis-dominant"),
    frac_tr_dom  = mean(driver == "Gtrans-dominant"),
    
    # mediated composition across links (absolute)
    frac_abs_E   = mean(abs(NIE_E)   / (abs(NIE_E)+abs(NIE_cis)+abs(NIE_tr) + 1e-12), na.rm = TRUE),
    frac_abs_cis = mean(abs(NIE_cis) / (abs(NIE_E)+abs(NIE_cis)+abs(NIE_tr) + 1e-12), na.rm = TRUE),
    frac_abs_tr  = mean(abs(NIE_tr)  / (abs(NIE_E)+abs(NIE_cis)+abs(NIE_tr) + 1e-12), na.rm = TRUE),
    
    PM_median = median(PM_abs, na.rm = TRUE),
    PM_mean   = mean(PM_abs, na.rm = TRUE),
    strength_mean = mean(strength, na.rm = TRUE),
    
    module = first(module),
    .groups = "drop"
  ) %>%
  mutate(
    driver_prot = case_when(
      frac_abs_E   >= frac_abs_cis & frac_abs_E   >= frac_abs_tr ~ "E-dominant",
      frac_abs_cis >= frac_abs_E   & frac_abs_cis >= frac_abs_tr ~ "Gcis-dominant",
      TRUE ~ "Gtrans-dominant"
    )
  ) %>%
  arrange(desc(n_links))

if (CFG$make_protein_hub_table) {
  write.csv(prot_summary, "protein_hub_summary.csv", row.names = FALSE)
}

# ============================================================
# 7B) Protein-centric Heatmap: Top proteins x disease categories
# ============================================================
if (CFG$make_heatmap) {
  
  topP <- prot_summary %>% slice_head(n = CFG$heatmap_topN) %>% pull(prot_term)
  
  heat_df <- DFm %>%
    transmute(
      prot_term = as.character(unlist(prot_term)),
      dz_cat    = as.character(unlist(dz_cat))
    )
  
  heat_df <- aggregate(
    x = list(n_links = rep(1, nrow(heat_df))),
    by = list(prot_term = heat_df$prot_term, dz_cat = heat_df$dz_cat),
    FUN = sum
  ) %>% as_tibble()
  
  # Order rows by total links
  prot_order <- heat_df %>%
    group_by(prot_term) %>% summarise(tw = sum(n_links), .groups = "drop") %>%
    arrange(desc(tw)) %>% pull(prot_term)
  
  dz_order <- heat_df %>%
    group_by(dz_cat) %>% summarise(tw = sum(n_links), .groups = "drop") %>%
    arrange(desc(tw)) %>% pull(dz_cat)
  
  heat_df <- heat_df %>%
    mutate(
      prot_term = factor(prot_term, levels = rev(prot_order)),
      dz_cat = factor(dz_cat, levels = dz_order),
      fillv = if (CFG$heatmap_log1p) log1p(n_links) else n_links
    )
  
  p_heat <- ggplot(heat_df, aes(x = dz_cat, y = prot_term, fill = fillv)) +
    geom_tile(color = "white", linewidth = 0.15) +
    theme_heap(12) +
    coord_cartesian(clip = "off") +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "right"
    ) +
    labs(
      title = paste0("Top ", CFG$heatmap_topN, " proteins by breadth: disease-category heatmap"),
      subtitle = if (CFG$heatmap_log1p) "Fill = log1p(# significant links)" else "Fill = # significant links",
      x = "Disease category", y = "Protein"
    )
  print(p_heat)
  ggsave("protein_hubs_heatmap_dzcat.png", p_heat, width = 12, height = 0.18*CFG$heatmap_topN + 4, dpi = 300, bg = "white")
}

# ============================================================
# 7C) Protein-centric Composition plot: mediated composition (E/cis/trans)
#     Simplex in 2D: x=frac_abs_E, y=frac_abs_cis, trans = 1-x-y
# ============================================================
# ============================================================
# 7C) Protein-centric Composition plot (TRUE ternary in an equilateral triangle)
# ============================================================
if (CFG$make_composition_plot) {
  
  comp_df <- prot_summary %>%
    slice_head(n = min(CFG$comp_topN, nrow(prot_summary))) %>%
    mutate(
      E = pmin(pmax(frac_abs_E, 0), 1),
      C = pmin(pmax(frac_abs_cis, 0), 1),
      T = pmax(0, 1 - E - C)
    ) %>%
    # renormalize just in case numerical jitter makes sums != 1
    mutate(s = E + C + T, E = E/s, C = C/s, T = T/s) %>%
    filter(is.finite(E), is.finite(C), is.finite(T))
  
  # ---- barycentric (E,C,T) -> equilateral triangle coordinates ----
  # corners: Trans=(0,0), E=(1,0), Cis=(0.5, sqrt(3)/2)
  b2xy <- function(E, C, T) {
    x <- E + 0.5 * C
    y <- (sqrt(3)/2) * C
    tibble(x = x, y = y)
  }
  
  tri_xy <- b2xy(comp_df$E, comp_df$C, comp_df$T)
  comp_df <- bind_cols(comp_df, tri_xy)
  
  # ---- helper to make grid lines for constant E/C/T ----
  make_grid <- function(vals = seq(0.2, 0.8, 0.2), n = 100) {
    lines <- list()
    
    for (v in vals) {
      # constant E = v, C from 0..(1-v)
      Cseq <- seq(0, 1 - v, length.out = n)
      Eseq <- rep(v, n)
      Tseq <- 1 - Eseq - Cseq
      lines[[length(lines) + 1]] <- b2xy(Eseq, Cseq, Tseq) %>% mutate(which = "E", val = v)
      
      # constant C = v, E from 0..(1-v)
      Eseq <- seq(0, 1 - v, length.out = n)
      Cseq <- rep(v, n)
      Tseq <- 1 - Eseq - Cseq
      lines[[length(lines) + 1]] <- b2xy(Eseq, Cseq, Tseq) %>% mutate(which = "C", val = v)
      
      # constant T = v, E from 0..(1-v)
      Eseq <- seq(0, 1 - v, length.out = n)
      Tseq <- rep(v, n)
      Cseq <- 1 - Eseq - Tseq
      lines[[length(lines) + 1]] <- b2xy(Eseq, Cseq, Tseq) %>% mutate(which = "T", val = v)
    }
    
    bind_rows(lines) %>% mutate(group = paste(which, val))
  }
  
  grid_df <- make_grid()
  
  # triangle boundary
  bound <- tibble(
    x = c(0, 1, 0.5, 0),
    y = c(0, 0, sqrt(3)/2, 0)
  )
  
  # nicer size scaling (avoid giant bubbles)
  comp_df <- comp_df %>%
    mutate(size_plot = n_links)
  
  p_comp <- ggplot() +
    geom_path(data = bound, aes(x = x, y = y), linewidth = 0.7) +
    geom_path(data = grid_df, aes(x = x, y = y, group = group),
              linewidth = 0.25, alpha = 0.25) +
    geom_point(
      data = comp_df,
      aes(x = x, y = y, color = driver_prot, size = size_plot),
      alpha = 0.55
    ) +
    scale_size_continuous(
      name   = "# links",
      range  = c(0.8, 6),
      breaks = c(5, 10, 20, 40, 60)
    ) +
    coord_equal(clip = "off") +
    theme_heap(14) +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank()
    ) +
    labs(
      title = "Protein-level mediated composition (ternary)",
      subtitle = paste0(
        "Top ", min(CFG$comp_topN, nrow(prot_summary)),
        " proteins by # links. Point size ~ sqrt(# links)."
      ),
      color = "Protein driver",
      size  = "# links"
    ) +
    guides(size = guide_legend(override.aes = list(alpha = 1))) +
    # corner labels (outside the triangle)
    annotate("text", x = 1.03, y = -0.03, label = "E = 1", hjust = 1, vjust = 1) +
    annotate("text", x = -0.03, y = -0.03, label = "trans = 1", hjust = 0, vjust = 1) +
    annotate("text", x = 0.50, y = sqrt(3)/2 + 0.03, label = "cis = 1", hjust = 0.5, vjust = 0)
  
  print(p_comp)
  ggsave("protein_composition_ternary.png", p_comp, width = 10, height = 8, dpi = 300, bg = "white")
}

message("Done. Files saved in working directory:
- alluvial_driver_module_dz_refined.png
- alluvial_driver_module_dz_refined_strong.png (if enabled)
- protein_hub_summary.csv (if enabled)
- protein_hubs_heatmap_dzcat.png (if enabled)
- protein_composition_simplex.png (if enabled)
- K_silhouette.png (if enabled)")



