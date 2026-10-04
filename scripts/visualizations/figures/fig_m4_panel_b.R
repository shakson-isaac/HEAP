#!/usr/bin/env Rscript
# ============================================================================
# fig_m4_panel_b.R  -- Fig6 composite panel b: interventional concordance heatmap
# ----------------------------------------------------------------------------
# Which exposures' HEAP proteomic signatures replicate proven interventions.
# Rows = top exposures by concordance breadth, cols = the 3 trials (HERITAGE
# exercise; GLP1RA STEP1/STEP2). Fill = weighted r; * = BH p<0.05 (n_eff>=8).
# A focused, cell-rendered cut of fig_intervention_compare for the composite.
#
# Renders two ways from this one file (per MULTIPANEL_FIGURE_GUIDE):
#   standalone (default)  -- full size for solo review
#   CELL (HEAP_CELL=1)    -- exact cell size, 5-7pt text, *_cell.{png,pdf}
#
# Input : support/intervention_compare/intervention_correlations.tsv
# Run   : HEAP_PATHS_FILE=.../00_paths.R HEAP_CELL=1 \
#           Rscript scripts/visualizations/figures/fig_m4_panel_b.R
# ============================================================================
local({
  cand <- c(file.path(getwd(), "scripts", "visualizations", "common"),
            "/n/groups/patel/shakson_ukb/HEAP/scripts/visualizations/common")
  cm <- cand[dir.exists(cand)][1]
  for (f in c("figure_paths.R","load_heap_results.R","plot_theme.R",
              "label_helpers.R","export_helpers.R")) source(file.path(cm, f))
})
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })

# --- render profile ---------------------------------------------------------
# Three render modes:
#   CELL  -- composite cell for Fig6 (exact size, 5-7pt)
#   FULL  -- one-PAGE supplementary landscape (HEAP_PANELB_MR=0, no cap): ~70 rows
#            must fit a 6.5 x 7.4 in portrait page, so use tight ~6.5pt text
#   else  -- MR-gated panel_b standalone (few rows, roomy large text)
CELL <- nzchar(Sys.getenv("HEAP_CELL"))
FULL <- !CELL && Sys.getenv("HEAP_PANELB_MR", "1") == "0" &&
        !nzchar(Sys.getenv("HEAP_PANELB_OUT", unset = ""))
BS   <- if (CELL) 7   else if (FULL) 7   else 10
STAR <- if (CELL) 3.0 else if (FULL) 2.6 else 4.5
TTL  <- if (CELL) 9   else if (FULL) 11  else 13
B_W  <- 4.50; B_H <- 3.5                          # widened 4.7 -> 5.5 so the exposure row labels fit (author, 2026-08-30)
sel_covar <- "base"

# Output basename: the standalone full landscape (HEAP_PANELB_MR=0, no per-
# category cap) is the supplementary figure fig_intervention_concordance_full;
# the MR-gated cell cut keeps the fig_m4_panel_b name for the composite. Allow an
# explicit HEAP_PANELB_OUT override for ad-hoc runs.
OUTNAME <- Sys.getenv("HEAP_PANELB_OUT", unset = "")

FIGDIR <- file.path(Sys.getenv("HEAP_PROJECT_ROOT", "/n/groups/patel/IGLOO/UKB/HEAP"), "figures", "exploratory/module4")
dir.create(FIGDIR, recursive = TRUE, showWarnings = FALSE)

.read <- function(file) {
  ov <- Sys.getenv("HEAP_INTERVENTION_DIR", unset = "")
  p <- if (nzchar(ov)) file.path(ov, file)
       else heap_project_output("support", "intervention_compare", file)
  if (!file.exists(p)) stop("Missing intervention input: ", p, call. = FALSE)
  fread(p)
}
# --- labels: curated short forms keyed on the stripped base stem -------------
SHORT <- c(
  number_of_days_week_of_vigorous_physical_activity_10_plus_minutes = "Vigorous activity",
  number_of_days_week_of_moderate_physical_activity_10_plus_minutes = "Moderate activity",
  types_of_physical_activity_in_last_4_weeks = "Strenuous sports",
  summed_days_activity = "Active days/week",
  met_minutes_per_week_for_vigorous_activity = "Vigorous MET-min/wk",
  at_or_above_moderate_vigorous_recommendation = "Activity guide",
  at_or_above_moderate_vigorous_walking_recommendation = "Walking guide",
  frequency_of_stair_climbing_in_last_4_weeks = "Stair climbing",
  ipaq_activity_group = "IPAQ group", usual_walking_pace = "Walking pace",
  types_of_transport_used_excluding_work = "Transport type",
  processed_meat_intake = "Processed meat", beef_intake = "Beef intake",
  pork_intake = "Pork intake", poultry_intake = "Poultry intake",
  lamb_mutton_intake = "Lamb/mutton intake",
  oily_fish_intake = "Oily fish", non_oily_fish_intake = "Non-oily fish",
  dried_fruit_intake = "Dried fruit", fresh_fruit_intake = "Fresh fruit",
  cooked_vegetable_intake = "Cooked vegetable", salad_raw_vegetable_intake = "Raw veg/salad",
  tea_intake = "Tea intake", coffee_intake = "Coffee intake", water_intake = "Water intake",
  bread_intake = "Bread intake", cereal_intake = "Cereal intake",
  milk_type_used = "Milk type", spread_type = "Spread type",
  major_dietary_changes_in_the_last_5_years = "Diet change",
  mineral_and_other_dietary_supplements = "Suppl.",
  salt_added_to_food = "Salt added to food", variation_in_diet = "Diet variation",
  never_eat_eggs_dairy_wheat_sugar = "Avoids eggs/dairy/wheat",
  smoking_status = "Smoking", pack_years_of_smoking = "Pack-years",
  current_tobacco_smoking = "Current smoking",
  alcohol_intake_frequency = "Alcohol freq.", alcohol_drinker_status = "Alcohol status",
  bread_type = "Bread type", cereal_type = "Cereal type",
  plays_computer_games = "Computer games", past_tobacco_smoking = "Past tobacco",
  snoring = "Snoring", nap_during_day = "Daytime napping",
  daytime_dozing_sleeping = "Daytime dozing", sleep_duration = "Sleep duration",
  age_first_had_sexual_intercourse = "Age at first sex",
  time_spend_outdoors_in_summer = "Time outdoors (summer)",
  time_spent_watching_television_tv = "TV time", time_spent_using_computer = "Computer time",
  index_of_multiple_deprivation_england = "Deprivation",
  index_of_multiple_deprivation = "Deprivation",
  education_score_england = "Education score", health_score_england = "Health score",
  alcohol_intake_versus_10_years_previously = "Alcohol vs 10y",
  use_of_sun_uv_protection = "Sun protection",
  number_of_days_week_walked_10_plus_minutes = "Walk days/wk")
# include the categorical level (.multi_X or trailing _Word) so distinct levels read right
prettify <- function(id) {
  lvl <- rep("", length(id))
  hm <- grepl("\\.multi_", id);                 lvl[hm] <- sub(".*\\.multi_", "", id[hm])
  hw <- !hm & grepl("_f[0-9]+_0_0_[A-Za-z]", id); lvl[hw] <- sub("^.*_f[0-9]+_0_0_", "", id[hw])
  lvl  <- gsub("[._]+", " ", lvl)
  lvl  <- sub("Fish oil.*", "Fish oil", lvl, ignore.case = TRUE)   # collapse long supplement level
  lvl  <- sub("None of the above", "none", lvl, ignore.case = TRUE)
  lvl  <- ifelse(tolower(lvl) == "yes", "", lvl)   # presence level -> show stem only
  base <- sub("\\.multi_.*$", "", id); base <- sub("_f[0-9].*$", "", base)
  stem <- ifelse(base %in% names(SHORT), SHORT[base], gsub("_", " ", heap_pretty_field(base)))
  # redundant only if every WORD of the level is already a word of the stem
  # (whole-word, not substring -- else "No" matches "s(no)ring" and gets dropped)
  red  <- mapply(function(s, l) {
            if (!nzchar(l)) return(FALSE)
            lw <- strsplit(tolower(l), "\\s+")[[1]]; sw <- strsplit(tolower(s), "\\s+")[[1]]
            all(lw %in% sw)
          }, stem, lvl)
  out  <- ifelse(nzchar(lvl) & !red, paste0(stem, ": ", lvl), stem)
  # CELL (composite) keeps tight 18-char labels; the standalone landscape has a
  # roomy left margin, so allow longer, more readable labels there.
  cap <- if (CELL) 17L else 34L
  ifelse(nchar(out) > cap, paste0(substr(out, 1, cap - 2L), "…"), out)
}

cor_dt <- .read("intervention_correlations.tsv")
for (cc in c("r","pval","n_eff","pval_BH")) cor_dt[[cc]] <- suppressWarnings(as.numeric(cor_dt[[cc]]))
cor_dt <- cor_dt[covarType == sel_covar]
# exposure -> fine category (from the scatter table)
catmap <- unique(.read("intervention_scatter.tsv")[, .(exposure_id, Category)])
cor_dt <- merge(cor_dt, catmap, by = "exposure_id", all.x = TRUE)
cor_dt[is.na(Category) | Category == "", Category := "Other"]
# MR-ANCHORING (validation): # of an exposure's concordant proteins that are GOLD
# (cis-pQTL colocalized OR Tier1) causal for cardiometabolic disease. Main figure
# keeps only MR-anchored exposures (so the concordance is disease-grounded);
# HEAP_PANELB_MR=0 -> full supplement heatmap (no MR gate).
MRGATE <- Sys.getenv("HEAP_PANELB_MR", "1") != "0"
DZCOLS <- c("T2D","Obesity","Lipids","Hypertension")            # unified cardiometabolic set
val <- .read("exposure_mr_validation.tsv")[, c("exposure_id","n_causal", DZCOLS), with = FALSE]
cor_dt <- merge(cor_dt, val, by = "exposure_id", all.x = TRUE)
for (cc in c("n_causal", DZCOLS)) cor_dt[is.na(get(cc)), (cc) := 0L]
CAT_LAB <- c(Diet_Weekly="Diet", Exercise_Freq="Exercise", Exercise_MET="Activity",
             Sleep="Sleep", Smoking="Smoking", Alcohol="Alcohol", Vitamins="Vitamins",
             Sun_Exposure="Sun", Internet_Usage="Screen", Sexual_Factors="Sexual",
             Deprivation_Indices="SES", Residential_Air_Pollution="Air poll.", Other="Other")
clab <- function(x) ifelse(x %in% names(CAT_LAB), CAT_LAB[x], gsub("_", " ", x))

# binary fields appear as mirror _Yes/_No levels (same info, opposite sign):
# keep the presence (_Yes) level, drop the redundant _No mirror.
cor_dt[, bfield := sub("(_f[0-9]+_0_0).*$", "\\1", exposure_id)]
yes_fields <- unique(cor_dt[grepl("_Yes$", exposure_id), bfield])
cor_dt <- cor_dt[!(grepl("_No$", exposure_id) & bfield %in% yes_fields)]

# --- SELECTION / DISPLAY RULES ----------------------------------------------
# (1) per cell: robust = n_eff>=MIN_NEFF; significant = robust & FDR<0.05. pval_BH
#     is NA below the n_eff floor (the FDR family was restricted to n_eff>=8 in
#     run_intervention_compare.R).
# (2) ONE ROW PER FIELD: collapse ordered-factor contrast levels to the single
#     MOST-SIGNIFICANT level; .multi_/categorical levels stay distinct.
# (3) eligibility: FDR-significant in >=1 trial.
# (4) GROUP BY exposure CATEGORY; within each category RANK BY SIGNIFICANCE
#     (smallest FDR first), keep top CAP per category.
MIN_NEFF <- 8
CAP <- as.integer(Sys.getenv("HEAP_PANELB_CATCAP", unset = "3"))
cor_dt[, robust := is.finite(n_eff) & n_eff >= MIN_NEFF]
cor_dt[, sigc   := robust & is.finite(pval_BH) & pval_BH < 0.05]
cor_dt[, gkey   := sub("(_f[0-9]+_0_0)[0-9]+$", "\\1", exposure_id)]

cscore <- cor_dt[, .(minp = suppressWarnings(min(pval_BH, na.rm = TRUE)),
                     nsig = sum(sigc, na.rm = TRUE), n_causal = n_causal[1]), by = .(gkey, exposure_id, Category)]
cscore[!is.finite(minp), minp := NA_real_]
elig <- cscore[nsig >= 1 & is.finite(minp)]                 # FDR-significant in >=1 trial
if (MRGATE) elig <- elig[n_causal >= 1]                     # AND >=1 Tier-1-or-above causal protein
best <- elig[order(gkey, minp)][, .SD[1], by = gkey]        # most-sig contrast/field
best <- best[order(Category, minp)][, crank := seq_len(.N), by = Category]
kept <- best[crank <= CAP]

cor_dt <- cor_dt[exposure_id %in% kept$exposure_id]
cor_dt[robust == FALSE, `:=`(r = NA_real_, pval_BH = NA_real_)]   # grey the sub-threshold cells
cor_dt[, exposure_lab := prettify(exposure_id)]

int_levels <- c("HERITAGE_effect","GLP1_effect1","GLP1_effect2")
int_labels <- c("HERITAGE","GLP1 STEP1","GLP1 STEP2")
cor_dt[, intervention_lab := factor(intervention, levels = int_levels, labels = int_labels)]
# GRADED significance stars (BH-FDR thresholds)
cor_dt[, star := fifelse(!is.finite(pval_BH), "",
                  fifelse(pval_BH < 0.001, "***",
                  fifelse(pval_BH < 0.01,  "**",
                  fifelse(pval_BH < 0.05,  "*", ""))))]

# category blocks ordered by their most-significant exposure; rows within a block
# ordered most-significant at TOP
kept <- merge(kept, unique(cor_dt[, .(exposure_id, exposure_lab)]), by = "exposure_id")
catord <- kept[, .(catmin = min(minp)), by = Category][order(catmin)]
kept[, Category := factor(Category, levels = catord$Category)]
labord <- kept[order(Category, -minp)]
cor_dt[, exposure_lab := factor(exposure_lab, levels = labord$exposure_lab)]
cor_dt[, Category_lab := factor(clab(Category), levels = clab(catord$Category))]
CAT_ORDER <- as.character(catord$Category)   # facet top->bottom order, for strip colours

DZLAB <- c(T2D="T2D", Obesity="Obesity", Lipids="Lipids", Hypertension="HTN")   # display
XLEV <- c(int_labels, unname(DZLAB))
plot_dt <- cor_dt[, .(exposure_id, exposure_lab, Category_lab,
                      intervention_lab = factor(as.character(intervention_lab), levels = XLEV), r, pval_BH, n_eff, star)]
# per-disease Tier-1-or-above CAUSAL counts (validation, specified per disease)
kd <- unique(cor_dt[, c("exposure_lab","Category_lab", DZCOLS), with = FALSE])
mrcols <- melt(kd, id.vars = c("exposure_lab","Category_lab"), measure.vars = DZCOLS, variable.name = "dz", value.name = "cnt")
mrcols[, intervention_lab := factor(DZLAB[as.character(dz)], levels = XLEV)]
mrcols[, cnt := as.integer(cnt)]
# The three trial columns carry a continuous r and need their width; the four
# disease columns carry a single digit, so they are set narrower. A discrete x
# gives every level the same slot, so positions are mapped numerically instead.
NXI  <- length(int_labels); WDZ <- 0.60
xpos <- c(seq_len(NXI), NXI + 0.5 + (seq_along(DZLAB) - 0.5) * WDZ)
names(xpos) <- XLEV
plot_dt[, `:=`(xn = xpos[as.character(intervention_lab)], wd = 0.98)]
mrcols[,  `:=`(xn = xpos[as.character(intervention_lab)], wd = WDZ * 0.96)]

MRBRK <- sort(unique(mrcols[cnt > 0]$cnt)); if (length(MRBRK) < 2) MRBRK <- c(1L, 2L)
rlim <- max(abs(plot_dt$r), na.rm = TRUE); rlim <- if (is.finite(rlim) && rlim > 0) rlim else 1
message(sprintf("panel b: %d fields across %d categories%s; +%d per-disease causal columns",
                nrow(kept), nrow(catord), if (MRGATE) " (MR-anchored)" else " (FULL/supplement)", length(DZCOLS)))

p <- ggplot(plot_dt, aes(xn, exposure_lab)) +
  geom_tile(aes(fill = r, width = wd), colour = "grey88", linewidth = 0.25) +
  geom_text(aes(label = star), size = STAR, vjust = 0.78, colour = "grey15") +
  geom_tile(data = mrcols, aes(width = wd), fill = "grey97", colour = "grey85", linewidth = 0.25) +
  geom_tile(data = mrcols[cnt > 0], aes(alpha = cnt, width = wd), fill = "#6A51A3", colour = "grey70", linewidth = 0.25) +
  geom_text(data = mrcols[cnt > 0], aes(label = cnt), size = if (CELL) 2.2 else 3.2, colour = "white", fontface = "bold") +
  scale_alpha_continuous(range = c(0.55, 1), breaks = MRBRK, limits = range(MRBRK),
                         name = if (CELL) "MR-causal proteins" else "Tier-1 MR-causal proteins",
                         guide = guide_legend(order = 2, direction = "horizontal",
                           override.aes = list(fill = "#6A51A3", colour = "grey70", linewidth = 0.25))) +
  scale_fill_gradient2(low = "#1B6CA8", mid = "white", high = "#B2182B", midpoint = 0,
                       limits = c(-rlim, rlim), name = "Weighted r", na.value = "grey90",
                       guide = guide_colourbar(order = 1)) +
  scale_x_continuous(breaks = unname(xpos), labels = XLEV, position = "top",
                     expand = expansion(add = c(0.03, 0.03))) +
  facet_grid(rows = vars(Category_lab), scales = "free_y", space = "free_y", switch = "y") +
  labs(title = if (CELL) "Interventional concordance" else "HEAP vs intervention proteomics",
       x = NULL, y = NULL, caption = NULL) +
  theme_heap(base_size = BS) +
  theme(axis.text.x = element_text(angle = 30, hjust = 0, vjust = 0,
                                   face = if (CELL || FULL) "plain" else "bold",
                                   size = if (CELL) 5.2 else if (FULL) 8.5 else NA,
                                   margin = margin(b = if (CELL) 1 else 4)),
        axis.text.y = element_text(size = if (CELL) 6.6 else if (FULL) 8 else 8.5,
                                   margin = margin(r = if (CELL) 1 else if (FULL) 2 else 3,
                                                   l = if (CELL) 2.5 else 3)),
        axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.spacing.y = grid::unit(if (CELL) 1.3 else if (FULL) 1.6 else 5, "pt"),
        strip.placement = "outside",
        strip.background = element_rect(fill = "grey92", colour = "grey75", linewidth = 0.3),
        strip.text.y.left = element_text(angle = 0, hjust = 0.5, face = "bold",
                                         size = if (CELL) 4.3 else if (FULL) 7.5 else 8.5, colour = "grey15",
                                         margin = margin(r = if (CELL) 2.5 else if (FULL) 2 else 4,
                                                         l = if (CELL) 1 else if (FULL) 2 else 4)),
        plot.title = element_text(size = if (CELL) TTL else if (FULL) TTL else 14, hjust = 0.5, face = "bold",
                                  margin = margin(b = if (CELL) 3 else if (FULL) 5 else 10)),
        plot.title.position = "panel",
        legend.position = if (CELL) "bottom" else "right",
        legend.key.height = grid::unit(if (CELL) 0.30 else if (FULL) 0.8 else 1.1, "lines"),
        legend.key.width  = grid::unit(if (CELL) 0.9 else 0.7, "lines"),
        legend.title = element_text(size = if (CELL) 6.8 else if (FULL) 7 else 9, vjust = 0.5),
        legend.text  = element_text(size = if (CELL) 6.3 else if (FULL) 6.5 else 8),
        legend.margin = margin(l = if (CELL) 0 else 6),
        legend.box = if (CELL) "horizontal" else "vertical",
        legend.title.position = "top",
        legend.box.spacing = grid::unit(if (CELL) 4 else 6, "pt"),
        legend.spacing.x = grid::unit(if (CELL) 10 else 6, "pt"),
        plot.caption = element_text(size = if (CELL) 5.0 else if (FULL) 6 else 8, colour = "grey35", hjust = 0,
                                    margin = margin(t = if (CELL) 2 else if (FULL) 3 else 8)),
        plot.margin = if (CELL) margin(3, 6, 2, 14) else if (FULL) margin(4, 6, 3, 4) else margin(8, 12, 8, 8))
if (CELL) p <- p + guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5))

# colour each category's facet strip with the canonical exposure-category palette
# (HEAP_ECAT_COLORS) by editing the built gtable — ggh4x/ggtext aren't installed.
strip_colored_grob <- function(plot, cat_keys) {
  g  <- ggplot2::ggplot_gtable(ggplot2::ggplot_build(plot))
  sl <- which(grepl("strip-l", g$layout$name)); sl <- sl[order(g$layout$t[sl])]
  fills <- unname(HEAP_ECAT_COLORS[cat_keys]); fills[is.na(fills)] <- "#BDBDBD"
  lum  <- colSums(grDevices::col2rgb(fills) * c(0.299, 0.587, 0.114))
  txtc <- ifelse(lum > 140, "grey10", "white")
  # strip grob nesting differs across ggplot2 versions (4.0 wraps it in a
  # strip.gTree); recurse the whole tree, recolouring every rect + text.
  recolor <- function(gr, fill, txt) {
    nm <- if (is.null(gr$name)) "" else gr$name
    if (inherits(gr, "rect") || grepl("background", nm)) {
      if (is.null(gr$gp)) gr$gp <- grid::gpar()
      gr$gp$fill <- fill; gr$gp$col <- "grey70"
    }
    if (inherits(gr, "text")) { if (is.null(gr$gp)) gr$gp <- grid::gpar(); gr$gp$col <- txt }
    if (!is.null(gr$children)) for (i in seq_along(gr$children)) gr$children[[i]] <- recolor(gr$children[[i]], fill, txt)
    if (!is.null(gr$grobs))    for (i in seq_along(gr$grobs))    gr$grobs[[i]]    <- recolor(gr$grobs[[i]], fill, txt)
    gr
  }
  for (i in seq_along(sl)) g$grobs[[sl[i]]] <- recolor(g$grobs[[sl[i]]], fills[i], txtc[i])
  g
}
emit_b <- function(png_f, pdf_f, W, H, dpi) {
  g <- tryCatch(strip_colored_grob(p, CAT_ORDER),
                error = function(e) { message("strip recolour skipped: ", conditionMessage(e)); ggplot2::ggplotGrob(p) })
  grDevices::png(png_f, width = W, height = H, units = "in", res = dpi, type = "cairo", bg = "white"); grid::grid.draw(g); grDevices::dev.off()
  grDevices::cairo_pdf(pdf_f, width = W, height = H, bg = "white"); grid::grid.draw(g); grDevices::dev.off()
}

if (CELL) {
  emit_b(file.path(FIGDIR, "fig_m4_panel_b_cell.png"), file.path(FIGDIR, "fig_m4_panel_b_cell.pdf"), B_W, B_H, 400)
  message("panel b CELL done (", nrow(kept), " fields)")
} else {
  # Standalone landscape: the FULL (no-MR-gate, no-cap) cut is the supplementary
  # figure fig_intervention_concordance_full; the MR-gated cut keeps the panel_b
  # name. Size scales with the row + category-strip count so every exposure gets
  # one readable row regardless of how many fields are BH-significant.
  n_rows  <- nrow(kept)
  n_strip <- nrow(catord)
  base    <- if (nzchar(OUTNAME)) OUTNAME
             else if (!MRGATE) "fig_intervention_concordance_full"
             else "fig_m4_panel_b"
  # COMPREHENSIVE REFERENCE figure: the full (no-MR-gate, no-cap) cut has ~70-80
  # rows. Rather than squeeze every row onto a fixed portrait page (which shrinks
  # the y-labels to an illegible smear), size the canvas to the row count
  # (~0.185 in/row) so each exposure gets one legible row when the PDF is zoomed.
  # The MR-gated panel_b cut (fewer rows) uses the same row-scaled rule.
  W <- 9.2
  H <- 1.9 + 0.185 * n_rows + 0.07 * n_strip
  H <- max(6.0, min(H, 18.0))
  emit_b(file.path(FIGDIR, paste0(base, ".png")), file.path(FIGDIR, paste0(base, ".pdf")), W, H, 200)
  message(sprintf("panel b standalone done -> %s (%d fields, %d categories, %.1f x %.1f in)",
                  base, n_rows, n_strip, W, H))
}
