#!/usr/bin/env Rscript

# ============================================================================
# fig_mediation_prioritization.R  [figure_id: fig_mediation_prioritization]
# ----------------------------------------------------------------------------
# CONSOLIDATED (2026-07-11). Merges the three Module-3 figures cited by ONE clause
# ("these proteins provided a prioritized set of exposure-responsive biology
# linked to modifiable disease risk"):
#   fig_mediation_decomp             direct vs through-protein decomposition
#   fig_mediation_targets            curated protein intermediaries
#   fig_mediation_exposure_anchoring mediation tracks Module-1 exposure R2
#
# Read together they say: the mediated share is small (a), yet the proteins it
# nominates are real, powered and disease-specific (b), and they are the SAME
# proteins Module 1 already flagged as exposure-responsive (c) -- so the
# prioritization is anchored in exposure biology rather than in noise.
#
#   a  decomposition  for the strongest lifestyle->protein->disease links, the
#                     total effect split into the part running through the protein
#                     (NIE) and the direct remainder (NDE), with % mediated
#   b  the shortlist  curated intermediaries for the diseases that actually carry
#                     strong mediation: mediated HR with 95% CI, colored by driver
#   c  the anchor     per exposure category, Spearman rho between a protein's
#                     Module-1 exposure R2 and the strength of its mediated effect.
#                     The major lifestyle domains anchor strongly; the categories
#                     with little proteomic reach (vitamins, sun) do not -- which
#                     is what an honest anchor looks like.
#
# Panel c replaces the 13-facet scatter of fig_mediation_exposure_anchoring: at
# 6.5in those facets were unreadable, and the per-category rho IS the claim.
#
# Input : module3 partitioned_categories + exploratory/module3/mediation_x_variance.tsv
# Output: figures/supplement/module3/fig_mediation_prioritization.{pdf,png} + data tsv
#
# Authored at 6.5in = the supplement's \textwidth (scale 1.0).
# ============================================================================

local({
  cand <- c(file.path(getwd(), "scripts", "visualizations", "common"),
            "/n/groups/patel/shakson_ukb/HEAP/scripts/visualizations/common")
  common <- cand[dir.exists(cand)][1]
  if (is.na(common)) stop("cannot locate common/ helpers")
  for (f in c("figure_paths", "load_heap_results", "plot_theme",
              "label_helpers", "export_helpers"))
    source(file.path(common, paste0(f, ".R")))
})
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork)
})

a <- commandArgs(trailingOnly = TRUE); a <- a[!startsWith(a, "--")]
a <- a[!a %in% c("fig_mediation_prioritization", "all_main", "all_supplement", "all", "website")]
covarType <- if (length(a) >= 1) a[1] else "base"
family    <- if (length(a) >= 2) a[2] else "lasso"
figure_id <- Sys.getenv("HEAP_FIGURE_ID", unset = "fig_mediation_prioritization")
BS <- 7.5
ALPHA <- 0.05

sel <- c("protID", "DZ_ID", "predictor", "predictor_class", "effect_type", "effect_logHR",
         "effect_HR", "delta_p", "delta_l95_HR", "delta_u95_HR", "instrument_present",
         "n", "n_cases", "protein_HR", "protein_p")
md <- load_module3_results(covarType = covarType, family = family,
                           mode = "partitioned_categories", select = sel)
md <- md[!(predictor_class %in% c("genetic_cis", "genetic_trans") & instrument_present == FALSE)]
pm <- heap_proportion_mediated(md)
pm[, category := heap_md_category(predictor)]
pm <- pm[!is.na(category) & is.finite(NIE_q) & NIE_q < ALPHA & pm_consistent == TRUE]
pm[, absNIE := abs(NIE_logHR)]

DRV <- c(Exercise_Freq = "Physical activity", Exercise_MET = "Physical activity",
         Diet_Weekly = "Diet", Alcohol = "Alcohol", Smoking = "Smoking", Sleep = "Sleep",
         Sun_Exposure = "Sun", Internet_Usage = "Internet", Vitamins = "Vitamins",
         Deprivation_Indices = "Deprivation", Sexual_Factors = "Sexual")
PAL <- c(`Physical activity` = "#1B7837", Diet = "#7FBF7B", Alcohol = "#E07B39",
         Smoking = "#525252", Sleep = "#6BAED6", Sun = "#FDB863", Internet = "#C51B7D",
         Vitamins = "#B8860B", Deprivation = "#80CDC1", Sexual = "#C994C7")
pm[, Lifestyle := fifelse(category %in% names(DRV), DRV[category], category)]
pm[, disease := heap_pretty_disease(DZ_ID)]
DZS <- c(`non insulin dependent diabetes mellitus` = "T2D", `chronic renal failure` = "CKD",
         obesity = "Obesity", `other diseases of liver` = "Liver disease",
         `heart failure` = "Heart failure", `acute renal failure` = "AKI")
pm[, dzL := fifelse(disease %in% names(DZS), DZS[disease], disease)]

# ---- a: direct vs through-protein --------------------------------------------
K <- 10L; NMIN_D <- 300L
d <- pm[is.finite(pm_display) & n_cases >= NMIN_D][order(-absNIE)][, .SD[1], by = .(protID, DZ_ID)]
d <- d[order(-absNIE)][seq_len(min(K, .N))]
d[, lab := sprintf("%s -> %s", protID, dzL)]
d[, lab := factor(lab, levels = rev(lab))]
dl <- rbind(
  d[, .(lab, logHR = NDE_logHR, fillv = "Direct (not via the protein)")],
  d[, .(lab, logHR = NIE_logHR, fillv = Lifestyle)])
message(sprintf("decomp: %d links | median %% mediated %.0f%%", nrow(d), 100 * median(d$pm_display)))

pa <- ggplot(dl, aes(logHR, lab, fill = fillv)) +
  geom_col(width = .66) +
  geom_text(data = d, aes(x = TE_logHR, y = lab,
                          label = sprintf("%d%%", round(100 * pm_display))),
            inherit.aes = FALSE, hjust = -0.25, size = 1.8, colour = "grey25") +
  scale_fill_manual(values = c(PAL, `Direct (not via the protein)` = "grey82"), name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0, .16))) +
  labs(x = "log HR of the exposure's total effect on disease", y = NULL) +
  theme_heap(base_size = BS) +
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(size = BS - 2),
        legend.position = "bottom", legend.key.size = unit(6, "pt"),
        legend.text = element_text(size = BS - 2.5),
        legend.margin = margin(0, 0, 0, 0),
        plot.margin = margin(10, 4, 2, 2)) +
  guides(fill = guide_legend(nrow = 2))

# ---- b: the curated shortlist -------------------------------------------------
MEFF <- 1.10; NMIN_T <- 500L; NDIS <- 3L; KPROT <- 6L
strong <- pm[absNIE >= log(MEFF) & n_cases >= NMIN_T]
keep_dz <- head(strong[, .(N = uniqueN(protID)), by = DZ_ID][order(-N)]$DZ_ID, NDIS)
strong <- strong[DZ_ID %in% keep_dz]
best <- strong[order(-absNIE)][, .SD[1], by = .(protID, DZ_ID)]
best <- best[order(DZ_ID, -absNIE)][, head(.SD, KPROT), by = DZ_ID]
nie <- md[effect_type == "NIE", .(protID, DZ_ID, predictor, nie_HR = effect_HR,
                                  lo = delta_l95_HR, hi = delta_u95_HR)]
best <- merge(best, nie, by = c("protID", "DZ_ID", "predictor"), all.x = TRUE)
best[, nie_HR := fifelse(is.finite(nie_HR), nie_HR, NIE_HR)]
shortdz <- function(x) { y <- heap_pretty_disease(x); fifelse(y %in% names(DZS), DZS[y], y) }
best[, dz := factor(shortdz(DZ_ID), levels = unique(shortdz(keep_dz)))]
best <- best[order(dz, nie_HR)]
best[, key := factor(paste(protID, DZ_ID, sep = "__"), levels = paste(protID, DZ_ID, sep = "__"))]
message(sprintf("targets: %d rows across %s", nrow(best), paste(levels(best$dz), collapse = " | ")))

pb <- ggplot(best, aes(nie_HR, key, colour = Lifestyle)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey55", linewidth = .3) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, linewidth = .4, alpha = .85) +
  geom_point(size = 1.2) +
  facet_grid(dz ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_colour_manual(values = PAL, name = NULL, guide = "none") +
  scale_y_discrete(labels = function(z) sub("__.*$", "", z)) +
  labs(x = "Mediated hazard ratio (NIE)", y = NULL) +
  theme_heap(base_size = BS) +
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(size = BS - 2.5),
        strip.text.y.left = element_text(angle = 0, size = BS - 2, face = "bold"),
        panel.spacing.y = unit(3, "pt"),
        strip.placement = "outside",
        plot.margin = margin(10, 4, 2, 2))

# ---- c: is the prioritization anchored in exposure-responsiveness? ------------
ANCH <- file.path(Sys.getenv("HEAP_PROJECT_ROOT", "/n/groups/patel/IGLOO/UKB/HEAP"), "figures", "exploratory/module3/mediation_x_variance.tsv")
j <- fread(ANCH)
j <- j[is.finite(m1_r2) & is.finite(NIE_HR) & is.finite(absNIE)]
j[, m1_r2 := pmax(m1_r2, 0)]
rc <- j[, .(rho = cor(absNIE, m1_r2, method = "spearman"), n = .N), by = category][n >= 50][order(rho)]
rc[, catf := factor(heap_category_pretty(category),
                    levels = heap_category_pretty(category))]
message(sprintf("anchoring: %d categories | rho %.2f (%s) to %.2f (%s)",
                nrow(rc), rc$rho[1], rc$catf[1], rc$rho[nrow(rc)], rc$catf[nrow(rc)]))

pc <- ggplot(rc, aes(rho, catf, fill = rho > 0)) +
  geom_vline(xintercept = 0, colour = "grey55", linewidth = .3) +
  geom_col(width = .66) +
  geom_text(aes(label = sprintf("%.2f  (n=%s)", rho, formatC(n, big.mark = ",", format = "d")),
                hjust = fifelse(rho > 0, -0.1, 1.1)),
            size = 1.8, colour = "grey25") +
  scale_fill_manual(values = c(`TRUE` = "#1B7837", `FALSE` = "#B0B0B0"), guide = "none") +
  scale_x_continuous(limits = c(-0.45, 1.05), expand = expansion(mult = c(.02, .02))) +
  labs(x = "Spearman rho:  Module-1 exposure R²  vs  strength of mediated effect",
       y = NULL) +
  theme_heap(base_size = BS) +
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(size = BS - 2),
        plot.margin = margin(10, 4, 2, 2))

p <- ((pa | pb) / pc) +
  plot_layout(heights = c(1.25, 1)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 9),
        plot.tag.position = c(0, 1))

out <- rbindlist(list(
  d[,    .(panel = "a", key = as.character(lab), value = pm_display)],
  best[, .(panel = "b", key = paste(protID, dz, sep = " | "), value = nie_HR)],
  rc[,   .(panel = "c", key = as.character(catf), value = rho)]), use.names = TRUE)

heap_emit_figure(p, figure_id, data = out, category = "supplement", subdir = "module3",
                 formats = c("pdf", "png"), width = 6.5, height = 6.4, website = TRUE)

message("fig_mediation_prioritization: done.")
