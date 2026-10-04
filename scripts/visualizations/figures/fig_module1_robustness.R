#!/usr/bin/env Rscript

# ============================================================================
# fig_module1_robustness.R  [figure_id: fig_module1_robustness]
# ----------------------------------------------------------------------------
# CONSOLIDATED Module-1 robustness figure. Replaces four figures that were all
# the same thing -- "the estimate barely moves when I perturb X", drawn as an
# agreement scatter with a correlation annotation -- differing only in X:
#   fig_traintest_stability          X = the train/test split (whole model)
#   fig_traintest_stability_categories X = the split, per exposure category (13 facets)
#   fig_spec_sensitivity             X = the covariate specification
#   fig_method_robustness            X = the penalty (lasso vs elastic-net)
#
#   a  Every perturbation on one axis. A forest of concordance coefficients
#      makes "robust to everything except BMI/clinical" a single glance instead
#      of three figures.
#   b  The one perturbation that DOES move the answer: adjusting for BMI (and
#      clinical covariates) roughly halves the exposure-responsive set and does
#      it specifically to the adiposity proteins (LEP, FABP4, IL1RN...). This is
#      a result, not a QC check.
#   c  Per-category generalization TRACKS REACH (Spearman 0.88): the categories
#      that fail to replicate train->test are exactly the ones that reach almost
#      no proteins (residential noise reaches 0 proteins and correlates at
#      -0.11 -- exactly what pure noise should do). The old 13-facet grid showed
#      13 scatters and buried this; one panel says it. The grid itself remains
#      available (cited from Methods) for readers who want the per-category detail.
#
# Input : module1_predictive_r2_score_partition/<experiment>/<covarType>/<method>/
#           predictive_r2_coarse_*  (needs r2_train)  +  predictive_r2_exposure_categories_*
# Output: figures/supplement/module1/fig_module1_robustness.{pdf,png} + data tsv
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
  library(data.table); library(ggplot2); library(patchwork); library(ggrepel)
})

figure_id <- Sys.getenv("HEAP_FIGURE_ID", unset = "fig_module1_robustness")
BS  <- 8
THR <- 0.01

PC  <- HEAP_PAL_COMPONENT
COL_OK  <- "#2E9E48"   # perturbation leaves the answer intact
COL_ATT <- "#D55E00"   # perturbation attenuates it
COL_BASE<- "#404040"

# ============================== data ========================================
# ---- covariate specifications (exposomic unique R2) -------------------------
# Prevalent disease appears as an EXCLUSION only. The +prevalent covariate run was
# computed but is not shown or released: conditioning on a variable downstream of
# exposure opens the collider path the reviews raised, and in the mediation deposit it
# collapsed 22,270 significant links to 1 -- the signature of adjusting away the
# estimand. `excl. prevalent` answers the same question without that risk.
SPECS <- list(
  base      = list(exp = "M1_base_lasso",           cov = "base",           lab = "base",            grp = "base"),
  draw      = list(exp = "M1_base_draw_lasso",      cov = "base_draw",      lab = "+ draw",          grp = "ok"),
  exclprev  = list(exp = "M1_base_exclprev_lasso",  cov = "base",           lab = "excl. prevalent", grp = "ok"),
  bmi       = list(exp = "M1_base_bmi_lasso",       cov = "base_bmi",       lab = "+ BMI",           grp = "att"),
  clinical  = list(exp = "M1_base_clinical_lasso",  cov = "base_clinical",  lab = "+ clinical",      grp = "att"))

getE <- function(s) {
  co <- load_module1_predictive_r2(s$cov, "lasso", "coarse", experiment = s$exp)
  co[get("method") == "score_unique_drop" & block == "E", .(E = pmax(0, mean(r2))), by = omic]
}
E <- lapply(SPECS, getE)
W <- Reduce(function(a, b) merge(a, b, by = "omic"),
            lapply(names(E), function(n) setnames(copy(E[[n]]), "E", n)))

spec_stat <- rbindlist(lapply(names(SPECS), function(n) data.table(
  block = "Covariate specification",
  lab   = SPECS[[n]]$lab,
  grp   = SPECS[[n]]$grp,
  n     = sum(E[[n]]$E >= THR),
  rho   = cor(W$base, W[[n]], method = "spearman"))))

# ---- penalty (lasso vs elastic-net), per component --------------------------
COMPS <- c(Exposomic = "E", Genetic = "G", GxE = "GxE")
D_l <- load_module1_predictive_r2("base", "lasso", "coarse", experiment = "M1_base_lasso")
D_e <- load_module1_predictive_r2("base", "enet",  "coarse", experiment = "M1_base_enet")
uniq <- function(D, blk) D[get("method") == "score_unique_drop" & block == blk,
                           .(val = mean(r2)), by = omic]
pen_stat <- rbindlist(lapply(names(COMPS), function(cn) {
  m <- merge(uniq(D_l, COMPS[[cn]]), uniq(D_e, COMPS[[cn]]), by = "omic",
             suffixes = c("_lasso", "_enet"))
  data.table(block = "Penalty (lasso → elastic-net)", lab = cn, grp = "ok",
             n = NA_integer_, rho = cor(m$val_lasso, m$val_enet, method = "spearman"))
}))

# ---- held-out generalization (whole model + per component) ------------------
pr <- load_module1_predictive_r2("base", "lasso", "coarse", experiment = "M1_base_lasso")
tot <- pr[get("method") == "score_model_total" & block == "C+G+E+GxE" &
            is.finite(r2) & is.finite(r2_train),
          .(train = mean(r2_train), test = mean(r2)), by = omic]
gen_stat <- data.table(block = "Held-out generalization", lab = "whole model", grp = "ok",
                       n = NA_integer_,
                       rho = cor(tot$train, tot$test, method = "spearman"))

fore <- rbindlist(list(gen_stat, pen_stat, spec_stat))
# short, wrapped strip labels: a long strip label makes facet_grid reserve a very
# wide strip column and starves the actual panel
fore[block == "Held-out generalization",       block := "Held-out\ngeneralization"]
fore[block == "Penalty (lasso → elastic-net)", block := "Penalty\n(lasso vs enet)"]
fore[block == "Covariate specification",       block := "Covariate\nspecification"]
BLKS <- c("Held-out\ngeneralization", "Penalty\n(lasso vs enet)", "Covariate\nspecification")
fore[, block := factor(block, levels = BLKS)]
setorder(fore, block, rho)
fore[, ylab := factor(seq_len(.N), labels = lab)]
fore[, ylab := factor(lab, levels = lab)]

# ---- per-category generalization vs reach ----------------------------------
ec <- load_module1_predictive_r2("base", "lasso", "exposure_categories",
                                 experiment = "M1_base_lasso")
u  <- ec[get("method") == "score_unique_drop" & is.finite(r2) & is.finite(r2_train)]
ca <- u[, .(train = mean(r2_train), test = mean(r2)), by = .(omic, block)]
cat_stat <- ca[, .(r = cor(train, test, use = "complete.obs"),
                   reach = sum(test >= 0.001, na.rm = TRUE)), by = block]
cat_stat[, category := heap_americanize(gsub("_", " ", sub("^PXS_", "", block)))]
RHO_CR <- cor(cat_stat$r, cat_stat$reach, method = "spearman")

# ============================== panels ======================================
# ---- a: concordance forest --------------------------------------------------
GC <- c(base = COL_BASE, ok = COL_OK, att = COL_ATT)
pa <- ggplot(fore, aes(rho, ylab, colour = grp)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey70", linewidth = .3) +
  geom_segment(aes(x = 0.8, xend = rho, y = ylab, yend = ylab), linewidth = .3, alpha = .5) +
  geom_point(size = 2) +
  geom_text(aes(label = sprintf("%.2f", rho)), hjust = -0.45, size = 2.2, colour = "grey20") +
  facet_grid(block ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_colour_manual(values = GC, guide = "none") +
  scale_x_continuous(limits = c(0.8, 1.045), breaks = c(0.8, 0.9, 1.0)) +
  labs(x = "Concordance with the base model (Spearman)", y = NULL) +
  theme_heap(base_size = BS) +
  theme(panel.grid = element_blank(),
        strip.placement = "outside",
        strip.background = element_blank(),
        strip.text.y.left = element_text(angle = 0, face = "bold", size = BS - 1.5,
                                         hjust = 1, lineheight = 1),
        plot.margin = margin(10, 6, 2, 4))

# ---- b: BMI attenuation -----------------------------------------------------
ADI <- c("LEP", "FABP4", "IL1RN", "CFH", "HGF", "IGSF9", "CA14", "MAMDC4")
bd  <- data.table(omic = W$omic, base = W$base, clinical = W$clinical)
bd[, flag := omic %in% ADI]
n_base <- spec_stat[lab == "base"]$n
n_clin <- spec_stat[lab == "+ clinical"]$n

pb <- ggplot(bd, aes(base, clinical)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey60", linewidth = .3) +
  geom_point(data = bd[!(flag)], colour = COL_ATT, size = .5, alpha = .3) +
  geom_point(data = bd[(flag)], colour = "grey10", size = 1.3) +
  geom_text_repel(data = bd[(flag)], aes(label = omic), size = 2.1, colour = "grey10",
                  segment.size = .2, min.segment.length = 0, max.overlaps = Inf,
                  box.padding = .3, seed = 1) +
  labs(x = expression("base exposomic"~R^2),
       y = expression("+ clinical exposomic"~R^2)) +
  theme_heap(base_size = BS) +
  theme(panel.grid = element_blank(), plot.margin = margin(10, 4, 2, 4))

# ---- c: generalization tracks reach -----------------------------------------
pc <- ggplot(cat_stat, aes(reach + 1, r)) +
  geom_hline(yintercept = 0, colour = "grey75", linewidth = .3) +
  geom_point(aes(colour = r), size = 2) +
  geom_text_repel(aes(label = category), size = 2, colour = "grey25",
                  segment.size = .2, min.segment.length = 0, max.overlaps = Inf,
                  box.padding = .32, seed = 2) +
  annotate("text", x = 1.4, y = 1.02, hjust = 0, size = 2.3, colour = "grey20",
           label = sprintf("Spearman = %.2f", RHO_CR)) +
  scale_x_log10(breaks = c(1, 11, 101, 1001), labels = c("0", "10", "100", "1000"),
                expand = expansion(mult = c(0.10, 0.16))) +   # room for right-hand labels
  scale_colour_gradient(low = COL_ATT, high = COL_OK, guide = "none") +
  scale_y_continuous(limits = c(-0.25, 1.08)) +
  labs(x = "# proteins reached",
       y = "train vs test correlation (r)") +
  theme_heap(base_size = BS) +
  theme(panel.grid = element_blank(), plot.margin = margin(10, 4, 2, 4))

# ============================== assemble ====================================
p <- (pa / (pb | pc)) +
  plot_layout(heights = c(1, 1.05)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 9),
        plot.tag.position = c(0, 1))

out <- rbindlist(list(
  fore[, .(panel = "a", key = paste(block, lab, sep = " | "), value1 = rho,
           value2 = as.numeric(n))],
  bd[, .(panel = "b", key = omic, value1 = base, value2 = clinical)],
  cat_stat[, .(panel = "c", key = category, value1 = r, value2 = as.numeric(reach))]),
  use.names = TRUE)

heap_emit_figure(p, figure_id, data = out, category = "supplement", subdir = "module1",
                 formats = c("pdf", "png"), width = 6.5, height = 6.6, website = TRUE)

message("fig_module1_robustness: done.")
print(fore[, .(block, lab, rho = round(rho, 3), n)])
cat(sprintf("\ncategory generalization vs reach: Spearman = %.3f\n", RHO_CR))
print(cat_stat[order(r), .(category, r = round(r, 2), reach)])
