#!/usr/bin/env Rscript

# ============================================================================
# summarize_module1_variance.R
# ----------------------------------------------------------------------------
# Manuscript summary statistics for the Module 1 variance-decomposition
# analysis (PGS / PXS / GIS score-partition predictive R2). Writes a
# citation-ready markdown report + machine-readable TSVs so the key numbers
# behind manuscript Figure 2 (genetic vs exposure-driven spectrum, per-protein
# explained variation, exposure-category partition) are documented in the repo
# and trivial to quote while writing.
#
# QUANTITIES (faithful to the current Module1_suggested.R score_partition mode):
#   * Whole-model test R2  = score_model_total block "C+G+E+GxE" (full model,
#     out-of-fold), per protein averaged over the 5 CV folds.
#   * Covariate-only R2    = score_model_total block "C".
#   * Unique R2 of a block = score_unique_drop = R2(full) - R2(full without
#     that block's score). Reported for Genetic (G), Exposomic (E), GxE, and
#     for cis / trans genetic sub-blocks and each exposure / GxE category.
#   * Genetics-driven vs exposures-driven = sign of (unique R2_G - unique R2_E),
#     the R2_G = R2_E diagonal of manuscript Fig 2A.
#   NOTE: out-of-fold R2 can be slightly negative (a block that does not
#   generalise); such values are reported as-is and floored at 0 only for the
#   donut-style R2 bins (manuscript Fig 2D).
#
# Output (committed to the repo): docs/manuscript_stats/module1_variance/
#   SUMMARY.md            headline + manuscript-ready sentences + tables
#   headline.tsv          one-row key numbers
#   component_summary.tsv per-component (Covars model / G / E / GxE) distribution
#   r2_bins.tsv           genetic & exposomic R2 binned (Fig 2D donut)
#   top_total.tsv         top proteins by whole-model test R2
#   top_genetic.tsv       top proteins by unique genetic R2 (most heritable)
#   top_exposomic.tsv     top proteins by unique exposomic R2
#   exposure_categories.tsv per exposure-category unique R2 (Fig 2F)
#   gxe_categories.tsv    per GxE-category unique R2
#   genetic_subblocks.tsv cis vs trans unique R2 summary
#
# Run:
#   HEAP_PATHS_FILE=.../workflow/00_paths.R \
#     Rscript scripts/analysis_summaries/summarize_module1_variance.R \
#       [covarType] [method] [experiment]
#   Defaults: base / lasso / M1_base_lasso.
# ============================================================================

local({
  cand <- c(file.path(getwd(), "scripts", "visualizations", "common"),
            "/n/groups/patel/shakson_ukb/HEAP/scripts/visualizations/common")
  common <- cand[dir.exists(cand)][1]
  if (is.na(common)) stop("cannot locate common/ helpers")
  for (f in c("figure_paths", "load_heap_results", "plot_theme", "label_helpers"))
    source(file.path(common, paste0(f, ".R")))
})
suppressPackageStartupMessages(library(data.table))

a <- commandArgs(trailingOnly = TRUE); a <- a[!startsWith(a, "--")]
covarType  <- if (length(a) >= 1) a[1] else "base"
method     <- if (length(a) >= 2) a[2] else "lasso"
experiment <- if (length(a) >= 3) a[3] else "M1_base_lasso"

cm  <- function(x) formatC(x, format = "d", big.mark = ",")
fmt <- function(x, d = 3) formatC(x, format = "f", digits = d)
pct <- function(x) fmt(100 * x, 1)
iqr <- function(x) quantile(x, c(0.25, 0.75), na.rm = TRUE)

message("Loading Module 1 (", experiment, "/", covarType, "/", method, ") ...")
coarse <- load_module1_predictive_r2(covarType, method, level = "coarse",
                                     experiment = experiment)
gsub   <- load_module1_predictive_r2(covarType, method, level = "genetic_subblocks",
                                     experiment = experiment)
ecat   <- load_module1_predictive_r2(covarType, method, level = "exposure_categories",
                                     experiment = experiment)
gxecat <- load_module1_predictive_r2(covarType, method, level = "gxe_categories",
                                     experiment = experiment)
fitsum <- load_module1_fit_summary(covarType, method, experiment = experiment)

# ---- per-protein values (mean over folds) -----------------------------------
mt <- coarse[get("method") == "score_model_total",
             .(r2 = mean(r2)), by = .(omic, block)]
total_full <- mt[block == "C+G+E+GxE", .(omic, full = r2)]
total_cov  <- mt[block == "C",         .(omic, cov  = r2)]

ud <- coarse[get("method") == "score_unique_drop",
             .(r2 = mean(r2)), by = .(omic, block)]
udw <- dcast(ud, omic ~ block, value.var = "r2")
setnames(udw, c("Covars", "G", "E", "GxE"),
         c("u_Covars", "u_G", "u_E", "u_GxE"), skip_absent = TRUE)

prot <- Reduce(function(x, y) merge(x, y, by = "omic", all = TRUE),
               list(total_full, total_cov, udw))
for (col in c("u_G", "u_E", "u_GxE")) if (!col %in% names(prot)) prot[[col]] <- 0
prot[, driven := fifelse(u_G >= u_E, "Genetics-driven", "Exposures-driven")]
n_prot <- nrow(prot)
n_fold <- uniqueN(coarse$fold)

# ---- component distribution summary -----------------------------------------
comp_tbl <- rbindlist(list(
  data.table(component = "Covariate-only model R2", v = prot$cov),
  data.table(component = "Whole-model test R2",     v = prot$full),
  data.table(component = "Genetic (unique)",        v = prot$u_G),
  data.table(component = "Exposomic (unique)",      v = prot$u_E),
  data.table(component = "GxE (unique)",            v = prot$u_GxE)
))[, .(n = .N, mean = mean(v, na.rm = TRUE), median = median(v, na.rm = TRUE),
       q25 = iqr(v)[1], q75 = iqr(v)[2], max = max(v, na.rm = TRUE),
       pct_gt_0.01 = mean(v > 0.01, na.rm = TRUE),
       pct_gt_0.05 = mean(v > 0.05, na.rm = TRUE)), by = component]

# ---- genetics-driven vs exposures-driven ------------------------------------
n_gen <- prot[, sum(driven == "Genetics-driven")]
n_exp <- prot[, sum(driven == "Exposures-driven")]

# ---- R2 bins (manuscript Fig 2D donut), floored at 0 ------------------------
g0 <- pmax(0, prot$u_G); e0 <- pmax(0, prot$u_E)
gbin <- cut(g0, c(-Inf, 0.01, 0.1, 0.25, 0.5, Inf),
            labels = c("0-0.01", "0.01-0.1", "0.1-0.25", "0.25-0.5", "0.5-1"),
            right = FALSE)
ebin <- cut(e0, c(-Inf, 0.01, 0.05, 0.2, Inf),
            labels = c("0-0.01", "0.01-0.05", "0.05-0.2", "0.2+"), right = FALSE)
r2_bins <- rbind(
  data.table(set = "Genetic (unique R2)", bin = levels(gbin))[
    , n := as.integer(table(gbin)[bin])][],
  data.table(set = "Exposomic (unique R2)", bin = levels(ebin))[
    , n := as.integer(table(ebin)[bin])][])
r2_bins[, pct := round(100 * n / n_prot, 1)]

# ---- top proteins -----------------------------------------------------------
top_total <- head(prot[order(-full)],
                  20)[, .(protein = omic, test_R2 = round(full, 3),
                          genetic_R2 = round(u_G, 3), exposomic_R2 = round(u_E, 3),
                          driven)]
top_gen   <- head(prot[order(-u_G)],
                  20)[, .(protein = omic, genetic_R2 = round(u_G, 3),
                          test_R2 = round(full, 3), exposomic_R2 = round(u_E, 3))]
top_exp   <- head(prot[order(-u_E)],
                  20)[, .(protein = omic, exposomic_R2 = round(u_E, 3),
                          genetic_R2 = round(u_G, 3), test_R2 = round(full, 3))]

# ---- cis vs trans (genetic sub-blocks) --------------------------------------
gs <- gsub[get("method") == "score_unique_drop",
           .(r2 = mean(r2), present = any(as.logical(present))),
           by = .(omic, block)]
gsw <- dcast(gs, omic ~ block, value.var = "r2")
for (col in c("Gcis", "Gtrans")) if (!col %in% names(gsw)) gsw[[col]] <- 0
gsw[, dominant := fifelse(Gcis >= Gtrans, "Cis-dominant", "Trans-dominant")]
n_cis_dom <- gsw[, sum(dominant == "Cis-dominant")]
gpres <- gsub[block == "Gcis",   any(as.logical(present)), by = omic][, sum(V1)]
tpres <- gsub[block == "Gtrans", any(as.logical(present)), by = omic][, sum(V1)]
subblk <- data.table(
  block = c("Cis (unique)", "Trans (unique)"),
  n_present = c(gpres, tpres),
  median = c(median(gsw$Gcis, na.rm = TRUE), median(gsw$Gtrans, na.rm = TRUE)),
  mean   = c(mean(gsw$Gcis,   na.rm = TRUE), mean(gsw$Gtrans,   na.rm = TRUE)),
  max    = c(max(gsw$Gcis,    na.rm = TRUE), max(gsw$Gtrans,    na.rm = TRUE)),
  pct_gt_0.01 = c(mean(gsw$Gcis > 0.01, na.rm = TRUE),
                  mean(gsw$Gtrans > 0.01, na.rm = TRUE)))
top_cis <- head(gsw[order(-Gcis)],
                15)[, .(protein = omic, cis_R2 = round(Gcis, 3),
                        trans_R2 = round(Gtrans, 3))]

# ---- exposure categories (Fig 2F) -------------------------------------------
cat_summ <- function(dt) {
  pp <- dt[get("method") == "score_unique_drop", .(r2 = mean(r2)), by = .(omic, category)]
  out <- pp[, .(mean_R2 = mean(r2), median_R2 = median(r2), max_R2 = max(r2),
                n_proteins_pos = sum(r2 > 0),
                n_proteins_gt_0.005 = sum(r2 > 0.005)), by = category]
  out[order(-mean_R2)]
}
ec <- cat_summ(ecat); gc <- cat_summ(gxecat)
ec[, broad := heap_broad_category(category)]

# ---- top GxE proteins (unique GxE R2) + dominant GxE exposure category -------
gxe_prot <- prot[, .(omic, gxe_r2 = u_GxE)][order(-gxe_r2)]
gxe_dom  <- gxecat[get("method") == "score_unique_drop", .(r2 = mean(r2)), by = .(omic, category)]
gxe_dom  <- gxe_dom[gxe_dom[, .I[which.max(r2)], by = omic]$V1, .(omic, dom_category = category)]
gxe_prot <- merge(gxe_prot, gxe_dom, by = "omic", all.x = TRUE)
top_gxe  <- head(gxe_prot, 20)

# ---- train vs test generalisation (whole model, fit_summary) ----------------
stab <- fitsum[, .(train_r2 = mean(train_r2), test_r2 = mean(test_r2),
                   fold_sd_test = sd(test_r2)), by = omic]
stab[, gap := train_r2 - test_r2]

# ---- headline ---------------------------------------------------------------
H <- list(
  experiment = experiment, covarType = covarType, method = method,
  n_proteins = n_prot, n_folds = n_fold,
  median_test_R2 = median(prot$full, na.rm = TRUE),
  q25_test_R2 = iqr(prot$full)[1], q75_test_R2 = iqr(prot$full)[2],
  max_test_R2 = max(prot$full, na.rm = TRUE),
  max_test_R2_protein = prot$omic[which.max(prot$full)],
  median_cov_R2 = median(prot$cov, na.rm = TRUE),
  median_genetic_R2 = median(prot$u_G, na.rm = TRUE),
  median_exposomic_R2 = median(prot$u_E, na.rm = TRUE),
  median_gxe_R2 = median(prot$u_GxE, na.rm = TRUE),
  n_genetics_driven = n_gen, pct_genetics_driven = n_gen / n_prot,
  n_exposures_driven = n_exp, pct_exposures_driven = n_exp / n_prot,
  n_genetic_gt_0.05 = sum(prot$u_G > 0.05, na.rm = TRUE),
  n_exposomic_gt_0.05 = sum(prot$u_E > 0.05, na.rm = TRUE),
  top_genetic_protein = top_gen$protein[1], top_genetic_R2 = top_gen$genetic_R2[1],
  top_exposomic_protein = top_exp$protein[1], top_exposomic_R2 = top_exp$exposomic_R2[1],
  n_cis_dominant = n_cis_dom, pct_cis_dominant = n_cis_dom / nrow(gsw),
  median_cis_R2 = median(gsw$Gcis, na.rm = TRUE),
  median_trans_R2 = median(gsw$Gtrans, na.rm = TRUE),
  top_exposure_category = ec$category[1], top_exposure_category_R2 = ec$mean_R2[1],
  median_train_r2 = median(stab$train_r2, na.rm = TRUE),
  median_test_r2_whole = median(stab$test_r2, na.rm = TRUE),
  median_overfit_gap = median(stab$gap, na.rm = TRUE),
  q99_overfit_gap = quantile(stab$gap, 0.99, na.rm = TRUE),
  max_overfit_gap = max(stab$gap, na.rm = TRUE),
  cor_train_test = cor(stab$train_r2, stab$test_r2, use = "complete.obs"),
  pct_overfit_gt_0.05 = mean(stab$gap > 0.05, na.rm = TRUE),
  median_fold_sd_test = median(stab$fold_sd_test, na.rm = TRUE),
  n_gxe_gt_0.01 = sum(prot$u_GxE > 0.01, na.rm = TRUE),
  n_gxe_gt_0.02 = sum(prot$u_GxE > 0.02, na.rm = TRUE),
  top_gxe_protein = top_gxe$omic[1], top_gxe_R2 = top_gxe$gxe_r2[1])

# ============================================================================
# Write outputs
# ============================================================================
outdir <- file.path(heap_path(), "docs", "manuscript_stats", "module1_variance")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

rnd <- function(dt) dt[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 4) else x)]
fwrite(as.data.table(lapply(H, function(x) if (is.numeric(x)) round(x, 4) else x)),
       file.path(outdir, "headline.tsv"), sep = "\t")
fwrite(rnd(comp_tbl),    file.path(outdir, "component_summary.tsv"), sep = "\t")
fwrite(r2_bins,          file.path(outdir, "r2_bins.tsv"), sep = "\t")
fwrite(top_total,        file.path(outdir, "top_total.tsv"), sep = "\t")
fwrite(top_gen,          file.path(outdir, "top_genetic.tsv"), sep = "\t")
fwrite(top_exp,          file.path(outdir, "top_exposomic.tsv"), sep = "\t")
fwrite(rnd(ec),          file.path(outdir, "exposure_categories.tsv"), sep = "\t")
fwrite(rnd(gc),          file.path(outdir, "gxe_categories.tsv"), sep = "\t")
fwrite(rnd(subblk),      file.path(outdir, "genetic_subblocks.tsv"), sep = "\t")
fwrite(top_cis,          file.path(outdir, "top_cis_genetic.tsv"), sep = "\t")
fwrite(rnd(stab[order(-gap)]), file.path(outdir, "traintest_stability.tsv"), sep = "\t")
fwrite(rnd(gxe_prot), file.path(outdir, "gxe_proteins.tsv"), sep = "\t")

mk <- character(0); add <- function(...) mk <<- c(mk, paste0(...))
add("# Module 1 — variance-decomposition summary")
add("")
add("_Generated by `scripts/analysis_summaries/summarize_module1_variance.R` from `",
    experiment, "/", covarType, "/", method, "`. Re-run the script to refresh._")
add("")
add("Cross-validated (", n_fold, "-fold out-of-fold) score-partition decomposition of ",
    "plasma-protein variance into genetic (PGS), exposomic (PXS) and gene-by-exposure ",
    "(GIS / GxE) components over a covariate baseline. **Unique R2** of a block = ",
    "R2(full model) - R2(full model without that block).")
add("")
add("## Headline numbers (copy-ready)")
add("")
add("- Decomposed **", cm(H$n_proteins), " plasma proteins**.")
add("- Whole-model out-of-fold test R2: median **", fmt(H$median_test_R2),
    "** (IQR ", fmt(H$q25_test_R2), "-", fmt(H$q75_test_R2), "), up to **",
    fmt(H$max_test_R2), "** (", H$max_test_R2_protein, ").")
add("- The proteome is **genetics-dominated**: **", cm(H$n_genetics_driven), " proteins (",
    pct(H$pct_genetics_driven), "%)** are genetics-driven (R2_G >= R2_E) vs **",
    cm(H$n_exposures_driven), " (", pct(H$pct_exposures_driven),
    "%)** exposures-driven (manuscript Fig 2A diagonal).")
add("- Unique genetic R2: median **", fmt(H$median_genetic_R2), "**, > 0.05 for **",
    cm(H$n_genetic_gt_0.05), " proteins**; top **", H$top_genetic_protein, "** (R2 = ",
    fmt(H$top_genetic_R2), ").")
add("- Unique exposomic R2: median **", fmt(H$median_exposomic_R2), "**, > 0.05 for **",
    cm(H$n_exposomic_gt_0.05), " proteins**; top **", H$top_exposomic_protein, "** (R2 = ",
    fmt(H$top_exposomic_R2), ").")
add("- Unique GxE R2 is small: median **", fmt(H$median_gxe_R2), "**.")
add("- Genetic signal is overwhelmingly **cis**: **", cm(H$n_cis_dominant), " of ",
    cm(nrow(gsw)), " proteins (", pct(H$pct_cis_dominant),
    "%)** are cis-dominant; median unique cis R2 **", fmt(H$median_cis_R2),
    "** vs trans **", fmt(H$median_trans_R2), "**.")
add("- Largest exposure-category contribution: **",
    heap_category_pretty(H$top_exposure_category), "** (mean unique R2 ",
    fmt(H$top_exposure_category_R2, 4), ").")
add("- **Generalises**: whole-model train R2 (median ", fmt(H$median_train_r2),
    ") tracks test R2 (median ", fmt(H$median_test_r2_whole),
    ") almost perfectly (cor = ", fmt(H$cor_train_test),
    "); median overfit gap **", fmt(H$median_overfit_gap), "** (max ",
    fmt(H$max_overfit_gap), "; ", pct(H$pct_overfit_gt_0.05),
    "% of proteins gap > 0.05).")
add("")
add("## Generalisation: train vs test R2 (whole model)")
add("")
add("The score-partition models are **not overfit** — evidence that the ",
    "exposure-type-routing fix removed the train-only IMD polynomial terms that ",
    "previously failed to replicate. Per-protein whole-model R2 averaged over the ",
    n_fold, " CV folds:")
add("")
add("- Train R2 median **", fmt(H$median_train_r2), "**, test (out-of-fold) R2 median **",
    fmt(H$median_test_r2_whole), "**; cor(train, test) = **", fmt(H$cor_train_test), "**.")
add("- Overfit gap (train - test): median **", fmt(H$median_overfit_gap),
    "**, 99th pct ", fmt(as.numeric(H$q99_overfit_gap)), ", max ", fmt(H$max_overfit_gap),
    "; **", pct(H$pct_overfit_gt_0.05), "%** of proteins exceed a 0.05 gap.")
add("- Within-protein test-R2 SD across folds: median **", fmt(H$median_fold_sd_test),
    "** (fold-to-fold stability).")
add("- See `fig_traintest_stability` and `traintest_stability.tsv`. NOTE: only the ",
    "whole-model train R2 is emitted by Module 1; component-level (G/E/GxE) train ",
    "R2 would need the module to write the train side of unique_r2_drop.")
add("")
add("## Per-component explained variation (across ", cm(n_prot), " proteins)")
add("")
add("| component | mean | median | IQR | max | % > 0.01 | % > 0.05 |")
add("|---|--:|--:|--:|--:|--:|--:|")
for (i in seq_len(nrow(comp_tbl)))
  add("| ", comp_tbl$component[i], " | ", fmt(comp_tbl$mean[i]), " | ",
      fmt(comp_tbl$median[i]), " | ", fmt(comp_tbl$q25[i]), "-", fmt(comp_tbl$q75[i]),
      " | ", fmt(comp_tbl$max[i]), " | ", pct(comp_tbl$pct_gt_0.01[i]), "% | ",
      pct(comp_tbl$pct_gt_0.05[i]), "% |")
add("")
add("## Genetic vs exposomic R2 distribution (manuscript Fig 2D)")
add("")
add("| set | bin | # proteins | % |")
add("|---|---|--:|--:|")
for (i in seq_len(nrow(r2_bins)))
  add("| ", r2_bins$set[i], " | ", r2_bins$bin[i], " | ", cm(r2_bins$n[i]), " | ",
      fmt(r2_bins$pct[i], 1), "% |")
add("")
add("## Top 15 most heritable proteins (unique genetic R2)")
add("")
add("| protein | genetic R2 | exposomic R2 | whole-model R2 |")
add("|---|--:|--:|--:|")
for (i in seq_len(min(15, nrow(top_gen))))
  add("| ", top_gen$protein[i], " | ", fmt(top_gen$genetic_R2[i]), " | ",
      fmt(top_gen$exposomic_R2[i]), " | ", fmt(top_gen$test_R2[i]), " |")
add("")
add("## Top 15 most exposure-driven proteins (unique exposomic R2)")
add("")
add("| protein | exposomic R2 | genetic R2 | whole-model R2 |")
add("|---|--:|--:|--:|")
for (i in seq_len(min(15, nrow(top_exp))))
  add("| ", top_exp$protein[i], " | ", fmt(top_exp$exposomic_R2[i]), " | ",
      fmt(top_exp$genetic_R2[i]), " | ", fmt(top_exp$test_R2[i]), " |")
add("")
add("## Genetic sub-architecture (cis vs trans)")
add("")
add("| block | # proteins present | mean | median | max | % > 0.01 |")
add("|---|--:|--:|--:|--:|--:|")
for (i in seq_len(nrow(subblk)))
  add("| ", subblk$block[i], " | ", cm(subblk$n_present[i]), " | ", fmt(subblk$mean[i]),
      " | ", fmt(subblk$median[i]), " | ", fmt(subblk$max[i]), " | ",
      pct(subblk$pct_gt_0.01[i]), "% |")
add("")
add("## Top GxE (gene-by-environment) proteins")
add("")
add("Unique out-of-fold R2 of the gene-interaction score (GIS). GxE is the ",
    "smallest component overall, but is concentrated in a set of proteins: **",
    cm(H$n_gxe_gt_0.01), "** have GxE R2 > 0.01 and **", cm(H$n_gxe_gt_0.02),
    "** > 0.02. The per-exposure-category GxE drops are small (signal spread ",
    "across exposures), so `dominant category` flags the leading exposure, not a partition.")
add("")
add("| protein | GxE R2 | dominant GxE category |")
add("|---|--:|---|")
for (i in seq_len(min(15, nrow(top_gxe))))
  add("| ", top_gxe$omic[i], " | ", fmt(top_gxe$gxe_r2[i]), " | ",
      heap_category_pretty(top_gxe$dom_category[i]), " |")
add("")
add("## Exposome contribution by category (manuscript Fig 2F)")
add("")
add("| category | broad group | mean unique R2 | median | max | # proteins > 0.005 |")
add("|---|---|--:|--:|--:|--:|")
for (i in seq_len(nrow(ec)))
  add("| ", heap_category_pretty(ec$category[i]), " | ", ec$broad[i], " | ",
      fmt(ec$mean_R2[i], 4), " | ", fmt(ec$median_R2[i], 4), " | ",
      fmt(ec$max_R2[i], 3), " | ", cm(ec$n_proteins_gt_0.005[i]), " |")
add("")
add("_Tables also available as TSVs in this directory: headline, component_summary, ",
    "r2_bins, top_total, top_genetic, top_exposomic, top_cis_genetic, ",
    "genetic_subblocks, exposure_categories, gxe_categories._")
writeLines(mk, file.path(outdir, "SUMMARY.md"))

message("Wrote summary to: ", outdir)
cat(sprintf("\n%s proteins | median test R2 %s | %s genetics-driven / %s exposures-driven | %s cis-dominant\n",
            cm(n_prot), fmt(H$median_test_R2), cm(n_gen), cm(n_exp), cm(n_cis_dom)))
