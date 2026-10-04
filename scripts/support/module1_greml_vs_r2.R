#!/usr/bin/env Rscript
# ============================================================================
# module1_greml_vs_r2.R  (support analysis for fig_greml_vs_r2)
# ----------------------------------------------------------------------------
# Cross-method concordance of the Module-1 variance decomposition: the GREML
# (multi-kernel REML) per-protein VARIANCE COMPONENTS (V/Vp for G, E, GxE) vs
# the HEAP per-protein PREDICTIVE R2 (unique drop-one, out-of-fold). These are
# two independent estimators of the same G/E/GxE partition:
#   * GREML  = how much variance EXISTS (the ceiling), variance-component REML.
#   * HEAP R2 = how much a cross-validated score PREDICTS out-of-sample.
#
# Inputs (canonical module output only):
#   GREML : population_architecture/<covar>/grm_cutoff_<cut>/<partition>/<P>_summary.tsv
#           (variance_G/E/GxE = V/Vp; converged flag; se_*)
#   HEAP  : module1_predictive_r2_score_partition/.../predictive_r2_coarse_*
#           (method == score_unique_drop, blocks G/E/GxE, mean over folds)
#
# Writes two cached tables the thin plotter reads (NO stats in the plotter):
#   <pop_arch>/<covar>/grm_cutoff_<cut>/concordance_greml_vs_heap_r2.tsv     (per protein x component)
#   <pop_arch>/<covar>/grm_cutoff_<cut>/concordance_greml_vs_heap_stats.tsv  (per component summary)
# ============================================================================
local({
  cm <- "/n/groups/patel/shakson_ukb/HEAP/scripts/visualizations/common"
  for (f in c("figure_paths.R","load_heap_results.R","label_helpers.R")) source(file.path(cm, f))
})
suppressPackageStartupMessages({ library(data.table) })

covarType <- Sys.getenv("PA_COVAR",     "base")
method    <- Sys.getenv("HEAP_METHOD",  "lasso")
grm_cut   <- Sys.getenv("PA_GRM_CUTOFF","0p025")
partition <- Sys.getenv("PA_PARTITION", "primary")

pa_dir <- file.path(heap_project_output("population_architecture"), covarType,
                    paste0("grm_cutoff_", grm_cut))
gdir   <- file.path(pa_dir, partition)
stopifnot(dir.exists(gdir))

# --- GREML variance components (converged proteins) -------------------------
fs <- list.files(gdir, pattern = "_summary\\.tsv$", full.names = TRUE)
message("GREML summaries: ", length(fs))
sel <- c("protein","variance_G","variance_E","variance_GxE","se_G","se_E","se_GxE","converged")
G <- rbindlist(lapply(fs, function(f) fread(f, select = sel, colClasses = "character")), fill = TRUE)
for (c in c("variance_G","variance_E","variance_GxE","se_G","se_E","se_GxE"))
  set(G, j = c, value = suppressWarnings(as.numeric(G[[c]])))
G[, protein := gsub('"', '', protein)]
G <- G[toupper(converged) == "TRUE"]
GL <- rbindlist(list(
  G[, .(protein, component = "G",   greml = variance_G,   greml_se = se_G)],
  G[, .(protein, component = "E",   greml = variance_E,   greml_se = se_E)],
  G[, .(protein, component = "GxE", greml = variance_GxE, greml_se = se_GxE)]))

# --- HEAP unique drop-one predictive R2 (per protein x block) ---------------
coarse <- load_module1_predictive_r2(covarType = covarType, method = method, level = "coarse")
ud <- coarse[get("method") == "score_unique_drop" & block %in% c("G","E","GxE"),
             .(r2 = mean(r2)), by = .(omic, block)]
setnames(ud, c("omic","block"), c("protein","component"))

# --- merge + write per-protein table ----------------------------------------
d <- merge(GL, ud, by = c("protein","component"))
d <- d[is.finite(greml) & is.finite(r2)]
d[, component := factor(component, levels = c("G","E","GxE"))]
setorder(d, component, -greml)

# --- per-component concordance stats (computed HERE, not in the plotter) -----
stats <- d[, {
  ok <- is.finite(greml) & is.finite(r2)
  .(n = sum(ok),
    pearson  = cor(greml[ok], r2[ok]),
    spearman = cor(greml[ok], r2[ok], method = "spearman"),
    median_greml = median(greml[ok]),
    median_r2    = median(r2[ok]),
    frac_r2_below_ceiling = mean(pmax(0, r2[ok]) <= greml[ok]))
}, by = component]

out_tab   <- file.path(pa_dir, "concordance_greml_vs_heap_r2.tsv")
out_stats <- file.path(pa_dir, "concordance_greml_vs_heap_stats.tsv")
fwrite(d,     out_tab,   sep = "\t")
fwrite(stats, out_stats, sep = "\t")
message("wrote ", out_tab, " (", nrow(d), " rows) and ", out_stats)
print(stats)
