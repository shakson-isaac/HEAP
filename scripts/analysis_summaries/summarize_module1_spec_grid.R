#!/usr/bin/env Rscript
# ============================================================================
# summarize_module1_spec_grid.R
# ----------------------------------------------------------------------------
# Supplementary DATA (Tier 2) for the Module 1 variance decomposition: the same
# per-protein component R2 as Supplementary Table S1, under EVERY alternative
# specification. See docs/SUPPLEMENT_DATA_STANDARD.md §2.2 -- the Table reports
# the base specification; the Data file carries the alternatives, stacked, with
# explicit axis columns.
#
# Emits two long files mirroring S1's two sheets:
#   varcomp_specs.tsv       coarse (C/G/E/GxE) + cumulative model totals
#   varcomp_specs_fine.tsv  Gcis/Gtrans + per-category E + per-category GxE
#
# Schema (both files):
#   experiment covariate_set estimator sample interactions | level protein
#   component | r2 ci_lo ci_hi K
#
# `experiment` alone is NOT decodable -- five runs share covarType=base,
# method=lasso and differ only by sample or interaction terms -- so each axis is
# an explicit column. Each alternative varies exactly ONE axis from base.
#
# r2 = mean over K out-of-fold estimates; ci = t-based 95% CI on that mean.
# Unique R2 may be negative (component hurts out-of-fold); NOT clipped.
#
#   Run:  Rscript scripts/analysis_summaries/summarize_module1_spec_grid.R
# ============================================================================
suppressPackageStartupMessages(library(data.table))
options(scipen = 999)

local({
  cand <- c(Sys.getenv("HEAP_PATHS_FILE", unset = ""),
            "/n/groups/patel/shakson_ukb/HEAP/workflow/00_paths.R")
  hit <- cand[nzchar(cand) & file.exists(cand)][1]; if (!is.na(hit)) source(hit)
})
source(file.path(heap_root, "scripts", "visualizations", "common", "load_heap_results.R"))

OUT <- heap_project_output("supp_deposit")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

## the one-at-a-time sensitivity grid around the base specification -----------
## NOTE: prevalent disease is handled by EXCLUSION (M1_base_exclprev_lasso), not by
## adjustment. Conditioning on prevalent disease -- which lies downstream of exposure --
## can open a collider path, the very bias the reviews raised, so the +prevalent
## covariate specification is deliberately not released even though it was computed.
GRID <- fread(text = '
experiment|covarType|method|covariate_set|estimator|sample|interactions
M1_base_lasso|base|lasso|base|lasso|all|none
M1_base_bmi_lasso|base_bmi|lasso|base_bmi|lasso|all|none
M1_base_draw_lasso|base_draw|lasso|base_draw|lasso|all|none
M1_base_clinical_lasso|base_clinical|lasso|base_clinical|lasso|all|none
M1_base_exclprev_lasso|base|lasso|base|lasso|excl_prevalent|none
M1_base_ridge|base|ridge|base|ridge|all|none
M1_base_enet|base|enet|base|enet|all|none
M1_base_lasso_GxC|base|lasso|base|lasso|all|GxC
M1_base_lasso_ExC|base|lasso|base|lasso|all|ExC
M1_base_lasso_GxC_ExC|base|lasso|base|lasso|all|GxC+ExC
', sep = "|", header = TRUE)

agg <- function(d, by) {
  d <- d[is.finite(r2)]
  d[, .(K = .N, r2 = mean(r2), sd = sd(r2)), by = by
    ][, `:=`(ci_lo = r2 - qt(0.975, K - 1) * sd / sqrt(K),
             ci_hi = r2 + qt(0.975, K - 1) * sd / sqrt(K))][]
}
tidy <- function(d, lev, comp_col)
  d[, .(level = lev, protein = omic, component = get(comp_col), r2, ci_lo, ci_hi, K)]

## Guard against partial runs: skip a level that is missing/empty and record it,
## rather than aborting the whole grid on one incomplete run. (M1_base_enet used to
## contribute `coarse` only, because its fine levels were written as 0-byte files;
## re-run 2026-08-13 with the category decompositions enabled, so it now has all three.)
coverage <- list()
safe_ld <- function(g, l) {
  d <- tryCatch(load_module1_predictive_r2(g$covarType, g$method, level = l,
                                           experiment = g$experiment),
                error = function(e) NULL)
  ok <- !is.null(d) && nrow(d) > 0L && "method" %in% names(d)
  coverage[[length(coverage) + 1L]] <<-
    data.table(experiment = g$experiment, level = l, rows = if (ok) nrow(d) else 0L)
  if (!ok) { message("    ! ", l, " empty -- skipped"); return(NULL) }
  d
}

coarse_l <- list(); fine_l <- list()
for (i in seq_len(nrow(GRID))) {
  g  <- GRID[i]
  message(sprintf("[%2d/%d] %s", i, nrow(GRID), g$experiment))
  ax <- function(d) if (is.null(d) || !nrow(d)) NULL else
    d[, `:=`(experiment = g$experiment, covariate_set = g$covariate_set,
             estimator = g$estimator, sample = g$sample, interactions = g$interactions)]

  co <- safe_ld(g, "coarse")
  if (!is.null(co)) coarse_l[[length(coarse_l) + 1L]] <- ax(rbindlist(list(
    tidy(agg(co[method == "score_unique_drop"], c("omic", "block")), "coarse",      "block"),
    tidy(agg(co[method == "score_model_total"], c("omic", "block")), "model_total", "block"))))

  gs <- safe_ld(g, "genetic_subblocks")
  ec <- safe_ld(g, "exposure_categories")
  gc <- safe_ld(g, "gxe_categories")
  parts <- list()
  if (!is.null(gs)) parts[[length(parts)+1L]] <-
    tidy(agg(gs[method == "score_unique_drop"], c("omic","block")), "genetic_subblocks", "block")
  if (!is.null(ec)) parts[[length(parts)+1L]] <-
    tidy(agg(ec[method == "score_unique_drop"], c("omic","category")), "exposure_categories", "category")
  if (!is.null(gc)) parts[[length(parts)+1L]] <-
    tidy(agg(gc[method == "score_unique_drop"], c("omic","category")), "gxe_categories", "category")
  if (length(parts)) fine_l[[length(fine_l) + 1L]] <- ax(rbindlist(parts))
}

finish <- function(l, path) {
  d <- rbindlist(l, fill = TRUE)
  setcolorder(d, c("experiment", "covariate_set", "estimator", "sample", "interactions",
                   "level", "protein", "component", "r2", "ci_lo", "ci_hi", "K"))
  for (c in c("r2", "ci_lo", "ci_hi")) set(d, j = c, value = round(d[[c]], 4))
  setorder(d, protein, level, component, experiment)
  fwrite(d, path, sep = "\t")
  message(sprintf("  %-24s %s rows x %d cols  (%.0f MB)", basename(path),
                  format(nrow(d), big.mark = ","), ncol(d), file.info(path)$size / 1e6))
  d
}
a <- finish(coarse_l, file.path(OUT, "varcomp_specs.tsv"))
b <- finish(fine_l,   file.path(OUT, "varcomp_specs_fine.tsv"))

## coverage report -- which run/level combinations made it in, and which did not
cov <- rbindlist(coverage)
fwrite(cov, file.path(OUT, "varcomp_specs_coverage.tsv"), sep = "\t")
gaps <- cov[rows == 0L]
message(sprintf("\n%d experiments requested | coarse: %d | fine: %d",
                nrow(GRID), uniqueN(a$experiment), uniqueN(b$experiment)))
message(sprintf("%d proteins | %d coarse + %d fine components per protein",
                uniqueN(a$protein), uniqueN(a[, .(level, component)]),
                uniqueN(b[, .(level, component)])))
if (nrow(gaps)) {
  message("\nGAPS (level absent/empty on disk -- excluded from the Data files):")
  for (i in seq_len(nrow(gaps)))
    message(sprintf("  %-24s %s", gaps$experiment[i], gaps$level[i]))
  message("These must be stated in the Supplementary Data legend.")
}
