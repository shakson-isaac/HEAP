#!/usr/bin/env Rscript
# ============================================================================
# module2_spec_sensitivity.R   (ANALYSIS -- support for fig_module2_spec_sensitivity)
# ----------------------------------------------------------------------------
# Robustness of the Module-2 exposure->protein associations across covariate
# specifications. Loads each spec's replicated statE (load_module2_replicated),
# and writes two tables consumed by the plotter:
#   spec_summary.tsv      -- per spec: n replicated, Spearman(beta) vs base,
#                            % of base-replicated retained, median |beta| ratio
#   attenuation.tsv       -- per (protein x exposure) base-replicated pair:
#                            beta in base vs each adjusted spec (for the scatter /
#                            most-attenuated panel; BMI-mediated drop e.g. LEP/FABP4)
#
# Run: module load gcc/14.2.0 R/4.4.2; export HEAP_PATHS_FILE=.../00_paths.R
#      Rscript scripts/analysis_summaries/module2_spec_sensitivity.R
# ============================================================================
suppressPackageStartupMessages({ library(data.table) }); source(Sys.getenv("HEAP_PATHS_FILE"))
local({ c <- "/n/groups/patel/shakson_ukb/HEAP/scripts/visualizations/common"
  for (f in c("figure_paths","load_heap_results","label_helpers")) source(file.path(c, paste0(f, ".R"))) })
OUT <- file.path(heap_project_output("module2"), "spec_sensitivity"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# NB "Deprivation as covariate" (M2_base_ses) is NOT a robustness check -- it is an
# ALTERNATIVE ESTIMAND, and the label says so. base_ses calls E_to_Covariate(), which
# DELETES the Deprivation_Indices category from the exposome and moves its 9 variables
# (household income, England IMD + its 7 domain sub-scores) into the covariate matrix.
# So this spec estimates the effect of a 12-category exposome conditional on deprivation,
# not the effect of the 13-category exposome. Because area deprivation plausibly lies
# UPSTREAM of lifestyle (deprivation -> smoking/diet/exercise -> protein), conditioning on
# it blocks a mediated path and removes real exposure signal: association counts fall ~57%
# while the surviving betas are essentially unchanged (Spearman ~0.99), which is the
# signature of over-adjustment, not of confounding. Excluded entirely from the GxE
# sensitivity figure (see module2_gxe_spec_sensitivity.R for why).
# Prevalent disease is shown as an EXCLUSION only -- conditioning on a variable
# downstream of exposure opens a collider path (the bias the reviews raised) and in
# the mediation deposit collapsed 22,270 significant links to 1. The deprivation-as-
# covariate estimand is likewise dropped: it moves a whole exposure domain into the
# covariate set, re-specifying the estimand rather than testing robustness.
SPECS <- data.table(
  exp = c("M2_base_main","M2_base_bmi","M2_base_clinical_main","M2_base_draw","M2_base_exclprev"),
  cov = c("base","base_bmi","base_clinical","base_draw","base"),
  lab = c("Base","+ BMI","+ Clinical\n(BMI+meds)","+ Blood draw","Exclude\nprevalent"),
  kind = c("primary","covariate","covariate","covariate","sample"))

load_spec <- function(i) {
  m <- tryCatch(load_module2_replicated(SPECS$cov[i], experiment = SPECS$exp[i]), error = function(e) NULL)
  if (is.null(m) || !nrow(m)) return(NULL)
  m <- m[is.finite(beta_train)]
  m[, .(ID, omicID, Category, beta = beta_train, repl = replicated == TRUE)]
}
specL <- lapply(seq_len(nrow(SPECS)), load_spec); names(specL) <- SPECS$exp
ok <- !vapply(specL, is.null, logical(1)); if (!ok[1]) stop("base spec failed to load")
base <- specL[["M2_base_main"]][, .(ID, omicID, Category, beta_base = beta, repl_base = repl)]

summ <- rbindlist(lapply(which(ok), function(i){ e <- SPECS$exp[i]; s <- specL[[e]]
  mg <- merge(base, s[, .(ID, omicID, beta_s = beta, repl_s = repl)], by = c("ID","omicID"))
  br <- mg[repl_base == TRUE]
  data.table(exp = e, lab = SPECS$lab[i], kind = SPECS$kind[i],
    n_repl = s[, sum(repl)], n_repl_base = base[, sum(repl_base)],
    spearman_all = suppressWarnings(cor(mg$beta_base, mg$beta_s, method = "spearman", use = "complete.obs")),
    spearman_repl = suppressWarnings(cor(br$beta_base, br$beta_s, method = "spearman", use = "complete.obs")),
    pct_retained = 100 * br[, mean(repl_s)],
    median_abs_ratio = br[, median(abs(beta_s)/abs(beta_base), na.rm = TRUE)]) }))
fwrite(summ, file.path(OUT, "spec_summary.tsv"), sep = "\t")

# attenuation: base-replicated pairs, beta in base vs the BMI & clinical specs
att_spec <- function(e){ if (!e %in% names(specL) || is.null(specL[[e]])) return(NULL)
  merge(base[repl_base == TRUE], specL[[e]][, .(ID, omicID, beta_s = beta, repl_s = repl)], by = c("ID","omicID"))[
    , .(ID, omicID, Category, beta_base, beta_s, repl_s, spec = e)] }
att <- rbindlist(lapply(c("M2_base_bmi","M2_base_clinical_main"), att_spec))
fwrite(att, file.path(OUT, "attenuation.tsv"), sep = "\t")

cat("specs loaded:", sum(ok), "/", nrow(SPECS), "\n"); print(summ[, .(lab=gsub("\n"," ",lab), n_repl, spearman_repl=round(spearman_repl,3), pct_retained=round(pct_retained), med_ratio=round(median_abs_ratio,2))])
cat("\nwrote", OUT, "/{spec_summary,attenuation}.tsv\n")
