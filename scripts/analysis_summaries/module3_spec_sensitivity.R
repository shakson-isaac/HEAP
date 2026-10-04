#!/usr/bin/env Rscript
# ============================================================================
# module3_spec_sensitivity.R   (ANALYSIS -- support for fig_module3_spec_sensitivity)
# ----------------------------------------------------------------------------
# Robustness of the Module-3 mediation (NIE) estimates across covariate
# specifications, model families and sample definitions. For each spec it loads
# the natural indirect effect (NIE) rows and writes three tables for the plotter:
#   spec_summary.tsv       -- per spec (aggregate exposomic NIE, predictor=PXS_total):
#                             # FDR-sig mediations, Spearman(NIE logHR) vs base,
#                             % of base-sig pairs retained, median |NIE| ratio
#   attenuation.tsv        -- base-sig exposomic NIE pairs: NIE logHR in base vs the
#                             +BMI / +clinical specs (for the attenuation scatter)
#   category_attenuation.tsv -- per exposure CATEGORY (partitioned): median |NIE|
#                             ratio +clinical/base over base-sig category pairs
#                             (which exposure axes are BMI/clinical-mediated)
#
# Run: module load gcc/14.2.0 R/4.4.2; export HEAP_PATHS_FILE=.../00_paths.R
#      Rscript scripts/analysis_summaries/module3_spec_sensitivity.R
# ============================================================================
suppressPackageStartupMessages({ library(data.table) }); source(Sys.getenv("HEAP_PATHS_FILE"))
local({ c <- "/n/groups/patel/shakson_ukb/HEAP/scripts/visualizations/common"
  for (f in c("figure_paths","load_heap_results","label_helpers")) source(file.path(c, paste0(f, ".R"))) })
OUT <- file.path(heap_project_output("module3"), "spec_sensitivity"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
SEL <- c("protID","DZ_ID","predictor","predictor_class","effect_type","effect_logHR","delta_p","n_cases")
FDRQ <- 0.05; NMIN <- 100

# --- spec grid (primary_total, aggregate exposomic + genetic NIE) -----------
# Prevalent disease is shown as an EXCLUSION only -- conditioning on a variable
# downstream of exposure opens a collider path (the bias the reviews raised) and in
# the mediation deposit collapsed 22,270 significant links to 1. The deprivation-as-
# covariate estimand is likewise dropped: it moves a whole exposure domain into the
# covariate set, re-specifying the estimand rather than testing robustness.
SPECS <- data.table(
  exp  = c("M3_base_lasso_primary","M3_base_bmi_lasso_primary","M3_base_clinical_lasso_primary",
           "M3_base_exclprev_lasso_primary","M3_base_draw_lasso_primary",
           "M3_base_ridge_primary","M3_base_enet_primary"),
  cov  = c("base","base_bmi","base_clinical","base","base_draw","base","base"),
  fam  = c("lasso","lasso","lasso","lasso","lasso","ridge","enet"),
  lab  = c("Base","+ BMI","+ Clinical\n(BMI+meds)","Exclude\nprevalent",
           "+ Blood draw","Ridge score","Elastic-net\nscore"),
  kind = c("primary","covariate","covariate","sample","covariate","model","model"))

# NUM-6, resolved 2026-08-08. This used to call a link significant on q<0.05
# alone, with BH recomputed on the n_cases>=NMIN subset, giving 23,549 at base.
# The deposit (summarize_module3_mediation.R, and the med_specs_* folders) also
# requires a sign-consistent proportion mediated and runs BH over every link in
# the driver class: 22,270. About 5% apart at every spec. The two never shipped
# together so nothing surfaced it -- but activating \Tref{med_spec_sensitivity}
# alongside the deposit would have put two counts for one quantity in front of a
# reader. Harmonized onto the DEPOSIT rule, which is stricter and is what the
# already-cited med_exposure_total.tsv uses, so no published number moved.
# The rule itself is med_link_rules.R -- do not reimplement it here.
source(file.path(heap_root, "scripts", "analysis_summaries", "med_link_rules.R"))
SEL_LINK <- c("protID","DZ_ID","predictor","predictor_class","effect_type",
              "effect_logHR","delta_se","delta_p","n_cases","instrument_present")
load_nie <- function(i) {
  m <- tryCatch(load_module3_results(covarType = SPECS$cov[i], family = SPECS$fam[i],
                  mode = "primary_total", experiment = SPECS$exp[i], select = SEL_LINK),
                error = function(e) NULL)
  if (is.null(m) || !nrow(m)) return(NULL)
  m <- m[predictor == "PXS_total"]
  nie <- m[effect_type == "NIE", .(protID, DZ_ID, predictor_class, n_cases,
             nie_logHR = effect_logHR, nie_se = delta_se, nie_p = delta_p,
             driver_present = instrument_present)]
  nde <- m[effect_type == "NDE", .(protID, DZ_ID, nde_logHR = effect_logHR)]
  d <- merge(nie, nde, by = c("protID", "DZ_ID"))
  if (!nrow(d)) return(NULL)
  d <- med_derive_links(d)
  d[is.finite(nie_logHR), .(protID, DZ_ID, logHR = nie_logHR, sig, n_cases)]
}
nieL <- lapply(seq_len(nrow(SPECS)), load_nie); names(nieL) <- SPECS$exp
ok <- !vapply(nieL, is.null, logical(1)); if (!ok[1]) stop("base spec failed to load")
base <- nieL[["M3_base_lasso_primary"]][, .(protID, DZ_ID, logHR_base = logHR, sig_base = sig)]

summ <- rbindlist(lapply(which(ok), function(i){ e <- SPECS$exp[i]; s <- nieL[[e]]
  mg <- merge(base, s[, .(protID, DZ_ID, logHR_s = logHR, sig_s = sig)], by = c("protID","DZ_ID"))
  bs <- mg[sig_base == TRUE]
  data.table(exp = e, lab = SPECS$lab[i], kind = SPECS$kind[i],
    n_sig = s[, sum(sig)], n_sig_base = base[, sum(sig_base)],
    n_diseases = s[, uniqueN(DZ_ID)],
    total_cases = s[, .(c = as.integer(median(n_cases))), by = DZ_ID][, sum(c)],
    spearman_all = suppressWarnings(cor(mg$logHR_base, mg$logHR_s, method = "spearman", use = "complete.obs")),
    spearman_sig = suppressWarnings(cor(bs$logHR_base, bs$logHR_s, method = "spearman", use = "complete.obs")),
    pct_retained = 100 * bs[, mean(sig_s)],
    median_abs_ratio = bs[, median(abs(logHR_s)/abs(logHR_base), na.rm = TRUE)]) }))
fwrite(summ, file.path(OUT, "spec_summary.tsv"), sep = "\t")

# --- attenuation: base-sig exposomic NIE, base vs +BMI / +clinical ----------
att_spec <- function(e){ if (!e %in% names(nieL) || is.null(nieL[[e]])) return(NULL)
  merge(base[sig_base == TRUE], nieL[[e]][, .(protID, DZ_ID, logHR_s = logHR, sig_s = sig)],
        by = c("protID","DZ_ID"))[, .(protID, DZ_ID, logHR_base, logHR_s, sig_s, spec = e)] }
att <- rbindlist(lapply(c("M3_base_bmi_lasso_primary","M3_base_clinical_lasso_primary"), att_spec))
fwrite(att, file.path(OUT, "attenuation.tsv"), sep = "\t")

# --- per-category attenuation (+clinical/base), partitioned -----------------
load_cat <- function(exp, cov){
  m <- tryCatch(load_module3_results(covarType = cov, family = "lasso",
                  mode = "partitioned_categories", experiment = exp, select = SEL),
                error = function(e) NULL)
  if (is.null(m) || !nrow(m)) return(NULL)
  m <- m[effect_type == "NIE" & predictor_class == "exposure_category" &
         is.finite(effect_logHR) & is.finite(delta_p) & n_cases >= NMIN]
  m[, category := sub("^PXS_", "", predictor)]
  m[, q := p.adjust(delta_p, "BH"), by = category]
  m[, .(protID, DZ_ID, category, logHR = effect_logHR, sig = q < FDRQ)]
}
cb <- load_cat("M3_base_lasso_partitioned", "base")
cc <- load_cat("M3_base_clinical_lasso_partitioned", "base_clinical")
cat_att <- NULL
if (!is.null(cb) && !is.null(cc)) {
  mg <- merge(cb[sig == TRUE, .(protID, DZ_ID, category, logHR_base = logHR)],
              cc[, .(protID, DZ_ID, category, logHR_s = logHR, sig_s = sig)],
              by = c("protID","DZ_ID","category"))
  cat_att <- mg[, .(n_base_sig = .N,
                    median_ratio = median(abs(logHR_s)/abs(logHR_base), na.rm = TRUE),
                    pct_retained = 100 * mean(sig_s)), by = category][order(median_ratio)]
  fwrite(cat_att, file.path(OUT, "category_attenuation.tsv"), sep = "\t")
}

cat("specs loaded:", sum(ok), "/", nrow(SPECS), "\n")
print(summ[, .(lab = gsub("\n"," ",lab), n_sig, n_diseases, total_cases,
               spearman_sig = round(spearman_sig,3),
               pct_retained = round(pct_retained), med_ratio = round(median_abs_ratio,2))])
if (!is.null(cat_att)) { cat("\ncategory attenuation (+clinical/base):\n"); print(cat_att) }
cat("\nwrote", OUT, "/{spec_summary,attenuation,category_attenuation}.tsv\n")
