#!/usr/bin/env Rscript
# ============================================================================
# module2_gxe_spec_sensitivity.R  (ANALYSIS -- support for fig_module2_gxe_spec_sensitivity)
# ----------------------------------------------------------------------------
# Robustness of the Module-2 polygenic GxE INTERACTIONS across covariate specs.
# GxE statistics are directionless F-tests (statFblock: p_GxE_joint / p_GcisxE /
# p_GtrxE) -- NO per-pair beta -- so robustness is p-value / replication based.
# A GxE pair is "replicated" when p_GxE_joint is Bonferroni-significant in BOTH
# 80/20 splits; cis / trans / joint-only by which component replicates. Writes:
#   gxe_spec_summary.tsv -- per spec: n replicated (+cis/trans/joint-only),
#                           Spearman of -log10 p_GxE_joint vs base, % retained
#   gxe_attenuation.tsv  -- base-replicated pairs: -log10 p_GxE_joint in base vs
#                           +BMI / +clinical (for the scatter; FOLR3/CCL3 hubs)
#
# Run: module load gcc/14.2.0 R/4.4.2; export HEAP_PATHS_FILE=.../00_paths.R
#      Rscript scripts/analysis_summaries/module2_gxe_spec_sensitivity.R
# ============================================================================
suppressPackageStartupMessages({ library(data.table) }); source(Sys.getenv("HEAP_PATHS_FILE"))
local({ c <- "/n/groups/patel/shakson_ukb/HEAP/scripts/visualizations/common"
  for (f in c("figure_paths","load_heap_results","label_helpers")) source(file.path(c, paste0(f, ".R"))) })
OUT <- file.path(heap_project_output("module2"), "spec_sensitivity"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# NB the +SES spec (M2_base_ses) is DELIBERATELY EXCLUDED from the GxE sensitivity
# analysis, for two reasons:
#   1. It is not a valid interaction test. A GxE test that conditions on a covariate C
#      must also carry the G x C and E x C terms, or the interaction estimate is biased.
#      base_ses adds SES MAIN EFFECTS only -- no G x SES, no E x SES. (We do apply exactly
#      this correction for age/sex: see fig_gxe_noise_floor panel b, which compares
#      base / +GxC / +ExC / +both.) So any GxE attrition under +SES is confounded with
#      model mis-specification and cannot be read as evidence of confounding.
#   2. It changes the estimand. base_ses calls E_to_Covariate(), which DELETES the
#      Deprivation_Indices category from the exposome and moves its 9 variables into the
#      covariate matrix. 16 of the 108 base-replicated GxE pairs are deprivation exposures
#      and are therefore untestable under +SES by construction, not "lost to adjustment".
# The spec is retained for the MAIN-EFFECT figure (fig_module2_spec_sensitivity), where
# reason 1 does not apply -- but it is labelled there as an alternative estimand
# (deprivation treated as covariate), not as a robustness check.
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

gxe_spec <- function(i) {
  tr <- tryCatch(load_module2_results(SPECS$cov[i], "train", experiment = SPECS$exp[i])$statFblock, error = function(e) NULL)
  te <- tryCatch(load_module2_results(SPECS$cov[i], "test",  experiment = SPECS$exp[i])$statFblock, error = function(e) NULL)
  if (is.null(tr) || is.null(te)) return(NULL)
  mg <- merge(tr[, .(ID, omicID, Category, jt_tr = p_GxE_joint, cis_tr = p_GcisxE, trn_tr = p_GtrxE)],
              te[, .(ID, omicID, jt_te = p_GxE_joint, cis_te = p_GcisxE, trn_te = p_GtrxE)], by = c("ID","omicID"))
  mg <- mg[is.finite(jt_tr) & is.finite(jt_te)]; thr <- 0.05 / nrow(mg)
  mg[, repl := jt_tr < thr & jt_te < thr]
  mg[, cisR := is.finite(cis_tr) & is.finite(cis_te) & cis_tr < thr & cis_te < thr]
  mg[, trnR := is.finite(trn_tr) & is.finite(trn_te) & trn_tr < thr & trn_te < thr]
  mg[, cls := fcase(repl & cisR, "cis", repl & trnR & !cisR, "trans", repl, "joint-only", default = "ns")]
  mg[, logp := -log10(pmin(jt_tr, jt_te))]; mg[, thr := thr]
  mg[]
}
specL <- lapply(seq_len(nrow(SPECS)), gxe_spec); names(specL) <- SPECS$exp
ok <- !vapply(specL, is.null, logical(1)); if (!ok[1]) stop("base GxE spec failed to load")
base <- specL[["M2_base_main"]][, .(ID, omicID, Category, repl_base = repl, cls_base = cls, logp_base = logp)]

summ <- rbindlist(lapply(which(ok), function(i){ e <- SPECS$exp[i]; s <- specL[[e]]
  mg <- merge(base, s[, .(ID, omicID, repl_s = repl, logp_s = logp)], by = c("ID","omicID"))
  br <- mg[repl_base == TRUE]
  data.table(exp = e, lab = SPECS$lab[i], kind = SPECS$kind[i],
    n_repl = s[, sum(repl)], n_cis = s[, sum(cls == "cis")], n_trans = s[, sum(cls == "trans")],
    n_jointonly = s[, sum(cls == "joint-only")], n_repl_base = base[, sum(repl_base)],
    spearman_all = suppressWarnings(cor(mg$logp_base, mg$logp_s, method = "spearman", use = "complete.obs")),
    spearman_repl = suppressWarnings(cor(br$logp_base, br$logp_s, method = "spearman", use = "complete.obs")),
    pct_retained = 100 * br[, mean(repl_s)]) }))
fwrite(summ, file.path(OUT, "gxe_spec_summary.tsv"), sep = "\t")

att_spec <- function(e){ if (!e %in% names(specL) || is.null(specL[[e]])) return(NULL)
  merge(base[repl_base == TRUE], specL[[e]][, .(ID, omicID, logp_s = logp, repl_s = repl)], by = c("ID","omicID"))[
    , .(ID, omicID, Category, cls_base, logp_base, logp_s, repl_s, spec = e)] }
att <- rbindlist(lapply(c("M2_base_bmi","M2_base_clinical_main"), att_spec))
fwrite(att, file.path(OUT, "gxe_attenuation.tsv"), sep = "\t")

cat("specs loaded:", sum(ok), "/", nrow(SPECS), "\n")
print(summ[, .(lab = gsub("\n"," ",lab), n_repl, n_cis, n_trans, n_jointonly, sp_repl = round(spearman_repl,3), pct_ret = round(pct_retained))])
cat("\nwrote", OUT, "/{gxe_spec_summary,gxe_attenuation}.tsv\n")
