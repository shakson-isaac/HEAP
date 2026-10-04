#!/usr/bin/env Rscript
# ============================================================================
# module6_holdout_ci.R
# ----------------------------------------------------------------------------
# Held-out generalization accuracy + a CI computed FROM the held-out data.
#
# For each exposure we read the person-level held-out scores (HoldoutScores.tsv:
# one row per eid x repeat-visit instance, with y_raw + pred_prot). We:
#   * point  = mean over the repeat visits of the per-visit held-out accuracy
#              (R^2 for continuous, AUC for binary)  -- the genuine out-of-sample
#              number, averaged over visits (instances 0/2/3).
#   * CI     = 95% bootstrap CI, resampling PEOPLE (eid clusters, so a person's
#              repeat visits move together), recomputing the same averaged statistic.
# Because the spread comes from resampling the ~3.4k held-out people (not the 50k
# training folds), the CI is the uncertainty of the held-out point itself, is
# computed identically for every exposure, and is visible (~+/-0.02-0.04).
#
# Writes one cached table the panel-b plotter reads:
#   <module6 out>/<covarType>/PESlong_<covarType>_HoldoutAccuracyCI.tsv
# ============================================================================
local({
  cm <- "/n/groups/patel/shakson_ukb/HEAP/scripts/visualizations/common"
  for (f in c("figure_paths.R","load_heap_results.R","label_helpers.R")) source(file.path(cm, f))
})
suppressPackageStartupMessages({ library(data.table) })

# covarType from argv so the reading metric can be produced for the alternative
# covariate specifications. Output stays inside the spec directory and carries the
# spec in its filename, so runs cannot collide.
args      <- commandArgs(trailingOnly = TRUE)
covarType <- if (length(args) && nzchar(args[1])) args[1] else "base"
B    <- as.integer(Sys.getenv("HOLDOUT_CI_B", "1000"))   # bootstrap reps
SEED <- 1L

outdir <- file.path(heap_project_output("module6_pes_longitudinal", covarType))
files  <- list.files(outdir, pattern = sprintf("^PESlong_%s_.*_HoldoutScores\\.tsv$", covarType), full.names = TRUE)
message("covarType = ", covarType)
message("found ", length(files), " HoldoutScores files; B=", B)

# --- metric helpers ---------------------------------------------------------
r2_fun <- function(y, p) {
  ok <- is.finite(y) & is.finite(p); y <- y[ok]; p <- p[ok]
  if (length(y) < 5 || var(y) == 0) return(NA_real_)
  1 - sum((y - p)^2) / sum((y - mean(y))^2)
}
# fast AUC (Mann-Whitney) on a binary y (0/1) scored by p
auc_fun <- function(y, p) {
  ok <- is.finite(y) & is.finite(p); y <- y[ok]; p <- p[ok]
  n1 <- sum(y == 1); n0 <- sum(y == 0)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  r <- rank(p)
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
# Average precision (area under the precision-recall curve), the step-wise
# estimator: AP = sum_k (R_k - R_{k-1}) * P_k. Reported alongside AUC because
# AUC is insensitive to class imbalance while AUPR is not -- several of these
# exposures have a positive class under 3%, where a high AUC can coexist with
# near-useless precision. The no-skill AUPR equals the prevalence, so AUPR must
# always be read against the prevalence column rather than against 0.5.
aupr_fun <- function(y, p) {
  ok <- is.finite(y) & is.finite(p); y <- y[ok]; p <- p[ok]
  n1 <- sum(y == 1)
  if (n1 == 0 || n1 == length(y)) return(NA_real_)
  o <- order(p, decreasing = TRUE); y <- y[o]
  tp <- cumsum(y == 1)
  prec <- tp / seq_along(y); rec <- tp / n1
  sum(diff(c(0, rec)) * prec)
}
# positive-class rate, averaged over repeat visits the same way the metrics are
prev_fun <- function(y, p) { y <- y[is.finite(y)]; if (!length(y)) return(NA_real_); mean(y == 1) }

# averaged-over-visits statistic for one (eid-subset of) the data
visit_mean_stat <- function(dt, metric, col = "pred_prot") {
  vals <- dt[, .(m = metric(y_raw, get(col))), by = instance]$m
  mean(vals, na.rm = TRUE)
}

one_exposure <- function(f) {
  x <- fread(f, select = c("eid","instance","exposure_id","exposure_type","y_raw","pred_prot","pred_cov","pred_full"))
  if (!nrow(x)) return(NULL)
  eid_chr <- x$exposure_id[1]; etype <- x$exposure_type[1]
  metric <- if (etype == "binary") auc_fun else r2_fun
  # Three nested models scored on the SAME held-out people, plus the increment
  # the proteome buys over covariates alone. The increment is bootstrapped as a
  # PAIRED difference -- both models recomputed inside each resample -- so its
  # interval accounts for the correlation between them. Taking the difference of
  # two separately-bootstrapped CIs would be far too wide.
  pt_prot <- visit_mean_stat(x, metric, "pred_prot")
  pt_cov  <- visit_mean_stat(x, metric, "pred_cov")
  pt_full <- visit_mean_stat(x, metric, "pred_full")
  pt_inc  <- pt_full - pt_cov
  # For binary exposures the SAME three nested models are also scored by average
  # precision, and the prevalence is carried so AUPR can be read against its
  # no-skill floor. Continuous exposures get NA here -- AUPR is undefined for them.
  is_bin  <- etype == "binary"
  pa_prot <- pa_cov <- pa_full <- pa_inc <- pa_prev <- NA_real_
  if (is_bin) {
    pa_prot <- visit_mean_stat(x, aupr_fun, "pred_prot")
    pa_cov  <- visit_mean_stat(x, aupr_fun, "pred_cov")
    pa_full <- visit_mean_stat(x, aupr_fun, "pred_full")
    pa_inc  <- pa_full - pa_cov
    pa_prev <- visit_mean_stat(x, prev_fun, "pred_prot")
  }
  eids <- unique(x$eid)
  setkey(x, eid)
  set.seed(SEED)
  bP <- bC <- bF <- bI <- numeric(B)
  aP <- aC <- aF <- aI <- aV <- rep(NA_real_, B)
  for (b in seq_len(B)) {
    samp <- sample(eids, length(eids), replace = TRUE)
    xb <- x[.(samp)]
    bP[b] <- visit_mean_stat(xb, metric, "pred_prot")
    bC[b] <- visit_mean_stat(xb, metric, "pred_cov")
    bF[b] <- visit_mean_stat(xb, metric, "pred_full")
    bI[b] <- bF[b] - bC[b]
    if (is_bin) {
      aP[b] <- visit_mean_stat(xb, aupr_fun, "pred_prot")
      aC[b] <- visit_mean_stat(xb, aupr_fun, "pred_cov")
      aF[b] <- visit_mean_stat(xb, aupr_fun, "pred_full")
      aI[b] <- aF[b] - aC[b]
      aV[b] <- visit_mean_stat(xb, prev_fun, "pred_prot")
    }
  }
  q <- function(v) { v <- v[is.finite(v)]
    if (length(v) > 10) quantile(v, c(.025,.975), names = FALSE) else c(NA_real_, NA_real_) }
  qP <- q(bP); qC <- q(bC); qF <- q(bF); qI <- q(bI)
  rP <- q(aP); rC <- q(aC); rF <- q(aF); rI <- q(aI); rV <- q(aV)
  data.table(exposure_id = eid_chr, exposure_type = etype,
             n_person = length(eids), n_obs = nrow(x),
             point = pt_prot, ci_lo = qP[1], ci_hi = qP[2], boot_sd = sd(bP, na.rm = TRUE),
             cov_point = pt_cov,  cov_lo  = qC[1], cov_hi  = qC[2],
             full_point = pt_full, full_lo = qF[1], full_hi = qF[2],
             increment = pt_inc,  increment_lo = qI[1], increment_hi = qI[2],
             increment_boot_sd = sd(bI, na.rm = TRUE),
             aupr_prot = pa_prot, aupr_prot_lo = rP[1], aupr_prot_hi = rP[2],
             aupr_cov  = pa_cov,  aupr_cov_lo  = rC[1], aupr_cov_hi  = rC[2],
             aupr_full = pa_full, aupr_full_lo = rF[1], aupr_full_hi = rF[2],
             aupr_increment = pa_inc, aupr_increment_lo = rI[1], aupr_increment_hi = rI[2],
             prevalence_holdout = pa_prev, prevalence_holdout_lo = rV[1],
             prevalence_holdout_hi = rV[2])
}

res <- rbindlist(lapply(files, one_exposure), fill = TRUE)
# attach category via the standard metadata helper
meta <- .heap_m6_exposure_meta(unique(res$exposure_id))
res <- merge(res, meta[, .(exposure_id, category)], by = "exposure_id", all.x = TRUE, sort = FALSE)

outfile <- file.path(outdir, sprintf("PESlong_%s_HoldoutAccuracyCI.tsv", covarType))
fwrite(res, outfile, sep = "\t")
message("wrote ", outfile, "  (", nrow(res), " exposures)")
print(res[order(-point)][1:8, .(exposure_id = substr(exposure_id,1,30), exposure_type,
        point = round(point,3), ci = sprintf("[%.3f, %.3f]", ci_lo, ci_hi), n_person)])
