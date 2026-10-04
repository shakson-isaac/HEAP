#!/usr/bin/env Rscript
# ============================================================================
# module6_within_ci.R
# ----------------------------------------------------------------------------
# Within-person Delta-correlation with a CI computed FROM the held-out people.
#
# For each (continuous) exposure we read the person-level held-out scores and pair
# each person's baseline visit (instance 0) with each repeat visit (2,3):
#   dY    = y_raw(follow-up)   - y_raw(baseline)            change in actual exposure
#   dProt = pred_prot(follow)  - pred_prot(baseline)        change in proteome score
#   dCov  = pred_cov(follow)   - pred_cov(baseline)         change in covariate score
# delta-correlation = cor(dY, dProt)  (and the covariate benchmark cor(dY, dCov)).
# CI = 95% bootstrap, resampling PEOPLE (eid clusters), recomputing both.
#
# Writes:  <module6 out>/<covarType>/PESlong_<covarType>_WithinDeltaCorCI.tsv
# Mirrors module6_holdout_ci.R so panel c gets the same rigor as panel b.
# ============================================================================
local({
  cm <- "/n/groups/patel/shakson_ukb/HEAP/scripts/visualizations/common"
  for (f in c("figure_paths.R","load_heap_results.R","label_helpers.R")) source(file.path(cm, f))
})
suppressPackageStartupMessages({ library(data.table) })

# covarType from argv so the panel-c tracking metric can be produced for the
# alternative covariate specifications (base_bmi, base_clinical, base_draw,
# base_exclprev). Output stays INSIDE the spec directory and carries the spec in
# its filename, so no run can overwrite another's result.
args      <- commandArgs(trailingOnly = TRUE)
covarType <- if (length(args) && nzchar(args[1])) args[1] else "base"
B    <- as.integer(Sys.getenv("WITHIN_CI_B", "1000"))
SEED <- 1L

outdir <- file.path(heap_project_output("module6_pes_longitudinal", covarType))
files  <- list.files(outdir, pattern = sprintf("^PESlong_%s_.*_HoldoutScores\\.tsv$", covarType), full.names = TRUE)
message("covarType = ", covarType)
message("found ", length(files), " HoldoutScores files; B=", B)

cor_safe <- function(a, b) {
  ok <- is.finite(a) & is.finite(b); a <- a[ok]; b <- b[ok]
  if (length(a) < 8 || sd(a) == 0 || sd(b) == 0) return(NA_real_)
  cor(a, b)
}

one_exposure <- function(f) {
  x <- fread(f, select = c("eid","instance","exposure_id","exposure_type","y_raw","pred_prot","pred_cov","pred_full"))
  if (!nrow(x)) return(NULL)
  etype <- x$exposure_type[1]                          # continuous AND binary -- one Delta-corr scale
  b0 <- x[instance == 0, .(eid, y0 = y_raw, p0 = pred_prot, c0 = pred_cov, f0 = pred_full)]
  fu <- x[instance %in% c(2,3), .(eid, y1 = y_raw, p1 = pred_prot, c1 = pred_cov, f1 = pred_full)]
  d  <- merge(fu, b0, by = "eid")
  d[, `:=`(dY = y1 - y0, dP = p1 - p0, dC = c1 - c0, dF = f1 - f0)]
  d  <- d[is.finite(dY) & is.finite(dP) & is.finite(dC) & is.finite(dF)]
  # for binary one-hots dY is in {-1,0,+1}; cor(dY, dProt) is a point-biserial-style
  # within-person Delta-correlation on the same scale as the continuous case. Require a
  # minimum number of people who actually changed state, else the estimate is noise.
  n_change <- sum(d$dY != 0)
  if (nrow(d) < 30 || sd(d$dY) == 0 || n_change < 15) return(NULL)
  dcor_prot <- cor_safe(d$dY, d$dP); dcor_cov <- cor_safe(d$dY, d$dC)
  dcor_full <- cor_safe(d$dY, d$dF)   # covariates + proteome, the analogue of prot_plus_cov
  eids <- unique(d$eid); setkey(d, eid)
  set.seed(SEED); bp <- numeric(B); bc <- numeric(B); bf <- numeric(B)
  for (b in seq_len(B)) {
    s  <- d[.(sample(eids, length(eids), replace = TRUE))]
    bp[b] <- cor_safe(s$dY, s$dP); bc[b] <- cor_safe(s$dY, s$dC); bf[b] <- cor_safe(s$dY, s$dF)
  }
  qp <- quantile(bp[is.finite(bp)], c(.025,.975), names = FALSE)
  qc <- quantile(bc[is.finite(bc)], c(.025,.975), names = FALSE)
  qf <- quantile(bf[is.finite(bf)], c(.025,.975), names = FALSE)
  data.table(exposure_id = x$exposure_id[1], exposure_type = etype,
             n_person = length(eids), n_pairs = nrow(d), n_change = n_change,
             dcor_prot = dcor_prot, prot_lo = qp[1], prot_hi = qp[2],
             dcor_cov  = dcor_cov,  cov_lo  = qc[1], cov_hi  = qc[2],
             dcor_full = dcor_full, full_lo = qf[1], full_hi = qf[2])
}

res <- rbindlist(lapply(files, one_exposure), fill = TRUE)
meta <- .heap_m6_exposure_meta(unique(res$exposure_id))
res <- merge(res, meta[, .(exposure_id, category)], by = "exposure_id", all.x = TRUE, sort = FALSE)

outfile <- file.path(outdir, sprintf("PESlong_%s_WithinDeltaCorCI.tsv", covarType))
fwrite(res, outfile, sep = "\t")
message("wrote ", outfile, "  (", nrow(res), " exposures)")
print(res[order(-dcor_prot)][1:8, .(exposure_id = substr(exposure_id,1,30),
        prot = round(dcor_prot,3), prot_ci = sprintf("[%.2f,%.2f]", prot_lo, prot_hi),
        cov = round(dcor_cov,3), n_person)])
