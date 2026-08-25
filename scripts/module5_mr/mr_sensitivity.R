#!/usr/bin/env Rscript
# mr_sensitivity.R — extra per-edge MR sensitivity analyses, shared by
# Module5.R (split-sample UKB) and Module5_deCODE.R (deCODE SomaScan).
#
# Called from run_mr_edge() after the main MR + heterogeneity/pleiotropy outputs.
# Each analysis is best-effort (tryCatch) and guarded by instrument count, so a
# failure or an under-powered edge never blocks the main result. Writes, next to
# the edge's other outputs:
#   <prefix>_steiger.tsv      Steiger directionality (reverse-causation guard)
#   <prefix>_singlesnp.tsv    per-SNP Wald ratios
#   <prefix>_leaveoneout.tsv  leave-one-out IVW
#   <prefix>_presso.tsv       MR-PRESSO global + outlier/distortion (nsnp >= min)

suppressPackageStartupMessages(library(data.table))

run_mr_sensitivity <- function(dat, outdir, prefix,
                               presso_min_snp = 4L,
                               presso_nperm   = 1000L) {
  .w <- function(x, suffix) {
    if (!is.null(x) && is.data.frame(x) && nrow(x) > 0)
      fwrite(as.data.table(x), file.path(outdir, paste0(prefix, suffix)), sep = "\t")
  }
  nsnp <- nrow(dat)

  # Steiger directionality — instrument should explain more variance in the
  # exposure than the outcome (guards against reverse causation / wrong arrow).
  tryCatch(.w(TwoSampleMR::directionality_test(dat), "_steiger.tsv"),
           error = function(e) NULL)

  # single-SNP Wald ratios + leave-one-out IVW (need >= 2 instruments)
  if (nsnp >= 2L) {
    tryCatch(.w(TwoSampleMR::mr_singlesnp(dat),   "_singlesnp.tsv"),   error = function(e) NULL)
    tryCatch(.w(TwoSampleMR::mr_leaveoneout(dat), "_leaveoneout.tsv"), error = function(e) NULL)
  }

  # MR-PRESSO — global pleiotropy + outlier-corrected estimate + distortion test
  # (needs several instruments; permutation-based, so guarded by nsnp).
  if (nsnp >= presso_min_snp && requireNamespace("MRPRESSO", quietly = TRUE)) {
    tryCatch({
      d  <- as.data.frame(dat)
      pr <- MRPRESSO::mr_presso(
        BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure",
        SdOutcome   = "se.outcome",   SdExposure   = "se.exposure",
        OUTLIERtest = TRUE, DISTORTIONtest = TRUE,
        data = d, NbDistribution = presso_nperm, SignifThreshold = 0.05)
      res <- as.data.table(pr$`Main MR results`)
      gt  <- pr$`MR-PRESSO results`$`Global Test`
      dd  <- pr$`MR-PRESSO results`$`Distortion Test`
      res[, presso_global_rssobs   := if (!is.null(gt)) gt$RSSobs  else NA_real_]
      res[, presso_global_pval     := if (!is.null(gt)) gt$Pvalue  else NA_real_]
      res[, presso_distortion_pval := if (!is.null(dd)) dd$Pvalue  else NA_real_]
      .w(res, "_presso.tsv")
    }, error = function(e) NULL)
  }
  invisible(NULL)
}
