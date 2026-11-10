#!/usr/bin/env Rscript
# ============================================================================
# collect_dc_holdout_ci.R
# ----------------------------------------------------------------------------
# Stacks the per-exposure outputs of module6_pes_disease_dc_holdout_ci.R into
# one table and checks it against the screen it is meant to annotate.
#
# Two things are asserted rather than assumed:
#   1. COVERAGE  -- every pair in pes_disease_dc_holdout.tsv gets a row here.
#   2. FIDELITY  -- the point estimates reproduce that file EXACTLY. They must,
#                   because the bootstrap resamples the test set of an unchanged
#                   fit; any drift means the 70/30 split was not reproduced and
#                   the intervals do not belong to the published numbers.
# A failure of either is fatal: shipping intervals attached to the wrong point
# estimate is worse than shipping no intervals at all.
#
#   Usage:  Rscript collect_dc_holdout_ci.R
#   Writes: multipes_disease/pes_disease_dc_holdout_ci.tsv
# ============================================================================
local({ cand <- c(Sys.getenv("HEAP_PATHS_FILE",""),"/n/groups/patel/shakson_ukb/HEAP/workflow/00_paths.R")
  h <- cand[nzchar(cand)&file.exists(cand)][1]; if(!is.na(h)) source(h) })
suppressPackageStartupMessages({ library(data.table) })
msg <- function(...) cat(format(Sys.time(),"[%H:%M:%S] "),...,"\n",sep="")

out_dir <- heap_project_output("module6_pes_longitudinal","multipes_disease")
ci_dir  <- file.path(out_dir,"dc_holdout_ci")
fs <- list.files(ci_dir, pattern="\\.tsv$", full.names=TRUE)
msg(length(fs), " per-exposure files")
if (!length(fs)) stop("no per-exposure outputs in ", ci_dir)

R <- rbindlist(lapply(fs, fread), fill=TRUE)
setorder(R, exposure_id, disease)
msg(nrow(R), " pairs | ", uniqueN(R$exposure_id), " exposures x ", uniqueN(R$disease), " diseases")

ref <- fread(file.path(out_dir,"pes_disease_dc_holdout.tsv"))
msg("screen has ", nrow(ref), " pairs")

key <- function(d) paste(d$exposure_id, d$disease)
miss <- setdiff(key(ref), key(R))
if (length(miss)) {
  msg("  ! MISSING ", length(miss), " pairs -- first few:"); print(head(miss, 5))
} else msg("  coverage: complete")

m <- merge(R[,.(exposure_id,disease,C0_ho,CP_ho,CE_ho,CEP_ho)],
           ref[,.(exposure_id,disease,r0=C0_ho,rP=CP_ho,rE=CE_ho,rEP=CEP_ho)],
           by=c("exposure_id","disease"))
dmax <- max(abs(c(m$C0_ho-m$r0, m$CP_ho-m$rP, m$CE_ho-m$rE, m$CEP_ho-m$rEP)), na.rm=TRUE)
msg(sprintf("  fidelity: %d matched | max|diff| = %.3e  %s", nrow(m), dmax,
            if (is.finite(dmax) && dmax < 1e-9) "EXACT" else "*** MISMATCH ***"))
if (!(is.finite(dmax) && dmax < 1e-9))
  stop("point estimates do not reproduce the screen -- the split was not reproduced; refusing to write")
if (length(miss)) stop("incomplete coverage -- rerun the missing array tasks before collecting")

fwrite(R, file.path(out_dir,"pes_disease_dc_holdout_ci.tsv"), sep="\t")
msg("wrote pes_disease_dc_holdout_ci.tsv (", nrow(R), " pairs, ", ncol(R), " cols)")

w <- R[is.finite(dC_overcov_ho_lo)]
msg(sprintf("  dC over covariates: %d of %d exclude 0 (%.0f%%)",
            sum(w$dC_overcov_ho_lo > 0), nrow(w), 100*mean(w$dC_overcov_ho_lo > 0)))
