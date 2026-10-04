#!/usr/bin/env Rscript
# ============================================================================
# module6_covariate_sensitivity.R  (support)
# ----------------------------------------------------------------------------
# Held-out READS accuracy of covariates-only vs covariates+PES for every
# exposure under each COVARIATE SPEC (the sensitivity runs). The proteome PES
# is identical across specs; what changes is the covariate model, so this
# isolates how much the PES ADDS beyond covariates and whether that survives
# adding BMI / clinical covariates.
#   For each exposure x spec (visit-averaged over held-out instances 0,2,3):
#     cov_only      = R2/AUC of covariate prediction of the exposure
#     cov_plus_pes  = R2/AUC of covariate+proteome prediction
#     incremental   = cov_plus_pes - cov_only
# Writes  module6_pes_longitudinal/covariate_sensitivity/covariate_sensitivity.tsv
# ============================================================================
local({ cand <- c(Sys.getenv("HEAP_PATHS_FILE",""),"/n/groups/patel/shakson_ukb/HEAP/workflow/00_paths.R")
  h <- cand[nzchar(cand)&file.exists(cand)][1]; if(!is.na(h)) source(h) })
suppressPackageStartupMessages({ library(data.table) })
for (f in c("figure_paths.R","load_heap_results.R","label_helpers.R"))
  source(file.path("/n/groups/patel/shakson_ukb/HEAP/scripts/visualizations/common", f))

root  <- heap_project_output("module6_pes_longitudinal", "base"); root <- dirname(root)
specs <- c("base","base_bmi","base_clinical","base_draw","base_exclprev")

r2  <- function(y,p){ ok<-is.finite(y)&is.finite(p); y<-y[ok]; p<-p[ok]
  if(length(y)<30 || sd(y)==0) return(NA_real_); 1 - sum((y-p)^2)/sum((y-mean(y))^2) }
auc <- function(y,p){ ok<-is.finite(y)&is.finite(p); y<-y[ok]; p<-p[ok]
  n1<-sum(y==1); n0<-sum(y==0); if(n1==0||n0==0) return(NA_real_)
  r<-rank(p); (sum(r[y==1]) - n1*(n1+1)/2)/(n1*n0) }
# average precision (AUPR estimator): mean of precision@k over the positives.
# Under the null (random score) its expectation is the positive-class prevalence.
aupr <- function(y,p){ ok<-is.finite(y)&is.finite(p); y<-y[ok]; p<-p[ok]
  n1<-sum(y==1); if(n1==0||n1==length(y)) return(NA_real_)
  o<-order(p, decreasing=TRUE); yy<-y[o]
  prec<-cumsum(yy==1)/seq_along(yy); sum(prec[yy==1])/n1 }
prev_pos <- function(y){ y<-y[is.finite(y)]; if(!length(y)) NA_real_ else mean(y==1) }

one <- function(f, spec) {
  x <- tryCatch(fread(f, select=c("instance","exposure_id","exposure_type","y_raw","pred_cov","pred_full")), error=function(e) NULL)
  if (is.null(x) || !nrow(x)) return(NULL)
  et <- x$exposure_type[1]; metric <- if (et=="binary") auc else r2; bin <- et=="binary"
  # visit-averaged (held-out instances 0,2,3), applying fn(y, pred_col)
  va  <- function(fn, col) mean(sapply(c(0,2,3), function(i){ d<-x[instance==i]; if(nrow(d)<30) NA_real_ else fn(d$y_raw, d[[col]]) }), na.rm=TRUE)
  prev <- if (bin) mean(sapply(c(0,2,3), function(i){ d<-x[instance==i]; if(nrow(d)<30) NA_real_ else prev_pos(d$y_raw) }), na.rm=TRUE) else NA_real_
  data.table(exposure_id=x$exposure_id[1], spec=spec, exposure_type=et,
             cov_only=va(metric,"pred_cov"), cov_plus_pes=va(metric,"pred_full"),
             prevalence=prev,
             aupr_cov  = if (bin) va(aupr,"pred_cov")  else NA_real_,
             aupr_full = if (bin) va(aupr,"pred_full") else NA_real_)
}

res <- rbindlist(lapply(specs, function(s){
  fs <- list.files(file.path(root, s), pattern="_HoldoutScores\\.tsv$", full.names=TRUE)
  if(!length(fs)) return(NULL)
  cat(sprintf("[%s] %d exposures\n", s, length(fs)))
  rbindlist(lapply(fs, one, spec=s), fill=TRUE)
}), fill=TRUE)
res[, incremental := cov_plus_pes - cov_only]
res[, incremental_aupr := aupr_full - aupr_cov]   # what the PES adds in precision-recall
res[, aupr_lift := aupr_full - prevalence]        # cov+PES AUPR above the null baseline (prevalence)
meta <- .heap_m6_exposure_meta(unique(res$exposure_id))
res <- merge(res, meta[, .(exposure_id, category)], by="exposure_id", all.x=TRUE, sort=FALSE)

outdir <- file.path(root, "covariate_sensitivity"); dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
fwrite(res, file.path(outdir, "covariate_sensitivity.tsv"), sep="\t")
cat("wrote ", file.path(outdir, "covariate_sensitivity.tsv"), " (", nrow(res), " rows)\n", sep="")
