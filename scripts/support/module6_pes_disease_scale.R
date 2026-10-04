#!/usr/bin/env Rscript
# ============================================================================
# module6_pes_disease_scale.R  (support: PES x disease prediction AT SCALE)
# ----------------------------------------------------------------------------
# For EVERY exposure PES (out-of-fold, pes_prot_z) and EVERY disease in the
# 181-disease >=100-case mediation set, fit FULL-COHORT Cox models and report
# the INFERENTIAL contribution of the PES:
#
#   base = age, sex, age^2, age*sex, age^2*sex, centre(factor), 20 genetic PCs
#   M0: ~ base                    M1: ~ base + PES
#   M2: ~ base + E                M3: ~ base + E + PES         (E = self-report)
#
# Per exposure x disease we report:
#   * HR per SD of PES over covariates        (M1 coef; Wald + LRT M0 vs M1)
#   * HR per SD of PES beyond covariates + E  (M3 coef; Wald + LRT M2 vs M3)
#   * apparent C-index for M0..M3             (flagged apparent; bootstrap is
#                                              run separately for headline pairs)
#
# Why full-cohort (no 70/30 re-split): the PES is ALREADY out-of-fold (read from
# *_TrainOOF.tsv), so the high-dimensional proteome->score model cannot leak. The
# remaining quantity is a single pre-computed covariate; testing its Cox
# coefficient by likelihood-ratio is valid and fully powered on all events. The
# old quadrant-scan 70/30 split only de-optimized the ~46-covariate Cox at the
# cost of 70% of events -> that is what capped the analysis at 16 diseases.
# C-index optimism (the one thing fit-on-all inflates) is handled downstream by
# an optimism-corrected bootstrap for the highlighted pairs only.
#
# BH-FDR is applied across the full exposure x disease family for each LRT.
# Writes multipes_disease/pes_disease_scale.tsv.
# ============================================================================
local({ cand <- c(Sys.getenv("HEAP_PATHS_FILE",""),"/n/groups/patel/shakson_ukb/HEAP/workflow/00_paths.R")
  h <- cand[nzchar(cand)&file.exists(cand)][1]; if(!is.na(h)) source(h) })
suppressPackageStartupMessages({ library(data.table); library(survival) })
msg <- function(...) cat(format(Sys.time(),"[%H:%M:%S] "),...,"\n",sep="")
RNGkind("L'Ecuyer-CMRG"); set.seed(42)

args      <- commandArgs(trailingOnly = TRUE)
covarType <- if (length(args) && nzchar(args[1])) args[1] else "base"
ncore  <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "4"))
MINEVT <- as.integer(Sys.getenv("HEAP_PES_MIN_EVENTS", unset = "100"))
# base_exclprev is NOT a covariate set. It is the base adjustment applied to the
# healthy-at-baseline subset, scored with the PES retrained on that subset -- two
# changes at once (sample and score), where the covariate specs change only the
# adjustment. Resolve the three axes explicitly rather than assuming covarType
# names a covariate set.
EXCLPREV  <- identical(covarType, "base_exclprev")
covarSet  <- if (EXCLPREV) "base" else covarType
scoreDir  <- if (EXCLPREV) "base_exclprev" else "base"
od      <- heap_project_output("module6_pes_longitudinal", scoreDir)
out_dir <- heap_project_output("module6_pes_longitudinal","multipes_disease")

# --- 181-disease set (DZ_ID == the age_*_first_reported_* column name) -------
dzset <- fread(file.path(out_dir,"disease_set_181.tsv"))$DZ_ID
NDZ  <- as.integer(Sys.getenv("HEAP_PES_NDZ",  unset = "0"))   # 0 = all (test knob)
if (NDZ > 0) dzset <- head(dzset, NDZ)
msg("disease set: ", length(dzset), " diseases")

# --- covariates (base spec, instance 0) -------------------------------------
msg("Loading HEAP.rds")
heap <- readRDS(heap_loader_rds)
dz   <- as.data.table(heap$disease$DZ_df)
cvl  <- as.data.table(heap$covars_long); rm(heap); gc()
# Cox adjustment is now driven by config/covariates/covariate_sets.yml rather
# than hardcoded here. The old literal was exactly the `base` set, and
# assert_base_formula_unchanged() enforces that it still is, so covarType="base"
# must reproduce every existing number.
source(file.path(dirname(Sys.getenv("HEAP_PATHS_FILE","/n/groups/patel/shakson_ukb/HEAP/workflow/00_paths.R")),
                 "config_helpers.R"))
source(file.path(Sys.getenv("HEAP_ROOT", "/n/groups/patel/shakson_ukb/HEAP"), "scripts/support/module6_disease_covars.R"))
zc    <- build_disease_covars(cvl, covar_set = covarSet)
# Use the UNFILTERED frame, matching the original: coxph drops incomplete rows at
# fit time via na.action. Pre-filtering to complete cases would be tidier but
# changes the sample the PES is standardised over, so `base` would no longer
# reproduce the published numbers. zc$n_complete records what the fits actually use.
cov0  <- zc$dt
if (EXCLPREV) {
  fcol <- "prevalent_major_disease"
  if (!fcol %in% names(cvl)) stop("sample filter needs `", fcol, "` in covars_long")
  keep <- cvl[[fcol]]; ids <- cvl$eid[cvl$instance == 0 & !is.na(keep) & keep == 0]
  n0 <- nrow(cov0); cov0 <- cov0[eid %in% ids]
  msg("  exclude_prevalent: ", n0, " -> ", nrow(cov0), " participants (healthy at baseline)")
}
basef <- zc$formula
if (identical(covarType, "base")) assert_base_formula_unchanged(basef)
msg("covariate set: ", covarType, " | ", length(zc$covars), " terms | ",
    nrow(cov0), " complete rows")
for(c in c("recode_age_of_death_0_0","age_of_removal_0_0","age_of_lastfollowup","recode_age_of_assessment_0_0"))
  if(c %in% names(dz)) dz[,(c):=as.numeric(get(c))]
dz[, censor_age := pmin(recode_age_of_death_0_0, age_of_removal_0_0, age_of_lastfollowup, na.rm=TRUE)]

# --- per-disease risk-set base tables (built once, shared by fork) ----------
base_by_dz <- list()
for(ac in dzset){
  if(!ac %in% names(dz)){ msg("MISSING disease col: ", ac); next }
  b <- merge(dz[, .(eid, event_age=as.numeric(get(ac)), asmt=recode_age_of_assessment_0_0, censor_age)], cov0, by="eid")
  b <- b[is.finite(asmt) & is.finite(censor_age) & censor_age>asmt]
  b <- b[is.na(event_age) | event_age>asmt]
  b[, status := as.integer(!is.na(event_age) & event_age<=censor_age)]
  b[, time := fifelse(status==1, event_age, censor_age)-asmt]
  b <- b[time>0]
  if(sum(b$status) < MINEVT) next
  # store ONLY (eid,status,time) per disease; covariates live once in cov0 and
  # are merged in per pair. (Storing the full covariate matrix x181 diseases
  # makes forked mclapply workers copy-on-write ~2GB each -> OOM.)
  base_by_dz[[ac]] <- b[, .(eid, status, time)]
}
msg("disease base tables passing >=", MINEVT, " events: ", length(base_by_dz), " of ", length(dzset))

appC <- function(time,status,lp){ ok<-is.finite(time)&is.finite(status)&is.finite(lp)
  if(sum(ok)<50 || sum(status[ok])<10) return(NA_real_)
  tryCatch(survival::concordance(Surv(time[ok],status[ok])~lp[ok], reverse=TRUE)$concordance, error=function(e) NA_real_) }

# pull the PES coefficient (term "PES") -> HR per SD, Wald p, and LRT vs reduced
pes_inf <- function(full, reduced){
  cf <- tryCatch(coef(full), error=function(e) NULL)
  if(is.null(cf) || !("PES" %in% names(cf)) || !is.finite(cf[["PES"]]))
    return(list(hr=NA_real_, lo=NA_real_, hi=NA_real_, wp=NA_real_, lrt=NA_real_))
  b  <- cf[["PES"]]; se <- sqrt(diag(vcov(full)))[["PES"]]
  wp <- 2*pnorm(-abs(b/se))
  lr <- 2*(full$loglik[2] - reduced$loglik[2])
  lrt <- if(is.finite(lr) && lr>=0) pchisq(lr, df=1, lower.tail=FALSE) else NA_real_
  list(hr=exp(b), lo=exp(b-1.96*se), hi=exp(b+1.96*se), wp=wp, lrt=lrt)
}

fs <- list.files(od, pattern="_TrainOOF\\.tsv$", full.names=TRUE)
NEXP <- as.integer(Sys.getenv("HEAP_PES_NEXP", unset = "0"))   # 0 = all (test knob)
if (NEXP > 0) fs <- head(fs, NEXP)
ge <- function(f) sub(sprintf("^PESlong_%s_(.*)_TrainOOF\\.tsv$", scoreDir), "\\1", basename(f))
msg(length(fs), " exposures x ", length(base_by_dz), " diseases = ",
    length(fs)*length(base_by_dz), " pairs | cores=", ncore)

run_exposure <- function(f){
  eid_exp <- ge(f)
  oof <- tryCatch(fread(f, select=c("eid","instance","y_raw","pes_prot_z"))[instance==0,
                    .(eid, E=as.numeric(y_raw), PES=as.numeric(pes_prot_z))], error=function(e) NULL)
  if(is.null(oof) || !nrow(oof)) return(NULL)
  out <- vector("list", length(base_by_dz)); k <- 0L
  for(nm in names(base_by_dz)){
    d <- merge(base_by_dz[[nm]], cov0, by="eid")          # covariates (shared) joined per pair
    d <- merge(d, oof, by="eid"); d <- d[is.finite(E) & is.finite(PES)]
    if(nrow(d)<2000 || sum(d$status)<MINEVT) next
    for (fc in names(d)) if (is.factor(d[[fc]])) set(d, j=fc, value=droplevels(d[[fc]]))
    sdp <- sd(d$PES); if(!is.finite(sdp) || sdp==0) next
    d[, PES := (PES-mean(PES))/sdp]                 # HR per 1 SD of PES in-sample
    fit <- function(extra) tryCatch(coxph(as.formula(paste0("Surv(time,status) ~ ", basef, extra)), d),
                                    error=function(e) NULL)
    m0 <- fit(""); mP <- fit(" + PES"); mE <- fit(" + E"); mEP <- fit(" + E + PES")
    if(is.null(m0) || is.null(mP) || is.null(mE) || is.null(mEP)) next
    aP  <- pes_inf(mP,  m0)          # PES over covariates
    aEP <- pes_inf(mEP, mE)          # PES beyond covariates + E (the headline)
    # apparent C-index SCREEN (in-sample; bootstrap refines candidates downstream)
    lp <- function(m) tryCatch(predict(m, type="lp"), error=function(e) rep(NA_real_, nrow(d)))
    C0 <- appC(d$time,d$status,lp(m0)); CP <- appC(d$time,d$status,lp(mP))
    CE <- appC(d$time,d$status,lp(mE)); CEP<- appC(d$time,d$status,lp(mEP))
    k <- k+1L
    out[[k]] <- data.table(
      exposure_id=eid_exp, disease=nm, n=nrow(d), events=sum(d$status),
      hr_pes=aP$hr,  hr_pes_lo=aP$lo,  hr_pes_hi=aP$hi,  wald_pes=aP$wp,  lrt_pes=aP$lrt,
      hr_pesE=aEP$hr, hr_pesE_lo=aEP$lo, hr_pesE_hi=aEP$hi, wald_pesE=aEP$wp, lrt_pesE=aEP$lrt,
      dC_overcov_app=CP-C0, dC_beyondE_app=CEP-CE)
  }
  if(!k) return(NULL)
  rbindlist(out[seq_len(k)], fill=TRUE)
}

res <- parallel::mclapply(fs, run_exposure, mc.cores=ncore)
ok  <- vapply(res, is.data.table, logical(1))
if(any(!ok)) msg("WARN: ", sum(!ok), " exposure(s) errored")
R <- rbindlist(res[ok], fill=TRUE)
R[, q_pes  := p.adjust(lrt_pes,  "BH")]
R[, q_pesE := p.adjust(lrt_pesE, "BH")]
OUTF <- if (identical(covarType,"base")) {
  file.path(out_dir,"pes_disease_scale_base_rerun.tsv")   # never clobber the shipped file
} else {
  file.path(out_dir, paste0("pes_disease_scale_", covarType, ".tsv"))
}
fwrite(R, OUTF, sep="\t")
if (identical(covarType,"base")) {
  ref <- tryCatch(fread(file.path(out_dir,"pes_disease_scale.tsv")), error=function(e) NULL)
  if (!is.null(ref)) {
    m <- merge(R[,.(exposure_id,disease,hr_pes,hr_pesE)],
               ref[,.(exposure_id,disease,r1=hr_pes,r2=hr_pesE)], by=c("exposure_id","disease"))
    dm <- max(abs(c(m$hr_pes-m$r1, m$hr_pesE-m$r2)), na.rm=TRUE)
    msg(sprintf("  BASE REPRODUCE CHECK: %d pairs | max|diff| in HR = %.3e  %s",
                nrow(m), dm, if (is.finite(dm) && dm < 1e-8) "EXACT" else "*** MISMATCH ***"))
  }
}
msg("wrote pes_disease_scale.tsv (", nrow(R), " pairs; ",
    uniqueN(R$exposure_id), " exposures x ", uniqueN(R$disease), " diseases)")
msg("PES beyond E: BH-sig pairs (q<0.05) = ", sum(R$q_pesE<0.05, na.rm=TRUE),
    " over ", sum(is.finite(R$q_pesE)), " tested")
