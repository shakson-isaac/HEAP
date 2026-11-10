#!/usr/bin/env Rscript
# ============================================================================
# module6_pes_disease_dc_holdout_ci.R  (support: CIs on the held-out dC screen)
# ----------------------------------------------------------------------------
# Adds 95% bootstrap intervals to every exposure x disease cell of
# pes_disease_dc_holdout.tsv WITHOUT changing any point estimate.
#
# WHY THIS ESTIMATOR RATHER THAN THE OUT-OF-BAG ONE
# module6_pes_disease_bootstrap.R refits all four Cox models inside every
# resample and evaluates out-of-bag. Measured cost: 7.20 s for the four refits
# vs 0.30 s for the four concordance evaluations -- the fits are 96% of the
# work, so B=400 costs ~3,000 s per pair and ~24,900 core-hours grid-wide.
# Here we instead fit ONCE on the same 70% split the screen already used, then
# resample the 30% TEST SET. That is ~126 s per pair (~1,045 core-hours,
# measured 23.8x cheaper) and, because the fit is untouched, the unresampled
# statistic reproduces the screen's published number exactly.
#
# WHAT THE INTERVAL MEANS -- state this in the legend, the two are NOT the same:
#   this script  -> uncertainty of the C-index GIVEN the fitted model
#                   (evaluation uncertainty only)
#   out-of-bag   -> uncertainty of the whole fit-and-evaluate procedure
# Increments stay valid either way: the test set is resampled ONCE per replicate
# and all four models are re-scored on that same resample, so the pairing
# between nested models -- and hence their correlation -- is preserved.
#
# REPRODUCING THE SPLIT. The screen draws its split with set.seed(1000+k) where
# k is a per-EXPOSURE counter over diseases, so the split depends on a disease's
# position in the loop. 165 exposures x 181 diseases == 29,865 rows in the
# screen output, i.e. no pair was ever skipped, so k is simply the 0-based index
# into names(base_by_dz). That makes the work parallelisable by exposure while
# still reproducing every split. VERIFY=1 checks this against the shipped file.
#
#   Usage:  Rscript module6_pes_disease_dc_holdout_ci.R <exposure_index>
#           HEAP_VERIFY=1 Rscript ... 1      # reproduce-only, no bootstrap
#   Writes: multipes_disease/dc_holdout_ci/<exposure_id>.tsv
# ============================================================================
local({ cand <- c(Sys.getenv("HEAP_PATHS_FILE",""),"/n/groups/patel/shakson_ukb/HEAP/workflow/00_paths.R")
  h <- cand[nzchar(cand)&file.exists(cand)][1]; if(!is.na(h)) source(h) })
suppressPackageStartupMessages({ library(data.table); library(survival) })
msg <- function(...) cat(format(Sys.time(),"[%H:%M:%S] "),...,"\n",sep="")
# MUST match module6_pes_disease_dc_holdout.R line 23. set.seed() draws from
# whichever RNG kind is current, so the same seed under Mersenne-Twister yields a
# completely different 70/30 split -- verified: omitting this put max|diff| at
# 6.1e-01 against the shipped screen instead of 0.
RNGkind("L'Ecuyer-CMRG"); set.seed(99)

args   <- commandArgs(trailingOnly=TRUE)
IDX    <- if (length(args) && nzchar(args[1])) as.integer(args[1]) else 1L
B      <- as.integer(Sys.getenv("HEAP_BOOT_B", unset="400"))
VERIFY <- nzchar(Sys.getenv("HEAP_VERIFY",""))
ncore  <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset="4"))
MINEVT <- as.integer(Sys.getenv("HEAP_PES_MIN_EVENTS", unset="100"))
od     <- heap_project_output("module6_pes_longitudinal","base")
out_dir<- heap_project_output("module6_pes_longitudinal","multipes_disease")
ci_dir <- file.path(out_dir,"dc_holdout_ci"); dir.create(ci_dir, recursive=TRUE, showWarnings=FALSE)
dzset  <- fread(file.path(out_dir,"disease_set_181.tsv"))$DZ_ID

msg("Loading HEAP.rds")
heap <- readRDS(heap_loader_rds); dz <- as.data.table(heap$disease$DZ_df)
cvl  <- as.data.table(heap$covars_long); rm(heap); gc()
sexcol<-grep("^sex_f31",names(cvl),value=T)[1]; agecol<-grep("age_when_attended_assessment_centre",names(cvl),value=T)[1]
PCS<-paste0("genetic_principal_components_f22009_0_",1:20); ctr<-"uk_biobank_assessment_centre_f54_0_0"
cov0<-cvl[instance==0,c("eid",agecol,sexcol,"age2","age_sex","age2_sex",ctr,PCS),with=F]
setnames(cov0,c(agecol,sexcol,ctr),c("age0","sex","centre")); numc<-c("age0","sex","age2","age_sex","age2_sex",PCS)
cov0[,(numc):=lapply(.SD,as.numeric),.SDcols=numc]; cov0[,centre:=factor(centre)]
rm(cvl); gc()
basef<-paste("age0+sex+age2+age_sex+age2_sex+centre",paste(PCS,collapse="+"),sep="+")
for(c in c("recode_age_of_death_0_0","age_of_removal_0_0","age_of_lastfollowup","recode_age_of_assessment_0_0")) if(c%in%names(dz)) dz[,(c):=as.numeric(get(c))]
dz[,censor_age:=pmin(recode_age_of_death_0_0,age_of_removal_0_0,age_of_lastfollowup,na.rm=T)]

# identical construction to the screen, so the analysis samples match row-for-row
base_by_dz<-list()
for(ac in dzset){ if(!ac%in%names(dz)) next
  b<-merge(dz[,.(eid,event_age=as.numeric(get(ac)),asmt=recode_age_of_assessment_0_0,censor_age)],cov0[,.(eid)],by="eid")
  b<-b[is.finite(asmt)&is.finite(censor_age)&censor_age>asmt]; b<-b[is.na(event_age)|event_age>asmt]
  b[,status:=as.integer(!is.na(event_age)&event_age<=censor_age)]; b[,time:=fifelse(status==1,event_age,censor_age)-asmt]; b<-b[time>0]
  if(sum(b$status)<MINEVT) next
  base_by_dz[[ac]]<-b[,.(eid,status,time)] }

cidx<-function(t,s,lp){ok<-is.finite(t)&is.finite(s)&is.finite(lp); if(sum(ok)<40||sum(s[ok])<8)return(NA_real_)
  tryCatch(concordance(Surv(t[ok],s[ok])~lp[ok],reverse=T)$concordance,error=function(e)NA_real_)}

fs<-sort(list.files(od,pattern="_TrainOOF\\.tsv$",full.names=TRUE))
ge<-function(f) sub("^PESlong_base_(.*)_TrainOOF\\.tsv$","\\1",basename(f))
if (IDX < 1 || IDX > length(fs)) stop("exposure index ", IDX, " outside 1..", length(fs))
f   <- fs[IDX]; eid <- ge(f)
msg("exposure ",IDX,"/",length(fs),": ",eid," | ",length(base_by_dz)," diseases | B=",B," | cores=",ncore)

oof <- fread(f,select=c("eid","instance","y_raw","pes_prot_z"))[instance==0,.(eid,E=as.numeric(y_raw),PES=as.numeric(pes_prot_z))]

one_disease <- function(j) {
  nm <- names(base_by_dz)[j]
  d<-merge(base_by_dz[[nm]],cov0,by="eid"); d<-merge(d,oof,by="eid"); d<-d[is.finite(E)&is.finite(PES)]
  if(nrow(d)<2000||sum(d$status)<MINEVT) return(NULL)
  sdp<-sd(d$PES); if(!is.finite(sdp)||sdp==0) return(NULL); d[,PES:=(PES-mean(PES))/sdp]
  # k is the 0-based position in the disease loop -- see header
  n<-nrow(d); set.seed(1000+(j-1L)); tr<-sample(n,floor(0.7*n)); te<-setdiff(seq_len(n),tr); dtr<-d[tr]
  fit<-function(x) tryCatch(coxph(as.formula(paste0("Surv(time,status)~",basef,x)),dtr),error=function(e)NULL)
  m0<-fit(""); mP<-fit("+PES"); mE<-fit("+E"); mEP<-fit("+E+PES")
  if(is.null(m0)||is.null(mP)||is.null(mE)||is.null(mEP)) return(NULL)
  pr<-function(m) tryCatch(predict(m,newdata=d[te],type="lp"),error=function(e)rep(NA_real_,length(te)))
  tt<-d$time[te]; ss<-d$status[te]
  l0<-pr(m0); lP<-pr(mP); lE<-pr(mE); lEP<-pr(mEP)
  C0<-cidx(tt,ss,l0); CP<-cidx(tt,ss,lP); CE<-cidx(tt,ss,lE); CEP<-cidx(tt,ss,lEP)
  o <- data.table(exposure_id=eid, disease=nm, n=n, events=sum(d$status),
                  C0_ho=C0, CP_ho=CP, CE_ho=CE, CEP_ho=CEP,
                  dC_overcov_ho=CP-C0, dC_beyondE_ho=CEP-CE, dC_e_ho=CE-C0)
  if (VERIFY) return(o)
  # --- bootstrap the TEST SET only; the fit is fixed --------------------------
  # One resample per replicate, re-scored by all four models, so nested models
  # stay paired and the increment interval keeps their correlation.
  m <- length(te); set.seed(7000L + j)
  bC0<-bCP<-bCE<-bCEP<-numeric(B)
  for (b in seq_len(B)) {
    ix <- sample.int(m, m, replace=TRUE)
    t2<-tt[ix]; s2<-ss[ix]
    bC0[b]<-cidx(t2,s2,l0[ix]); bCP[b]<-cidx(t2,s2,lP[ix])
    bCE[b]<-cidx(t2,s2,lE[ix]); bCEP[b]<-cidx(t2,s2,lEP[ix])
  }
  q<-function(v){v<-v[is.finite(v)]; if(length(v)>10) quantile(v,c(.025,.975),names=FALSE) else c(NA_real_,NA_real_)}
  qo<-q(bCP-bC0); qb<-q(bCEP-bCE)
  o[, `:=`(nboot_ho = B,
           C0_ho_lo=q(bC0)[1],  C0_ho_hi=q(bC0)[2],
           CP_ho_lo=q(bCP)[1],  CP_ho_hi=q(bCP)[2],
           CE_ho_lo=q(bCE)[1],  CE_ho_hi=q(bCE)[2],
           CEP_ho_lo=q(bCEP)[1],CEP_ho_hi=q(bCEP)[2],
           dC_overcov_ho_lo=qo[1], dC_overcov_ho_hi=qo[2],
           dC_overcov_ho_pgt0=mean((bCP-bC0)>0, na.rm=TRUE),
           dC_beyondE_ho_lo=qb[1], dC_beyondE_ho_hi=qb[2],
           dC_beyondE_ho_pgt0=mean((bCEP-bCE)>0, na.rm=TRUE))]
  o[]
}

res <- parallel::mclapply(seq_along(base_by_dz), one_disease, mc.cores=ncore)
R <- rbindlist(res[vapply(res,is.data.table,logical(1))], fill=TRUE)
msg("  ",nrow(R)," pairs")

# --- reproduce check against the shipped screen ------------------------------
ref <- fread(file.path(out_dir,"pes_disease_dc_holdout.tsv"))[exposure_id==eid]
if (nrow(ref)) {
  m <- merge(R[,.(disease,C0_ho,CP_ho,CE_ho,CEP_ho)],
             ref[,.(disease,r0=C0_ho,rP=CP_ho,rE=CE_ho,rEP=CEP_ho)], by="disease")
  dmax <- max(abs(c(m$C0_ho-m$r0, m$CP_ho-m$rP, m$CE_ho-m$rE, m$CEP_ho-m$rEP)), na.rm=TRUE)
  exact <- is.finite(dmax) && dmax < 1e-9
  msg(sprintf("  REPRODUCE CHECK: %d matched pairs | max|diff| = %.3e  %s",
              nrow(m), dmax, if (exact) "EXACT" else "*** MISMATCH ***"))
  # Fatal, so this can gate the array without a human reading the log. A mismatch
  # means the 70/30 split was not reproduced, and intervals would then be attached
  # to point estimates from a DIFFERENT split -- plausible-looking and wrong.
  if (!exact) stop("split not reproduced (max|diff| = ", signif(dmax,3), "); refusing to continue")
} else msg("  ! no reference rows for this exposure")

if (!VERIFY) { fwrite(R, file.path(ci_dir, paste0(eid,".tsv")), sep="\t"); msg("  wrote ",file.path(ci_dir,paste0(eid,".tsv"))) }
