#!/usr/bin/env Rscript
# ============================================================================
# module6_pes_disease_bootstrap.R  (support: bootstrap CIs for HIGHLIGHT pairs)
# ----------------------------------------------------------------------------
# For the handful of exposure->disease pairs we actually highlight in the main
# figure / text, quantify whether the PES C-index gain is significant. The
# C-index difference of two correlated nested Cox models has no clean closed
# form, so we bootstrap: in each of B resamples we refit the models on the
# resample and evaluate the C-index on the OUT-OF-BAG participants (honest, no
# in-sample optimism). We report TWO gains per pair:
#   dC_overE0 = C(cov+PES)   - C(cov)        # PES over covariates (Fig 6d metric)
#   dC_beyond = C(cov+E+PES) - C(cov+E)      # PES beyond covariates + self-report
# A gain is "significant" if its 95% percentile CI excludes 0.
# Writes multipes_disease/pes_disease_bootstrap.tsv.
# ============================================================================
local({ cand <- c(Sys.getenv("HEAP_PATHS_FILE",""),"/n/groups/patel/shakson_ukb/HEAP/workflow/00_paths.R")
  h <- cand[nzchar(cand)&file.exists(cand)][1]; if(!is.na(h)) source(h) })
suppressPackageStartupMessages({ library(data.table); library(survival) })
msg <- function(...) cat(format(Sys.time(),"[%H:%M:%S] "),...,"\n",sep="")
RNGkind("L'Ecuyer-CMRG"); set.seed(2026)
ncore <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset="8"))
B     <- as.integer(Sys.getenv("HEAP_BOOT_B", unset="400"))
od    <- heap_project_output("module6_pes_longitudinal","base")
out_dir <- heap_project_output("module6_pes_longitudinal","multipes_disease")

# highlight pairs: exposure_id, disease grep pattern, readable label.
# If HEAP_BOOT_PAIRS points to a TSV (cols: exposure_id, disease, label), use it
# (the screen-selected candidate set); otherwise fall back to the named defaults.
pairs_file <- Sys.getenv("HEAP_BOOT_PAIRS", unset="")
if (nzchar(pairs_file) && file.exists(pairs_file)) {
  pf <- fread(pairs_file)
  P <- lapply(seq_len(nrow(pf)), function(i) list(pf$exposure_id[i], pf$disease[i], pf$label[i]))
} else {
P <- list(
  list("pack_years_of_smoking_f20161_0_0",                 "j43_first_reported_emphysema",            "Pack-years smoking -> Emphysema"),
  list("smoking_status_f20116_0_0_Current",                "j43_first_reported_emphysema",            "Current smoking -> Emphysema"),
  list("pack_years_of_smoking_f20161_0_0",                 "j44_first_reported_other_chronic",        "Pack-years smoking -> COPD"),
  list("met_minutes_per_week_for_vigorous_activity_f22039_0_0","e11_first_reported_non_insulin",       "Vigorous activity -> Type-2 diabetes"),
  list("alcohol_intake_frequency_f1558_0_0",               "k70_first_reported_alcoholic_liver",      "Alcohol frequency -> Alcoholic liver disease"),
  list("usual_walking_pace_f924_0_0",                      "e66_first_reported_obesity",              "Walking pace -> Obesity (broad control)"),
  list("usual_walking_pace_f924_0_0",                      "e11_first_reported_non_insulin",          "Walking pace -> Type-2 diabetes (broad control)"))
}

msg("Loading HEAP.rds")
heap <- readRDS(heap_loader_rds); dz <- as.data.table(heap$disease$DZ_df)
cvl  <- as.data.table(heap$covars_long); rm(heap); gc()
sexcol<-grep("^sex_f31",names(cvl),value=T)[1]; agecol<-grep("age_when_attended_assessment_centre",names(cvl),value=T)[1]
PCS<-paste0("genetic_principal_components_f22009_0_",1:20); ctr<-"uk_biobank_assessment_centre_f54_0_0"
cov0<-cvl[instance==0,c("eid",agecol,sexcol,"age2","age_sex","age2_sex",ctr,PCS),with=F]
setnames(cov0,c(agecol,sexcol,ctr),c("age0","sex","centre")); numc<-c("age0","sex","age2","age_sex","age2_sex",PCS)
cov0[,(numc):=lapply(.SD,as.numeric),.SDcols=numc]; cov0[,centre:=factor(centre)]
rm(cvl); gc()   # free the full covars_long (~2.8GB) so per-pair mclapply forks stay small
basef<-paste("age0+sex+age2+age_sex+age2_sex+centre",paste(PCS,collapse="+"),sep="+")
for(c in c("recode_age_of_death_0_0","age_of_removal_0_0","age_of_lastfollowup","recode_age_of_assessment_0_0")) if(c%in%names(dz)) dz[,(c):=as.numeric(get(c))]
dz[,censor_age:=pmin(recode_age_of_death_0_0,age_of_removal_0_0,age_of_lastfollowup,na.rm=T)]

cidx<-function(t,s,lp){ok<-is.finite(t)&is.finite(s)&is.finite(lp); if(sum(ok)<40||sum(s[ok])<8)return(NA_real_)
  tryCatch(concordance(Surv(t[ok],s[ok])~lp[ok],reverse=T)$concordance,error=function(e)NA_real_)}

build<-function(eid, dzpat){
  ac<-grep(dzpat,names(dz),value=T,ignore.case=T); ac<-ac[grepl("^age_",ac)][1]; if(is.na(ac))return(NULL)
  f<-file.path(od,paste0("PESlong_base_",eid,"_TrainOOF.tsv"))
  oof<-tryCatch(fread(f,select=c("eid","instance","y_raw","pes_prot_z"))[instance==0,.(eid,E=as.numeric(y_raw),PES=as.numeric(pes_prot_z))],error=function(e)NULL); if(is.null(oof))return(NULL)
  b<-merge(dz[,.(eid,event_age=as.numeric(get(ac)),asmt=recode_age_of_assessment_0_0,censor_age)],cov0,by="eid")
  b<-b[is.finite(asmt)&is.finite(censor_age)&censor_age>asmt]; b<-b[is.na(event_age)|event_age>asmt]
  b[,status:=as.integer(!is.na(event_age)&event_age<=censor_age)]; b[,time:=fifelse(status==1,event_age,censor_age)-asmt]; b<-b[time>0]
  d<-merge(b,oof,by="eid"); d<-d[is.finite(E)&is.finite(PES)]; d[,PES:=as.numeric(scale(PES))]; d
}

run_pair<-function(p){
  eid<-p[[1]]; d<-build(eid,p[[2]]); if(is.null(d)||nrow(d)<2000||sum(d$status)<100) return(NULL)
  n<-nrow(d); fml<-function(x) as.formula(paste0("Surv(time,status)~",basef,x))
  boot<-parallel::mclapply(1:B, function(b){
    set.seed(2026+b)
    ix<-sample.int(n,n,replace=TRUE); oob<-setdiff(seq_len(n),unique(ix))
    if(length(oob)<200) return(NULL); tr<-d[ix]; te<-d[oob]
    fit<-function(x) tryCatch(coxph(fml(x),tr),error=function(e)NULL)
    m0<-fit(""); mP<-fit("+PES"); mE<-fit("+E"); mEP<-fit("+E+PES")
    if(is.null(m0)||is.null(mP)||is.null(mE)||is.null(mEP)) return(NULL)
    pr<-function(m) tryCatch(predict(m,newdata=te,type="lp"),error=function(e)rep(NA,nrow(te)))
    C0<-cidx(te$time,te$status,pr(m0)); CP<-cidx(te$time,te$status,pr(mP))
    CE<-cidx(te$time,te$status,pr(mE)); CEP<-cidx(te$time,te$status,pr(mEP))
    # absolute held-out C-index at each nested rung (covariates -> +E -> +E+PES)
    # plus the two PES increments, all on the same OOB sample
    c(C_cov=C0, C_covP=CP, C_covE=CE, C_covEP=CEP, dC_overE0=CP-C0, dC_beyond=CEP-CE)
  }, mc.cores=ncore)
  boot<-do.call(rbind, boot[!vapply(boot,is.null,logical(1))])
  ci<-function(v){v<-v[is.finite(v)]; c(mean(v), quantile(v,.025), quantile(v,.975), mean(v>0))}
  a<-ci(boot[,"dC_overE0"]); e<-ci(boot[,"dC_beyond"])
  cc<-ci(boot[,"C_cov"]); cp<-ci(boot[,"C_covP"]); ce<-ci(boot[,"C_covE"]); cep<-ci(boot[,"C_covEP"])
  data.table(label=p[[3]], exposure_id=eid, disease=p[[2]], n=n, events=sum(d$status), nboot=nrow(boot),
    C_cov=cc[1], C_cov_lo=cc[2], C_cov_hi=cc[3],
    C_covP=cp[1], C_covP_lo=cp[2], C_covP_hi=cp[3],
    C_covE=ce[1], C_covE_lo=ce[2], C_covE_hi=ce[3],
    C_covEP=cep[1], C_covEP_lo=cep[2], C_covEP_hi=cep[3],
    dC_overCov=a[1], overCov_lo=a[2], overCov_hi=a[3], overCov_pgt0=a[4],
    dC_beyondE=e[1], beyondE_lo=e[2], beyondE_hi=e[3], beyondE_pgt0=e[4])
}
msg("bootstrapping ", length(P), " highlight pairs (B=",B,", cores=",ncore,")")
R<-rbindlist(lapply(P, function(p){ msg("  ",p[[3]]); run_pair(p) }), fill=TRUE)
outfile <- Sys.getenv("HEAP_BOOT_OUT", unset="pes_disease_bootstrap.tsv")
fwrite(R, file.path(out_dir,outfile), sep="\t")
msg("wrote ",outfile," (",nrow(R)," pairs)")
print(R[, .(label, events,
  overCov=sprintf("%.3f [%.3f,%.3f]",dC_overCov,overCov_lo,overCov_hi),
  beyondE=sprintf("%.3f [%.3f,%.3f]",dC_beyondE,beyondE_lo,beyondE_hi))])
