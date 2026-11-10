#!/usr/bin/env Rscript
# ============================================================================
# module6_pes_disease_dc_holdout.R  (support: HELD-OUT dC screen, all pairs)
# ----------------------------------------------------------------------------
# Apparent (in-sample) C-index gain is biased toward 0 here: the ~46-covariate
# baseline (age/sex terms + 20 PCs + centre) overfits in-sample, inflating the
# baseline C and leaving no apparent room for the PES. Only an OUT-OF-SAMPLE
# C-index reveals the true PES gain (held-out and bootstrap agree; apparent does
# not). So for every exposure x disease we fit on a 70% split and evaluate the
# C-index on the held-out 30%:
#   dC_overcov_ho = C(cov+PES)   - C(cov)        # PES over covariates (Fig 6d)
#   dC_beyondE_ho = C(cov+E+PES) - C(cov+E)      # PES beyond cov + self-report
#   dC_e_ho       = C(cov+E)     - C(cov)        # self-report over covariates
# This is a SCREEN (single split, noisy per cell) used to rank candidates; the
# highlighted pairs get an out-of-bag bootstrap CI downstream. HR/LRT come from
# the full-cohort scan (pes_disease_scale.tsv) and are merged separately.
# Writes multipes_disease/pes_disease_dc_holdout.tsv.
# ============================================================================
local({ cand <- c(Sys.getenv("HEAP_PATHS_FILE",""),"/n/groups/patel/shakson_ukb/HEAP/workflow/00_paths.R")
  h <- cand[nzchar(cand)&file.exists(cand)][1]; if(!is.na(h)) source(h) })
suppressPackageStartupMessages({ library(data.table); library(survival) })
msg <- function(...) cat(format(Sys.time(),"[%H:%M:%S] "),...,"\n",sep="")
RNGkind("L'Ecuyer-CMRG"); set.seed(99)
ncore  <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset="4"))
MINEVT <- as.integer(Sys.getenv("HEAP_PES_MIN_EVENTS", unset="100"))
# Optional spec argument, resolved on the SAME three axes as
# module6_pes_disease_scale.R. base_exclprev is not a covariate set: it keeps the
# base adjustment, restricts to participants free of major disease at baseline,
# and reads the PES retrained on that subset. Default "" reproduces the original
# base run byte-for-byte (same paths, same seeds, same output name).
argv     <- commandArgs(trailingOnly = TRUE)
covarType<- if (length(argv)) argv[1] else ""
EXCLPREV <- identical(covarType, "base_exclprev")
if (nzchar(covarType) && !EXCLPREV)
  stop("module6_pes_disease_dc_holdout.R: unsupported spec '", covarType,
       "' (only base_exclprev; the held-out screen fixes the adjustment)")
scoreDir <- if (EXCLPREV) "base_exclprev" else "base"
od     <- heap_project_output("module6_pes_longitudinal", scoreDir)
out_dir<- heap_project_output("module6_pes_longitudinal","multipes_disease")
dzset  <- fread(file.path(out_dir,"disease_set_181.tsv"))$DZ_ID
msg("disease set: ", length(dzset), " diseases")

msg("Loading HEAP.rds")
heap <- readRDS(heap_loader_rds); dz <- as.data.table(heap$disease$DZ_df)
cvl  <- as.data.table(heap$covars_long); rm(heap); gc()
sexcol<-grep("^sex_f31",names(cvl),value=T)[1]; agecol<-grep("age_when_attended_assessment_centre",names(cvl),value=T)[1]
PCS<-paste0("genetic_principal_components_f22009_0_",1:20); ctr<-"uk_biobank_assessment_centre_f54_0_0"
cov0<-cvl[instance==0,c("eid",agecol,sexcol,"age2","age_sex","age2_sex",ctr,PCS),with=F]
setnames(cov0,c(agecol,sexcol,ctr),c("age0","sex","centre")); numc<-c("age0","sex","age2","age_sex","age2_sex",PCS)
cov0[,(numc):=lapply(.SD,as.numeric),.SDcols=numc]; cov0[,centre:=factor(centre)]
if (EXCLPREV) {
  # must be taken BEFORE cvl is freed below
  fcol <- "prevalent_major_disease"
  if (!fcol %in% names(cvl)) stop("covars_long has no `", fcol, "` column")
  keep <- cvl[[fcol]]
  ids  <- cvl$eid[cvl$instance == 0 & !is.na(keep) & keep == 0]
  n0 <- nrow(cov0); cov0 <- cov0[eid %in% ids]
  msg("  exclude_prevalent: ", n0, " -> ", nrow(cov0), " participants (healthy at baseline)")
  if (!nrow(cov0)) stop("exclude_prevalent left no participants")
}
rm(cvl); gc()  # free covars_long to keep mclapply forks small
basef<-paste("age0+sex+age2+age_sex+age2_sex+centre",paste(PCS,collapse="+"),sep="+")
for(c in c("recode_age_of_death_0_0","age_of_removal_0_0","age_of_lastfollowup","recode_age_of_assessment_0_0")) if(c%in%names(dz)) dz[,(c):=as.numeric(get(c))]
dz[,censor_age:=pmin(recode_age_of_death_0_0,age_of_removal_0_0,age_of_lastfollowup,na.rm=T)]
base_by_dz<-list()
for(ac in dzset){ if(!ac%in%names(dz)) next
  b<-merge(dz[,.(eid,event_age=as.numeric(get(ac)),asmt=recode_age_of_assessment_0_0,censor_age)],cov0[,.(eid)],by="eid")
  b<-b[is.finite(asmt)&is.finite(censor_age)&censor_age>asmt]; b<-b[is.na(event_age)|event_age>asmt]
  b[,status:=as.integer(!is.na(event_age)&event_age<=censor_age)]; b[,time:=fifelse(status==1,event_age,censor_age)-asmt]; b<-b[time>0]
  if(sum(b$status)<MINEVT) next
  base_by_dz[[ac]]<-b[,.(eid,status,time)] }
msg("disease tables: ", length(base_by_dz))

cidx<-function(t,s,lp){ok<-is.finite(t)&is.finite(s)&is.finite(lp); if(sum(ok)<40||sum(s[ok])<8)return(NA_real_)
  tryCatch(concordance(Surv(t[ok],s[ok])~lp[ok],reverse=T)$concordance,error=function(e)NA_real_)}
fs<-list.files(od,pattern="_TrainOOF\\.tsv$",full.names=TRUE)
ge<-function(f) sub(sprintf("^PESlong_%s_(.*)_TrainOOF\\.tsv$", scoreDir),"\\1",basename(f))
msg(length(fs)," exposures x ",length(base_by_dz)," diseases | cores=",ncore)

run_exposure<-function(f){
  eid<-ge(f)
  oof<-tryCatch(fread(f,select=c("eid","instance","y_raw","pes_prot_z"))[instance==0,.(eid,E=as.numeric(y_raw),PES=as.numeric(pes_prot_z))],error=function(e)NULL)
  if(is.null(oof)||!nrow(oof)) return(NULL)
  out<-vector("list",length(base_by_dz)); k<-0L
  for(nm in names(base_by_dz)){
    d<-merge(base_by_dz[[nm]],cov0,by="eid"); d<-merge(d,oof,by="eid"); d<-d[is.finite(E)&is.finite(PES)]
    if(nrow(d)<2000||sum(d$status)<MINEVT) next
    sdp<-sd(d$PES); if(!is.finite(sdp)||sdp==0) next; d[,PES:=(PES-mean(PES))/sdp]
    n<-nrow(d); set.seed(1000+k); tr<-sample(n,floor(0.7*n)); te<-setdiff(seq_len(n),tr); dtr<-d[tr]
    fit<-function(x) tryCatch(coxph(as.formula(paste0("Surv(time,status)~",basef,x)),dtr),error=function(e)NULL)
    m0<-fit(""); mP<-fit("+PES"); mE<-fit("+E"); mEP<-fit("+E+PES")
    if(is.null(m0)||is.null(mP)||is.null(mE)||is.null(mEP)) next
    pr<-function(m) tryCatch(predict(m,newdata=d[te],type="lp"),error=function(e)rep(NA_real_,length(te)))
    C0<-cidx(d$time[te],d$status[te],pr(m0)); CP<-cidx(d$time[te],d$status[te],pr(mP))
    CE<-cidx(d$time[te],d$status[te],pr(mE)); CEP<-cidx(d$time[te],d$status[te],pr(mEP))
    k<-k+1L
    out[[k]]<-data.table(exposure_id=eid,disease=nm,n=n,events=sum(d$status),
      C0_ho=C0,CP_ho=CP,CE_ho=CE,CEP_ho=CEP,
      dC_overcov_ho=CP-C0, dC_beyondE_ho=CEP-CE, dC_e_ho=CE-C0)
  }
  if(!k) return(NULL); rbindlist(out[seq_len(k)],fill=TRUE)
}
res<-parallel::mclapply(fs,run_exposure,mc.cores=ncore)
ok<-vapply(res,is.data.table,logical(1)); if(any(!ok)) msg("WARN ",sum(!ok)," errored")
R<-rbindlist(res[ok],fill=TRUE)
OUTF <- if (EXCLPREV) "pes_disease_dc_holdout_base_exclprev.tsv" else "pes_disease_dc_holdout.tsv"
fwrite(R,file.path(out_dir,OUTF),sep="\t")   # never clobber the shipped base screen
msg("wrote ",OUTF," (",nrow(R)," pairs)")
msg("held-out dC_beyondE: >0.02 ",sum(R$dC_beyondE_ho>0.02,na.rm=T),
    " | >0.05 ",sum(R$dC_beyondE_ho>0.05,na.rm=T),
    " | dC_overcov >0.05 ",sum(R$dC_overcov_ho>0.05,na.rm=T))
