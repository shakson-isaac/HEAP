#!/usr/bin/env Rscript
# ============================================================================
# module6_quadrant_scan.R  (support: ALL-exposure read x disease quadrant)
# ----------------------------------------------------------------------------
# For every exposure, compute the HONEST held-out (70/30, OOF i0) Cox C-index
# gain of its PES over a CURATED panel of well-powered, organ-diverse diseases:
#   M0 covariates(BASE: age,sex,age^2,age*sex,age^2*sex,centre,20 PCs) -> M1 +PES(pes_prot_z) -> M2 +E(y_raw).
# Per exposure we keep the BEST disease (max held-out PES gain) so the quadrant
# y-axis is robust (no rare-disease apparent-C flukes). x-axis = panel-b read.
# Writes multipes_disease/quadrant_scan.tsv (per exposure x disease) +
#        multipes_disease/quadrant_scan_best.tsv (per exposure best disease).
# ============================================================================
local({ cand <- c(Sys.getenv("HEAP_PATHS_FILE",""),"/n/groups/patel/shakson_ukb/HEAP/workflow/00_paths.R")
  h <- cand[nzchar(cand)&file.exists(cand)][1]; if(!is.na(h)) source(h) })
suppressPackageStartupMessages({ library(data.table); library(survival) })
msg <- function(...) cat(format(Sys.time(),"[%H:%M:%S] "),...,"\n",sep=""); set.seed(42)
# covarType from argv. Outputs are SUFFIXED by spec so an alternative
# specification can never overwrite the base results that Fig 6 panel d reads;
# `base` deliberately keeps the original unsuffixed filenames so nothing existing
# moves. Verify a base re-run is byte-identical before running any other spec.
args      <- commandArgs(trailingOnly = TRUE)
covarType <- if (length(args) && nzchar(args[1])) args[1] else "base"
SFX       <- if (covarType == "base") "" else paste0("_", covarType)
od <- heap_project_output("module6_pes_longitudinal", covarType)
out_dir <- heap_project_output("module6_pes_longitudinal","multipes_disease")

DZ <- c("Type-2 diabetes"="e11_first_reported_non_insulin","Ischaemic heart disease"="i25_first_reported_chronic_ischaemic",
  "Hypertension"="i10_first_reported_essential","Heart failure"="i50_first_reported_heart_failure",
  "Atrial fibrillation"="i48_first_reported","Stroke"="i63_first_reported_cerebral_infarction",
  "COPD"="j44_first_reported_other_chronic_obstructive","Emphysema"="j43_first_reported_emphysema",
  "Asthma"="j45_first_reported_asthma","Alcoholic liver disease"="k70_first_reported_alcoholic_liver",
  "Chronic kidney disease"="n18_first_reported_chronic_renal","Obesity"="e66_first_reported_obesity",
  "Lipid disorder"="e78_first_reported_disorders_of_lipoprotein","Depression"="f32_first_reported_depressive",
  "Lung cancer"="c34_first_reported","Knee osteoarthritis"="m17_first_reported_gonarthrosis",
  "Dementia"="g30_first_reported_alzheimer")

msg("Loading HEAP.rds")
heap <- readRDS(heap_loader_rds); dz <- as.data.table(heap$disease$DZ_df); cvl <- as.data.table(heap$covars_long); rm(heap); gc()
sexcol<-grep("^sex_f31",names(cvl),value=TRUE)[1]; agecol<-grep("age_when_attended_assessment_centre",names(cvl),value=TRUE)[1]
# BASE covariate set (config/covariates/covariate_sets.yml == frozen-risk CovarSpec$base):
# age, sex, age^2, age*sex, age^2*sex, assessment centre (factor), 20 genetic PCs.
PCS <- paste0("genetic_principal_components_f22009_0_", 1:20); ctr <- "uk_biobank_assessment_centre_f54_0_0"
cov0 <- cvl[instance==0, c("eid", agecol, sexcol, "age2","age_sex","age2_sex", ctr, PCS), with=FALSE]
setnames(cov0, c(agecol, sexcol, ctr), c("age0","sex","centre"))
numc <- c("age0","sex","age2","age_sex","age2_sex", PCS); cov0[, (numc) := lapply(.SD, as.numeric), .SDcols=numc]
cov0[, centre := factor(centre)]
basef <- paste("age0 + sex + age2 + age_sex + age2_sex + centre", paste(PCS, collapse=" + "), sep=" + ")
for(c in c("recode_age_of_death_0_0","age_of_removal_0_0","age_of_lastfollowup","recode_age_of_assessment_0_0"))
  if(c %in% names(dz)) dz[,(c):=as.numeric(get(c))]
dz[, censor_age := pmin(recode_age_of_death_0_0, age_of_removal_0_0, age_of_lastfollowup, na.rm=TRUE)]

# precompute per-disease base table (eid, covars, status, time)
base_by_dz <- list()
for(nm in names(DZ)){
  ac <- grep(DZ[[nm]], names(dz), value=TRUE, ignore.case=TRUE); ac <- ac[grepl("^age_",ac)][1]
  if(is.na(ac)){ msg("MISSING disease col: ",nm); next }
  b <- merge(dz[, .(eid, event_age=as.numeric(get(ac)), asmt=recode_age_of_assessment_0_0, censor_age)], cov0, by="eid")
  b <- b[is.finite(asmt)&is.finite(censor_age)&censor_age>asmt]; b <- b[is.na(event_age)|event_age>asmt]
  b[, status := as.integer(!is.na(event_age)&event_age<=censor_age)]
  b[, time := fifelse(status==1, event_age, censor_age)-asmt]; b <- b[time>0]
  base_by_dz[[nm]] <- b[, c("eid","status","time", numc, "centre"), with=FALSE]
}
msg("disease base tables: ", length(base_by_dz))

ci <- function(time,status,lp){ ok<-is.finite(time)&is.finite(status)&is.finite(lp)
  if(sum(ok)<50||sum(status[ok])<10) return(NA_real_)
  tryCatch(concordance(Surv(time[ok],status[ok])~lp[ok], reverse=TRUE)$concordance, error=function(e) NA_real_) }

fs <- list.files(od, pattern="_TrainOOF\\.tsv$", full.names=TRUE)
ge <- function(f) sub(sprintf("^PESlong_%s_(.*)_TrainOOF\\.tsv$", covarType), "\\1", basename(f))
msg(length(fs)," exposures x ",length(base_by_dz)," diseases")
# base model is ~45 covariates -> parallelise over exposures (mclapply, fork shares base_by_dz)
RNGkind("L'Ecuyer-CMRG"); set.seed(42)
run_exposure <- function(f){
  eid_exp <- ge(f)
  oof <- tryCatch(fread(f, select=c("eid","instance","y_raw","pes_prot_z"))[instance==0, .(eid, E=as.numeric(y_raw), PES=as.numeric(pes_prot_z))], error=function(e) NULL)
  if(is.null(oof)||!nrow(oof)) return(NULL)
  out <- list()
  for(nm in names(base_by_dz)){
    d <- merge(base_by_dz[[nm]], oof, by="eid"); d <- d[is.finite(E)&is.finite(PES)]
    if(nrow(d)<2000 || sum(d$status)<60) next
    n<-nrow(d); tr<-sample(n,floor(0.7*n)); te<-setdiff(seq_len(n),tr); dtr<-d[tr]
    lp<-function(extra) tryCatch(predict(coxph(as.formula(paste0("Surv(time,status) ~ ", basef, extra)), dtr), d[te]), error=function(e) rep(NA_real_, length(te)))
    C0<-ci(d$time[te],d$status[te],lp(""))
    C1<-ci(d$time[te],d$status[te],lp(" + PES"))
    C2<-ci(d$time[te],d$status[te],lp(" + E"))
    out[[nm]]<-data.table(exposure_id=eid_exp, disease=nm, events=sum(d$status), C0=C0, C1_PES=C1, C2_E=C2, dC_pes=C1-C0, dC_e=C2-C0)
  }
  rbindlist(out, fill=TRUE)
}
res <- parallel::mclapply(fs, run_exposure, mc.cores=4L)
ok <- vapply(res, is.data.table, logical(1)); if(any(!ok)) msg("WARN: ", sum(!ok), " exposure(s) errored in parallel")
R <- rbindlist(res[ok], fill=TRUE)
fwrite(R, file.path(out_dir, sprintf("quadrant_scan%s.tsv", SFX)), sep="\t")
best <- R[is.finite(dC_pes)][order(exposure_id,-dC_pes)][, .SD[1], by=exposure_id]
fwrite(best, file.path(out_dir, sprintf("quadrant_scan_best%s.tsv", SFX)), sep="\t")
msg("wrote quadrant_scan.tsv (",nrow(R),") + quadrant_scan_best.tsv (",nrow(best),")")
