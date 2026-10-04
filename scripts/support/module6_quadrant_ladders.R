#!/usr/bin/env Rscript
# ============================================================================
# module6_quadrant_ladders.R   (support analysis for Q3 read x disease quadrant)
# ----------------------------------------------------------------------------
# For a curated slate of (exposure -> disease) pairs, compute the HONEST held-out
# (70/30 split, OOF i0) Cox C-index ladder:
#   M0 covariates(BASE: age,sex,age^2,age*sex,age^2*sex,centre,20 PCs) -> M1 +PES -> M2 +E -> M3 +E+PES.
# This is the clean, well-powered version (not the leaky apparent C-index): the
# PES is already out-of-fold and we still evaluate on a held-out 30% test split.
# Used to place exemplars in the read x disease quadrant figure (Score A/B/C/D).
# Writes multipes_disease/quadrant_ladders.tsv.
# ============================================================================
local({ cand <- c(Sys.getenv("HEAP_PATHS_FILE",""),"/n/groups/patel/shakson_ukb/HEAP/workflow/00_paths.R")
  h <- cand[nzchar(cand)&file.exists(cand)][1]; if(!is.na(h)) source(h) })
suppressPackageStartupMessages({ library(data.table); library(survival) })
msg <- function(...) cat(format(Sys.time(),"[%H:%M:%S] "),...,"\n",sep="")
set.seed(42)
# covarType from argv. Outputs are SUFFIXED by spec so an alternative
# specification can never overwrite the base results that Fig 6 panel d reads;
# `base` deliberately keeps the original unsuffixed filenames so nothing existing
# moves. Verify a base re-run is byte-identical before running any other spec.
args      <- commandArgs(trailingOnly = TRUE)
covarType <- if (length(args) && nzchar(args[1])) args[1] else "base"
SFX       <- if (covarType == "base") "" else paste0("_", covarType)
od <- heap_project_output("module6_pes_longitudinal", covarType)
out_dir <- heap_project_output("module6_pes_longitudinal","multipes_disease"); dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)

# ---- curated slate: exposure_id | disease regex | disease label | quad guess ----
SL <- fread(text='exposure_id\tdz_re\tdz_lab\tquad
current_tobacco_smoking_f1239_0_0_Yes._on_most_or_all_days\tj44_first_reported_other_chronic_obstructive\tCOPD\tC
current_tobacco_smoking_f1239_0_0_Yes._on_most_or_all_days\tj43_first_reported_emphysema\tEmphysema\tC
alcohol_intake_frequency_f1558_0_0\tf10_first_reported_mental_and_behavioural_disorders_due_to_use_of_alcohol\tAlcohol-use disorder\tC
number_of_days_week_of_vigorous_physical_activity_10_plus_minutes_f904_0_0\te11_first_reported_non_insulin\tType-2 diabetes\tB
summed_days_activity_f22033_0_0\te11_first_reported_non_insulin\tType-2 diabetes\tB
coffee_intake_f1498_0_0\te11_first_reported_non_insulin\tType-2 diabetes\tA
oily_fish_intake_f1329_0_0\ti25_first_reported_chronic_ischaemic\tIschaemic heart\tA
oily_fish_intake_f1329_0_0\te78_first_reported_disorders_of_lipoprotein\tLipid disorder\tA
nap_during_day_f1190_0_0\tf32_first_reported_depressive\tDepression\tD
usual_walking_pace_f924_0_0\te11_first_reported_non_insulin\tType-2 diabetes\tD
pm2_5_mean\te11_first_reported_non_insulin\tType-2 diabetes\tD', sep='\t', header=TRUE)

msg("Loading HEAP.rds (disease + covars)")
heap <- readRDS(heap_loader_rds)
dz <- as.data.table(heap$disease$DZ_df)
cvl <- as.data.table(heap$covars_long); rm(heap); gc()
sexcol <- grep("^sex_f31",names(cvl),value=TRUE)[1]; agecol <- grep("age_when_attended_assessment_centre",names(cvl),value=TRUE)[1]
# BASE covariate set (matches frozen-risk CovarSpec$base / module6_quadrant_scan.R)
PCS <- paste0("genetic_principal_components_f22009_0_", 1:20); ctr <- "uk_biobank_assessment_centre_f54_0_0"
cov0 <- cvl[instance==0, c("eid", agecol, sexcol, "age2","age_sex","age2_sex", ctr, PCS), with=FALSE]
setnames(cov0, c(agecol, sexcol, ctr), c("age0","sex","centre"))
numc <- c("age0","sex","age2","age_sex","age2_sex", PCS); cov0[, (numc) := lapply(.SD, as.numeric), .SDcols=numc]
cov0[, centre := factor(centre)]
basef <- paste("age0 + sex + age2 + age_sex + age2_sex + centre", paste(PCS, collapse=" + "), sep=" + ")
for(c in c("recode_age_of_death_0_0","age_of_removal_0_0","age_of_lastfollowup","recode_age_of_assessment_0_0"))
  if(c %in% names(dz)) dz[,(c):=as.numeric(get(c))]
dz[, censor_age := pmin(recode_age_of_death_0_0, age_of_removal_0_0, age_of_lastfollowup, na.rm=TRUE)]

ci <- function(time,status,lp){ ok<-is.finite(time)&is.finite(status)&is.finite(lp)
  if(sum(ok)<50||sum(status[ok])<10) return(NA_real_)
  tryCatch(concordance(Surv(time[ok],status[ok])~lp[ok], reverse=TRUE)$concordance, error=function(e) NA_real_) }

oof_cache <- list()
get_oof <- function(eid_exp){
  if(is.null(oof_cache[[eid_exp]])){
    f <- file.path(od, sprintf("PESlong_%s_%s_TrainOOF.tsv", covarType, eid_exp))
    d <- fread(f, select=c("eid","instance","y_raw","pes_prot_z"))[instance==0, .(eid, E=as.numeric(y_raw), PES=as.numeric(pes_prot_z))]
    oof_cache[[eid_exp]] <<- d
  }
  oof_cache[[eid_exp]]
}

res <- list()
for(i in seq_len(nrow(SL))){
  acol <- grep(SL$dz_re[i], names(dz), value=TRUE, ignore.case=TRUE); acol <- acol[grepl("^age_",acol)][1]
  if(is.na(acol)){ msg("no disease col for ", SL$dz_lab[i]); next }
  oof <- get_oof(SL$exposure_id[i])
  d <- merge(dz[, .(eid, event_age=as.numeric(get(acol)), asmt=recode_age_of_assessment_0_0, censor_age)], cov0, by="eid")
  d <- merge(d, oof, by="eid")
  d <- d[is.finite(asmt)&is.finite(censor_age)&censor_age>asmt&is.finite(E)&is.finite(PES)]
  d <- d[is.na(event_age)|event_age>asmt]
  d[, status := as.integer(!is.na(event_age)&event_age<=censor_age)]
  d[, time := fifelse(status==1, event_age, censor_age)-asmt]; d <- d[time>0]
  n <- nrow(d); tr <- sample(n, floor(0.7*n)); te <- setdiff(seq_len(n), tr); dtr<-d[tr]; dte<-d[te]
  lp <- function(extra) tryCatch(predict(coxph(as.formula(paste0("Surv(time,status) ~ ", basef, extra)), dtr), dte), error=function(e) rep(NA_real_, length(te)))
  C0 <- ci(d$time[te],d$status[te], lp(""))
  C1 <- ci(d$time[te],d$status[te], lp(" + PES"))
  C2 <- ci(d$time[te],d$status[te], lp(" + E"))
  C3 <- ci(d$time[te],d$status[te], lp(" + E + PES"))
  res[[i]] <- data.table(exposure_id=SL$exposure_id[i], disease=SL$dz_lab[i], quad_guess=SL$quad[i],
     n=n, events=sum(d$status), C0=C0, C1_PES=C1, C2_E=C2, C3_both=C3, dC_pes=C1-C0, dC_e=C2-C0)
  msg(sprintf("%-26s -> %-22s n=%d ev=%d | C0=%.3f +PES=%.3f +E=%.3f", substr(SL$exposure_id[i],1,26), SL$dz_lab[i], n, sum(d$status), C0, C1, C2))
}
out <- rbindlist(res, fill=TRUE)
fwrite(out, file.path(out_dir, sprintf("quadrant_ladders%s.tsv", SFX)), sep="\t")
msg("wrote ", file.path(out_dir, sprintf("quadrant_ladders%s.tsv", SFX)))
print(out[, .(exposure_id=substr(exposure_id,1,24), disease, quad_guess, events, C0=round(C0,3), C1_PES=round(C1_PES,3), C2_E=round(C2_E,3), dC_pes=round(dC_pes,3), dC_e=round(dC_e,3))])
