#!/usr/bin/env Rscript
# ============================================================================
# module3_disease_specificity.R  (pleiotropic vs disease-SPECIFIC mediators)
# ----------------------------------------------------------------------------
# For a given disease, which mediating proteins are PLEIOTROPIC (mediate many
# diseases = shared reporters) vs DISEASE-SPECIFIC (mediate few = candidate causal
# intermediaries, the kind MR validates, e.g. ASGR1)? Per (protein,disease) FDR-sig
# dominant link (n>=100): dominant exposure, NIE, protein disease-PLEIOTROPY (#
# diseases it mediates), and PROPORTION MEDIATED (PM = aggregate-exposome NIE/TE).
# Focus on the cardiometabolic diseases; flag ASGR1. Output: disease_mediators.tsv.
# ============================================================================
local({
  cm <- c(file.path(getwd(),"scripts","visualizations","common"),
          "/n/groups/patel/shakson_ukb/HEAP/scripts/visualizations/common")
  cm <- cm[dir.exists(cm)][1]
  for (f in c("figure_paths","load_heap_results","plot_theme","label_helpers")) source(file.path(cm,paste0(f,".R")))
})
suppressPackageStartupMessages(library(data.table))
OUT <- file.path(Sys.getenv("HEAP_PROJECT_ROOT", "/n/groups/patel/IGLOO/UKB/HEAP"), "figures", "exploratory/module3")
a <- commandArgs(trailingOnly=TRUE); covarType <- if(length(a)>=1) a[1] else "base"; family <- if(length(a)>=2) a[2] else "lasso"

# SPECIFICATION -> (experiment, covariate directory). These are not the same
# thing and for one specification they disagree: base_exclprev is a SAMPLE
# specification, so its covariate_set stays `base` and its runs live under
# .../M3_base_exclprev_lasso_partitioned/base/lasso/..., not under a
# `base_exclprev` directory. load_module3_results() resolves a path from
# covarType alone, so calling it with "base_exclprev" looks for a folder that
# was never created and the run dies after the summariser has already succeeded.
# Passing `experiment` explicitly is what the loader provides for exactly this.
SPEC_DIRS <- list(
  base          = list(exp="M3_base_lasso",          cov="base"),
  base_bmi      = list(exp="M3_base_bmi_lasso",      cov="base_bmi"),
  base_clinical = list(exp="M3_base_clinical_lasso", cov="base_clinical"),
  base_draw     = list(exp="M3_base_draw_lasso",     cov="base_draw"),
  base_exclprev = list(exp="M3_base_exclprev_lasso", cov="base")
)
if (!covarType %in% names(SPEC_DIRS))
  stop("unknown specification '", covarType, "'. One of: ", paste(names(SPEC_DIRS), collapse=", "))
SD <- SPEC_DIRS[[covarType]]
NMIN <- 100
sel <- c("protID","DZ_ID","predictor","predictor_class","effect_type","effect_logHR","effect_HR",
         "delta_p","instrument_present","n","n_cases","protein_HR","protein_p")
md <- load_module3_results(covarType=SD$cov, family=family, mode="partitioned_categories",
                           experiment=paste0(SD$exp, "_partitioned"), select=sel)
md <- md[!(predictor_class %in% c("genetic_cis","genetic_trans") & instrument_present==FALSE)]
pm <- heap_proportion_mediated(md); pm[, category := heap_md_category(predictor)]
e <- pm[!is.na(category) & is.finite(NIE_logHR) & n_cases>=NMIN]
e[, aa := abs(NIE_logHR)][, sig := is.finite(NIE_q) & NIE_q<0.05 & pm_consistent==TRUE]
ple <- e[sig==TRUE, .(pleiotropy=uniqueN(DZ_ID)), by=protID]                 # disease-breadth
dom <- e[order(protID,DZ_ID,-aa)][, .(dom_cat=category[1], dom_NIE=NIE_HR[1], dom_q=NIE_q[1],
  dom_sig=sig[1], n_sig_cat=sum(sig), n_cases=n_cases[1], protein_HR=protein_HR[1]), by=.(protID,DZ_ID)]
# proportion mediated (aggregate exposome) from primary_total
pri <- load_module3_results(covarType=SD$cov, family=family, mode="primary_total",
                            experiment=paste0(SD$exp, "_primary"))
pmt <- heap_proportion_mediated(pri); pmt[, driver := heap_md_driver_component(predictor_class)]
PMx <- pmt[driver=="Exposome" & is.finite(pm), .(protID, DZ_ID, PM=pm)]
j <- merge(dom[dom_sig==TRUE], ple, by="protID")
j <- merge(j, PMx, by=c("protID","DZ_ID"), all.x=TRUE)
j[, aa := abs(log(dom_NIE))]
sysmap <- function(dz){ L<-toupper(sub("^age_([a-zA-Z]).*","\\1",dz))
  m<-c(E="Endocrine/metabolic",I="Circulatory",N="Renal/GU",J="Respiratory",K="Digestive",
       D="Blood/immune",M="Musculoskeletal",F="Mental/neuro"); o<-m[L]; o[is.na(o)]<-"Other"; unname(o) }
j[, disease := heap_pretty_disease(DZ_ID)][, system := sysmap(DZ_ID)]
j[, specific := pleiotropy<=3]                                              # disease-specific protein
# Non-base specifications write their own file. disease_mediators.tsv is what
# fig_mediation_pleiotropy.R reads for the printed panel (325 disease-specific,
# 303 pleiotropic hubs), and this script already accepted a covarType argument
# while always writing that one name -- so running it for a sensitivity silently
# replaced the base result with the sensitivity's. It did exactly that once.
DEST <- if (identical(covarType, "base")) "disease_mediators.tsv" else
  sprintf("disease_mediators_%s.tsv", covarType)
fwrite(j[order(DZ_ID,pleiotropy,-aa)], file.path(OUT, DEST), sep="\t")

cat(sprintf("pleiotropy: median=%d  specific(<=3 dz)=%d proteins  hubs(>=20 dz)=%d\n",
            as.integer(median(ple$pleiotropy)), ple[pleiotropy<=3,.N], ple[pleiotropy>=20,.N]))
cat("\n=== ASGR1 profile (the case example) ===\n")
print(j[protID=="ASGR1"][order(-aa), .(disease=substr(disease,1,28), system, dom_cat, dom_NIE=round(dom_NIE,3),
        PM=round(PM,2), pleiotropy, n_cases)])
CM <- c(e11="type-2 diabetes", e66="obesity", e78="lipoprotein disorder", n18="chronic kidney",
        n17="acute kidney", i25="ischaemic heart dis", i48="atrial fibrillation", i50="heart failure")
for (cd in names(CM)) {
  dz <- j[grepl(paste0("_",cd,"_"), DZ_ID)]
  if (!nrow(dz)) next
  cat(sprintf("\n=== %s : most DISEASE-SPECIFIC mediators (low pleiotropy) ===\n", CM[cd]))
  print(head(dz[order(pleiotropy,-aa), .(protID, dom_cat, dom_NIE=round(dom_NIE,3), PM=round(PM,2),
              pleiotropy, protein_HR=round(protein_HR,2), n_cases)], 8))
}
cat("\nwrote", DEST, "\n")
