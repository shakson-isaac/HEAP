#!/usr/bin/env Rscript
# ============================================================================
# module3_mediation_x_variance.R  (support for fig_mediation_exposure_anchoring)
# ----------------------------------------------------------------------------
# Joins the two arms of the mediation decomposition per (exposure category, protein):
#   exposure->protein arm = Module 1 unique R2 the category explains for the protein
#                           (fold-averaged, score_unique_drop, exposure_categories)
#   protein->disease arm  = Module 3 protein->disease Cox HR
#   result                = Module 3 mediated NIE of the category through the protein
# NIE ~ (exposure->protein) x (protein->disease); this table lets the figure show
# whether a category's top mediators are exposure-anchored (high M1 R2) or
# disease-anchored (high protein HR). Output: exploratory/module3/mediation_x_variance.tsv
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
NMIN <- 300

# --- Module 3: per (category, protein) strongest mediated link ---
sel <- c("protID","DZ_ID","predictor","predictor_class","effect_type","effect_logHR","effect_HR",
         "delta_p","instrument_present","n","n_cases","protein_HR","protein_p")
md <- load_module3_results(covarType=covarType, family=family, mode="partitioned_categories", select=sel)
pm <- heap_proportion_mediated(md); pm[, category := heap_md_category(predictor)]
e <- pm[!is.na(category) & is.finite(NIE_q) & NIE_q<0.05 & pm_consistent==TRUE & n_cases>=NMIN]
e[, absNIE := abs(NIE_logHR)]
m3 <- e[, { k <- which.max(absNIE)
  .(NIE_HR=NIE_HR[k], absNIE=max(absNIE), top_disease=heap_pretty_disease(DZ_ID[k]),
    n_diseases=uniqueN(DZ_ID), n_cases=n_cases[k], protein_dz_HR=protein_HR[k]) }, by=.(category, protID)]

# --- Module 1: fold-averaged unique R2 per (protein, category) ---
m1 <- load_module1_predictive_r2(level="exposure_categories", covarType=covarType,
                                 method=family, experiment=paste0("M1_",covarType,"_",family))
m1a <- m1[, .(m1_r2=mean(r2, na.rm=TRUE)), by=.(protID=omic, category)]

j <- merge(m3, m1a, by=c("protID","category"), all.x=TRUE)
fwrite(j[order(category,-absNIE)], file.path(OUT,"mediation_x_variance.tsv"), sep="\t")
cat(sprintf("wrote mediation_x_variance.tsv | %d (category,protein) rows | %d categories | M1 matched %d\n",
            nrow(j), uniqueN(j$category), j[is.finite(m1_r2),.N]))
