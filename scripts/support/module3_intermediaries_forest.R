#!/usr/bin/env Rscript
# module3_intermediaries_forest.R  -> intermediaries_forest.tsv
# Dominant-exposure NIE + 95% CI per (protein, disease) for the forest-plot panel d.
# Uses RAW partitioned NIE rows (which carry delta_l95_HR/u95_HR), n>=100.
local({ cm <- c(file.path(getwd(),"scripts","visualizations","common"),"/n/groups/patel/shakson_ukb/HEAP/scripts/visualizations/common"); cm <- cm[dir.exists(cm)][1]
  for (f in c("figure_paths","load_heap_results","plot_theme","label_helpers")) source(file.path(cm,paste0(f,".R"))) })
suppressPackageStartupMessages(library(data.table))
OUT <- file.path(Sys.getenv("HEAP_PROJECT_ROOT", "/n/groups/patel/IGLOO/UKB/HEAP"), "figures", "exploratory/module3"); NMIN <- 100
sel <- c("protID","DZ_ID","predictor","predictor_class","effect_type","effect_logHR","effect_HR",
         "delta_se","delta_l95_HR","delta_u95_HR","delta_p","n_cases","protein_HR")
md <- load_module3_results(covarType="base", family="lasso", mode="partitioned_categories", select=sel)
nie <- md[effect_type=="NIE" & predictor_class=="exposure_category" & is.finite(effect_logHR) & n_cases>=NMIN]
heap_md_fdr(nie, "delta_p", "delta_q")
nie[, category := heap_md_category(predictor)][, a := abs(effect_logHR)]
nie <- nie[!is.na(category)]
dom <- nie[order(protID,DZ_ID,-a)][, .(dom_cat=category[1], NIE_HR=effect_HR[1], l95=delta_l95_HR[1],
  u95=delta_u95_HR[1], dom_q=delta_q[1], n_cases=n_cases[1], protein_HR=protein_HR[1]), by=.(protID,DZ_ID)]
dom <- dom[dom_q<0.05 & is.finite(l95) & is.finite(u95)]
sysmap <- function(dz){ L<-toupper(sub("^age_([a-zA-Z]).*","\\1",dz))
  m<-c(E="Endocrine/metabolic",I="Circulatory",N="Renal/GU",J="Respiratory",K="Digestive",
       D="Blood/immune",M="Musculoskeletal",F="Mental/neuro"); o<-m[L]; o[is.na(o)]<-"Other"; unname(o) }
dom[, disease := heap_pretty_disease(DZ_ID)][, system := sysmap(DZ_ID)][, aa := abs(log(NIE_HR))]
# join PM from disease_mediators
pm <- fread(file.path(OUT,"disease_mediators.tsv"))[, .(protID, DZ_ID, PM)]
dom <- merge(dom, pm, by=c("protID","DZ_ID"), all.x=TRUE)
fwrite(dom[order(-aa)], file.path(OUT,"intermediaries_forest.tsv"), sep="\t")
cat(sprintf("wrote intermediaries_forest.tsv | %d (protein,disease) dominant-exposure NIE+CI rows\n", nrow(dom)))
