#!/usr/bin/env Rscript
# ============================================================================
# module3_category_disease_full.R  (support: full exposure x disease axis structure)
# ----------------------------------------------------------------------------
# Per (exposure category, disease): how many proteins carry an FDR-significant,
# sign-consistent exposure-mediated effect (NIE) -- using the PARTITIONED NIE for
# THAT category (not the dominant-exposure collapse), so every axis is visible:
# smoking->respiratory, alcohol->liver/mental, diet->digestive, etc., not just the
# cardiometabolic-renal core. Output: cat_disease_full.tsv.
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
MEFF <- log(1.10); NMIN <- 100
sel <- c("protID","DZ_ID","predictor","predictor_class","effect_type","effect_logHR","effect_HR",
         "delta_p","instrument_present","n","n_cases","protein_HR","protein_p")
md <- load_module3_results(covarType=covarType, family=family, mode="partitioned_categories", select=sel)
md <- md[!(predictor_class %in% c("genetic_cis","genetic_trans") & instrument_present==FALSE)]
pm <- heap_proportion_mediated(md); pm[, category := heap_md_category(predictor)]
e <- pm[!is.na(category) & is.finite(NIE_q) & NIE_q<0.05 & pm_consistent==TRUE & n_cases>=NMIN]
e[, disease := heap_pretty_disease(DZ_ID)]
sysmap <- function(dz){ L<-toupper(sub("^age_([a-zA-Z]).*","\\1",dz))
  m<-c(A="Infection",B="Infection",C="Neoplasm",D="Blood/immune",E="Endocrine/metabolic",F="Mental/neuro",
       G="Mental/neuro",H="Eye/ear",I="Circulatory",J="Respiratory",K="Digestive",L="Skin",
       M="Musculoskeletal",N="Renal/GU",R="Symptoms"); o<-m[L]; o[is.na(o)]<-"Other"; unname(o) }
e[, system := sysmap(DZ_ID)]
g <- e[, .(n_prot=uniqueN(protID), n_strong=uniqueN(protID[abs(NIE_logHR)>=MEFF]),
           dir_pos=mean(NIE_logHR>0)), by=.(category, DZ_ID, disease, system)]
fwrite(g[order(category,-n_prot)], file.path(OUT,"cat_disease_full.tsv"), sep="\t")

cat(sprintf("\n(category, disease) FDR-sig cells: %d | categories %d | diseases %d\n",
            nrow(g), uniqueN(g$category), uniqueN(g$DZ_ID)))
cat("\n=== axes by organ system (total mediator-protein count across categories) ===\n")
print(e[, .(n_links=.N, n_prot=uniqueN(protID)), by=system][order(-n_links)])
for (cc in c("Smoking","Alcohol","Diet_Weekly","Exercise_Freq","Sleep","Sun_Exposure")) {
  cat(sprintf("\n=== %s : top mediated diseases (FDR-sig) ===\n", cc))
  print(head(g[category==cc][order(-n_prot), .(disease, system, n_prot, n_strong, dir=round(dir_pos,2))], 10))
}
cat("\nwrote cat_disease_full.tsv\n")
