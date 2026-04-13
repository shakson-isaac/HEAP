suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(ggplot2)
})

# LOAD DATA:
MRres <- fread(file = "/n/groups/patel/shakson_ukb/UK_Biobank/Output/MR_edges/summary/MRmotifs.csv")
HEAPassoc <- fread("/n/groups/patel/shakson_ukb/UK_Biobank/Output/App/Tables/ReplicatedEassoc.csv")
HEAPdz <- fread("/n/groups/patel/shakson_ukb/UK_Biobank/Output/Mediation/Results/Type5.csv")

MDsig <- function(DF){
  alpha = 0.05
  bonf_corr <- nrow(DF)
  
  res <- DF %>%
    filter(
      pxs_NIE_logHR_delta_p      < alpha/bonf_corr |
        gcis_NIE_logHR_delta_p   < alpha/bonf_corr |
        gtrans_NIE_logHR_delta_p < alpha/bonf_corr
    )
  
  return(res)
}
HEAPdz <- MDsig(HEAPdz)

#Compare MRres and HEAPassoc:
MR_HEAP_EP <- merge(MRres, HEAPassoc, by.x = c("Exposure","Protein"),
      by.y = c("Eid_train","omicID"), allow.cartesian = TRUE)





#Get MR hit values for EP:
EP_viz <- MR_HEAP_EP %>%
                select(c("Exposure","Protein",
                         "beta_EP","padj_EP",
                         "Estimate_train",
                         "Pr(>|t|)_train")) %>%
                unique()

colnames(EP_viz) <- c("Exposure","Protein",
                      "MR_Beta","MR_P",
                      "HEAP_Beta","HEAP_P")

plot(EP_viz$MR_Beta, EP_viz$HEAP_Beta)
cor(EP_viz$MR_Beta, EP_viz$HEAP_Beta, use = "complete.obs")
  
EP_viz_MRsig <- EP_viz %>% filter(MR_P < 0.05)
plot(EP_viz_MRsig$MR_Beta, EP_viz_MRsig$HEAP_Beta)
cor(EP_viz_MRsig$MR_Beta, EP_viz_MRsig$HEAP_Beta, use = "complete.obs")

# Get MRhit and compare to Indirect Effects from Mediation Analysis
#Compare MRres and HEAPassoc:

#Get FinnGen code for HEAPdz
# Load the UKB/FinnGen Disease Mapping:
UKBfinngen <- fread("/n/groups/patel/IGLOO/FinnGen/UKBFinnGenDisease.csv")

HEAPdz <- merge(HEAPdz, UKBfinngen, 
                by.x = c("DZ_ID"),
                by.y = c("Disease"))

HEAPdz <- HEAPdz %>%
  rename(Disease = "FinnGen") %>%
  separate_rows(Disease, sep = ";\\s*") %>%
  mutate(Disease = paste0("finngen_R12_", Disease))



MR_HEAP_PD <- merge(MRres, HEAPdz, 
                    by.x = c("Protein","Disease"),
                    by.y = c("prot_term","Disease"), 
                    allow.cartesian = TRUE)



plot(MR_HEAP_PD$beta_PDcis, MR_HEAP_PD$gcis_NIE_HR)
cor(MR_HEAP_PD$beta_PDcis, MR_HEAP_PD$gcis_NIE_HR, use = "complete.obs")

# Mediation Plot Stat:
plot(MR_HEAP_PD$beta_EP * MR_HEAP_PD$beta_PDcis, MR_HEAP_PD$pxs_NIE_logHR)
cor(MR_HEAP_PD$beta_EP * MR_HEAP_PD$beta_PDcis, MR_HEAP_PD$pxs_NIE_logHR, use = "complete.obs")


cor(MR_HEAP_PD$beta_PDcis, MR_HEAP_PD$pxs_NIE_logHR, use = "complete.obs")
cor(MR_HEAP_PD$beta_PDtrans, MR_HEAP_PD$pxs_NIE_logHR, use = "complete.obs")
cor(MR_HEAP_PD$beta_DP, MR_HEAP_PD$pxs_NIE_logHR, use = "complete.obs")


#PD association consistency:
PD_check <- MR_HEAP_PD %>%
                filter(padj_PDcis < 0.05 &
                        gcis_NIE_logHR_delta_p < (0.05/713076))

cor.test(PD_check$beta_PDcis, PD_check$gcis_NIE_logHR, use = "complete.obs")
cor.test(PD_check$beta_PDcis, PD_check$gcis_NDE_logHR, use = "complete.obs")
#cor(PD_check$beta_PDcis, PD_check$gcis_NDE_HR, use = "complete.obs")

#cor(PD_check$beta_DP, PD_check$gcis_NIE_logHR, use = "complete.obs")


# Indirect Effect Check:
PXS_IEcis <- MR_HEAP_PD %>%
              filter(padj_PDcis < 0.05 &
                       padj_EP < 0.05 &
                       pxs_NIE_logHR_delta_p < (0.05/713076))

cor(PXS_IEcis$beta_EP, PXS_IEcis$pxs_NIE_logHR, use = "complete.obs")
cor(PXS_IEcis$beta_PDcis, PXS_IEcis$pxs_NIE_logHR, use = "complete.obs")
cor.test(PXS_IEcis$beta_PDcis, PXS_IEcis$pxs_NIE_logHR, use = "complete.obs")


PXS_IEcisv2 <- MR_HEAP_PD %>%
  filter(padj_PDcis < 0.05 &
           pxs_NIE_logHR_delta_p < (0.05/713076))
cor.test(PXS_IEcisv2$beta_PDcis, PXS_IEcisv2$pxs_NIE_logHR, use = "complete.obs")



PXS_IEtrans <- MR_HEAP_PD %>%
  filter(padj_PDtrans < 0.05 &
           padj_EP < 0.05 &
           pxs_NIE_logHR_delta_p < (0.05/713076))
cor.test(PXS_IEtrans$beta_PDtrans, PXS_IEtrans$pxs_NIE_logHR, use = "complete.obs")



### Maybe Module version?
MotifA <- MR_HEAP_PD %>% filter(motif_label == "A")

