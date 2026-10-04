# ============================================================================
# program_clusters.R  -- canonical biological-program cluster assignment for the
# Module-2 exposure -> program -> tissue analysis (curated 2026-06-23).
#
# heap_program_cluster(desc): map a Reactome pathway Description to a coarse
#   program cluster. Order = most specific first (immune > ECM > neuronal >
#   vascular > growth-factor/RTK > glyco/lipid metabolism > muscle > other).
# HEAP_PROGRAM_TISSUES: GTEx tissue -> organ-system grouping (Brain split out).
# HEAP_PROGRAM_LEVELS / HEAP_TISSUE_LEVELS: display order for the tripartite.
# HEAP_TRIPARTITE_EXEMPLARS: curated exemplar exposures (one true direction each).
#
# Provenance: cluster membership audited in
#   figures/exploratory/module2/curation/cluster_membership.tsv
# ============================================================================

heap_program_cluster <- function(desc) {
  x <- tolower(desc)
  data.table::fcase(
    grepl("adaptive immune", x), "Adaptive immune",
    grepl("neutrophil|innate immune|interleukin|immunoregulatory|tnf|nf-kb|cytokine|chemokine|complement|interferon|toll|antigen|nod1|nod2|clec7a|dectin|dap12|inflammasome|nlrp|\\(nlr\\)|leucine rich repeat", x), "Innate immune",
    grepl("proteoglycan|glycosaminoglycan|hs-gag|chondroitin|dermatan|heparan|syndecan|b4galt7|b3gat3|b3galt6|tetrasaccharide|cs/ds|collagen|matrix|integrin|elastic|metalloprotein|keratin|cornified|laminin|adherens|cell-cell junction|cell junction organization", x), "ECM / proteoglycan",
    grepl("neuronal|synapse|synaptic|axon|l1cam|nervous system|chemical synapses|neuro", x), "Neuronal / synaptic",
    grepl("vascular wall|platelet|hemostas|coagul|fibrin|clot|scavenger|angiotensin|kallikrein|kinin|contact activation", x), "Vascular / hemostasis / RAAS",
    grepl("insulin-like growth|igfbp|tgfb|growth factor|peptide hormone|peptide ligand|eph-ephrin|ephrin|erbb|fgfr|ntrk|trka|notch", x), "Growth-factor / RTK",
    grepl("glycosphingolipid|sphingolipid|glycosylation|lipid|adipo|cholesterol|fatty|sterol|triglycer|biological oxidation|carbohydrate metabolism|amino acid|nucleotide|^metabolism|metabolism of|diseases of metabolism", x), "Glycan / lipid metabolism",
    grepl("muscle|contraction", x), "Muscle",
    rep(TRUE, length(x)), "Other")
}

HEAP_PROGRAM_LEVELS <- c("Innate immune","Adaptive immune","ECM / proteoglycan","Neuronal / synaptic",
                         "Growth-factor / RTK","Glycan / lipid metabolism","Vascular / hemostasis / RAAS","Muscle")
HEAP_TISSUE_LEVELS  <- c("Blood","Lung","Liver","Adipose","Skin","Vascular","Brain","Muscle")

HEAP_PROGRAM_TISSUES <- list(
  Lung="lung",
  Liver=c("liver","stomach","small_intestine_terminal_ileum","colon_transverse","minor_salivary_gland","esophagus_mucosa"),
  Blood=c("spleen","whole_blood","cells_ebv-transformed_lymphocytes"),
  Adipose=c("adipose_subcutaneous","adipose_visceral_omentum"),
  Vascular=c("artery_coronary","artery_aorta","artery_tibial"),
  Skin=c("skin_not_sun_exposed_suprapubic","skin_sun_exposed_lower_leg","cells_cultured_fibroblasts"),
  Brain=c("brain_amygdala","brain_hippocampus","brain_anterior_cingulate_cortex_ba24","brain_caudate_basal_ganglia",
          "brain_cortex","brain_frontal_cortex_ba9","brain_nucleus_accumbens_basal_ganglia","brain_hypothalamus",
          "nerve_tibial","brain_putamen_basal_ganglia","pituitary","brain_substantia_nigra",
          "brain_cerebellar_hemisphere","brain_cerebellum"),
  Muscle=c("muscle_skeletal","heart_left_ventricle","heart_atrial_appendage"))

# curated exemplar exposures: harmful (red) then protective (blue); each has a
# single coherent direction, so exposure->program edges are unambiguous.
HEAP_TRIPARTITE_EXEMPLARS <- data.frame(
  e   = c("pack_years_of_smoking_f20161_0_0","no2_mean","index_of_multiple_deprivation_england_f26410_0_0",
          "processed_meat_intake_f1349_0_0","beef_intake_f1369_0_0",
          "oily_fish_intake_f1329_0_0","usual_walking_pace_f924_0_0","summed_days_activity_f22033_0_0",
          "cereal_intake_f1458_0_0","dried_fruit_intake_f1319_0_0"),
  lab = c("Smoking","Air pollution (NO₂)","Deprivation","Processed meat","Red meat",
          "Oily fish","Walking pace","Active days/wk","Cereal / fiber","Dried fruit"),
  grp = c(rep("Harmful",5), rep("Protective",5)),
  stringsAsFactors = FALSE)
