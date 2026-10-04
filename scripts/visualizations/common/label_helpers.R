#!/usr/bin/env Rscript

# ============================================================================
# label_helpers.R — shared exposure / protein / disease / component labels
# ----------------------------------------------------------------------------
# Pretty-printing and category lookups that were duplicated across the legacy
# Module1/Module2/Module3/ModulePred plotting scripts. Backed by the canonical
# config files so labels stay in sync with the analysis configuration.
# ============================================================================

local({
  if (exists("heap_config", mode = "function")) return(invisible())
  cand <- c(file.path(getwd(), "scripts", "visualizations", "common", "figure_paths.R"),
            "/n/groups/patel/shakson_ukb/HEAP/scripts/visualizations/common/figure_paths.R")
  hit <- cand[file.exists(cand)][1]
  if (!is.na(hit)) source(hit)
})

suppressPackageStartupMessages({ library(data.table) })

# ---------------------------------------------------------------------------
# Exposure variable -> fine category -> broad category
# ---------------------------------------------------------------------------

#' Exposure category-group map (fine -> broad), from the canonical config TSV.
#' @return data.table or NULL if the config is absent
heap_exposure_category_map <- function() {
  f <- heap_config("exposure_sets", "analysis_exposure_category_groups.tsv")
  if (!file.exists(f)) return(NULL)
  fread(f)
}

#' Map fine exposure categories (e.g. "Diet_Weekly", "Exercise_MET") to broad
#' group labels (e.g. "Lifestyle: Diet and Activity") via the canonical config.
#' Unmapped categories are returned unchanged. Vectorised.
heap_broad_category <- function(category, label = TRUE) {
  m <- heap_exposure_category_map()
  if (is.null(m)) return(as.character(category))
  col <- if (label && "broad_group_label" %in% names(m)) "broad_group_label" else "broad_group"
  lut <- setNames(m[[col]], m$category)
  out <- unname(lut[as.character(category)])
  out[is.na(out)] <- as.character(category)[is.na(out)]
  out
}

#' Active analysis exposures with include flags (canonical config TSV).
heap_exposure_table <- function() {
  f <- heap_config("exposure_sets", "analysis_exposures.tsv")
  if (!file.exists(f)) return(NULL)
  fread(f)
}

# ---------------------------------------------------------------------------
# British -> American spelling normaliser for display labels.
# ----------------------------------------------------------------------------
# Disease names (ICD-10/UKB), Reactome pathway names and a few exposure labels
# carry British spellings ("behavioural", "fibre", "ischaemic", "oesophageal",
# "haemorrhage", "anaemia", "oedema", "diarrhoea"). The manuscript is written in
# American English, so every label that reaches a figure routes through this
# normaliser. Case is preserved for the lowercase / Title / UPPER variants of
# each word (so "Ischaemic", "ISCHAEMIC" and "ischaemic" all normalise correctly).
# ---------------------------------------------------------------------------

# British substring -> American substring (lowercase canonical forms). Each is a
# fixed substring match, so morphological variants are covered (e.g. "behaviour"
# also fixes "behavioural"; "haem" fixes haemo-/haema-/-haemorrhage/-haematuria).
HEAP_BRITISH_AMERICAN <- c(
  behaviour = "behavior", fibre = "fiber", oedema = "edema",
  oesophag  = "esophag",  oestrog = "estrog", ischaem = "ischem",
  haem = "hem", anaem = "anem", rrhoea = "rrhea", pnoea = "pnea",
  tumour = "tumor", leukaem = "leukem", coeliac = "celiac",
  paediatr = "pediatr", gynaec = "gynec", anaesth = "anesth",
  foetal = "fetal", foetus = "fetus", caesar = "cesar",
  orthopaed = "orthoped", amoeb = "ameb", colour = "color",
  favour = "favor", licence = "license", defence = "defense",
  catalogue = "catalog",
  # NB: substitution is SUBSTRING-based (gsub fixed=TRUE), so every entry must be safe
  # as a substring of a correct American word. Deliberately EXCLUDED for that reason:
  #   analyse->analyze      corrupts "analyses"        -> "analyzes"
  #   characteris->characteriz corrupts "characteristics" -> "characteriztics"
  #   organis->organiz      corrupts "organism"        -> "organizm"
  signalling = "signaling", modelling = "modeling",
  labelled = "labeled", labelling = "labeling",
  normalis = "normaliz", standardis = "standardiz", summaris = "summariz",
  utilis = "utiliz", recognis = "recogniz", minimis = "minimiz", maximis = "maximiz",
  categoris = "categoriz", prioritis = "prioritiz",
  neighbour = "neighbor", centre = "center", grey = "gray")

#' Normalise British spellings to American in display labels. Vectorised; applies
#' each substitution in lowercase, Title-case and UPPER-case so the surrounding
#' casing of the matched word is preserved.
heap_americanize <- function(x) {
  x <- as.character(x)
  ttl <- function(s) { substr(s, 1, 1) <- toupper(substr(s, 1, 1)); s }
  for (i in seq_along(HEAP_BRITISH_AMERICAN)) {
    b <- names(HEAP_BRITISH_AMERICAN)[i]; a <- unname(HEAP_BRITISH_AMERICAN[i])
    x <- gsub(b, a, x, fixed = TRUE)                       # lower
    x <- gsub(ttl(b), ttl(a), x, fixed = TRUE)             # Title
    x <- gsub(toupper(b), toupper(a), x, fixed = TRUE)     # UPPER
  }
  x
}

# ---------------------------------------------------------------------------
# UKB-field name prettifier (shared helper, previously inlined many times)
# ---------------------------------------------------------------------------

#' Canonical SHORT exposure labels (config/exposure_sets/exposure_labels.tsv) —
#' the labeling analogue of the canonical category palette. Single source of
#' truth for short exposure names across every figure/table in the study.
#' Generated/edited via scripts/visualizations/generate_exposure_labels.R.
heap_exposure_label_map <- function() {
  f <- heap_config("exposure_sets", "exposure_labels.tsv")
  if (!file.exists(f)) return(NULL)
  fread(f)
}

#' Short display label for one or more exposure variable IDs. Falls back to the
#' generic prettifier for any variable not in the registry. Vectorised.
heap_exposure_label <- function(x) {
  m <- heap_exposure_label_map()
  x <- as.character(x)
  if (is.null(m)) return(heap_pretty_field(x))
  lut <- setNames(m$short_label, m$variable)
  out <- unname(lut[x])
  miss <- is.na(out)
  if (any(miss)) out[miss] <- heap_pretty_field(x[miss])
  out
}

#' Strip UK Biobank field suffixes and humanize a variable name.
#' e.g. "bread_intake_f1438_0_0" -> "Bread intake"
heap_pretty_field <- function(x) {
  x <- gsub("_f[0-9]+_[0-9]+_[0-9]+$", "", x)   # drop _fXXXX_i_a
  x <- gsub("_f[0-9]+$", "", x)
  x <- gsub("_+", " ", x)
  x <- trimws(x)
  # Title-case first letter only (keep acronyms intact-ish).
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  heap_americanize(x)
}

#' Humanize a FinnGen/UKB disease ID into a short readable label.
#' e.g. "age_e78_first_reported_disorders_of_lipoprotein_metabolism..." ->
#'      "Disorders of lipoprotein metabolism..."
#'      "finngen_R12_E4_OBESITYCAL" -> "Obesitycal"
#'      "finngen_R12_AB1_GASTROENTERITIS_NOS" -> "Gastroenteritis NOS"
#' Strips the UKB age_*_first_reported_ wrapper and the FinnGen
#' finngen_R<NN>_ + leading ICD chapter code, then humanizes. FinnGen endpoint
#' names keep their own (upper-case) casing so acronyms (T2D, AF, NOS) stay
#' intact rather than being mangled by sentence-casing. Vectorised.
heap_pretty_disease <- function(x) {
  x <- as.character(x)
  is_fg <- grepl("^finngen_R[0-9]+_", x)
  # --- UKB first-reported form ---
  x <- gsub("^age_[a-z][0-9]+_first_reported_", "", x)
  x <- gsub("_f[0-9]+_[0-9]+_[0-9]+$", "", x)
  # --- FinnGen form: drop release prefix + leading ICD chapter code token ---
  if (any(is_fg)) {
    x[is_fg] <- sub("^finngen_R[0-9]+_", "", x[is_fg])
    # leading chapter code: letter(s) + digit(s) (E4_, AB1_, I9_, M13_); leave
    # acronym-style codes without a trailing letter-block (T2D_, AUD) alone.
    x[is_fg] <- sub("^[A-Z]+[0-9]+_", "", x[is_fg])
  }
  x <- gsub("_+", " ", x)
  heap_americanize(trimws(x))
}

# ---------------------------------------------------------------------------
# GTEx tissue display names (canonical short labels for the enrichment figures).
# Raw GTEx v10 keys are underscored/lowercase (e.g. brain_putamen_basal_ganglia);
# this maps them to clean "Organ (region)" labels. Fallback: drop a leading
# "cells_", underscores -> spaces, sentence case.
# ---------------------------------------------------------------------------
HEAP_TISSUE_LABELS <- c(
  adipose_subcutaneous = "Adipose (subcut.)", adipose_visceral_omentum = "Adipose (visceral)",
  adrenal_gland = "Adrenal gland", artery_aorta = "Artery (aorta)",
  artery_coronary = "Artery (coronary)", artery_tibial = "Artery (tibial)",
  brain_amygdala = "Brain (amygdala)", brain_anterior_cingulate_cortex_ba24 = "Brain (cingulate ctx)",
  brain_caudate_basal_ganglia = "Brain (caudate)", brain_cerebellar_hemisphere = "Brain (cerebellar hem.)",
  brain_cerebellum = "Brain (cerebellum)", brain_cortex = "Brain (cortex)",
  brain_frontal_cortex_ba9 = "Brain (frontal ctx)", brain_hippocampus = "Brain (hippocampus)",
  brain_hypothalamus = "Brain (hypothalamus)", brain_nucleus_accumbens_basal_ganglia = "Brain (nuc. accumbens)",
  brain_putamen_basal_ganglia = "Brain (putamen)", brain_spinal_cord_cervical_c1 = "Spinal cord",
  brain_substantia_nigra = "Brain (substantia nigra)", breast_mammary_tissue = "Breast (mammary)",
  cells_cultured_fibroblasts = "Fibroblasts (cultured)", cells_ebv_transformed_lymphocytes = "Lymphocytes (EBV)",
  cervix_ectocervix = "Cervix (ecto)", cervix_endocervix = "Cervix (endo)",
  colon_sigmoid = "Colon (sigmoid)", colon_transverse = "Colon (transverse)",
  esophagus_gastroesophageal_junction = "Esophagus (GE junction)", esophagus_mucosa = "Esophagus (mucosa)",
  esophagus_muscularis = "Esophagus (muscularis)", fallopian_tube = "Fallopian tube",
  heart_atrial_appendage = "Heart (atrial)", heart_left_ventricle = "Heart (LV)",
  kidney_cortex = "Kidney (cortex)", kidney_medulla = "Kidney (medulla)",
  liver = "Liver", lung = "Lung", minor_salivary_gland = "Salivary gland",
  muscle_skeletal = "Muscle (skeletal)", nerve_tibial = "Nerve (tibial)",
  ovary = "Ovary", pancreas = "Pancreas", pituitary = "Pituitary", prostate = "Prostate",
  skin_not_sun_exposed_suprapubic = "Skin (not sun-exp.)", skin_sun_exposed_lower_leg = "Skin (sun-exp.)",
  small_intestine_terminal_ileum = "Small intestine", spleen = "Spleen", stomach = "Stomach",
  testis = "Testis", thyroid = "Thyroid", uterus = "Uterus", vagina = "Vagina",
  whole_blood = "Whole blood")

#' Clean GTEx tissue display label(s). Vectorised; unknown keys fall back to a
#' general prettifier (drop leading "cells_", underscores -> spaces, sentence case).
heap_pretty_tissue <- function(x) {
  x <- as.character(x)
  out <- unname(HEAP_TISSUE_LABELS[x])
  miss <- is.na(out)
  if (any(miss)) {
    g <- gsub("^cells_", "", x[miss]); g <- trimws(gsub("_+", " ", g))
    substr(g, 1, 1) <- toupper(substr(g, 1, 1)); out[miss] <- g
  }
  heap_americanize(out)
}

#' Wrap long labels to multiple lines at a word boundary (for axis tick labels
#' such as Reactome pathway names). Vectorised.
heap_wrap_label <- function(x, width = 40L) {
  x <- heap_americanize(x)
  vapply(as.character(x), function(s) paste(strwrap(s, width = width), collapse = "\n"),
         character(1), USE.NAMES = FALSE)
}

#' ICD-10 chapter (organ system) for a UKB first-reported disease ID.
#'
#' The DZ_ID embeds the ICD-10 code as "age_<letter><nn>_first_reported_...".
#' Maps the leading letter (with the C/D neoplasm-vs-blood and H eye-vs-ear
#' splits handled by the two-digit code) to a short organ-system label used to
#' group/colour diseases in the mediation figures. Vectorised.
#' @param x disease IDs
#' @return character vector of system labels (NA-safe; "Other" for unmatched)
heap_icd_chapter <- function(x) {
  x <- tolower(as.character(x))
  m  <- regmatches(x, regexpr("^age_([a-z])([0-9]{2})", x))
  let <- toupper(substr(sub("^age_", "", m), 1, 1))
  num <- suppressWarnings(as.integer(substr(sub("^age_.", "", m), 1, 2)))
  let[!nzchar(m)] <- NA; num[!nzchar(m)] <- NA
  base <- c(A="Infectious", B="Infectious",
            C="Neoplasms",
            E="Endocrine/metabolic", F="Mental & behavioral", G="Nervous",
            I="Circulatory", J="Respiratory", K="Digestive", L="Skin",
            M="Musculoskeletal", N="Genitourinary",
            O="Pregnancy", P="Perinatal", Q="Congenital",
            R="Symptoms/signs", S="Injury", T="Injury")
  out <- unname(base[let])
  out[let == "D" & !is.na(num) & num <= 48] <- "Neoplasms"
  out[let == "D" & !is.na(num) & num >= 50] <- "Blood & immune"
  out[let == "H" & !is.na(num) & num <= 59] <- "Eye"
  out[let == "H" & !is.na(num) & num >= 60] <- "Ear"
  out[is.na(out)] <- "Other"
  out
}

#' Collapse Module 1 R2groups fine group labels into coarse variance components.
#'
#' R2groups `group` values look like: "Covars", "Gcis", "Gtrans", "E_<category>",
#' "GxEcis_<category>", "GxEtrans_<category>". This maps each to one of:
#' Covariates / Genetic / Exposome / GxE (factor, ordered for stacking).
#' @param group character vector of R2groups group labels
heap_r2_coarse_component <- function(group) {
  g <- as.character(group)
  out <- ifelse(g == "Covars", "Covariates",
         ifelse(g %in% c("Gcis", "Gtrans"), "Genetic",
         ifelse(grepl("^GxE", g), "GxE",
         ifelse(grepl("^E_", g), "Exposome", "Other"))))
  factor(out, levels = c("Covariates", "Genetic", "Exposome", "GxE", "Other"))
}

#' Extract the exposure category from an E_/GxEcis_/GxEtrans_ group label.
#' e.g. "GxEtrans_Diet_Weekly" -> "Diet_Weekly"; "Covars"/"Gcis" -> NA.
heap_r2_exposure_category <- function(group) {
  g <- as.character(group)
  cat <- sub("^(E|GxEcis|GxEtrans)_", "", g)
  ifelse(grepl("^(E_|GxEcis_|GxEtrans_)", g), cat, NA_character_)
}

#' Canonical pretty names for variance components / model blocks.
heap_pretty_component <- function(x) {
  map <- c(C = "Covariates", G = "Genetic", E = "Exposome",
           PGS = "Genetic (PGS)", PXS = "Exposome (PXS)",
           GxE = "Gene x Exposure", Covars = "Covariates",
           GStrans = "Trans genetic", GScis = "Cis genetic")
  out <- map[as.character(x)]
  out[is.na(out)] <- as.character(x)[is.na(out)]
  unname(out)
}

# ---------------------------------------------------------------------------
# Module 3 (GEM mediation) driver labels — shared across the mediation figures.
# The Module 3 `predictor` column is one of: G_raw / PXS_total (primary_total),
# Gcis_raw / Gtrans_raw / PXS_<Category> (partitioned_categories), or
# Gcis_raw / Gtrans_raw / PXSgrp_<Group> (partitioned_grouped_categories). These
# helpers translate that vocabulary into the canonical exposure-category codes
# (so scale_*_exposure() colours apply) and readable driver labels.
# ---------------------------------------------------------------------------

#' Coarse driver component (Genetic vs Exposome) from a Module 3 predictor_class.
#' Handles primary (genetic_total/exposure_total) and partitioned
#' (genetic_cis/genetic_trans/exposure_category/exposure_group) vocabularies.
heap_md_driver_component <- function(predictor_class) {
  g <- as.character(predictor_class)
  out <- ifelse(grepl("genetic", g), "Genetic",
         ifelse(grepl("exposure", g), "Exposome", "Other"))
  factor(out, levels = c("Genetic", "Exposome", "Other"))
}

#' Driver GROUP that keeps the cis/trans genetic split, for figures that compare
#' genetic-cis vs genetic-trans vs exposomic proportion-mediated / effect. Maps a
#' Module 3 predictor_class to a display group. The partitioned model supplies
#' genetic_cis/genetic_trans; the primary model's pooled genetic_total and the
#' exposomic classes collapse to their natural group.
heap_md_driver_group <- function(predictor_class) {
  g <- as.character(predictor_class)
  out <- c(genetic_cis       = "Genetic (cis)",
           genetic_trans     = "Genetic (trans)",
           genetic_total     = "Genetic (total)",
           exposure_total    = "Exposomic (PXS)",
           exposure_category = "Exposomic (PXS)",
           exposure_group    = "Exposomic (PXS)")[g]
  out[is.na(out)] <- g[is.na(out)]
  factor(unname(out), levels = c("Genetic (cis)", "Genetic (trans)",
                                 "Genetic (total)", "Exposomic (PXS)"))
}

#' Canonical colours for the cis/trans/exposomic driver groups (matches the
#' cis=blue, trans=light-blue, exposome=orange convention used elsewhere).
heap_md_driver_group_colors <- function() {
  c("Genetic (cis)"   = unname(HEAP_PAL_COMPONENT[["Genetic"]]),
    "Genetic (trans)" = "#7FB3D5",
    "Genetic (total)" = "#3A6EA5",
    "Exposomic (PXS)" = unname(HEAP_PAL_COMPONENT[["Exposome"]]))
}

#' Fine exposure-category CODE from a Module 3 predictor (canonical key for the
#' HEAP exposure palette). Strips the PXS_/PXSgrp_ prefix; genetic predictors and
#' PXS_total map to NA. e.g. "PXS_Diet_Weekly" -> "Diet_Weekly".
heap_md_category <- function(predictor) {
  p <- as.character(predictor)
  cat <- sub("^PXS(grp)?_", "", p)
  ifelse(grepl("^PXS(grp)?_", p) & p != "PXS_total", cat, NA_character_)
}

#' Readable driver label for legends/axes. Genetic components get explicit
#' cis/trans/total labels; exposure categories get the prettified category name.
#' e.g. "Gcis_raw"->"Genetic (cis)", "PXS_Smoking"->"Smoking",
#'      "PXS_total"->"Exposome (total)".
heap_md_predictor_label <- function(predictor) {
  p <- as.character(predictor)
  gmap <- c(G_raw = "Genetic (total)", Gcis_raw = "Genetic (cis)",
            Gtrans_raw = "Genetic (trans)", PXS_total = "Exposome (total)")
  out <- gmap[p]
  cat <- heap_md_category(p)
  out[is.na(out) & !is.na(cat)] <- gsub("_", " ", cat[is.na(out) & !is.na(cat)])
  out[is.na(out)] <- gsub("_", " ", p[is.na(out)])
  unname(out)
}
