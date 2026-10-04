#!/usr/bin/env Rscript

# ============================================================================
# run_intervention_compare.R  —  HEAP support analysis
# ----------------------------------------------------------------------------
# Compares HEAP exposure->protein effects against external intervention
# proteomics (GLP1 STEP1/STEP2 trials; HERITAGE exercise trial), weighting each
# protein by Olink<->SomaScan cross-platform reliability.
#
# Ported from the analysis half of the legacy
#   scripts/visualizations/Visualizations/analyHEAP_interventionv2.R
# but: (1) reads canonical IGLOO inputs, (2) reads Module 2 univar_assoc outputs
# instead of the 1.1 GB HEAPassoc.qs, (3) emits flat figure-ready TSVs instead of
# the S4 INTconstruct/.qs object. Plotting stays in Visualizations/Module4/*.
#
# Inputs (all IGLOO-canonical via workflow/00_paths.R):
#   Module 2 associations  heap_project_output("module2", <covarType>)/univar_assoc_*.rds
#   GLP1 trial             heap_interventions_or_legacy("GLP1_proteomics.xlsx")
#   HERITAGE trial         heap_interventions_or_legacy("jciinsight_prot.xlsx")
#   Olink<->SomaScan rel.  heap_olinksoma()
#
# Outputs -> heap_project_output("support", "intervention_compare"):
#   intervention_correlations.tsv  covarType,exposure_id,intervention,r,pval,pval_BH,n_eff,sig_any
#   intervention_scatter.tsv       covarType,exposure_id,protein,beta_HEAP,se_HEAP,
#                                  HERITAGE_effect,HERITAGE_se,GLP1_effect1,GLP1_se1,
#                                  GLP1_effect2,GLP1_se2,olink_soma_r
#   (both also mirrored to figures/data/)
#
# Usage:
#   module load gcc/14.2.0 R/4.4.2
#   HEAP_PATHS_FILE=/n/groups/patel/shakson_ukb/HEAP/workflow/00_paths.R \
#     Rscript scripts/support/intervention_compare/run_intervention_compare.R [covarType ...]
#   (no args => every covarType found under module2/; e.g. `... Type5`)
#
# Design notes (intentional differences from legacy):
#   - Outputs ALL exposure x intervention correlations plus a `sig_any` flag
#     (BH<0.05 in >=1 intervention) instead of pre-filtering, so the plotting
#     layer decides what to show.
#   - HEAP side: test split, Bonferroni 0.05/n over all tested exposure-protein
#     pairs in that covarType (matches legacy `Pr(>|t|) < 0.05/n()`).
#   - Missing reliability r_cross is mean-imputed; weight w_rel = pmax(r_cross,0)
#     (matches legacy v2 behavior).
# ============================================================================

local({
  cand <- c(Sys.getenv("HEAP_PATHS_FILE", unset = ""),
            file.path(getwd(), "workflow", "00_paths.R"),
            "/n/groups/patel/shakson_ukb/HEAP/workflow/00_paths.R")
  hit <- cand[nzchar(cand) & file.exists(cand)][1]
  if (is.na(hit)) stop("Could not locate workflow/00_paths.R (set HEAP_PATHS_FILE).")
  source(hit)
})

suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(tidyr); library(purrr)
  library(tibble); library(readxl); library(weights)
})

SPLIT <- "test"   # match legacy: correlate on the held-out split
MIN_NEFF <- 8     # effective-N floor: corr below this is excluded from the FDR family

# --- resolve which covarTypes to run ----------------------------------------
# Canonical Module 2 output is experiment-nested: module2/<experiment>/<covarType>/
# (the primary full-proteome base run is M2_base_main/base). Older Type* pilots
# used the flat layout module2/<covarType>/. We resolve both: prefer the nested
# experiment subdir, fall back to flat. The experiment defaults to M2_base_main
# and can be overridden with HEAP_M2_EXPERIMENT (empty -> flat legacy layout).
m2_root_igloo <- heap_project_output("module2")
m2_root_local <- heap_output("module2")
m2_root <- if (dir.exists(m2_root_igloo)) m2_root_igloo else m2_root_local
if (!dir.exists(m2_root))
  stop("No Module 2 output dir found (looked at IGLOO + local). Run Module 2 first.")

EXPERIMENT <- Sys.getenv("HEAP_M2_EXPERIMENT", unset = "M2_base_main")
exp_root   <- if (nzchar(EXPERIMENT)) file.path(m2_root, EXPERIMENT) else m2_root

# the directory holding univar_assoc_*.rds for one covarType (nested or flat)
.statE_dir <- function(covarType) {
  cand <- unique(c(if (nzchar(EXPERIMENT)) file.path(m2_root, EXPERIMENT, covarType),
                   file.path(m2_root, covarType)))
  has_rds <- function(d) dir.exists(d) &&
    length(list.files(d, pattern = "^univar_assoc_.*\\.rds$")) > 0
  hit <- cand[vapply(cand, has_rds, logical(1))]
  if (!length(hit))
    stop("No univar_assoc_*.rds for covarType '", covarType, "' (looked at: ",
         paste(cand, collapse = ", "), ") — run Module 2 first.")
  hit[1]
}

args <- commandArgs(trailingOnly = TRUE)
covarTypes <- if (length(args)) args else {
  src <- if (nzchar(EXPERIMENT) && dir.exists(exp_root)) exp_root else m2_root
  list.dirs(src, recursive = FALSE, full.names = FALSE)
}
covarTypes <- covarTypes[nzchar(covarTypes)]
if (!length(covarTypes))
  stop("No covarType subdirectories under ", exp_root)
message("Intervention compare | source: ", exp_root,
        " | experiment: ", if (nzchar(EXPERIMENT)) EXPERIMENT else "(flat)",
        " | covarTypes: ", paste(covarTypes, collapse = ", "))

# --- aggregate Module 2 statE (element 1) across proteins for one covarType --
load_statE <- function(covarType, split = SPLIT) {
  d <- .statE_dir(covarType)
  files <- list.files(d, pattern = "^univar_assoc_.*\\.rds$", full.names = TRUE)
  # unique(): statE rows (term x protein) are unique by construction; guards
  # against exact-duplicate rows from the fixed Module2.R batch-accumulator bug.
  unique(rbindlist(lapply(files, function(f) {
    obj <- readRDS(f); as.data.table(obj[[split]][[1]])  # [[1]] = statE
  }), fill = TRUE))
}

# --- external intervention inputs (read once) -------------------------------
heritage_path <- heap_interventions_or_legacy("jciinsight_prot.xlsx")
glp1_path     <- heap_interventions_or_legacy("GLP1_proteomics.xlsx")
olinksoma     <- heap_olinksoma()
for (p in c(heritage_path, glp1_path, olinksoma))
  if (!file.exists(p)) stop("Missing intervention input: ", p)

heritage_raw <- read_excel(heritage_path, sheet = excel_sheets(heritage_path)[1], skip = 2)
glp1_sheets  <- excel_sheets(glp1_path)
STEP1 <- read_excel(glp1_path, sheet = glp1_sheets[2])
STEP2 <- read_excel(glp1_path, sheet = glp1_sheets[3])

# Olink<->SomaScan reliability (skip=3 header offset, matches legacy)
prot_rel <- fread(olinksoma, skip = 3)
prot_rel <- prot_rel[, c("gene_name", "olink_nonnorm_corr", "olink_smpnorm_corr")]
setnames(prot_rel, c("EntrezGeneSymbol", "r_cross", "r_crossv2"))
prot_rel <- na.omit(prot_rel, cols = c("EntrezGeneSymbol", "r_cross"))
# OlinkSoma has duplicate gene symbols (multiple assays/panels); collapse to one
# row per gene so reliability joins stay 1:1 (the legacy left these duplicated,
# which silently double-counted those proteins in the weighted correlation).
prot_rel <- prot_rel[, .(r_cross = mean(r_cross), r_crossv2 = mean(r_crossv2)),
                     by = EntrezGeneSymbol]

# intervention effect tables (q<0.05), aggregated per gene symbol.
# HERITAGE, like GLP1, has duplicate gene symbols (multiple aptamers/probes);
# collapse to one row per gene (mean effect, max se) so the protein joins stay
# 1:1 — otherwise the scatter left_join fans out into a many-to-many that
# double-counts those proteins (and trips a dplyr warning).
heritage_eff <- heritage_raw %>%
  filter(`False Discovery Rate (q-value)` < 0.05) %>%
  mutate(HERITAGE_se = `log(10) Fold Change` / `t-statistic`) %>%
  rename(HERITAGE_effect = `log(10) Fold Change`) %>%
  group_by(EntrezGeneSymbol) %>%
  summarise(HERITAGE_effect = mean(HERITAGE_effect),
            HERITAGE_se     = max(HERITAGE_se), .groups = "drop")

glp1_eff <- function(STEP, eff_nm, se_nm) {
  STEP %>% filter(qvalue < 0.05) %>%
    group_by(EntrezGeneSymbol) %>%
    summarise(!!eff_nm := mean(effect_size), !!se_nm := max(std_error), .groups = "drop")
}
glp1_step1_eff <- glp1_eff(STEP1, "GLP1_effect1", "GLP1_se1")
glp1_step2_eff <- glp1_eff(STEP2, "GLP1_effect2", "GLP1_se2")

INTERVENTIONS <- c("HERITAGE_effect", "GLP1_effect1", "GLP1_effect2")

# --- weighted Pearson r + p via effective N (from legacy) --------------------
wtd_cor_and_p <- function(x, y, w) {
  ok <- is.finite(x) & is.finite(y) & is.finite(w) & (w > 0)
  x <- x[ok]; y <- y[ok]; w <- w[ok]
  if (length(x) < 3) return(c(r = NA_real_, p = NA_real_, neff = NA_real_))
  r <- suppressWarnings(wtd.cor(x, y, weight = w)[1, 1])
  neff <- (sum(w)^2) / sum(w^2)
  if (!is.finite(r) || !is.finite(neff) || neff <= 2)
    return(c(r = NA_real_, p = NA_real_, neff = neff))
  tval <- r * sqrt((neff - 2) / pmax(1e-12, 1 - r^2))
  c(r = r, p = 2 * pt(-abs(tval), df = neff - 2), neff = neff)
}

# --- per-covarType pipeline --------------------------------------------------
run_one <- function(covarType) {
  # REPLICATED associations: a (term ID x protein) pair must pass Bonferroni in
  # BOTH the train and test splits — the manuscript's significance criterion.
  # This also drops the ordered-factor high-degree polynomial-contrast terms
  # whose raw coefficients explode in one split but do not reproduce (e.g.
  # |beta|>>1 diet contrasts) and would otherwise corrupt the weighted
  # correlation. The displayed effect/SE is from the held-out TEST split.
  te <- load_statE(covarType, "test")
  tr <- load_statE(covarType, "train")
  te2 <- te[, .(ID, omicID, Eid, Category,
                Estimate, `Std. Error`, p_test = get("Pr(>|t|)"))]
  tr2 <- tr[, .(ID, omicID, p_train = get("Pr(>|t|)"))]
  m   <- merge(te2, tr2, by = c("ID", "omicID"))
  thr <- 0.05 / nrow(m)                             # Bonferroni over all tested pairs
  repE <- m[is.finite(p_train) & is.finite(p_test) & p_train < thr & p_test < thr]
  if (!nrow(repE)) { message("  ", covarType, ": no replicated pairs; skipping"); return(NULL) }

  sigE <- repE[, .(ID, EntrezGeneSymbol = omicID, Estimate)]
  wide <- sigE %>%
    pivot_wider(names_from = ID, values_from = Estimate) %>%
    select(where(~ sum(!is.na(.)) >= 3))           # keep exposures with >=3 proteins
  exposure_cols <- setdiff(names(wide), "EntrezGeneSymbol")
  if (!length(exposure_cols)) { message("  ", covarType, ": no exposures with >=3 proteins"); return(NULL) }

  merged <- list(wide,
                 heritage_eff[, c("EntrezGeneSymbol", "HERITAGE_effect")],
                 glp1_step1_eff[, c("EntrezGeneSymbol", "GLP1_effect1")],
                 glp1_step2_eff[, c("EntrezGeneSymbol", "GLP1_effect2")]) %>%
    reduce(full_join, by = "EntrezGeneSymbol") %>%
    left_join(prot_rel[, c("EntrezGeneSymbol", "r_cross")], by = "EntrezGeneSymbol") %>%
    mutate(r_cross = ifelse(is.na(r_cross), mean(r_cross, na.rm = TRUE), r_cross),
           w_rel   = pmax(r_cross, 0))

  cor_rows <- map_dfr(exposure_cols, function(ec) map_dfr(INTERVENTIONS, function(ic) {
    res <- wtd_cor_and_p(merged[[ec]], merged[[ic]], merged$w_rel)
    tibble(covarType = covarType, exposure_id = ec, intervention = ic,
           r = unname(res["r"]), pval = unname(res["p"]), n_eff = unname(res["neff"]))
  }))
  # FDR (Benjamini-Hochberg) PER TRIAL, over the INTERPRETABLE family only: a
  # correlation with effective N < MIN_NEFF is unstable (driven by 2-3 proteins),
  # so we exclude it from the multiple-testing family rather than let it inflate
  # the family size. pval_BH stays NA for those (and for non-finite p).
  cor_rows <- cor_rows %>%
    group_by(intervention) %>%
    mutate(pval_BH = {
      pb   <- rep(NA_real_, dplyr::n())
      keep <- is.finite(n_eff) & n_eff >= MIN_NEFF & is.finite(pval)
      if (any(keep)) pb[keep] <- p.adjust(pval[keep], method = "BH")
      pb
    }) %>% ungroup() %>%
    group_by(exposure_id) %>% mutate(sig_any = any(pval_BH < 0.05, na.rm = TRUE)) %>% ungroup()

  # carry Eid + Category alongside the model-term ID: the MR annotation
  # (annotate_mr.R) keys on MRmotifs `Exposure`, which matches the per-level term
  # (ID) for categorical exposures but the base Eid for ordered-factor exposures,
  # so the downstream join needs both keys.
  scatter <- repE[, .(covarType = covarType, exposure_id = ID, Eid = Eid,
                      Category = Category,
                      protein = omicID, beta_HEAP = Estimate, se_HEAP = `Std. Error`)] %>%
    left_join(heritage_eff,   by = c("protein" = "EntrezGeneSymbol")) %>%
    left_join(glp1_step1_eff, by = c("protein" = "EntrezGeneSymbol")) %>%
    left_join(glp1_step2_eff, by = c("protein" = "EntrezGeneSymbol")) %>%
    left_join(prot_rel[, .(protein = EntrezGeneSymbol, olink_soma_r = r_cross)], by = "protein")

  message(sprintf("  %s: %d exposures x 3 interventions; %d sig-any exposures; scatter rows=%d",
                  covarType, length(exposure_cols),
                  length(unique(cor_rows$exposure_id[cor_rows$sig_any])), nrow(scatter)))
  list(cor = cor_rows, scatter = scatter)
}

results <- lapply(covarTypes, run_one)
results <- results[!vapply(results, is.null, logical(1))]
if (!length(results)) stop("No covarType produced results.")

cor_all     <- rbindlist(lapply(results, `[[`, "cor"), fill = TRUE)
scatter_all <- rbindlist(lapply(results, `[[`, "scatter"), fill = TRUE)

# --- write outputs -----------------------------------------------------------
out_dir  <- heap_project_output("support", "intervention_compare")
fig_data <- file.path(heap_project_root("figures"), "data")
dir.create(out_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(fig_data, recursive = TRUE, showWarnings = FALSE)
fwrite(cor_all,     file.path(out_dir,  "intervention_correlations.tsv"), sep = "\t")
fwrite(scatter_all, file.path(out_dir,  "intervention_scatter.tsv"),      sep = "\t")
fwrite(cor_all,     file.path(fig_data, "intervention_correlations.tsv"), sep = "\t")
fwrite(scatter_all, file.path(fig_data, "intervention_scatter.tsv"),      sep = "\t")

message("\nWrote:\n  ", file.path(out_dir, "intervention_correlations.tsv"),
        "\n  ", file.path(out_dir, "intervention_scatter.tsv"),
        "\n  (+ copies under figures/data/)\nDONE.")
