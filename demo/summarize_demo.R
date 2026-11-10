#!/usr/bin/env Rscript
# demo/summarize_demo.R -- compare each demo stage's estimates with the effects
# planted by make_synthetic_heap.R. Stages that were not run are skipped.
# Usage: Rscript demo/summarize_demo.R <demo_dir>
suppressPackageStartupMessages(library(data.table))
D   <- commandArgs(trailingOnly = TRUE)[1]
out <- file.path(D, "igloo", "UKB", "HEAP", "output")
pl  <- fread(file.path(D, "planted_effects.tsv"))
x1  <- pl$exposure_x1[1]; x2 <- pl$exposure_x2[1]
show <- function(x) print(as.data.frame(x), row.names = FALSE)
section <- function(title, path, fn) {
  if (!file.exists(path)) return(invisible())
  cat("\n=== ", title, " ===\n", sep = ""); fn(path)
}

# ---- Module 1: held-out R2 unique to each block, mean over folds -----------
section("Module 1 -- held-out R2 unique to each block (mean of 5 folds)",
        file.path(out, "module1_predictive_r2_score_partition", "base", "lasso",
                  "predictive_r2_coarse_1.txt"), function(p) {
  m1 <- rbindlist(lapply(list.files(dirname(p), "^predictive_r2_coarse_.*\\.txt$", full.names = TRUE), fread), use.names = TRUE, fill = TRUE)
  u  <- m1[method == "score_unique_drop", .(r2 = mean(r2)), by = .(omic, block)]
  u  <- as.data.frame(dcast(u, omic ~ block, value.var = "r2"))
  u[-1] <- round(u[-1], 3)
  cat("Planted: genetics ~0.09, exposures ~0.07, GxE ~0.02 of protein variance.\n\n")
  show(u[, c("omic", "G", "E", "GxE", "Covars")])
})

# ---- Module 1 (GREML): variance components with GCTA ------------------------
gp <- file.path(D, "greml", "igloo", "UKB", "HEAP", "output", "population_architecture",
                "base", "grm_cutoff_0p025")
section("Module 1 -- GREML variance components (GCTA REML, G + E + GxE kernels)",
        file.path(D, "greml", "planted_greml.tsv"), function(p) {
  fs <- list.files(file.path(gp, "primary"), "^protein_[0-9]+_summary\\.tsv$", full.names = TRUE)
  if (!length(fs)) { cat("No GREML results (stage skipped or failed).\n"); return(invisible()) }
  pg <- fread(p)[1]
  cat(sprintf("Planted: G %.2f, E %.2f, GxE %.2f of protein variance.\n", pg$planted_G, pg$planted_E, pg$planted_GxE))
  rs <- file.path(gp, "..", "relatedness", "grm_cutoff_0p025", "master_relatedness_summary.tsv")
  if (file.exists(rs)) { r <- fread(rs); cat("Relatedness cutoff 0.025:\n"); show(r) }
  g <- rbindlist(lapply(fs, fread), use.names = TRUE, fill = TRUE)
  f <- function(v, se) sprintf("%.2f (%.2f)", v, se)
  cat("\nEstimate (SE). GxE is weakly identified at this sample size.\n\n")
  show(g[, .(protein, n, G = f(variance_G, se_G), E = f(variance_E, se_E),
             GxE = f(variance_GxE, se_GxE), converged)])
})

# ---- Module 2: planted exposures vs all others ------------------------------
section("Module 2 -- associations (held-out test set)",
        file.path(out, "module2", "base", "univar_assoc_1.rds"), function(p) {
  fs <- list.files(dirname(p), "^univar_assoc_.*\\.rds$", full.names = TRUE)
  m2 <- lapply(fs, function(f) readRDS(f)$test)
  flat <- function(d) {        # Module 2 can store numeric columns as list columns
    d <- as.data.table(d)
    for (c in names(d)) if (is.list(d[[c]]))
      set(d, j = c, value = vapply(d[[c]], function(z) if (length(z)) as.numeric(z[[1]]) else NA_real_, 0))
    d }
  E  <- rbindlist(lapply(m2, function(x) flat(x[[1]])), use.names = TRUE, fill = TRUE)
  setnames(E, "Pr(>|t|)", "p")
  F  <- rbindlist(lapply(m2, function(x) flat(x[[4]])), use.names = TRUE, fill = TRUE)
  pick <- function(v) E[ID == v, .(protein = omicID, exposure = v,
                                   estimate = round(Estimate, 3), p = signif(p, 2))]
  cat("Planted: ", x1, " 0.25 SD; ", x2, " 0.15 SD.\n\n", sep = "")
  show(rbind(pick(x1), pick(x2)))
  null <- E[!ID %in% c(x1, x2)]
  cat(sprintf("\nAll other %d exposure terms x %d proteins: %.1f%% at p < 0.05 (~5%% expected by chance).\n",
              uniqueN(null$ID), uniqueN(null$omicID), 100 * mean(null$p < 0.05, na.rm = TRUE)))
  cat("\nPlanted cis-score x ", x1, " interaction (0.15 SD):\n\n", sep = "")
  show(F[ID == x1, .(protein = omicID, p_GcisxE = signif(p_GcisxE, 2))])
})

# ---- Module 3: indirect (mediated) effect of the exposome score -------------
section("Module 3 -- mediation of exposure -> disease through each protein",
        file.path(out, "module3", "base", "lasso", "primary_total", "MDres_1.txt"), function(p) {
  m3 <- rbindlist(lapply(list.files(dirname(p), "^MDres_.*\\.txt$", full.names = TRUE), fread), use.names = TRUE, fill = TRUE)
  m3 <- m3[predictor == "PXS_total" & effect_type == "NIE"]
  m3[, disease := ifelse(DZ_ID == pl$disease_mediated[1], "disease 1 (mediated)", "disease 2 (null)")]
  cat("Planted: only ", pl$mediator_protein[1], " lies on the path to disease 1 (log-HR ",
      pl$planted_protein_logHR[1], " per SD).\n\n", sep = "")
  show(m3[order(disease, delta_p), .(disease, protein = protID, NIE_logHR = round(effect_logHR, 3),
                                     p = signif(delta_p, 2), cases = n_cases)])
})

# ---- MR: tiers of the planted edges and the canonical motif counts ---------
section("MR -- tiers of the planted edges (Module5.R on synthetic summary statistics)",
        file.path(out, "mr_edges", "summary", "mr_tiered_edges.tsv"), function(p) {
  planted <- c(paste0("exposure_", 1:4), paste0("protein_", 1:4))
  te <- fread(p)[dataset == "UKB" & src_id %in% planted &
                 edge_dir %in% c("E_to_P", "Pcis_to_D", "Pcis_to_E", "E_to_D")]
  cat("Triad k = exposure_k -> protein_k -> disease_k.  1: true mediator; 2: reverse\n",
      "causation; 3: Steiger direction unresolved; 4: fails colocalization.\n\n", sep = "")
  show(te[order(edge_dir, src_id), .(edge = edge_dir, from = src_id, to = tgt_id,
                                     tier = mr_tier_final, reason = tier_reason)])
  mc <- file.path(D, "repo", "docs", "manuscript_stats", "module5", "mr_motif_counts.tsv")
  if (file.exists(mc)) {
    cat("\nTriad motifs (Tier 1 edges in either pQTL arm; Tier 1+ = both arms):\n\n")
    show(fread(mc)[, .(motif, tier1_triads, tier1plus_triads)])
  }
})

# ---- Module 6: proteome-based exposure score on held-out repeat visits -----
m6p <- file.path(out, "module6_pes_longitudinal", "base", paste0("PESlong_base_", x1, "_"))
section("Module 6 -- proteome-based exposure score (held-out participants)",
        paste0(m6p, "HoldoutByVisitMetrics.tsv"), function(p) {
  cat("Trained on baseline-only participants; scored on ", pl$n_holdout[1],
      " others at baseline (0) and repeat visits (2, 3).

", sep = "")
  h <- fread(p)[model %in% c("prot_only", "cov_only")]
  show(h[, .(visit = instance, model, n, r2 = round(r2, 3))])
  w <- fread(paste0(m6p, "WithinPersonChangeMetrics.tsv"))[model == "prot_only"]
  cat("\nWithin-person change: does the score move with the exposure between visits?\n\n")
  show(w[, .(visits = paste0("0 -> ", followup_instance), n, r = round(delta_cor, 2),
             p = signif(p_delta_pes, 2))])
})
