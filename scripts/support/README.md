# Supporting analyses

Analyses that run on the outputs of the main modules and produce the tables behind specific
figure panels and supplementary items. Each script's header describes its inputs and
outputs. Where an order is shown (→), run the scripts in that order.

## Fig. 1: variance decomposition

| Script | Does |
|---|---|
| `module1_greml_vs_r2.R` | Compares GREML variance components with HEAP predictive R², protein by protein |

## Fig. 3: mediation

| Script | Does |
|---|---|
| `module3_disease_specificity.R` → `module3_category_disease_full.R` → `module3_intermediaries_forest.R` | Separates proteins that mediate many diseases from disease-specific ones; counts mediating proteins per exposure category and disease; indirect effects with 95% CIs for the forest plot |
| `module3_mediation_x_variance.R` | Joins, per exposure category and protein, the variance the category explains (Module 1) with the indirect effect on disease (Module 3) |

## Fig. 4: Mendelian randomization

| Directory | Does |
|---|---|
| [`mr_tables/`](mr_tables) | `build_mr_tables.R` (per-arm tables, tiers, triad motifs) → `compare_arms.R` (replication across the UK Biobank and deCODE arms) |
| [`coloc/`](coloc/README.md) | Colocalization of cis-pQTL edges, and the tables for the regional plots |

## Fig. 5: intervention trials

| Script | Does |
|---|---|
| `intervention_compare/run_intervention_compare.R` → `build_mr_tiered_pd.R` → `build_mr_protein_class.R` → `annotate_mr.R` | Compares exposure-associated proteins with exercise and semaglutide trial proteomics, and labels them with their MR evidence |
| `build_shared_language_network.R` | Node and edge tables for the network of exposures, trials and shared proteins |

## Fig. 6: proteome-based exposure scores (PES)

| Script | Does |
|---|---|
| `module6_holdout_ci.R`, `module6_within_ci.R` | Held-out accuracy and within-person change, with confidence intervals from the held-out participants |
| `module6_visit_timing.R` | Years between repeat visits |
| `module6_covariate_sensitivity.R` | Held-out accuracy of covariates alone vs. covariates + PES under each covariate set |
| `module6_pes_score_correlation.R` | Correlations between the PES of different exposures |
| `module6_pes_disease_scale.R` | Every PES against every disease (Cox models, full cohort) |
| `module6_pes_disease_dc_holdout.R` → `module6_pes_disease_dc_holdout_ci.R` → `collect_dc_holdout_ci.R` | Held-out gain in C-index from adding each PES, for every exposure–disease pair, with bootstrap CIs |
| `module6_quadrant_scan.R` → `module6_quadrant_ladders.R` → `module6_pes_disease_ladder_slate.R` → `module6_pes_disease_bootstrap.R` → `module6_pes_disease_ladder_annotate.R` | Held-out C-index ladders (covariates, + PES, + exposure) for selected pairs, with bootstrap CIs |
| `module6_disease_covars.R` | Shared: the Cox covariate frame used by the disease scripts |

## Other

| Script | Does |
|---|---|
| `build_exposure_level_labels.R` | Readable labels for each level of the ordinal exposures (dose–response figures) |

The SLURM launchers for the longer Module 6 jobs are in [`slurm/module6/`](../../slurm/module6).
