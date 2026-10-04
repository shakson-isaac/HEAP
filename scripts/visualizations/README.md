# Figure code

The code that draws every main and supplementary figure in the HEAP manuscript.

The figure scripts read **completed pipeline outputs**: module results plus a few
summary tables built from them. Without UK Biobank access you cannot regenerate those
outputs, but the summary-level results behind every figure are in the paper's
Supplementary Tables and Supplementary Data, also downloadable from
[heap.bio](https://heap.bio).

## Rebuilding the figures

Run the full pipeline first ([`docs/REPRODUCIBILITY.md`](../../docs/REPRODUCIBILITY.md)),
then the summary scripts below, then:

```bash
export HEAP_PROJECT_ROOT=/path/to/output/root   # the directory holding output/
bash scripts/visualizations/make_figures.sh          # Figs 1-6 + supplement
bash scripts/visualizations/make_figures.sh fig3     # a single main figure
bash scripts/visualizations/make_figures.sh supp     # supplement only
```

Figures are written to `$HEAP_PROJECT_ROOT/figures/`, and the six finished main figures
are also copied to `figures/manuscript/Fig1.pdf` … `Fig6.pdf`. Each supplementary
script writes the exact data it plotted to `figures/data/<module>/<figure>.tsv`.

Besides R, the main figures need [tectonic](https://tectonic-typesetting.github.io) to
compile the TikZ schematics and the vector layouts of Figs 3–4, and `pdftocairo`
(poppler) to rasterize them. Set `TECTONIC=/path/to/tectonic` if it is not on `PATH`.

## Main figures

Each main figure is assembled from panels drawn separately. `make_figures.sh` runs
them in this order.

| Figure | Panels (in build order) | Assembled by |
|---|---|---|
| **Fig. 1** Genetic vs. exposomic spectrum | `figures/fig1_exp_S1.tex`, `figures/fig1_partition_block.tex` → `build_panelA_S1.R`; `figures/cmp_greml_heap_spectrum.R`; `figures/fig_expo_signatures.R` | `build_module1_fig1_horizontal.R` |
| **Fig. 2** Exposure–protein associations | drawn within the builder | `build_module2_fig3_composite.R` |
| **Fig. 3** Mediation | `figures/fig_module3_schematic_compact.tex`; `figures/fig_mediation_scale_main.R`; `figures/fig_mediation_pleiotropy.R`; `figures/fig_mediation_forest.R` | `build_module3_composite_vector.R` (+ tectonic) |
| **Fig. 4** Mendelian randomization | `build_mr_tier_ladder.R`; `mr_schematics/mr_evidence_dag.tex`; `build_mr_panelb_folded.R`; `build_mr_panelc.R`; `build_mr_paneld.R`; `figures/fig_mr_main_mediators_dag.R` | `build_mr_composite_vector.R` (+ tectonic) |
| **Fig. 5** Intervention trials | `figures/fig_m4_schematic_compact.tex`; `figures/fig_m4_panel_b.R`; `figures/fig_m4_panel_d.R`; `figures/fig_m4_shared_network.R` | `build_module4_composite_v4.R` |
| **Fig. 6** Proteome-based exposure scores | `figures/fig_module6_schematic_compact.tex`; `figures/fig_m6_panel_b.R`; `figures/fig_m6_panel_c.R`; `figures/fig_m6_panel_d.R` | `build_module6_composite_landscape.R` |

Panel scripts run with `HEAP_CELL=1`, which writes the panel at its print size for the
composite. The `fig_m4_*` and `module4` names refer to the paper's Module 5
(interventions); see the module-number table in the [top-level README](../../README.md).

## Supplementary figures

One script per figure, `figures/<figure>.R`. `config/figures/figure_registry.tsv`
records each figure's inputs and output location.

| Figure | Script |
|---|---|
| HEAP vs. GREML variance-component concordance | `fig_greml_vs_r2.R` |
| Variance partition per protein | `fig_variance_architecture.R` |
| Module 1 robustness across perturbations | `fig_module1_robustness.R` |
| Train vs. test stability by exposure category | `fig_traintest_stability_categories.R` |
| GxE at the noise floor | `fig_gxe_noise_floor.R` |
| Pathway enrichment by variance component | `fig_pathways_by_component.R` |
| Per-category biology of exposure-responsive proteins | `fig_category_biology.R` |
| Exposure correlation matrix | `fig_exposure_cormatrix.R` |
| Signed ExWAS (Miami) plot | `fig_exwas_miami.R` |
| Dose–response of exposure–protein effects | `fig_dose_response.R` |
| Covariate sensitivity, associations | `fig_module2_spec_sensitivity.R` |
| Covariate sensitivity, GxE associations | `fig_module2_gxe_spec_sensitivity.R` |
| Polygenic GxE architecture | `fig_gxe_architecture.R` |
| Tissue enrichment (per exposure; themes) | `fig_tissue_enrichment.R`, `fig_tissue_themes.R` |
| Pathway enrichment (per exposure; themes) | `fig_pathway_enrichment.R`, `fig_pathway_themes.R` |
| Mediation hubs | `fig_mediation_flows.R` |
| Mediation prioritization | `fig_mediation_prioritization.R` |
| Proportion mediated | `fig_mediation_proportion.R` |
| Covariate sensitivity, mediation | `fig_module3_spec_sensitivity.R` |
| Exposure-GWAS instrument diagnostics | `fig_instrument_diagnostics.R` |
| Exposure-GWAS exemplars | `fig_gwas_exemplars.R` |
| Genetic correlation between exposures (LDSC) | `fig_ldsc_rg.R` |
| MR motif distribution | `fig_mr_motif_overview.R` |
| MR edge attrition | `fig_mr_attrition.R` |
| MR rigor and replication | `fig_mr_rigor.R` |
| MR refines mediation | `fig_mr_refines_mediation.R` |
| Colocalization locus plot | `fig_mr_coloc.R` |
| PES incremental value over covariates | `fig_pes_incremental_value.R` |
| PES robustness to covariate specification | `fig_pes_covariate_robustness.R` |
| PES accuracy vs. panel size | `fig_pes_panel_size.R` |
| PES exposure specificity | `fig_pes_specificity.R` |
| PES within-person tracking design | `fig_pes_tracking_design.R` |
| PES tracking between imaging visits | `fig_pes_imaging_tracking.R` |

## Summary tables the figures read

Some figures read summary tables computed from the module outputs rather than the raw
outputs. Run these once, after the pipeline and before `make_figures.sh`. Each script's
header documents its inputs and arguments.

| Used by | Scripts (in order) |
|---|---|
| Fig. 1, GREML concordance | `scripts/support/module1_greml_vs_r2.R` |
| Fig. 2 | `scripts/analysis_summaries/module2_program_tissue_edges.R` |
| Fig. 3, mediation supp. | `scripts/support/module3_disease_specificity.R` → `module3_category_disease_full.R` → `module3_intermediaries_forest.R`; `module3_mediation_x_variance.R` |
| Fig. 4, MR supp. | `scripts/support/mr_tables/build_mr_tables.R` → `compare_arms.R`; `scripts/support/coloc/` (see its README); `scripts/analysis_summaries/summarize_mr_triads.R` |
| Fig. 5 | `scripts/support/intervention_compare/run_intervention_compare.R` → `build_mr_tiered_pd.R` → `build_mr_protein_class.R` → `annotate_mr.R`; `scripts/support/build_shared_language_network.R` |
| Fig. 6, PES supp. | `scripts/support/module6_holdout_ci.R`, `module6_within_ci.R`, `module6_visit_timing.R`, `module6_covariate_sensitivity.R`, `module6_pes_disease_scale.R`, `module6_pes_score_correlation.R`, `module6_quadrant_scan.R` → `module6_quadrant_ladders.R` → `module6_pes_disease_ladder_slate.R` → `module6_pes_disease_ladder_annotate.R` |
| Sensitivity supp. | `scripts/analysis_summaries/module2_spec_sensitivity.R`, `module2_gxe_spec_sensitivity.R`, `module3_spec_sensitivity.R` |
| GWAS supp. | `scripts/analysis_summaries/summarize_gwas_loci.R` |
| Dose–response | `scripts/support/build_exposure_level_labels.R` |

Fig. 5 also reads three published supplementary tables: HERITAGE exercise proteomics
(Robbins *et al.*, *JCI Insight* 2023), STEP 1/STEP 2 semaglutide proteomics
(Maretty *et al.*, *Nat Med* 2025), and Olink–SomaScan correlations (Eldjarn *et al.*,
*Nature* 2023). Sources are in [`docs/DATA_ACCESS.md`](../../docs/DATA_ACCESS.md).

## Shared code

`common/` holds what every figure script uses: paths (`figure_paths.R`), result loaders
(`load_heap_results.R`), the house theme and exposure-category palette (`plot_theme.R`,
`config/figures/exposure_category_palette.tsv`), exposure labels (`label_helpers.R`),
and figure output (`export_helpers.R`, `figure_registry.R`).
