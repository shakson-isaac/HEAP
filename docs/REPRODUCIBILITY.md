# Reproducing the analyses

Re-running HEAP needs an approved UK Biobank application and a cluster. To try the code
without either, run the synthetic-data [demo](../demo/README.md). To use the published
results, go to [heap.bio](https://heap.bio) or the paper's Supplementary Tables and Data;
they are outputs of this code and cannot be used to re-run it.

Software: [`ENVIRONMENT.md`](ENVIRONMENT.md). Data sources: [`DATA_ACCESS.md`](DATA_ACCESS.md).
Stage inputs and outputs, machine-readable: `config/io_map.yml`.

## The analysis dataset

Every module reads one object, `HEAP.rds`, plus per-protein genetic scores. How you
build it from your UK Biobank download depends on how your lab stores those data; the
loader (`scripts/loaders/HEAP_loader.R`) shows how we built ours. All column names come
from `config/`, so the same object can be rebuilt under any application.

| Element | Contents | Used by |
|---|---|---|
| `prot_baseline` | `eid` + one column per protein (Olink NPX, baseline visit) | all modules |
| `E_baseline` | Named list, one data frame per exposure category: `eid` + one column per exposure, named as in `config/exposure_sets/analysis_exposures.tsv` | all modules |
| `covars_baseline` | `eid` + covariates, named as in `config/covariates/covariate_sets.yml` | all modules |
| `protIDs`, `Elist_names`, `Eid_cat`, `ordinalIDs`, `covars_list` | Protein list, category names, exposure → category map, ordinal exposures, covariate names | all modules |
| `prot_long`, `E_long`, `covars_long`, `split_df` | The same with an `instance` (visit) column; visits per participant | Module 6 |
| `disease`, `prevalent_disease` | Age at first report of each incident disease; prevalent-disease flags | Modules 3, 6 |

Genetic scores are one PLINK `.sscore` file per protein (`#IID`, `SCORE1_AVG`) in
`$HEAP_IGLOO_ROOT/UKB/ProtGScis/` and `ProtGStrans/`, named by OmicsPred ID;
`UKB/OMICSPRED/UKB_Olink_multi_ancestry_models_val_results_portal.csv` maps protein
(`Gene`) to `OMICSPRED_ID`. `demo/make_synthetic_heap.R` writes all of these in this format.

## Setup

```bash
export HEAP_ROOT=/path/to/this/repo
export HEAP_SCRATCH_ROOT=/path/to/scratch
export HEAP_IGLOO_ROOT=/path/to/output/root     # outputs: $HEAP_IGLOO_ROOT/UKB/HEAP/output
```

The `base` covariate set is the primary analysis; the other sets in
`config/covariates/covariate_sets.yml` and the filters in `config/samples/sample_filters.yml`
are sensitivity analyses.

## Run order

```
Loader (HEAP.rds)
 ├─ Genetic scores ──┐
 ├─ Exposure GWAS ───┼─ LD score regression
 │                   ├─ Module 1 ─→ Module 3
 │                   ├─ Module 2 ─→ enrichment
 │                   └─ MR (also needs exposure GWAS + pQTL)
 ├─ GREML
 └─ Module 6
```

**1. Protein genetic scores**

```bash
Rscript scripts/genetic_scores/OMICSPredCisTrans.R
sbatch slurm/genetic_scores/UKBgen_ProtGS.sh
sbatch slurm/genetic_scores/UKBgen_ProtGScis.sh
sbatch slurm/genetic_scores/UKBgen_ProtGStrans.sh
```

**2. Analysis dataset.** `sbatch slurm/loaders/run_HEAP_loader.sh`, then
`Rscript workflow/regenerate_exposure_manifest.R` to refresh the exposure list.

**3. Analysis modules.** Each named experiment in `config/modules/<module>_experiments.yml`
fixes the covariate set, model family and sample filter; the launcher builds a job
manifest and submits the array:

```bash
EXPERIMENT=M1_base_lasso bash slurm/module1/HEAPmodule1_manifest.sh
```

| Analysis (code module) | Launcher | Needs |
|---|---|---|
| Variance decomposition (module 1) | `slurm/module1/HEAPmodule1_manifest.sh` | `HEAP.rds` |
| Associations (module 2) | `slurm/module2/HEAPmodule2_manifest.sh`; then `scripts/module2_associations/summarize_replicated_associations.R` | `HEAP.rds` |
| Mediation (module 3) | `slurm/module3/HEAPmodule3_manifest.sh` | the matching Module 1 experiment |
| Enrichment (module 4) | `slurm/module4_enrichment/run_gsea.sh` | Module 2 |
| MR (module 5) | `slurm/module5/HEAPmodule5_manifest.sh` | exposure GWAS |
| PES (module 6) | see below | `HEAP.rds` |

Module 6:

```bash
Rscript scripts/module6_prediction/Module6_config.R
bash slurm/module6/submit_longitudinal_workflow.sh pilot   # one exposure, smoke test
bash slurm/module6/submit_longitudinal_workflow.sh all
sbatch slurm/module6/Module6_compact_pes_array.sh
sbatch slurm/module6/Module6_longitudinal_pes_frozenrisk_array.sh
```

More detail: `scripts/module4_enrichment/README.md`, `scripts/module5_mr/README.md`,
`slurm/module5/README.md`.

**4. Exposure GWAS** (excludes the proteomics cohort):

```bash
Rscript scripts/gwas_regenie/prepare_gwas_exposures.R
sbatch --array=1-$(wc -l < slurm/gwas_regenie/evars_continuous_heap.txt) \
       slurm/gwas_regenie/gwas_regenie_exposures_continuous_v2.sh
sbatch --array=1-$(wc -l < slurm/gwas_regenie/evars_binary_heap.txt) \
       slurm/gwas_regenie/gwas_regenie_exposures_binary_v2.sh
```

Re-run `prepare_gwas_exposures.R` whenever `analysis_exposures.tsv` changes. Detail:
`slurm/gwas_regenie/README_gwas_exposure_workflow.md`.

**5. LD score regression.** `slurm/ldsc/README_ldsc_workflow.md`; collect with
`scripts/ldsc/collect_ldsc_h2.R` and `collect_ldsc_rg.R`.

**6. GREML.** `RUN_GROUP=heap_v1 bash slurm/greml/submit_greml_workflow.sh`

## Figures

```bash
bash scripts/visualizations/make_figures.sh          # all main + supplementary figures
bash scripts/visualizations/make_figures.sh fig4     # one main figure
```

[`scripts/visualizations/README.md`](../scripts/visualizations/README.md) maps each panel
to its script.

## Supplementary Tables and Data

Each item is computed from the module outputs by the scripts below (paths relative to
`scripts/`). Formatting into the workbook and zip is done in the manuscript build and
changes no value.

| Results | Script(s) |
|---|---|
| Variance decomposition, train/test stability, specification grid | `analysis_summaries/summarize_module1_component_r2.R`, `summarize_module1_variance.R`, `summarize_module1_spec_grid.R` |
| GREML | `analysis_summaries/summarize_greml_base.R` |
| Associations and GxE, per-specification deposits | `analysis_summaries/summarize_module2_supp_tables.R`, `export_assoc_spec_deposit.R`, `module2_spec_sensitivity.R`, `module2_gxe_spec_sensitivity.R` |
| Mediation, disease list, per-specification deposits | `analysis_summaries/summarize_module3_mediation.R`, `summarize_module3_disease_list.R`, `export_med_spec_deposit.R`, `module3_spec_sensitivity.R` |
| Enrichment | `module4_enrichment/04_enrichment_tables.R`, `run_module1_*_gsea.R` |
| MR and colocalization | `analysis_summaries/build_mr_supp_tables.R`, `summarize_mr_triads.R`, `module5_mr/Module5_deCODE.R`, `support/coloc/run_coloc_systematic.R` |
| Intervention comparison | `support/intervention_compare/run_intervention_compare.R` |
| PES accuracy and disease prediction | `support/module6_holdout_ci.R`; `support/module6_pes_disease_ladder_slate.R` → `module6_pes_disease_bootstrap.R` → `module6_pes_disease_ladder_annotate.R`; `support/module6_pes_disease_dc_holdout.R`, `module6_pes_disease_dc_holdout_ci.R`, `collect_dc_holdout_ci.R`; `analysis_summaries/export_pes_lookup_deposits.R` |
| PES weights | `module6_prediction/Module6_export_deployable_lasso_weights.R` |

## Runtime

Planning estimates derived mostly from the SLURM time limits each stage requests, not
benchmarked runs. Scale: 2,686 proteins, 169 exposures; GWAS over 91 binary and 81
continuous exposures.

| Stage | CPU-hours (≈) | Wall-clock (≈) |
|---|---|---|
| Loader | 2 | 2 h (one job, ~100 GB) |
| Genetic scores | 28,000 | 2–4 days |
| Exposure GWAS | 10,000 | 7 days |
| Module 1 | 6,000 | 0.5–1 day |
| Module 2 (+ sensitivity) | 10,000 | 0.5–1 day |
| Module 3 | 5,000 | 0.5–1 day |
| MR | 3,000 | 0.5 day |
| Module 6 | 3,000 | 1–2 days |
| GREML | 35,000 | 2–4 days |
| Enrichment | < 5 | minutes |

With 100–300 concurrent job slots, an end-to-end run takes about 1–2 weeks.
