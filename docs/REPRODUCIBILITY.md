# Reproducibility

How to reproduce the HEAP analyses from UK Biobank source data. This is the single
run guide: prerequisites, stage order, entry points, and expected runtimes.

- **Getting the data** → [`DATA_ACCESS.md`](DATA_ACCESS.md)
- **Software and cluster environment** → [`ENVIRONMENT.md`](ENVIRONMENT.md)
- **Machine-readable stage inputs/outputs and completion checks** → `config/io_map.yml`

## 0. Configure paths

No analysis path is hardcoded. Set these, then source `workflow/00_paths.R` from R:

```bash
export HEAP_ROOT=/path/to/this/repo
export HEAP_SCRATCH_ROOT=/path/to/scratch
export HEAP_IGLOO_ROOT=/path/to/output/root
```

`00_paths.R` also prepends a shared R library to `.libPaths()` and sets `umask 002`
so outputs are group-writable. Canonical outputs land under
`$HEAP_IGLOO_ROOT/UKB/HEAP/output`, referred to below as `$OUT`.

Install R dependencies with `Rscript scripts/setup/install_r_packages.R`.

## Covariate scheme

`base` is the primary adjustment for all main results. The remaining sets
(`base_bmi`, `base_draw`, `base_clinical`, `base_ses`, `base_prevalent`) and the
`exclude_prevalent` sample filter are sensitivity analyses. Sets are defined in
`config/covariates/covariate_sets.yml`, sample filters in
`config/samples/sample_filters.yml`.

## 1. Protein genetic scores

Split OmicsPred score files into cis/trans, then build the score inputs read by
Modules 1 and 2:

```bash
Rscript scripts/genetic_scores/OMICSPredCisTrans.R
sbatch slurm/genetic_scores/UKBgen_ProtGS.sh
sbatch slurm/genetic_scores/UKBgen_ProtGScis.sh
sbatch slurm/genetic_scores/UKBgen_ProtGStrans.sh
```

## 2. Foundation loader — everything depends on this

`scripts/loaders/HEAP_loader.R` is the single unified loader. It reads raw UK Biobank
data and writes the canonical `HEAP.rds` consumed by every module; no module reads raw
UK Biobank data directly.

```bash
sbatch slurm/loaders/run_HEAP_loader.sh
```

After a loader rerun, refresh the exposure manifest with
`Rscript workflow/regenerate_exposure_manifest.R`, which updates
`config/exposure_sets/analysis_exposures.tsv` — the single source of truth for the
exposure set and their types.

Modules obtain data from `HEAP.rds` via `as_pxs_baseline(heap)` (Modules 1, 2, 3) or
`as_pxs_longitudinal(heap)` (Module 6), both defined in `00_paths.R`.

## 3. Analysis modules

Modules 1, 2, 3 and 5 share one submission pattern. Each named experiment is defined in
`config/modules/<module>_experiments.yml`; the launcher generates a manifest from it,
sizes the job array from the manifest's row count, selects wall time, memory and cores
from the model family recorded there, and submits:

```bash
EXPERIMENT=M1_base_lasso bash slurm/module1/HEAPmodule1_manifest.sh
```

The experiment name encodes the covariate set, model family and sample filter. Outputs
land in `$OUT/<output_subdir>/<experiment>/`.

| Stage | R entry point | Launcher | Depends on |
|---|---|---|---|
| Module 1 — variance decomposition | `module1_variance_decomposition/Module1_suggested.R` | `slurm/module1/HEAPmodule1_manifest.sh` | `HEAP.rds` |
| Module 2 — E/G/GxE associations | `module2_associations/Module2.R` | `slurm/module2/HEAPmodule2_manifest.sh` | `HEAP.rds` |
| Module 2 — replicated summary | `module2_associations/summarize_replicated_associations.R` | — | Module 2 |
| Module 3 — mediation | `module3_mediation/Module3.R` | `slurm/module3/HEAPmodule3_manifest.sh` | **Module 1** |
| Module 4 — enrichment | `module4_enrichment/` (`00`–`04`) | `slurm/module4_enrichment/run_gsea.sh` | Module 2 |
| Module 5 — Mendelian randomization | `module5_mr/Module5.R` | `slurm/module5/HEAPmodule5_manifest.sh` | Exposure GWAS (§4) |
| Module 6 — proteomic exposure scores | `module6_prediction/Module6_prod_longitudinal.R` | `slurm/module6/submit_longitudinal_workflow.sh` | `HEAP.rds` |

R entry points are relative to `scripts/`.

Module 6 does not use the manifest pattern. Build its exposure specs first, then run the
pilot before committing to the full array:

```bash
Rscript scripts/module6_prediction/Module6_config.R
bash slurm/module6/submit_longitudinal_workflow.sh pilot   # one exposure, smoke test
bash slurm/module6/submit_longitudinal_workflow.sh all     # full array
sbatch slurm/module6/Module6_compact_pes_array.sh
sbatch slurm/module6/Module6_longitudinal_pes_frozenrisk_array.sh
```

Module 3 guards its upstream `mediation_scores` and fails cleanly if the matching
Module 1 experiment (same covariate set and family) has not been run.

Stage detail: `scripts/module4_enrichment/README.md`, `scripts/module5_mr/README.md`,
and `slurm/module5/README.md`.

## 4. Exposure GWAS (regenie)

Two-sample design that excludes the proteomics cohort; covariate set `base`.

```bash
Rscript scripts/gwas_regenie/prepare_gwas_exposures.R
sbatch --array=1-$(wc -l < slurm/gwas_regenie/evars_continuous_heap.txt) \
       slurm/gwas_regenie/gwas_regenie_exposures_continuous_v2.sh
sbatch --array=1-$(wc -l < slurm/gwas_regenie/evars_binary_heap.txt) \
       slurm/gwas_regenie/gwas_regenie_exposures_binary_v2.sh
```

Re-run `prepare_gwas_exposures.R` whenever `analysis_exposures.tsv` changes — the
exposure lists and array sizes derive from it. Full detail:
`slurm/gwas_regenie/README_gwas_exposure_workflow.md`.

## 5. LD score regression

SNP-heritability and genetic correlation on the exposure GWAS. See
`slurm/ldsc/README_ldsc_workflow.md`; results are collected with
`scripts/ldsc/collect_ldsc_h2.R` and `collect_ldsc_rg.R`.

## 6. Population architecture (GREML)

Self-contained three-stage chain: build the master unrelated set, run the per-protein
GREML array, then summarize.

```bash
export RUN_GROUP=heap_v1
bash slurm/greml/submit_greml_workflow.sh     # SPEC=base
```

Outputs: `$OUT/population_architecture/<SPEC>/grm_cutoff_<label>/`.

## Dependency order

```
Loader (HEAP.rds)
 ├─ Genetic scores ──┐
 ├─ Exposure GWAS ───┼─ LD score regression
 │                   ├─ Module 1 ─→ Module 3
 │                   ├─ Module 2 ─→ Module 4 enrichment
 │                   └─ Module 5 (also needs exposure GWAS + pQTL)
 ├─ GRM master ─→ Population architecture (GREML)
 └─ Module 6
```

## Scale and runtime

Planning estimates, not guarantees. Scale: 2,686 proteins and 169 exposomic features
in the main analysis; the exposure GWAS runs over 91 binary and 81 continuous exposures
(array sizes are derived from the `evars_*` lists at submission time). Per-protein
modules run as batched arrays sized from the generated manifest.

| Stage | CPU-hours (≈) | Wall-clock at moderate concurrency (≈) |
|---|---|---|
| Loader | 2 | ~2 h (single job, ~100 GB) |
| Genetic scores (cis + trans + total) | ~28,000 | 2–4 days |
| Exposure GWAS | ~10,000 | ~7 days (continuous is the long pole) |
| Module 1 | ~6,000 | 0.5–1 day |
| Module 2 (+ sensitivity) | ~10,000 | 0.5–1 day |
| Module 3 | ~5,000 | 0.5–1 day |
| Module 5 (MR) | ~3,000 | ~0.5 day |
| Module 6 | ~3,000 | 1–2 days |
| Population architecture (GREML) | ~35,000 | 2–4 days (heaviest array) |
| Enrichment | < 5 | minutes |

Wall-clock ≈ (per-task time) × (array size) / (concurrent slots). On a cluster with
100–300 concurrent slots, an end-to-end run is on the order of **1–2 weeks**,
dominated by the genetic scores, exposure GWAS, and GREML arrays.

## Contact

Questions: shakson_isaac@g.harvard.edu
