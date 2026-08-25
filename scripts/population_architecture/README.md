# Population Architecture Pipeline

This module is separate from the existing predictive `Prot_ExPGS` / mediation workflow. It estimates population architecture quantities for proteins using a multi-kernel REML model and keeps predictive out-of-fold scores strictly as descriptive comparison targets.

## What It Reports

- `variance_G`: variance proportion from the genotype GRM.
- `variance_E`: variance proportion from the exposure similarity kernel.
- `variance_GxE`: variance proportion from the Hadamard interaction kernel `K_G ∘ K_E`.
- `variance_Covars_fixed`: variance share of fitted fixed covariates, labeled separately from variance components.
- `variance_CovarsxG_sensitivity` and `variance_CovarsxE_sensitivity`: optional sensitivity kernels, only in the sensitivity model.

## Key Inputs

- `UKB_PGS_PXS_load.rds` for staged proteins, exposures, and covariates.
- `UKBallchr.pgen` for the full genotype panel with real UKB IDs.
- `ukb_nonimputed_snps.pvar` as the LD-pruned SNP list used to build the genotype GRM.
- Existing OOF component files under `Data/Parallel/Module1/<Type>` for descriptive predictive `delta-R2` comparison.

## Main Scripts

- `scripts/export_architecture_inputs.R`
- `scripts/build_genotype_grm.R`
- `scripts/build_environment_kernels.R`
- `scripts/run_population_architecture.R`
- `scripts/summarize_population_architecture.R`

## Cohort Policy

- Primary runs use one maximal proteomics-eligible master cohort defined by genotype availability plus complete fixed covariates for the selected covariate specification.
- `G`, `E`, and `GxE` are built once on that master cohort and reused across proteins.
- Each protein is then fit on its own non-missing phenotype subset inside the shared cohort.
- This keeps proteins comparable and avoids rebuilding a new genotype GRM for every outcome.
- The default recommendation is to build the master GRM on the proteomics-eligible cohort, not the full `~500k` UKB.

## Example Pilot Commands

```bash
module load gcc/14.2.0
module load R/4.5.2

cd /n/groups/patel/shakson_ukb/UK_Biobank

Rscript population_architecture/scripts/export_architecture_inputs.R \
  population_architecture/config/default_config.R pilot_type5 Type5 \
  --proteins=DKKL1,TLR3,LILRB5 \
  --max-samples=500 \
  --seed=1 \
  --force=true

Rscript population_architecture/scripts/build_genotype_grm.R \
  population_architecture/config/default_config.R pilot_type5 Type5

Rscript population_architecture/scripts/build_environment_kernels.R \
  population_architecture/config/default_config.R pilot_type5 Type5 \
  --center-exposures=true \
  --include-covar-kernels=true

Rscript population_architecture/scripts/run_population_architecture.R \
  population_architecture/config/default_config.R pilot_type5 Type5 \
  --model=primary \
  --proteins=DKKL1,TLR3,LILRB5 \
  --center-exposures=true \
  --min-protein-n=450 \
  --continue-on-error=true \
  --force=true

Rscript population_architecture/scripts/run_population_architecture.R \
  population_architecture/config/default_config.R pilot_type5 Type5 \
  --model=sensitivity \
  --proteins=TLR3 \
  --center-exposures=true \
  --min-protein-n=450 \
  --continue-on-error=true

Rscript population_architecture/scripts/summarize_population_architecture.R \
  population_architecture/config/default_config.R pilot_type5 Type5 \
  --model=primary \
  --center-exposures=true \
  --min-protein-n=450
```

## Example Full Run Commands

Single-command Slurm launch:

```bash
cd /n/groups/patel/shakson_ukb/UK_Biobank

sbatch \
  --export=ALL,RUN_ID=full_type3_centered,SPEC=Type3,CENTER_EXPOSURES=true,MIN_PROTEIN_N=2000,RUN_SENSITIVITY=true \
  population_architecture/slurm/population_architecture_full.sh
```

Focused run for the requested 7 proteins plus 2 benchmark extras (`TLR3`, `DKKL1`):

```bash
cd /n/groups/patel/shakson_ukb/UK_Biobank

sbatch population_architecture/slurm/population_architecture_focus9.sh
```

Preferred array-based focused run:

```bash
cd /n/groups/patel/shakson_ukb/UK_Biobank

bash population_architecture/slurm/submit_population_architecture_focus9_array.sh
```

This submits:

- one prep job to export the master cohort and build shared `G`, `E`, and `GxE` kernels
- one job array with one task per protein
- one finalize job to rebuild summaries from the per-protein row outputs

Manual step-by-step launch:

```bash
module load gcc/14.2.0
module load R/4.5.2

cd /n/groups/patel/shakson_ukb/UK_Biobank

Rscript population_architecture/scripts/export_architecture_inputs.R \
  population_architecture/config/default_config.R full_type3_centered Type3

Rscript population_architecture/scripts/build_genotype_grm.R \
  population_architecture/config/default_config.R full_type3_centered Type3

Rscript population_architecture/scripts/build_environment_kernels.R \
  population_architecture/config/default_config.R full_type3_centered Type3 \
  --center-exposures=true \
  --include-covar-kernels=true

Rscript population_architecture/scripts/run_population_architecture.R \
  population_architecture/config/default_config.R full_type3_centered Type3 \
  --model=primary \
  --center-exposures=true

Rscript population_architecture/scripts/run_population_architecture.R \
  population_architecture/config/default_config.R full_type3_centered Type3 \
  --model=sensitivity \
  --center-exposures=true \
  --continue-on-error=true

Rscript population_architecture/scripts/summarize_population_architecture.R \
  population_architecture/config/default_config.R full_type3_centered Type3 \
  --model=primary \
  --center-exposures=true
```

## Notes

- `Type3` is recommended as the primary full-run covariate specification, with `Type5` retained as a richer fixed-effect sensitivity specification.
- The focused protein list for the requested run is stored in `population_architecture/config/protein_sets/focus9_master.txt`.
- The 2 benchmark extras in that focused set are `TLR3` and `DKKL1`, carried forward from earlier smoke tests to give one previously more stable target and one stress-test target.
- The array workflow is the preferred operational mode for focused runs because failures are isolated to individual proteins instead of one monolithic protein loop.
- The `E` kernel now follows the core `Module1.R` exposure preprocessing more closely: columns with missingness above `cfg$exposure_missing_rate_max` are dropped before kernel construction, and `ordinalIDs` are one-hot encoded as categorical features even if stored numerically.
- For the centered vs non-centered exposure sensitivity, rerun the kernel build, REML, and summary steps with `--center-exposures=false` and a different run ID.
- The genotype GRM is built once per exported cohort and then reused across proteins.
