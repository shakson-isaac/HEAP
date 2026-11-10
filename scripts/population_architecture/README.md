# GREML variance decomposition

Estimates, for each protein, the share of variance explained by genetics (G), the
exposome (E) and their interaction (GxE) with multi-kernel REML in GCTA. This is the
GREML estimate in Fig. 1, next to HEAP's predictive R² from
[`module1_variance_decomposition`](../module1_variance_decomposition).

| Component | Kernel |
|---|---|
| `variance_G` | Genomic relationship matrix from LD-pruned genotyped SNPs |
| `variance_E` | Exposure similarity matrix (exposures standardized; categorical exposures one-hot encoded) |
| `variance_GxE` | Element-wise product of the G and E matrices |

G, E and GxE are fitted jointly. The fixed covariates come from the chosen covariate set
(default `base`); their share is reported separately as `variance_Covars_fixed`.

## Run

```bash
bash slurm/greml/submit_greml_workflow.sh
```

This submits three dependent jobs:

1. `BuildMasterUnrelated.sh`: builds one genomic relationship matrix for every
   participant with proteomics, genotypes and complete covariates, and keeps one member
   of each pair related above `GRM_CUTOFF` (default 0.025).
2. `ProtGremlGrmCutoff.sh`: an array with one task per protein. Each task builds the
   kernels on the unrelated set and fits the model on participants with that protein
   measured (at least `MIN_PROTEIN_N`, default 2,000).
3. `SummaryGremlGrmCutoff.sh`: combines the per-protein results into one table.

Overrides: `SPEC` (covariate set), `GRM_CUTOFF`, `MIN_PROTEIN_N`, `PROTEIN_SET`
(a file of protein names), `PERM_OUTPUT_ROOT`. The jobs call the R scripts in
`scripts/` in this order: `export_architecture_inputs.R`, `filter_prep_inputs_by_keep.R`,
`build_genotype_grm.R`, `build_environment_kernels.R`, `run_population_architecture.R`,
`summarize_population_architecture.R`. Settings are in `config/default_config.R`.

## Inputs

The analysis dataset (`HEAP.rds`), UK Biobank genotypes in PLINK 2 format
(`$HEAP_IGLOO_ROOT/UKB/gwas/UKBallchr.{pgen,pvar,psam}`) with an LD-pruned SNP list
(`ukb_nonimputed_snps.pvar`), and the GCTA and PLINK 2 binaries
([`docs/ENVIRONMENT.md`](../../docs/ENVIRONMENT.md)). The
[demo](../../demo/README.md) runs these three jobs on synthetic genotypes.
