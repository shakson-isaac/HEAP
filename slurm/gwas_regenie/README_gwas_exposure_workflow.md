# Exposure GWAS (regenie)

Genome-wide association studies of the lifestyle exposures, run with regenie. They supply
the exposure instruments for Mendelian randomization
([`scripts/module5_mr`](../../scripts/module5_mr/README.md)).

**Sample.** Everyone in the proteomics cohort is removed before any phenotype file is
written, so the exposure GWAS and the pQTLs come from non-overlapping UK Biobank
participants, as two-sample MR requires.

**Exposures.** Those in `config/exposure_sets/analysis_exposures.tsv` with `include == 1`
and under 20% missing in the proteomics cohort at baseline.

| Exposure type | regenie model |
|---|---|
| Continuous, ordinal | Quantitative, rank-inverse normal transformed (`--apply-rint`) |
| Binary | `--bt`, with Firth correction in step 2 (`--firth --approx --pThresh 0.01`) |

**Covariates.** Age, sex, age², age × sex, age² × sex and genetic principal components 1–20.

**Genotypes.** Step 1 uses the genotyped SNPs (`ukb_nonimputed_snps`), step 2 the full panel
(`UKBallchr`), both in PLINK 2 format. Their directory is set as `GWAS_PATH` in the two
array scripts.

## Run

```bash
conda env create -f slurm/gwas_regenie/environment.yml   # regenie 4.1, once
bash slurm/gwas_regenie/submit_gwas_exposures.sh          # --dry-run to preview
```

The launcher runs `scripts/gwas_regenie/prepare_gwas_exposures.R`, which writes the
exposure lists (`evars_continuous_heap.txt`, `evars_binary_heap.txt`) and QC tables, then
submits one array task per exposure: `gwas_regenie_exposures_continuous_v2.sh` and
`gwas_regenie_exposures_binary_v2.sh`. Each task builds its own phenotype and covariate
files from `HEAP.rds` and runs regenie steps 1 and 2.

## Outputs

| File | Contents |
|---|---|
| `output/gwas/regenie_step2/<exposure>/<exposure>.regenie` | Summary statistics (MR input) |
| `output/gwas_regenie/sample_counts.tsv` | Participants at each filtering step |
| `output/gwas_regenie/exposure_gwas_qc.tsv` | Per exposure: type, sample size, cases and controls, whether `--apply-rint` was used |

Paths are under `$HEAP_IGLOO_ROOT/UKB/HEAP/`.
`scripts/gwas_regenie/stage_gwas_deposit.R` converts the summary statistics into
bgzipped, indexed files for public deposit.
