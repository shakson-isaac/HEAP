# HEAP Exposure GWAS Workflow (REGENIE)

## Overview

This workflow runs genome-wide association studies (GWAS) for lifestyle
exposures using REGENIE, supporting downstream two-sample Mendelian
Randomization (MR) analyses.

The workflow is split into two phases:
1. **Preparation** (R script, batch mode) – generates the per-type exposure
   lists and QC tables from HEAP.rds (no per-exposure pheno/covar staging).
2. **REGENIE jobs** (SLURM array) – one REGENIE job per exposure; each job
   regenerates its own pheno/covar input (single-exposure mode) from HEAP.rds.

---

## Two-sample MR design

Exposure GWAS summary statistics must come from a sample that does **not**
overlap with the proteomics cohort used in HEAP protein association analyses.
This is required so that the two-sample MR instruments (exposure → protein)
are derived from independent samples.

`prepare_gwas_exposures.R` enforces this by identifying proteomics cohort
participants from `heap$prot_baseline$eid` and removing them before writing
any phenotype/covariate files.

---

## Exposure eligibility

Eligibility is defined by `HEAP/config/exposure_sets/analysis_exposures.tsv`:

| Column | Meaning |
|--------|---------|
| `variable` | Exposure column name in HEAP.rds |
| `variable_type` | `continuous`, `ordinal`, or `binary` |
| `miss_rate_prot_i0` | Missingness rate in the proteomics cohort at instance 0 |
| `include` | 1 = candidate for inclusion, 0 = excluded |

An exposure is eligible if **`include == 1` AND `miss_rate_prot_i0 < 0.20`**.

The 20% missingness threshold is consistent with the HEAP-wide 80%
complete-case design for exposomic measurements.

---

## Quantitative vs binary trait normalization

| `variable_type` | REGENIE model | `--apply-rint` |
|-----------------|--------------|----------------|
| `continuous` | quantitative | YES |
| `ordinal` | quantitative (treated as continuous scale) | YES |
| `binary` | `--bt` + Firth logistic | NO |

`--apply-rint` applies a rank-inverse normal transformation to quantitative
phenotypes inside REGENIE (step 2 only). Binary traits use `--bt` in both
steps, with Firth correction (`--firth --approx --pThresh 0.01`) in step 2.

---

## Quick start (recommended)

One launcher does prep + sizing + submission. Run it from an O2 login node:

```bash
/n/groups/patel/shakson_ukb/HEAP/slurm/gwas_regenie/submit_gwas_exposures.sh
```

It runs `prepare_gwas_exposures.R` (via `srun`, since prep loads HEAP.rds and is
too heavy for the login node), reads the line counts of the two generated
`evars_*_heap.txt` lists to size each array automatically, then submits the
continuous and binary REGENIE arrays. No one computes `N_CONT`/`N_BIN` by hand.

```bash
submit_gwas_exposures.sh --dry-run     # print every srun/sbatch, run nothing
submit_gwas_exposures.sh --skip-prep   # lists already current; just submit
submit_gwas_exposures.sh --help        # full option list
```

Before the first production run, make sure HEAP.rds is current
(`Rscript scripts/loaders/HEAP_loader.R`), then inspect the QC tables it writes
(`output/gwas_regenie/sample_counts.tsv`, `exposure_gwas_qc.tsv`) — verify
`apply_rint` is TRUE only for quantitative traits.

---

## Manual step-by-step (what the launcher automates)

### 1. Ensure HEAP.rds is up to date

```bash
cd /n/groups/patel/shakson_ukb/HEAP
Rscript scripts/loaders/HEAP_loader.R
```

### 2. Run the exposure preparation script

```bash
cd /n/groups/patel/shakson_ukb/HEAP
Rscript scripts/gwas_regenie/prepare_gwas_exposures.R
```

This script (in **batch mode**) will:
- Read `analysis_exposures.tsv`, filter to eligible exposures
- Load HEAP.rds, stack all baseline exposure data
- Exclude proteomics participants
- Write `slurm/gwas_regenie/evars_continuous_heap.txt`
- Write `slurm/gwas_regenie/evars_binary_heap.txt`
- Write `output/gwas_regenie/exposure_gwas_qc.tsv`
- Write `output/gwas_regenie/sample_counts.tsv`

Batch mode does **not** write per-exposure `pheno.txt`/`covar.txt` files: each
REGENIE array job regenerates its own (single-exposure mode, into job-local tmp),
so staging all of them here would be wasted I/O.

### 3. Check the QC table

```bash
cat /n/groups/patel/shakson_ukb/HEAP/output/gwas_regenie/sample_counts.tsv
head /n/groups/patel/shakson_ukb/HEAP/output/gwas_regenie/exposure_gwas_qc.tsv
```

Verify sample counts and that `apply_rint` is TRUE only for quantitative traits.

### 4. Submit SLURM array jobs

Get the array sizes from the evar list files:
```bash
N_CONT=$(wc -l < slurm/gwas_regenie/evars_continuous_heap.txt)
N_BIN=$(wc -l  < slurm/gwas_regenie/evars_binary_heap.txt)
echo "Continuous: ${N_CONT}, Binary: ${N_BIN}"
```

Submit:
```bash
cd /n/groups/patel/shakson_ukb/HEAP

# Quantitative exposures (continuous + ordinal)
sbatch --array=1-${N_CONT} slurm/gwas_regenie/gwas_regenie_exposures_continuous_v2.sh

# Binary exposures
sbatch --array=1-${N_BIN} slurm/gwas_regenie/gwas_regenie_exposures_binary_v2.sh
```

---

## File and directory layout

```
HEAP/
├── config/exposure_sets/
│   └── analysis_exposures.tsv          ← exposure eligibility list (source of truth)
├── scripts/gwas_regenie/
│   └── prepare_gwas_exposures.R        ← phenotype/covariate prep script
├── slurm/gwas_regenie/
│   ├── submit_gwas_exposures.sh        ← one-command launcher (prep → size → submit)
│   ├── evars_continuous_heap.txt       ← generated by prep script (regenerated each run)
│   ├── evars_binary_heap.txt           ← generated by prep script (regenerated each run)
│   ├── gwas_regenie_exposures_continuous_v2.sh   ← SLURM array (quantitative)
│   ├── gwas_regenie_exposures_binary_v2.sh       ← SLURM array (binary)
│   ├── environment.yml                 ← regenie conda env spec
│   └── README_gwas_exposure_workflow.md
└── output/gwas_regenie/
    ├── exposure_gwas_qc.tsv            ← per-exposure QC metrics (prep batch mode)
    └── sample_counts.tsv               ← global sample flow counts (prep batch mode)

/n/groups/patel/IGLOO/UKB/HEAP/                       ← canonical outputs (group-shared)
├── output/gwas/regenie_step2/<exposure_name>/
│   ├── <exposure_name>.regenie         ← final step-2 summary stats (MR input)
│   └── <exposure_name>.log             ← REGENIE step-2 log
└── logs/gwas/
    ├── gwas_exposure_cont_<user>_<JID>_<AID>.out
    └── gwas_exposure_bin_<user>_<JID>_<AID>.out

$SLURM_TMPDIR/heap_gwas_<jobid>_<taskid>/             ← per-job temp, node-local, auto-cleaned
├── exposure_input/<exposure_name>/
│   ├── pheno.txt                       ← FID IID <trait>   (generated per job)
│   └── covar.txt                       ← FID IID age sex age2 ... PC20
├── regenie_step1/                      ← step-1 null model (ephemeral)
└── regenie_step1_lowmem/               ← step-1 lowmem scratch

${HEAP_SCRATCH_ROOT}/regenie_step2/<exposure_name>/   ← step-2 working dir (fresh each run)
```

---

## Covariates used

| Covariate | Description |
|-----------|-------------|
| `age_when_attended_assessment_centre_f21003_0_0` | Age at baseline assessment |
| `sex_f31_0_0` | Sex (recoded: Male=1, Female=0) |
| `age2` | Age squared |
| `age_sex` | Age × sex interaction |
| `age2_sex` | Age² × sex interaction |
| `genetic_principal_components_f22009_0_1..20` | Genetic PCs 1–20 |

---

## QC table columns (exposure_gwas_qc.tsv)

| Column | Description |
|--------|-------------|
| `exposure_id` | Variable name |
| `category` | Exposure category (e.g., Smoking, Diet_Weekly) |
| `exposure_type` | continuous / ordinal / binary |
| `miss_rate_from_tsv` | Missingness rate from analysis_exposures.tsv |
| `n_gwas_base` | Non-proteomics participants with complete covariates |
| `n_with_phenotype` | Non-missing phenotype values in GWAS base |
| `n_missing_phenotype` | Missing phenotype values in GWAS base |
| `n_proteomics_excluded` | Proteomics cohort participants removed |
| `n_final_gwas` | Final GWAS sample size for this exposure |
| `apply_rint` | TRUE if `--apply-rint` is used in REGENIE |
| `n_cases_binary` | Number of cases (binary traits only) |
| `n_controls_binary` | Number of controls (binary traits only) |

---

## Versioning

The `_v2` scripts are the current, canonical workflow. They replaced the
original `_v1` scripts (without version suffix), which used the old
`UK_Biobank/BScripts/GWAS/REGENIE/DATA/` input paths and have since been
removed. The `_v2` scripts generate per-exposure pheno/covar files inside a
job-local temp dir and copy step-2 summary stats to the canonical IGLOO
location (`output/gwas/regenie_step2/<exposure>/<exposure>.regenie`).
