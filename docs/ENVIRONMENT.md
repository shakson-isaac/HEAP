# Environment

## System

Developed and run on Linux (Harvard Medical School O2 cluster, x86-64, Red Hat Enterprise
Linux 9.8); not tested on macOS or Windows. No GPU or other special hardware. The full
analysis runs on a SLURM cluster: per-task needs range from 1 core / 4 GB (genetic scores)
to 1 core / 60 GB (GREML) and 8 cores / 50 GB (regenie), as declared in `slurm/`. The
[demo](../demo/README.md) needs one core and under 2 GB; its GREML step also needs GCTA and PLINK 2.

## R

R 4.4.2 (built with GCC 14.2.0) and 46 packages, listed with the versions used in
[`config/r_packages.tsv`](../config/r_packages.tsv). From the repository root:

```bash
Rscript scripts/setup/install_r_packages.R
```

This installs the missing packages into your default R library, or into `HEAP_RLIB` if
set: CRAN packages with `install.packages`, Bioconductor packages with `BiocManager`, and
three GitHub packages with `remotes` (`MRCIEU/TwoSampleMR`, `MRCIEU/genetics.binaRies`,
and `rmgpanw/ukbwranglr` pinned to commit `50d7ae4`). Current releases are installed,
except two pinned versions: ggplot2 3.5.2 (4.x breaks the `ggtree → enrichplot →
clusterProfiler` stack used for enrichment) and qs 0.27.3 (archived on CRAN, installed
from the archive). On Linux, packages are built from source and need C, C++ and Fortran
compilers. Rerunning the command installs only what is still missing.

## Genetics tools

| Tool | Version | Used for |
|---|---|---|
| [regenie](https://rgcgithub.github.io/regenie/) | 4.1 | Exposure GWAS (`conda env create -f slurm/gwas_regenie/environment.yml`) |
| [GCTA](https://yanglab.westlake.edu.cn/software/gcta/) | 1.94.1 | GREML |
| [PLINK 2](https://www.cog-genomics.org/plink/2.0/) | 2.0.20220814 | Protein genetic scores |
| [bgenix](https://enkre.net/cgi-bin/code/bgen) | 1.1.7 | Extracting variants from UK Biobank `.bgen` files |
| [LDSC](https://github.com/bulik/ldsc) | `master` (no commit pinned) | LD score regression (Python 2.7; `slurm/ldsc/environment.yml`, or `environment_pinned.yml` if that fails to resolve) |

## Environment variables

The analysis scripts resolve their paths through these variables. When one is unset,
they fall back to the original cluster's paths (`/n/groups/patel/...`), which will not exist
elsewhere. Some helper scripts, mainly the colocalization scripts in
`scripts/support/coloc/` and a few figure and summary scripts, still set those paths
directly; edit the path lines near the top of each to run them elsewhere. `IGLOO` is the
name of the lab's shared data directory on the original cluster.

| Variable | Points at |
|---|---|
| `HEAP_ROOT` | This repository |
| `HEAP_IGLOO_ROOT` | Root for outputs and reference data |
| `HEAP_SCRATCH_ROOT` | Scratch space for job staging |
| `HEAP_LOADER_RDS` | The analysis dataset (`HEAP.rds`) |
| `HEAP_RLIB` | R library to use (optional) |
| `HEAP_REGENIE_ENV` | The regenie conda environment |
| `GCTA_BIN` | The `gcta64` binary |
| `LDSC_HOME`, `LDSC_LD_DIR` | LDSC checkout; LD score directory (`eur_w_ld_chr/`) |
