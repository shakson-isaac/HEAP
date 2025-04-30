# Environment

All analyses were performed using **R version 4.2.1**. The following key R packages were used:

- `glmnet` — for penalized regression and model fitting  
- `survival` (specifically `coxph`) — for Cox proportional hazards models  
- `clusterProfiler` — for functional enrichment analysis

Additional R packages used for data manipulation, visualization, and parallelization are listed in the corresponding scripts and can be installed via `install.packages()` or `BiocManager::install()` as needed.

## HPC Environment

We utilized a compute cluster running the **SLURM** job scheduling system to parallelize proteome-wide analyses. Scripts such as `runProt_*.R` were executed across multiple nodes and cores using custom `bash` scripts and job arrays.

Typical compute requirements per job:
- **Memory**: 10–20 GB RAM  
- **Cores**: 1–8 per job  
- **Runtime**: up to 6 hours depending on module

## Modules and External Tools

The following software modules were loaded as part of the environment (using `module load`):

- `gcc/9.2.0` — compiler for R and system dependencies  
- `cmake/3.22.2` — for C/C++-based package compilation  
- `fftw/3.3.10` — required for packages using fast Fourier transforms  
- `plink2/2.0.20220814` — for genomic data processing  
- `sqlite3/3.43.1` — for lightweight SQL database access

## Contact

For further information or questions, please contact:  shakson_isaac@g.harvard.edu
