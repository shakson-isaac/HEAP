# Environment

What you need to install to run HEAP, and how the scripts find it. Nothing below is
hardcoded — every tool and reference dataset is located through an environment
variable, so HEAP does not assume the cluster it was developed on.

## R

**R 4.4.2**, compiled under **GCC 14.2.0**. The analysis scripts need ~46 packages,
tracked in [`config/r_packages.tsv`](../config/r_packages.tsv) (36 CRAN, 8 Bioconductor,
2 from GitHub). Install them all with:

```bash
Rscript scripts/setup/install_r_packages.R --only-missing
```

`workflow/00_paths.R` prepends `$HEAP_RLIB` to `.libPaths()` when set, so a shared or
project-local library can be used without touching your personal one.

> **Version pin — keep ggplot2 on 3.5.x.** ggplot2 4.0.0 removed the internal
> `check_linewidth()` that `ggtree` still calls, which breaks the
> `ggtree → enrichplot → clusterProfiler` stack used by the enrichment module with a
> lazy-load error. HEAP was run against ggplot2 3.5.2.

Two packages come from GitHub rather than CRAN:
[`MRCIEU/TwoSampleMR`](https://github.com/MRCIEU/TwoSampleMR) and
[`MRCIEU/genetics.binaRies`](https://github.com/MRCIEU/genetics.binaRies) (the latter
ships the PLINK binary that TwoSampleMR's clumping calls).

## Genetics toolchain

| Tool | Version used | Where to get it | Used by |
|---|---|---|---|
| **regenie** | 4.1 | [rgcgithub.github.io/regenie](https://rgcgithub.github.io/regenie/) — or the conda spec below | Exposure GWAS |
| **GCTA** | 1.94.1 | [yanglab.westlake.edu.cn/software/gcta](https://yanglab.westlake.edu.cn/software/gcta/) | Population architecture (GREML) |
| **PLINK 2** | 2.0.20220814 | [cog-genomics.org/plink/2.0](https://www.cog-genomics.org/plink/2.0/) | Protein genetic scores |
| **bgenix** | — | [BGEN library](https://enkre.net/cgi-bin/code/bgen) | Extracting variants from UKB `.bgen` |
| **LDSC** | upstream `master` (no commit pinned) | [github.com/bulik/ldsc](https://github.com/bulik/ldsc) — plus the conda spec below | SNP-heritability, genetic correlation |

### Conda environments (checked in)

Two environments are specified in the repository, so you do not have to reconstruct them:

```bash
# regenie (exposure GWAS)
conda env create -f slurm/gwas_regenie/environment.yml -p /path/to/regenie_env
export HEAP_REGENIE_ENV=/path/to/regenie_env

# LDSC (Python 2.7 — LDSC has not been ported to Python 3)
conda env create -f slurm/ldsc/environment.yml -p /path/to/ldsc_env
```

`slurm/ldsc/environment_pinned.yml` holds a fully solved version of the LDSC
environment if the loose spec fails to resolve.

## Reference data

These are third-party datasets HEAP reads but does not redistribute. Point the
corresponding variable at your own copy.

| Dataset | Purpose | Source |
|---|---|---|
| UK Biobank phenotypes, proteomics, genotypes | Everything | [UK Biobank application](https://www.ukbiobank.ac.uk/enable-your-research/apply-for-access) — see [`DATA_ACCESS.md`](DATA_ACCESS.md) |
| OmicsPred protein score weights | cis/trans protein genetic scores | [omicspred.org](https://www.omicspred.org/) |
| LDSC EUR LD scores (`eur_w_ld_chr/`) and the HapMap3 `w_hm3.snplist` | LD score regression | Distributed with [LDSC](https://github.com/bulik/ldsc) — see its documentation |
| 1000 Genomes EUR LD reference (`.bed/.bim/.fam`) | MR clumping | Any EUR PLINK reference panel |
| UK Biobank pQTL summary statistics | Mendelian randomization | UK Biobank Pharma Proteomics Project |
| deCODE plasma pQTL | MR replication arm | [decode.com/summarydata](https://www.decode.com/summarydata/) |
| FinnGen disease GWAS | MR outcomes | [finngen.fi](https://www.finngen.fi/en/access_results) |

If your GWAS cohort is not European, swap `LDSC_LD_DIR` and the MR LD reference for an
ancestry-matched panel.

## Environment variables

Set these to point HEAP at your own installations and data. None have a portable
default — the values baked into the scripts refer to the original cluster.

| Variable | Points at |
|---|---|
| `HEAP_ROOT` | This repository |
| `HEAP_SCRATCH_ROOT` | Fast scratch space for job staging |
| `HEAP_IGLOO_ROOT` | Root for outputs and shared reference data |
| `HEAP_RLIB` | R library to prepend to `.libPaths()` (optional) |
| `HEAP_REGENIE_ENV` | The regenie conda environment |
| `GCTA_BIN` | The `gcta64` binary |
| `LDSC_HOME` | Your LDSC checkout (the directory holding `ldsc.py`) |
| `LDSC_LD_DIR` | LD score directory (`eur_w_ld_chr/`) |
| `HEAP_LOADER_RDS` | The `HEAP.rds` object built by the loader |

## Compute

HEAP was run on a SLURM cluster. Per-task requirements vary widely by stage — from
1 core / 4 GB for genetic scores to 1 core / 60 GB for the GREML array and
8 cores / 50 GB for regenie. Per-stage resources are declared in the `slurm/` scripts,
and expected runtimes are tabulated in [`REPRODUCIBILITY.md`](REPRODUCIBILITY.md).

Nothing requires SLURM in principle — the R entry points can be driven by any scheduler,
or run serially, given enough time.

## Contact

Questions: shakson_isaac@g.harvard.edu
