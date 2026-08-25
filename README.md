# HEAP

**H**uman **E**xposomic **A**rchitecture of the **P**roteome — analysis code for the HEAP
manuscript.

![HEAP overview](HEAP.png)

HEAP partitions variation in 2,686 plasma proteins across 53,014 UK Biobank participants
into genetic, exposomic, and gene-by-environment components against 169 exposomic features,
then traces which proteins act as causal intermediates between exposure and disease and which
are downstream reporters of exposure.

Interactive results: **[heap.bio](https://heap.bio)**

---

## Repository layout

```
workflow/      Path configuration, config validation, and manifest generation
config/        YAML/TSV configuration - covariate sets, sample filters, exposure sets,
               protein sets, and per-module experiment definitions. Nothing that
               governs an analysis is hardcoded in a script.
scripts/       Analysis code, one directory per module (see below)
slurm/         Job scripts for the genetics toolchain (GWAS, LDSC, GREML,
               protein genetic scores) -- see Release scope below
docs/          How to reproduce, how to get the data, and the software environment
```

### Analysis modules

| Directory | What it does |
|---|---|
| `scripts/loaders` | Assembles the analysis matrices from UK Biobank source data |
| `scripts/genetic_scores` | Splits OmicsPred protein scores into cis / trans / total |
| `scripts/module1_variance_decomposition` | Partitions protein variance into genetic, exposomic, GxE, and covariate components |
| `scripts/module2_associations` | Univariate exposure, genetic, and GxE association models |
| `scripts/module3_mediation` | Mediation of exposure→disease effects through the proteome |
| `scripts/module4_enrichment` | Tissue and pathway GSEA of the association results |
| `scripts/module5_mr` | Bidirectional Mendelian randomization and colocalization |
| `scripts/module6_prediction` | Proteomic exposure scores and longitudinal disease prediction |
| `scripts/population_architecture` | GREML variance-component estimation |
| `scripts/gwas_regenie` | Genome-wide association analysis of the exposures |
| `scripts/ldsc` | LD score regression on the exposure GWAS |
| `scripts/setup` | Environment and dependency setup |

### Manuscript vs. code module numbers

These deliberately differ. Use manuscript numbers when reading the paper, code numbers
when navigating this repository.

| Manuscript | Concept | Code |
|---|---|---|
| Module 1 | Variance / exposure-responsive spectrum | `scripts/module1_variance_decomposition` |
| Module 2 | Exposure–protein association | `scripts/module2_associations` |
| Module 3 | Mediation | `scripts/module3_mediation` |
| Module 4 | **Mendelian randomization** | `scripts/module5_mr` |
| Module 5 | Interventional comparison | figure code — not in this release |
| Module 6 | Proteomic exposure scores | `scripts/module6_prediction` |
| *(un-numbered)* | Tissue / pathway enrichment | `scripts/module4_enrichment` |

Note that the code prefix `module4` refers to **enrichment**, an un-numbered supporting
analysis — not to the manuscript's Module 4, which is Mendelian randomization and lives
under `scripts/module5_mr`. In the manuscript, MR precedes the interventional comparison
because the intervention analysis annotates proteins by their MR causal edge.

## Getting started

All paths are configured through environment variables — no analysis path is hardcoded:

```bash
export HEAP_ROOT=/path/to/this/repo
export HEAP_SCRATCH_ROOT=/path/to/scratch
export HEAP_IGLOO_ROOT=/path/to/output/root
```

Then source `workflow/00_paths.R` from R, which also sets the shared library path and
output permissions. The full stage-by-stage run guide is
[`docs/REPRODUCIBILITY.md`](docs/REPRODUCIBILITY.md), with the machine-readable
dependency map in `config/io_map.yml`.

Analyses were run on a SLURM cluster (Harvard O2) with R 4.4.2 under GCC 14.2.0. See
[`docs/ENVIRONMENT.md`](docs/ENVIRONMENT.md) for the full environment specification.

## Data availability

This repository contains **code only**. UK Biobank individual-level data are
access-restricted and cannot be redistributed; they must be obtained directly from
UK Biobank under an approved application. See
[`docs/DATA_ACCESS.md`](docs/DATA_ACCESS.md).

To reproduce the analyses with your own UK Biobank instance, follow
[`docs/REPRODUCIBILITY.md`](docs/REPRODUCIBILITY.md).

## Release scope

This is a staged release. It contains the analysis code, the configuration that
drives it, and the job scripts for the genetics toolchain — exposure GWAS
(regenie), LD score regression, GREML variance components, and protein genetic
scores.

Not included in this release:

- **Figure generation** — plotting code, the figure registry, and figure legends.
  Visualizations are published separately alongside the
  [HEAP website](https://heap.bio).
- **Job wrappers for Modules 1–6** — thin `sbatch` array wrappers. The analysis
  they submit is in `scripts/`, which is included here.
- **Derived summary statistics** — these accompany the manuscript as
  supplementary data.

`config/io_map.yml` maps the complete pipeline, so it refers to some stages not
included above.

## Citation

<!-- TODO: fill in once the manuscript has a DOI / preprint URL -->

## License

See [LICENSE](LICENSE).
