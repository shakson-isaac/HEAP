# HEAP: Human Exposomic Architecture of the Proteome

Analysis code for **"Human Plasma Proteomics Links Modifiable Lifestyle Exposome to Disease
Risk"** (Isaac *et al.*, in review).

[**Interactive results — heap.bio**](https://heap.bio) ·
[Reproduce the analysis](docs/REPRODUCIBILITY.md) ·
[Reproduce the figures](#reproducing-the-figures) ·
[Data access](docs/DATA_ACCESS.md) ·
[Software environment](docs/ENVIRONMENT.md)

![HEAP overview](HEAP.png)

HEAP links lifestyle exposures, plasma proteins and disease risk in more than 50,000
UK Biobank participants (2,686 Olink proteins, 169 exposomic features). It asks how much
of each protein's variation is explained by genetics, the exposome, and their
interaction; which exposure–protein associations replicate; which proteins lie on the
path from an exposure to disease (mediation, Mendelian randomization, colocalization);
whether observational signatures match the proteomic response to exercise training and
GLP-1 receptor agonists; and whether proteome-based exposure scores track behavior and
predict disease.

---

## Contents

- [What is in this repository](#what-is-in-this-repository)
- [Quick start](#quick-start)
- [Reproducing the analysis](#reproducing-the-analysis)
- [Reproducing the figures](#reproducing-the-figures)
- [Finding more detail](#finding-more-detail)
- [Citation](#citation) · [License](#license) · [Contact](#contact)

## What is in this repository

```
scripts/     Analysis code, one directory per analysis (table below)
  visualizations/   Code that draws every main and supplementary figure
  support/, analysis_summaries/   Summary tables computed from module outputs for the figures
config/      Every analysis setting: covariate sets, sample filters, exposure and
             protein sets, per-module experiment definitions. Scripts hardcode none.
slurm/       Job scripts that submit each analysis to a SLURM cluster
workflow/    Path configuration, config validation, manifest generation
docs/        Run guide, data access, software environment
```

The paper is organized into six analyses. The code directories are numbered in the
order they were written, which differs from the paper in two places (rows marked †):

| Paper | Analysis | Main figure | Code |
|---|---|---|---|
| Module 1 | Genetic vs. exposomic variance of each protein | Fig. 1 | [`scripts/module1_variance_decomposition`](scripts/module1_variance_decomposition) |
| Module 2 | Exposure–protein associations | Fig. 2 | [`scripts/module2_associations`](scripts/module2_associations) |
| Module 3 | Mediation of exposure → disease through proteins | Fig. 3 | [`scripts/module3_mediation`](scripts/module3_mediation) |
| Module 4 † | Mendelian randomization and colocalization | Fig. 4 | [`scripts/module5_mr`](scripts/module5_mr) |
| Module 5 | Comparison with intervention trials | Fig. 5 | [`scripts/support/intervention_compare`](scripts/support/intervention_compare) |
| Module 6 | Proteome-based exposure scores (PES) | Fig. 6 | [`scripts/module6_prediction`](scripts/module6_prediction) |
| — † | Tissue and pathway enrichment (supports Figs. 2–3) | | [`scripts/module4_enrichment`](scripts/module4_enrichment) |

Supporting analyses: [`loaders`](scripts/loaders) (builds the analysis dataset from
UK Biobank), [`genetic_scores`](scripts/genetic_scores) (cis/trans protein genetic
scores), [`gwas_regenie`](scripts/gwas_regenie) and [`ldsc`](scripts/ldsc) (exposure
GWAS and LD score regression, which supply the MR instruments), and
[`population_architecture`](scripts/population_architecture) (GREML variance
components).

## Quick start

**1. Get the code**

```bash
git clone https://github.com/shakson-isaac/HEAP.git
cd HEAP
```

**2. Install the software.** R 4.4.2 and ~46 packages:

```bash
Rscript scripts/setup/install_r_packages.R --only-missing
```

The genetics steps also need regenie, GCTA, PLINK 2 and LDSC; conda specifications are
checked in. Versions and sources: [`docs/ENVIRONMENT.md`](docs/ENVIRONMENT.md).

**3. Get the data.** UK Biobank individual-level data need an approved
[UK Biobank application](https://www.ukbiobank.ac.uk/enable-your-research/apply-for-access)
and cannot be redistributed. The external summary statistics (OmicsPred, pQTL, FinnGen,
HERITAGE, semaglutide trials) are public. Sources: [`docs/DATA_ACCESS.md`](docs/DATA_ACCESS.md).

*No UK Biobank access?* The summary-level results behind every figure are in the
paper's Supplementary Tables and Data, downloadable from [heap.bio](https://heap.bio).

**4. Point HEAP at your directories.** No path is hardcoded:

```bash
export HEAP_ROOT=/path/to/this/repo
export HEAP_SCRATCH_ROOT=/path/to/scratch
export HEAP_IGLOO_ROOT=/path/to/output/root
```

## Reproducing the analysis

The full run guide, with commands, dependencies and runtimes for each stage, is
[`docs/REPRODUCIBILITY.md`](docs/REPRODUCIBILITY.md). In outline:

```
UK Biobank ─→ loader (HEAP.rds) ─┬─→ Module 1 ─→ Module 3
                                 ├─→ Module 2 ─→ enrichment
                                 ├─→ Module 6
  genetic scores, exposure GWAS ─┴─→ Module 5 MR (+ pQTL, FinnGen)
                                     GREML, LD score regression
```

1. Build the analysis dataset: `sbatch slurm/loaders/run_HEAP_loader.sh`.
2. Run the genetics steps (protein genetic scores, exposure GWAS, LDSC, GREML).
3. Run each module. Experiments are named in `config/modules/<module>_experiments.yml`,
   for example `EXPERIMENT=M1_base_lasso bash slurm/module1/HEAPmodule1_manifest.sh`.

The primary analysis uses the `base` covariate set. The other sets in
`config/covariates/covariate_sets.yml` are sensitivity analyses. An end-to-end run
takes about 1–2 weeks on a cluster with 100–300 concurrent job slots, most of it in the
genetic scores, exposure GWAS and GREML.

## Reproducing the figures

Every main and supplementary figure is drawn by code in
[`scripts/visualizations`](scripts/visualizations). After the pipeline has run:

```bash
bash scripts/visualizations/make_figures.sh          # Figs 1-6 + supplementary
bash scripts/visualizations/make_figures.sh fig4     # a single main figure
```

The finished main figures are written to `$HEAP_PROJECT_ROOT/figures/manuscript/`.
[`scripts/visualizations/README.md`](scripts/visualizations/README.md) maps each figure
and panel to its script and lists the summary scripts that build the figures' inputs.

**Without UK Biobank access.** The figure scripts need pipeline outputs, which require
UK Biobank data. The summary-level results behind every figure are released with the
paper as Supplementary Tables and Supplementary Data, downloadable from
[heap.bio](https://heap.bio).

## Finding more detail

| I want to know… | Look in |
|---|---|
| Why each method was chosen and how it was specified | The paper's Methods and Supplementary Notes |
| The exact command and order to run each stage | [`docs/REPRODUCIBILITY.md`](docs/REPRODUCIBILITY.md) |
| What each stage reads and writes | [`config/io_map.yml`](config/io_map.yml) |
| The covariates, sample filters or exposure set used | [`config/`](config) |
| Software versions and reference datasets | [`docs/ENVIRONMENT.md`](docs/ENVIRONMENT.md) |
| Details of a specific analysis | the README in that analysis's directory, e.g. [`scripts/module5_mr`](scripts/module5_mr/README.md), [`slurm/gwas_regenie`](slurm/gwas_regenie/README_gwas_exposure_workflow.md) |
| The results themselves | Supplementary Tables and Data accompanying the paper, and [heap.bio](https://heap.bio) |

## Citation

Isaac S, Ellis RJ, Jee YH, Murthy VL, Udler MS, Neale BM, Sunyaev S, Gusev A,
Martin AR, Patel CJ. *Human Plasma Proteomics Links Modifiable Lifestyle Exposome to
Disease Risk.* In review.

<!-- TODO: add DOI and BibTeX once the preprint/article is public -->

## License

Code is released under the [MIT License](LICENSE). UK Biobank data are subject to the
UK Biobank access conditions and are not covered by this license.

## Contact

Shakson Isaac — shakson_isaac@g.harvard.edu. Bug reports and questions:
[GitHub issues](https://github.com/shakson-isaac/HEAP/issues).
