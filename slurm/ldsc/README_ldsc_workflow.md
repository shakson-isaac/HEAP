# LD score regression

Runs LDSC (Bulik-Sullivan *et al.* 2015) on the exposure GWAS
([`slurm/gwas_regenie`](../gwas_regenie/README_gwas_exposure_workflow.md)) to estimate, per
exposure:

- **SNP heritability** (`h2`, observed scale);
- **the LDSC intercept**, a test of genomic inflation: close to 1 when inflation reflects
  polygenic signal, above 1 when it reflects confounding;
- **genetic correlations** (`rg`) between well-powered exposures (h2 z ≥ 4 and h2 ≥ 0.02,
  listed in `ldsc_rg_exposures.txt`).

Results appear in the GWAS supplementary figures (`fig_ldsc_h2`, `fig_ldsc_intercept`,
`fig_ldsc_rg`).

## Run

```bash
conda env create -f slurm/ldsc/environment.yml   # Python 2.7; or environment_pinned.yml
bash slurm/ldsc/submit_ldsc.sh                   # h2 for every finished GWAS; --rg adds rg
Rscript scripts/ldsc/collect_ldsc_h2.R           # after the jobs finish
Rscript scripts/ldsc/collect_ldsc_rg.R           # if rg was run
```

`submit_ldsc.sh` lists the finished exposure GWAS (`ldsc_exposures.txt`) and submits
`ldsc_h2_array.sh`, one task per exposure, each running `scripts/ldsc/run_ldsc_h2.sh`:
convert regenie output (LOG10P → P, effect allele `ALLELE1`), munge to HapMap3 SNPs, and run
`ldsc.py --h2`. `ldsc_rg.sh` computes the pairwise correlations. Exposures that already have
a result are skipped (`LDSC_FORCE=1` reruns them). One task needs 1 core, 16 GB, under 2 hours.

## Inputs and outputs

Needs the LDSC code and the EUR LD scores (`eur_w_ld_chr/`, with `w_hm3.snplist`), set with
`LDSC_HOME`, `LDSC_LD_DIR`, `LDSC_SNPLIST` and `LDSC_ENV`
([`docs/ENVIRONMENT.md`](../../docs/ENVIRONMENT.md)). Outputs, under
`$HEAP_IGLOO_ROOT/UKB/HEAP/output/gwas/ldsc/`: `munged/`, `h2/<exposure>.log`, `rg/`, and the
tables `ldsc_h2_summary.tsv` and `ldsc_rg_summary.tsv`.

The LD reference is European, matching the GWAS sample. For binary exposures `h2` is on the
observed scale; no liability-scale conversion is applied.
