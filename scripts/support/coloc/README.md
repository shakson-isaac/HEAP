# Colocalization

Tests whether a protein's cis-pQTL signal and an outcome GWAS signal share one causal
variant, for the cis-pQTL edges that pass Mendelian randomization. A cis edge counts
as colocalized at PP.H4 ≥ 0.8. A cis Tier 1 edge that was tested and did not
colocalize is demoted to Tier 2.

Every pass uses `coloc::coloc.abf` with its default priors (p1 = p2 = 1e-4,
p12 = 1e-5) over a ±500 kb window around the protein's lead cis instrument. coloc.abf
assumes a single causal variant and needs no LD reference.

Outcomes are FinnGen R12 diseases for protein → disease edges and the HEAP exposure
GWAS for protein → exposure edges. Both pQTL arms are covered: UK Biobank (UKB-PPP)
and deCODE.

## Scripts, in run order

Run these after the MR tables exist (`scripts/support/mr_tables/`). The two
colocalization passes feed the MR tier table.

| Script | Does | Writes (under `output/support/coloc/`) |
|---|---|---|
| `run_coloc_shortlist.R` | Colocalizes the cis Tier 1 edges; `compare_arms.R` folds the result into `mr_tiered_edges.tsv` | `coloc_shortlist_results.tsv` |
| `run_coloc_systematic.R` | Colocalizes every cis edge still marked `pending`; `build_mr_tables.R` reads the result back | `coloc_results.tsv` |

`run_coloc_systematic.R` accepts `--nchunks N --chunk i` to split the work across jobs.
Chunks write per-locus files only; a final run without `--nchunks` reuses them and
writes `coloc_results.tsv`.

## Regional plots

These three steps rebuild the per-SNP tables behind the colocalization locus plot
(supplementary figure, `scripts/visualizations/figures/fig_mr_coloc.R`) and the
regional view on heap.bio. They do not change any colocalization result.

| Script | Does | Writes |
|---|---|---|
| `build_coloc_manifest.R` | One row per cis locus, with its pQTL file, outcome file and lead SNP resolved | `coloc_manifest.tsv` |
| `run_coloc_locus.R` (via `slurm/coloc/HEAPcoloc_locus.sh`) | Reruns coloc.abf for each manifest row and keeps the harmonized SNP table | `per_locus/<arm>__<protein>__<outcome>_*.tsv` |
| `export_coloc_web.R` | Adds LD (r²) to the lead variant with PLINK and the 1000G EUR panel | `web/<locus>_plot_table.tsv`, `_genes.tsv`, `_meta.tsv` |

The LD reference is GRCh37 while the summary statistics are GRCh38, so variants are
matched to the panel by rsID, not by position.

## Inputs

UKB-PPP and deCODE pQTL summary statistics, FinnGen R12 summary statistics and
manifest, the exposure GWAS (`scripts/gwas_regenie`), and a 1000G EUR PLINK panel
(regional plots only). Sources and environment variables:
[`docs/ENVIRONMENT.md`](../../../docs/ENVIRONMENT.md).
