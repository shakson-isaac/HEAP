# Tissue and pathway enrichment

Tests whether the proteins associated with each exposure are enriched for tissue-specific
or pathway gene sets. Results feed Fig. 2b and the enrichment supplementary figures.

The directory is named `module4` for historical reasons. Fig. 4 (MR) is in
[`module5_mr/`](../module5_mr/), and the Fig. 5 panels are built by
`build_module4_composite*`. Full map: [Figures and code](../../README.md#figures-and-code).

## Run order

| Step | Script | Does |
|---|---|---|
| 1 | `00_build_genesets.R` | GTEx tissue (Tau and fold-change), Human Protein Atlas and KEGG gene sets |
| 2 | `01_entrez_map.R` | Gene symbol → Entrez ID map for the Olink proteins |
| 3 | `00b_prepare_gsea_input.R` | Ranked exposure → protein associations from the Module 2 output |
| 4 | `slurm/module4_enrichment/run_gsea.sh` | `03_run_gsea.R` (tissue and Reactome GSEA per exposure) → `04_enrichment_tables.R` (flat tables) → tissue and pathway heatmaps |
| 5 | `slurm/module4_enrichment/run_category_gsea.sh` | `05_category_enrichment.R`: GSEA per exposure category |
| 6 | `run_module1_r2_gsea.R`, `run_module1_greml_gsea.R`, `run_module1_category_gsea.R` | GSEA of proteins ranked by their Module 1 variance components (supplementary) |

`02_enrichment_core.R` holds the shared functions. `03b_run_ora.R` and
`run_module1_category_ora.R` run over-representation analysis as an alternative to GSEA;
the paper uses GSEA.

## Inputs and outputs

Inputs: the Module 2 associations, GTEx v10 median TPM, Human Protein Atlas tables and the
OmicsPred protein list (sources in [`docs/DATA_ACCESS.md`](../../docs/DATA_ACCESS.md)).
Outputs go to `$HEAP_IGLOO_ROOT/UKB/HEAP/output/module4_enrichment/`: the gene sets, the
GSEA result object `HEAPgsea.qs`, and `tissue_enrichment.csv` / `pathway_enrichment.csv`
(one row per exposure × tissue or pathway).
