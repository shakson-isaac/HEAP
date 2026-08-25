# Module 4 — Tissue / Pathway Enrichment

> **📍 Numbering (the code `module4` label is overloaded):** *this* dir is the
> **un-numbered** enrichment analysis (feeds Fig 3 biology + FigS8/S9). The manuscript's
> **Module 4 / Fig 4 is MR** (= code [`module5_mr/`](../module5_mr/)), and the intervention
> figures built under the `module4` namespace (`build_module4_composite*`) are the
> manuscript's **Module 5 / Fig 5**. Full map:
> [Manuscript vs. code module numbers](../../README.md#manuscript-vs-code-module-numbers).

A dedicated **analysis module** (not a visualization step) that tests whether the
proteins associated with each exposure (Module 2 output) are enriched for
tissue-specific or pathway gene sets. Outputs are canonical IGLOO tables that
visualization scripts then plot.

> Status: **scaffold** (2026-06-03). Stubs define the I/O contract and name the
> legacy source to port from. Not yet implemented — see
> `docs/SUPPORT_ANALYSIS_PLAN.md §2`.

## Inputs (all IGLOO-canonical via `workflow/00_paths.R` helpers)
- Module 2 significant exposures → `heap_project_output("module2", "ReplicatedEassoc.csv")`
  (and/or per-batch `module2/<covarType>/univar_assoc_*.rds`)
- GTEx v10 median-TPM GCTs → `heap_gtex_rna_dir()`  (`/n/groups/patel/IGLOO/UKB/GTEX/RNA`)
- HPA tables → `heap_hpa_or_legacy("normal_ihc_data.tsv" | "subcellular_location.tsv")`
- OmicsPred protein map (universe) → `heap_omicspred_or_legacy(...)`

## Outputs → `heap_project_output("module4_enrichment", ...)`
```
genesets/{GTEX_tissue,HPA_specific,HPA_enriched,HPA_expressed,HPA_secretome,KEGG_T2G,KEGG_T2N}.txt
GTEX_tau_scores.csv
OlinkEntrezConv.txt          # cached BioMart map
HEAPgsea.qs                  # GSEA result object
tissue_enrichment.csv        # figure-input (NES, p.adjust per exposure x tissue)
pathway_enrichment.csv       # figure-input (NES, p.adjust per exposure x pathway)
ora_tissue.csv / ora_pathway.csv   # ORA results (separate functionality)
```

## Scripts
| Script | Role | Ports from |
|---|---|---|
| `00_build_genesets.R` | GTEx Tau + FC>4 tissue Term2Gene; HPA sets; KEGG T2G/T2N | `ModuleExt/GTEX_ident.R`, `ModuleExt/HPA_ident.R` |
| `01_entrez_map.R` | one BioMart call → cached `OlinkEntrezConv.txt` | `ModuleExt/ObtainEntrezMapping.R` |
| `02_enrichment_core.R` | sourced fn library: `CreateGenelists`, `convertEntrez`, GSEA/ORA wrappers | `Module2/HEAPassoc_pathway.R`, `TissueSpec/TissueSpec_RunExample.R` |
| `03_run_gsea.R` | **GSEA** (primary): GTEx-tissue + Reactome `gsePathway` → `HEAPgsea.qs` | `Module2/HEAPassoc_pathway.R` |
| `03b_run_ora.R` | **ORA** (separate functionality): `enricher`/`enrichPathway`/`enrichGO` over sig protein lists → `ora_*.csv` | `TissueSpec/TissueSpec_RunExample.R` (ORA_tissue/ORA_paths) |
| `04_enrichment_tables.R` | flatten `HEAPgsea.qs` → `tissue_enrichment.csv`, `pathway_enrichment.csv` | `Module2/HEAPassoc_pathwaytables.R` |
| `05_enrichment_figures.R` | ComplexHeatmap tissue/pathway figures + GTEx Tau density (plotting only) | `Module2/HEAPassoc_pathwayviz.R` |

Scope decision (2026-06-03): **GSEA is the primary/manuscript path** (`03_run_gsea.R`);
ORA is provided as a **separate, optional** script (`03b_run_ora.R`) for
functionality, not the default figure path.
