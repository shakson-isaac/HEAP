# Mendelian randomization

Bidirectional two-sample MR across exposure → protein → disease triads (Fig. 4). The
directory is numbered 5 for historical reasons. Full map:
[Figures and code](../../README.md#figures-and-code).

## Two arms

| Arm | Experiment | Protein instruments | Runner |
|---|---|---|---|
| UK Biobank | `MR_UKB_primary` | UK Biobank Olink pQTL | `Module5.R` |
| deCODE replication | `MR_deCODE_replication` | deCODE SomaScan pQTL (no UK Biobank overlap) | `Module5_deCODE.R` |

Both arms use the same triads, the same exposure GWAS (regenie, run on UK Biobank
participants without proteomics, so exposure and pQTL samples do not overlap) and FinnGen
R12 disease GWAS. Each pair is tested in both directions, `EP PE PD DP ED DE`, with cis and
trans protein instruments tested separately.

## Run order

| Step | Command | Does |
|---|---|---|
| 1 | `bash scripts/module5_mr/rebuild_triads.sh` | Builds the triads: replicated exposure → protein associations from Module 2 (`summarize_replicated_associations.R`), joined to the protein → disease pairs with a significant indirect effect in Module 3 (`Module5_load.R`); then writes the job manifests (`generate_module5_manifest()` in `workflow/generate_manifests.R`) |
| 2 (optional) | `bash slurm/module5/HEAPmodule5_prewarm.sh` | Selects and clumps each instrument once (`prewarm_instruments.R`) so the array jobs start from a cache; without it the runners build the cache as they go |
| 3 | `EXPERIMENT=MR_UKB_primary bash slurm/module5/HEAPmodule5_manifest.sh`, then the same with `MR_deCODE_replication` | One array task per manifest row: `<runner> <edge_type> <chunk> <n_chunks>`. Pass `DEPENDENCY=afterok:<job id>` to wait for step 2 |
| 4 | `Rscript scripts/support/mr_tables/build_mr_tables.R UKB` (and `DECODE`) | Per-arm tables: Benjamini–Hochberg q, sensitivity flags, tiers, triad motifs |
| 5 | [`scripts/support/coloc/`](../support/coloc/README.md) | Colocalization of cis edges; rerun step 4 afterwards so the tiers use it |
| 6 | `Rscript scripts/support/mr_tables/compare_arms.R`, then `scripts/analysis_summaries/summarize_mr_triads.R` | Replication across arms (Tier 1+) and the triad motif counts |

The [demo](../../demo/README.md) runs steps 3–6 on synthetic summary statistics.

## Per-edge output

`output/mr_edges/MR_UKB_primary/<edge>/<from>/<to>/` (deCODE: `output/mr_edges_decode/`),
one file per analysis: the summary, IVW / MR-Egger / weighted median / weighted mode, Cochran's Q,
Egger intercept, and from `mr_sensitivity.R` Steiger directionality, single-SNP,
leave-one-out and MR-PRESSO results. Instrument caches are kept under `output/mr/` and are
reused by reruns.

Instruments are selected at p < 5 × 10⁻⁸; the threshold is set in the runners.
The tier and motif definitions are in the paper's Methods.
