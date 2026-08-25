# Module 5 — Mendelian Randomization (HEAP MR framework)

> **📍 Numbering:** This is **code** Module 5. In the **manuscript** it is **Module 4 /
> Figure 4** — MR is presented *before* interventions. Full code↔paper map:
> [Manuscript vs. code module numbers](../../README.md#manuscript-vs-code-module-numbers).

Bidirectional MR across the HEAP triad **E → P → D** (Exposure → Protein →
Disease), run with two protein-instrument arms that share an identical triad
set so results are directly comparable.

## The two arms

| Arm | Experiment | Protein instruments | Exposure instruments | Disease outcomes | Runner |
|-----|-----------|---------------------|----------------------|------------------|--------|
| **Split-sample UKB** | `MR_UKB_primary` | UKB Olink pQTL (proteomics cohort) | UKB REGENIE GWAS trained on **non-proteomics** individuals | FinnGen R12 | `Module5.R` |
| **deCODE replication** | `MR_deCODE_replication` | deCODE SomaScan pQTL (external, zero UKB overlap) | same UKB non-proteomics REGENIE GWAS | FinnGen R12 | `Module5_deCODE.R` |

"Split sample" = the exposure GWAS and the protein pQTLs come from
**non-overlapping UKB subsamples** (exposure GWAS excludes everyone with
proteomics; pQTLs are measured only in the proteomics cohort), so E↔P MR is
two-sample valid. deCODE is the fully external pQTL replication. Both arms read
the **same** `global_edges/edges_<TYPE>.tsv` triad lists.

Sensitivity: `MR_UKB_relaxed_p` (p<5e-6 instruments; status `pending`).

## Edge types (full bidirectional network)

`EP, PE, PD, DP, ED, DE` — each undirected triad pair is tested in both
directions (forward causal + reverse-causation check). cis and trans protein
instruments are tested separately for the P→* edges.

## Pipeline

```
Module 3 (partitioned mediation, MDres_*.txt)
        │  significant NIE per predictor (Bonferroni per predictor per experiment)
        ▼
Module5_load.R  ──►  global_edges/edges_<TYPE>.tsv   (the organized triad set)
        │
        ▼
generate_module5_manifest()  ──►  manifests/module5/<EXPERIMENT>.tsv
        │  (per-edge-type chunking from actual edge counts; `runner` column)
        ▼
HEAPmodule5_manifest.sh  ──►  Module5.R / Module5_deCODE.R per row
        │
        ▼
output/mr_edges/MR_UKB_primary/<edge>/<from>/<to>/<prefix>_summary.tsv  (+ mr_methods, harmonised, heterogeneity, pleiotropy, steiger, singlesnp, leaveoneout, presso)
output/mr_edges_decode/MR_deCODE_replication/<edge>/...                  (deCODE arm)
   ^ edge-output base = the manifest output_path (HEAP_MR_OUTDIR), so it is experiment-keyed
```

> **One-shot rebuild:** `bash scripts/module5_mr/rebuild_triads.sh` runs the
> three connector steps in order (refresh `ReplicatedEassoc.csv` from base →
> `Module5_load.R` → `generate_module5_manifest`). The steps below document each
> stage individually.

### 1. Build the triad/edge set — `Module5_load.R`

Reads the **main partitioned** Module 3 experiments resolved from
`config/modules/module3_experiments.yml` (priority `main`, status `ready`,
`mediation_mode: partitioned_categories`) at their experiment-nested output
paths `module3/<experiment>/<covar>/<family>/<mode>/`. Applies a Bonferroni
threshold per predictor per experiment to the NIE estimates, unions
significant hits across covariate specs, maps UKB diseases → FinnGen, and
writes edge lists to `output/mr_edges/global_edges/`.

Run:
```bash
module load gcc/14.2.0 R/4.4.2
HEAP_ROOT=/n/groups/patel/shakson_ukb/HEAP \
HEAP_PATHS_FILE=$HEAP_ROOT/workflow/00_paths.R \
Rscript scripts/module5_mr/Module5_load.R
```

**Canonical runner-facing files** (what the runners read; FinnGen disease ids):

| File(s) | Columns | Definition |
|---------|---------|------------|
| `edges_PD.tsv` / `edges_DP.tsv` | `Protein, Disease` | protein↔disease pairs lying on a significant **E→P→D** chain, i.e. a **PXS_\<category\> (exposome) NIE** is significant for the pair. The MR validates the **exposome** mediation (E,P,D), not genetic mediation. |
| `edges_EP.tsv` / `edges_PE.tsv` | `Exposure, ExposureCategory, Protein` | replicated exposure→protein associations (`ReplicatedEassoc.csv`) within the priority set |
| `edges_ED.tsv` / `edges_DE.tsv` | `Exposure, Disease` | category-level E→P→D triples expanded to actual exposure variables via the replicated E→P map |

**Triad inventory:** `mr_triads.tsv` — columns `Exposure, ExposureCategory,
Protein, Disease(FinnGen), Disease_UKB, ICD10`. One row per fully-specified
E→P→D triad being tested (E→P replicated in Module 2 AND PXS_\<category\> NIE
significant in Module 3). This is the documented "what is being tested" table.
See [`TRIAD_CONNECTIONS.md`](TRIAD_CONNECTIONS.md) for the full provenance/funnel.

Provenance/analytic files (not read by the runners): `MR_priority_table.tsv`,
`edges_PD_anyG{,_finngen}.tsv` (the **genetic** Gcis/Gtrans-NIE-anchored P→D set, an
alternative selection), `edges_PD_Gcis/Gtrans{,_finngen}.tsv`,
`edges_EPD_categories{,_finngen}.tsv`, `edges_PD_all_priority_ukb.tsv`,
`HEAPres.tsv` (legacy alias of the triad inventory, read by old viz scripts).

Current counts: **PD/DP = 5,261** (exposome-mediated), **EP/PE = 18,794**,
**ED/DE = 3,687**; **mr_triads = 77,063** fully-specified triads (114 exposures ×
665 proteins × 71 diseases). NB: these depend on `ReplicatedEassoc.csv` being
current — it must be regenerated from the **base** run (not the stale Type3 copy).
UKB diseases that map to several FinnGen phenotypes (`;`-separated) are split into
separate outcomes. Instrument coverage (post-GWAS) is verified in `TRIAD_CONNECTIONS.md`.

### How `ReplicatedEassoc.csv` is built
Producer: [`summarize_replicated_associations.R`](../module2_associations/summarize_replicated_associations.R).
Module 2 fits univariate exposure→protein associations **separately in a train
and a test split** (per-protein `univar_assoc_<idx>.rds`, each with `$train`/`$test`).
The producer aggregates the `statE` component across all proteins, outer-joins
train×test, sets a Bonferroni threshold `0.05 / (#tested associations)`, and keeps
rows significant in **both** splits (`Pr(>|t|)_train` and `Pr(>|t|)_test` < thr) —
"replicated = sigBOTH". Output is covarType-specific (default `base`); re-run as
`Rscript .../summarize_replicated_associations.R base` after any Module 2 rerun.

### 2. Generate manifests

```r
source("workflow/00_paths.R"); source("workflow/config_helpers.R")
source("workflow/generate_manifests.R")
generate_module5_manifest("MR_UKB_primary")
generate_module5_manifest("MR_deCODE_replication")
```
Chunks are sized **per edge type** from the edge-list row counts
(`pairs_per_chunk`, default 25), so no array task gets an empty chunk. Each row
carries the `runner` to dispatch.

### 3a. (Recommended) Pre-warm the instrument caches

```bash
bash slurm/module5/HEAPmodule5_prewarm.sh
```
Clumps every **unique** instrument once (1,067 proteins × cis+trans for each arm,
~120 exposures, 71 diseases) into the shared `output/mr/*` caches, so the main
array doesn't redundantly clump the same protein in ~8 concurrent tasks. One
array per kind; `WARM_ONLY=1` makes the runners expose their functions without
running edges (so the cache is identical to what they read). Wall-clock ≈ **25–40
min** at the default parallelism (clumping ≈ 54 s each, gated by the ~13-protein
`ukb_protein` tasks); bump `NSLICES_*` for faster. Diseases/exposures with no
genome-wide-significant variant self-skip (no instrument).

### 3b. Submit the MR array (gate on the pre-warm)

```bash
# capture the prewarm job id(s), then:
EXPERIMENT=MR_UKB_primary        DEPENDENCY=afterok:<prewarm_jobid> bash slurm/module5/HEAPmodule5_manifest.sh
EXPERIMENT=MR_deCODE_replication DEPENDENCY=afterok:<prewarm_jobid> bash slurm/module5/HEAPmodule5_manifest.sh
```
(Without the pre-warm, just drop `DEPENDENCY` — the runners still build caches on
demand, only less efficiently.)
The launcher self-submits an array (`1..N_rows`), resolves each row, and runs
`<runner> <edge_type> <chunk_id> <n_chunks>`. Logs → `IGLOO/UKB/HEAP/logs/module5/`.
Gate on upstream with `DEPENDENCY=afterok:<jobid>` (O2 ignores `SBATCH_DEPENDENCY`).

## Sensitivity analyses (per edge)

`run_mr_edge()` emits, per edge: **IVW + MR-Egger + weighted median + weighted
mode** (`_mr_methods`), **Cochran-Q heterogeneity** (`_heterogeneity`), **Egger
intercept** (`_pleiotropy`). `mr_sensitivity.R` (shared by both runners) adds:

| File | Analysis | Guard |
|------|----------|-------|
| `_steiger.tsv` | Steiger directionality (reverse-causation) | needs outcome N |
| `_singlesnp.tsv` | per-SNP Wald ratios | nsnp ≥ 2 |
| `_leaveoneout.tsv` | leave-one-out IVW | nsnp ≥ 2 |
| `_presso.tsv` | MR-PRESSO global + outlier + distortion | nsnp ≥ 4 |

Each is best-effort (`tryCatch`) so it never blocks the main result. Steiger
needs a sample size for the binary FinnGen outcome (the `.gz` files have none),
so the disease reader injects case/control N from `FinnGen/finngen_R12_manifest.tsv`.
Design-level robustness also comes from cis/trans separation, the two independent
instrument arms (UKB↔deCODE), and the bidirectional reverse edges.

## Runtime (measured on a compute node, 4 cores)

Cold-cache (first run) is dominated by PLINK reloading the EUR LD reference on each
clump call: a P↔D pair (cis + trans clump + disease load + MR + all sensitivity)
runs ~150–250 s. So a 25-pair **PD/DP chunk ≈ 1–2 h**; **EP/PE/ED/DE chunks are
lighter** (the outcome side isn't clumped). The launcher's **`-t 0-08:00`, `32G`,
`4c`, `short`** comfortably covers this with margin and fits the 12 h short cap.
Reruns reuse the on-disk instrument caches and finish in minutes.

> Cold-run efficiency note: with the array un-throttled, many tasks clump the
> same protein concurrently (caches are shared but computed on demand), so the
> first run does redundant clumping. To avoid it, optionally pre-warm caches in a
> low-parallelism pass (clump each unique protein/disease/exposure once) before
> the full array — then every task starts warm.

## Status / gotchas

- **GWAS dependency:** the runners need the exposure REGENIE step-2 outputs,
  UKB/deCODE pQTL, and FinnGen sumstats. The triad set + manifests + launcher
  are ready now; MR results require the GWAS to finish.
- **Instruments are disk-cached** under `output/mr/{clumps,protein_inst,protein_inst_decode,disease_inst}/`
  — the first full run is the expensive one; reruns reuse caches.
- **`p_threshold` is hardcoded (5e-8)** in both runners; the manifest column is
  informational. `MR_UKB_relaxed_p` (5e-6) needs the runner parameterized
  before it takes effect.
- `slurm/module5/HEAPmodule5v2.sh` is **superseded** by the manifest launcher.
