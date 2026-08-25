# HEAP MR triad connections & provenance

How the **E → P → D** triads tested by Module 5 are assembled, which file makes
each connection, and how the triad count is derived. The goal: every triad is
auditable back to a named Module 2 / Module 3 output and the statistic that
qualified it.

A triad has three legs. Each leg is a separate causal claim, instrumented and
tested independently by MR, but **selected** from a specific upstream output:

```
        Module 2 (associations)                Module 3 (partitioned mediation)
        ───────────────────────                ────────────────────────────────
  E ───► P     replicated E→P assoc      P ───► D   PXS_<category> NIE significant
  (which exposures act on which          (which proteins mediate which disease,
   protein, sigBOTH train+test)           for an exposome category)
        │                                          │
        └──────────────── joined on (category, protein) ───────────────┘
                                   ▼
                         mr_triads.tsv  (E, category, P, D)
```

## The two connector files

### Leg E→P  —  `output/module2/ReplicatedEassoc.csv`  (the Module 2 connector)
The **dedicated Module 2 output that connects to MR.** Produced by
[`summarize_replicated_associations.R`](../module2_associations/summarize_replicated_associations.R):

- Module 2 fits univariate exposure→protein associations **separately in a train
  and a test split** (`output/module2/<experiment>/<covar>/univar_assoc_*.rds`,
  each with `$train`/`$test`).
- The producer aggregates the `statE` component over all proteins, outer-joins
  train×test, sets a Bonferroni threshold `0.05 / (#tested associations)`, and
  keeps rows significant in **both** splits → "replicated = sigBOTH".
- Key columns used downstream: `omicID` (protein), `Eid_train` (exposure
  variable = REGENIE id), `Category_train` (exposome category).

> **Freshness matters.** This file is covarType-specific and gets overwritten by
> whichever run produced it. It must reflect the canonical **base** run. A stale
> `Type3` copy (29 rows, from a near-empty partial run) once throttled the triad
> count to ~20; the current base run (400 proteins) yields **25,403** replicated
> E→P associations. Regenerate after any Module 2 rerun:
> `Rscript scripts/module2_associations/summarize_replicated_associations.R base`

### Legs P→D and E→P→D  —  Module 3 partitioned mediation (the Module 3 connector)
`output/module3/<experiment>/<covar>/lasso/partitioned_categories/MDres_*.txt`
(main experiments resolved from `config/modules/module3_experiments.yml`).

Each row is one `(protID, DZ_ID, predictor)` mediation estimate. The predictors
that matter:

| predictor | meaning | selects |
|-----------|---------|---------|
| `PXS_<category>` | the exposome category's proteomic-score component | **P→D pairs on a real E→P→D chain** (exposome NIE) → `edges_PD` |
| `Gcis_raw`, `Gtrans_raw` | the protein's genetic component | alternative genetic-anchored P→D set → `edges_PD_anyG` (provenance only) |

A protein-disease pair enters the canonical **PD/DP** set if a `PXS_<category>`
**NIE** (natural indirect effect) is significant for it — Bonferroni per
predictor per experiment, unioned across covariate specs. We use the **exposome**
NIE (not the genetic NIE) so the MR validates the exposome→protein→disease
mediation that HEAP discovered.

## Assembly (`Module5_load.R`) → canonical edge files

| Output | Built from | Definition |
|--------|-----------|------------|
| `edges_PD.tsv` / `edges_DP.tsv` | M3 PXS NIE | protein↔disease on a significant E→P→D chain |
| `edges_EP.tsv` / `edges_PE.tsv` | M2 `ReplicatedEassoc` ∩ priority proteins/categories | exposure↔protein |
| `edges_ED.tsv` / `edges_DE.tsv` | M3 triples × M2 E→P map | exposure↔disease |
| `mr_triads.tsv` | join of the above on (category, protein) | fully-specified E,P,D triads |

Disease ids are FinnGen R12 (`finngen_R12_*`, mapped via
`FinnGen/UKBFinnGenDisease.csv`); exposures are REGENIE ids; proteins are gene
symbols (mapped to Olink / deCODE SomaScan inside the runners).

## The triad funnel (current counts)

```
Module 2 base run ........... 400 proteins, 702,552 E→P tests
  └─ replicated (sigBOTH) ... 25,403 E→P associations         (ReplicatedEassoc.csv)
Module 3 partitioned ........ 13,545,600 mediation rows
  └─ PXS NIE significant ..... 8,129 E(category)→P→D triples  (edges_EPD_categories)
       └─ × FinnGen map ...... 5,261 protein↔disease           (edges_PD / DP)
Join E→P × E(cat)→P→D on (category, protein)
  └─ fully-specified triads .. 77,063                          (mr_triads.tsv)
       spanning 114 exposures × 665 proteins × 71 diseases
Edges actually MR-tested:
  EP/PE = 18,794   PD/DP = 5,261   ED/DE = 3,687
```

**FinnGen multi-phenotype mapping.** `UKBFinnGenDisease.csv` maps some UKB
diseases to several FinnGen phenotypes (`;`-separated, e.g. `T2D; T2D_WIDE`,
`E4_OBESITY; E4_OBESITYCAL; E4_OBESITYNAS` — 13 of 135 rows). `Module5_load.R`
splits these into one row per phenotype, so each becomes a distinct,
file-resolvable `<finngen_R12_*>.gz` outcome (this is why P–D = 5,261, not 3,665,
and diseases = 71, not 60).

## Instrument coverage (verified after GWAS completion)

Every triad entity was checked against the actual GWAS files the runners read:

| Source | Coverage | Notes |
|--------|----------|-------|
| Exposure GWAS (REGENIE) | 119/120 | 1 degenerate multi-select activity one-hot self-skips |
| Disease GWAS (FinnGen)  | 71/71 ✓ | all resolve after the multi-id split |
| UKB Olink pQTL          | 678/681 | IL6, TNF map to 4 Olink assays each (runner needs 1:1); HLA_E not on the panel |
| deCODE SomaScan pQTL    | 525/681 (77%) | cross-platform limit — uncovered proteins simply don't replicate in the deCODE arm |

The handful of unmapped entities fail gracefully per-edge (`safe_edge` logs and
continues); they do not block the run.

Every count above is reproduced by re-running, in order:
1. `summarize_replicated_associations.R base`  (refresh the Module 2 connector)
2. `Module5_load.R`                            (rebuild edges + mr_triads)
3. `generate_module5_manifest(<experiment>)`   (resize the array)
