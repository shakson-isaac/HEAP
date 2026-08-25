# slurm/module5 — Module 5 (Mendelian Randomization) job launchers

SLURM launchers for the HEAP MR framework. The analysis/pipeline details live in
[`scripts/module5_mr/README.md`](../../scripts/module5_mr/README.md) and
[`TRIAD_CONNECTIONS.md`](../../scripts/module5_mr/TRIAD_CONNECTIONS.md); this file
is the **how-to-run-the-jobs** reference.

## Scripts

| Script | Purpose |
|--------|---------|
| `HEAPmodule5_prewarm.sh` | Pre-compute (clump) every unique instrument once into the shared caches, so the main array starts warm. Run **first**. |
| `HEAPmodule5_manifest.sh` | The main MR array. Dispatches the per-row `runner` (UKB `Module5.R` vs deCODE `Module5_deCODE.R`) over the manifest. |
| `HEAPmodule5v2.sh` | **SUPERSEDED** — hardcoded UKB-only, fixed n_chunks=5000, not config-driven. Kept for reference; do not use. |

Both active scripts are **self-submitting** (run them on a login/compute node; they
`sbatch` themselves as array jobs) and write logs to
`/n/groups/patel/IGLOO/UKB/HEAP/logs/module5/`. Outputs and the umask are
group-writable (`umask 002`).

## Run order

```bash
cd /n/groups/patel/shakson_ukb/HEAP

# 0. (if triads/manifests are stale) rebuild the edge set + manifests
bash scripts/module5_mr/rebuild_triads.sh

# 1. Pre-warm the instrument caches (run FIRST). Note the 4 job ids it prints.
bash slurm/module5/HEAPmodule5_prewarm.sh

# 2. Submit the two MR arms, gated on the pre-warm completing.
#    Use afterok with all 4 prewarm ids (colon-separated), or afterany.
EXPERIMENT=MR_UKB_primary \
  DEPENDENCY=afterok:<ukb_protein_id>:<decode_protein_id>:<exposure_id>:<disease_id> \
  bash slurm/module5/HEAPmodule5_manifest.sh
EXPERIMENT=MR_deCODE_replication \
  DEPENDENCY=afterok:<...same ids...> \
  bash slurm/module5/HEAPmodule5_manifest.sh
```

> O2 ignores `SBATCH_DEPENDENCY`, so dependencies **must** be passed via the
> `DEPENDENCY` env var (the launchers translate it to `sbatch --dependency`).
> Skipping the pre-warm is allowed (drop `DEPENDENCY`): the runners still build
> caches on demand, just less efficiently.

## `HEAPmodule5_prewarm.sh`

Submits **one array per kind** — `ukb_protein`, `decode_protein`, `exposure`,
`disease` — each running `prewarm_instruments.R <kind> <slice_idx> <n_slices>`.
The runners are sourced with `HEAP_MR_WARM_ONLY=1`, which exposes their `CFG` +
`get_*_instruments()` but skips the edge run, so the cache is identical to what
the main runners read. `exposure`/`disease` caches are shared by both arms (warmed
via the UKB runner); only proteins are arm-specific (`protein_inst` vs
`protein_inst_decode`).

Defaults: `-t 0-03:00`, `--mem=24G`, `-c 2`, `-p short`; parallelism
`NSLICES_ukb_protein=80`, `_decode_protein=80`, `_exposure=12`, `_disease=8`
(180 tasks). Instrument universe ≈ 1,067 proteins (PD ∪ PE), ~120 exposures,
71 diseases; clumping ≈ 54 s each (PLINK reloads the 8.5M-variant EUR ref per
call). **Wall-clock ≈ 25–40 min** at default parallelism — bump `NSLICES_*` for
faster. Entities with no genome-wide-significant variant self-skip (no instrument).

Override examples:
```bash
NSLICES_ukb_protein=120 NSLICES_decode_protein=120 bash slurm/module5/HEAPmodule5_prewarm.sh
TIME=0-02:00 MEM=32G bash slurm/module5/HEAPmodule5_prewarm.sh
```

## `HEAPmodule5_manifest.sh`

```bash
EXPERIMENT=MR_UKB_primary bash slurm/module5/HEAPmodule5_manifest.sh
```
Reads `manifests/module5/<EXPERIMENT>.tsv`, self-submits `--array=1-<N_rows>`, and
each task resolves its row (`runner`, `edge_type`, `chunk_id`, `n_chunks`) and runs
`<runner> <edge_type> <chunk_id> <n_chunks>`. Defaults `-t 0-08:00`, `--mem=32G`,
`-c 4`, `-p short` (env-overridable: `TIME/MEM/CPUS/PARTITION`; `MAKE_PLOTS=1` for
per-edge scatter plots). With warm caches, tasks finish in minutes; cold PD/DP
chunks run ~1–2 h (covered by the 8 h ceiling).

Experiments (from `config/modules/module5_experiments.yml`): `MR_UKB_primary`
(split-sample UKB), `MR_deCODE_replication` (deCODE SomaScan), `MR_UKB_relaxed_p`
(sensitivity, pending — needs the runner to read `p_threshold`).

## Logs & monitoring

```bash
squeue -u $USER -o "%.12i %.20j %.8T %R" | grep -E "M5warm|M5_"
ls /n/groups/patel/IGLOO/UKB/HEAP/logs/module5/
```
- Pre-warm: `logs/module5/prewarm_<kind>_<jobid>_<task>.{out,err}`
- Main:     `logs/module5/<EXPERIMENT>_<jobid>_<task>.{out,err}`
- Per-chunk runner logs: `output/mr_edges{,_decode}/<EXPERIMENT>/logs/<edge_type>/chunk_<idx>.log`
  (the launcher exports `HEAP_MR_OUTDIR` = the manifest `output_path`, so edge
  outputs are experiment-keyed: `output/mr_edges/MR_UKB_primary/...`,
  `output/mr_edges_decode/MR_deCODE_replication/...`)
