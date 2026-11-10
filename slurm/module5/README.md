# MR job launchers

| Script | Does |
|---|---|
| `HEAPmodule5_prewarm.sh` | Optional. Builds the instrument caches, one array per kind (UK Biobank proteins, deCODE proteins, exposures, diseases) |
| `HEAPmodule5_manifest.sh` | The MR array for one experiment (`EXPERIMENT=MR_UKB_primary` or `MR_deCODE_replication`) |

Both submit themselves as SLURM arrays when run with `bash`. Resources can be overridden
with `TIME`, `MEM`, `CPUS` and `PARTITION`, and a job dependency with
`DEPENDENCY=afterok:<job id>`. Logs go to `$HEAP_IGLOO_ROOT/UKB/HEAP/logs/module5/`.

Run order and outputs: [`scripts/module5_mr/README.md`](../../scripts/module5_mr/README.md).
