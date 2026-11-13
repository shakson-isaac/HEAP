#!/bin/bash
#SBATCH -c 4
#SBATCH -t 0-02:00
#SBATCH --mem=24G
#SBATCH -p short
#
# Convert the REGENIE exposure GWAS into the public deposit: one bgzipped file
# per exposure plus a tabix index.
#
#   sbatch --array=1-169 stage_gwas_deposit_array.sh
#
# Array size is the number of .regenie files, which is NOT the number of
# directories: three exposure directories are empty (former_alcohol_drinker,
# chronotype "definitely evening", previously_smoked_most_days) because those
# GWAS produced no output. Count before submitting:
#
#   find /n/groups/patel/IGLOO/UKB/HEAP/output/gwas/regenie_step2 \
#        -name '*.regenie' | wc -l
#
# One exposure takes ~35 s and 302 MB, measured on
# age_first_had_sexual_intercourse. The whole set is ~51 GB.
#
# Each task converts ONE exposure, so a failure costs one file rather than the
# run. The manifest and README are written once at the end by --finalize, after
# the array drains:
#
#   sbatch --dependency=afterok:<ARRAY_JOBID> stage_gwas_deposit_array.sh --finalize
#
set -euo pipefail

export PATH=/n/app/htslib/1.21-gcc-14.2.0/bin:/n/app/R/4.4.2-gcc-14.2.0/bin:$PATH
export LD_LIBRARY_PATH=/n/app/gcc/14.2.0/lib64:${LD_LIBRARY_PATH:-}

SCRIPT=/n/groups/patel/shakson_ukb/HEAP/scripts/gwas_regenie/stage_gwas_deposit.R
OUT=${HEAP_DEPOSIT_OUT:-/n/groups/patel/IGLOO/UKB/HEAP/output/gwas_deposit}
mkdir -p "$OUT"

# Per-task scratch. SLURM_TMPDIR is unset on O2, so the fallback is a path that
# is unique per user AND per job -- a bare /tmp/heap would collide across the
# array and two tasks would overwrite each other's intermediate TSV.
export TMPDIR="${SLURM_TMPDIR:-/tmp/$USER-${SLURM_JOB_ID:-$$}-${SLURM_ARRAY_TASK_ID:-0}}"
mkdir -p "$TMPDIR"
trap 'rm -rf "$TMPDIR"' EXIT

if [[ "${1:-}" == "--finalize" ]]; then
  Rscript "$SCRIPT" --out "$OUT" --finalize
  exit 0
fi

Rscript "$SCRIPT" --out "$OUT" --index "${SLURM_ARRAY_TASK_ID:?set --array}"
