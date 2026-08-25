#!/bin/bash
#SBATCH -c 1
#SBATCH -t 0-02:00
#SBATCH --mem=16G
#SBATCH -p short
###############################################################################
# ldsc_h2_array.sh — LD Score Regression heritability + intercept, one exposure
#                    per array task.
#
# Array size is set at submission time:
#   sbatch --array=1-<N> ldsc_h2_array.sh
# where <N> = number of lines in slurm/ldsc/ldsc_exposures.txt. The launcher
# slurm/ldsc/submit_ldsc.sh regenerates that list (completed exposure GWAS) and
# sizes the array automatically.
#
# Each task just resolves its exposure id and calls the shared per-exposure
# pipeline scripts/ldsc/run_ldsc_h2.sh (preprocess -> munge -> ldsc.py --h2).
# Outputs: ${IGLOO}/UKB/HEAP/output/gwas/ldsc/{munged,h2}/<exposure>.*
#
# Uses the LDSC package + LD reference + conda env staged in IGLOO so any
# hpc_patel member can run it. Override any of those with the LDSC_* env vars
# documented in run_ldsc_h2.sh.
###############################################################################
set -euo pipefail
umask 0002   # group-writable outputs for hpc_patel team runs

HEAP_ROOT="${HEAP_ROOT:-/n/groups/patel/shakson_ukb/HEAP}"
export HEAP_ROOT
export HEAP_PATHS_FILE="${HEAP_ROOT}/workflow/00_paths.R"

mkdir -p "/n/groups/patel/IGLOO/UKB/HEAP/logs/ldsc"

EXP_LIST="${LDSC_EXP_LIST:-${HEAP_ROOT}/slurm/ldsc/ldsc_exposures.txt}"
if [[ ! -f "${EXP_LIST}" ]]; then
  echo "ERROR: exposure list not found: ${EXP_LIST}" >&2
  echo "Run: ${HEAP_ROOT}/slurm/ldsc/submit_ldsc.sh   (regenerates the list)" >&2
  exit 1
fi

Ename=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${EXP_LIST}" | tr -d '[:space:]')
if [[ -z "${Ename}" ]]; then
  echo "No exposure at array index ${SLURM_ARRAY_TASK_ID} in ${EXP_LIST}; exiting." >&2
  exit 0
fi

echo "[ldsc] task ${SLURM_ARRAY_TASK_ID}: ${Ename}"
exec bash "${HEAP_ROOT}/scripts/ldsc/run_ldsc_h2.sh" "${Ename}"
