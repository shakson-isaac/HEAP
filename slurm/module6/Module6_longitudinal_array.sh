#!/bin/bash
#SBATCH -t 1-00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=40G
#SBATCH -p medium
#SBATCH -J M6long_array
#SBATCH --array=1-172
umask 002  # group-writable outputs for hpc_patel multi-user runs


HEAP_ROOT="${HEAP_ROOT:-/n/groups/patel/shakson_ukb/HEAP}"
export HEAP_ROOT
export HEAP_PATHS_FILE="${HEAP_ROOT}/workflow/00_paths.R"
mkdir -p "${HEAP_IGLOO_ROOT:-/n/groups/patel/IGLOO}/UKB/HEAP/logs/slurm"

set -euo pipefail

module load gcc/14.2.0
module load R/4.4.2

ulimit -n 10000

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export R_THREADS=1

BASE_DIR="${HEAP_ROOT}"
R_SCRIPT="${HEAP_ROOT}/scripts/module6_prediction/Module6_prod_longitudinal.R"
IGLOO_HEAP="${HEAP_IGLOO_ROOT:-/n/groups/patel/IGLOO}/UKB/HEAP"
# HEAP.rds (module6 derives the longitudinal view in-process via as_pxs_longitudinal)
LOADER_RDS="${IGLOO_HEAP}/intermediate/HEAP.rds"

EXPOSURES_TXT="${EXPOSURES_TXT:-${HEAP_ROOT}/config/exposure_sets/module6_exposures.txt}"
COVAR_TYPE="${COVAR_TYPE:-base}"

mkdir -p "${HEAP_IGLOO_ROOT:-/n/groups/patel/IGLOO}/UKB/HEAP/logs/slurm/module6"
mkdir -p "${HEAP_IGLOO_ROOT:-/n/groups/patel/IGLOO}/UKB/HEAP/logs/status/module6"

if [[ ! -s "${LOADER_RDS}" ]]; then
  echo "ERROR: missing HEAP.rds (run run_HEAP_loader.sh first): ${LOADER_RDS}" >&2
  exit 2
fi

if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
  echo "ERROR: this script must be submitted as a Slurm array job." >&2
  exit 2
fi

EXPOSURE_ID="$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${EXPOSURES_TXT}" | tr -d '\r' | sed 's/[[:space:]]*$//')"
if [[ -z "${EXPOSURE_ID}" ]]; then
  echo "ERROR: exposure_id empty for task ${SLURM_ARRAY_TASK_ID}; check ${EXPOSURES_TXT}" >&2
  exit 2
fi

echo "[$(date)] Task ${SLURM_ARRAY_TASK_ID}: covar=${COVAR_TYPE}, exposure=${EXPOSURE_ID}, host=${HOSTNAME}"
Rscript "${R_SCRIPT}" "${COVAR_TYPE}" "${EXPOSURE_ID}"
echo "[$(date)] Task ${SLURM_ARRAY_TASK_ID} completed"
