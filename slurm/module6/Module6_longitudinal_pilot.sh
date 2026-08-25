#!/bin/bash
#SBATCH -t 1-00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=180G
#SBATCH -p medium
#SBATCH -J M6long_pilot
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
LOADER_RDS="${IGLOO_HEAP}/intermediate/HEAP.rds"

COVAR_TYPE="${1:-${COVAR_TYPE:-base}}"
EXPOSURE_ID="${2:-${EXPOSURE_ID:-smoking_status_f20116_0_0_Current}}"

mkdir -p "${HEAP_IGLOO_ROOT:-/n/groups/patel/IGLOO}/UKB/HEAP/logs/slurm/module6"
mkdir -p "${HEAP_IGLOO_ROOT:-/n/groups/patel/IGLOO}/UKB/HEAP/logs/status/module6"

if [[ ! -s "${LOADER_RDS}" ]]; then
  echo "ERROR: missing HEAP.rds (run run_HEAP_loader.sh first): ${LOADER_RDS}" >&2
  exit 2
fi

echo "[$(date)] Starting Module 6 longitudinal pilot on ${HOSTNAME}"
echo "Covariate type: ${COVAR_TYPE}"
echo "Exposure: ${EXPOSURE_ID}"
echo "R script: ${R_SCRIPT}"

Rscript "${R_SCRIPT}" "${COVAR_TYPE}" "${EXPOSURE_ID}"

echo "[$(date)] Pilot completed successfully"
