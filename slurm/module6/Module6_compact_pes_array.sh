#!/bin/bash
#SBATCH -t 0-00:30
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=30G
#SBATCH -p short
#SBATCH -J M6compact
#SBATCH --array=1-169
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
R_SCRIPT="${HEAP_ROOT}/scripts/module6_prediction/Module6_compact_pes_deployability.R"
EXPOSURES_TXT="${EXPOSURES_TXT:-${HEAP_ROOT}/config/exposure_sets/module6_exposures.txt}"
COVAR_TYPE="${COVAR_TYPE:-base}"

KS="${KS:-10,25,50,100,200,500,all}"
SELECTION_MODES="${SELECTION_MODES:-lasso,portable,portable_weighted}"
PORTABILITY_FILE="${PORTABILITY_FILE:-/n/groups/patel/IGLOO/UKB/OlinkSoma/OlinkSoma.csv}"
PORTABILITY_CORR_COL="${PORTABILITY_CORR_COL:-olink_smpnorm_corr}"
PORTABILITY_THRESHOLD="${PORTABILITY_THRESHOLD:-0.70}"
PORTABILITY_THRESHOLDS="${PORTABILITY_THRESHOLDS:-0.50,0.60,0.70,0.80,0.90}"
OUTPUT_TAG="${OUTPUT_TAG:-CompactPESThresholdSweep}"

mkdir -p "${HEAP_IGLOO_ROOT:-/n/groups/patel/IGLOO}/UKB/HEAP/logs/slurm/module6"
mkdir -p "${HEAP_IGLOO_ROOT:-/n/groups/patel/IGLOO}/UKB/HEAP/logs/status/module6"

if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
  echo "ERROR: submit this script as a Slurm array job. Default subset: --array=1,5,9,10,38,40,50,81,82,83,90,100,101,102,106,110,112,120,131%8 ${0}" >&2
  exit 2
fi

EXPOSURE_ID="$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${EXPOSURES_TXT}" | tr -d '\r' | sed 's/[[:space:]]*$//')"
if [[ -z "${EXPOSURE_ID}" ]]; then
  echo "ERROR: exposure_id empty for task ${SLURM_ARRAY_TASK_ID}; check ${EXPOSURES_TXT}" >&2
  exit 2
fi

ARTIFACT="${HEAP_IGLOO_ROOT:-/n/groups/patel/IGLOO}/UKB/HEAP/output/module6_pes_longitudinal/${COVAR_TYPE}/PESlong_${COVAR_TYPE}_${EXPOSURE_ID//[^A-Za-z0-9_.-]/_}_FinalModelArtifact.rds"
if [[ ! -s "${ARTIFACT}" ]]; then
  echo "ERROR: missing required final PES artifact for ${EXPOSURE_ID}: ${ARTIFACT}" >&2
  echo "This compact workflow compresses existing trained PES models; run Module6_prod_longitudinal.R first for this exposure." >&2
  exit 3
fi

echo "[$(date)] Compact PES task ${SLURM_ARRAY_TASK_ID}: covar=${COVAR_TYPE}, exposure=${EXPOSURE_ID}, host=${HOSTNAME}"
echo "[$(date)] ks=${KS}; selection_modes=${SELECTION_MODES}; portability=${PORTABILITY_FILE}; thresholds=${PORTABILITY_THRESHOLDS}"

Rscript "${R_SCRIPT}" "${COVAR_TYPE}" "${EXPOSURE_ID}" \
  --ks "${KS}" \
  --selection-modes "${SELECTION_MODES}" \
  --model-types prot \
  --portability-file "${PORTABILITY_FILE}" \
  --portability-corr-col "${PORTABILITY_CORR_COL}" \
  --portability-threshold "${PORTABILITY_THRESHOLD}" \
  --portability-thresholds "${PORTABILITY_THRESHOLDS}" \
  --output-tag "${OUTPUT_TAG}"

echo "[$(date)] Compact PES task ${SLURM_ARRAY_TASK_ID} completed"
