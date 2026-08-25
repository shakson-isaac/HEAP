#!/bin/bash
# Combined Module 6 longitudinal PES + frozen Cox risk-calculator array.
# Each array task:
#   1. Re-runs Module6_prod_longitudinal.R for one exposure to refresh PES artifacts.
#   2. Runs Module6_longitudinal_cox_frozenrisk.R on the same exposure.

#SBATCH -t 1-00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=40G
#SBATCH -p medium
#SBATCH -J M6pes_frz
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
WORK_DIR="${HEAP_ROOT}/slurm/module6"
PES_SCRIPT="${HEAP_ROOT}/scripts/module6_prediction/Module6_prod_longitudinal.R"
COX_SCRIPT="${HEAP_ROOT}/scripts/module6_prediction/Module6_longitudinal_cox_frozenrisk.R"
IGLOO_HEAP="${HEAP_IGLOO_ROOT:-/n/groups/patel/IGLOO}/UKB/HEAP"
LOADER_RDS="${IGLOO_HEAP}/intermediate/HEAP.rds"

EXPOSURES_TXT="${EXPOSURES_TXT:-${HEAP_ROOT}/config/exposure_sets/module6_exposures.txt}"
COVAR_TYPE="${COVAR_TYPE:-base}"
RESULT_DIR="${RESULT_DIR:-${HEAP_IGLOO_ROOT:-/n/groups/patel/IGLOO}/UKB/HEAP/output/module6_pes_longitudinal/${COVAR_TYPE}}"
STATUS_DIR="${HEAP_IGLOO_ROOT:-/n/groups/patel/IGLOO}/UKB/HEAP/logs/status/module6"

SCORE_TYPES="${SCORE_TYPES:-prot,full}"
INSTANCES="${INSTANCES:-0,2,3}"
MIN_N="${MIN_N:-200}"
MIN_EVENTS="${MIN_EVENTS:-20}"
EVAL_MIN_EVENTS="${EVAL_MIN_EVENTS:-5}"
RISK_HORIZONS="${RISK_HORIZONS:-1,3,5,10}"
PES_SCALE="${PES_SCALE:-oof}"
# Empty by default so output stems are FrozenCoxEval / RiskCalculatorMetrics, which is
# the tag the fig_pes_cox_* figures + load_module6_frozenrisk() expect (default tag="").
OUTPUT_TAG="${OUTPUT_TAG:-}"

MAX_DISEASES="${MAX_DISEASES:-}"
DISEASE_REGEX="${DISEASE_REGEX:-}"
# Off by default: per-person RiskCalculatorScores is ~0.14 GB/exposure even at 20
# diseases (much larger over all diseases) and the figures use eval/fits/deciles, not it.
SAVE_RISK_SCORES="${SAVE_RISK_SCORES:-0}"
SAVE_PERSON_SCORES="${SAVE_PERSON_SCORES:-0}"
SAVE_MODELS="${SAVE_MODELS:-0}"
# Reuse the already-fit prod PES (skip the ~24h refit) when HoldoutScores.tsv exists.
SKIP_PES_IF_HOLDOUT_DONE="${SKIP_PES_IF_HOLDOUT_DONE:-1}"

mkdir -p "${STATUS_DIR}" "${RESULT_DIR}"

if [[ ! -s "${LOADER_RDS}" ]]; then
  echo "ERROR: missing HEAP.rds (run run_HEAP_loader.sh first): ${LOADER_RDS}" >&2
  exit 2
fi

if [[ ! -s "${PES_SCRIPT}" ]]; then
  echo "ERROR: missing PES script: ${PES_SCRIPT}" >&2
  exit 2
fi

if [[ ! -s "${COX_SCRIPT}" ]]; then
  echo "ERROR: missing frozen Cox script: ${COX_SCRIPT}" >&2
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

SAFE_EXPOSURE="$(printf '%s' "${EXPOSURE_ID}" | sed 's/[^A-Za-z0-9_.-]/_/g')"
PREFIX="${RESULT_DIR}/PESlong_${COVAR_TYPE}_${SAFE_EXPOSURE}"
HOLDOUT_TSV="${PREFIX}_HoldoutScores.tsv"
STATUS_TSV="${STATUS_DIR}/pes_frozenrisk_array_status.tsv"

echo "[$(date)] Combined task ${SLURM_ARRAY_TASK_ID} started on ${HOSTNAME}"
echo "[$(date)] covar=${COVAR_TYPE}"
echo "[$(date)] exposure=${EXPOSURE_ID}"
echo "[$(date)] result_dir=${RESULT_DIR}"
echo "[$(date)] pes_scale=${PES_SCALE}, risk_horizons=${RISK_HORIZONS}, output_tag=${OUTPUT_TAG}"

if [[ "${SKIP_PES_IF_HOLDOUT_DONE}" == "1" && -s "${HOLDOUT_TSV}" ]]; then
  echo "[$(date)] Skipping PES refresh because ${HOLDOUT_TSV} already exists"
else
  echo "[$(date)] Running PES longitudinal model"
  Rscript "${PES_SCRIPT}" "${COVAR_TYPE}" "${EXPOSURE_ID}"
  echo "[$(date)] PES longitudinal model completed"
fi

COX_ARGS=(
  "${COVAR_TYPE}"
  "${EXPOSURE_ID}"
  "--result-dir" "${RESULT_DIR}"
  "--score-types" "${SCORE_TYPES}"
  "--instances" "${INSTANCES}"
  "--min-n" "${MIN_N}"
  "--min-events" "${MIN_EVENTS}"
  "--eval-min-events" "${EVAL_MIN_EVENTS}"
  "--pes-scale" "${PES_SCALE}"
  "--horizons" "${RISK_HORIZONS}"
  "--output-tag" "${OUTPUT_TAG}"
)

if [[ -n "${MAX_DISEASES}" ]]; then
  COX_ARGS+=("--max-diseases" "${MAX_DISEASES}")
fi

if [[ -n "${DISEASE_REGEX}" ]]; then
  COX_ARGS+=("--disease-regex" "${DISEASE_REGEX}")
fi

if [[ "${SAVE_RISK_SCORES}" == "1" ]]; then
  COX_ARGS+=("--save-risk-scores")
fi

if [[ "${SAVE_PERSON_SCORES}" == "1" ]]; then
  COX_ARGS+=("--save-person-scores")
fi

if [[ "${SAVE_MODELS}" == "1" ]]; then
  COX_ARGS+=("--save-models")
else
  COX_ARGS+=("--no-save-models")
fi

echo "[$(date)] Running frozen Cox risk calculator"
echo "[$(date)] Rscript ${COX_SCRIPT} ${COX_ARGS[*]}"
Rscript "${COX_SCRIPT}" "${COX_ARGS[@]}"
echo "[$(date)] Frozen Cox risk calculator completed"

{
  flock 9
  if [[ ! -s "${STATUS_TSV}" ]]; then
    printf "timestamp\tjob_id\tarray_task_id\tcovar_type\texposure_id\toutput_tag\tstatus\n"
  fi
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$(date --iso-8601=seconds)" "${SLURM_ARRAY_JOB_ID:-NA}" "${SLURM_ARRAY_TASK_ID}" "${COVAR_TYPE}" "${EXPOSURE_ID}" "${OUTPUT_TAG}" "completed"
} 9>>"${STATUS_TSV}"

echo "[$(date)] Combined task ${SLURM_ARRAY_TASK_ID} completed"
