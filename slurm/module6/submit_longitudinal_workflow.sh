#!/bin/bash
umask 002  # group-writable outputs for hpc_patel multi-user runs

HEAP_ROOT="${HEAP_ROOT:-/n/groups/patel/shakson_ukb/HEAP}"
export HEAP_ROOT
export HEAP_PATHS_FILE="${HEAP_ROOT}/workflow/00_paths.R"
mkdir -p "${HEAP_IGLOO_ROOT:-/n/groups/patel/IGLOO}/UKB/HEAP/logs/slurm"

set -euo pipefail

# Module 6 reads HEAP.rds directly and derives the longitudinal PXS in-process via
# as_pxs_longitudinal() (00_paths.R) -- like modules 1/2/3. There is NO separate
# longitudinal loader step or intermediate RDS. Prerequisites: run
# slurm/loaders/run_HEAP_loader.sh (HEAP.rds) and Module6_config.R (exposure_specs.tsv +
# module6_exposures.txt) first.
BASE_DIR="${HEAP_ROOT}"
WORK_DIR="${HEAP_ROOT}/slurm/module6"
STATUS_DIR="${HEAP_IGLOO_ROOT:-/n/groups/patel/IGLOO}/UKB/HEAP/logs/status/module6"
SLURM_DIR="${HEAP_IGLOO_ROOT:-/n/groups/patel/IGLOO}/UKB/HEAP/logs/slurm/module6"
IGLOO_HEAP="${HEAP_IGLOO_ROOT:-/n/groups/patel/IGLOO}/UKB/HEAP"
HEAP_RDS="${IGLOO_HEAP}/intermediate/HEAP.rds"

MODE="${1:-pilot}"
COVAR_TYPES="${COVAR_TYPES:-base}"
PILOT_COVAR_TYPE="${PILOT_COVAR_TYPE:-base}"
PILOT_EXPOSURE_ID="${PILOT_EXPOSURE_ID:-smoking_status_f20116_0_0_Current}"
EXPOSURES_TXT="${EXPOSURES_TXT:-${HEAP_ROOT}/config/exposure_sets/module6_exposures.txt}"

if [[ ! -s "${HEAP_RDS}" ]]; then
  echo "ERROR: HEAP.rds not found: ${HEAP_RDS} (run slurm/loaders/run_HEAP_loader.sh first)" >&2
  exit 2
fi
if [[ ! -s "${EXPOSURES_TXT}" ]]; then
  echo "ERROR: exposure list not found: ${EXPOSURES_TXT} (run Module6_config.R first)" >&2
  exit 2
fi

mkdir -p "${STATUS_DIR}" "${SLURM_DIR}"
JOB_IDS_FILE="${STATUS_DIR}/job_ids.tsv"
echo -e "step\tjob_id\tnote" > "${JOB_IDS_FILE}"

cd "${WORK_DIR}"

submit_pilot() {
  local pilot_job
  pilot_job="$(sbatch --parsable "${WORK_DIR}/Module6_longitudinal_pilot.sh" "${PILOT_COVAR_TYPE}" "${PILOT_EXPOSURE_ID}")"
  echo -e "pilot\t${pilot_job}\t${PILOT_COVAR_TYPE} ${PILOT_EXPOSURE_ID}" | tee -a "${JOB_IDS_FILE}"
}

submit_arrays() {
  local n_exp
  n_exp="$(wc -l < "${EXPOSURES_TXT}" | tr -d ' ')"
  for covar_type in ${COVAR_TYPES}; do
    local array_job
    array_job="$(sbatch --parsable --array="1-${n_exp}" \
      --export=ALL,COVAR_TYPE="${covar_type}",EXPOSURES_TXT="${EXPOSURES_TXT}" \
      "${WORK_DIR}/Module6_longitudinal_array.sh")"
    echo -e "array\t${array_job}\t${covar_type} ${EXPOSURES_TXT}" | tee -a "${JOB_IDS_FILE}"
  done
}

case "${MODE}" in
  pilot)
    submit_pilot
    ;;
  all)
    submit_arrays
    ;;
  *)
    echo "Usage: $0 {pilot|all}" >&2
    exit 2
    ;;
esac

echo
echo "Submitted workflow mode: ${MODE}"
echo "Job IDs: ${JOB_IDS_FILE}"
