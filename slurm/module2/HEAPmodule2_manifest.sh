#!/bin/bash
umask 002  # group-writable outputs for hpc_patel multi-user runs
# Manifest-driven Module 2 array job.
#
# USAGE
# -----
#   # Step 1: Generate the manifest (run once per experiment)
#   Rscript -e "
#     source('workflow/00_paths.R')
#     source('workflow/config_helpers.R')
#     source('workflow/generate_manifests.R')
#     generate_module2_manifest('M2_base_main')
#   "
#
#   # Step 2: Submit
#   EXPERIMENT=M2_base_main bash slurm/module2/HEAPmodule2_manifest.sh
#
# REQUIRED ENV VARS
# -----------------
#   EXPERIMENT    Name of experiment (must match config/modules/module2_experiments.yml)
#
# OPTIONAL ENV VARS
# -----------------
#   HEAP_ROOT     (default: /n/groups/patel/shakson_ukb/HEAP)
#   MANIFEST      Override manifest path
#   DEPENDENCY    sbatch --dependency spec to gate this run, e.g. afterok:<jobid>
#                 (SBATCH_DEPENDENCY is ignored on this cluster, so use this)

set -euo pipefail

HEAP_ROOT="${HEAP_ROOT:-/n/groups/patel/shakson_ukb/HEAP}"
export HEAP_ROOT
export HEAP_PATHS_FILE="${HEAP_ROOT}/workflow/00_paths.R"

EXPERIMENT="${EXPERIMENT:-}"
if [[ -z "${EXPERIMENT}" ]]; then
  echo "ERROR: EXPERIMENT env var is required." >&2
  echo "  Example: EXPERIMENT=M2_base_main bash ${0}" >&2
  exit 1
fi

IGLOO_HEAP="${HEAP_IGLOO_ROOT:-/n/groups/patel/IGLOO}/UKB/HEAP"
MANIFEST="${MANIFEST:-${IGLOO_HEAP}/manifests/module2/${EXPERIMENT}.tsv}"

if [[ ! -f "${MANIFEST}" ]]; then
  echo "ERROR: Manifest not found: ${MANIFEST}" >&2
  echo "Generate it first:" >&2
  echo "  Rscript -e \"source('workflow/00_paths.R'); source('workflow/config_helpers.R'); source('workflow/generate_manifests.R'); generate_module2_manifest('${EXPERIMENT}')\"" >&2
  exit 1
fi

N_ROWS=$(tail -n +2 "${MANIFEST}" | wc -l)
if [[ "${N_ROWS}" -eq 0 ]]; then
  echo "ERROR: Manifest has no data rows: ${MANIFEST}" >&2
  exit 1
fi

RSCRIPT="${HEAP_ROOT}/scripts/module2_associations/Module2.R"
LOG_DIR="${IGLOO_HEAP}/logs/module2"
mkdir -p "${LOG_DIR}"

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Submitting Module 2 manifest job"
echo "  Experiment : ${EXPERIMENT}"
echo "  Manifest   : ${MANIFEST}"
echo "  N rows     : ${N_ROWS}"

# Array spec: ARRAY_TASKS overrides with an explicit list/range (e.g. failed
# chunk indices to rerun: "289" or "12,40,103-110"); otherwise full 1-N_ROWS.
ARRAY_SPEC="${ARRAY_TASKS:-1-${N_ROWS}}"
echo "  Array spec : ${ARRAY_SPEC}"

# Submitter vs array task: SLURM_ARRAY_TASK_ID is set ONLY in the array element.
# (SLURM_JOB_ID is also set in interactive/OnDemand allocations, e.g. a VS Code
# tunnel, so guarding on it misfires when launching from inside one.)
if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
  # Optional dependency gate. SBATCH_DEPENDENCY is ignored on this cluster, so
  # pass --dependency explicitly via the DEPENDENCY env var, e.g.
  #   DEPENDENCY=afterok:<jobid> EXPERIMENT=... bash <this script>
  dep_arg=()
  [[ -n "${DEPENDENCY:-}" ]] && dep_arg=(--dependency="${DEPENDENCY}")
  exec sbatch \
    -t "${TIME:-0-02:00}" \
    --ntasks=1 \
    --cpus-per-task=1 \
    --mem=20G \
    -p short \
    --array="${ARRAY_SPEC}" \
    --export=ALL \
    "${dep_arg[@]}" \
    -J "M2_${EXPERIMENT}" \
    -o "${LOG_DIR}/${EXPERIMENT}_%A_%a.out" \
    -e "${LOG_DIR}/${EXPERIMENT}_%A_%a.err" \
    --mail-type=FAIL \
    "$0"
fi

# -----------------------------------------------------------------------
# Job body
# -----------------------------------------------------------------------

module load gcc/14.2.0
module load R/4.4.2
ulimit -n 10000

IDX="${SLURM_ARRAY_TASK_ID}"
echo "[$(date '+%H:%M:%S')] START  experiment=${EXPERIMENT}  array_index=${IDX}"

Rscript "${RSCRIPT}" \
  --manifest "${MANIFEST}" \
  --array-index "${IDX}"

STATUS=$?
echo "[$(date '+%H:%M:%S')] FINISH  experiment=${EXPERIMENT}  array_index=${IDX}  exit=${STATUS}"
exit "${STATUS}"
