#!/bin/bash
umask 002  # group-writable outputs for hpc_patel multi-user runs
# Manifest-driven Module 3 (Mediation) array job.
#
# USAGE
# -----
#   # Step 1: Generate the manifest (run once per experiment)
#   Rscript -e "
#     source('workflow/00_paths.R')
#     source('workflow/config_helpers.R')
#     source('workflow/generate_manifests.R')
#     generate_module3_manifest('M3_base_lasso_primary')
#   "
#
#   # Step 2: Submit the array job
#   EXPERIMENT=M3_base_lasso_primary bash slurm/module3/HEAPmodule3_manifest.sh
#
# REQUIRED ENV VARS
# -----------------
#   EXPERIMENT    Name of the experiment (must match a key in
#                 config/modules/module3_experiments.yml and have a
#                 generated manifest)
#
# OPTIONAL ENV VARS
# -----------------
#   HEAP_ROOT     (default: /n/groups/patel/shakson_ukb/HEAP)
#   MANIFEST      Override manifest path
#   RESUME_FROM   Array index to start from (for partial reruns)
#   DEPENDENCY    sbatch --dependency spec to gate this run, e.g. afterok:<jobid>
#                 (SBATCH_DEPENDENCY is ignored on this cluster, so use this)

set -euo pipefail

HEAP_ROOT="${HEAP_ROOT:-/n/groups/patel/shakson_ukb/HEAP}"
export HEAP_ROOT
export HEAP_PATHS_FILE="${HEAP_ROOT}/workflow/00_paths.R"

EXPERIMENT="${EXPERIMENT:-}"
if [[ -z "${EXPERIMENT}" ]]; then
  echo "ERROR: EXPERIMENT env var is required." >&2
  echo "  Example: EXPERIMENT=M3_base_lasso_primary bash ${0}" >&2
  exit 1
fi

IGLOO_HEAP="${HEAP_IGLOO_ROOT:-/n/groups/patel/IGLOO}/UKB/HEAP"
MANIFEST="${MANIFEST:-${IGLOO_HEAP}/manifests/module3/${EXPERIMENT}.tsv}"

if [[ ! -f "${MANIFEST}" ]]; then
  echo "ERROR: Manifest not found: ${MANIFEST}" >&2
  echo "Generate it first:" >&2
  echo "  Rscript -e \"" >&2
  echo "    source('workflow/00_paths.R')" >&2
  echo "    source('workflow/config_helpers.R')" >&2
  echo "    source('workflow/generate_manifests.R')" >&2
  echo "    generate_module3_manifest('${EXPERIMENT}')" >&2
  echo "  \"" >&2
  exit 1
fi

N_ROWS=$(tail -n +2 "${MANIFEST}" | wc -l)
if [[ "${N_ROWS}" -eq 0 ]]; then
  echo "ERROR: Manifest has no data rows: ${MANIFEST}" >&2
  exit 1
fi

# Time limit scales with mediation_mode (column 6 of the manifest).
# partitioned_categories is the heaviest (cis + trans + every PXS_<category>
# predictor x all diseases per task) and overran the old flat 2:30 limit.
MODE=$(awk -F'\t' 'NR==2{print $6}' "${MANIFEST}")
case "${MODE}" in
  partitioned_categories)         TIME="0-06:00" ;;
  partitioned_grouped_categories) TIME="0-04:00" ;;
  *)                              TIME="0-02:30" ;;   # primary_total etc.
esac

# Array spec: ARRAY_TASKS overrides with an explicit list (e.g. "68,80,85" to
# rerun specific timed-out tasks); otherwise RESUME_FROM..N_ROWS (full range).
RESUME_FROM="${RESUME_FROM:-1}"
ARRAY_SPEC="${ARRAY_TASKS:-${RESUME_FROM}-${N_ROWS}}"

RSCRIPT="${HEAP_ROOT}/scripts/module3_mediation/Module3.R"
LOG_DIR="${IGLOO_HEAP}/logs/module3"
mkdir -p "${LOG_DIR}"

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Submitting Module 3 manifest job"
echo "  Experiment : ${EXPERIMENT}"
echo "  Manifest   : ${MANIFEST}"
echo "  N rows     : ${N_ROWS}"
echo "  Mode       : ${MODE}"
echo "  Array spec : ${ARRAY_SPEC}"
echo "  Resources  : time=${TIME}, mem=40G, partition=short"

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
    -t "${TIME}" \
    --ntasks=1 \
    --cpus-per-task=1 \
    --mem=40G \
    -p short \
    --array="${ARRAY_SPEC}" \
    --export=ALL \
    "${dep_arg[@]}" \
    -J "M3_${EXPERIMENT}" \
    -o "${LOG_DIR}/${EXPERIMENT}_%A_%a.out" \
    -e "${LOG_DIR}/${EXPERIMENT}_%A_%a.err" \
    --mail-type=FAIL \
    "$0"
fi

# -----------------------------------------------------------------------
# Job body — runs inside Slurm
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
