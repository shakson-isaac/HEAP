#!/bin/bash
umask 002  # group-writable outputs for hpc_patel multi-user runs
# Manifest-driven Module 1 array job.
#
# USAGE
# -----
#   # Step 1: Generate the manifest (run once per experiment)
#   Rscript -e "
#     source('workflow/00_paths.R')
#     source('workflow/config_helpers.R')
#     source('workflow/generate_manifests.R')
#     generate_module1_manifest('M1_base_lasso')
#   "
#
#   # Step 2: Submit the array job
#   EXPERIMENT=M1_base_lasso bash slurm/module1/HEAPmodule1_manifest.sh
#
# REQUIRED ENV VARS
# -----------------
#   EXPERIMENT    Name of the experiment (must match a key in
#                 config/modules/module1_experiments.yml and have a
#                 generated manifest at
#                 /n/groups/patel/IGLOO/UKB/HEAP/manifests/module1/<EXPERIMENT>.tsv)
#
# OPTIONAL ENV VARS
# -----------------
#   HEAP_ROOT     (default: /n/groups/patel/shakson_ukb/HEAP)
#   MANIFEST      Override manifest path
#   DEPENDENCY    sbatch --dependency spec to gate this run, e.g. afterok:<jobid>
#                 (SBATCH_DEPENDENCY is ignored on this cluster, so use this)
#
# RESOURCE SCALING
# ----------------
#   Family is read from the manifest row, so resource allocation is
#   determined here from the FAMILY env var derived from the manifest,
#   or by reading a sample row before submission.
#   Safe default: 8h / 30G / short handles lasso, ridge, enet.
#   For rf: increase to 24h / 50G / medium and set RF_NUM_THREADS=4.

set -euo pipefail

HEAP_ROOT="${HEAP_ROOT:-/n/groups/patel/shakson_ukb/HEAP}"
export HEAP_ROOT
export HEAP_PATHS_FILE="${HEAP_ROOT}/workflow/00_paths.R"

EXPERIMENT="${EXPERIMENT:-}"
if [[ -z "${EXPERIMENT}" ]]; then
  echo "ERROR: EXPERIMENT env var is required." >&2
  echo "  Example: EXPERIMENT=M1_base_lasso bash ${0}" >&2
  exit 1
fi

# Default manifest path (override via MANIFEST env var)
IGLOO_HEAP="${HEAP_IGLOO_ROOT:-/n/groups/patel/IGLOO}/UKB/HEAP"
MANIFEST="${MANIFEST:-${IGLOO_HEAP}/manifests/module1/${EXPERIMENT}.tsv}"

if [[ ! -f "${MANIFEST}" ]]; then
  echo "ERROR: Manifest not found: ${MANIFEST}" >&2
  echo "Generate it first:" >&2
  echo "  Rscript -e \"" >&2
  echo "    source('workflow/00_paths.R')" >&2
  echo "    source('workflow/config_helpers.R')" >&2
  echo "    source('workflow/generate_manifests.R')" >&2
  echo "    generate_module1_manifest('${EXPERIMENT}')" >&2
  echo "  \"" >&2
  exit 1
fi

N_ROWS=$(tail -n +2 "${MANIFEST}" | wc -l)
if [[ "${N_ROWS}" -eq 0 ]]; then
  echo "ERROR: Manifest has no data rows: ${MANIFEST}" >&2
  exit 1
fi

# Determine resource class from the experiment family (read from manifest header)
FAMILY=$(awk -F'\t' 'NR==2{print $5}' "${MANIFEST}")
case "${FAMILY}" in
  lasso|ridge|enet)
    TIME="0-05:00"; MEM="20G"; PARTITION="short"; THREADS=1
    ;;
  rf)
    TIME="1-00:00"; MEM="50G"; PARTITION="medium"; THREADS=4
    ;;
  *)
    # Unknown family: use a safe default
    TIME="0-05:00"; MEM="20G"; PARTITION="short"; THREADS=1
    echo "WARNING: Unknown family '${FAMILY}'; using default resources." >&2
    ;;
esac

# Resource overrides for reruns: e.g. enet chunks that TIMEOUT on the 5h short
# partition need a longer wall on medium. Override any of these via env var:
#   TIME=0-12:00 PARTITION=medium MEM=32G THREADS=1 ARRAY_TASKS=31,65,121 ...
TIME="${TIME_OVERRIDE:-${TIME}}"
MEM="${MEM_OVERRIDE:-${MEM}}"
PARTITION="${PARTITION_OVERRIDE:-${PARTITION}}"
THREADS="${THREADS_OVERRIDE:-${THREADS}}"

# Array spec: ARRAY_TASKS overrides with an explicit list/range (e.g. failed
# chunk indices to rerun: "31,65,121-130"); otherwise full 1-N_ROWS.
ARRAY_SPEC="${ARRAY_TASKS:-1-${N_ROWS}}"

RSCRIPT="${HEAP_ROOT}/scripts/module1_variance_decomposition/Module1_suggested.R"
LOG_DIR="${IGLOO_HEAP}/logs/module1"
mkdir -p "${LOG_DIR}"

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Submitting Module 1 manifest job"
echo "  Experiment : ${EXPERIMENT}"
echo "  Manifest   : ${MANIFEST}"
echo "  N rows     : ${N_ROWS}"
echo "  Family     : ${FAMILY}"
echo "  Resources  : time=${TIME}, mem=${MEM}, partition=${PARTITION}, threads=${THREADS}"
echo "  Array spec : ${ARRAY_SPEC}"

# If not already inside a Slurm job, submit self
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
    --cpus-per-task="${THREADS}" \
    --mem="${MEM}" \
    -p "${PARTITION}" \
    --array="${ARRAY_SPEC}" \
    --export=ALL \
    "${dep_arg[@]}" \
    -J "M1_${EXPERIMENT}" \
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

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

IDX="${SLURM_ARRAY_TASK_ID}"

echo "[$(date '+%H:%M:%S')] START  experiment=${EXPERIMENT}  array_index=${IDX}"

# Pass manifest path and array index; the R script reads its own row.
Rscript "${RSCRIPT}" \
  --manifest "${MANIFEST}" \
  --array-index "${IDX}"

STATUS=$?
echo "[$(date '+%H:%M:%S')] FINISH  experiment=${EXPERIMENT}  array_index=${IDX}  exit=${STATUS}"
exit "${STATUS}"
