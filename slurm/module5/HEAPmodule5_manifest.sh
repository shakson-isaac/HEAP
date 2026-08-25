#!/bin/bash
umask 002  # group-writable outputs for hpc_patel multi-user runs
# Manifest-driven Module 5 (Mendelian Randomization) array job.
#
# Runs the SAME triad/edge lists (global_edges/edges_<TYPE>.tsv) through either
# MR arm, selected by EXPERIMENT:
#   MR_UKB_primary        -> split-sample UKB  (runner Module5.R)
#   MR_deCODE_replication -> deCODE SomaScan   (runner Module5_deCODE.R)
# The runner for each row is read from the manifest `runner` column, so one
# launcher serves both arms.
#
# USAGE
# -----
#   # Step 0 (once): build the triad/edge lists
#   Rscript scripts/module5_mr/Module5_load.R
#
#   # Step 1 (once per experiment): generate the manifest
#   Rscript -e "
#     source('workflow/00_paths.R')
#     source('workflow/config_helpers.R')
#     source('workflow/generate_manifests.R')
#     generate_module5_manifest('MR_UKB_primary')
#   "
#
#   # Step 2: submit the array job
#   EXPERIMENT=MR_UKB_primary        bash slurm/module5/HEAPmodule5_manifest.sh
#   EXPERIMENT=MR_deCODE_replication bash slurm/module5/HEAPmodule5_manifest.sh
#
# REQUIRED ENV VARS
# -----------------
#   EXPERIMENT    Key in config/modules/module5_experiments.yml with a manifest
#                 at .../manifests/module5/<EXPERIMENT>.tsv
#
# OPTIONAL ENV VARS
# -----------------
#   HEAP_ROOT     (default: /n/groups/patel/shakson_ukb/HEAP)
#   MANIFEST      Override manifest path
#   TIME/MEM/CPUS/PARTITION  Override resource request
#   MAKE_PLOTS    1 to emit per-edge scatter plots (default 0; viz layer handles plots)
#   DEPENDENCY    sbatch --dependency spec to gate this run (SBATCH_DEPENDENCY is
#                 ignored on O2, so pass it here), e.g. afterok:<jobid>

set -euo pipefail

HEAP_ROOT="${HEAP_ROOT:-/n/groups/patel/shakson_ukb/HEAP}"
export HEAP_ROOT
export HEAP_PATHS_FILE="${HEAP_ROOT}/workflow/00_paths.R"

EXPERIMENT="${EXPERIMENT:-}"
if [[ -z "${EXPERIMENT}" ]]; then
  echo "ERROR: EXPERIMENT env var is required." >&2
  echo "  Example: EXPERIMENT=MR_UKB_primary bash ${0}" >&2
  exit 1
fi

IGLOO_HEAP="${HEAP_IGLOO_ROOT:-/n/groups/patel/IGLOO}/UKB/HEAP"
MANIFEST="${MANIFEST:-${IGLOO_HEAP}/manifests/module5/${EXPERIMENT}.tsv}"

if [[ ! -f "${MANIFEST}" ]]; then
  echo "ERROR: Manifest not found: ${MANIFEST}" >&2
  echo "Generate it first:" >&2
  echo "  Rscript -e \"source('workflow/00_paths.R'); source('workflow/config_helpers.R'); source('workflow/generate_manifests.R'); generate_module5_manifest('${EXPERIMENT}')\"" >&2
  exit 1
fi

N_ROWS=$(tail -n +2 "${MANIFEST}" | wc -l)
if [[ "${N_ROWS}" -eq 0 ]]; then
  echo "ERROR: Manifest has no data rows: ${MANIFEST}" >&2
  exit 1
fi

# Resources (override via env). MR = clumping (plink on EUR LD ref) + loading
# disease/protein/exposure GWAS. Instruments are disk-cached, so the first full
# run is the expensive one; reruns reuse caches. Sized for the heavy PD/DP rows.
TIME="${TIME:-0-08:00}"
MEM="${MEM:-32G}"
CPUS="${CPUS:-4}"
PARTITION="${PARTITION:-short}"
MAKE_PLOTS="${MAKE_PLOTS:-0}"

# Array spec: ARRAY_TASKS overrides with an explicit list/range (e.g. failed
# task indices to rerun: "1002,1881-1905,..."); otherwise the full 1-N_ROWS range.
ARRAY_SPEC="${ARRAY_TASKS:-1-${N_ROWS}}"

LOG_DIR="${IGLOO_HEAP}/logs/module5"
mkdir -p "${LOG_DIR}"

# -----------------------------------------------------------------------
# Submitter branch: submit self as an array job, then exit.
# (Guard on SLURM_ARRAY_TASK_ID, NOT SLURM_JOB_ID — the latter is also set
#  inside interactive/VS Code allocations.)
# -----------------------------------------------------------------------
if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] Submitting Module 5 MR manifest job"
  echo "  Experiment : ${EXPERIMENT}"
  echo "  Manifest   : ${MANIFEST}"
  echo "  N rows     : ${N_ROWS}"
  echo "  Array spec : ${ARRAY_SPEC}"
  echo "  Resources  : time=${TIME}, mem=${MEM}, partition=${PARTITION}, cpus=${CPUS}"

  dep_arg=()
  [[ -n "${DEPENDENCY:-}" ]] && dep_arg=(--dependency="${DEPENDENCY}")
  exec sbatch \
    -t "${TIME}" \
    --ntasks=1 \
    --cpus-per-task="${CPUS}" \
    --mem="${MEM}" \
    -p "${PARTITION}" \
    --array="${ARRAY_SPEC}" \
    --export=ALL \
    "${dep_arg[@]}" \
    -J "M5_${EXPERIMENT}" \
    -o "${LOG_DIR}/${EXPERIMENT}_%A_%a.out" \
    -e "${LOG_DIR}/${EXPERIMENT}_%A_%a.err" \
    --mail-type=FAIL \
    "$0"
fi

# -----------------------------------------------------------------------
# Job body — runs inside the Slurm array element.
# -----------------------------------------------------------------------
module load gcc/14.2.0
module load R/4.4.2
ulimit -n 10000 2>/dev/null || true  # best-effort; some nodes cap the hard limit

export OMP_NUM_THREADS="${CPUS}"
export OPENBLAS_NUM_THREADS="${CPUS}"

IDX="${SLURM_ARRAY_TASK_ID}"

# Resolve this task's manifest row by the array_index column (header-driven, so
# robust to column reordering).
get_col() {
  awk -F'\t' -v idx="${IDX}" -v want="$1" '
    NR==1 { for (i=1;i<=NF;i++) h[$i]=i; next }
    $(h["array_index"]) == idx { print $(h[want]) }
  ' "${MANIFEST}"
}

RUNNER="$(get_col runner)"
EDGE_TYPE="$(get_col edge_type)"
CHUNK_ID="$(get_col chunk_id)"
N_CHUNKS="$(get_col n_chunks)"
OUTBASE="$(get_col output_path)"

if [[ -z "${RUNNER}" || -z "${EDGE_TYPE}" || -z "${CHUNK_ID}" || -z "${N_CHUNKS}" ]]; then
  echo "ERROR: could not resolve manifest row for array_index=${IDX} in ${MANIFEST}" >&2
  exit 1
fi

# Drive the runner's edge-output base from the manifest's (IGLOO-canonical,
# experiment-keyed) output_path, so config = where outputs actually land.
[[ -n "${OUTBASE}" ]] && export HEAP_MR_OUTDIR="${OUTBASE}"

RSCRIPT="${HEAP_ROOT}/scripts/module5_mr/${RUNNER}"
if [[ ! -f "${RSCRIPT}" ]]; then
  echo "ERROR: runner script not found: ${RSCRIPT}" >&2
  exit 1
fi

echo "[$(date '+%H:%M:%S')] START  ${EXPERIMENT}  idx=${IDX}  runner=${RUNNER}  edge=${EDGE_TYPE}  chunk=${CHUNK_ID}/${N_CHUNKS}"

# Module5.R / Module5_deCODE.R take positional args: <edge_type> <chunk_idx> <n_chunks> [make_plots]
Rscript "${RSCRIPT}" "${EDGE_TYPE}" "${CHUNK_ID}" "${N_CHUNKS}" "${MAKE_PLOTS}"

STATUS=$?
echo "[$(date '+%H:%M:%S')] FINISH ${EXPERIMENT}  idx=${IDX}  edge=${EDGE_TYPE}  chunk=${CHUNK_ID}  exit=${STATUS}"
exit "${STATUS}"
