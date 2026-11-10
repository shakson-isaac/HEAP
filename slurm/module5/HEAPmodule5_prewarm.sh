#!/bin/bash
umask 002  # group-writable outputs
# Pre-warm the Module 5 instrument caches: clump every UNIQUE protein / exposure /
# disease ONCE before the main MR array, so the 2,222-task array starts warm and
# doesn't redundantly clump the same instruments concurrently (~5-8x saving on the
# cold run). Run this FIRST; then submit HEAPmodule5_manifest.sh.
#
# USAGE
#   bash slurm/module5/HEAPmodule5_prewarm.sh
#   # gate the main run on it:  DEPENDENCY=afterok:<prewarm_jobid> EXPERIMENT=... bash HEAPmodule5_manifest.sh
#
# Submits one small array per kind (caches are shared: exposure/disease warmed via
# the UKB runner are reused by the deCODE arm; only proteins are arm-specific).
#
# OPTIONAL ENV: HEAP_ROOT, TIME, MEM, CPUS, PARTITION, DEPENDENCY,
#               NSLICES_ukb_protein / _decode_protein / _exposure / _disease

set -euo pipefail

HEAP_ROOT="${HEAP_ROOT:-/n/groups/patel/shakson_ukb/HEAP}"
export HEAP_ROOT
export HEAP_PATHS_FILE="${HEAP_ROOT}/workflow/00_paths.R"
IGLOO_HEAP="${HEAP_IGLOO_ROOT:-/n/groups/patel/IGLOO}/UKB/HEAP"
LOG_DIR="${IGLOO_HEAP}/logs/module5"
mkdir -p "${LOG_DIR}"

TIME="${TIME:-0-03:00}"
MEM="${MEM:-24G}"
CPUS="${CPUS:-2}"
PARTITION="${PARTITION:-short}"

# Parallelism per kind (tasks). More slices = shorter wall (no array throttle).
# Defaults sized so each task handles ~13 proteins (~108 s each → ~25 min/task).
# ~1067 unique proteins (PD ∪ PE), ~120 exposures, 71 diseases; ~54 s per clump.
NSLICES_ukb_protein="${NSLICES_ukb_protein:-80}"
NSLICES_decode_protein="${NSLICES_decode_protein:-80}"
NSLICES_exposure="${NSLICES_exposure:-12}"
NSLICES_disease="${NSLICES_disease:-8}"

# -----------------------------------------------------------------------
# Worker branch: inside an array task (KIND + NSLICES exported by submitter)
# -----------------------------------------------------------------------
if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
  module load gcc/14.2.0
  module load R/4.4.2
  ulimit -n 10000 2>/dev/null || true  # best-effort; some nodes cap the hard limit
  export OMP_NUM_THREADS="${CPUS}" OPENBLAS_NUM_THREADS="${CPUS}"
  echo "[$(date '+%H:%M:%S')] prewarm kind=${KIND} slice=${SLURM_ARRAY_TASK_ID}/${NSLICES}"
  Rscript "${HEAP_ROOT}/scripts/module5_mr/prewarm_instruments.R" \
    "${KIND}" "${SLURM_ARRAY_TASK_ID}" "${NSLICES}"
  exit $?
fi

# -----------------------------------------------------------------------
# Submitter branch: one array job per kind
# -----------------------------------------------------------------------
dep_arg=()
[[ -n "${DEPENDENCY:-}" ]] && dep_arg=(--dependency="${DEPENDENCY}")

submit_kind() {
  local k="$1" ns="$2"
  echo "Submitting prewarm: kind=${k}  array=1-${ns}  (t=${TIME}, mem=${MEM}, p=${PARTITION})"
  sbatch \
    -t "${TIME}" --ntasks=1 --cpus-per-task="${CPUS}" --mem="${MEM}" -p "${PARTITION}" \
    --array="1-${ns}" \
    --export=ALL,KIND="${k}",NSLICES="${ns}" \
    "${dep_arg[@]}" \
    -J "M5warm_${k}" \
    -o "${LOG_DIR}/prewarm_${k}_%A_%a.out" \
    -e "${LOG_DIR}/prewarm_${k}_%A_%a.err" \
    "$0"
}

# KINDS env selects which kinds to submit (default all). Idempotent — already-cached
# instruments are skipped — so a partial re-run (e.g. KINDS=decode_protein) just fills gaps.
KINDS="${KINDS:-ukb_protein decode_protein exposure disease}"
for k in ${KINDS}; do
  case "${k}" in
    ukb_protein)    submit_kind ukb_protein    "${NSLICES_ukb_protein}" ;;
    decode_protein) submit_kind decode_protein "${NSLICES_decode_protein}" ;;
    exposure)       submit_kind exposure       "${NSLICES_exposure}" ;;
    disease)        submit_kind disease        "${NSLICES_disease}" ;;
    *) echo "WARNING: unknown kind '${k}' — skipping" >&2 ;;
  esac
done

echo "Pre-warm arrays submitted (kinds: ${KINDS}). Gate the main run with DEPENDENCY=afterok:<jobid> if desired."
