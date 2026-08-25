#!/bin/bash
###############################################################################
# submit_gwas_exposures.sh
#
# One-command launcher for the HEAP exposure GWAS (REGENIE).
# Run this from an O2 login node:
#
#     /n/groups/patel/shakson_ukb/HEAP/slurm/gwas_regenie/submit_gwas_exposures.sh
#
# It does, in order:
#   1. Runs prepare_gwas_exposures.R (batch mode) to (re)generate the per-type
#      exposure lists (evars_continuous_heap.txt / evars_binary_heap.txt) and the
#      QC tables. Dispatched with `srun` (not the login node) only because it loads
#      HEAP.rds and stacks all exposures. Batch prep no longer stages per-exposure
#      pheno/covar files -- each array job regenerates its own in single-exposure
#      mode -- so this step is now lists + QC only: a single-threaded, ~1-2 minute,
#      ~8 GB job (HEAP.rds is 1.3 GB on disk, ~8 GB peak RSS when loaded + joined).
#      Defaults below are sized to that (16 GB / 1 core / 30 min) with headroom; the
#      memory-heavy work (regenie step 1) lives in the array jobs, not here.
#   2. Reads the line counts of the two generated lists -- this is how N_CONT and
#      N_BIN are determined; no one has to count by hand.
#   3. Submits one REGENIE array job per exposure for each type, sized to match.
#
# Options:
#   --skip-prep   Lists already exist and are current; skip step 1, just submit.
#   --dry-run     Print every srun/sbatch command but execute nothing.
#   -h, --help    Show this header.
#
# Resource overrides for the prep step (env vars, all optional):
#   PREP_PART (short)  PREP_TIME (0-00:30)  PREP_MEM (16G)  PREP_CPUS (1)
#
# Override the HEAP checkout with HEAP_ROOT=/path/to/HEAP if not the default.
###############################################################################

set -euo pipefail
umask 0002   # group-writable outputs for hpc_patel team runs (files 664, dirs 775)

# --- Resolve HEAP root (from this script's location, overridable) ------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
HEAP_ROOT="${HEAP_ROOT:-$(cd "${SCRIPT_DIR}/../.." && pwd)}"
export HEAP_ROOT
export HEAP_PATHS_FILE="${HEAP_ROOT}/workflow/00_paths.R"

GWAS_DIR="${HEAP_ROOT}/slurm/gwas_regenie"
PREP_R="${HEAP_ROOT}/scripts/gwas_regenie/prepare_gwas_exposures.R"
EVARS_CONT="${GWAS_DIR}/evars_continuous_heap.txt"
EVARS_BIN="${GWAS_DIR}/evars_binary_heap.txt"
ARRAY_CONT="${GWAS_DIR}/gwas_regenie_exposures_continuous_v2.sh"
ARRAY_BIN="${GWAS_DIR}/gwas_regenie_exposures_binary_v2.sh"

# --- Prep step resources (overridable via env) -------------------------------
# Right-sized from real runs: HEAP.rds is 1.3 GB on disk and peaks at ~8 GB RSS once
# loaded + joined; the batch loop is single-threaded and finishes in ~1-2 min. 16 GB /
# 1 core / 30 min gives ~2x headroom. (Was 64 GB / 4 / 3h, based on a wrong ~50 GB
# HEAP.rds assumption -- that figure actually belongs to regenie step 1 in the arrays.)
PREP_PART="${PREP_PART:-short}"
PREP_TIME="${PREP_TIME:-0-00:30}"
PREP_MEM="${PREP_MEM:-16G}"
PREP_CPUS="${PREP_CPUS:-1}"

# --- Parse options -----------------------------------------------------------
SKIP_PREP=0
DRY_RUN=0
usage() { sed -n '3,27p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; }
for arg in "$@"; do
  case "$arg" in
    --skip-prep) SKIP_PREP=1 ;;
    --dry-run)   DRY_RUN=1 ;;
    -h|--help)   usage; exit 0 ;;
    *) echo "ERROR: unknown argument '$arg' (try --help)" >&2; exit 2 ;;
  esac
done

run() { echo "+ $*"; [[ "$DRY_RUN" == 1 ]] || "$@"; }

# --- Sanity: required files present ------------------------------------------
for f in "$PREP_R" "$ARRAY_CONT" "$ARRAY_BIN"; do
  [[ -f "$f" ]] || { echo "ERROR: missing required file: $f" >&2; exit 1; }
done

# --- Step 1: prep (generates the evar lists + QC tables) ---------------------
if [[ "$SKIP_PREP" == 0 ]]; then
  echo "[1/3] prep: prepare_gwas_exposures.R via srun (mem=${PREP_MEM}, part=${PREP_PART})"
  run srun --partition="$PREP_PART" --time="$PREP_TIME" --mem="$PREP_MEM" \
      --cpus-per-task="$PREP_CPUS" --job-name=heap_gwas_prep \
      bash -lc "umask 0002 && module load gcc/14.2.0 R/4.4.2 && Rscript '$PREP_R'"
else
  echo "[1/3] prep: --skip-prep set; reusing existing evar lists"
fi

# --- Step 2: determine N_CONT / N_BIN from the generated lists ---------------
# awk NR matches how the array scripts index the file (sed -n "Np"), and counts a
# final line even if it lacks a trailing newline.
count_lines() { [[ -f "$1" ]] && awk 'END{print NR+0}' "$1" || echo 0; }
N_CONT=$(count_lines "$EVARS_CONT")
N_BIN=$(count_lines "$EVARS_BIN")
echo "[2/3] list sizes: continuous=${N_CONT}  binary=${N_BIN}"

if [[ "$N_CONT" -eq 0 && "$N_BIN" -eq 0 ]]; then
  echo "ERROR: both exposure lists are empty/missing. Did prep succeed?" >&2
  echo "       Expected: ${EVARS_CONT}" >&2
  echo "                 ${EVARS_BIN}" >&2
  exit 1
fi

# --- Step 3: submit one array job per type, sized to its list ----------------
echo "[3/3] submitting REGENIE arrays"
if [[ "$N_CONT" -gt 0 ]]; then
  run sbatch --array="1-${N_CONT}" "$ARRAY_CONT"
else
  echo "  (no continuous exposures -- skipping continuous array)"
fi
if [[ "$N_BIN" -gt 0 ]]; then
  run sbatch --array="1-${N_BIN}" "$ARRAY_BIN"
else
  echo "  (no binary exposures -- skipping binary array)"
fi

echo "[done] track progress with:  squeue -u \$USER"
