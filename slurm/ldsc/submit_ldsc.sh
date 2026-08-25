#!/bin/bash
###############################################################################
# submit_ldsc.sh — one-command launcher for the HEAP LDSC heritability stage.
#
# Run from an O2 login node:
#     /n/groups/patel/shakson_ukb/HEAP/slurm/ldsc/submit_ldsc.sh
#
# It does, in order:
#   1. (Re)generates slurm/ldsc/ldsc_exposures.txt = every exposure that has a
#      COMPLETED GWAS (a <exp>.regenie under output/gwas/regenie_step2/<exp>/).
#      LDSC can only run on finished GWAS, so the list is derived from them.
#   2. Reads the line count to size the array (no manual N).
#   3. sbatch-submits ldsc_h2_array.sh as an array (one exposure per task).
#
# Options:
#   --rg          After the h2 array, also submit ldsc_rg.sh (genetic
#                 correlations) with an afterok dependency on the h2 array.
#   --dry-run     Print the sbatch/find commands but execute nothing.
#   -h, --help    Show this header.
#
# Override the HEAP checkout with HEAP_ROOT=/path/to/HEAP and the exposure GWAS
# root with HEAP_REGENIE_STEP2_DIR if they are not the defaults.
###############################################################################
set -euo pipefail
umask 0002

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
HEAP_ROOT="${HEAP_ROOT:-$(cd "${SCRIPT_DIR}/../.." && pwd)}"
IGLOO="${HEAP_IGLOO_ROOT:-/n/groups/patel/IGLOO}"
STEP2_DIR="${HEAP_REGENIE_STEP2_DIR:-${IGLOO}/UKB/HEAP/output/gwas/regenie_step2}"

EXP_LIST="${SCRIPT_DIR}/ldsc_exposures.txt"
ARRAY="${SCRIPT_DIR}/ldsc_h2_array.sh"
RG="${SCRIPT_DIR}/ldsc_rg.sh"

DO_RG=0; DRY_RUN=0
usage() { sed -n '3,24p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; }
for arg in "$@"; do
  case "$arg" in
    --rg)      DO_RG=1 ;;
    --dry-run) DRY_RUN=1 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown argument '$arg' (try --help)" >&2; exit 2 ;;
  esac
done
run() { echo "+ $*"; [[ "$DRY_RUN" == 1 ]] || "$@"; }

[[ -f "$ARRAY" ]] || { echo "ERROR: missing $ARRAY" >&2; exit 1; }
[[ -d "$STEP2_DIR" ]] || { echo "ERROR: exposure GWAS dir not found: $STEP2_DIR" >&2; exit 1; }

# --- Step 1: list exposures with a completed GWAS ----------------------------
# The scan is read-only and just (re)writes a tracked text file, so it always
# runs -- even under --dry-run -- so the reported count reflects reality.
echo "[1/3] scanning completed exposure GWAS under ${STEP2_DIR}"
: > "${EXP_LIST}"
for d in "${STEP2_DIR}"/*/; do
  e=$(basename "${d}")
  if [[ -f "${d}/${e}.regenie" || -f "${d}/regenie_step2_${e}_${e}.regenie" ]]; then
    echo "${e}" >> "${EXP_LIST}"
  fi
done
sort -o "${EXP_LIST}" "${EXP_LIST}"

count_lines() { [[ -f "$1" ]] && awk 'END{print NR+0}' "$1" || echo 0; }
N=$(count_lines "${EXP_LIST}")
echo "[2/3] ${N} completed exposures -> ${EXP_LIST}"
if [[ "${N}" -eq 0 ]]; then
  echo "ERROR: no completed exposure GWAS found; nothing to submit." >&2
  exit 1
fi

# --- Step 3: submit the h2 array (+ optional rg with afterok dependency) ------
echo "[3/3] submitting LDSC h2 array (1-${N})"
if [[ "$DRY_RUN" == 1 ]]; then
  echo "+ sbatch --array=1-${N} ${ARRAY}"
  [[ "$DO_RG" == 1 ]] && echo "+ sbatch --dependency=afterok:<h2_jobid> ${RG}"
else
  H2_SUB=$(sbatch --array="1-${N}" "${ARRAY}")
  echo "${H2_SUB}"
  H2_ID=$(echo "${H2_SUB}" | awk '{print $NF}')
  if [[ "$DO_RG" == 1 ]]; then
    echo "    submitting ldsc_rg.sh with afterok:${H2_ID}"
    sbatch --dependency="afterok:${H2_ID}" "${RG}"
  fi
fi

echo "[done] track with:  squeue -u \$USER"
echo "       collect results after completion:"
echo "         module load gcc/14.2.0 R/4.4.2"
echo "         HEAP_PATHS_FILE=${HEAP_ROOT}/workflow/00_paths.R Rscript ${HEAP_ROOT}/scripts/ldsc/collect_ldsc_h2.R"
