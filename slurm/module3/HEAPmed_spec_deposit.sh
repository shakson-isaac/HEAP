#!/bin/bash
umask 002  # group-writable outputs for hpc_patel multi-user runs
# Per-specification mediation deposit: aggregate the 8 Module-3 covariate
# specifications into mediation/exposome/<spec>.tsv + mediation/genetic/<spec>.tsv
# for HEAP_Supplementary_Data.zip.
#
# Two stages, because one spec is ~800 MDres files and several GB:
#   stage 1  array 1-8, one spec each -> intermediate/med_spec_cache/<exp>.rds
#   stage 2  single job gated on stage 1 -> writes the TSVs, verifies that the
#            regenerated BASE reproduces the shipped med_exposure_total.tsv, and
#            refuses to ship if it does not (rule 5: that is an author decision).
#
# USAGE
# -----
#   bash slurm/module3/HEAPmed_spec_deposit.sh            # submits both stages
#   STAGE=finalize bash slurm/module3/HEAPmed_spec_deposit.sh   # re-run stage 2 only
#
# OPTIONAL ENV VARS
# -----------------
#   HEAP_ROOT     (default: /n/groups/patel/shakson_ukb/HEAP)
#   ARRAY_TASKS   explicit list to rerun specific specs, e.g. "3,7"
#   DEPENDENCY    sbatch --dependency spec to gate stage 1
#                 (SBATCH_DEPENDENCY is ignored on this cluster, so use this)
set -euo pipefail

HEAP_ROOT="${HEAP_ROOT:-/n/groups/patel/shakson_ukb/HEAP}"
export HEAP_ROOT
export HEAP_PATHS_FILE="${HEAP_ROOT}/workflow/00_paths.R"
IGLOO_HEAP="${HEAP_IGLOO_ROOT:-/n/groups/patel/IGLOO}/UKB/HEAP"
RSCRIPT="${HEAP_ROOT}/scripts/analysis_summaries/export_med_spec_deposit.R"
LOG_DIR="${IGLOO_HEAP}/logs/module3"      # IGLOO, not the source tree (not group-writable)
mkdir -p "${LOG_DIR}"

# the spec grid, in the same order as SPECS in the R script
SPEC_LIST=(
  M3_base_lasso_primary
  M3_base_bmi_lasso_primary
  M3_base_clinical_lasso_primary
  M3_base_prevalent_lasso_primary
  M3_base_exclprev_lasso_primary
  M3_base_draw_lasso_primary
  M3_base_ridge_primary
  M3_base_enet_primary
)
N=${#SPEC_LIST[@]}

# ---------------------------------------------------------------------------
# Submitter. SLURM_ARRAY_TASK_ID is set ONLY inside an array element; guarding on
# SLURM_JOB_ID would misfire when launching from an OnDemand/VS Code allocation.
# ---------------------------------------------------------------------------
if [[ -z "${SLURM_ARRAY_TASK_ID:-}" && "${STAGE:-}" != "finalize_body" ]]; then
  if [[ "${STAGE:-}" == "finalize" ]]; then
    exec sbatch -t 0-01:00 --ntasks=1 --cpus-per-task=1 --mem=64G -p short \
      --export=ALL,STAGE=finalize_body \
      -J "M3_medspec_finalize" \
      -o "${LOG_DIR}/med_spec_finalize_%j.out" -e "${LOG_DIR}/med_spec_finalize_%j.err" \
 "$0"
  fi

  dep_arg=(); [[ -n "${DEPENDENCY:-}" ]] && dep_arg=(--dependency="${DEPENDENCY}")
  ARRAY_SPEC="${ARRAY_TASKS:-1-${N}}"    # no %N throttle, by house rule
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] Submitting per-specification mediation deposit"
  echo "  Specs      : ${N} (${ARRAY_SPEC})"
  echo "  Resources  : stage1 time=0-03:00 mem=48G | stage2 time=0-01:00 mem=64G"

  JID=$(sbatch --parsable \
    -t 0-03:00 --ntasks=1 --cpus-per-task=1 --mem=48G -p short \
    --array="${ARRAY_SPEC}" --export=ALL "${dep_arg[@]}" \
    -J "M3_medspec" \
    -o "${LOG_DIR}/med_spec_%A_%a.out" -e "${LOG_DIR}/med_spec_%A_%a.err" \
 "$0")
  echo "  stage 1 job: ${JID}"

  FID=$(sbatch --parsable \
    -t 0-01:00 --ntasks=1 --cpus-per-task=1 --mem=64G -p short \
    --dependency="afterok:${JID}" --export=ALL,STAGE=finalize_body \
    -J "M3_medspec_finalize" \
    -o "${LOG_DIR}/med_spec_finalize_%j.out" -e "${LOG_DIR}/med_spec_finalize_%j.err" \
 "$0")
  echo "  stage 2 job: ${FID}  (afterok:${JID})"
  echo
  echo "  watch:  squeue -u \$USER -j ${JID},${FID}"
  echo "  logs :  ${LOG_DIR}/med_spec_*"
  exit 0
fi

# ---------------------------------------------------------------------------
# Job body
# ---------------------------------------------------------------------------
module load gcc/14.2.0
module load R/4.4.2
ulimit -n 10000

if [[ "${STAGE:-}" == "finalize_body" ]]; then
  echo "[$(date '+%H:%M:%S')] START  finalize"
  Rscript "${RSCRIPT}" --finalize
  STATUS=$?
  echo "[$(date '+%H:%M:%S')] FINISH finalize  exit=${STATUS}"
  exit "${STATUS}"
fi

IDX="${SLURM_ARRAY_TASK_ID}"
EXPERIMENT="${SPEC_LIST[$((IDX - 1))]}"
echo "[$(date '+%H:%M:%S')] START  spec=${EXPERIMENT}  array_index=${IDX}"
Rscript "${RSCRIPT}" --spec "${EXPERIMENT}"
STATUS=$?
echo "[$(date '+%H:%M:%S')] FINISH spec=${EXPERIMENT}  exit=${STATUS}"
exit "${STATUS}"
