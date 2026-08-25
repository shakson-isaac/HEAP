#!/bin/bash
#SBATCH -c 8
#SBATCH -t 2-00:00
#SBATCH --mem=50G
#SBATCH -p medium
# Array size is set at submission time:
#   sbatch --array=1-<N_BIN> gwas_regenie_exposures_binary_v2.sh
# where N_BIN is the number of lines in evars_binary_heap.txt.
#
# Step 0 (one-time setup before submitting):
#   Rscript scripts/gwas_regenie/prepare_gwas_exposures.R
#   (generates evars_binary_heap.txt and evars_continuous_heap.txt)
#
# Each array job generates its OWN per-exposure pheno/covar files inside a
# job-local temporary directory -- no canonical scratch dependency.
#
# GWAS for binary exposure phenotypes.
# Uses --bt (logistic regression) in both steps and Firth correction in step 2.
# Does NOT apply --apply-rint (rank-inverse normal is for quantitative traits only).
#
# Two-sample MR design:
#   prepare_gwas_exposures.R excludes all proteomics cohort participants.

set -euo pipefail

# Resolve HEAP root
HEAP_ROOT="${HEAP_ROOT:-/n/groups/patel/shakson_ukb/HEAP}"
umask 0002   # group-writable outputs for hpc_patel team runs
export HEAP_ROOT
export HEAP_PATHS_FILE="${HEAP_ROOT}/workflow/00_paths.R"

# Exposure name from HEAP-managed evar list
EVARS_FILE="${HEAP_ROOT}/slurm/gwas_regenie/evars_binary_heap.txt"
if [[ ! -f "${EVARS_FILE}" ]]; then
  echo "ERROR: ${EVARS_FILE} not found." >&2
  echo "Run: Rscript scripts/gwas_regenie/prepare_gwas_exposures.R" >&2
  exit 1
fi
Ename=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${EVARS_FILE}")
if [[ -z "${Ename}" ]]; then
  echo "ERROR: No exposure at array index ${SLURM_ARRAY_TASK_ID} in ${EVARS_FILE}" >&2
  exit 1
fi
echo "[${Ename}] Task ${SLURM_ARRAY_TASK_ID}: binary exposure GWAS"

# Paths
# regenie conda env: HEAP_REGENIE_ENV override -> shared IGLOO copy -> personal home env.
CONDA_ENV="${HEAP_REGENIE_ENV:-${HEAP_IGLOO_ROOT:-/n/groups/patel/IGLOO}/envs/regenie_env}"
[ -d "${CONDA_ENV}" ] || CONDA_ENV="${CONDA_ENVS:-$HOME/.conda/envs}/regenie_env"
GWAS_PATH="/n/groups/patel/IGLOO/UKB/gwas"
SCRATCH_ROOT="${HEAP_SCRATCH_ROOT:-/n/scratch/users/${USER:0:1}/${USER}}"

# Per-job temporary directories (all node-local or per-job; never canonical)
JOB_TMPDIR="${SLURM_TMPDIR:-/tmp}/heap_gwas_${SLURM_JOB_ID:-$$}_${SLURM_ARRAY_TASK_ID:-0}"

# Exposure pheno/covar input: generated fresh per job (ephemeral; not canonical scratch).
# prepare_gwas_exposures.R ALWAYS nests its output under a per-exposure subdir
# (<output-dir>/<exposure>/pheno.txt), so hand it the PARENT dir and let it create the
# ${Ename} subdir itself. TMP_EXPOSURE_INPUT_DIR is that final nested dir (what step2 reads).
TMP_EXPOSURE_INPUT_PARENT="${JOB_TMPDIR}/exposure_input"
TMP_EXPOSURE_INPUT_DIR="${TMP_EXPOSURE_INPUT_PARENT}/${Ename}"

# Regenie step1 null model: node-local (large; ephemeral)
OUTDIR_STEP1="${JOB_TMPDIR}/regenie_step1"
TMPDIR_STEP1="${JOB_TMPDIR}/regenie_step1_lowmem"

# Regenie step2 working directory: scratch (cleaned after IGLOO copy)
OUTDIR_STEP2_SCRATCH="${SCRATCH_ROOT}/regenie_step2/${Ename}"

# Canonical IGLOO final destination for step2 summary stats
IGLOO_HEAP="${HEAP_IGLOO_ROOT:-/n/groups/patel/IGLOO}/UKB/HEAP"
OUTDIR_STEP2_IGLOO="${HEAP_REGENIE_STEP2_DIR:-${IGLOO_HEAP}/output/gwas/regenie_step2}/${Ename}"

# Optional: keep job-local tmp for debugging (default: clean on exit)
KEEP_GWAS_TMP="${KEEP_GWAS_TMP:-0}"

# Start the scratch step2 working dir FRESH so stale files from an older run of
# this exposure can never be picked up or copied to the canonical location.
rm -rf "${OUTDIR_STEP2_SCRATCH}"
mkdir -p "${TMP_EXPOSURE_INPUT_DIR}" "${OUTDIR_STEP1}" "${TMPDIR_STEP1}" \
         "${OUTDIR_STEP2_SCRATCH}" "${OUTDIR_STEP2_IGLOO}"
mkdir -p "${HEAP_IGLOO_ROOT:-/n/groups/patel/IGLOO}/UKB/HEAP/logs/gwas"

cleanup_tmp() {
  if [[ "${KEEP_GWAS_TMP}" != "1" && -d "${JOB_TMPDIR}" ]]; then
    rm -rf "${JOB_TMPDIR}"
  fi
}
trap cleanup_tmp EXIT

# Step 0: Generate pheno/covar files for this exposure using R
# This is a per-job ephemeral step -- output goes to TMP_EXPOSURE_INPUT_DIR only.
echo "[${Ename}] Generating pheno/covar files (single-exposure mode)"
module load gcc/14.2.0
module load R/4.4.2
Rscript "${HEAP_ROOT}/scripts/gwas_regenie/prepare_gwas_exposures.R" \
  --exposure-id "${Ename}" \
  --output-dir  "${TMP_EXPOSURE_INPUT_PARENT}"

if [[ ! -f "${TMP_EXPOSURE_INPUT_DIR}/pheno.txt" ]]; then
  echo "ERROR: pheno.txt not created for ${Ename} in ${TMP_EXPOSURE_INPUT_DIR}" >&2
  exit 1
fi
echo "[${Ename}] Pheno/covar files ready"

module load conda/miniforge3/24.11.3-0
eval "$(conda shell.bash hook)"
conda activate "${CONDA_ENV}"

# Step 1: whole-genome regression (null model, binary)
echo "[${Ename}] Running REGENIE step 1 (binary, --bt)"
regenie \
  --step 1 \
  --pgen "${GWAS_PATH}/ukb_nonimputed_snps" \
  --phenoFile "${TMP_EXPOSURE_INPUT_DIR}/pheno.txt" \
  --covarFile "${TMP_EXPOSURE_INPUT_DIR}/covar.txt" \
  --bt \
  --bsize 1000 \
  --threads 8 \
  --lowmem \
  --lowmem-prefix "${TMPDIR_STEP1}/regenie_tmp_preds" \
  --out "${OUTDIR_STEP1}/regenie_step1_${Ename}"

# Step 2: association testing with Firth logistic regression
echo "[${Ename}] Running REGENIE step 2 (binary, --bt --firth)"
regenie \
  --step 2 \
  --pgen "${GWAS_PATH}/UKBallchr" \
  --phenoFile "${TMP_EXPOSURE_INPUT_DIR}/pheno.txt" \
  --covarFile "${TMP_EXPOSURE_INPUT_DIR}/covar.txt" \
  --bt \
  --firth --approx --pThresh 0.01 \
  --pred "${OUTDIR_STEP1}/regenie_step1_${Ename}_pred.list" \
  --bsize 400 \
  --threads 8 \
  --out "${OUTDIR_STEP2_SCRATCH}/regenie_step2"
  # NOTE: do NOT put ${Ename} in the --out basename. REGENIE appends
  # _<phenotype>.regenie (phenotype == Ename), so embedding Ename here DOUBLES it
  # and a long exposure (e.g. the 117-char physical-activity one) overflows the
  # 255-byte filename limit -> "cannot write file". The scratch dir already
  # namespaces per-exposure, and the copy step below globs *.regenie.

echo "[${Ename}] REGENIE binary GWAS complete"

# Copy final step2 summary stats to the canonical IGLOO location with a SIMPLE
# name: <exposure>.regenie (the per-exposure directory already namespaces it, so
# the old regenie_step2_<exp>_<exp>.regenie double-naming is dropped). rm before
# cp so any hpc_patel member can overwrite a previous run's file.
echo "[${Ename}] Copying step2 outputs to IGLOO: ${OUTDIR_STEP2_IGLOO}"
PRODUCED_REGENIE=$(ls "${OUTDIR_STEP2_SCRATCH}"/*.regenie 2>/dev/null | head -1)
PRODUCED_LOG=$(ls "${OUTDIR_STEP2_SCRATCH}"/*.log 2>/dev/null | head -1)
CANONICAL_FILE="${OUTDIR_STEP2_IGLOO}/${Ename}.regenie"

if [[ -n "${PRODUCED_REGENIE}" ]]; then
  rm -f "${CANONICAL_FILE}"
  cp -f "${PRODUCED_REGENIE}" "${CANONICAL_FILE}"
  chmod g+w "${CANONICAL_FILE}" 2>/dev/null || true
fi
if [[ -n "${PRODUCED_LOG}" ]]; then
  rm -f "${OUTDIR_STEP2_IGLOO}/${Ename}.log"
  cp -f "${PRODUCED_LOG}" "${OUTDIR_STEP2_IGLOO}/${Ename}.log"
  chmod g+w "${OUTDIR_STEP2_IGLOO}/${Ename}.log" 2>/dev/null || true
fi

if [[ ! -f "${CANONICAL_FILE}" ]]; then
  echo "WARNING: Expected canonical file not found after copy: ${CANONICAL_FILE}" >&2
  echo "  Check scratch: ${OUTDIR_STEP2_SCRATCH}/" >&2
else
  echo "[${Ename}] Canonical IGLOO file confirmed: ${CANONICAL_FILE}"
fi
