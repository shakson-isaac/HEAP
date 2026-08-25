#!/bin/bash
#SBATCH -t 1-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH -p medium
#SBATCH -J poparch_master_rel

set -euo pipefail

HEAP_ROOT="${HEAP_ROOT:-/n/groups/patel/shakson_ukb/HEAP}"
umask 0002   # group-writable outputs for hpc_patel team runs
export HEAP_ROOT
export HEAP_PATHS_FILE="${HEAP_ROOT}/workflow/00_paths.R"

SCRIPTS_ROOT="${HEAP_ROOT}/scripts/population_architecture/scripts"
LOG_DIR="${HEAP_IGLOO_ROOT:-/n/groups/patel/IGLOO}/UKB/HEAP/logs/greml"
mkdir -p "${LOG_DIR}"

module load gcc/14.2.0
module load R/4.4.2

CONFIG="${CONFIG:-${HEAP_ROOT}/scripts/population_architecture/config/default_config.R}"
# Canonical IGLOO-rooted output (was repo-local HEAP/output).
PERM_OUTPUT_ROOT="${PERM_OUTPUT_ROOT:-${HEAP_IGLOO_ROOT:-/n/groups/patel/IGLOO}/UKB/HEAP/output/population_architecture}"
# Per-USER node-local base: bare /tmp is not setgid, so a hardcoded shared name
# (e.g. shi872_poparch_tmp) created by one member is unwritable by another member
# whose task lands on the same node -> mktemp -d fails with "Permission denied".
# Namespacing by ${USER} gives each member their own top-level dir (auto-cleaned).
SCRATCH_ROOT_BASE="${SCRATCH_ROOT_BASE:-${SLURM_TMPDIR:-${TMPDIR:-/tmp}}/${USER}_poparch_tmp}"
# GCTA: prefer shared IGLOO install, fall back to legacy UK_Biobank bin.
GCTA_BIN="${GCTA_BIN:-${HEAP_IGLOO_ROOT:-/n/groups/patel/IGLOO}/GCTA/gcta64}"
[ -x "${GCTA_BIN}" ] || GCTA_BIN="/n/groups/patel/shakson_ukb/UK_Biobank/bin/gcta-1.94.1-linux-kernel-3-x86_64/gcta64"

GRM_CUTOFF="${GRM_CUTOFF:-0.025}"
CUTOFF_LABEL="${CUTOFF_LABEL:-${GRM_CUTOFF//./p}}"
SPEC="${SPEC:-base}"
RUN_ID="${RUN_ID:-master_relatedness_${SPEC}_grmcutoff_${CUTOFF_LABEL}}"
FORCE="${FORCE:-true}"
PROFILE_PREP="${PROFILE_PREP:-true}"

scratch_group_dir="${SCRATCH_ROOT_BASE}/${RUN_ID}"
mkdir -p "${SCRATCH_ROOT_BASE}"
task_root="$(mktemp -d "${scratch_group_dir}_XXXXXX")"
task_output_root="${task_root}/output"
perm_dir="${PERM_OUTPUT_ROOT}/${SPEC}/relatedness/grm_cutoff_${CUTOFF_LABEL}"
center_dir="centered"

mkdir -p "${task_output_root}"
mkdir -p "${perm_dir}"

cleanup() {
  if [[ -n "${task_root:-}" && -d "${task_root}" ]]; then
    rm -rf "${task_root}"
  fi
}
trap cleanup EXIT

echo "[$(date '+%Y-%m-%d %H:%M:%S')] RUN_ID=${RUN_ID}"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] GRM cutoff=${GRM_CUTOFF}"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Scratch root=${task_root}"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Permanent relatedness dir=${perm_dir}"
df -h "${task_root}" || true

run_step() {
  local step_name="$1"
  shift
  if [[ "${PROFILE_PREP}" == "true" ]]; then
    /usr/bin/time -v -o "${LOG_DIR}/${RUN_ID}_${step_name}_time.log" "$@"
  else
    "$@"
  fi
}

run_step export \
  env "POPARCH_OUTPUT_ROOT=${task_output_root}" \
  Rscript "${SCRIPTS_ROOT}/export_architecture_inputs.R" \
  "${CONFIG}" "${RUN_ID}" "${SPEC}" \
  "--complete-case-exposures=false" \
  "--protein-specific-prep=false" \
  "--force=${FORCE}"

run_step genotype_grm \
  env "POPARCH_OUTPUT_ROOT=${task_output_root}" \
  Rscript "${SCRIPTS_ROOT}/build_genotype_grm.R" \
  "${CONFIG}" "${RUN_ID}" "${SPEC}" \
  "--force=${FORCE}"

prep_root="${task_output_root}/${RUN_ID}/${SPEC}/${center_dir}"
full_grm_prefix="${prep_root}/kernels/geno_ld_pruned"
singleton_prefix="${prep_root}/logs/master_unrelated_${CUTOFF_LABEL}"

/usr/bin/time -v -o "${LOG_DIR}/${RUN_ID}_singleton_${CUTOFF_LABEL}_time.log" \
  "${GCTA_BIN}" \
  --grm "${full_grm_prefix}" \
  --grm-singleton "${GRM_CUTOFF}" \
  --out "${singleton_prefix}"

singleton_keep="${singleton_prefix}.singleton.txt"
if [[ ! -f "${singleton_keep}" ]]; then
  echo "Expected singleton keep file was not created: ${singleton_keep}" >&2
  exit 1
fi

cp -f "${singleton_keep}" "${perm_dir}/master_unrelated.singleton.txt"
cp -f "${prep_root}/inputs/sample_manifest.tsv" "${perm_dir}/master_sample_manifest.tsv"
cp -f "${prep_root}/inputs/metadata.tsv" "${perm_dir}/master_metadata.tsv"
cp -f "${prep_root}/kernels/geno_ld_pruned.grm.id" "${perm_dir}/master_geno_ld_pruned.grm.id"

env KEEP_FILE="${singleton_keep}" \
    MANIFEST_FILE="${prep_root}/inputs/sample_manifest.tsv" \
    OUT_FILE="${perm_dir}/master_relatedness_summary.tsv" \
    CUTOFF_VALUE="${GRM_CUTOFF}" \
  Rscript - <<'EOF'
manifest <- read.delim(Sys.getenv("MANIFEST_FILE"), sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
keep <- read.delim(Sys.getenv("KEEP_FILE"), sep = "", header = FALSE, stringsAsFactors = FALSE)
if (ncol(keep) < 2L) stop("Keep file must have at least two columns.")
names(keep)[1:2] <- c("FID", "IID")
out <- data.frame(
  metric = c("n_master_samples", "n_unrelated_samples", "grm_cutoff"),
  value  = c(nrow(manifest), nrow(unique(keep[, 1:2, drop = FALSE])), Sys.getenv("CUTOFF_VALUE")),
  stringsAsFactors = FALSE
)
write.table(out, file = Sys.getenv("OUT_FILE"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
EOF

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Saved unrelated keep file to ${perm_dir}/master_unrelated.singleton.txt"
