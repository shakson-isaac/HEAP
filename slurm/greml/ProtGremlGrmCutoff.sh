#!/bin/bash
#SBATCH -t 0-12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=60G
#SBATCH -p short
#SBATCH -J poparch_cutoff
#SBATCH --array=1-2923

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
PROTEIN_SET="${PROTEIN_SET:-${HEAP_ROOT}/scripts/population_architecture/config/protein_sets/all_proteins_from_loader.txt}"
IGLOO_HEAP="${HEAP_IGLOO_ROOT:-/n/groups/patel/IGLOO}/UKB/HEAP"
# Canonical permanent output root: IGLOO-rooted for reproducibility.
# Per-protein summaries are first written to node-local scratch, then
# copied to PERM_OUTPUT_ROOT before job exit (see cp command below).
PERM_OUTPUT_ROOT="${PERM_OUTPUT_ROOT:-${IGLOO_HEAP}/output/population_architecture}"
export POPARCH_OUTPUT_ROOT="${PERM_OUTPUT_ROOT}"
# Node-local scratch for fast per-job working directory (auto-cleaned).
# Per-USER node-local base: bare /tmp is not setgid, so a hardcoded shared name
# (e.g. shi872_poparch_tmp) created by one member is unwritable by another member
# whose task lands on the same node -> mktemp -d fails with "Permission denied".
# Namespacing by ${USER} gives each member their own top-level dir (auto-cleaned).
SCRATCH_ROOT_BASE="${SCRATCH_ROOT_BASE:-${SLURM_TMPDIR:-${TMPDIR:-/tmp}}/${USER}_poparch_tmp}"

: "${RUN_GROUP:?Set RUN_GROUP at submission time}"

# If EXPERIMENT is set, read GRM_CUTOFF, SPEC, and MIN_PROTEIN_N from
# config/modules/population_architecture_experiments.yml.
# This makes reviewer-facing parameters fully config-driven.
# EXPERIMENT and individual env vars can still override each other:
#   - EXPERIMENT provides defaults from the named config entry
#   - GRM_CUTOFF / SPEC env vars override those defaults if also set
EXPERIMENT="${EXPERIMENT:-}"
if [[ -n "${EXPERIMENT}" ]]; then
  _EXP_JSON=$(Rscript --vanilla -e "
    source('${HEAP_PATHS_FILE:-${HEAP_ROOT}/workflow/00_paths.R}')
    h <- tryCatch(source('${HEAP_ROOT}/workflow/config_helpers.R'), error=function(e) NULL)
    if (!exists('load_experiment_config', mode='function')) quit(status=0)
    exp <- tryCatch(
      load_experiment_config('population_architecture', '${EXPERIMENT}'),
      error = function(e) { message(e); NULL }
    )
    if (is.null(exp)) quit(status=0)
    cat('GRM_CUTOFF_CFG=', exp\$grm_cutoff %||% 0.025, '\n', sep='')
    cat('SPEC_CFG=',        exp\$covariate_set %||% 'base', '\n', sep='')
    cat('MIN_PROTEIN_N_CFG=', exp\$min_protein_n %||% 2000, '\n', sep='')
  " 2>/dev/null || true)

  # Source the key=value pairs emitted by Rscript
  eval "${_EXP_JSON}" 2>/dev/null || true

  # Apply config values only if the caller did not explicitly set the env var
  GRM_CUTOFF="${GRM_CUTOFF:-${GRM_CUTOFF_CFG:-0.025}}"
  SPEC="${SPEC:-${SPEC_CFG:-base}}"
  MIN_PROTEIN_N="${MIN_PROTEIN_N:-${MIN_PROTEIN_N_CFG:-2000}}"

  echo "[config] Experiment '${EXPERIMENT}': GRM_CUTOFF=${GRM_CUTOFF} SPEC=${SPEC} MIN_PROTEIN_N=${MIN_PROTEIN_N}"
fi

GRM_CUTOFF="${GRM_CUTOFF:-0.025}"
CUTOFF_LABEL="${CUTOFF_LABEL:-${GRM_CUTOFF//./p}}"
GREML_RUN_ID="${GREML_RUN_ID:-${RUN_GROUP}_grmcutoff_${CUTOFF_LABEL}}"

SPEC="${SPEC:-base}"
CENTER_EXPOSURES="${CENTER_EXPOSURES:-true}"
INCLUDE_COVAR_KERNELS="${INCLUDE_COVAR_KERNELS:-true}"
MODEL="${MODEL:-primary}"
MIN_PROTEIN_N="${MIN_PROTEIN_N:-2000}"
PREP_FORCE="${PREP_FORCE:-false}"
GREML_FORCE="${GREML_FORCE:-true}"
PROFILE_PREP="${PROFILE_PREP:-true}"
CONTINUE_ON_ERROR="${CONTINUE_ON_ERROR:-true}"
MASTER_KEEP_FILE="${MASTER_KEEP_FILE:-${PERM_OUTPUT_ROOT}/${SPEC}/relatedness/grm_cutoff_${CUTOFF_LABEL}/master_unrelated.singleton.txt}"

if [[ ! -f "${PROTEIN_SET}" ]]; then
  echo "Protein set file not found: ${PROTEIN_SET}" >&2
  exit 1
fi

mapfile -t proteins < <(grep -v '^[[:space:]]*$' "${PROTEIN_SET}")
if [[ "${#proteins[@]}" -eq 0 ]]; then
  echo "No proteins found in ${PROTEIN_SET}" >&2
  exit 1
fi

if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
  echo "SLURM_ARRAY_TASK_ID is required." >&2
  exit 1
fi

task_index=$((SLURM_ARRAY_TASK_ID - 1))
if (( task_index < 0 || task_index >= ${#proteins[@]} )); then
  echo "SLURM_ARRAY_TASK_ID ${SLURM_ARRAY_TASK_ID} is out of range for ${#proteins[@]} proteins." >&2
  exit 1
fi

protein="${proteins[$task_index]}"
protein_slug="$(echo "${protein}" | tr '[:upper:]' '[:lower:]')"
prep_run_id="${RUN_GROUP}_${protein_slug}_prep"
profile_label="${GREML_RUN_ID}_${MODEL}_${protein}"
center_dir="$( [[ "${CENTER_EXPOSURES}" == "true" ]] && echo centered || echo uncentered )"
perm_model_dir="${PERM_OUTPUT_ROOT}/${SPEC}/grm_cutoff_${CUTOFF_LABEL}/${MODEL}"
scratch_group_dir="${SCRATCH_ROOT_BASE}/${RUN_GROUP}"

mkdir -p "${scratch_group_dir}"
task_root="$(mktemp -d "${scratch_group_dir}/${SLURM_ARRAY_TASK_ID}_${protein_slug}_XXXXXX")"
task_output_root="${task_root}/output"
mkdir -p "${task_output_root}"
mkdir -p "${perm_model_dir}"
prep_root="${task_output_root}/${prep_run_id}/${SPEC}/${center_dir}"

cleanup() {
  if [[ -n "${task_root:-}" && -d "${task_root}" ]]; then
    rm -rf "${task_root}"
  fi
}
trap cleanup EXIT

# Sanity-check: PERM_OUTPUT_ROOT must be under IGLOO, not local or scratch
if [[ "${PERM_OUTPUT_ROOT}" == *"/n/scratch/"* || "${PERM_OUTPUT_ROOT}" == *"${HEAP_ROOT}/output"* ]]; then
  echo "ERROR: PERM_OUTPUT_ROOT must be IGLOO-rooted, not scratch or local repo output." >&2
  echo "  Current value: ${PERM_OUTPUT_ROOT}" >&2
  echo "  Expected prefix: ${IGLOO_HEAP}/output/..." >&2
  exit 1
fi

echo "[$(date '+%Y-%m-%d %H:%M:%S')] RUN_GROUP=${RUN_GROUP}"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] GREML_RUN_ID=${GREML_RUN_ID}"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Task ${SLURM_ARRAY_TASK_ID}/${#proteins[@]} protein=${protein}"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Prep run: ${prep_run_id}"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] GRM cutoff=${GRM_CUTOFF}"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Master keep file=${MASTER_KEEP_FILE}"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Node-local scratch root=${task_root}"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Permanent (IGLOO) model dir=${perm_model_dir}"
# Note: POPARCH_OUTPUT_ROOT is passed per R subprocess as an explicit env override
# (env "POPARCH_OUTPUT_ROOT=${task_output_root}") so computation uses the node-local
# tmp dir. The final summary TSV is cp-ed to perm_model_dir (IGLOO) before job exit.
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Export script mtime=$(stat -c '%y' "${SCRIPTS_ROOT}/export_architecture_inputs.R" 2>/dev/null)"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Filter script mtime=$(stat -c '%y' "${SCRIPTS_ROOT}/filter_prep_inputs_by_keep.R" 2>/dev/null)"
df -h "${task_root}" || true

if [[ ! -f "${MASTER_KEEP_FILE}" ]]; then
  echo "Master relatedness keep file not found: ${MASTER_KEEP_FILE}" >&2
  echo "Run BuildMasterUnrelated.sh first." >&2
  exit 1
fi

COMPLETE_CASE_EXPOSURES="true"
PROTEIN_SPECIFIC_PREP="true"
MAX_SAMPLES="${MAX_SAMPLES:-}"
SEED="${SEED:-1}"

run_step() {
  local step_name="$1"
  shift
  if [[ "${PROFILE_PREP}" == "true" ]]; then
    /usr/bin/time -v -o "${LOG_DIR}/${prep_run_id}_${step_name}_time.log" "$@"
  else
    "$@"
  fi
}

export_args=(
  "${CONFIG}" "${prep_run_id}" "${SPEC}"
  "--proteins=${protein}"
  "--complete-case-exposures=${COMPLETE_CASE_EXPOSURES}"
  "--protein-specific-prep=${PROTEIN_SPECIFIC_PREP}"
  "--force=${PREP_FORCE}"
)

if [[ -n "${MAX_SAMPLES}" ]]; then
  export_args+=("--max-samples=${MAX_SAMPLES}" "--seed=${SEED}")
fi

run_step export \
  env "POPARCH_OUTPUT_ROOT=${task_output_root}" \
  Rscript "${SCRIPTS_ROOT}/export_architecture_inputs.R" \
  "${export_args[@]}"

metadata_path="${prep_root}/inputs/metadata.tsv"
if [[ -f "${metadata_path}" ]]; then
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] Export metadata snapshot:"
  grep -E '^(n_samples|n_samples_before_exposure_complete|n_samples_after_exposure_complete|n_exposures_retained_complete_case|n_exposures_dropped_pre_complete_case)' "${metadata_path}" || true
fi

run_step filter_relatedness \
  env "POPARCH_OUTPUT_ROOT=${task_output_root}" \
  Rscript "${SCRIPTS_ROOT}/filter_prep_inputs_by_keep.R" \
  "${CONFIG}" "${prep_run_id}" "${SPEC}" \
  "--keep-file=${MASTER_KEEP_FILE}" \
  "--center-exposures=${CENTER_EXPOSURES}" \
  "--cutoff=${GRM_CUTOFF}"

if [[ -f "${metadata_path}" ]]; then
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] Post-relatedness metadata snapshot:"
  grep -E '^(n_samples|n_samples_before_relatedness_cutoff|n_samples_after_relatedness_cutoff|grm_relatedness_cutoff)' "${metadata_path}" || true
fi

run_step genotype_grm \
  env "POPARCH_OUTPUT_ROOT=${task_output_root}" \
  Rscript "${SCRIPTS_ROOT}/build_genotype_grm.R" \
  "${CONFIG}" "${prep_run_id}" "${SPEC}" \
  "--force=${PREP_FORCE}"

run_step environment_kernels \
  env "POPARCH_OUTPUT_ROOT=${task_output_root}" \
  Rscript "${SCRIPTS_ROOT}/build_environment_kernels.R" \
  "${CONFIG}" "${prep_run_id}" "${SPEC}" \
  "--center-exposures=${CENTER_EXPOSURES}" \
  "--include-covar-kernels=${INCLUDE_COVAR_KERNELS}" \
  "--force=true"

/usr/bin/time -v -o "${LOG_DIR}/${profile_label}_cutoff_${CUTOFF_LABEL}_time.log" \
  env "POPARCH_OUTPUT_ROOT=${task_output_root}" \
  Rscript "${SCRIPTS_ROOT}/run_population_architecture.R" \
    "${CONFIG}" "${GREML_RUN_ID}" "${SPEC}" \
    "--prep-run-id=${prep_run_id}" \
    "--model=${MODEL}" \
    "--proteins=${protein}" \
    "--center-exposures=${CENTER_EXPOSURES}" \
    "--min-protein-n=${MIN_PROTEIN_N}" \
    "--continue-on-error=${CONTINUE_ON_ERROR}" \
    "--keep-artifact-paths=false" \
    "--write-combined-summary=false" \
    "--force=${GREML_FORCE}"

scratch_row_path="${task_output_root}/${GREML_RUN_ID}/${SPEC}/${center_dir}/models/${MODEL}/${protein}_summary.tsv"
perm_row_path="${perm_model_dir}/${protein}_summary.tsv"
if [[ ! -f "${scratch_row_path}" ]]; then
  echo "Expected summary row was not created: ${scratch_row_path}" >&2
  exit 1
fi

env CUTOFF_VALUE="${GRM_CUTOFF}" ROW_PATH="${scratch_row_path}" Rscript - <<'EOF'
df <- read.delim(Sys.getenv("ROW_PATH"), sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
df$grm_cutoff <- as.numeric(Sys.getenv("CUTOFF_VALUE"))
write.table(df, file = Sys.getenv("ROW_PATH"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
EOF

cp -f "${scratch_row_path}" "${perm_row_path}"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Saved summary row to ${perm_row_path}"
