#!/bin/bash
#SBATCH -J poparch_cutoff_fin
#SBATCH -p short
#SBATCH -t 0-00:20:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G

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
PERM_OUTPUT_ROOT="${PERM_OUTPUT_ROOT:-${HEAP_IGLOO_ROOT:-/n/groups/patel/IGLOO}/UKB/HEAP/output/population_architecture}"

: "${RUN_GROUP:?Set RUN_GROUP at submission time}"
GRM_CUTOFF="${GRM_CUTOFF:-0.025}"
CUTOFF_LABEL="${CUTOFF_LABEL:-${GRM_CUTOFF//./p}}"
GREML_RUN_ID="${GREML_RUN_ID:-${RUN_GROUP}_grmcutoff_${CUTOFF_LABEL}}"

SPEC="${SPEC:-base}"
CENTER_EXPOSURES="${CENTER_EXPOSURES:-true}"
MODEL="${MODEL:-primary}"
MIN_PROTEIN_N="${MIN_PROTEIN_N:-2000}"

model_dir="${PERM_OUTPUT_ROOT}/${SPEC}/grm_cutoff_${CUTOFF_LABEL}/${MODEL}"
summary_dir="${model_dir}"
plots_dir="${model_dir}/plots"

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Aggregating GREML cutoff run ${GREML_RUN_ID}"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Looking for per-protein rows in ${model_dir}"

if compgen -G "${model_dir}/*_summary.tsv" > /dev/null; then
  Rscript "${SCRIPTS_ROOT}/summarize_population_architecture.R" \
    "${CONFIG}" "${GREML_RUN_ID}" "${SPEC}" \
    "--model=${MODEL}" \
    "--center-exposures=${CENTER_EXPOSURES}" \
    "--model-dir=${model_dir}" \
    "--summary-dir=${summary_dir}" \
    "--plots-dir=${plots_dir}" \
    "--min-protein-n=${MIN_PROTEIN_N}"
else
  echo "No per-protein summary rows found; skipping."
  exit 0
fi
