#!/bin/bash
#SBATCH --job-name=exposure_coding_check
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=25G
#SBATCH -t 0-00:20
#SBATCH -p short
umask 002  # group-writable outputs for hpc_patel multi-user runs

set -euo pipefail

HEAP_ROOT="${HEAP_ROOT:-/n/groups/patel/shakson_ukb/HEAP}"
export HEAP_ROOT

module load gcc/14.2.0
module load R/4.4.2

echo "[$(date '+%Y-%m-%d %H:%M:%S')] exposure_coding_check starting"
echo "  HEAP_ROOT : ${HEAP_ROOT}"
echo "  Job ID    : ${SLURM_JOB_ID}"

Rscript "${HEAP_ROOT}/scripts/visualizations/exposure_coding_check.R"

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Done"
echo "PDFs: ${HEAP_ROOT}/output/exposure_coding_check/"
