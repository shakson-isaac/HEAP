#!/bin/bash
#SBATCH -t 0-02:30
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH -p short
#SBATCH -J M6dzboot
umask 002

HEAP_ROOT="${HEAP_ROOT:-/n/groups/patel/shakson_ukb/HEAP}"
export HEAP_ROOT HEAP_PATHS_FILE="${HEAP_ROOT}/workflow/00_paths.R"
mkdir -p /n/groups/patel/IGLOO/UKB/HEAP/logs/module6
set -euo pipefail
module load gcc/14.2.0 R/4.4.2
ulimit -n 10000
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1 R_THREADS=1
export HEAP_BOOT_B="${HEAP_BOOT_B:-400}"
echo "[$(date)] M6 PES disease bootstrap: cpus=${SLURM_CPUS_PER_TASK}, B=${HEAP_BOOT_B}, host=${HOSTNAME}"
Rscript "${HEAP_ROOT}/scripts/support/module6_pes_disease_bootstrap.R"
echo "[$(date)] done"
