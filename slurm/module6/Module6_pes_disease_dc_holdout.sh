#!/bin/bash
#SBATCH -t 0-09:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=64G
#SBATCH -p short
#SBATCH -J M6dcho
umask 002  # group-writable outputs for hpc_patel multi-user runs

HEAP_ROOT="${HEAP_ROOT:-/n/groups/patel/shakson_ukb/HEAP}"
export HEAP_ROOT
export HEAP_PATHS_FILE="${HEAP_ROOT}/workflow/00_paths.R"
mkdir -p /n/groups/patel/IGLOO/UKB/HEAP/logs/module6

set -euo pipefail
module load gcc/14.2.0
module load R/4.4.2
ulimit -n 10000

# keep BLAS single-threaded; parallelism is across exposures via mclapply
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1 R_THREADS=1
export HEAP_PES_MIN_EVENTS="${HEAP_PES_MIN_EVENTS:-100}"

# Optional specification. Unset (the default) runs base exactly as before and
# writes pes_disease_dc_holdout.tsv; SPEC=base_exclprev restricts to participants
# healthy at baseline, reads the score retrained on that subset, and writes to a
# separate file so the shipped base screen is never clobbered.
#   sbatch --export=ALL,SPEC=base_exclprev Module6_pes_disease_dc_holdout.sh
SPEC="${SPEC:-}"

echo "[$(date)] M6 PES dC holdout: spec=${SPEC:-base}, cpus=${SLURM_CPUS_PER_TASK}, host=${HOSTNAME}, min_events=${HEAP_PES_MIN_EVENTS}"
Rscript "${HEAP_ROOT}/scripts/support/module6_pes_disease_dc_holdout.R" ${SPEC}
echo "[$(date)] M6 PES dC holdout completed"
