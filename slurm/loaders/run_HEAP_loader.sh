#!/bin/bash
#SBATCH --job-name=HEAP_loader
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH -t 0-01:30
#SBATCH -p short
umask 002  # group-writable outputs for hpc_patel multi-user runs

# Usage:
#   sbatch run_HEAP_loader.sh   # full run (~45 min, always includes disease)
#
# Output (canonical IGLOO location):
#   /n/groups/patel/IGLOO/UKB/HEAP/intermediate/HEAP.rds
#   /n/groups/patel/IGLOO/UKB/HEAP/intermediate/HEAP_audit_*.tsv
#
# Override output path (temporary/testing only):
#   HEAP_LOADER_RDS=/path/to/custom.rds sbatch run_HEAP_loader.sh

set -euo pipefail

HEAP_ROOT="${HEAP_ROOT:-/n/groups/patel/shakson_ukb/HEAP}"
export HEAP_ROOT

module load gcc/14.2.0
module load R/4.4.2

ulimit -n 10000

echo "[$(date '+%Y-%m-%d %H:%M:%S')] HEAP_loader starting"
echo "  HEAP_ROOT         : ${HEAP_ROOT}"
echo "  HEAP_SKIP_DISEASE : ${HEAP_SKIP_DISEASE:-FALSE}"
echo "  Job ID            : ${SLURM_JOB_ID}"
echo "  Node              : $(hostname)"
echo "  Memory            : ${SLURM_MEM_PER_NODE}MB"

START=$(date +%s)

Rscript "${HEAP_ROOT}/scripts/loaders/HEAP_loader.R"
STATUS=$?

END=$(date +%s)
ELAPSED=$(( END - START ))
ELAPSED_MIN=$(( ELAPSED / 60 ))

echo ""
echo "[$(date '+%Y-%m-%d %H:%M:%S')] HEAP_loader finished"
echo "  Exit status : ${STATUS}"
echo "  Elapsed     : ${ELAPSED}s (${ELAPSED_MIN} min)"

if [ "$STATUS" -eq 0 ]; then
  IGLOO_HEAP="${HEAP_IGLOO_ROOT:-/n/groups/patel/IGLOO}/UKB/HEAP"
  RDS="${HEAP_LOADER_RDS:-${IGLOO_HEAP}/intermediate/HEAP.rds}"
  if [ -f "$RDS" ]; then
    SIZE=$(du -sh "$RDS" | cut -f1)
    echo "  Output size : ${SIZE}  -> ${RDS}"
    echo ""
    echo "Audit files written:"
    ls -lh "${IGLOO_HEAP}/intermediate/HEAP_audit_"*.tsv 2>/dev/null || echo "  (none found)"
  fi
fi

exit $STATUS
