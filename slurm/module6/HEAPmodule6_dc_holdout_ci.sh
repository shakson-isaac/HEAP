#!/bin/bash
#SBATCH -J m6_dcci
#SBATCH -t 0-05:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH -p short
#SBATCH --array=1-165

# ============================================================
# 95% bootstrap intervals for every exposure x disease cell of the
# held-out Cox C-index screen (29,865 pairs), one exposure per task.
#
# Fits ONCE on the same 70% split the screen used, then resamples the
# 30% test set B=400 times. Measured 23.8x cheaper than refitting
# inside every replicate (7.20 s of refit vs 0.30 s of concordance
# per replicate), and the point estimates reproduce the shipped
# screen exactly because the fit is untouched.
#
#   ~181 diseases x ~126 s / 4 cores = ~1.6 h per task
#   ~1,045 core-hours total
#
# Submit:  sbatch HEAPmodule6_dc_holdout_ci.sh
# Collect: Rscript scripts/support/collect_dc_holdout_ci.R
# Feeds:   export_pes_lookup_deposits.R -> \Dataref{pes_disease_lookup}
# ============================================================

HEAP_ROOT="${HEAP_ROOT:-/n/groups/patel/shakson_ukb/HEAP}"
export HEAP_ROOT
export HEAP_PATHS_FILE="${HEAP_ROOT}/workflow/00_paths.R"
mkdir -p /n/groups/patel/IGLOO/UKB/HEAP/logs/module6

module load gcc/14.2.0 R/4.4.2

echo "=== dc holdout CI | exposure ${SLURM_ARRAY_TASK_ID} | $(date) ==="
Rscript "${HEAP_ROOT}/scripts/support/module6_pes_disease_dc_holdout_ci.R" "${SLURM_ARRAY_TASK_ID}"
echo "=== done ${SLURM_ARRAY_TASK_ID} | exit=$? | $(date) ==="
