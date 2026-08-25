#!/bin/bash
umask 002  # group-writable outputs for hpc_patel multi-user runs
# Submit the full GREML GRM-cutoff workflow.
#
# Usage (no env vars required — just run it):
#   bash submit_greml_workflow.sh
#
# Optional overrides (all have sensible defaults):
#   RUN_GROUP=heap_v1  (run label; namespaces scratch/prep dirs + time-log names,
#                       NOT the output path — final outputs are SPEC/cutoff/MODEL-keyed)
#   GRM_CUTOFF=0.025   SPEC=base   MODEL=primary   MIN_PROTEIN_N=2000
#   PERM_OUTPUT_ROOT=...   PROTEIN_SET=...
#
# Workflow steps (each depends on the previous):
#   1. BuildMasterUnrelated  — build whole-cohort GRM, derive singleton keep list
#   2. ProtGremlGrmCutoff    — per-protein GREML array job (2923 tasks by default)
#   3. SummaryGremlGrmCutoff — aggregate per-protein rows into combined TSV

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
HEAP_ROOT="${HEAP_ROOT:-/n/groups/patel/shakson_ukb/HEAP}"

# Run label, baked in with a default so the workflow runs as a one-liner.
# Override by exporting RUN_GROUP before calling, e.g. for a distinct concurrent run.
# Note: this does NOT affect dashboard detection (heap_status.R keys on the
# `poparch*` SLURM job names + `*_summary.tsv` output files, not RUN_GROUP).
RUN_GROUP="${RUN_GROUP:-heap_v1}"

export RUN_GROUP HEAP_ROOT

GRM_CUTOFF="${GRM_CUTOFF:-0.025}"
CUTOFF_LABEL="${GRM_CUTOFF//./p}"
SPEC="${SPEC:-base}"
MODEL="${MODEL:-primary}"
PERM_OUTPUT_ROOT="${PERM_OUTPUT_ROOT:-${HEAP_IGLOO_ROOT:-/n/groups/patel/IGLOO}/UKB/HEAP/output/population_architecture}"
PROTEIN_SET="${PROTEIN_SET:-${HEAP_ROOT}/scripts/population_architecture/config/protein_sets/all_proteins_from_loader.txt}"

if [[ ! -f "${PROTEIN_SET}" ]]; then
  echo "Protein set not found: ${PROTEIN_SET}" >&2
  exit 1
fi
N_PROTEINS="$(grep -cv '^[[:space:]]*$' "${PROTEIN_SET}")"

echo "=== GREML GRM-cutoff workflow ==="
echo "  RUN_GROUP          = ${RUN_GROUP}"
echo "  GRM_CUTOFF         = ${GRM_CUTOFF}"
echo "  SPEC               = ${SPEC}"
echo "  MODEL              = ${MODEL}"
echo "  PERM_OUTPUT_ROOT   = ${PERM_OUTPUT_ROOT}"
echo "  Proteins           = ${N_PROTEINS}"
echo ""

# Validate a job id is a real Slurm number before chaining afterok on it.
require_jid() { [[ "${1:-}" =~ ^[0-9]+$ ]] || { echo "ERROR: sbatch failed for ${2} (got JID='${1:-}')." >&2; exit 1; }; }

# Step 1: Build master unrelated keep file
STEP1_JID=$(sbatch \
  --export=ALL \
  --parsable \
  "${SCRIPT_DIR}/BuildMasterUnrelated.sh")
require_jid "${STEP1_JID}" "BuildMasterUnrelated"
echo "Submitted BuildMasterUnrelated: job ${STEP1_JID}"

# Step 2: Per-protein GREML array (depends on step 1)
STEP2_JID=$(sbatch \
  --export=ALL \
  --array="1-${N_PROTEINS}" \
  --dependency="afterok:${STEP1_JID}" \
  --parsable \
  "${SCRIPT_DIR}/ProtGremlGrmCutoff.sh")
require_jid "${STEP2_JID}" "ProtGremlGrmCutoff"
echo "Submitted ProtGremlGrmCutoff array (${N_PROTEINS} tasks): job ${STEP2_JID}"

# Step 3: Aggregate per-protein summary rows (depends on all of step 2)
STEP3_JID=$(sbatch \
  --export=ALL \
  --dependency="afterok:${STEP2_JID}" \
  --parsable \
  "${SCRIPT_DIR}/SummaryGremlGrmCutoff.sh")
require_jid "${STEP3_JID}" "SummaryGremlGrmCutoff"
echo "Submitted SummaryGremlGrmCutoff: job ${STEP3_JID}"

echo ""
echo "Output will be written to:"
echo "  ${PERM_OUTPUT_ROOT}/${SPEC}/grm_cutoff_${CUTOFF_LABEL}/"
echo "Logs: ${HEAP_ROOT}/logs/greml/"
