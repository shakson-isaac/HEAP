#!/bin/bash
# Rebuild the HEAP MR triad set end-to-end, in dependency order.
#
#   1. Refresh the Module 2 connector  (ReplicatedEassoc.csv, from the base run)
#   2. Rebuild the triad/edge lists     (Module5_load.R -> global_edges/*, mr_triads.tsv)
#   3. Regenerate the Module 5 manifests (per-edge-type chunking)
#
# See TRIAD_CONNECTIONS.md for what each step connects and the resulting counts.
#
# Usage:
#   bash scripts/module5_mr/rebuild_triads.sh
#   COVARTYPE=base EXPERIMENTS="MR_UKB_primary MR_deCODE_replication" \
#     bash scripts/module5_mr/rebuild_triads.sh
set -euo pipefail
umask 002

HEAP_ROOT="${HEAP_ROOT:-/n/groups/patel/shakson_ukb/HEAP}"
export HEAP_ROOT
export HEAP_PATHS_FILE="${HEAP_ROOT}/workflow/00_paths.R"
COVARTYPE="${COVARTYPE:-base}"
EXPERIMENTS="${EXPERIMENTS:-MR_UKB_primary MR_deCODE_replication}"

module load gcc/14.2.0 2>/dev/null || true
module load R/4.4.2   2>/dev/null || true

echo "==> [1/3] Refresh Module 2 connector (ReplicatedEassoc.csv, covarType=${COVARTYPE})"
Rscript "${HEAP_ROOT}/scripts/module2_associations/summarize_replicated_associations.R" "${COVARTYPE}"

echo "==> [2/3] Rebuild triad/edge lists (Module5_load.R)"
Rscript "${HEAP_ROOT}/scripts/module5_mr/Module5_load.R"

echo "==> [3/3] Regenerate Module 5 manifests: ${EXPERIMENTS}"
Rscript -e "
  source('${HEAP_ROOT}/workflow/00_paths.R')
  source('${HEAP_ROOT}/workflow/config_helpers.R')
  source('${HEAP_ROOT}/workflow/generate_manifests.R')
  for (e in strsplit('${EXPERIMENTS}', ' +')[[1]]) if (nzchar(e)) generate_module5_manifest(e)
"

echo "==> DONE. Triad set rebuilt. Inventory: \$(global_edges)/mr_triads.tsv"
