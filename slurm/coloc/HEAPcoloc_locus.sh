#!/bin/bash
umask 002  # group-writable outputs for hpc_patel multi-user runs
# Per-locus colocalization array job.
#
# Re-runs scripts/support/coloc/run_coloc_locus.R over every row of the coloc
# manifest so that each locus RETAINS its <out>_harmonized_snps.tsv -- the
# harmonized pQTL x disease SNP table that coloc.abf actually consumed.
#
# Why: the summaries shipped, the per-SNP tables did not. Only one locus
# (ASGR1 x E4_LIPOPROT) still had its harmonized table, so the regional
# colocalization view on heap.bio could be drawn for exactly one of 65 pairs.
# The runner already writes the file (run_coloc_locus.R:295); it just needs to
# be run again with the output kept.
#
# Downstream of this:
#   scripts/support/coloc/export_coloc_web.R
#     -> web/<locus>_plot_table.tsv (snp, chr, pos, p_trait1, p_trait2, r2),
#        read by the colocalization locus plot and the heap.bio regional view.
#
# USAGE
# -----
#   bash slurm/coloc/HEAPcoloc_locus.sh
#   ARRAY_TASKS=3,7,11 bash slurm/coloc/HEAPcoloc_locus.sh   # rerun failures
#
# ENV
# ---
#   MANIFEST     override the manifest path
#   ARRAY_TASKS  explicit array spec instead of 1-N
#   TIME/MEM/CPUS/PARTITION
#   DEPENDENCY   sbatch --dependency spec (O2 ignores SBATCH_DEPENDENCY)
set -euo pipefail

REPO="${HEAP_REPO:-/n/groups/patel/shakson_ukb/HEAP}"
IGLOO_HEAP="${HEAP_IGLOO_ROOT:-/n/groups/patel/IGLOO}/UKB/HEAP"
MANIFEST="${MANIFEST:-${IGLOO_HEAP}/output/support/coloc/coloc_manifest.tsv}"
OUT_DIR="${IGLOO_HEAP}/output/support/coloc/per_locus"

if [[ ! -f "${MANIFEST}" ]]; then
  echo "ERROR: coloc manifest not found: ${MANIFEST}" >&2
  echo "Build it first: Rscript scripts/support/coloc/build_coloc_manifest.R" >&2
  exit 1
fi

N_ROWS=$(tail -n +2 "${MANIFEST}" | wc -l)
[[ "${N_ROWS}" -eq 0 ]] && { echo "ERROR: manifest has no rows" >&2; exit 1; }

# One 500kb window per task: reads a slice of the pQTL and FinnGen sumstats,
# harmonizes, runs coloc.abf. I/O-bound and short.
TIME="${TIME:-0-02:00}"
MEM="${MEM:-16G}"
CPUS="${CPUS:-2}"
PARTITION="${PARTITION:-short}"
ARRAY_SPEC="${ARRAY_TASKS:-1-${N_ROWS}}"

LOG_DIR="${IGLOO_HEAP}/logs/coloc"
mkdir -p "${LOG_DIR}" "${OUT_DIR}"

# ---- submitter branch: submit self as an array, then exit -------------------
# Guard on SLURM_ARRAY_TASK_ID, not SLURM_JOB_ID -- the latter is also set
# inside interactive/VS Code allocations.
if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] Submitting coloc per-locus array"
  echo "  Manifest   : ${MANIFEST}"
  echo "  N rows     : ${N_ROWS}"
  echo "  Array spec : ${ARRAY_SPEC}"
  echo "  Out dir    : ${OUT_DIR}"
  echo "  Resources  : time=${TIME}, mem=${MEM}, cpus=${CPUS}, partition=${PARTITION}"

  dep_arg=()
  [[ -n "${DEPENDENCY:-}" ]] && dep_arg=(--dependency="${DEPENDENCY}")
  exec sbatch \
    -t "${TIME}" \
    --ntasks=1 \
    --cpus-per-task="${CPUS}" \
    --mem="${MEM}" \
    -p "${PARTITION}" \
    --array="${ARRAY_SPEC}" \
    --export=ALL \
    "${dep_arg[@]}" \
    -J "coloc_locus" \
    -o "${LOG_DIR}/coloc_%A_%a.out" \
    -e "${LOG_DIR}/coloc_%A_%a.err" \
    --mail-type=FAIL \
    "$0"
fi

# ---- job body --------------------------------------------------------------
module load gcc/14.2.0
module load R/4.4.2

ROW="${SLURM_ARRAY_TASK_ID}"
cd "${REPO}"

# SLURM_TMPDIR is unset on O2 and /tmp is shared, so scratch is per user+job.
SCRATCH="/tmp/${USER}/coloc_${SLURM_ARRAY_JOB_ID}_${ROW}"
mkdir -p "${SCRATCH}"
export TMPDIR="${SCRATCH}"
trap 'rm -rf "${SCRATCH}"' EXIT

echo "[$(date '+%Y-%m-%d %H:%M:%S')] coloc row ${ROW}/${N_ROWS}"
awk -F'\t' -v r="${ROW}" 'NR==1{for(i=1;i<=NF;i++)h[i]=$i} NR==r+1{
  for(i=1;i<=NF;i++) if (h[i]=="arm"||h[i]=="protID"||h[i]=="disease_or_exposure_id"||h[i]=="edge_dir"||h[i]=="lead_snp")
    printf "  %-24s %s\n", h[i], $i }' "${MANIFEST}"

# No --out: that argument is a full PREFIX, not a directory, so passing a
# directory would write "<dir>_coloc_summary.tsv". The runner's own default
# (support/coloc/per_locus/<arm>__<protein>__<disease>) is what we want -- it
# names the arm explicitly, matches the summary naming, and lands beside rather
# than on top of the legacy flat-directory outputs.
Rscript scripts/support/coloc/run_coloc_locus.R \
  --manifest "${MANIFEST}" \
  --row "${ROW}"

echo "[$(date '+%Y-%m-%d %H:%M:%S')] row ${ROW} done"
