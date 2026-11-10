#!/usr/bin/env bash
# ============================================================================
# demo/run_demo.sh -- run the unmodified HEAP analysis code on small synthetic
# datasets with planted effects. No UK Biobank data, cluster or SLURM needed.
#
#   bash demo/run_demo.sh                 # every stage
#   bash demo/run_demo.sh m2 m1           # only these stages (in pipeline order)
#   DEMO_DIR=/path bash demo/run_demo.sh  # output directory (default ./heap_demo_run)
#
# Stages, in the paper's order: m1 (variance decomposition), m2 (associations),
#         m3 (mediation; needs m1), mr (Mendelian randomization), m6 (proteome-based
#         exposure score).
#
# The code runs from $DEMO_DIR/repo: links to this repository's scripts/,
# workflow/, slurm/ and demo/, plus a copy of config/ in which the exposure and protein
# lists carry the demo's generic labels. Every data path is redirected through
# the environment variables a full run uses (docs/ENVIRONMENT.md), so nothing
# outside DEMO_DIR is read or written. Set HEAP_RLIB if your R packages live in
# a non-default library.
# ============================================================================
set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DEMO="${DEMO_DIR:-heap_demo_run}"; mkdir -p "$DEMO"; DEMO="$(cd "$DEMO" && pwd)"
N_PARTICIPANTS="${N_PARTICIPANTS:-2000}"
N_CHUNKS=2                       # 6 proteins, run as 2 array tasks of 3
ALL_STAGES=(m1 m2 m3 mr m6)
STAGES=("$@"); [ ${#STAGES[@]} -eq 0 ] && STAGES=("${ALL_STAGES[@]}")

# ---- the repository as the demo sees it ------------------------------------
ROOT="$DEMO/repo"
if [ ! -d "$ROOT" ]; then
  mkdir -p "$ROOT"
  for d in scripts workflow demo slurm; do ln -s "$REPO/$d" "$ROOT/$d"; done
  cp -r "$REPO/config" "$ROOT/config"        # demo lists overwrite two files below
fi

export HEAP_REPO="$REPO" HEAP_ROOT="$ROOT"
export HEAP_PATHS_FILE="$ROOT/workflow/00_paths.R"
export HEAP_IGLOO_ROOT="$DEMO/igloo"
export HEAP_SCRATCH_ROOT="$DEMO/scratch"
export HEAP_PROJECT_ROOT="$DEMO/igloo/UKB/HEAP"
export HEAP_LOADER_RDS="$DEMO/HEAP_synthetic.rds"
OUT="$HEAP_PROJECT_ROOT/output"
mkdir -p "$HEAP_SCRATCH_ROOT"
cd "$ROOT"

say() { echo "[$(date '+%H:%M:%S')] $*"; }
run() {  # run <log-name> <command...>: log to DEMO, show the tail on failure
  local log="$DEMO/$1.log"; shift
  "$@" > "$log" 2>&1 || { echo "FAILED: $*  (log: $log)"; tail -20 "$log"; exit 1; }
}

stage_m2() {
  say "Module 2 -- exposure, genetic and GxE associations (base)"
  for i in $(seq 1 "$N_CHUNKS"); do
    run "module2_$i" Rscript scripts/module2_associations/Module2.R "$i" "$N_CHUNKS" base
  done
}
stage_m1() {
  say "Module 1 -- variance decomposition, HEAP predictive R2 (base, lasso)"
  for i in $(seq 1 "$N_CHUNKS"); do
    run "module1_$i" Rscript scripts/module1_variance_decomposition/Module1_suggested.R \
      "$i" "$N_CHUNKS" base lasso
  done
  greml
}
greml() {   # GREML variance components with GCTA (needs gcta64 + plink2)
  local gcta plink2
  gcta="${GCTA_BIN:-$(command -v gcta64 || true)}"; plink2="${PLINK2_BIN:-$(command -v plink2 || true)}"
  if [ ! -x "$gcta" ] || [ ! -x "$plink2" ]; then
    say "Module 1 -- GREML skipped: set GCTA_BIN and PLINK2_BIN (or put gcta64/plink2 on PATH)"
    return 0
  fi
  say "Module 1 -- variance decomposition, GREML (GCTA; 4,000 people, 50,000 SNPs)"
  local G="$DEMO/greml"; mkdir -p "$G/scratch"
  (
    export HEAP_IGLOO_ROOT="$G/igloo" HEAP_SCRATCH_ROOT="$G/scratch" \
           HEAP_PROJECT_ROOT="$G/igloo/UKB/HEAP" HEAP_LOADER_RDS="$G/HEAP_synthetic.rds"
    mkdir -p "$HEAP_IGLOO_ROOT/GCTA" "$HEAP_IGLOO_ROOT/UKB/gwas"
    ln -sf "$gcta" "$HEAP_IGLOO_ROOT/GCTA/gcta64"; ln -sf "$plink2" "$HEAP_IGLOO_ROOT/UKB/gwas/plink2"
    export PLINK2_BIN="$plink2" GCTA_BIN="$gcta"
    if [ ! -f "$HEAP_LOADER_RDS" ]; then
      run greml_heap  Rscript demo/make_synthetic_heap.R "$G" "$G/config_unused" 4000 3
      run greml_geno  Rscript demo/make_synthetic_greml.R "$G" "$HEAP_IGLOO_ROOT" 50000 3
    fi
    module() { :; }; export -f module          # the job scripts call `module load`
    printf 'protein_%s
' 1 2 3 > "$G/proteins.txt"
    export RUN_GROUP=demo CONFIG="$ROOT/demo/greml_demo_config.R" PROTEIN_SET="$G/proteins.txt" \
           PERM_OUTPUT_ROOT="$HEAP_PROJECT_ROOT/output/population_architecture" \
           SCRATCH_ROOT_BASE="$G/scratch/poparch" MIN_PROTEIN_N=500 PROFILE_PREP=false \
           INCLUDE_COVAR_KERNELS=false SLURM_CPUS_PER_TASK="${SLURM_CPUS_PER_TASK:-2}"
    run greml_unrelated bash slurm/greml/BuildMasterUnrelated.sh
    for i in 1 2 3; do run "greml_protein_$i" env SLURM_ARRAY_TASK_ID="$i" bash slurm/greml/ProtGremlGrmCutoff.sh; done
    run greml_summary bash slurm/greml/SummaryGremlGrmCutoff.sh
  )
}
stage_m3() {
  say "Module 3 -- mediation of exposure -> disease through proteins (base, lasso)"
  local mf="$DEMO/module3_manifest.tsv"
  printf 'array_index\tchunk_id\tn_chunks\tcovariate_set\tfamily\tmediation_mode\texperiment_name\tsample_filter\tupstream_score_dir\toutput_path\n' > "$mf"
  for i in 1 2; do
    printf '%s\t%s\t2\tbase\tlasso\tprimary_total\tdemo_M3\tnone\t%s\t%s\n' "$i" "$i" \
      "$OUT/module1_predictive_r2_score_partition/base/lasso/mediation_scores" "$OUT/module3" >> "$mf"
  done
  for i in 1 2; do
    run "module3_$i" Rscript scripts/module3_mediation/Module3.R --manifest "$mf" --array-index "$i"
  done
}
stage_m6() {
  local x1; x1=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if($i=="exposure_x1") c=i} NR==2{print $c}' "$DEMO/planted_effects.tsv")
  say "Module 6 -- proteome-based exposure score for $x1 (base)"
  run module6 Rscript scripts/module6_prediction/Module6_prod_longitudinal.R base "$x1"
}
stage_mr() {
  say "MR -- Module5.R on synthetic GWAS summary statistics (UKB and deCODE pQTL arms), then tiering and triad motifs"
  rm -rf "$OUT"/mr "$OUT"/mr_edges "$OUT"/mr_edges_decode   # Module5.R skips edges that already have results
  run mr_inputs Rscript demo/make_synthetic_gwas.R "$HEAP_IGLOO_ROOT" "$HEAP_PROJECT_ROOT"
  for e in EP PE PD DP ED DE; do
    run "mr_ukb_$e"    env HEAP_MR_OUTDIR="$OUT/mr_edges/MR_UKB_primary" \
                         Rscript scripts/module5_mr/Module5.R "$e" 1 1
    run "mr_decode_$e" env HEAP_MR_OUTDIR="$OUT/mr_edges_decode/MR_deCODE_replication" \
                         Rscript scripts/module5_mr/Module5_deCODE.R "$e" 1 1
  done
  run mr_tables_ukb    env CPUS=1 Rscript scripts/support/mr_tables/build_mr_tables.R UKB
  run mr_tables_decode env CPUS=1 Rscript scripts/support/mr_tables/build_mr_tables.R DECODE
  run mr_arms          Rscript scripts/support/mr_tables/compare_arms.R
  run mr_triads        Rscript scripts/analysis_summaries/summarize_mr_triads.R
}

t0=$(date +%s)
if [ ! -f "$HEAP_LOADER_RDS" ]; then
  say "Synthetic data -> $DEMO"
  Rscript demo/make_synthetic_heap.R "$DEMO" "$ROOT/config" "$N_PARTICIPANTS" 6
fi
for s in "${ALL_STAGES[@]}"; do
  for want in "${STAGES[@]}"; do [ "$s" = "$want" ] && "stage_$s"; done
done
say "Summary"
Rscript demo/summarize_demo.R "$DEMO"
say "Done in $(( $(date +%s) - t0 )) s. Outputs: $OUT"
