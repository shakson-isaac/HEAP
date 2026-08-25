#!/bin/bash
#SBATCH -J m1_catgsea
#SBATCH -t 0-02:00
#SBATCH -c 2
#SBATCH --mem=24G
#SBATCH -p short
umask 002
module load gcc/14.2.0 R/4.4.2
export R_LIBS=/n/groups/patel/shakson_ukb/Rlib/gsea_env   # group-shared GSEA stack
HEAP_ROOT=/n/groups/patel/shakson_ukb/HEAP
export HEAP_PATHS_FILE="$HEAP_ROOT/workflow/00_paths.R"
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_NUM_THREADS=1
S=$HEAP_ROOT/scripts
# per-category GSEA on Module-1 exposomic R2 rankings -> per-category biology figure
echo "[$(date '+%H:%M:%S')] per-category GSEA ..." && Rscript $S/module4_enrichment/run_module1_category_gsea.R base lasso && \
echo "[$(date '+%H:%M:%S')] category-biology figure ..." && Rscript $S/visualizations/figures/fig_category_biology.R
echo "[$(date '+%H:%M:%S')] DONE exit=$?"
