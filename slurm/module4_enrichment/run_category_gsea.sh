#!/bin/bash
#SBATCH -J m4_catgsea
#SBATCH -t 0-02:00
#SBATCH -c 2
#SBATCH --mem=20G
#SBATCH -p short
umask 002
module load gcc/14.2.0 R/4.4.2
export R_LIBS=/n/groups/patel/shakson_ukb/Rlib/gsea_env   # group-shared GSEA stack
HEAP_ROOT=/n/groups/patel/shakson_ukb/HEAP
export HEAP_PATHS_FILE="$HEAP_ROOT/workflow/00_paths.R"
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_NUM_THREADS=1
S=$HEAP_ROOT/scripts
echo "[$(date '+%H:%M:%S')] 05 category GSEA ..."  && Rscript $S/module4_enrichment/05_category_enrichment.R && \
echo "[$(date '+%H:%M:%S')] category heatmap ..."  && Rscript $S/visualizations/figures/fig_enrichment_category.R
echo "[$(date '+%H:%M:%S')] DONE exit=$?"
