#!/bin/bash
#SBATCH -J m4_gsea
#SBATCH -t 0-04:00
#SBATCH -c 6
#SBATCH --mem=32G
#SBATCH -p short
umask 002
module load gcc/14.2.0 R/4.4.2
# GSEA stack (clusterProfiler/org.Hs.eg.db/ReactomePA) lives in the group-shared
# lib (mirror of the orphaned ~/R-4.1.1 set; loads fine under R 4.4.2).
export R_LIBS=/n/groups/patel/shakson_ukb/Rlib/gsea_env
HEAP_ROOT=/n/groups/patel/shakson_ukb/HEAP
export HEAP_PATHS_FILE="$HEAP_ROOT/workflow/00_paths.R"
export HEAP_GSEA_CORES=6
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_NUM_THREADS=1
INPUT=/n/groups/patel/IGLOO/UKB/HEAP/output/module4_enrichment/gsea_input_M2_base_main.csv
S=$HEAP_ROOT/scripts
echo "[$(date '+%H:%M:%S')] 03 GSEA ..."        && Rscript $S/module4_enrichment/03_run_gsea.R --input "$INPUT" && \
echo "[$(date '+%H:%M:%S')] 04 flatten ..."     && Rscript $S/module4_enrichment/04_enrichment_tables.R && \
echo "[$(date '+%H:%M:%S')] tissue heatmap ..." && Rscript $S/visualizations/figures/fig_tissue_enrichment.R && \
echo "[$(date '+%H:%M:%S')] pathway heatmap ..."&& Rscript $S/visualizations/figures/fig_pathway_enrichment.R
echo "[$(date '+%H:%M:%S')] DONE exit=$?"
