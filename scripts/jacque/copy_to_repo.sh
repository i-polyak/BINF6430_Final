#!/bin/bash
# Copy final (seed 6430) results into the group repo
set -e
SRC=/scratch/caldwell.ja/final_project/results
REPO=/scratch/caldwell.ja/final_project/src/BINF6430_Final/results

# SV1 tables
mkdir -p $REPO/deseq2/jacque/sva
cp $SRC/sampleTable_SVs_*.txt $REPO/deseq2/jacque/sva/

# DESeq2, gene-level comparison, and Hallmark GSEA for each run
for run in matched_sv1 matched_none matched_batch allGEO_sv1 allGEO_none allGEO_batch authorsSV1_sv1; do
  d=$SRC/DE_$run
  mkdir -p $REPO/deseq2/jacque/$run $REPO/comparisons/jacque/$run
  cp $d/DEtable_*.txt $REPO/deseq2/jacque/$run/
  cp $d/comparison_to_paper.txt $REPO/comparisons/jacque/$run/
  for r in stat lfc; do
    mkdir -p $REPO/fgsea/jacque/$run/GSEA_$r
    cp $d/GSEA_$r/GSEA_comparison_summary.txt $d/GSEA_$r/GSEA_vs_paper_all.txt $REPO/fgsea/jacque/$run/GSEA_$r/
  done
done

# GSEA on the authors' Table S2, four rankings
for r in lfc stat lfc_bm10 signp; do
  mkdir -p $REPO/fgsea/jacque/paperS2/$r
  cp $SRC/GSEA_paperS2_$r/GSEA_comparison_summary.txt $SRC/GSEA_paperS2_$r/GSEA_vs_paper_all.txt $REPO/fgsea/jacque/paperS2/$r/
done

# Combined Hallmark table, if collect_gsea.R has been run
if [ -f $SRC/GSEA_all_runs_stat.txt ]; then cp $SRC/GSEA_all_runs_stat.txt $REPO/fgsea/jacque/; fi

# Heatmaps (PDF and PNG)
mkdir -p $REPO/fgsea/jacque/figures
cp $SRC/figures/NES_heatmap_* $REPO/fgsea/jacque/figures/

du -sh $REPO/*/jacque

