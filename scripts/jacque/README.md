# Scripts (Jacque)

Reproduction of the Bloom et al. (JCI Insight) mouse OVX hypothalamus analysis,
from kallisto output through DESeq2 and Hallmark GSEA.

## Environment
- before running scripts --MUST `conda deactivate` !!!!!!!
- R 4.2.1 / Bioconductor 3.15 (DESeq2 1.36.0, tximport 1.24.0, fgsea 1.22.0), matching the paper
- Run through `run_R421.sh`, which defines `rscript_r421` (container in `final_project/sif/`,
  R library `final_project/Rlib_bioc315/`). Deactivate conda before sourcing it.
- All analysis steps use `set.seed(6430)`, the group's common seed.
- pull_down_r421.sh was the setup step to pull down the .sif file.
- Container image: `final_project/sif/`, pulled with `pull_down_r421`

## Pipeline order
| Step | Script | Submit with | Output (under /scratch/caldwell.ja/final_project/results/) |
|---|---|---|---|
| tximport | tximport_mouse.R | tximport_mouse.sh | txi.rds, unadjusted_TPMs.txt, gene_annot_M30.rds || SVA (SmartSVA) | sva_mouse.R | sva_mouse.sh | sampleTable_SVs_<set>.txt |
| Authors' SV1 | authors_sv1.R | run directly | sampleTable_SVs_authorsSV1.txt |
| DESeq2 | deseq2_mouse.R | deseq2_mouse.sh | DE_<set>_<design>/DEtable_<contrast>.txt |
| Gene-level comparison to paper | compare_to_paper.R | compare_to_paper.sh | DE_<set>_<design>/comparison_to_paper.txt |
| Hallmark GSEA | fgsea_mouse.R | fgsea_mouse.sh | DE_<set>_<design>/GSEA_<rank>/, GSEA_paperS2_<rank>/ |
| Combined GSEA table | collect_gsea.R | run directly | GSEA_all_runs_stat.txt |
| NES heatmaps | nes_heatmap.R | run directly | figures/NES_heatmap_<set>_<design>_<rank>.pdf/.png |
| Copy results into repo | copy_to_repo.sh | bash | results/{deseq2,fgsea,comparisons}/jacque/ |

Each step reads the previous step's output, so after any upstream change, rerun everything downstream.
Chain jobs with `sbatch --dependency=afterok:<jobid>`.

## Arguments
- `<set>`: `matched` (30 samples), `allGEO` (33 samples), `authorsSV1` (30 samples, authors' SV1)
- `<design>`: `sv1` (~ SV1 + condition), `none` (~ condition), `batch` (~ batch + condition)
- `<rank>`: GSEA ranking metric, `stat` (DESeq2 Wald statistic) or `lfc` (log2FoldChange);
  `paperS2` runs also use `lfc_bm10` and `signp`

## Other scripts
- mapping_rate.sh: collects kallisto mapping rates for each sample
- pull_down_r421.sh: pulls the R 4.2.1 container into `final_project/sif/` (environment setup)
- check_basemean.R: one-off test of whether the authors' samples match the GEO 33 (baseMean comparison)
- check_versions.R: records package versions to session_info.txt

