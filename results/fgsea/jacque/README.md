# Hallmark GSEA (Jacque)

Made by fgsea_mouse.R: fgsea against the 50 MSigDB mouse Hallmark gene sets
(mh.all.v2023.1.Mm.symbols.gmt, from the authors' repo), protein-coding genes,
eps = 0, minSize = 15, maxSize = 500, seed 6430. Results are compared with the
authors' Table S1 (jci.insight.207270.sdt1.xlsx).

## Layout
- `<set>_<design>/GSEA_stat/` and `GSEA_lfc/`: GSEA on our DESeq2 results, genes ranked by
  the Wald statistic (`stat`) or log2FoldChange (`lfc`). Runs: matched and allGEO × sv1, none,
  batch, plus authorsSV1_sv1.
- `paperS2/<rank>/`: GSEA on the authors' own DE table (Table S2), testing their GSEA step
  independently of our DE results. Rankings: `lfc`, `stat`, `lfc_bm10` (log2FoldChange,
  genes with baseMean ≥ 10 only), `signp` (sign(log2FC) × −log10(p)).
- `GSEA_all_runs_stat.txt`: every run × contrast, `stat` ranking (collect_gsea.R)
- `figures/`: NES heatmaps, paper vs ours, sv1 design, `stat` ranking (nes_heatmap.R),
  PDF and PNG:
  - `NES_heatmap_matched_sv1_stat`
  - `NES_heatmap_allGEO_sv1_stat`
  - `NES_heatmap_authorsSV1_sv1_stat`

  Rows are the 50 Hallmark pathways, ordered by the paper's NES in PH 4m vs intact.
  Color is NES; * padj < 0.05, ** < 0.01, *** < 0.001.

## Files in each run folder
- `GSEA_vs_paper_all.txt`: every pathway × contrast, our NES and padj next to the paper's
- `GSEA_comparison_summary.txt`: per contrast, NES Spearman vs paper, % same sign,
  and the number of significant pathways (padj < 0.05) in ours, in the paper, and in both

## Notes
- Table S1 pathway names were matched to Hallmark names by edit distance. The paper's
  "mTORC2 signaling" corresponds to HALLMARK_MTORC1_SIGNALING.
- `stat` is the main ranking, because it reproduces Table S1 from Table S2; `lfc` is the
  ranking in the authors' posted code.
