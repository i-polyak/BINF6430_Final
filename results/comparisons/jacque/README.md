# Gene-level comparison to the paper (Jacque)

Made by compare_to_paper.R: compares our DESeq2 results with the authors' Table S2
(jci.insight.207270.sdt2.xlsx), contrast by contrast. Genes are joined on Ensembl ID
(version suffix removed).

## Layout
One subfolder per DESeq2 run, each with `comparison_to_paper.txt`:
- `matched_sv1`, `matched_none`, `matched_batch` (30 samples)
- `allGEO_sv1`, `allGEO_none`, `allGEO_batch` (33 samples)
- `authorsSV1_sv1` (30 samples, authors' SV1)

Contrasts (one row each): PH 14dOVX vs intact, PH 4mOVX vs intact, PH 4mOVX vs 14dOVX,
POA 14dOVX vs intact, POA 4mOVX vs intact, POA 4mOVX vs 14dOVX.

## Columns
- `ours`: our DEGs (protein-coding, padj < 0.05)
- `paper_padj`, `shared_padj`: the paper's DEGs by padj < 0.05, and how many we share
- `paper_yes`, `shared_yes`: the paper's DEGs by its `significant == "yes"` column, and how many we share
- `genes_shared`: genes present in both tables
- `baseMean_spearman`: Spearman correlation of baseMean, all shared genes
- `LFC_spearman`: Spearman correlation of log2FoldChange, all shared genes
- `LFC_pearson_bm10`: Pearson correlation of log2FoldChange, genes with paper baseMean ≥ 10
- `stat_spearman`, `stat_pearson`: Spearman and Pearson correlations of the Wald statistic, all shared genes

The DE tables these were computed from are not in the repo (see `deseq2/jacque/README.md`).
