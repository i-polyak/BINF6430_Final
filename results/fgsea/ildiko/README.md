# fgsea results, 30 samples (Ildiko)

Input: DESeq2 results in results/deseq2/ildiko/ (30-sample run).
Gene sets: authors' mh.all.v2023.1.Mm.symbols.gmt (MSigDB mouse Hallmark 2023.1, from their repo).
Settings as in the authors' GeneSetEnrichmentAnalysis.Rmd: protein_coding genes only,
fgsea(eps = 0, minSize = 15, maxSize = 500). Seed 6430 added (not in authors' code).
fgsea 1.22.0, ComplexHeatmap 2.12.1.

Two rankings:
- log2FC: the authors' published method (raw log2FoldChange)
- stat: sensitivity check (DESeq2 Wald statistic)

Files
GSEA_hallmark_log2FC.tsv            all 50 pathways x 4 comparisons, authors' ranking
GSEA_hallmark_stat.tsv              same, stat ranking
GSEA_summary.tsv                    significant pathways per comparison for each ranking (incl. baseMean cutoff tests)
Heatmap_GSEA_hallmark_log2FC.pdf    26 Fig 1D pathways, color scale -4/0/+4 (published legend)
Heatmap_GSEA_hallmark_stat.pdf      same, stat ranking
