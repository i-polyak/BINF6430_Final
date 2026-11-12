# Comparisons with the authors' supplemental tables (Ildiko)

Authors' supplement: insight.jci.org/articles/view/207270 (S1 = mouse GSEA, S2 = mouse DESeq2, S3 = "normalized counts").

Files
GSEA_vs_S1.tsv                  significant pathways and NES correlation with S1 for: authors' code on their S2 (log2FC),
                                their S2 ranked by stat, and our 30-sample log2FC and stat rankings
GSEA_S1_vs_S2stat_matrix.tsv    6 x 6 NES correlations, S1 comparisons (rows) vs S2 stat-ranked (columns)
DE30_vs_S2.tsv                  our 30-sample DESeq2 vs S2: shared genes, stat and log2FC correlation, DEG overlap
DE33_vs_S2.tsv                  our 33-sample DESeq2 vs S2, with and without our SV1
baseMean_vs_S2.tsv              quantification check: our 33-sample baseMean vs S2 baseMean

Key findings
- Authors' published GSEA code (raw log2FC ranking) on their own S2 gives 0 significant pathways and
  NES correlation 0.50-0.75 with S1; ranking S2 by stat reproduces S1 closely for 4 of 6 comparisons.
- S1 "PH 14d vs intact" matches no S2 comparison.
- Quantification reproduced: baseMean correlation 0.999.
- Remaining DESeq2 differences come from SV1, which cannot be reproduced from published materials.
- S3 is SV1-adjusted, back-transformed VST values, not plain normalized counts.
