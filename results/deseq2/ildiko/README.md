# DESeq2 results (Ildiko)

Produced by scripts/ildiko/05_tximport_deseq2.R
R 4.2.3, Bioconductor 3.15, DESeq2 1.36.0, tximport 1.24.0 (see sessionInfo.txt)

Samples: the 30 libraries shared by GEO GSE288245 and the authors' analysis (sample_table_30.tsv).
Design: ~ SV1 + condition, one model with both regions, authors' SV1 values.
Comparisons: each time point vs intact, within region.

Files
DE_<region>_<time>_vs_intact.tsv   all tested genes, sorted by padj (gene type from GENCODE vM30)
DEG_summary.tsv                    DEG counts at padj < 0.05 vs the paper
sample_table_30.tsv                samples, conditions, batch, SV1
sessionInfo.txt                    R and package versions

Not included: dds.rds (fitted model, recreated by the script)

## Two versions

This folder: 30-sample run, the libraries shared by GEO and the authors' GitHub sampleTable_SVs.txt,
using the authors' published SV1 values.

33_samples/: the 33 GEO libraries, which Supplemental Table S3 shows is the sample set of the
published analysis, with SV1 recomputed. This is the primary reproduction; the 30-sample run is kept
for comparison. See 33_samples/README.md and results/comparisons/ildiko/.
