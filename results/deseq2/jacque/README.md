# DESeq2 results (Jacque)

## In this folder
`sva/`: sample tables with surrogate variables, one row per sample
- `sampleTable_SVs_matched.txt`, `sampleTable_SVs_allGEO.txt`: from sva_mouse.R (SmartSVA with
  standard model matrices ~condition vs ~1, on log2(counts + 1) of genes with total counts > 1,
  seed 6430)
- `sampleTable_SVs_authorsSV1.txt`: the authors' published SV1 (from `sampleTable_SVs.txt` in
  their GitHub repo), joined to our samples by SRR; 30 samples (authors_sv1.R)

SV1 diagnostics (R² of SV1 against batch, region, and condition) are in
`logs/jacque/sva_mouse_11875470.out`.

## Not in this folder: the full DE tables
The DE tables (7 runs × 6 contrasts, ~230 MB) are too large for GitHub. They are on Explorer at
`/scratch/caldwell.ja/final_project/results/DE_<set>_<design>/DEtable_<contrast>.txt`
and can be regenerated exactly

