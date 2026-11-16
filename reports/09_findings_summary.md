# BINF 6430 Reproducibility Project — Findings Summary (Anaita Lashley)
**Scope:** 33-sample GEO primary reproduction, built entirely on Explorer, R 4.2.1 / Bioconductor 3.15.

---

## 1. Environment reproducibility — a genuine, citable finding on its own

Building an exact R 4.2.1 / Bioconductor 3.15 environment on Explorer (no laptop involved) surfaced three separate instances of the same underlying problem: **a package's "latest" release is not the same as the release that matches a pinned older Bioconductor stack**, and installing blindly silently breaks reproducibility rather than erroring cleanly.

- **`sva` contamination:** `BiocManager::install()` silently resolved to Bioconductor 3.16 at one point mid-session, installing `sva` and ~20 other packages at the wrong version even though `BiocManager::version()` still reported 3.15. Caught only because `renv::snapshot()` flagged a pre-flight validation failure. Fixed by force-reinstalling the full dependency list pinned explicitly to 3.15.
- **`ashr` dropped, never used:** considered for fold-change shrinkage, required a CRAN install the rest of the group didn't have; removed in favor of `type="normal"` to keep comparability across the group.
- **`GlobalOptions` packaging bug:** the current CRAN release (0.1.4, published very recently) references `utils::findMatches`, a function that has never actually been exported by base R's `utils` — a genuine bug in the package itself, not a version-compatibility issue in the usual sense. Fixed by installing the older, stable 0.1.2 release from CRAN's archive instead of "latest."

**Takeaway for the report:** three separate times, the default "just install the current version" instinct would have either silently broken reproducibility (sva) or failed outright (GlobalOptions). This is worth a paragraph in the methods/pitfalls section as its own observation about environment reproducibility being harder in practice than version-pinning alone suggests.

## 2. Reconstructing the kallisto index/quant pipeline's missing documentation

The original handoff log described building the kallisto index and running quantification, but only the quant step existed as a tracked script (`kallisto_quant.sbatch`); the index-building command had only ever been run interactively and was never captured. Reconstructed as `01_build_kallisto_index.sbatch` for the record — worth noting in the report that the reproducibility exercise caught a gap in its own documentation, not just the original paper's.

## 3. The tximport / GENCODE pipe-header pitfall

kallisto's `.h5` output stores the **entire** GENCODE pipe-delimited transcript header as `target_id` (confirmed by direct `rhdf5::h5read()` inspection), but **tximport's own kallisto importer silently strips to the first pipe field and drops the version suffix internally before comparing to `tx2gene`** — a behavior not documented in tximport's own error message, which only shows the already-stripped file-side IDs. Neither a full-header `tx2gene` nor a versioned-bare-ID `tx2gene` matched; only a pipe-stripped, version-stripped `tx2gene$TXNAME` worked. This cost significant debugging time and is worth flagging for Jacqueline/Ildiko if they build tx2gene independently.

## 4. The gene-biotype collapse bug — a real, quantified artifact

GENCODE's transcript fasta headers carry **transcript**-level biotype in their final field, not gene-level biotype. A gene with multiple transcripts of different biotypes (e.g. a protein-coding isoform and a retained-intron isoform of the same gene) produces multiple rows in a naive `gene_info` table, one per distinct biotype string. Directly measured: **21,640 genes** (out of 56,691) in this dataset have more than one distinct transcript biotype. Left uncorrected, this silently inflates any DEG count that doesn't explicitly filter by biotype, since merge() duplicates rows for affected genes. Fixed by collapsing to one row per gene, calling a gene protein-coding if any of its transcripts are (the standard practical convention, since no true gene-level GENCODE annotation was used here).

## 5. Matching the authors' literal method vs. "better practice" — three deliberate, documented divergences

Two independent sources (the full methods text, and the authors' actual Rmd scripts on GitHub) were used to establish exactly what the authors did, as distinct from what Project 07/08 (the course's own teaching modules) recommend as better practice:

- **No manual DESeq2 prefilter.** Methods text states DESeq2 ran "on the entire transcriptome annotation" — not even Project 07's own recommended `rowSums(counts>=10)>=group size` filter. Matched literally; Project 07's stronger filter is flagged as a documented, not-yet-run sensitivity check.
- **No fold-change shrinkage mentioned anywhere in their methods.** Neither `apeglm`, `ashr`, nor `lfcShrink`'s "normal" method appears in their text — they almost certainly report raw, unshrunken Wald-test log2FoldChange. This project's primary DEG counts and saved results use the unshrunken object to match; a shrunken column is computed separately and explicitly labeled as an addition, not a reproduction.
- **Ranking method for GSEA.** Their own `GeneSetEnrichmentAnalysis.Rmd` ranks by **raw log2FoldChange only** — there is no "stat" ranking anywhere in their code. "Stat ranking" throughout this project is the **group's own diagnostic hypothesis** (that their published figure may have used a different ranking than their shared code), not something the authors themselves ran. This distinction is easy to blur and should stay explicit in the write-up.

## 6. The authors' own documented inconsistencies (useful, citable material)

Independently of anything in the pipeline code, direct inspection of the authors' published supplemental materials surfaced several small but real errors:

- **A genuine mislabel, not just a typo:** Supplemental Table S1 contains a pathway labeled `MTORC2_SIGNALING` — there is no Hallmark "MTORC2" gene set, only MTORC1. This is a mislabeling in their published supplement, not a formatting quirk.
- **A formatting typo in their own `Comparison` column:** `"POA 14d vs  intact"` contains a literal double space, confirmed by direct inspection — a naive exact-string join would have silently failed to match this comparison with no error.
- **Several spelling typos in their own gene-set labels:** `ANGIOGENSIS`, `BILD_ACID_METABOLISM`, `INFLAMMATORY_REPONSE`, `P53_PATHWAYS`, `PI3K_ATK_MTOR_SIGNALING`, `TNF_ALPHA_SIGNALING_VIA_NK_B`.
- **A typo in their own hosted Hallmark gmt file** (separate from the table typos above): `REACTIVE_OXIGEN_SPECIES_PATHWAY` ("OXIGEN," not "OXYGEN") — confirms the error is in the file they distribute, not something introduced downstream.
- **A plotting bug in their own GSEA heatmap script:** their `colorRamp2(c(0, 1.25, 2.5), c("blue","white","red"))` color scale has no provision for negative NES values — any negative enrichment would render as solid blue, visually indistinguishable from NES = 0, since `circlize` clamps out-of-range values to the nearest breakpoint color. Their `Heatmap()` call also labels `row_title = "Condition"` on a matrix already transposed to pathways-as-rows — the axis titles don't match the actual orientation of their own transposed matrix.

These are all worth a consolidated "errata observed in the authors' own materials" subsection — they're exactly the kind of finding a reproducibility study is meant to surface, independent of whether your own numbers match theirs.

## 7. Numeric reproduction results (33-sample primary run)

**DESeq2 vs. Supplemental Table S2** (gene-level, all genes, not restricted to protein-coding):

| Comparison | stat correlation | log2FC correlation |
|---|---|---|
| PH 14d vs intact | 0.839 | 0.507 |
| PH 4m vs intact | 0.832 | 0.423 |
| POA 14d vs intact | 0.807 | 0.371 |
| POA 4m vs intact | 0.753 | 0.367 |

**fgsea vs. Supplemental Table S1** (pathway-level, stat ranking, after correcting both sides' naming inconsistencies):

| Comparison | NES correlation | Pathways significant in both |
|---|---|---|
| PH 14d vs intact | 0.330 | 3 of 46 matched |
| PH 4m vs intact | **0.915** | 21 of 46 matched |
| POA 14d vs intact | 0.636 | 7 of 46 matched |
| POA 4m vs intact | 0.799 | 12 of 46 matched |

**A genuinely interesting, specific observation worth its own sentence in the discussion:** PH 14d has the *highest* gene-level stat correlation (0.839) of the four comparisons, but the *lowest* GSEA NES correlation (0.330) and weakest significance overlap in both tables. This isn't a contradiction — continuous-ranking agreement (stat correlation) can be reasonably strong even while binary significance calls are unstable, which is consistent with 14 days being the weakest, earliest timepoint biologically (smallest effect sizes, most threshold-sensitive).

**Independent confirmation of the group's log2FC-ranking diagnostic:** running the authors' literal ranking method (raw log2FoldChange) on this independently-built 33-sample pipeline reproduces the same near-total collapse (1, 0, 1, 0 significant pathways) that Jacqueline found on the 30-sample run with the authors' own DESeq2 table. This is now confirmed on two independent sample sets and two independent environment builds — a much stronger claim than a single demonstration.

## 8. Seed-dependent SmartSVA crash (Jacqueline's finding, your own run as a data point)

Jacqueline found the authors' as-written SmartSVA call (with its double-wrapped intercept) crashes with `TridiagEigen: eigen decomposition failed` under seed 6430 on the 30-sample matched set, but runs under seed 1 — and gives different results between seeds even where it doesn't crash on the all-GEO set. Your own 33-sample run, using the corrected (non-double-wrapped) SmartSVA call, completed cleanly under seed 6430 with no crash. This is consistent with — though doesn't fully prove — the hypothesis that the instability is specific to the double-wrapping bug combined with the smaller 30-sample subset, rather than a general SmartSVA/seed-6430 incompatibility.

## 9. Open items / suggested follow-ups, not yet done

- `isva` version not explicitly checked against the group's comparison table
- stat correlation "without SV1" row not run (only the with-SV1 model was fit)
- Supplemental Table S3 (expression/counts) not compared — has a broken multi-row header requiring a different read strategy; set aside as lower priority
- Shrunken-log2FC ranking (Project 08's suggested third GSEA ranking method) not run, only the authors' two methods
- Project 07's independent prefilter, as a documented sensitivity check against the no-prefilter primary run, left as commented-out code, not executed
