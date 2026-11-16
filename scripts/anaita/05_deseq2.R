# 05_deseq2.R
# BINF 6430 reproducibility project — DESeq2 on the 33-sample GEO run
#
# Run as an sbatch job -- see 05_submit_deseq2.sbatch
#
# Comparison to Project 07 (Reagan/Salmon data), same stage:
#   - Project 07 builds its DESeqDataSet with DESeqDataSetFromTximport(),
#     which attaches tximport's average-transcript-length offset per
#     sample/gene automatically (the sidebar on "tximport offset" in that
#     workbook explains why: it protects gene-level tests against isoform
#     switching, per Soneson et al. 2015). This script does the same --
#     loading txi.rds (the FULL tximport object, offsets included) rather
#     than the already-rounded, offset-free counts_filtered.rds that 04
#     built purely for the SVA step. Reusing counts_filtered.rds here
#     instead would silently drop that correction -- worth flagging if
#     you ever see this done differently elsewhere in the group's scripts.
#   - Project 07's filter is rowSums(counts(dds) >= 10) >= <smallest group
#     size>. The filter used in 04 (rowSums(counts) > 1) was deliberately
#     weak, because it exists only to match the authors' literal SVA
#     prefiltering code -- it is NOT meant to double as the DESeq2
#     filtering step. This script applies Project 07's stronger,
#     depth-aware filter fresh, independently of the SVA-stage one.
#   - Project 07 fits two SEPARATE DESeqDataSet objects (BMAT, MM) because
#     the two cell types have very different transcriptomes and shouldn't
#     share one dispersion trend. This project's design keeps both brain
#     regions (PH, POA) in one combined model instead, with SV1 absorbing
#     much of the region-level difference (see 04's Stage 6 note) -- the
#     opposite design choice, worth re-reading Project 07's sidebar on
#     "one model with an interaction instead of two" when writing the
#     report's design-critique section.
#   - This project's four comparisons (PH 14d vs intact, PH 4m vs intact,
#     POA 14d vs intact, POA 4m vs intact) are NOT simple coefficients
#     against one shared baseline, because POA_intact is itself a non-
#     reference level in a single combined `condition` factor (PH_intact
#     is the baseline). Project 07's `co_cultured` is a straightforward
#     two-level Wald test against one reference; here, each POA comparison
#     needs the character-vector contrast = c("condition", "level",
#     "reference") form so DESeq2 computes the right linear combination
#     regardless of which level happens to be the design's baseline.
#
# WHAT THIS SCRIPT ACTUALLY MATCHES, vs. what it flags as a suggested
# improvement instead:
#   - PREFILTERING: the authors' methods text states DESeq2 was run "on the
#     entire transcriptome annotation," with protein-coding restriction
#     applied only to DOWNSTREAM REPORTING, not before fitting. That means
#     they did NOT apply a manual prefilter at all -- not even something as
#     basic as Project 07's rowSums(counts >= 10) >= group size. This
#     script therefore runs DESeq2 on every gene that survived tximport,
#     with no manual filter, to match that literally. (DESeq2's own
#     built-in independent filtering inside results() still runs
#     automatically either way -- that's part of the method itself, not an
#     extra filter a user adds, so it's present regardless.)
#     SUGGESTED IMPROVEMENT (Project 07): Project 07 teaches applying an
#     independent prefilter like rowSums(counts >= 10) >= group size BEFORE
#     fitting, arguing it increases power by not even considering genes
#     that could never reach significance. Worth running as a documented
#     sensitivity check alongside the authors'-method run, not as the
#     primary number.
#   - SHRINKAGE: the authors' methods mention no shrinkage method at all
#     (no lfcShrink, apeglm, ashr, or "normal") -- they most likely report
#     raw Wald-test log2FoldChange directly from results(). This script's
#     primary reported numbers (DEG counts, saved res_list) use the
#     UNSHRUNKEN res object, matching that.
#     SUGGESTED IMPROVEMENT (Project 07 + Project 08): both modules
#     recommend shrinking fold changes before ranking or plotting, since
#     raw estimates from genes with few reads are unstable (Project 07's
#     "Gene B, 5 reads" example; Project 08 makes the same point for
#     choosing a GSEA ranking metric). This script still computes
#     res_shrunk_list as that suggested addition, clearly separated from
#     the authors'-method primary results below -- not used for DEG counts
#     or saved as if it were the authors' number.

suppressPackageStartupMessages({
  library(DESeq2)
  library(Biostrings)
  library(dplyr)
})

## --- load tximport object (with offsets) and SV1/condition table ---
txi <- readRDS("~/binf6430_repro_env/txi.rds")
coldata_final <- readRDS("~/binf6430_repro_env/coldata_final.rds")

stopifnot(all(rownames(coldata_final) == colnames(txi$counts)))

## --- build gene_id -> gene_type lookup from the reference fasta ---
# GENCODE pipe-delimited headers carry biotype as the second-to-last field,
# e.g. "...|Xkr4-201|Xkr4|3634|protein_coding|" -- same source Project 07
# gets from its separately-provided tx_table$description column.
ref_fasta <- "~/binf6430_backup/final_project/reference/gencode.vM30.transcripts.fa.gz"
headers <- names(readDNAStringSet(ref_fasta))
header_fields <- strsplit(headers, "\\|")

gene_info_raw <- data.frame(
  gene_id   = sapply(header_fields, `[`, 2),
  gene_name = sapply(header_fields, `[`, 6),
  gene_type = sapply(header_fields, function(x) x[length(x)])  # last field is biotype; R's strsplit drops the trailing empty string after the final "|", so this is NOT length(x)-1
)

# IMPORTANT: the biotype field above is the TRANSCRIPT's biotype, not the
# gene's. A gene with multiple transcripts can have multiple different
# biotype strings (e.g. one protein_coding isoform, one NMD/retained-intron
# isoform), which means gene_id is NOT unique in gene_info_raw as built --
# merging on it downstream would silently duplicate rows for any such gene
# and inflate DEG counts that don't filter by biotype. Collapse to one row
# per gene_id here: call a gene protein_coding if ANY of its transcripts
# are, which is the standard practical convention.
n_dup_genes <- gene_info_raw %>% dplyr::count(gene_id) %>% dplyr::filter(n > 1) %>% nrow()
cat(sprintf("Genes with more than one transcript biotype (pre-collapse): %d\n", n_dup_genes))

gene_info <- gene_info_raw %>%
  group_by(gene_id, gene_name) %>%
  summarise(
    gene_type = if (any(gene_type == "protein_coding")) "protein_coding" else gene_type[1],
    .groups = "drop"
  ) %>%
  as.data.frame()

stopifnot(!any(duplicated(gene_info$gene_id)))  # confirm the collapse actually worked

cat("gene_info preview:\n")
print(head(gene_info))

## --- build DESeqDataSet from tximport (preserves length offsets) ---
dds <- DESeqDataSetFromTximport(txi, colData = coldata_final, design = ~ SV1 + condition)
cat(sprintf("\nGenes entering DESeq2 (entire transcriptome annotation, no manual prefilter, matching the authors' stated method): %d\n",
            nrow(dds)))

## --- run DESeq2 directly on all genes, no manual prefilter ---
# DESeq2's own built-in independent filtering (inside results(), based on
# mean count across samples) still applies automatically -- that's part of
# the method itself, not a user-added filter, so genes with essentially no
# signal will still show padj = NA downstream regardless of this choice.
dds <- DESeq(dds)

## --- SUGGESTED IMPROVEMENT (Project 07), NOT the primary run: ---
## An independent prefilter before fitting, which Project 07 argues
## increases power. Left here as a documented sensitivity check rather
## than deleted, so it's easy to re-enable and compare:
# min_group_size <- min(table(coldata_final$condition))
# dds_prefiltered <- dds[rowSums(counts(dds, normalized = FALSE) >= 10) >= min_group_size, ]
# dds_prefiltered <- DESeq(dds_prefiltered)
cat("\nresultsNames(dds):\n")
print(resultsNames(dds))

## --- the four comparisons, matching the paper's Fig 2 contrasts ---
comparisons <- list(
  PH_14d_vs_intact  = c("condition", "PH_14dOVX",  "PH_intact"),
  PH_4m_vs_intact   = c("condition", "PH_4mOVX",   "PH_intact"),
  POA_14d_vs_intact = c("condition", "POA_14dOVX", "POA_intact"),
  POA_4m_vs_intact  = c("condition", "POA_4mOVX",  "POA_intact")
)

res_list       <- list()
res_shrunk_list <- list()

for (name in names(comparisons)) {
  ctr <- comparisons[[name]]
  res <- results(dds, contrast = ctr, alpha = 0.05)
  # IMPORTANT: the authors' published methods text describes kallisto ->
  # tximport -> DESeq2 -> SV1 as covariate -> fgsea on preranked lists, and
  # says NOTHING about fold-change shrinkage at all -- no lfcShrink, no
  # apeglm, no ashr, no "normal". So neither shrinkage option here is "what
  # the authors do"; they most likely report raw, unshrunken Wald-test
  # log2FoldChange values directly from results(). This shrunken column is
  # therefore a deliberate addition on top of the authors' literal method,
  # not a reproduction of it -- kept because Project 07/Project 08 both
  # recommend it for honest ranking/plotting with this few replicates, but
  # it should NOT be presented as matching the authors' approach, and it
  # has no effect on the padj-based DEG counts above either way (shrinkage
  # only changes log2FoldChange, not padj).
  #
  # apeglm needs coef, not contrast, so it can't be used for the POA
  # comparisons here (POA_intact isn't the model's baseline). ashr accepts
  # contrast but needs an extra CRAN install the rest of the group doesn't
  # have. Using type="normal": DESeq2's original built-in shrinkage, no
  # extra install, accepts contrast directly -- the most practical choice
  # among the three, not the most "authentic" one, since none of them are.
  res_shrunk <- lfcShrink(dds, contrast = ctr, type = "normal", res = res)

  res_list[[name]]        <- res
  res_shrunk_list[[name]] <- res_shrunk

  res_df <- as.data.frame(res) 
  res_df$gene_id <- rownames(res_df)
  res_df <- merge(res_df, gene_info, by = "gene_id", all.x = TRUE)

  n_total  <- sum(res_df$padj < 0.05, na.rm = TRUE)
  n_coding <- sum(res_df$padj < 0.05 & res_df$gene_type == "protein_coding", na.rm = TRUE)

  cat(sprintf("\n%s: %d DEGs total (padj<0.05), %d protein-coding\n",
              name, n_total, n_coding))
}

## --- save everything for the fgsea stage ---
saveRDS(dds, "~/binf6430_repro_env/dds.rds")
saveRDS(res_list, "~/binf6430_repro_env/res_list.rds")
saveRDS(res_shrunk_list, "~/binf6430_repro_env/res_shrunk_list.rds")
saveRDS(gene_info, "~/binf6430_repro_env/gene_info.rds")

cat("\nSaved dds.rds, res_list.rds, res_shrunk_list.rds, gene_info.rds — ready for fgsea.\n")
