# 06_fgsea.R
# BINF 6430 reproducibility project — fgsea on the 33-sample DESeq2 results
#
# Run as an sbatch job -- see 06_submit_fgsea.sbatch
#
# Comparison to Project 08 (Reagan/clusterProfiler GSEA), same stage:
#   - Project 08 runs GSEA through clusterProfiler's gseGO()/GSEA() wrappers,
#     which call the fgsea algorithm internally for the preranked case --
#     same underlying method as calling fgsea() directly here, just a
#     different interface (gseGO() also handles GO-graph navigation, which
#     isn't relevant to a Hallmark-gene-set run like this one).
#   - Project 08 names THREE candidate ranking metrics: raw log2FC, shrunken
#     log2FC, and the Wald statistic -- and argues shrunken LFC is often the
#     better middle ground, since raw LFC lets poorly-measured genes (e.g.
#     "Gene B, 5 reads") dominate the ranking the same way Project 07 warned
#     about for individual gene calls. This script only replicates the
#     authors' own TWO methods (raw log2FC, stat) because that's what the
#     group's diagnostic specifically tested -- but Project 08's framing is
#     a reason a shrunken-LFC ranking (res_shrunk_list from 05) could be
#     worth adding as a third comparison in the report, not just the two
#     the authors used.
#   - Project 08 explicitly warns: "never pass a thresholded list to
#     gseGO(); the result looks like GSEA and means nothing." This script
#     already does it the right way -- ranking every tested gene, not just
#     the padj<0.05 subset -- but that warning is worth keeping in mind if
#     the fgsea call is ever modified later.
#   - Project 08 sets the same seed (6430) for the same stated reason:
#     "GSEA permutations are random." Confirms this project's seed choice
#     lines up with how the course itself treats GSEA's randomness, not
#     just an arbitrary number the group picked.
#   - Project 08 cautions that GSEA p-values (like ORA's) are best read as
#     "a ranking of sets... not as calibrated error rates," and notes the
#     Broad Institute's own GSEA software suggests FDR < 0.25 for
#     exploratory, sample-permutation runs -- looser than the 0.05 used
#     here and in the paper. Worth a line in the report acknowledging 0.05
#     is a stricter bar than GSEA's own documentation recommends for this
#     kind of preranked analysis.
#   - Project 08's duplicate-gene-name handling before gseGO() is a plain
#     !duplicated() call -- keeping whichever entry appears first,
#     arbitrarily. This script instead keeps the entry with the largest
#     |stat|/|log2FoldChange| for each symbol (see build_ranks() below) --
#     a deliberate tie-break rather than an arbitrary one, which is a
#     real difference worth naming if asked to justify the choice.
#
# CONFIRMED FROM THE AUTHORS' ACTUAL GSEA SCRIPT (GeneSetEnrichmentAnalysis.Rmd,
# jcbloom/hypo_menopauseNavarro), not just their methods paragraph:
#   - Their ranking is RAW log2FoldChange ONLY. There is no "stat" ranking
#     anywhere in their code. "Ranking by stat" is the GROUP'S OWN diagnostic
#     hypothesis (that their published FIGURE may have used a different
#     ranking than their shared CODE), not something the authors themselves
#     ever ran. Keep this distinction explicit in the write-up: the log2FC
#     ranking below is a literal reproduction; the stat ranking is the
#     group's own test of that hypothesis, run for comparison only.
#   - They filter to protein_coding genes BEFORE building the ranked vector
#     (gene.type == "protein_coding" on the log2FC table, before ranking) --
#     not after, and not as a separate reporting step. This script now does
#     the same for both rankings, so the stat-ranking comparison is apples-
#     to-apples with what their own code structurally does to the gene list.
#   - Their fgsea() call sets minSize = 15, maxSize = 500 explicitly. Using
#     fgsea's defaults instead (as this script did previously) silently
#     changes which pathways are even tested, not just the p-values -- fixed
#     below to match their literal call.
#   - No set.seed() anywhere in their script, consistent with the group's own
#     comparison table ("Random seed: none" for the authors). The seed set
#     here (6430) remains the GROUP'S agreed convention, not a reproduction.

suppressPackageStartupMessages({
  library(fgsea)
  library(dplyr)
})

set.seed(6430)  # per the group's Friday-meeting agreement

## --- load DESeq2 results and gene annotation from stage 5 ---
res_list       <- readRDS("~/binf6430_repro_env/res_list.rds")
gene_info      <- readRDS("~/binf6430_repro_env/gene_info.rds")

## --- load the Hallmark mouse gene set file ---
# Must be the EXACT file the authors used, not msigdbr (current msigdbr
# pulls a newer MSigDB release with different gene set contents -- the
# mistake Ildiko caught Jacqueline about). Pull this file from the authors'
# repo (jcbloom/hypo_menopauseNavarro, "Gene Set Enrichment Analysis" folder)
# and place it at the path below before running this script.
gmt_path <- "~/binf6430_repro_env/mh.all.v2023.1.Mm.symbols.gmt"
stopifnot(file.exists(gmt_path))
hallmark <- gmtPathways(gmt_path)
cat(sprintf("Loaded %d Hallmark pathways from %s\n", length(hallmark), gmt_path))

## --- build a ranking vector for one comparison, one method ---
# fgsea needs gene SYMBOLS (to match the gmt file), not Ensembl IDs, and no
# duplicate names or NAs. Multiple Ensembl IDs can map to one symbol; where
# that happens, keep the entry with the largest |stat|/|log2FoldChange| (the
# most confident call for that symbol) -- the authors' own script shows no
# explicit handling for this collision at all, so this is a documented
# assumption on our end, not a gap in matching them.
#
# protein_coding restriction added here (not just at DEG-counting time in
# 05_deseq2.R) specifically BECAUSE the authors' own script filters to
# gene.type == "protein_coding" before ranking, not after.
build_ranks <- function(res, gene_info, rank_by = c("stat", "log2FoldChange")) {
  rank_by <- match.arg(rank_by)
  df <- as.data.frame(res)
  df$gene_id <- rownames(df)
  df <- merge(df, gene_info[, c("gene_id", "gene_name", "gene_type")], by = "gene_id")
  df <- df[df$gene_type == "protein_coding", ]
  df <- df[!is.na(df[[rank_by]]) & !is.na(df$gene_name) & df$gene_name != "", ]

  df <- df %>%
    group_by(gene_name) %>%
    slice_max(order_by = abs(.data[[rank_by]]), n = 1, with_ties = FALSE) %>%
    ungroup()

  ranks <- setNames(df[[rank_by]], df$gene_name)
  sort(ranks, decreasing = TRUE)
}

## --- run fgsea for one comparison, both ranking methods ---
run_fgsea_both <- function(res, gene_info, hallmark, label) {
  ranks_stat   <- build_ranks(res, gene_info, "stat")
  ranks_log2fc <- build_ranks(res, gene_info, "log2FoldChange")

  # minSize = 15, maxSize = 500: the authors' literal parameters, from their
  # own GeneSetEnrichmentAnalysis.Rmd. Using fgsea's defaults instead changes
  # which Hallmark pathways are even tested, not just their p-values.
  fgsea_stat   <- fgsea(pathways = hallmark, stats = ranks_stat,   eps = 0, minSize = 15, maxSize = 500)
  fgsea_log2fc <- fgsea(pathways = hallmark, stats = ranks_log2fc, eps = 0, minSize = 15, maxSize = 500)

  n_sig_stat   <- sum(fgsea_stat$padj   < 0.05, na.rm = TRUE)
  n_sig_log2fc <- sum(fgsea_log2fc$padj < 0.05, na.rm = TRUE)

  cat(sprintf("\n%s:\n  stat ranking:   %d significant pathways (padj<0.05)\n  log2FC ranking: %d significant pathways (padj<0.05)\n",
              label, n_sig_stat, n_sig_log2fc))

  list(stat = fgsea_stat, log2fc = fgsea_log2fc,
       ranks_stat = ranks_stat, ranks_log2fc = ranks_log2fc)
}

## --- run all four comparisons ---
fgsea_results <- list()
for (name in names(res_list)) {
  fgsea_results[[name]] <- run_fgsea_both(res_list[[name]], gene_info, hallmark, name)
}

## --- reference point: paper's published Fig 1D counts (stat-ranking target) ---
cat("\nPaper's published significant-pathway counts (Fig 1D / S1), for comparison:\n")
cat("  PH 14d vs intact:  4\n  PH 4m vs intact:  25\n  POA 14d vs intact: 10\n  POA 4m vs intact: 16\n")

## --- save everything ---
saveRDS(fgsea_results, "~/binf6430_repro_env/fgsea_results.rds")
cat("\nSaved fgsea_results.rds.\n")
