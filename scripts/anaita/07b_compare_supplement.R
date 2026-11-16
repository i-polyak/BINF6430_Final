# 07b_compare_supplement.R
# BINF 6430 reproducibility project — compare our 33-sample results directly
# against the authors' published supplemental tables S1 (GSEA) and S2 (DESeq2),
# not just the headline counts quoted in Fig 1D / Fig 2.
#
# Run as an sbatch job -- see 07b_submit_compare.sbatch

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
})

## =========================================================================
## PART A: GSEA comparison (our fgsea_results.rds vs. Supplemental Table S1)
## =========================================================================

fgsea_results <- readRDS("~/binf6430_repro_env/fgsea_results.rds")

s1 <- read_excel("~/binf6430_repro_env/supplement_S1_GSEA.xlsx", sheet = "TableS1")
cat(sprintf("S1: %d rows, %d columns\n", nrow(s1), ncol(s1)))

# IMPORTANT: the authors' own Comparison column has inconsistent whitespace
# -- e.g. "POA 14d vs  intact" with a double space, confirmed by inspection.
# Collapse repeated whitespace on both sides of the match so a literal typo
# in their file doesn't silently break (or worse, silently fail to match
# without complaint) our comparison.
s1$Comparison_norm <- gsub("\\s+", " ", trimws(s1$Comparison))
cat("Unique Comparison values in S1 (after whitespace normalization):\n")
print(unique(s1$Comparison_norm))

comparison_map <- c(
  PH_14d_vs_intact  = "PH 14d vs intact",
  PH_4m_vs_intact   = "PH 4m vs intact",
  POA_14d_vs_intact = "POA 14d vs intact",
  POA_4m_vs_intact  = "POA 4m vs intact"
)
stopifnot(all(comparison_map %in% s1$Comparison_norm))  # fails loudly if still wrong after normalizing

# Pathway names differ in format: our gmt uses MSigDB style (HALLMARK_ADIPOGENESIS),
# the authors' table uses title case with spaces (Adipogenesis). Normalize both
# to the same canonical form (uppercase, underscores/spaces collapsed) to match.
#
# Beyond formatting, the authors' S1 table has genuine typos in several
# pathway names -- confirmed by inspection, not guessed: ANGIOGENSIS,
# BILD_ACID_METABOLISM, INFLAMMATORY_REPONSE, P53_PATHWAYS,
# PI3K_ATK_MTOR_SIGNALING, TNF_ALPHA_SIGNALING_VIA_NK_B, and most notably
# MTORC2_SIGNALING -- there is no Hallmark "MTORC2" gene set, only MTORC1, so
# this one is a real mislabel in their published supplement, not just a
# formatting difference. These are corrected here ONLY so the numeric
# comparison can proceed; they are NOT silently fixed and forgotten -- see
# the typo_pathways list below, which is reported separately in the output
# so the mislabeling stays visible and citable rather than disappearing into
# a quietly-higher match rate.
normalize_pathway <- function(x) {
  x <- toupper(x)
  x <- sub("^HALLMARK_", "", x)
  x <- gsub("-", "", x)                    # "ESTROGEN_RESPONSE_-_EARLY" style
  x <- gsub("/", "", x)                    # "G2/M_CHECKPOINT" style
  x <- gsub("_DOWN$", "_DN", x)            # MSigDB's official short form
  x <- gsub("[ _]+", "_", x)
  x
}

# Confirmed typos in the authors' S1 file (left as their literal strings on
# the right, normalized form on the left, for the report's pitfalls section):
typo_pathways <- c(
  "ANGIOGENSIS"                  = "ANGIOGENESIS",
  "BILD_ACID_METABOLISM"         = "BILE_ACID_METABOLISM",
  "INFLAMMATORY_REPONSE"         = "INFLAMMATORY_RESPONSE",
  "P53_PATHWAYS"                 = "P53_PATHWAY",
  "PI3K_ATK_MTOR_SIGNALING"      = "PI3K_AKT_MTOR_SIGNALING",
  "TNF_ALPHA_SIGNALING_VIA_NK_B" = "TNFA_SIGNALING_VIA_NFKB",
  "MTORC2_SIGNALING"             = "MTORC1_SIGNALING"  # genuine mislabel, not just a typo -- there is no Hallmark MTORC2 set
)
cat("\nAuthors' S1 pathway-name typos being corrected for matching purposes (reported, not hidden):\n")
print(typo_pathways)

normalize_pathway_s1 <- function(x) {
  x_norm <- normalize_pathway(x)
  ifelse(x_norm %in% names(typo_pathways), typo_pathways[x_norm], x_norm)
}

gsea_comparison_results <- list()

for (our_name in names(comparison_map)) {
  their_name <- comparison_map[[our_name]]

  ours <- fgsea_results[[our_name]]$stat  # stat ranking: the one that approximates their figure
  ours_df <- data.frame(
    pathway_norm = normalize_pathway(ours$pathway),
    our_NES  = ours$NES,
    our_padj = ours$padj
  )

  theirs <- s1 %>%
    filter(Comparison_norm == their_name) %>%
    transmute(
      pathway_norm = normalize_pathway_s1(Pathway),  # typo-corrected, for matching
      pathway_raw  = Pathway,                        # kept as-written, for reporting
      their_NES  = `Normalized Enrichment Score`,
      their_padj = adjusted_pvalue
    )

  merged <- inner_join(ours_df, theirs, by = "pathway_norm")
  cat(sprintf("\n%s: %d of our %d pathways matched to %d of their %d pathways\n",
              our_name, nrow(merged), nrow(ours_df), nrow(merged), nrow(theirs)))

  if (nrow(merged) < nrow(ours_df) || nrow(merged) < nrow(theirs)) {
    cat("  NOTE: not all pathways matched -- check for naming mismatches:\n")
    cat("  Ours not matched:   ", paste(setdiff(ours_df$pathway_norm, merged$pathway_norm), collapse = ", "), "\n")
    cat("  Theirs not matched: ", paste(setdiff(theirs$pathway_norm, merged$pathway_norm), collapse = ", "), "\n")
  }

  nes_cor <- cor(merged$our_NES, merged$their_NES)
  our_sig   <- merged$our_padj   < 0.05
  their_sig <- merged$their_padj < 0.05
  agree_table <- table(ours_significant = our_sig, theirs_significant = their_sig)

  cat(sprintf("  NES correlation: %.4f\n", nes_cor))
  cat("  Agreement on significance (padj<0.05):\n")
  print(agree_table)

  gsea_comparison_results[[our_name]] <- list(merged = merged, nes_cor = nes_cor, agree_table = agree_table)
}

## =========================================================================
## PART B: DESeq2 comparison (our res_list.rds vs. Supplemental Table S2)
## =========================================================================

res_list <- readRDS("~/binf6430_repro_env/res_list.rds")

s2_sheet_map <- c(
  PH_14d_vs_intact  = "PH14dOVXvsPHintact",
  PH_4m_vs_intact   = "PH4mOVXvsPHintact",
  POA_14d_vs_intact = "POA14dOVXvsPOAintact",
  POA_4m_vs_intact  = "POA4mOVXvsPOAintact"
)

deseq2_comparison_results <- list()

for (our_name in names(s2_sheet_map)) {
  sheet <- s2_sheet_map[[our_name]]
  theirs <- read_excel("~/binf6430_repro_env/supplement_S2_DESeq2.xlsx", sheet = sheet)

  # Check whether 'ensembl' carries a version suffix like our rownames do --
  # printed once so this can be caught rather than silently failing to join.
  if (our_name == names(s2_sheet_map)[1]) {
    cat("\nExample 'ensembl' values from S2 (checking for version suffix):\n")
    print(head(theirs$ensembl))
    cat("Example rownames from our res_list (same check):\n")
    print(head(rownames(res_list[[our_name]])))
  }

  ours <- as.data.frame(res_list[[our_name]])
  ours$ensembl_bare <- sub("\\..*", "", rownames(ours))  # strip version, just in case
  theirs$ensembl_bare <- sub("\\..*", "", theirs$ensembl)

  merged <- merge(
    ours[, c("ensembl_bare", "log2FoldChange", "stat", "padj")],
    theirs[, c("ensembl_bare", "log2FoldChange", "stat", "padj", "gene_type")],
    by = "ensembl_bare", suffixes = c("_ours", "_theirs")
  )
  cat(sprintf("\n%s: %d genes matched (ours: %d genes, theirs: %d genes)\n",
              our_name, nrow(merged), nrow(ours), nrow(theirs)))

  stat_cor   <- cor(merged$stat_ours, merged$stat_theirs, use = "complete.obs")
  log2fc_cor <- cor(merged$log2FoldChange_ours, merged$log2FoldChange_theirs, use = "complete.obs")

  our_sig   <- merged$padj_ours   < 0.05
  their_sig <- merged$padj_theirs < 0.05
  agree_table <- table(ours_significant = our_sig, theirs_significant = their_sig)

  cat(sprintf("  stat correlation:   %.4f\n", stat_cor))
  cat(sprintf("  log2FC correlation: %.4f\n", log2fc_cor))
  cat("  Agreement on significance (padj<0.05):\n")
  print(agree_table)

  deseq2_comparison_results[[our_name]] <- list(
    merged = merged, stat_cor = stat_cor, log2fc_cor = log2fc_cor, agree_table = agree_table
  )
}

## --- save everything ---
saveRDS(gsea_comparison_results, "~/binf6430_repro_env/gsea_vs_supplement.rds")
saveRDS(deseq2_comparison_results, "~/binf6430_repro_env/deseq2_vs_supplement.rds")
cat("\nSaved gsea_vs_supplement.rds and deseq2_vs_supplement.rds.\n")
