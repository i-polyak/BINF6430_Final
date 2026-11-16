# 08_heatmap.R
# BINF 6430 reproducibility project — recreate the Fig 1D-style NES heatmap
# from our own 33-sample fgsea results (stat ranking).
#
# Run as an sbatch job -- see 08_submit_heatmap.sbatch
#
# Built from the authors' own plotting code (GeneSetEnrichmentAnalysis.Rmd,
# jcbloom/hypo_menopauseNavarro), with one deliberate correction:
#   Their colorRamp2(c(0, 1.25, 2.5), c("blue", "white", "red")) only covers
#   NES = 0 to 2.5. NES is signed (can be negative), and circlize clamps
#   out-of-range values to the nearest breakpoint color -- so any negative
#   enrichment in their own heatmap would render as solid blue, visually
#   identical to NES = 0. That's a real flaw in their own figure, not
#   something to replicate. This script uses a symmetric diverging scale
#   centered at 0, sized to the actual range of our NES values, so negative
#   and positive enrichment are both visible.
#
#   Their Heatmap() call also labels row_title = "Condition" on a matrix
#   they had already transposed to pathways-as-rows -- the axis titles don't
#   match the transposed orientation in their own code. Not replicated here;
#   labels below match the actual matrix orientation.

suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
})

fgsea_results <- readRDS("~/binf6430_repro_env/fgsea_results.rds")

comparisons <- c("PH_14d_vs_intact", "PH_4m_vs_intact", "POA_14d_vs_intact", "POA_4m_vs_intact")
comparison_labels <- c("PH 14d", "PH 4m", "POA 14d", "POA 4m")  # for column headers

# Build the NES and padj matrices: rows = pathways, columns = comparisons.
# Using the STAT ranking results, since that's the one that approximates
# the authors' own published figure (see 06_fgsea.R notes: their own
# log2FC-ranked code does not reproduce their own figure).
pathways <- fgsea_results[[1]]$stat$pathway
nes_mat  <- matrix(NA, nrow = length(pathways), ncol = length(comparisons),
                    dimnames = list(pathways, comparison_labels))
padj_mat <- matrix(NA, nrow = length(pathways), ncol = length(comparisons),
                    dimnames = list(pathways, comparison_labels))

for (i in seq_along(comparisons)) {
  res <- fgsea_results[[comparisons[i]]]$stat
  res <- res[match(pathways, res$pathway), ]
  nes_mat[, i]  <- res$NES
  padj_mat[, i] <- res$padj
}

# Clean up pathway names for display: drop the HALLMARK_ prefix, underscores to spaces
rownames(nes_mat)  <- gsub("_", " ", sub("^HALLMARK_", "", rownames(nes_mat)))
rownames(padj_mat) <- rownames(nes_mat)

cat(sprintf("Matrix: %d pathways x %d comparisons\n", nrow(nes_mat), ncol(nes_mat)))
cat(sprintf("NES range: %.2f to %.2f\n", min(nes_mat, na.rm = TRUE), max(nes_mat, na.rm = TRUE)))

# Symmetric diverging color scale centered at 0, sized to the actual data --
# the deliberate fix for the authors' asymmetric, negative-blind scale noted above.
max_abs_nes <- max(abs(nes_mat), na.rm = TRUE)
col_fun <- colorRamp2(c(-max_abs_nes, 0, max_abs_nes), c("blue", "white", "red"))

fontsize <- 10

pdf("~/binf6430_repro_env/Heatmap_Mouse_GSEA_Hallmark_33sample.pdf", width = 8, height = 11)
Heatmap(nes_mat,
        name = "NES",
        col = col_fun,
        column_title = "Comparison",      # matches the matrix's actual columns
        row_title = "Hallmark Pathway",   # matches the matrix's actual rows
        column_names_rot = 45,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        row_names_gp = gpar(fontsize = fontsize),
        column_names_gp = gpar(fontsize = fontsize),
        cell_fun = function(j, i, x, y, w, h, fill) {
          if (!is.na(padj_mat[i, j])) {
            if (padj_mat[i, j] < 0.001) {
              grid.text("***", x, y, gp = gpar(fontsize = fontsize))
            } else if (padj_mat[i, j] < 0.01) {
              grid.text("**", x, y, gp = gpar(fontsize = fontsize))
            } else if (padj_mat[i, j] < 0.05) {
              grid.text("*", x, y, gp = gpar(fontsize = fontsize))
            }
          }
        })
dev.off()

cat("\nSaved Heatmap_Mouse_GSEA_Hallmark_33sample.pdf\n")
