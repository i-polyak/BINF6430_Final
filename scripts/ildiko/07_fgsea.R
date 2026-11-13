# =============================================================================
# 07_fgsea.R
# BINF6430 Final Project, Group 1: reproduction of Bloom et al. (2026), mouse arm
# Author: Ildiko Polyak
# Date:   2026-11-13
#
# Purpose
#   GSEA (fgsea) on Hallmark gene sets for both DESeq2 runs:
#     - 30-sample run (05_tximport_deseq2.R)  -> fgsea_results/
#     - 33-sample run (06_deseq2_33samples.R) -> fgsea_results_33/
#   Two rankings each: log2FoldChange (the authors' published method) and the DESeq2
#   Wald statistic "stat" (sensitivity check). For the 30-sample run, the authors'
#   ranking is also repeated after baseMean cutoffs (1, 5, 10, 20; our choice, not
#   from the paper) to test whether low-count genes explain the lack of significance.
#   Heatmaps show the 26 pathways in the paper's Fig 1D.
#
# Settings, following the authors' GeneSetEnrichmentAnalysis.Rmd
#   Gene sets: mh.all.v2023.1.Mm.symbols.gmt (from the authors' repo)
#   protein_coding genes only, ranked high to low, named by gene symbol
#   fgsea(eps = 0, minSize = 15, maxSize = 500)
#   Added (not in authors' code): set.seed(6430) before each fgsea call
#   Heatmap color scale -4/0/+4 = the published Fig 1D legend (authors' repo code
#   uses 0/1.25/2.5); asterisks *, **, *** = padj < 0.05, 0.01, 0.001
#
# Environment
#   R 4.2.3, Bioconductor 3.15: fgsea 1.22.0, ComplexHeatmap 2.12.1, circlize 0.4.16
#
# Inputs
#   Phase3/inputs/mh.all.v2023.1.Mm.symbols.gmt
#   Phase3/DESeq2_results/DE_<comparison>_vs_intact.tsv      (30 samples)
#   Phase3/DESeq2_results_33/DE_<comparison>_vs_intact.tsv   (33 samples)
#
# Outputs (in each fgsea folder; folders must already exist)
#   GSEA_hallmark_log2FC.tsv, GSEA_hallmark_stat.tsv   all 50 pathways x 4 comparisons
#   GSEA_summary.tsv                                   significant pathways per ranking
#   Heatmap_GSEA_hallmark_log2FC.pdf, Heatmap_GSEA_hallmark_stat.pdf
# =============================================================================


# ---- 1. Paths and packages --------------------------------------------------
phase3 <- "C:/Users/ipoly/OneDrive/Desktop/BINF6430_Trans/Final Project/Phase3"
inputs <- file.path(phase3, "inputs")

library(fgsea)
library(ComplexHeatmap)
library(circlize)

pathways    <- gmtPathways(file.path(inputs, "mh.all.v2023.1.Mm.symbols.gmt"))
stopifnot(length(pathways) == 50)
comparisons <- c("PH_14d", "PH_4m", "POA_14d", "POA_4m")


# ---- 2. Functions -----------------------------------------------------------
# Run fgsea on one DESeq2 results folder, ranking by rank_col, optionally keeping
# only genes with baseMean >= min_baseMean.
run_fgsea <- function(de_dir, rank_col, min_baseMean = NULL) {
  out <- lapply(comparisons, function(n) {
    de <- read.delim(file.path(de_dir, paste0("DE_", n, "_vs_intact.tsv")))
    de <- de[de$gene_type == "protein_coding" & !is.na(de[[rank_col]]), ]
    if (!is.null(min_baseMean)) de <- de[de$baseMean >= min_baseMean, ]
    de <- de[order(de[[rank_col]], decreasing = TRUE), ]
    ranks <- setNames(de[[rank_col]], de$gene_name)
    set.seed(6430)
    res <- fgsea(pathways = pathways, stats = ranks, eps = 0.0, minSize = 15, maxSize = 500)
    res$comparison <- n
    res
  })
  names(out) <- comparisons
  out
}

# Stack 4 comparisons into one table; leading-edge genes as comma-separated text
# (as in the authors' script).
to_table <- function(lst) {
  tab <- data.table::rbindlist(lst)
  tab$leadingEdge <- sapply(tab$leadingEdge, paste, collapse = ",")
  data.table::setcolorder(tab, c("comparison", "pathway"))
  tab
}

sig <- function(lst) sapply(lst, function(r) sum(r$padj < 0.05))

# Reproducibility check: compare new NES with a previously saved table, if present
check_prev <- function(lst, file) {
  if (file.exists(file)) {
    old <- read.delim(file)
    new <- as.data.frame(to_table(lst))
    key_old <- paste(old$comparison, old$pathway)
    key_new <- paste(new$comparison, new$pathway)
    cat(basename(dirname(file)), basename(file), "- max NES difference from saved:",
        signif(max(abs(old$NES[match(key_new, key_old)] - new$NES)), 3), "\n")
  }
}

# Heatmap: the 26 Fig 1D pathways in the paper's order
fig1d <- c("HALLMARK_ADIPOGENESIS", "HALLMARK_ALLOGRAFT_REJECTION", "HALLMARK_ANGIOGENESIS",
           "HALLMARK_APOPTOSIS", "HALLMARK_CHOLESTEROL_HOMEOSTASIS", "HALLMARK_COAGULATION",
           "HALLMARK_COMPLEMENT", "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
           "HALLMARK_ESTROGEN_RESPONSE_EARLY", "HALLMARK_ESTROGEN_RESPONSE_LATE",
           "HALLMARK_FATTY_ACID_METABOLISM", "HALLMARK_GLYCOLYSIS", "HALLMARK_HYPOXIA",
           "HALLMARK_IL2_STAT5_SIGNALING", "HALLMARK_IL6_JAK_STAT3_SIGNALING",
           "HALLMARK_INFLAMMATORY_RESPONSE", "HALLMARK_INTERFERON_ALPHA_RESPONSE",
           "HALLMARK_INTERFERON_GAMMA_RESPONSE", "HALLMARK_KRAS_SIGNALING_UP",
           "HALLMARK_MTORC1_SIGNALING", "HALLMARK_MYOGENESIS", "HALLMARK_P53_PATHWAY",
           "HALLMARK_TGF_BETA_SIGNALING", "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
           "HALLMARK_UV_RESPONSE_UP", "HALLMARK_XENOBIOTIC_METABOLISM")
stopifnot(length(fig1d) == 26, all(fig1d %in% names(pathways)))

make_mats <- function(lst) {
  nes  <- sapply(lst, function(r) r$NES[match(fig1d, r$pathway)])
  padj <- sapply(lst, function(r) r$padj[match(fig1d, r$pathway)])
  nes[is.na(nes)]   <- 0       # as in the authors' heatmap code
  padj[is.na(padj)] <- 0.99
  rownames(nes) <- rownames(padj) <- gsub("_", " ", sub("HALLMARK_", "", fig1d))
  colnames(nes) <- colnames(padj) <- c("PH 14d", "PH 4m", "POA 14d", "POA 4m")
  list(nes = nes, padj = padj)
}

col_fun <- colorRamp2(c(-4, 0, 4), c("blue", "white", "red"))

draw_hm <- function(m, title, file) {
  ht <- Heatmap(m$nes, name = "NES", column_title = title, col = col_fun,
                column_names_rot = 45, cluster_rows = FALSE, cluster_columns = FALSE,
                row_names_gp = gpar(fontsize = 10), column_names_gp = gpar(fontsize = 10),
                cell_fun = function(j, i, x, y, w, h, fill) {
                  p <- m$padj[i, j]
                  if (p < 0.001)     grid.text("***", x, y, gp = gpar(fontsize = 10))
                  else if (p < 0.01) grid.text("**",  x, y, gp = gpar(fontsize = 10))
                  else if (p < 0.05) grid.text("*",   x, y, gp = gpar(fontsize = 10))
                })
  pdf(file); draw(ht); dev.off()
  draw(ht)
}

# Run everything for one DESeq2 run and save to out_dir
gsea_run <- function(de_dir, out_dir, title_prefix, cutoffs = NULL) {
  g_lfc  <- run_fgsea(de_dir, "log2FoldChange")
  g_stat <- run_fgsea(de_dir, "stat")

  check_prev(g_lfc,  file.path(out_dir, "GSEA_hallmark_log2FC.tsv"))
  check_prev(g_stat, file.path(out_dir, "GSEA_hallmark_stat.tsv"))

  data.table::fwrite(to_table(g_lfc),  file.path(out_dir, "GSEA_hallmark_log2FC.tsv"), sep = "\t")
  data.table::fwrite(to_table(g_stat), file.path(out_dir, "GSEA_hallmark_stat.tsv"),   sep = "\t")

  summ <- rbind(log2FC_all = sig(g_lfc), stat_all = sig(g_stat))
  if (!is.null(cutoffs)) {
    for (cut in cutoffs) {
      summ <- rbind(summ, sig(run_fgsea(de_dir, "log2FoldChange", min_baseMean = cut)))
      rownames(summ)[nrow(summ)] <- paste0("log2FC_baseMean_ge_", cut)
    }
  }
  print(summ)
  write.table(data.frame(ranking = rownames(summ), summ),
              file.path(out_dir, "GSEA_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

  draw_hm(make_mats(g_lfc),  paste0(title_prefix, ": log2FC ranking (authors' method)"),
          file.path(out_dir, "Heatmap_GSEA_hallmark_log2FC.pdf"))
  draw_hm(make_mats(g_stat), paste0(title_prefix, ": stat ranking (sensitivity)"),
          file.path(out_dir, "Heatmap_GSEA_hallmark_stat.pdf"))

  print(list.files(out_dir))
  invisible(list(log2FC = g_lfc, stat = g_stat))
}


# ---- 3. 30-sample run -------------------------------------------------------
# 2026-11-12 result: significant log2FC 0/0/0/0; stat 13/31/10/20;
# baseMean cutoffs 1/5/10/20 gave at most 0-4 per comparison
gsea30 <- gsea_run(file.path(phase3, "DESeq2_results"),
                   file.path(phase3, "fgsea_results"),
                   "Hallmark Pathways", cutoffs = c(1, 5, 10, 20))


# ---- 4. 33-sample run (primary) ---------------------------------------------
# 2026-11-13 result: significant log2FC 0/0/0/0; stat 5/27/5/17 (S1: 4/25/10/16)
gsea33 <- gsea_run(file.path(phase3, "DESeq2_results_33"),
                   file.path(phase3, "fgsea_results_33"),
                   "Hallmark (33 samples)")
