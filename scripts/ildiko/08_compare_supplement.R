# =============================================================================
# 08_compare_supplement.R
# BINF6430 Final Project, Group 1: reproduction of Bloom et al. (2026), mouse arm
# Author: Ildiko Polyak
# Date:   2026-11-13
#
# Purpose
#   Compare our results with the authors' supplemental tables
#   (insight.jci.org/articles/view/207270):
#     S1 = mouse GSEA results (50 pathways x 6 comparisons, NES + padj)
#     S2 = mouse DESeq2 output (6 sheets)
#     S3 = "normalized counts" (actually SV1-adjusted, back-transformed VST values)
#   Checks:
#     A. Which samples the published analysis used (S3 libraries vs GEO and GitHub files)
#     B. Paper's DEG counts vs S2 (the "significant" column vs padj < 0.05)
#     C. Authors' published GSEA code run on their own S2 (log2FC ranking), and S2
#        ranked by stat, vs S1; plus our 30- and 33-sample GSEA vs S1
#     D. S1 vs S2-stat GSEA for all 6 comparisons (6 x 6 NES correlation matrix)
#     E. Pathway-by-pathway agreement for the 26 Fig 1D pathways (S1 vs our 33-sample stat)
#     F. Our 30-sample DESeq2 vs S2
#     G. Our 33-sample DESeq2 vs S2, with our SV1 and without SV1
#     H. Quantification check: our 33-sample baseMean vs S2 baseMean
#
# Environment
#   R 4.2.3, Bioconductor 3.15: DESeq2 1.36.0, fgsea 1.22.0, readxl 1.4.3
#
# Inputs
#   Phase3/inputs/: supplemental table 1/2/3.xlsx, samplekey.txt,
#                   sample_mapping_GSE288245.tsv, mh.all.v2023.1.Mm.symbols.gmt
#   Phase3/DESeq2_results/, DESeq2_results_33/ (DE tables, dds33.rds)
#   Phase3/fgsea_results/, fgsea_results_33/ (GSEA tables)
#   Run 05, 06, and 07 first.
#
# Outputs (Phase3/supplement_comparisons/, folder must already exist)
#   GSEA_vs_S1.tsv, GSEA_S1_vs_S2stat_matrix.tsv, Fig1D_pathways_vs_S1.tsv,
#   DE30_vs_S2.tsv, DE33_vs_S2.tsv, baseMean_vs_S2.tsv
# =============================================================================


# ---- 1. Paths and packages --------------------------------------------------
phase3 <- "C:/Users/ipoly/OneDrive/Desktop/BINF6430_Trans/Final Project/Phase3"
inputs <- file.path(phase3, "inputs")
d30    <- file.path(phase3, "DESeq2_results")
d33    <- file.path(phase3, "DESeq2_results_33")
sc     <- file.path(phase3, "supplement_comparisons")

library(readxl)
library(fgsea)
library(DESeq2)

pathways    <- gmtPathways(file.path(inputs, "mh.all.v2023.1.Mm.symbols.gmt"))
comparisons <- c("PH_14d", "PH_4m", "POA_14d", "POA_4m")
strip_ver   <- function(x) sub("\\.\\d+$", "", x)
wt <- function(x, f) write.table(x, f, sep = "\t", quote = FALSE, row.names = FALSE)

# Reproducibility check: compare numeric columns with a previously saved file
check_prev <- function(new, file) {
  if (!file.exists(file)) { cat(basename(file), "- new file\n"); return(invisible()) }
  old <- read.delim(file, check.names = FALSE)
  cols <- intersect(names(new)[sapply(new, is.numeric)], names(old))
  d <- max(sapply(cols, function(cl) max(abs(as.numeric(old[[cl]]) - new[[cl]]), na.rm = TRUE)))
  cat(basename(file), "- max difference from saved (shared columns):", signif(d, 3), "\n")
}


# ---- 2. Read the authors' tables --------------------------------------------
s2file <- file.path(inputs, "supplemental table 2.xlsx")
sheets6 <- c(PH_14d = "PH14dOVXvsPHintact", PH_4m = "PH4mOVXvsPHintact",
             POA_14d = "POA14dOVXvsPOAintact", POA_4m = "POA4mOVXvsPOAintact",
             PH_4m_vs_14d = "PH4mOVXvsPH14dOVX", POA_4m_vs_14d = "POA4mOVXvsPOA14dOVX")
s2_list <- lapply(sheets6, function(s) as.data.frame(read_excel(s2file, sheet = s)))

# S1 pathway names contain typos; match to the GMT names after removing spaces and
# punctuation. "mTORC2 signaling" is not a Hallmark set; mapped to MTORC1 (assumption).
# The MSigDB GMT itself spells REACTIVE_OXIGEN_SPECIES_PATHWAY.
s1 <- as.data.frame(read_excel(file.path(inputs, "supplemental table 1.xlsx")))
s1_key <- function(x) toupper(gsub("[^A-Za-z0-9]", "", x))
fix <- c(ANGIOGENSIS = "ANGIOGENESIS", BILDACIDMETABOLISM = "BILEACIDMETABOLISM",
         INFLAMMATORYREPONSE = "INFLAMMATORYRESPONSE", TNFALPHASIGNALINGVIANKB = "TNFASIGNALINGVIANFKB",
         MTORC2SIGNALING = "MTORC1SIGNALING", P53PATHWAYS = "P53PATHWAY",
         PI3KATKMTORSIGNALING = "PI3KAKTMTORSIGNALING", REACTIVEOXYGENSPECIES = "REACTIVEOXIGENSPECIESPATHWAY",
         KRASSIGNALINGDOWN = "KRASSIGNALINGDN", UVRESPONSEDOWN = "UVRESPONSEDN")
k <- s1_key(s1$Pathway); k <- ifelse(k %in% names(fix), fix[k], k)
gmt_key <- setNames(names(pathways), s1_key(sub("HALLMARK_", "", names(pathways))))
s1$hallmark <- gmt_key[k]
stopifnot(!any(is.na(s1$hallmark)))
s1_comp6 <- c(PH_14d = "PH 14d vs intact", PH_4m = "PH 4m vs intact",
              POA_14d = "POA 14d vs  intact", POA_4m = "POA 4m vs intact",
              PH_4m_vs_14d = "PH 4m vs 14d", POA_4m_vs_14d = "POA 4m vs 14d")
s1_comp <- s1_comp6[comparisons]


# ---- A. Which samples did the published analysis use? -----------------------
s3raw <- suppressMessages(as.data.frame(read_excel(file.path(inputs, "supplemental table 3.xlsx"),
                                                   col_names = FALSE)))
s3_libs <- unlist(s3raw[2, -1])
map <- read.delim(file.path(inputs, "sample_mapping_GSE288245.tsv"), colClasses = "character")
geo_libs <- map$library_name[map$status %in% c("in_GEO_and_analysis", "GEO_only")]
gh_libs  <- read.delim(file.path(inputs, "samplekey.txt"))$sample
cat("\nA. S3 libraries:", length(s3_libs),
    "| identical to the 33 GEO libraries:", setequal(s3_libs, geo_libs), "\n")
cat("   In authors' GitHub samplekey but not in S3:", setdiff(gh_libs, s3_libs), "\n")
cat("   In S3 but not in authors' GitHub samplekey:", setdiff(s3_libs, gh_libs), "\n")
# 2026-11-12: S3 = the 33 GEO libraries; GitHub-only: L16_5443_S1, L16_5452_S10, L534_3_S3


# ---- B. Paper's DEG counts vs S2 --------------------------------------------
cat("\nB. Protein-coding DEGs in S2:\n")
print(t(sapply(comparisons, function(n) {
  d <- s2_list[[n]]; pc <- d$gene_type %in% "protein_coding"
  c(significant_yes = sum(pc & d$significant == "yes"),
    padj_below_0.05 = sum(pc & d$padj < 0.05, na.rm = TRUE))
})))
# 2026-11-12: significant_yes 13/82/14/70 (= paper); padj < 0.05 26/101/27/77


# ---- C. GSEA vs S1 -----------------------------------------------------------
# Authors' GSEA code (GeneSetEnrichmentAnalysis.Rmd) on their own S2, two rankings
gsea_s2 <- function(rank_col, sheet_names) {
  out <- lapply(sheet_names, function(n) {
    d <- s2_list[[n]]
    d <- d[d$gene_type %in% "protein_coding" & !is.na(d[[rank_col]]), ]
    d <- d[order(d[[rank_col]], decreasing = TRUE), ]
    ranks <- setNames(d[[rank_col]], d$gene_name)
    set.seed(6430)
    fgsea(pathways = pathways, stats = ranks, eps = 0.0, minSize = 15, maxSize = 500)
  })
  names(out) <- sheet_names
  out
}
g_s2_lfc  <- gsea_s2("log2FoldChange", comparisons)
g_s2_stat <- gsea_s2("stat", names(sheets6))

# Our GSEA results, read from 07_fgsea.R output
read_gsea <- function(dir, file) {
  tab <- read.delim(file.path(phase3, dir, file))
  split(tab, tab$comparison)[comparisons]
}
g30_lfc  <- read_gsea("fgsea_results",    "GSEA_hallmark_log2FC.tsv")
g30_stat <- read_gsea("fgsea_results",    "GSEA_hallmark_stat.tsv")
g33_lfc  <- read_gsea("fgsea_results_33", "GSEA_hallmark_log2FC.tsv")
g33_stat <- read_gsea("fgsea_results_33", "GSEA_hallmark_stat.tsv")

compare_to_s1 <- function(lst) t(sapply(comparisons, function(n) {
  a <- s1[s1$Comparison == s1_comp[[n]], ]
  b <- as.data.frame(lst[[n]])[match(a$hallmark, lst[[n]]$pathway), ]
  c(sig = sum(b$padj < 0.05), NES_cor_with_S1 = round(cor(a$`Normalized Enrichment Score`, b$NES), 3))
}))

gsea_vs_s1 <- data.frame(
  S1_sig             = sapply(comparisons, function(n) sum(s1$adjusted_pvalue[s1$Comparison == s1_comp[[n]]] < 0.05)),
  authors_code_on_S2 = compare_to_s1(g_s2_lfc),
  S2_ranked_by_stat  = compare_to_s1(g_s2_stat),
  ours30_log2FC      = compare_to_s1(g30_lfc),
  ours30_stat        = compare_to_s1(g30_stat),
  ours33_log2FC      = compare_to_s1(g33_lfc),
  ours33_stat        = compare_to_s1(g33_stat))
gsea_vs_s1 <- data.frame(comparison = rownames(gsea_vs_s1), gsea_vs_s1)
cat("\nC. GSEA vs S1:\n"); print(gsea_vs_s1, row.names = FALSE)
check_prev(gsea_vs_s1, file.path(sc, "GSEA_vs_S1.tsv"))
wt(gsea_vs_s1, file.path(sc, "GSEA_vs_S1.tsv"))
# 2026-11-12: authors' code on S2: 0 significant, NES cor 0.52/0.75/0.50/0.72;
# S2 by stat: 6/27/11/14, cor 0.15/0.97/0.90/0.93


# ---- D. S1 vs S2 ranked by stat, all 6 comparisons --------------------------
cor_matrix <- sapply(names(g_s2_stat), function(s2n) {
  sapply(names(s1_comp6), function(s1n) {
    a <- s1[s1$Comparison == s1_comp6[[s1n]], ]
    b <- g_s2_stat[[s2n]][match(a$hallmark, g_s2_stat[[s2n]]$pathway), ]
    round(cor(a$`Normalized Enrichment Score`, b$NES), 2)
  })
})
cor_df <- data.frame(S1_comparison = rownames(cor_matrix), cor_matrix)
cat("\nD. S1 (rows) vs S2 ranked by stat (columns), NES correlation:\n"); print(cor_df, row.names = FALSE)
check_prev(cor_df, file.path(sc, "GSEA_S1_vs_S2stat_matrix.tsv"))
wt(cor_df, file.path(sc, "GSEA_S1_vs_S2stat_matrix.tsv"))


# ---- E. Fig 1D pathways: S1 vs our 33-sample stat GSEA ----------------------
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
fig1d_tab <- do.call(rbind, lapply(comparisons, function(n) {
  a <- s1[s1$Comparison == s1_comp[[n]], ]
  a <- a[match(fig1d, a$hallmark), ]
  b <- g33_stat[[n]][match(fig1d, g33_stat[[n]]$pathway), ]
  sig_s1 <- a$adjusted_pvalue < 0.05
  sig_us <- b$padj < 0.05
  data.frame(comparison = n, pathway = sub("HALLMARK_", "", fig1d),
             S1_NES = round(a$`Normalized Enrichment Score`, 3), S1_padj = signif(a$adjusted_pvalue, 3),
             ours33_stat_NES = round(b$NES, 3), ours33_stat_padj = signif(b$padj, 3),
             agreement = ifelse(sig_s1 & sig_us, "both significant",
                         ifelse(sig_s1, "S1 only", ifelse(sig_us, "ours only", "neither"))),
             same_direction = sign(a$`Normalized Enrichment Score`) == sign(b$NES))
}))
cat("\nE. Fig 1D pathways (26), S1 vs our 33-sample stat GSEA:\n")
print(t(sapply(comparisons, function(n) {
  x <- fig1d_tab[fig1d_tab$comparison == n, ]
  c(S1_significant = sum(x$agreement %in% c("both significant", "S1 only")),
    both_significant = sum(x$agreement == "both significant"),
    S1_only = sum(x$agreement == "S1 only"), ours_only = sum(x$agreement == "ours only"),
    same_direction = sum(x$same_direction))
})))
print(fig1d_tab[fig1d_tab$agreement %in% c("S1 only", "ours only"), ], row.names = FALSE)
check_prev(fig1d_tab, file.path(sc, "Fig1D_pathways_vs_S1.tsv"))
wt(fig1d_tab, file.path(sc, "Fig1D_pathways_vs_S1.tsv"))
# 2026-11-13: both significant 3/18/4/10; S1 only 1/0/2/1; ours only 0/2/0/2;
# same direction 19/25/19/23 of 26


# ---- F. Our 30-sample DESeq2 vs S2 ------------------------------------------
de30 <- t(sapply(comparisons, function(n) {
  ours <- read.delim(file.path(d30, paste0("DE_", n, "_vs_intact.tsv")))
  ours$ens <- strip_ver(ours$gene_id)
  theirs <- s2_list[[n]]
  m <- merge(ours, theirs, by.x = "ens", by.y = "ensembl", suffixes = c(".ours", ".theirs"))
  ours_deg   <- m$padj.ours < 0.05 & m$gene_type.ours == "protein_coding"
  theirs_deg <- m$significant == "yes" & m$gene_type.theirs %in% "protein_coding"
  c(ours_genes = nrow(ours), theirs_genes = nrow(theirs), shared_genes = nrow(m),
    stat_cor = round(cor(m$stat.ours, m$stat.theirs, use = "complete.obs"), 3),
    lfc_spearman = round(cor(m$log2FoldChange.ours, m$log2FoldChange.theirs,
                             method = "spearman", use = "complete.obs"), 3),
    ours_DEG = sum(ours_deg, na.rm = TRUE), theirs_DEG = sum(theirs_deg, na.rm = TRUE),
    DEG_overlap = sum(ours_deg & theirs_deg, na.rm = TRUE),
    pc_type_disagree = sum(xor(m$gene_type.ours == "protein_coding",
                               m$gene_type.theirs %in% "protein_coding"), na.rm = TRUE))
}))
de30 <- data.frame(comparison = rownames(de30), de30)
cat("\nF. 30-sample DESeq2 vs S2 (theirs_DEG = S2 'significant' column):\n"); print(de30, row.names = FALSE)
check_prev(de30, file.path(sc, "DE30_vs_S2.tsv"))
wt(de30, file.path(sc, "DE30_vs_S2.tsv"))


# ---- G. Our 33-sample DESeq2 vs S2, with and without SV1 --------------------
# Without SV1: refit the saved 33-sample model with design ~ condition
dds33  <- readRDS(file.path(d33, "dds33.rds"))
dds33n <- dds33
design(dds33n) <- ~ condition
dds33n <- DESeq(dds33n)
contrasts <- list(PH_14d  = c("condition", "PH_14dOVX",  "PH_intact"),
                  PH_4m   = c("condition", "PH_4mOVX",   "PH_intact"),
                  POA_14d = c("condition", "POA_14dOVX", "POA_intact"),
                  POA_4m  = c("condition", "POA_4mOVX",  "POA_intact"))

de33 <- t(sapply(comparisons, function(n) {
  ours <- read.delim(file.path(d33, paste0("DE_", n, "_vs_intact.tsv")))
  ours$ens <- strip_ver(ours$gene_id)
  theirs <- s2_list[[n]]
  m <- merge(ours, theirs, by.x = "ens", by.y = "ensembl", suffixes = c(".ours", ".theirs"))
  ours_deg   <- m$padj.ours < 0.05 & m$gene_type.ours == "protein_coding"
  theirs_deg <- m$padj.theirs < 0.05 & m$gene_type.theirs %in% "protein_coding"
  nosv <- as.data.frame(results(dds33n, contrast = contrasts[[n]], alpha = 0.05))
  nosv$ens <- strip_ver(rownames(nosv))
  mn <- merge(nosv, theirs, by.x = "ens", by.y = "ensembl", suffixes = c(".o", ".t"))
  c(ours_pc_DEG = sum(ours$padj < 0.05 & ours$gene_type == "protein_coding", na.rm = TRUE),
    S2_pc_DEG   = sum(theirs$padj < 0.05 & theirs$gene_type %in% "protein_coding", na.rm = TRUE),
    paper       = c(PH_14d = 13, PH_4m = 82, POA_14d = 14, POA_4m = 70)[[n]],
    DEG_overlap = sum(ours_deg & theirs_deg, na.rm = TRUE),
    stat_cor    = round(cor(m$stat.ours, m$stat.theirs, use = "complete.obs"), 3),
    lfc_spearman = round(cor(m$log2FoldChange.ours, m$log2FoldChange.theirs,
                             method = "spearman", use = "complete.obs"), 3),
    stat_cor_without_SV1 = round(cor(mn$stat.o, mn$stat.t, use = "complete.obs"), 3))
}))
de33 <- data.frame(comparison = rownames(de33), de33)
cat("\nG. 33-sample DESeq2 vs S2 (both padj < 0.05), with and without SV1:\n"); print(de33, row.names = FALSE)
check_prev(de33, file.path(sc, "DE33_vs_S2.tsv"))
wt(de33, file.path(sc, "DE33_vs_S2.tsv"))
# 2026-11-12: stat_cor 0.839/0.833/0.812/0.757; without SV1 0.888/0.868/0.881/0.852


# ---- H. Quantification: our 33-sample baseMean vs S2 baseMean ---------------
# baseMean = mean of DESeq2 normalized counts; does not depend on SV1.
# (S3 is not used here: it is VST + removeBatchEffect(SV1) + 2^x, not counts.)
o <- read.delim(file.path(d33, "DE_PH_4m_vs_intact.tsv")); o$ens <- strip_ver(o$gene_id)
m <- merge(o[, c("ens", "baseMean")], s2_list[["PH_4m"]][, c("ensembl", "baseMean")],
           by.x = "ens", by.y = "ensembl", suffixes = c(".ours", ".theirs"))
e <- m$baseMean.theirs > 10
q <- quantile(m$baseMean.ours[e] / m$baseMean.theirs[e], c(0.05, 0.25, 0.5, 0.75, 0.95))
bm <- data.frame(genes_compared = nrow(m), expressed_genes = sum(e),
                 spearman_all = round(cor(m$baseMean.ours, m$baseMean.theirs, method = "spearman"), 4),
                 pearson_log2_expressed = round(cor(log2(m$baseMean.ours[e]), log2(m$baseMean.theirs[e])), 4),
                 ratio_q05 = q[[1]], ratio_q25 = q[[2]], ratio_median = q[[3]], ratio_q75 = q[[4]], ratio_q95 = q[[5]])
cat("\nH. Quantification, 33-sample baseMean vs S2:\n"); print(bm, row.names = FALSE)
check_prev(bm, file.path(sc, "baseMean_vs_S2.tsv"))
wt(bm, file.path(sc, "baseMean_vs_S2.tsv"))
# 2026-11-12: Spearman 0.9993; log2 Pearson (expressed) 0.9991; ratio median 0.999

print(list.files(sc))
