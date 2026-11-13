# =============================================================================
# 06_deseq2_33samples.R
# BINF6430 Final Project, Group 1: reproduction of Bloom et al. (2026), mouse arm
# Author: Ildiko Polyak
# Date:   2026-11-13
#
# Purpose
#   DESeq2 on the 33 GEO samples (GSE288245), the sample set of the published
#   analysis (authors' Supplemental Table S3 lists exactly these 33 libraries; their
#   GitHub samplekey.txt / sampleTable_SVs.txt list a different 33). SV1 is
#   recomputed with SmartSVA, following the authors' human script
#   (DEanalysis_human.Rmd), because their mouse SVA code and their SV1 values for
#   these 33 samples were not published. This is the primary reproduction.
#
# Environment
#   R 4.2.3, Bioconductor 3.15: DESeq2 1.36.0, tximport 1.24.0, rtracklayer 1.56.1,
#   SmartSVA 0.1.3 (with sva 3.44.0, qvalue 2.28.0), readxl 1.4.3.
#   Full versions are written to DESeq2_results_33/sessionInfo.txt
#
# Inputs (Phase3/inputs/ unless noted)
#   sample_mapping_GSE288245.tsv     SRR, library, GEO title, status (built in Phase 1)
#   gencode.vM30.annotation.gtf.gz   GENCODE vM30 annotation
#   sampleTable_SVs.txt              authors' GitHub SV1 (used only for a sanity check)
#   supplemental table 2.xlsx        authors' DESeq2 output (for the DEG summary)
#   Phase3/ExplorerAbundanceFiles/SRRxxxxxxxx/abundance.h5   kallisto output, 33 samples
#
# Outputs (Phase3/DESeq2_results_33/, folder must already exist)
#   DE_<comparison>_vs_intact.tsv    full results per comparison, sorted by padj
#   DEG_summary.tsv                  DEG counts vs S2 and paper, overlap, agreement with S2
#   sample_table_33.tsv              samples, conditions, GEO batch, recomputed SV1
#   dds33.rds                        fitted DESeq2 object (for GSEA)
#   sessionInfo.txt                  R and package versions
# =============================================================================


# ---- 1. Paths and packages --------------------------------------------------
phase3 <- "C:/Users/ipoly/OneDrive/Desktop/BINF6430_Trans/Final Project/Phase3"
inputs <- file.path(phase3, "inputs")
abund  <- file.path(phase3, "ExplorerAbundanceFiles")
d33    <- file.path(phase3, "DESeq2_results_33")

library(tximport)
library(DESeq2)
library(SmartSVA)
library(readxl)


# ---- 2. Sample table: the 33 GEO samples ------------------------------------
# Conditions are built from the GEO titles (e.g. "PH_sham_replicate3"); sham = intact.
map <- read.delim(file.path(inputs, "sample_mapping_GSE288245.tsv"), colClasses = "character")
geo <- map[map$status %in% c("in_GEO_and_analysis", "GEO_only"), ]
stopifnot(nrow(geo) == 33)

region <- sub("_.*", "", geo$GEO_title)
tp     <- sapply(strsplit(geo$GEO_title, "_"), `[`, 2)
tp[tp == "sham"] <- "intact"
cd33 <- data.frame(SRR = geo$SRR, library = geo$library_name,
                   condition = paste0(region, "_", tp), GEO_batch = geo$GEO_batch)
rownames(cd33) <- cd33$SRR

# Check: matches the authors' labels for the 30 samples also in their GitHub files
shared_rows <- geo$status == "in_GEO_and_analysis"
stopifnot(all(cd33$condition[shared_rows] == geo$authors_group[shared_rows]))

cd33$condition <- factor(cd33$condition,
  levels = c("PH_intact", "PH_14dOVX", "PH_4mOVX", "POA_intact", "POA_14dOVX", "POA_4mOVX"))
print(table(cd33$condition))   # expect 4, 6, 6, 5, 6, 6 (same as S3)


# ---- 3. tx2gene and gene annotation from GENCODE vM30 -----------------------
gtf <- rtracklayer::import(file.path(inputs, "gencode.vM30.annotation.gtf.gz"))
tx <- gtf[gtf$type == "transcript"]
tx2gene <- data.frame(TXNAME = tx$transcript_id, GENEID = tx$gene_id)
g <- gtf[gtf$type == "gene"]
geneinfo <- data.frame(gene_id = g$gene_id, gene_name = g$gene_name, gene_type = g$gene_type)
rm(gtf, tx, g)


# ---- 4. tximport (same settings as the authors) -----------------------------
files33 <- file.path(abund, cd33$SRR, "abundance.h5")
names(files33) <- cd33$SRR
stopifnot(all(file.exists(files33)))
txi33 <- tximport(files33, type = "kallisto", tx2gene = tx2gene, ignoreAfterBar = TRUE)

dds33 <- DESeqDataSetFromTximport(txi33, cd33, ~ condition)
dds33 <- dds33[rowSums(counts(dds33)) > 1, ]   # authors' prefilter; 42,561 genes on 2026-11-12
print(nrow(dds33))


# ---- 5. SV1, following the authors' human SVA code --------------------------
# Copied from DEanalysis_human.Rmd, with condition in place of AGE:
#   - log2 of RAW counts (their code computes a size-factor-scaled matrix but does not use it)
#   - models re-wrapped with model.matrix(), duplicating the intercept, as in their code
#   - seed 6430 added (not in their code) so results are identical across the group
Y    <- as.matrix(log2(counts(dds33) + 1))
mod1 <- model.matrix(~ condition, data = cd33)
mod0 <- model.matrix(~ 1, data = cd33)
Y.r  <- t(resid(lm(t(Y) ~ mod1, data = cd33)))
n.sv <- EstDimRMT(Y.r, FALSE)$dim + 1
mod1 <- model.matrix(~ mod1, cd33)
mod0 <- model.matrix(~ mod0, cd33)
cat("Generating", n.sv, "SVs\n")   # 19 on 2026-11-12
set.seed(6430)
sv.obj <- smartsva.cpp(Y, mod = mod1, mod0 = mod0, B = 1000, alpha = 1, n.sv = n.sv)
cd33$SV1 <- sv.obj$sv[, 1]

# Sanity check vs the authors' GitHub SV1 (30 shared samples; sign is arbitrary)
sv_gh <- read.delim(file.path(inputs, "sampleTable_SVs.txt"), row.names = 1)
both  <- cd33$library %in% rownames(sv_gh)
cat("Correlation with authors' GitHub SV1:",
    round(cor(cd33$SV1[both], sv_gh[cd33$library[both], "SV1"]), 3), "\n")   # -0.91 on 2026-11-12

# Reproducibility check vs the previously saved SV1, if present
prev <- file.path(d33, "sample_table_33.tsv")
if (file.exists(prev)) {
  old <- read.delim(prev)
  cat("Max difference from previously saved SV1:",
      signif(max(abs(old$SV1[match(cd33$SRR, old$SRR)] - cd33$SV1)), 3), "\n")
}


# ---- 6. DESeq2: design ~ SV1 + condition ------------------------------------
colData(dds33)$SV1 <- cd33$SV1
design(dds33) <- ~ SV1 + condition
dds33 <- DESeq(dds33)


# ---- 7. The 4 comparisons vs intact, with gene names and types --------------
contrasts <- list(
  PH_14d  = c("condition", "PH_14dOVX",  "PH_intact"),
  PH_4m   = c("condition", "PH_4mOVX",   "PH_intact"),
  POA_14d = c("condition", "POA_14dOVX", "POA_intact"),
  POA_4m  = c("condition", "POA_4mOVX",  "POA_intact")
)
res33 <- lapply(contrasts, function(ct) {
  r <- as.data.frame(results(dds33, contrast = ct, alpha = 0.05))
  idx <- match(rownames(r), geneinfo$gene_id)
  r$gene_name <- geneinfo$gene_name[idx]
  r$gene_type <- geneinfo$gene_type[idx]
  r
})


# ---- 8. DEG summary vs the authors' S2 and the paper ------------------------
# Like-for-like: our padj < 0.05 vs their padj < 0.05, protein-coding genes.
# Paper's 13/82/14/70 use S2's "significant" column, which applies an unstated extra filter.
s2file <- file.path(inputs, "supplemental table 2.xlsx")
sheets <- c(PH_14d = "PH14dOVXvsPHintact", PH_4m = "PH4mOVXvsPHintact",
            POA_14d = "POA14dOVXvsPOAintact", POA_4m = "POA4mOVXvsPOAintact")
strip_ver <- function(x) sub("\\.\\d+$", "", x)

de33_compare <- t(sapply(names(sheets), function(n) {
  ours <- res33[[n]]
  ours$ens <- strip_ver(rownames(ours))
  theirs <- as.data.frame(read_excel(s2file, sheet = sheets[[n]]))
  m <- merge(ours, theirs, by.x = "ens", by.y = "ensembl", suffixes = c(".ours", ".theirs"))
  ours_deg   <- m$padj.ours < 0.05 & m$gene_type.ours == "protein_coding"
  theirs_deg <- m$padj.theirs < 0.05 & m$gene_type.theirs %in% "protein_coding"
  c(ours_pc_DEG  = sum(res33[[n]]$padj < 0.05 & res33[[n]]$gene_type == "protein_coding", na.rm = TRUE),
    S2_pc_DEG    = sum(theirs$padj < 0.05 & theirs$gene_type %in% "protein_coding", na.rm = TRUE),
    paper        = c(PH_14d = 13, PH_4m = 82, POA_14d = 14, POA_4m = 70)[[n]],
    DEG_overlap  = sum(ours_deg & theirs_deg, na.rm = TRUE),
    stat_cor     = round(cor(m$stat.ours, m$stat.theirs, use = "complete.obs"), 3),
    lfc_spearman = round(cor(m$log2FoldChange.ours, m$log2FoldChange.theirs,
                             method = "spearman", use = "complete.obs"), 3))
}))
print(de33_compare)
# 2026-11-12 result: ours_pc_DEG 17/61/44/65; S2 26/101/27/77; stat_cor 0.839/0.833/0.812/0.757


# ---- 9. Save outputs --------------------------------------------------------
wt <- function(x, f) write.table(x, f, sep = "\t", quote = FALSE, row.names = FALSE)

for (n in names(res33)) {
  r <- res33[[n]]
  out <- data.frame(gene_id = rownames(r), gene_name = r$gene_name, gene_type = r$gene_type,
                    r[, c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")])
  wt(out[order(out$padj), ], file.path(d33, paste0("DE_", n, "_vs_intact.tsv")))
}
wt(data.frame(comparison = rownames(de33_compare), de33_compare), file.path(d33, "DEG_summary.tsv"))
wt(cd33, file.path(d33, "sample_table_33.tsv"))
saveRDS(dds33, file.path(d33, "dds33.rds"))
writeLines(capture.output(sessionInfo()), file.path(d33, "sessionInfo.txt"))

print(list.files(d33))
