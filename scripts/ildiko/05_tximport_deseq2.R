# =============================================================================
# 05_tximport_deseq2.R
# BINF6430 Final Project, Group 1: reproduction of Bloom et al. (2026), mouse arm
# Author: Ildiko Polyak
# Date:   2026-11-12
#
# Purpose
#   Import kallisto abundances (30 samples shared between GEO GSE288245 and the
#   authors' analysis), fit DESeq2 with the authors' design (~ SV1 + condition),
#   extract the 4 time point vs intact comparisons, and compare DEG counts and
#   labeled genes to the paper's Figure 2.
#
# Environment (must match to reproduce)
#   R 4.2.3, Bioconductor 3.15
#   DESeq2 1.36.0, tximport 1.24.0, rhdf5 2.40.0, rtracklayer 1.56.1
#   Full versions are written to DESeq2_results/sessionInfo.txt
#
# Inputs (Phase3/inputs/)
#   gencode.vM30.annotation.gtf.gz   GENCODE vM30 annotation (same release as kallisto index)
#   sample_mapping_GSE288245.tsv     SRR to library name to condition (built in Phase 1)
#   sampleTable_SVs.txt              authors' condition, batch, and SV1 per library
#   ExplorerAbundanceFiles/SRRxxxxxxxx/abundance.h5   kallisto 0.48.0 output, 33 samples
#
# Outputs (Phase3/DESeq2_results/)
#   DE_<comparison>_vs_intact.tsv    full results per comparison, sorted by padj
#   DEG_summary.tsv                  DEG counts vs paper
#   sample_table_30.tsv              samples used in the model
#   dds.rds                          fitted DESeq2 object (for fgsea)
#   sessionInfo.txt                  R and package versions
# =============================================================================


# ---- 1. Paths ---------------------------------------------------------------
phase3 <- "C:/Users/ipoly/OneDrive/Desktop/BINF6430_Trans/Final Project/Phase3"
inputs <- file.path(phase3, "inputs")
abund  <- file.path(phase3, "ExplorerAbundanceFiles")
outdir <- file.path(phase3, "DESeq2_results")
dir.create(outdir, showWarnings = FALSE)

library(tximport)
library(DESeq2)


# ---- 2. Sample table: the 30 samples shared by GEO and the authors' analysis --
# The authors analyzed 33 libraries. 30 of them are on GEO. The 3 missing are all
# intact controls (L16_5443_S1 PH, L16_5452_S10 POA, L534_3_S3 PH).
map <- read.delim(file.path(inputs, "sample_mapping_GSE288245.tsv"), colClasses = "character")
sv  <- read.delim(file.path(inputs, "sampleTable_SVs.txt"), row.names = 1)

shared <- map[map$status == "in_GEO_and_analysis", ]
stopifnot(nrow(shared) == 30)

coldata <- data.frame(
  SRR       = shared$SRR,
  library   = shared$library_name,
  condition = sv[shared$library_name, "condition"],
  batch     = sv[shared$library_name, "batch"],
  SV1       = sv[shared$library_name, "SV1"]   # authors' precomputed values, estimated on 33 samples
)
rownames(coldata) <- coldata$SRR
stopifnot(all(coldata$condition == shared$authors_group))

coldata$condition <- factor(coldata$condition,
  levels = c("PH_intact", "PH_14dOVX", "PH_4mOVX", "POA_intact", "POA_14dOVX", "POA_4mOVX"))
print(table(coldata$condition))   # expect 3, 5, 6, 4, 6, 6

files <- file.path(abund, coldata$SRR, "abundance.h5")
names(files) <- coldata$SRR
stopifnot(all(file.exists(files)))


# ---- 3. tx2gene and gene annotation from GENCODE vM30 ------------------------
# Same GTF the authors used. Versioned IDs are kept so they match kallisto's IDs.
gtf <- rtracklayer::import(file.path(inputs, "gencode.vM30.annotation.gtf.gz"))

tx <- gtf[gtf$type == "transcript"]
tx2gene <- data.frame(TXNAME = tx$transcript_id, GENEID = tx$gene_id)   # 149,419 transcripts

g <- gtf[gtf$type == "gene"]
geneinfo <- data.frame(gene_id = g$gene_id, gene_name = g$gene_name, gene_type = g$gene_type)
# 56,691 genes, 21,919 protein_coding. Note: the authors took gene types from
# current Ensembl via biomaRt instead, which may change the protein-coding set.

rm(gtf, tx, g)


# ---- 4. tximport (same settings as the authors) -----------------------------
txi <- tximport(files, type = "kallisto", tx2gene = tx2gene, ignoreAfterBar = TRUE)
stopifnot(all(colnames(txi$counts) == rownames(coldata)))
print(round(colSums(txi$counts) / 1e6, 1))   # counted reads per sample, millions


# ---- 5. DESeq2: one model, both regions, design ~ SV1 + condition -----------
dds <- DESeqDataSetFromTximport(txi, coldata, ~ SV1 + condition)
dds <- dds[rowSums(counts(dds)) > 1, ]   # authors' prefilter; keeps 42,166 genes
dds <- DESeq(dds)


# ---- 6. The 4 comparisons vs intact, with gene names and types --------------
contrasts <- list(
  PH_14d  = c("condition", "PH_14dOVX",  "PH_intact"),
  PH_4m   = c("condition", "PH_4mOVX",   "PH_intact"),
  POA_14d = c("condition", "POA_14dOVX", "POA_intact"),
  POA_4m  = c("condition", "POA_4mOVX",  "POA_intact")
)

res_list <- lapply(contrasts, function(ct) {
  r <- as.data.frame(results(dds, contrast = ct, alpha = 0.05))
  idx <- match(rownames(r), geneinfo$gene_id)
  r$gene_name <- geneinfo$gene_name[idx]
  r$gene_type <- geneinfo$gene_type[idx]
  r
})


# ---- 7. DEG counts vs the paper (Figure 2) ----------------------------------
deg_summary <- data.frame(
  all_DEGs       = sapply(res_list, function(r) sum(r$padj < 0.05, na.rm = TRUE)),
  protein_coding = sapply(res_list, function(r) sum(r$padj < 0.05 & r$gene_type == "protein_coding", na.rm = TRUE)),
  paper          = c(13, 82, 14, 70)
)
print(deg_summary)
# 2026-11-12 result: protein-coding 19, 56, 30, 55 vs paper 13, 82, 14, 70


# ---- 8. Paper's top gene, Gm28040 -------------------------------------------
# GENCODE vM30 has two genes named Gm28040. ENSMUSG00000115958 matches the paper
# (up in PH, down in POA); ENSMUSG00000098306 is unchanged.
gm <- do.call(rbind, lapply(names(res_list), function(n) {
  x <- res_list[[n]][which(res_list[[n]]$gene_name == "Gm28040"), c("gene_name", "log2FoldChange", "padj")]
  cbind(comparison = n, gene_id = rownames(x), x, row.names = NULL)
}))
print(gm)


# ---- 9. Genes labeled in the paper's Figure 2 volcano plots -----------------
# Names read from the figure images; all are significant in the paper.
paper_genes <- list(
  PH_14d  = c("Egr2", "Sgk1", "Gm34066", "Ulk4", "Vps13d", "Abca8a", "Eno1b", "Ltbp2", "Ccdc93"),
  PH_4m   = c("Gm34066", "Alas2", "C4b", "Serpina3n", "Btg2", "Flt1", "Aox3", "Slc4a1", "Prl", "Gm43517"),
  POA_14d = c("Pomc", "Irx3", "Tnnt1", "1700012B09Rik", "Gca", "Sgk1", "Endou", "Erich5", "Egr1", "Ermard"),
  POA_4m  = c("Btg2", "Irx3", "Cyp2e1", "Kdr", "Fbln5", "Tnnt1", "Mpeg1", "Mndal", "Zc3hav1", "Adgrg2",
              "Erich5", "Flt1", "Ctla2a", "Adm", "Lmod2", "Spred3", "Six1", "Ppp1r3a", "Shox2", "Prl")
)

gene_check <- do.call(rbind, lapply(names(paper_genes), function(n) {
  r <- res_list[[n]]
  do.call(rbind, lapply(paper_genes[[n]], function(g) {
    hits <- r[which(r$gene_name == g), ]
    if (nrow(hits) == 0) return(data.frame(comparison = n, gene = g, log2FC = NA, padj = NA))
    best <- hits[order(hits$padj), ][1, ]   # if a name has 2 IDs, keep the more significant
    data.frame(comparison = n, gene = g, log2FC = round(best$log2FoldChange, 2), padj = signif(best$padj, 2))
  }))
}))
gene_check$significant <- !is.na(gene_check$padj) & gene_check$padj < 0.05
print(gene_check)
print(table(gene_check$comparison, gene_check$significant))
# 2026-11-12 result: significant PH 14d 1/9, PH 4m 7/10, POA 14d 6/10, POA 4m 16/20


# ---- 10. Save outputs -------------------------------------------------------
for (n in names(res_list)) {
  r <- res_list[[n]]
  out <- data.frame(gene_id = rownames(r), gene_name = r$gene_name, gene_type = r$gene_type,
                    r[, c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")])
  out <- out[order(out$padj), ]
  write.table(out, file.path(outdir, paste0("DE_", n, "_vs_intact.tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
}

write.table(data.frame(comparison = rownames(deg_summary), deg_summary),
            file.path(outdir, "DEG_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(coldata, file.path(outdir, "sample_table_30.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
saveRDS(dds, file.path(outdir, "dds.rds"))
writeLines(capture.output(sessionInfo()), file.path(outdir, "sessionInfo.txt"))

print(list.files(outdir))
