suppressPackageStartupMessages(library(DESeq2))

args <- commandArgs(trailingOnly = TRUE)
set_name    <- if (length(args) > 0) args[1] else "matched"   # "matched" or "allGEO"
design_name <- if (length(args) > 1) args[2] else "sv1"       # "sv1", "none", or "batch"

proj   <- "/scratch/caldwell.ja/final_project"
outdir <- file.path(proj, "results", paste0("DE_", set_name, "_", design_name))
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

txi        <- readRDS(file.path(proj, "results/txi.rds"))
gene_annot <- readRDS(file.path(proj, "results/gene_annot_M30.rds"))
coldata    <- read.delim(file.path(proj, "results",
                         paste0("sampleTable_SVs_", set_name, ".txt")))
rownames(coldata)  <- coldata$sample
coldata$condition  <- factor(coldata$condition)
coldata$batch      <- factor(coldata$batch)
coldata$SV1        <- as.numeric(coldata$SV1)

# subset txi to this sample set
keep <- coldata$sample
txi$counts    <- txi$counts[, keep]
txi$abundance <- txi$abundance[, keep]
txi$length    <- txi$length[, keep]
txi$infReps   <- NULL
stopifnot(identical(colnames(txi$counts), rownames(coldata)))

design_formula <- switch(design_name,
                         sv1   = ~ SV1 + condition,
                         none  = ~ condition,
                         batch = ~ batch + condition)
cat("Sample set:", set_name, "| design:", deparse(design_formula), "|",
    ncol(txi$counts), "samples\n")

# authors' steps
dds <- DESeqDataSetFromTximport(txi, coldata, design_formula)
dds <- dds[rowSums(counts(dds)) > 1, ]
cat("Genes after filter:", nrow(dds), "\n")
dds <- DESeq(dds)
saveRDS(dds, file.path(outdir, "dds.rds"))

contrasts <- list(
  PH_14dOVX_vs_intact  = c("condition", "PH_14dOVX",  "PH_intact"),
  PH_4mOVX_vs_intact   = c("condition", "PH_4mOVX",   "PH_intact"),
  POA_14dOVX_vs_intact = c("condition", "POA_14dOVX", "POA_intact"),
  POA_4mOVX_vs_intact  = c("condition", "POA_4mOVX",  "POA_intact"),
  PH_4mOVX_vs_14dOVX   = c("condition", "PH_4mOVX",   "PH_14dOVX"),
  POA_4mOVX_vs_14dOVX  = c("condition", "POA_4mOVX",  "POA_14dOVX"))
paper_pc <- c(13, 82, 14, 70, NA, NA)   # Figure 2 protein-coding DEGs

summ <- data.frame()
for (i in seq_along(contrasts)) {
  nm  <- names(contrasts)[i]
  res <- results(dds, contrast = contrasts[[i]], alpha = 0.05)
  dat <- as.data.frame(res)

  idx <- match(rownames(dat), gene_annot$gene_id)
  dat$gene_name <- gene_annot$gene_name[idx]
  dat$gene_type <- gene_annot$gene_type[idx]
  dat$significant <- ifelse(!is.na(dat$padj) & dat$padj < 0.05, "2", "1")

  write.table(dat, file.path(outdir, paste0("DEtable_", nm, ".txt")),
              sep = "\t", quote = FALSE, col.names = NA)

  sig <- dat$significant == "2"
  pc  <- dat$gene_type == "protein_coding"
  summ <- rbind(summ, data.frame(contrast = nm,
                                 DEGs_all = sum(sig),
                                 DEGs_protein_coding = sum(sig & pc, na.rm = TRUE),
                                 paper_protein_coding = paper_pc[i]))
}

cat("\n"); print(summ, row.names = FALSE)
write.table(summ, file.path(outdir, "DEG_summary.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)


