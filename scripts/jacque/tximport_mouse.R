library(tximport)

proj <- "/scratch/caldwell.ja/final_project"

# samples
samples <- read.delim(file.path(proj, "samplekey.txt"), header = TRUE)
files <- file.path(proj, "kallisto/quant", samples$sample, "abundance.h5")
names(files) <- samples$sample
stopifnot(all(file.exists(files)))
cat("Samples:", length(files), "\n")

# tx2gene from GENCODE vM30 GTF (transcript lines only)
gtf <- read.delim(file.path(proj, "refs/gencode.vM30.annotation.gtf.gz"),
                  header = FALSE, comment.char = "#", quote = "",
                  stringsAsFactors = FALSE)
tx <- gtf[gtf$V3 == "transcript", ]

get_attr <- function(x, key) sub(paste0('.*', key, ' "([^"]+)".*'), "\\1", x)

tx2gene <- data.frame(TXNAME = get_attr(tx$V9, "transcript_id"),
                      GENEID = get_attr(tx$V9, "gene_id"))

gene_annot <- unique(data.frame(gene_id   = get_attr(tx$V9, "gene_id"),
                                gene_name = get_attr(tx$V9, "gene_name"),
                                gene_type = get_attr(tx$V9, "gene_type")))
cat("Transcripts in tx2gene:", nrow(tx2gene), "\n")
cat("Genes:", nrow(gene_annot), "\n")

# tximport (same call as authors)
txi <- tximport(files, type = "kallisto", tx2gene = tx2gene, ignoreAfterBar = TRUE)
cat("Gene x sample count matrix:", dim(txi$counts), "\n")

# save
saveRDS(txi, file.path(proj, "results/txi.rds"))
saveRDS(gene_annot, file.path(proj, "results/gene_annot_M30.rds"))
write.table(txi$abundance, file.path(proj, "results/unadjusted_TPMs.txt"),
            sep = "\t", quote = FALSE)

