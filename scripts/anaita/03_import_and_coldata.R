# 03_import_and_coldata.R
# BINF 6430 reproducibility project — build tx2gene, run tximport, save txi.rds
#
# Run as an sbatch job -- see 03_submit_import.sbatch
#
# Comparison to Project 07 (Reagan/Salmon data), same stage:
#   - Project 07 is handed a ready-made tx2gene_gencode_v46.tsv, with
#     transcript ID and version already split into separate columns, and
#     just pastes them back together. Here, tx2gene has to be built from
#     the raw GENCODE vM30 fasta headers, which are fully pipe-delimited
#     with no separate version column -- hence the strsplit/sub steps below.
#   - Project 07 imports Salmon's quant.sf output; this imports kallisto's
#     abundance.h5. Different quantifier, same tximport() call -- the
#     type = "kallisto" argument is the only functional difference in the
#     import call itself.
#   - Project 07's workbook warns that an annotation-release mismatch
#     between tx2gene and the index can silently drop transcripts with no
#     error. The same risk applies here, compounded by kallisto's internal
#     ID-stripping behavior (see the tx2gene note below) -- which is exactly
#     the problem worked through earlier in this project.

suppressPackageStartupMessages({
  library(Biostrings)
  library(tximport)
  library(rhdf5)
})

## --- paths ---
ref_fasta <- "~/binf6430_backup/final_project/reference/gencode.vM30.transcripts.fa.gz"
quant_dir <- "~/binf6430_backup/final_project/quant/GSE288245"

## --- build tx2gene ---
# NOTE: kallisto's .h5 output strips GENCODE's pipe-delimited header down to
# the first field (transcript ID) and drops the version suffix internally
# before matching against tx2gene. tx2gene$TXNAME must therefore be
# pre-stripped to bare Ensembl IDs (no version, no pipe) or tximport reports
# a 100% transcript mismatch even though the IDs look compatible.
#
# Project 07 comparison: its tx2gene is built with
#   transcript_id = paste(<FILL_IN>, <FILL_IN>, sep = ".")
# i.e. ID and version are already separate columns in the provided TSV, so
# the student just pastes them back into one string. There is no such
# pre-split file here -- the version has to be peeled off the single pipe
# string instead, which is what sub("\\..*", ...) below does.
headers <- names(readDNAStringSet(ref_fasta))

tx2gene <- data.frame(
  TXNAME = sub("\\..*", "", sapply(strsplit(headers, "\\|"), `[`, 1)),
  GENEID = sapply(strsplit(headers, "\\|"), `[`, 2)
)

cat("tx2gene preview:\n")
print(head(tx2gene))

## --- build file list ---
samples <- list.dirs(quant_dir, full.names = FALSE, recursive = FALSE)
files <- setNames(file.path(quant_dir, samples, "abundance.h5"), samples)

stopifnot(length(files) == 33)
cat(sprintf("\nFound %d sample files.\n", length(files)))

## --- tximport ---
# Project 07 comparison: same function, type = "salmon" there instead of
# "kallisto" here, reading quant.sf instead of abundance.h5. The console
# message tximport prints about transcripts missing from tx2gene is the
# same diagnostic Project 07's write-up prompt 1 asks you to read and
# report -- worth checking that line in this job's .out log the same way.
txi <- tximport(files, type = "kallisto", tx2gene = tx2gene, ignoreTxVersion = TRUE)

cat("\ntximport result:\n")
cat("  dim(counts):", paste(dim(txi$counts), collapse = " x "), "\n")
cat("  sample names:", paste(head(colnames(txi$counts)), collapse = ", "), "...\n")

## --- save for next stage ---
saveRDS(txi, "~/binf6430_repro_env/txi.rds")
cat("\nSaved txi.rds — ready for 04_sv1_computation.R.\n")
