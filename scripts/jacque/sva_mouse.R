library(SmartSVA)
library(isva)

args <- commandArgs(trailingOnly = TRUE)
set_name <- if (length(args) > 0) args[1] else "matched"   # "matched" or "allGEO"

proj <- "/scratch/caldwell.ja/final_project"

txi <- readRDS(file.path(proj, "results/txi.rds"))
samples <- read.delim(file.path(proj, "samplekey.txt"), header = TRUE)
stopifnot(identical(colnames(txi$counts), samples$sample))

# choose sample set
if (set_name == "matched") {
  keep <- samples$in_authors_analysis
} else {
  keep <- rep(TRUE, nrow(samples))
}
samples <- samples[keep, ]
counts  <- txi$counts[, keep]
cat("Sample set:", set_name, "|", ncol(counts), "samples\n")

# filter, then log2 of unscaled counts (authors' recipe)
counts <- counts[rowSums(counts) > 1, ]
log2c  <- as.matrix(log2(counts + 1))
cat("Genes used for SVA:", nrow(log2c), "\n")

samples$condition <- factor(samples$condition)
mod1 <- model.matrix(~condition, data = samples)
mod0 <- model.matrix(~1, data = samples)

# number of SVs (authors' recipe: EstDimRMT + 1)
Y.r  <- t(resid(lm(t(log2c) ~ mod1)))
n.sv <- EstDimRMT(Y.r, FALSE)$dim + 1
cat("Generating", n.sv, "SVs\n")

set.seed(6430)
sv.obj <- smartsva.cpp(log2c, mod = mod1, mod0 = mod0, B = 1000, alpha = 1, n.sv = n.sv)

sv <- as.data.frame(sv.obj$sv)
colnames(sv) <- paste0("SV", 1:n.sv)
out <- cbind(samples, sv)

write.table(out, file.path(proj, "results", paste0("sampleTable_SVs_", set_name, ".txt")),
            sep = "\t", quote = FALSE, row.names = FALSE)

# diagnostics: what does SV1 track?
sink(file.path(proj, "results", paste0("SV1_diagnostics_", set_name, ".txt")), split = TRUE)
cat("Seed: 6430 | n.sv:", n.sv, "\n")
r2 <- function(f) round(summary(lm(f, data = out))$r.squared, 3)
cat("\nR-squared of SV1 explained by:\n")
cat("  batch:    ", r2(SV1 ~ batch), "\n")
cat("  region:   ", r2(SV1 ~ region), "\n")
cat("  condition:", r2(SV1 ~ condition), "\n")
cat("\nSV1 mean by condition:\n"); print(round(tapply(out$SV1, out$condition, mean), 3))
sink()

