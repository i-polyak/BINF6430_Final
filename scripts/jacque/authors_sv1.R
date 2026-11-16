proj <- "/scratch/caldwell.ja/final_project"

# authors' SV1 (row names are library names)
auth <- read.delim(file.path(proj, "src/hypo_menopauseNavarro",
                             "Differential Gene Expression Analysis/Mouse/sampleTable_SVs.txt"),
                   row.names = 1)
auth$library <- rownames(auth)
cat("Authors' samples:", nrow(auth), "\n")

# Ildiko's mapping: library name -> SRR
map <- read.delim(file.path(proj, "src/BINF6430_Final/metadata/sample_mapping_GSE288245.tsv"))
map <- map[!is.na(map$SRR) & map$SRR != "", c("SRR", "library_name")]

j <- merge(auth, map, by.x = "library", by.y = "library_name")
cat("Authors' samples matched to an SRR:", nrow(j), "\n")
cat("Authors' samples with no SRR:", paste(setdiff(auth$library, j$library), collapse = ", "), "\n\n")

# correlate with our SV1
for (v in c("matched", "allGEO")) {
  ours <- read.delim(file.path(proj, "results", paste0("sampleTable_SVs_", v, ".txt")))
  x <- merge(j[, c("SRR", "condition", "SV1")], ours[, c("sample", "condition", "SV1")],
             by.x = "SRR", by.y = "sample", suffixes = c(".authors", ".ours"))
  cat(sprintf("%-17s shared: %d | condition labels agree: %d | Pearson r: %.3f | Spearman: %.3f\n",
              v, nrow(x), sum(x$condition.authors == x$condition.ours),
              cor(x$SV1.authors, x$SV1.ours),
              cor(x$SV1.authors, x$SV1.ours, method = "spearman")))
}

# new SV table: matched 30 samples, authors' SV1
base <- read.delim(file.path(proj, "results", "sampleTable_SVs_matched.txt"))
base <- base[, c("sample", "condition", "region", "batch", "animal", "in_authors_analysis")]
base$SV1 <- j$SV1[match(base$sample, j$SRR)]
stopifnot(!any(is.na(base$SV1)))

out <- file.path(proj, "results", "sampleTable_SVs_authorsSV1.txt")
write.table(base, out, sep = "\t", quote = FALSE, row.names = FALSE)
cat("\nWrote:", out, "with", nrow(base), "samples\n")

