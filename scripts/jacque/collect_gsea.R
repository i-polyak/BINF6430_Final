args   <- commandArgs(trailingOnly = TRUE)
resdir <- if (length(args) > 0) args[1] else "/scratch/caldwell.ja/final_project/results"
rank   <- if (length(args) > 1) args[2] else "stat"

runs <- c("matched_sv1", "matched_none", "matched_batch",
          "allGEO_sv1", "allGEO_none", "allGEO_batch", "authorsSV1_sv1")

all <- data.frame()
for (r in runs) {
  f <- file.path(resdir, paste0("DE_", r), paste0("GSEA_", rank), "GSEA_comparison_summary.txt")
  if (!file.exists(f)) { cat("missing:", f, "\n"); next }
  x <- read.delim(f)
  x$run <- r
  x$modified <- format(file.mtime(f), "%Y-%m-%d %H:%M")
  all <- rbind(all, x[, c("run", "contrast", "NES_spearman", "sig_ours",
                          "sig_paper", "sig_both", "modified")])
}

cat("Results folder:", resdir, "| ranking:", rank, "\n\nFile timestamps:\n")
print(unique(all[, c("run", "modified")]), row.names = FALSE)

first <- function(v) v[1]
cols  <- intersect(runs, unique(all$run))

cat("\nNES Spearman vs Table S1:\n")
nes <- tapply(all$NES_spearman, list(all$contrast, all$run), first)
print(nes[, cols])

cat("\nPaper's significant pathways recovered, sig_both (of sig_ours):\n")
rec <- tapply(paste0(all$sig_both, " (of ", all$sig_ours, ")"),
              list(all$contrast, all$run), first)
print(noquote(rec[, cols]))

cat("\nTotals, sig_both / sig_paper:\n")
print(tapply(all$sig_both, all$run, sum)[cols])
print(tapply(all$sig_paper, all$run, sum)[cols])

write.table(all, file.path(resdir, paste0("GSEA_all_runs_", rank, ".txt")),
            sep = "\t", quote = FALSE, row.names = FALSE)

