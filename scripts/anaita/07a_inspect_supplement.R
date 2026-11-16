# 07a_inspect_supplement.R
# BINF 6430 reproducibility project — inspect the authors' supplemental tables
# before writing the actual comparison, so the join/merge logic is built on
# real column names rather than a guess (same discipline as checking
# abundance.tsv's target_id format before building tx2gene).

suppressPackageStartupMessages({
  library(readxl)
})

files <- list(
  S1_GSEA   = "~/binf6430_repro_env/supplement_S1_GSEA.xlsx",
  S2_DESeq2 = "~/binf6430_repro_env/supplement_S2_DESeq2.xlsx",
  S3_counts = "~/binf6430_repro_env/supplement_S3_counts.xlsx"
)

for (label in names(files)) {
  path <- files[[label]]
  if (!file.exists(path)) {
    cat(sprintf("\n%s: FILE NOT FOUND at %s\n", label, path))
    next
  }
  cat(sprintf("\n=== %s ===\n", label))
  sheets <- excel_sheets(path)
  cat("Sheets:", paste(sheets, collapse = ", "), "\n")
  for (sh in sheets) {
    cat(sprintf("\n-- sheet: %s --\n", sh))
    df <- read_excel(path, sheet = sh, n_max = 5)
    print(dim(df))
    print(colnames(df))
    print(head(df, 3))
  }
}
