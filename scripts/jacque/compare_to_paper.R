suppressPackageStartupMessages(library(readxl))

args <- commandArgs(trailingOnly = TRUE)
set_name    <- if (length(args) > 0) args[1] else "matched"
design_name <- if (length(args) > 1) args[2] else "sv1"

proj  <- "/scratch/caldwell.ja/final_project"
indir <- file.path(proj, "results", paste0("DE_", set_name, "_", design_name))
xl    <- file.path(proj, "paper", "jci.insight.207270.sdt2.xlsx")

sheet_map <- c(PH_14dOVX_vs_intact  = "PH14dOVXvsPHintact",
               PH_4mOVX_vs_intact   = "PH4mOVXvsPHintact",
               POA_14dOVX_vs_intact = "POA14dOVXvsPOAintact",
               POA_4mOVX_vs_intact  = "POA4mOVXvsPOAintact",
               PH_4mOVX_vs_14dOVX   = "PH4mOVXvsPH14dOVX",
               POA_4mOVX_vs_14dOVX  = "POA4mOVXvsPOA14dOVX")

sp <- function(a, b) round(cor(a, b, method = "spearman", use = "complete.obs"), 3)
pe <- function(a, b) round(cor(a, b, use = "complete.obs"), 3)

summ <- data.frame()
for (nm in names(sheet_map)) {
  ours <- read.delim(file.path(indir, paste0("DEtable_", nm, ".txt")), row.names = 1)
  ours$ensembl <- sub("\\..*", "", rownames(ours))      # strip version
  theirs <- as.data.frame(read_excel(xl, sheet = sheet_map[[nm]]))

  # DEG sets (protein-coding)
  our_deg  <- ours$ensembl[which(ours$padj < 0.05 & ours$gene_type == "protein_coding")]
  padj_deg <- theirs$ensembl[which(theirs$padj < 0.05 & theirs$gene_type == "protein_coding")]
  yes_deg  <- theirs$ensembl[which(theirs$significant == "yes" & theirs$gene_type == "protein_coding")]

  # gene-level agreement on shared genes
  m <- merge(ours, theirs, by = "ensembl", suffixes = c(".ours", ".paper"))
  expr <- m$baseMean.paper >= 10

  summ <- rbind(summ, data.frame(
    contrast          = nm,
    ours              = length(our_deg),
    paper_padj        = length(padj_deg),
    shared_padj       = length(intersect(our_deg, padj_deg)),
    paper_yes         = length(yes_deg),
    shared_yes        = length(intersect(our_deg, yes_deg)),
    genes_shared      = nrow(m),
    baseMean_spearman = sp(m$baseMean.ours, m$baseMean.paper),
    LFC_spearman      = sp(m$log2FoldChange.ours, m$log2FoldChange.paper),
    LFC_pearson_bm10  = pe(m$log2FoldChange.ours[expr], m$log2FoldChange.paper[expr]),
    stat_spearman     = sp(m$stat.ours, m$stat.paper),
    stat_pearson      = pe(m$stat.ours, m$stat.paper)))
}

cat("Sample set:", set_name, "| design:", design_name, "\n")
print(summ, row.names = FALSE)
write.table(summ, file.path(indir, "comparison_to_paper.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)
