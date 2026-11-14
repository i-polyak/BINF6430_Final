suppressPackageStartupMessages({
  library(fgsea)
  library(data.table)
  library(readxl)
})

args <- commandArgs(trailingOnly = TRUE)
set_name    <- if (length(args) > 0) args[1] else "matched"   # "matched", "allGEO", or "paperS2"
design_name <- if (length(args) > 1) args[2] else "sv1"       # ignored for paperS2
rank_by     <- if (length(args) > 2) args[3] else "lfc"       # "lfc", "stat", "lfc_bm10", "signp"

proj     <- "/scratch/caldwell.ja/final_project"
paperdir <- file.path(proj, "paper")
gmt_file <- file.path(paperdir, "mh.all.v2023.1.Mm.symbols.gmt")
s1_file  <- file.path(paperdir, "jci.insight.207270.sdt1.xlsx")
s2_file  <- file.path(paperdir, "jci.insight.207270.sdt2.xlsx")

if (set_name == "paperS2") {
  outdir    <- file.path(proj, "results", paste0("GSEA_paperS2_", rank_by))
  run_label <- paste("paperS2", rank_by)
} else {
  indir     <- file.path(proj, "results", paste0("DE_", set_name, "_", design_name))
  outdir    <- file.path(indir, paste0("GSEA_", rank_by))
  run_label <- paste(set_name, design_name, rank_by)
}
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

contrasts <- c(PH_14dOVX_vs_intact  = "PH14dOVXvsPHintact",
               PH_4mOVX_vs_intact   = "PH4mOVXvsPHintact",
               POA_14dOVX_vs_intact = "POA14dOVXvsPOAintact",
               POA_4mOVX_vs_intact  = "POA4mOVXvsPOAintact",
               PH_4mOVX_vs_14dOVX   = "PH4mOVXvsPH14dOVX",
               POA_4mOVX_vs_14dOVX  = "POA4mOVXvsPOA14dOVX")

pathways <- gmtPathways(gmt_file)
cat("Run:", run_label, "| gene sets:", length(pathways), "\n")

# ---- fgsea for each contrast ----
all_res <- list()
for (nm in names(contrasts)) {
  if (set_name == "paperS2") {
    de <- as.data.frame(read_excel(s2_file, sheet = contrasts[[nm]]))
  } else {
    de <- read.delim(file.path(indir, paste0("DEtable_", nm, ".txt")), row.names = 1)
  }
  de <- de[which(de$gene_type == "protein_coding"), ]
  de <- de[!is.na(de$gene_name), ]

  score <- switch(rank_by,
                  lfc      = de$log2FoldChange,
                  stat     = de$stat,
                  lfc_bm10 = ifelse(de$baseMean >= 10, de$log2FoldChange, NA),
                  signp    = sign(de$log2FoldChange) * -log10(de$pvalue))
  keep  <- !is.na(score) & is.finite(score)
  ranks <- sort(setNames(score[keep], de$gene_name[keep]), decreasing = TRUE)

  set.seed(6430)
  res <- fgsea(pathways = pathways, stats = ranks,
               eps = 0.0, minSize = 15, maxSize = 500)
  res$contrast <- nm
  res$leadingEdge <- sapply(res$leadingEdge, paste, collapse = ",")
  fwrite(res, file.path(outdir, paste0("GSEA_", nm, ".txt")), sep = "\t")
  all_res[[nm]] <- res
}
all_res <- rbindlist(all_res)

# ---- match to Table S1 ----
s1 <- as.data.frame(read_excel(s1_file))
names(s1)[names(s1) == "Normalized Enrichment Score"] <- "NES_paper"
names(s1)[names(s1) == "adjusted_pvalue"] <- "padj_paper"

key <- function(x) tolower(gsub("[^A-Za-z0-9]", "", x))

gmt_names   <- sub("^HALLMARK_", "", names(pathways))
paper_names <- unique(s1$Pathway)
d    <- adist(key(paper_names), key(gmt_names))
best <- apply(d, 1, which.min)
pmap <- data.frame(Pathway = paper_names,
                   pathway = paste0("HALLMARK_", gmt_names[best]),
                   edit_distance = d[cbind(seq_along(paper_names), best)])
if (any(duplicated(pmap$pathway))) cat("WARNING: two paper names matched the same gene set\n")

lab <- function(nm) {
  p <- strsplit(nm, "_")[[1]]
  paste(p[1], sub("OVX", "", p[2]), "vs", sub("OVX", "", p[4]))
}
all_res$Comparison <- sapply(all_res$contrast, lab)
missing_lab <- setdiff(key(unique(all_res$Comparison)), key(s1$Comparison))
if (length(missing_lab) > 0) cat("WARNING: no paper match for:", missing_lab, "\n")

s1 <- merge(s1, pmap[, c("Pathway", "pathway")], by = "Pathway")
s1$ckey <- key(s1$Comparison)
all_res$ckey <- key(all_res$Comparison)

m <- merge(as.data.frame(all_res), s1[, c("pathway", "ckey", "NES_paper", "padj_paper")],
           by = c("pathway", "ckey"))
fwrite(m, file.path(outdir, "GSEA_vs_paper_all.txt"), sep = "\t")

summ <- do.call(rbind, lapply(split(m, m$contrast), function(x) data.frame(
  contrast      = x$contrast[1],
  NES_spearman  = round(cor(x$NES, x$NES_paper, method = "spearman", use = "complete.obs"), 3),
  same_sign_pct = round(100 * mean(sign(x$NES) == sign(x$NES_paper), na.rm = TRUE), 1),
  sig_ours      = sum(x$padj < 0.05, na.rm = TRUE),
  sig_paper     = sum(x$padj_paper < 0.05, na.rm = TRUE),
  sig_both      = sum(x$padj < 0.05 & x$padj_paper < 0.05, na.rm = TRUE))))

cat("\nGSEA comparison to Table S1 (", run_label, "):\n")
print(summ, row.names = FALSE)
write.table(summ, file.path(outdir, "GSEA_comparison_summary.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

