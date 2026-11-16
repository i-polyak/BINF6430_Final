suppressPackageStartupMessages(library(ggplot2))

args <- commandArgs(trailingOnly = TRUE)
set_name    <- if (length(args) > 0) args[1] else "matched"
design_name <- if (length(args) > 1) args[2] else "sv1"
rank_by     <- if (length(args) > 2) args[3] else "stat"

proj   <- "/scratch/caldwell.ja/final_project"
infile <- file.path(proj, "results", paste0("DE_", set_name, "_", design_name),
                    paste0("GSEA_", rank_by), "GSEA_vs_paper_all.txt")
outdir <- file.path(proj, "results", "figures")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

m <- read.delim(infile)

long <- rbind(
  data.frame(pathway = m$pathway, Comparison = m$Comparison, source = "Paper",
             NES = m$NES_paper, padj = m$padj_paper),
  data.frame(pathway = m$pathway, Comparison = m$Comparison, source = "Ours",
             NES = m$NES, padj = m$padj))

long$stars <- ifelse(long$padj < 0.001, "***",
              ifelse(long$padj < 0.01,  "**",
              ifelse(long$padj < 0.05,  "*", "")))

# order pathways by the paper's NES in PH 4m vs intact
ref <- m[m$Comparison == "PH 4m vs intact", ]
path_order <- ref$pathway[order(ref$NES_paper)]
long$pathway <- factor(long$pathway, levels = path_order,
                       labels = gsub("_", " ", sub("^HALLMARK_", "", path_order)))

long$Comparison <- factor(long$Comparison,
                          levels = c("PH 14d vs intact", "PH 4m vs intact", "PH 4m vs 14d",
                                     "POA 14d vs intact", "POA 4m vs intact", "POA 4m vs 14d"))
long$source <- factor(long$source, levels = c("Paper", "Ours"))

p <- ggplot(long, aes(x = source, y = pathway, fill = NES)) +
  geom_tile(color = "white") +
  geom_text(aes(label = stars), size = 2.2, vjust = 0.75) +
  facet_grid(. ~ Comparison) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                       limits = c(-3, 3), oob = scales::squish) +
  theme_minimal(base_size = 8) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid  = element_blank(),
        strip.text  = element_text(face = "bold")) +
  labs(x = NULL, y = NULL,
       title = "Hallmark GSEA: published (Table S1) vs. reproduction",
       subtitle = paste0("Ours: ", set_name, " samples, design ", design_name,
                         ", ranked by ", rank_by,
                         "  |  * padj < 0.05, ** < 0.01, *** < 0.001"))

base <- file.path(outdir, paste0("NES_heatmap_", set_name, "_", design_name, "_", rank_by))
ggsave(paste0(base, ".pdf"), p, width = 11, height = 9)
ggsave(paste0(base, ".png"), p, width = 11, height = 9, dpi = 200)
cat("Saved:", paste0(base, ".pdf"), "and .png\n")

