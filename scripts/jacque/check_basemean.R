library(readxl)

proj   <- "/scratch/caldwell.ja/final_project"
theirs <- as.data.frame(read_excel(file.path(proj, "paper", "jci.insight.207270.sdt2.xlsx"),
                                   sheet = "PH4mOVXvsPHintact"))

for (s in c("matched", "allGEO")) {
  ours <- read.delim(file.path(proj, "results", paste0("DE_", s, "_none"),
                               "DEtable_PH_4mOVX_vs_intact.txt"), row.names = 1)
  ours$ensembl <- sub("\\..*", "", rownames(ours))
  m <- merge(ours[, c("ensembl", "baseMean")], theirs[, c("ensembl", "baseMean")],
             by = "ensembl", suffixes = c(".ours", ".paper"))
  m <- m[m$baseMean.paper >= 10, ]
  cat(sprintf("%-8s genes: %d   Pearson(log baseMean): %.5f   median ratio ours/paper: %.4f\n",
              s, nrow(m), cor(log(m$baseMean.ours), log(m$baseMean.paper)),
              median(m$baseMean.ours / m$baseMean.paper)))
}

