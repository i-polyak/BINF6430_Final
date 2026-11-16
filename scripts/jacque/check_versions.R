pkgs <- c("tximport", "DESeq2", "fgsea", "SmartSVA", "sva", "isva", "msigdbr")

cat (R.version.string, "\n")
cat("Bioconductor", as.character(BiocManager::version()), "\n")
for (p in pkgs) {
	v <- tryCatch(as.character(packageVersion(p)), error=function(e) "NOT IBSTALLED")
	cat(sprintf("%-10s %s\n", p, v))
}

sink("session_info.txt")
sessionInfo()
sink()

