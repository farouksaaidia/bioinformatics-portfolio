#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(dplyr)
  library(optparse)
})

option_list <- list(
  make_option(c("-d","--deg_dir"), type="character", help="Directory containing DEG CSV files"),
  make_option(c("-o","--output"), type="character", help="Output directory for enrichment results")
)
opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$deg_dir) || is.null(opt$output)) stop("❌ Provide --deg_dir and --output")

dir.create(opt$output, showWarnings=FALSE, recursive=TRUE)
files <- list.files(opt$deg_dir, pattern="\\.csv$", full.names=TRUE)
for (f in files) {
  name <- tools::file_path_sans_ext(basename(f))
  degs <- read.csv(f)
  genes <- degs$gene[degs$p_val_adj < 0.05 & abs(degs$avg_log2FC) > 0.25]
  entrez <- bitr(genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")$ENTREZID
  if (length(entrez) > 0) {
    go <- enrichGO(entrez, OrgDb="org.Hs.eg.db", ont="BP", readable=TRUE)
    kegg <- enrichKEGG(entrez, organism="hsa")
    write.csv(as.data.frame(go), file.path(opt$output, paste0(name, "_GO.csv")))
    write.csv(as.data.frame(kegg), file.path(opt$output, paste0(name, "_KEGG.csv")))
    cat("✅ Enrichment done for", name, "\n")
  }
}
cat("🎯 Functional enrichment complete.\n")
