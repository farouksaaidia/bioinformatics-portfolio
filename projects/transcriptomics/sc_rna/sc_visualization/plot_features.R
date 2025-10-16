#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(ggplot2)
})

option_list <- list(
  make_option(c("-i","--input"), type="character", help="Input Seurat .rds"),
  make_option(c("-g","--genes"), type="character", help="Comma-separated list of genes"),
  make_option(c("-o","--output"), type="character", help="Output directory")
)
opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$input) || is.null(opt$genes) || is.null(opt$output)) stop("❌ Missing arguments")

dir.create(opt$output, showWarnings=FALSE, recursive=TRUE)
so <- readRDS(opt$input)
genes <- unlist(strsplit(opt$genes, ","))

for (gene in genes) {
  p <- FeaturePlot(so, features=gene) + ggtitle(paste("Expression of", gene))
  ggsave(file.path(opt$output, paste0(gene, "_feature.png")), p, width=6, height=5)
}
cat("✅ Feature plots generated for", length(genes), "genes in", opt$output, "\n")
