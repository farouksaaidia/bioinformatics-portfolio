#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(optparse)
})

option_list <- list(
  make_option(c("-i", "--input"), type="character", help="Input Seurat RDS with PCA/UMAP/tSNE"),
  make_option(c("-o", "--output_dir"), type="character", help="Output directory for plots")
)
opt <- parse_args(OptionParser(option_list=option_list))

dir.create(opt$output_dir, showWarnings=FALSE, recursive=TRUE)
obj <- readRDS(opt$input)

plots <- list(
  PCA = DimPlot(obj, reduction="pca") + ggtitle("PCA Plot"),
  UMAP = DimPlot(obj, reduction="umap") + ggtitle("UMAP Plot"),
  TSNE = DimPlot(obj, reduction="tsne") + ggtitle("t-SNE Plot")
)

for (name in names(plots)) {
  outfile <- file.path(opt$output_dir, paste0(name, "_plot.png"))
  ggsave(outfile, plots[[name]], width=6, height=5)
  cat("✅ Saved", outfile, "\n")
}
