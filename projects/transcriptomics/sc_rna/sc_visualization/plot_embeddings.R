#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(ggplot2)
  library(patchwork)
})

option_list <- list(
  make_option(c("-i", "--input"), type="character", help="Input Seurat .rds"),
  make_option(c("-m", "--metadata"), type="character", default="cell_type", help="Metadata column to color by"),
  make_option(c("-o", "--output"), type="character", help="Output directory")
)
opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$input) || is.null(opt$output)) stop("❌ Missing input/output")

dir.create(opt$output, showWarnings=FALSE, recursive=TRUE)
so <- readRDS(opt$input)

p1 <- DimPlot(so, reduction="umap", group.by=opt$metadata, label=TRUE) + ggtitle(paste("UMAP -", opt$metadata))
p2 <- DimPlot(so, reduction="tsne", group.by=opt$metadata, label=TRUE) + ggtitle(paste("tSNE -", opt$metadata))

pdf(file.path(opt$output, "embedding_plots.pdf"))
print(p1 + p2)
dev.off()

ggsave(file.path(opt$output, "umap.png"), p1, width=7, height=6)
ggsave(file.path(opt$output, "tsne.png"), p2, width=7, height=6)
cat("✅ Embedding plots generated in", opt$output, "\n")
