#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(ggplot2)
  library(patchwork)
})

option_list <- list(
  make_option(c("-i","--input"), type="character", help="Input Seurat .rds with clustering results"),
  make_option(c("-o","--output_pdf"), type="character", help="Output PDF report path")
)
opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$input) || is.null(opt$output_pdf)) stop("❌ Provide input and output paths")

so <- readRDS(opt$input)
pdf(opt$output_pdf, width=12, height=8)

if ("seurat_clusters" %in% colnames(so@meta.data)) {
  print(SpatialDimPlot(so, group.by="seurat_clusters") + ggtitle("Seurat Clusters"))
}
if ("bayes_domains" %in% colnames(so@meta.data)) {
  print(SpatialDimPlot(so, group.by="bayes_domains") + ggtitle("BayesSpace Domains"))
}

dev.off()
message("✅ Clustering report generated at ", opt$output_pdf)
