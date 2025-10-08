#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Seurat)
  library(optparse)
})

option_list <- list(
  make_option(c("-i", "--input"), type="character", help="Input RDS (PCA Seurat object)"),
  make_option(c("-o", "--output"), type="character", help="Output RDS with t-SNE embedding")
)
opt <- parse_args(OptionParser(option_list=option_list))

obj <- readRDS(opt$input)
obj <- RunTSNE(obj, dims = 1:30)
saveRDS(obj, opt$output)
cat("✅ t-SNE completed and saved to", opt$output, "\n")
