#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Seurat)
  library(optparse)
})

option_list <- list(
  make_option(c("-i", "--input"), type="character", help="Input RDS (PCA-computed Seurat object)"),
  make_option(c("-o", "--output"), type="character", help="Output RDS with UMAP embedding")
)
opt <- parse_args(OptionParser(option_list=option_list))

if (is.null(opt$input) || is.null(opt$output)) {
  stop("Usage: Rscript run_umap_seurat.R -i <input> -o <output>")
}

obj <- readRDS(opt$input)
obj <- RunUMAP(obj, dims = 1:30)
saveRDS(obj, opt$output)
cat("✅ UMAP completed and saved to", opt$output, "\n")
