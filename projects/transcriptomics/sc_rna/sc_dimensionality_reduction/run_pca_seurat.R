#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Seurat)
  library(optparse)
})

option_list <- list(
  make_option(c("-i", "--input"), type="character", help="Input RDS (normalized Seurat object)"),
  make_option(c("-o", "--output"), type="character", help="Output RDS with PCA results")
)
opt <- parse_args(OptionParser(option_list=option_list))

if (is.null(opt$input) || is.null(opt$output)) {
  stop("Usage: Rscript run_pca_seurat.R -i <input> -o <output>")
}

seurat_obj <- readRDS(opt$input)
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(seurat_obj))
saveRDS(seurat_obj, opt$output)
cat("✅ PCA completed and saved to", opt$output, "\n")
