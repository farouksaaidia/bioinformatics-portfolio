#!/usr/bin/env Rscript
# Compute per-spot QC metrics for spatial transcriptomics data (Seurat Spatial)
suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(FNN)
})

option_list <- list(
  make_option(c("-i","--input"), type="character", help="Input Seurat .rds file"),
  make_option(c("-o","--output"), type="character", help="Output prefix for CSV and .rds"),
  make_option(c("-k","--neighbors"), type="integer", default=6, help="Number of nearest spatial neighbors (default: 6)")
)
opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$input) || is.null(opt$output)) stop("❌ Provide --input and --output")

seurat_obj <- readRDS(opt$input)
message("📊 Calculating spatial QC metrics...")

# Basic QC
seurat_obj$log_counts <- log10(Matrix::rowSums(GetAssayData(seurat_obj, "counts")) + 1)
seurat_obj$n_features <- Matrix::rowSums(GetAssayData(seurat_obj, "counts") > 0)
seurat_obj$percent_mt <- PercentageFeatureSet(seurat_obj, pattern="^MT-")

# Spatial neighbor metrics
coords <- NULL
if (length(seurat_obj@images) > 0) coords <- seurat_obj@images[[1]]@coordinates
if (!is.null(coords)) {
  coords <- as.matrix(coords)
  knn <- FNN::get.knn(coords, k=opt$neighbors)
  seurat_obj$neighbor_mean_counts <- sapply(seq_len(nrow(coords)), function(i) {
    mean(seurat_obj$log_counts[knn$nn.index[i,]])
  })
  seurat_obj$neighbor_mean_features <- sapply(seq_len(nrow(coords)), function(i) {
    mean(seurat_obj$n_features[knn$nn.index[i,]])
  })
}

write.csv(seurat_obj@meta.data, paste0(opt$output, "_qc_metrics.csv"), row.names=TRUE)
saveRDS(seurat_obj, paste0(opt$output, "_qc.rds"))
message("✅ QC metrics computed and saved to ", opt$output)
