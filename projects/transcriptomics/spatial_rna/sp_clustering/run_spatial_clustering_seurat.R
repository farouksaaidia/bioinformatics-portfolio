#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(dplyr)
})

option_list <- list(
  make_option(c("-i","--input"), type="character", help="Input normalized Seurat .rds"),
  make_option(c("-o","--output"), type="character", help="Output clustered Seurat .rds"),
  make_option(c("--resolution"), type="double", default=0.6, help="Resolution for clustering (default: 0.6)"),
  make_option(c("--use_spatial"), action="store_true", default=FALSE, help="Include spatial coordinates in clustering")
)
opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$input) || is.null(opt$output)) stop("❌ Provide both --input and --output")

so <- readRDS(opt$input)
DefaultAssay(so) <- DefaultAssay(so)
message("📥 Running clustering with resolution=", opt$resolution)

# Build graph and cluster
so <- FindNeighbors(so, dims=1:30, verbose=FALSE)
so <- FindClusters(so, resolution=opt$resolution, verbose=FALSE)

# Optionally incorporate spatial coordinates (simple joint embedding)
if (opt$use_spatial && "imagecol" %in% colnames(so@meta.data)) {
  message("🌍 Adding spatial coordinates to clustering")
  coords <- scale(so@meta.data[, c("imagecol","imagerow")])
  combined <- cbind(Embeddings(so, "pca")[,1:20], coords)
  g <- FindNeighbors(so, reduction="pca", graph.name="spatial_nn", dims=1:20)
  so <- FindClusters(so, graph.name="spatial_nn", resolution=opt$resolution, algorithm=1)
}

saveRDS(so, opt$output)
message("✅ Clustering complete -> ", opt$output)
