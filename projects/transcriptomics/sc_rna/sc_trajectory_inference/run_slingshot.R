#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(slingshot)
  library(RColorBrewer)
})

option_list <- list(
  make_option(c("-i","--input"), type="character", help="Input Seurat .rds with embeddings and cluster labels"),
  make_option(c("-e","--embedding"), type="character", default="umap", help="Embedding to use (umap or pca)"),
  make_option(c("-c","--cluster_col"), type="character", default="seurat_clusters", help="Metadata cluster column"),
  make_option(c("-o","--output"), type="character", help="Output Seurat .rds with slingshot pseudotime (list of curves)")
)
opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$input) || is.null(opt$output)) stop("❌ Provide --input and --output")

so <- readRDS(opt$input)
emb_name <- tolower(opt$embedding)
if (! (emb_name %in% names(so@reductions))) stop(paste0("❌ Embedding '", emb_name, "' not found in Seurat object."))

emb <- Embeddings(so, emb_name)
if (!(opt$cluster_col %in% colnames(so@meta.data))) stop("❌ cluster_col not found in metadata")

clusters <- as.character(so@meta.data[[opt$cluster_col]])
sds <- slingshot(emb, clusterLabels = clusters, start.clus = NULL, stretch = 1)
# extract pseudotime for each lineage and attach first lineage as default
pt_list <- slingPseudotime(sds)
# add first pseudotime column to meta
so$sling_pseudotime <- as.numeric(pt_list[,1])
# store Slingshot object in misc for downstream plotting
so@misc$slingshot <- sds

saveRDS(so, opt$output)
message("✅ Slingshot pseudotime saved to ", opt$output)
