#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(gridExtra)
})

option_list <- list(
  make_option(c("-i","--input"), type="character", help="Input Seurat .rds with dimensionality reductions"),
  make_option(c("-o","--output_pdf"), type="character", help="Output PDF report path"),
  make_option(c("-f","--features"), type="character", default="nFeature_RNA,percent_mt", help="Comma-separated features to overlay on spatial plots")
)
opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$input) || is.null(opt$output_pdf)) stop("❌ Provide --input and --output_pdf")

so <- readRDS(opt$input)
pdf(opt$output_pdf, width=12, height=8)

# List available reductions
reds <- names(so@reductions)
if (length(reds) == 0) {
  cat("⚠️ No dimensionality reductions found in object.\n")
} else {
  for (r in reds) {
    # attempt to plot first two dims
    try({
      emb <- Embeddings(so, r)
      df <- as.data.frame(emb[,1:2])
      colnames(df) <- c("Dim1","Dim2")
      df$label <- if ("seurat_clusters" %in% colnames(so@meta.data)) so$seurat_clusters else so$sample_id
      p <- ggplot(df, aes(x=Dim1, y=Dim2, color=label)) + geom_point(size=0.6) + ggtitle(paste("Embedding:", r))
      print(p)
    }, silent = TRUE)
  }
}

# Spatial overlays for requested features
features <- unlist(strsplit(opt$features, ","))
for (feat in features) {
  if (feat %in% colnames(so@meta.data)) {
    p <- SpatialFeaturePlot(so, features = feat) + ggtitle(paste("Spatial:", feat))
    print(p)
  }
}

dev.off()
message("✅ Dimensionality reduction report generated at ", opt$output_pdf)
