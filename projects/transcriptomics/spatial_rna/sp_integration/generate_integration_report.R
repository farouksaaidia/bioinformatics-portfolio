#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(ggplot2)
  library(cowplot)
})

option_list <- list(
  make_option(c("-i","--input"), type="character", help="Input integrated Seurat .rds"),
  make_option(c("-b","--batch_col"), type="character", default="sample_id", help="Batch/sample metadata column"),
  make_option(c("-o","--output_pdf"), type="character", help="Output PDF report path")
)
opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$input) || is.null(opt$output_pdf)) stop("❌ Provide --input and --output_pdf")

so <- readRDS(opt$input)
pdf(opt$output_pdf, width=12, height=9)

# UMAP colored by batch and by predicted cell type if available
if ("umap" %in% names(so@reductions)) {
  print(DimPlot(so, reduction="umap", group.by = opt$batch_col) + ggtitle("UMAP by batch"))
}
if ("predicted_cell_type" %in% colnames(so@meta.data)) {
  print(DimPlot(so, reduction="umap", group.by = "predicted_cell_type", label=TRUE) + ggtitle("UMAP by predicted_cell_type"))
}

# Check batch mixing by silhouette-like visual: proportion per cluster per batch (stacked bar)
if ("seurat_clusters" %in% colnames(so@meta.data)) {
  tab <- table(so$seurat_clusters, so@meta.data[[opt$batch_col]])
  prop <- prop.table(tab, margin = 1)
  # quick stacked bar
  df <- as.data.frame(prop)
  colnames(df) <- c("cluster","batch","prop")
  p <- ggplot(df, aes(x=cluster, y=prop, fill=batch)) + geom_bar(stat="identity") + ggtitle("Cluster composition by batch")
  print(p)
}

dev.off()
message("✅ Integration report written to ", opt$output_pdf)
