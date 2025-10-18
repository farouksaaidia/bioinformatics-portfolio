#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(ggplot2)
  library(patchwork)
})

option_list <- list(
  make_option(c("-i","--input"), type="character", help="Input Seurat .rds file"),
  make_option(c("-o","--output_dir"), type="character", help="Output directory for plots")
)
opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$input) || is.null(opt$output_dir)) stop("❌ Provide input and output_dir")

dir.create(opt$output_dir, showWarnings=FALSE, recursive=TRUE)
so <- readRDS(opt$input)

message("🧭 Generating spatial QC plots...")
p1 <- SpatialFeaturePlot(so, features="log_counts") + ggtitle("Log10 Total Counts")
p2 <- SpatialFeaturePlot(so, features="n_features") + ggtitle("Detected Features")
p3 <- SpatialFeaturePlot(so, features="percent_mt") + ggtitle("Mitochondrial %")

pdf(file.path(opt$output_dir, "spatial_qc_plots.pdf"), width=10, height=8)
print(p1 + p2 + p3)
dev.off()

message("✅ QC plots saved to ", opt$output_dir)
