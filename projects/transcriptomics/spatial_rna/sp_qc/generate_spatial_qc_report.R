#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(ggplot2)
  library(Seurat)
  library(knitr)
  library(rmarkdown)
})

option_list <- list(
  make_option(c("-i","--input"), type="character", help="Input Seurat .rds file"),
  make_option(c("-q","--qc_csv"), type="character", help="QC metrics CSV"),
  make_option(c("-o","--output"), type="character", help="Output PDF report path")
)
opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$input) || is.null(opt$qc_csv) || is.null(opt$output))
  stop("❌ Provide --input, --qc_csv, and --output")

seurat_obj <- readRDS(opt$input)
qc <- read.csv(opt$qc_csv, row.names=1)
summary_stats <- data.frame(
  metric=c("Mean counts","Mean features","Mean %MT"),
  value=c(mean(qc$spots_total_counts, na.rm=T), mean(qc$spots_n_features, na.rm=T), mean(qc$percent_mt, na.rm=T))
)

pdf(opt$output, width=10, height=8)
gridExtra::grid.table(summary_stats)
print(SpatialFeaturePlot(seurat_obj, features="percent_mt") + ggtitle("Mitochondrial Percentage Map"))
dev.off()

cat("✅ Spatial QC report generated →", opt$output, "\n")
