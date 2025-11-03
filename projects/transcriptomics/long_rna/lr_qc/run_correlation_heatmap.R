#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(pheatmap)
  library(readr)
})

option_list <- list(
  make_option(c("-c", "--counts"), type="character", help="Normalized counts matrix (CSV/TSV)"),
  make_option(c("-o", "--outdir"), type="character", help="Output directory")
)
opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$counts) || is.null(opt$outdir)) stop("Provide --counts and --outdir")

dir.create(opt$outdir, showWarnings=FALSE, recursive=TRUE)
counts <- read.delim(opt$counts, row.names=1)
corr <- cor(counts, method="pearson")

pheatmap(corr, color=colorRampPalette(c("white", "steelblue"))(50),
         main="Sample Correlation Heatmap", fontsize=10, filename=file.path(opt$outdir, "sample_correlation_heatmap.pdf"))

cat("✅ Correlation heatmap generated at", opt$outdir, "\n")
