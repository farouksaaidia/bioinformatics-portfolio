#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(ggplot2)
  library(gridExtra)
})

option_list <- list(
  make_option(c("-i", "--input"), type="character", help="Directory with generated plots"),
  make_option(c("-o", "--output"), type="character", help="Output PDF report")
)
opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$input) || is.null(opt$output)) stop("❌ Provide input and output")

files <- list.files(opt$input, pattern="\\.png$", full.names=TRUE)
pdf(opt$output, width=10, height=8)
for (f in files) {
  img <- png::readPNG(f)
  grid::grid.raster(img)
}
dev.off()
cat("✅ Visualization report generated →", opt$output, "\n")
