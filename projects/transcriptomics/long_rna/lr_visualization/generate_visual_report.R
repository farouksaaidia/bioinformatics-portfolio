#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(gridExtra)
  library(ggplot2)
})

option_list <- list(
  make_option(c("-i", "--indir"), type="character", help="Directory with visualization outputs (plots, PDFs)"),
  make_option(c("-o", "--output"), type="character", help="Output PDF report path")
)
opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$indir) || is.null(opt$output)) stop("Provide --indir and --output")

plots <- list.files(opt$indir, pattern="\\.pdf$", full.names=TRUE)
pdf(opt$output, width=8, height=10)
for (p in plots) {
  grid::grid.newpage()
  grid::grid.text(basename(p), y=0.95, gp=grid::gpar(fontsize=14))
  grid::grid.raster(png::readPNG(p), interpolate=TRUE)
}
dev.off()
cat("✅ Visualization report compiled at", opt$output, "\n")
