#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(ggplot2)
  library(gridExtra)
})

option_list <- list(
  make_option(c("-i", "--indir"), type="character", help="Directory containing QC plots and tables"),
  make_option(c("-o", "--output"), type="character", help="Output PDF file path")
)
opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$indir) || is.null(opt$output)) stop("Provide --indir and --output")

plots <- list.files(opt$indir, pattern="\\.pdf$", full.names=TRUE)
pdf(opt$output, width=8, height=10)
for (p in plots) {
  grid::grid.newpage()
  grid::grid.text(paste("QC Plot:", basename(p)), y=0.95, gp=grid::gpar(fontsize=14))
  grid::grid.raster(png::readPNG(p), interpolate=TRUE)
}
dev.off()
cat("✅ QC report generated:", opt$output, "\n")
