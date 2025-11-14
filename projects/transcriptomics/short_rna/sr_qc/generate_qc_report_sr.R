#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
})

option_list <- list(
  make_option(c("-i","--indir"), type="character", help="Directory containing QC metrics + plots"),
  make_option(c("-o","--output"), type="character", help="Output PDF path")
)

opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$indir) || is.null(opt$output)) stop("Required: --indir --output")

pdf(opt$output)

imgs <- list.files(opt$indir, pattern="pdf$", full.names=TRUE)
for (f in imgs) {
  tryCatch({
    grid::grid.newpage()
    grid::grid.raster(png::readPNG(f))
  }, error=function(e) {})
}

dev.off()
cat("QC summary report saved.\n")
