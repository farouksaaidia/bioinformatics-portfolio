#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(Giotto)
})

option_list <- list(
  make_option(c("-i","--input"), type="character", help="Input Giotto object .rds"),
  make_option(c("-o","--output"), type="character", help="Output CSV file for neighborhood enrichment")
)
opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$input) || is.null(opt$output)) stop("❌ Provide --input and --output")

gobj <- readRDS(opt$input)

cat("📊 Running Giotto spatial neighborhood enrichment analysis...\n")
results <- cellProximityEnrichment(gobject=gobj,
                                   cluster_column="cell_type",
                                   spatial_network_name="Delaunay_network",
                                   adjust_method="fdr")

write.csv(results, opt$output, row.names=FALSE)
cat("✅ Neighborhood enrichment analysis complete ->", opt$output, "\n")
