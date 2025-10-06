#!/usr/bin/env Rscript
# Purpose: Detect empty droplets using DropletUtils

suppressMessages(library(DropletUtils))
suppressMessages(library(Matrix))
suppressMessages(library(optparse))

option_list <- list(
  make_option(c("-i","--input"), type="character", help="Input count matrix directory, comma-separated for multiple samples"),
  make_option(c("-o","--output"), type="character", help="Output directory for filtered counts")
)
opt <- parse_args(OptionParser(option_list=option_list))
inputs <- strsplit(opt$input, ",")[[1]]
dir.create(opt$output, showWarnings = FALSE, recursive = TRUE)

for(dir in inputs){
    cat("Processing:", dir, "\n")
    counts <- read10xCounts(dir)
    e.out <- emptyDrops(counts)
    filtered <- counts[, !is.na(e.out$FDR) & e.out$FDR <= 0.01]
    saveRDS(filtered, file=file.path(opt$output, paste0(basename(dir), "_filtered.rds")))
}

cat("Empty droplet detection finished for all samples.\n")
