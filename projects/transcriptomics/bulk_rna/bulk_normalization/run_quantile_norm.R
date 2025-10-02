#!/usr/bin/env Rscript

# Quantile Normalization using preprocessCore
suppressPackageStartupMessages({
  library(preprocessCore)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: run_quantile_norm.R <input_counts.txt> <output_norm.txt>")
}

input <- args[1]
output <- args[2]

counts <- read.table(input, header = TRUE, row.names = 1)
norm_counts <- normalize.quantiles(as.matrix(counts))
colnames(norm_counts) <- colnames(counts)
rownames(norm_counts) <- rownames(counts)

write.table(norm_counts, file = output, sep = "\t", quote = FALSE)
