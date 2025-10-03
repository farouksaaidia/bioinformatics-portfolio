#!/usr/bin/env Rscript

# Differential expression with Sleuth (Kallisto input)
suppressMessages(library(sleuth))

args <- commandArgs(trailingOnly = TRUE)
sample_info <- read.table(args[1], header=TRUE, stringsAsFactors=FALSE)

so <- sleuth_prep(sample_info, ~ condition)
so <- sleuth_fit(so, ~ condition, 'full')
so <- sleuth_fit(so, ~1, 'reduced')
so <- sleuth_lrt(so, 'reduced', 'full')

res <- sleuth_results(so, 'reduced:full', 'lrt')
write.csv(res, file=args[2])
