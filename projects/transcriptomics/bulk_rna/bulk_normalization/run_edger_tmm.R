#!/usr/bin/env Rscript
# edgeR TMM normalization
# Input: counts.txt
# Output: edger_tmm_normalized.txt

suppressPackageStartupMessages(library(edgeR))

args <- commandArgs(trailingOnly=TRUE)
input <- args[1]
output <- args[2]

counts <- read.table(input, header=TRUE, row.names=1)
dge <- DGEList(counts=counts)
dge <- calcNormFactors(dge, method="TMM")
norm_counts <- cpm(dge, normalized.lib.sizes=TRUE)
write.table(norm_counts, file=output, sep="\t", quote=FALSE)
