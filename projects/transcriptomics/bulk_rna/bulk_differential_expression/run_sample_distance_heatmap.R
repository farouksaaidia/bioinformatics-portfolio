#!/usr/bin/env Rscript
# Sample distance heatmap
suppressMessages(library(DESeq2))
suppressMessages(library(pheatmap))

args <- commandArgs(trailingOnly = TRUE)
counts <- read.csv(args[1], row.names=1)
condition <- read.csv(args[2])

dds <- DESeqDataSetFromMatrix(counts, DataFrame(condition), ~ condition)
vsd <- vst(dds, blind=TRUE)
dists <- dist(t(assay(vsd)))
pheatmap(as.matrix(dists), filename="sample_distance_heatmap.png")
