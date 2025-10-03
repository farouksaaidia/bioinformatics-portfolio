#!/usr/bin/env Rscript
# Heatmap of top variable genes
suppressMessages(library(DESeq2))
suppressMessages(library(pheatmap))

args <- commandArgs(trailingOnly = TRUE)
counts <- read.csv(args[1], row.names=1)

topVarGenes <- head(order(rowVars(as.matrix(counts)), decreasing=TRUE), 50)
mat <- counts[topVarGenes, ]
pheatmap(mat, scale="row", filename="heatmap_topgenes.png")
