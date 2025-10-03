#!/usr/bin/env Rscript
# PCA Plot of samples
suppressMessages(library(DESeq2))
suppressMessages(library(ggplot2))

args <- commandArgs(trailingOnly = TRUE)
counts <- read.csv(args[1], row.names=1)   # normalized counts
condition <- read.csv(args[2])             # sample metadata

dds <- DESeqDataSetFromMatrix(counts, DataFrame(condition), ~ condition)
vsd <- vst(dds, blind=TRUE)
pca <- plotPCA(vsd, intgroup="condition")
ggsave("pca_plot.png", plot=pca)
