#!/usr/bin/env Rscript

# DESeq2 Likelihood Ratio Test
suppressMessages(library(DESeq2))

args <- commandArgs(trailingOnly = TRUE)
counts <- read.table(args[1], header=TRUE, row.names=1)
colData <- read.table(args[2], header=TRUE, row.names=1)

dds <- DESeqDataSetFromMatrix(countData = counts, colData = colData, design = ~ condition)
dds <- DESeq(dds, test="LRT", reduced=~1)

res <- results(dds)
write.csv(as.data.frame(res), paste0(args[3], "_DESeq2_LRT.csv"))
