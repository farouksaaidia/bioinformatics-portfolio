#!/usr/bin/env Rscript

# Differential expression with DESeq2 (Wald test)
suppressMessages(library(DESeq2))

# Args: counts_file, sample_info_file, output_prefix
args <- commandArgs(trailingOnly = TRUE)
counts <- read.table(args[1], header=TRUE, row.names=1)
colData <- read.table(args[2], header=TRUE, row.names=1)

dds <- DESeqDataSetFromMatrix(countData = counts, colData = colData, design = ~ condition)
dds <- DESeq(dds)

res <- results(dds)
write.csv(as.data.frame(res), paste0(args[3], "_DESeq2_results.csv"))
