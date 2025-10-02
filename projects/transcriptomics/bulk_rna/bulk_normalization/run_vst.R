#!/usr/bin/env Rscript

# Variance Stabilizing Transformation (VST) with DESeq2
suppressPackageStartupMessages({
  library(DESeq2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: run_vst.R <input_counts.txt> <sample_info.txt> <output_norm.txt>")
}

counts_file <- args[1]
sample_file <- args[2]
output <- args[3]

counts <- read.table(counts_file, header = TRUE, row.names = 1)
coldata <- read.table(sample_file, header = TRUE, row.names = 1)

dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~ 1)
vsd <- vst(dds, blind = TRUE)

write.table(assay(vsd), file = output, sep = "\t", quote = FALSE)
