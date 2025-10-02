#!/usr/bin/env Rscript
# DESeq2 normalization
# Input: counts.txt (gene_id as rownames, samples as columns)
# Output: deseq2_normalized.txt

suppressPackageStartupMessages(library(DESeq2))

args <- commandArgs(trailingOnly=TRUE)
input <- args[1]
output <- args[2]

counts <- read.table(input, header=TRUE, row.names=1)
dds <- DESeqDataSetFromMatrix(countData=counts, colData=DataFrame(condition=factor(rep("A", ncol(counts)))), design=~1)
dds <- estimateSizeFactors(dds)
norm_counts <- counts(dds, normalized=TRUE)
write.table(norm_counts, file=output, sep="\t", quote=FALSE)
