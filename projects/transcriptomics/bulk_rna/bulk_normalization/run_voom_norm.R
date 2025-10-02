#!/usr/bin/env Rscript
# Limma-voom normalization
# Input: counts.txt
# Output: voom_normalized.txt

suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(edgeR))

args <- commandArgs(trailingOnly=TRUE)
input <- args[1]
output <- args[2]

counts <- read.table(input, header=TRUE, row.names=1)
dge <- DGEList(counts=counts)
dge <- calcNormFactors(dge)
v <- voom(dge, plot=FALSE)
write.table(v$E, file=output, sep="\t", quote=FALSE)
