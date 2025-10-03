#!/usr/bin/env Rscript

# Differential expression with NOISeq
suppressMessages(library(NOISeq))

args <- commandArgs(trailingOnly = TRUE)
counts <- read.table(args[1], header=TRUE, row.names=1)
colData <- read.table(args[2], header=TRUE, row.names=1)

# Build NOISeq object
mydata <- readData(data=counts, factors=colData)

# Perform NOISeqBIO (biological replicates)
res <- noiseqbio(mydata, k=0.5, norm="tmm", factor="condition", lc=1, replicates="biological")

write.table(degenes(res, q=0.95, M=NULL), file=args[3], sep=",", quote=FALSE)
