#!/usr/bin/env Rscript

# Differential expression with baySeq
suppressMessages(library(baySeq))

args <- commandArgs(trailingOnly = TRUE)
counts <- as.matrix(read.table(args[1], header=TRUE, row.names=1))
colData <- read.table(args[2], header=TRUE, row.names=1)

# Setup groupings
groups <- list(NDE=rep(1, ncol(counts)), DE=c(rep(1, sum(colData$condition==1)), rep(2, sum(colData$condition==2))))

CD <- new("countData", data=counts, replicates=colData$condition, groups=groups)
libsizes(CD) <- getLibsizes(CD)
CD <- getPriors.NB(CD, samplesize=1000, estimation="QL")
CD <- getLikelihoods(CD, pET="BIC")

res <- topCounts(CD, group="DE", number=nrow(counts))
write.csv(res, file=args[3])
