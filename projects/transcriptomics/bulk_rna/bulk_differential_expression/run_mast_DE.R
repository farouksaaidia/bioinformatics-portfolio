#!/usr/bin/env Rscript

# Differential expression with MAST
suppressMessages(library(MAST))

args <- commandArgs(trailingOnly = TRUE)
counts <- as.matrix(read.table(args[1], header=TRUE, row.names=1))
colData <- read.table(args[2], header=TRUE, row.names=1)

sca <- FromMatrix(exprsArray = counts, cData = colData)

zlmCond <- zlm(~condition, sca)
summaryCond <- summary(zlmCond, doLRT='condition')
res <- summaryCond$datatable

write.csv(res, file=args[3])
