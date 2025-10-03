#!/usr/bin/env Rscript

# Differential expression with edgeR
suppressMessages(library(edgeR))

args <- commandArgs(trailingOnly = TRUE)
counts <- read.table(args[1], header=TRUE, row.names=1)
colData <- read.table(args[2], header=TRUE, row.names=1)

group <- factor(colData$condition)
y <- DGEList(counts=counts, group=group)
y <- calcNormFactors(y)
design <- model.matrix(~group)
y <- estimateDisp(y, design)

fit <- glmQLFit(y, design)
qlf <- glmQLFTest(fit)
res <- topTags(qlf, n=Inf)

write.csv(res, file=args[3])
