#!/usr/bin/env Rscript

# Differential expression with limma-voom
suppressMessages(library(limma))
suppressMessages(library(edgeR))

args <- commandArgs(trailingOnly = TRUE)
counts <- read.table(args[1], header=TRUE, row.names=1)
colData <- read.table(args[2], header=TRUE, row.names=1)

group <- factor(colData$condition)
y <- DGEList(counts=counts, group=group)
y <- calcNormFactors(y)

design <- model.matrix(~group)
v <- voom(y, design, plot=FALSE)
fit <- lmFit(v, design)
fit <- eBayes(fit)

res <- topTable(fit, number=Inf, sort.by="P")
write.csv(res, file=args[3])
