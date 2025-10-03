#!/usr/bin/env Rscript
# MA Plot for DE results
suppressMessages(library(ggplot2))

args <- commandArgs(trailingOnly = TRUE)
res <- read.csv(args[1])   # DE results with log2FoldChange, baseMean, padj

png("ma_plot.png")
with(res, plot(baseMean, log2FoldChange, pch=20, cex=0.5, log="x",
     col=ifelse(padj<0.05, "red", "black")))
abline(h=c(-1,1), col="blue")
dev.off()
