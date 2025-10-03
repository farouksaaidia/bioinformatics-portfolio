#!/usr/bin/env Rscript
# Volcano Plot
suppressMessages(library(ggplot2))

args <- commandArgs(trailingOnly = TRUE)
res <- read.csv(args[1])

res$threshold <- as.factor(res$padj < 0.05 & abs(res$log2FoldChange) > 1)
ggplot(res, aes(x=log2FoldChange, y=-log10(padj), col=threshold)) +
  geom_point(alpha=0.5) +
  theme_minimal()
ggsave("volcano_plot.png")
