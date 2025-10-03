#!/usr/bin/env Rscript
# Barplot of top DEGs
suppressMessages(library(ggplot2))

args <- commandArgs(trailingOnly = TRUE)
res <- read.csv(args[1])

topgenes <- head(res[order(res$padj), ], 10)
ggplot(topgenes, aes(x=reorder(gene, -log2FoldChange), y=log2FoldChange)) +
  geom_bar(stat="identity", fill="steelblue") +
  coord_flip() + theme_minimal()
ggsave("barplot_topgenes.png")
