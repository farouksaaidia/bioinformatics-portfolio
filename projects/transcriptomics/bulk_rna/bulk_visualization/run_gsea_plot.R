#!/usr/bin/env Rscript
# GSEA enrichment plot
suppressMessages(library(fgsea))
suppressMessages(library(ggplot2))

args <- commandArgs(trailingOnly = TRUE)
ranks <- read.csv(args[1], row.names=1)   # ranked gene list
pathways <- gmtPathways(args[2])          # gmt file

fgseaRes <- fgsea(pathways, stats=ranks[,1], minSize=15, maxSize=500)
plotEnrichment(pathways[[1]], stats=ranks[,1])
ggsave("gsea_plot.png")
