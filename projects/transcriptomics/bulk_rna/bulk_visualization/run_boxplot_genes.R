#!/usr/bin/env Rscript
# Boxplot for selected genes
suppressMessages(library(ggplot2))

args <- commandArgs(trailingOnly = TRUE)
counts <- read.csv(args[1], row.names=1)
condition <- read.csv(args[2])
gene <- args[3]

data <- data.frame(expr=counts[gene, ], condition=condition$condition)
ggplot(data, aes(x=condition, y=expr, fill=condition)) +
  geom_boxplot() + theme_minimal() +
  ggtitle(paste("Expression of", gene))
ggsave(paste0("boxplot_", gene, ".png"))
