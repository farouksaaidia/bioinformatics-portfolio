#!/usr/bin/env Rscript
# Pathway enrichment plot
suppressMessages(library(ggplot2))

args <- commandArgs(trailingOnly = TRUE)
enrich <- read.csv(args[1])   # GO/KEGG enrichment results

ggplot(enrich, aes(x=reorder(Term, -log10(pvalue)), y=-log10(pvalue))) +
  geom_bar(stat="identity", fill="darkgreen") +
  coord_flip() + theme_minimal()
ggsave("pathway_enrichment_plot.png")
