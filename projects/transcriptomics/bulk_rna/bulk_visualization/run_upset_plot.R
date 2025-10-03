#!/usr/bin/env Rscript
# UpSet plot for overlaps
suppressMessages(library(UpSetR))

args <- commandArgs(trailingOnly = TRUE)
input <- read.csv(args[1])   # binary matrix of DEGs per method

png("upset_plot.png")
upset(input, sets=colnames(input), order.by="freq")
dev.off()
