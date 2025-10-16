#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(ggplot2)
  library(dplyr)
})

option_list <- list(
  make_option(c("-i", "--input"), type="character", help="Input Seurat .rds"),
  make_option(c("-g", "--group"), type="character", default="cell_type", help="Grouping variable (e.g. cell_type)"),
  make_option(c("-c", "--condition"), type="character", default="condition", help="Condition/sample variable"),
  make_option(c("-o", "--output"), type="character", help="Output directory")
)
opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$input) || is.null(opt$output)) stop("❌ Missing input/output")

dir.create(opt$output, showWarnings=FALSE, recursive=TRUE)
so <- readRDS(opt$input)

df <- so@meta.data %>%
  count(!!sym(opt$group), !!sym(opt$condition)) %>%
  group_by(!!sym(opt$condition)) %>%
  mutate(freq = n / sum(n))

p <- ggplot(df, aes_string(x=opt$condition, y="freq", fill=opt$group)) +
  geom_bar(stat="identity", position="fill") +
  theme_minimal() +
  ylab("Fraction") + xlab(opt$condition) +
  ggtitle(paste("Cluster Composition by", opt$condition))

ggsave(file.path(opt$output, "cluster_composition.png"), p, width=8, height=6)
cat("✅ Cluster composition plot saved in", opt$output, "\n")
