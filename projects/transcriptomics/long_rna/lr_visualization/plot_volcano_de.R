#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(ggplot2)
  library(optparse)
  library(dplyr)
  library(readr)
})

option_list <- list(
  make_option(c("-i", "--input"), type="character", help="Input DE results file (CSV with gene_id, log2FoldChange, padj)"),
  make_option(c("-o", "--outdir"), type="character", help="Output directory"),
  make_option(c("-p", "--padj"), type="double", default=0.05, help="Adjusted p-value threshold (default 0.05)")
)
opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$input) || is.null(opt$outdir)) stop("Provide --input and --outdir")

dir.create(opt$outdir, showWarnings=FALSE, recursive=TRUE)
de <- read.csv(opt$input)
de$status <- ifelse(de$padj < opt$padj & abs(de$log2FoldChange) > 1,
                    ifelse(de$log2FoldChange > 0, "Up", "Down"), "NS")

p <- ggplot(de, aes(x=log2FoldChange, y=-log10(padj), color=status)) +
  geom_point(alpha=0.6, size=1.5) +
  scale_color_manual(values=c("Up"="#e41a1c", "Down"="#377eb8", "NS"="grey70")) +
  geom_vline(xintercept=c(-1,1), linetype="dashed", color="black") +
  geom_hline(yintercept=-log10(opt$padj), linetype="dotted", color="black") +
  theme_minimal(base_size=12) +
  labs(title="Volcano Plot", x="log2 Fold Change", y="-log10 Adjusted p-value")
ggsave(file.path(opt$outdir, "volcano_plot.pdf"), p, width=7, height=6)
cat("✅ Volcano plot generated at", opt$outdir, "\n")
