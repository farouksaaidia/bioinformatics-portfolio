#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(optparse)
  library(dplyr)
})

option_list <- list(
  make_option(c("-i", "--input"), type="character", help="Input DE results CSV (must contain gene_id and padj)"),
  make_option(c("-o", "--outdir"), type="character", help="Output directory"),
  make_option(c("-p", "--padj"), type="double", default=0.05, help="Adjusted p-value cutoff")
)
opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$input) || is.null(opt$outdir)) stop("Provide --input and --outdir")

dir.create(opt$outdir, showWarnings=FALSE, recursive=TRUE)
de <- read.csv(opt$input)
sig_genes <- subset(de, padj < opt$padj)$gene_id

if (length(sig_genes) < 5) stop("Not enough significant genes for KEGG analysis")

ekegg <- enrichKEGG(gene = sig_genes, organism = "hsa", pAdjustMethod = "BH")
write.csv(as.data.frame(ekegg), file.path(opt$outdir, "kegg_enrichment_results.csv"))
cat("✅ KEGG enrichment complete. Results saved to", opt$outdir, "\n")
