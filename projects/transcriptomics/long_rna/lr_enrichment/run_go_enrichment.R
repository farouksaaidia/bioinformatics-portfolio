#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(optparse)
  library(dplyr)
  library(readr)
})

option_list <- list(
  make_option(c("-i", "--input"), type="character", help="Input DE results file (CSV with gene_id and padj)"),
  make_option(c("-o", "--outdir"), type="character", help="Output directory for enrichment results"),
  make_option(c("-p", "--padj"), type="double", default=0.05, help="Adjusted p-value threshold (default: 0.05)")
)
opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$input) || is.null(opt$outdir)) stop("Provide --input and --outdir")

dir.create(opt$outdir, showWarnings=FALSE, recursive=TRUE)
de <- read.csv(opt$input)
sig_genes <- subset(de, padj < opt$padj)$gene_id
if (length(sig_genes) < 5) stop("Not enough significant genes for enrichment")

ego <- enrichGO(gene = sig_genes, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP", pAdjustMethod = "BH", readable = TRUE)
write.csv(as.data.frame(ego), file.path(opt$outdir, "go_enrichment_results.csv"))
cat("✅ GO enrichment complete. Results saved to", opt$outdir, "\n")
