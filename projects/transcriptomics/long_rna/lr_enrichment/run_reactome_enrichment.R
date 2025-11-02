#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(ReactomePA)
  library(optparse)
  library(dplyr)
  library(readr)
})

option_list <- list(
  make_option(c("-i", "--input"), type="character", help="Input DE results (CSV with gene_id, padj)"),
  make_option(c("-o", "--outdir"), type="character", help="Output directory"),
  make_option(c("-p", "--padj"), type="double", default=0.05, help="Adjusted p-value cutoff")
)
opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$input) || is.null(opt$outdir)) stop("Provide --input and --outdir")

dir.create(opt$outdir, showWarnings=FALSE, recursive=TRUE)
de <- read.csv(opt$input)
sig_genes <- subset(de, padj < opt$padj)$gene_id

if (length(sig_genes) < 5) stop("Not enough significant genes for Reactome enrichment")

reactome <- enrichPathway(gene = sig_genes, organism = "human", pAdjustMethod = "BH", readable = TRUE)
write.csv(as.data.frame(reactome), file.path(opt$outdir, "reactome_enrichment_results.csv"))
cat("✅ Reactome enrichment complete. Results saved to", opt$outdir, "\n")
