#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(pheatmap)
  library(dplyr)
  library(readr)
  library(optparse)
})

option_list <- list(
  make_option(c("-c", "--counts"), type="character", help="Normalized counts matrix (CSV/TSV)"),
  make_option(c("-d", "--de"), type="character", help="DE results CSV with gene_id, log2FoldChange, padj"),
  make_option(c("-o", "--outdir"), type="character", help="Output directory"),
  make_option(c("-n", "--topn"), type="integer", default=50, help="Top N genes to display")
)
opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$counts) || is.null(opt$de) || is.null(opt$outdir)) stop("Provide --counts, --de, and --outdir")

dir.create(opt$outdir, showWarnings=FALSE, recursive=TRUE)
counts <- read.delim(opt$counts, row.names=1)
de <- read.csv(opt$de)
top_genes <- head(de[order(de$padj), "gene_id"], opt$topn)
sub_counts <- counts[top_genes, ]
pheatmap(sub_counts, scale="row", clustering_distance_rows="euclidean",
         clustering_distance_cols="correlation", show_rownames=TRUE, show_colnames=TRUE,
         filename=file.path(opt$outdir, "heatmap_top_genes.pdf"))
cat("✅ Heatmap of top DE genes saved to", opt$outdir, "\n")
