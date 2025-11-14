#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(ggplot2)
  library(data.table)
  library(dplyr)
  library(FactoMineR)
  library(factoextra)
})

option_list <- list(
  make_option(c("-c","--counts"), type="character", help="Normalized counts CSV/TSV"),
  make_option(c("-o","--outdir"), type="character", help="Output directory")
)

opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$counts) || is.null(opt$outdir)) stop("Required: --counts --outdir")

dir.create(opt$outdir, showWarnings=FALSE, recursive=TRUE)

counts <- fread(opt$counts)
colnames(counts)[1] <- "gene_id"

mat <- as.matrix(counts[, -1, with=FALSE])
mat <- t(mat)

pca <- PCA(mat, graph=FALSE)
pdf(file.path(opt$outdir, "pca_samples.pdf"))
fviz_pca_ind(pca)
dev.off()

pdf(file.path(opt$outdir, "sample_distributions.pdf"))
boxplot(mat, las=2, outline=FALSE)
dev.off()

cat("QC plots saved in output directory.\n")
