#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(DESeq2)
  library(optparse)
  library(readr)
  library(dplyr)
})

option_list <- list(
  make_option(c("-c","--counts"), type="character", help="Raw counts CSV/TSV (genes x samples)"),
  make_option(c("-m","--metadata"), type="character", help="Sample metadata TSV with condition column"),
  make_option(c("-o","--outdir"), type="character", help="Output directory")
)

opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$counts) || is.null(opt$metadata) || is.null(opt$outdir)) stop("Missing required arguments.")

dir.create(opt$outdir, showWarnings=FALSE, recursive=TRUE)

counts <- read.table(opt$counts, header=TRUE, sep=",", row.names=1, check.names=FALSE)
meta <- read.table(opt$metadata, header=TRUE, sep="\t")

dds <- DESeqDataSetFromMatrix(countData=counts, colData=meta, design=~condition)
dds <- estimateSizeFactors(dds)

norm_counts <- counts(dds, normalized=TRUE)
write.csv(norm_counts, file=file.path(opt$outdir, "normalized_counts_deseq2.csv"))

cat("Normalization complete. Output: normalized_counts_deseq2.csv\n")
