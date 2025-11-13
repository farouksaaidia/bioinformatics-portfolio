#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(edgeR)
  library(optparse)
  library(readr)
})

option_list <- list(
  make_option(c("-c","--counts"), type="character", help="Raw counts CSV/TSV"),
  make_option(c("-o","--outdir"), type="character", help="Output directory")
)

opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$counts) || is.null(opt$outdir)) stop("Missing required arguments.")

dir.create(opt$outdir, showWarnings=FALSE, recursive=TRUE)

counts <- read.table(opt$counts, header=TRUE, row.names=1, sep=",", check.names=FALSE)
dge <- DGEList(counts=counts)
dge <- calcNormFactors(dge)

cpm_mat <- cpm(dge, normalized.lib.sizes=TRUE)
write.csv(cpm_mat, file=file.path(opt$outdir, "normalized_counts_cpm.csv"))

cat("CPM normalization done: normalized_counts_cpm.csv\n")
