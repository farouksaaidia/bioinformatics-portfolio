#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(dplyr)
})

option_list <- list(
  make_option(c("-b","--bamdir"), type="character", help="Directory containing BAM files"),
  make_option(c("-c","--counts"), type="character", help="Raw or normalized counts CSV/TSV"),
  make_option(c("-o","--outdir"), type="character", help="Output directory for QC metrics")
)

opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$bamdir) || is.null(opt$counts) || is.null(opt$outdir)) {
  stop("Required: --bamdir --counts --outdir")
}

dir.create(opt$outdir, showWarnings=FALSE, recursive=TRUE)

counts <- fread(opt$counts)
colnames(counts)[1] <- "gene_id"

samples <- setdiff(colnames(counts), "gene_id")

qc <- data.frame(
  sample = samples,
  detected_genes = apply(counts[ , ..samples], 2, function(x) sum(x > 0))
)

write.csv(qc, file.path(opt$outdir, "qc_metrics.csv"), row.names=FALSE)

cat("QC metrics saved to qc_metrics.csv\n")
