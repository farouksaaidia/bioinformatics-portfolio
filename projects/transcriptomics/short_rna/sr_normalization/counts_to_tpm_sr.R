#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
})

option_list <- list(
  make_option(c("-c","--counts"), type="character", help="Counts matrix (genes x samples)"),
  make_option(c("-l","--lengths"), type="character", help="File with gene_id and length (bp)"),
  make_option(c("-o","--outdir"), type="character", help="Output directory")
)

opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$counts) || is.null(opt$lengths) || is.null(opt$outdir)) stop("Missing arguments.")

dir.create(opt$outdir, showWarnings=FALSE, recursive=TRUE)

counts <- read.table(opt$counts, header=TRUE, sep=",", row.names=1, check.names=FALSE)
len <- read.table(opt$lengths, header=TRUE, sep="\t")

if (!"length" %in% colnames(len)) stop("Lengths file must contain 'length' column.")

len <- len[match(rownames(counts), len$gene_id), ]
gene_lengths_kb <- len$length / 1000

rpk <- sweep(counts, 1, gene_lengths_kb, "/")
tpm <- sweep(rpk, 2, colSums(rpk) / 1e6, "/")

write.csv(tpm, file=file.path(opt$outdir, "counts_tpm.csv"))
cat("TPM file saved: counts_tpm.csv\n")
