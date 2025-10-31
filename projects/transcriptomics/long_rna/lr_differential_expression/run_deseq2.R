#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(DESeq2)
  library(optparse)
  library(readr)
})

option_list <- list(
  make_option(c("-c", "--counts"), type="character", help="Input normalized count matrix (CSV/TSV)"),
  make_option(c("-m", "--metadata"), type="character", help="Sample metadata file (TSV with columns: sample, condition)"),
  make_option(c("-o", "--outdir"), type="character", help="Output directory"),
  make_option(c("-f", "--format"), type="character", default="csv", help="Input format: csv or tsv")
)
opt <- parse_args(OptionParser(option_list=option_list))

if (is.null(opt$counts) || is.null(opt$metadata) || is.null(opt$outdir))
  stop("Provide --counts, --metadata, and --outdir")

dir.create(opt$outdir, showWarnings=FALSE, recursive=TRUE)

read_fun <- ifelse(opt$format == "csv", read.csv, read.delim)
counts <- read_fun(opt$counts, row.names=1)
meta <- read.delim(opt$metadata, header=TRUE, stringsAsFactors=TRUE)

if (!all(colnames(counts) %in% meta$sample))
  stop("Sample names in counts do not match metadata file")

dds <- DESeqDataSetFromMatrix(countData=round(counts), colData=meta, design=~condition)
dds <- DESeq(dds)
res <- results(dds)
write.csv(as.data.frame(res), file.path(opt$outdir, "deseq2_results.csv"), row.names=TRUE)

# Save normalized counts
norm_counts <- counts(dds, normalized=TRUE)
write.csv(norm_counts, file.path(opt$outdir, "normalized_counts.csv"), row.names=TRUE)

cat("✅ DESeq2 analysis complete. Results saved in", opt$outdir, "\n")
