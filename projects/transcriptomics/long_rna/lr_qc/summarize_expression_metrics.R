#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(ggplot2)
  library(readr)
})

option_list <- list(
  make_option(c("-c", "--counts"), type="character", help="Input normalized counts matrix (CSV/TSV)"),
  make_option(c("-o", "--outdir"), type="character", help="Output directory")
)
opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$counts) || is.null(opt$outdir)) stop("Provide --counts and --outdir")

dir.create(opt$outdir, showWarnings=FALSE, recursive=TRUE)
counts <- read.delim(opt$counts, header=TRUE, row.names=1, sep="\t", check.names=FALSE)

gene_sums <- data.frame(
  Gene=row.names(counts),
  TotalCounts=rowSums(counts),
  DetectedSamples=rowSums(counts > 0)
)
write.csv(gene_sums, file.path(opt$outdir, "gene_expression_summary.csv"), row.names=FALSE)

p <- ggplot(gene_sums, aes(x=log10(TotalCounts+1))) +
  geom_histogram(bins=50, fill="#3182bd", alpha=0.7) +
  labs(title="Gene Count Distribution (log10 scale)", x="log10(Total Counts)", y="Frequency") +
  theme_minimal(base_size=12)
ggsave(file.path(opt$outdir, "gene_count_distribution.pdf"), p, width=7, height=5)

cat("✅ Expression metrics summary complete. Files saved to", opt$outdir, "\n")
