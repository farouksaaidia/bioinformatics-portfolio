#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(ggplot2)
  library(ggrepel)
  library(readr)
})

option_list <- list(
  make_option(c("-c", "--counts"), type="character", help="Normalized counts matrix (CSV/TSV)"),
  make_option(c("-m", "--metadata"), type="character", help="Sample metadata (TSV with columns: sample, condition)"),
  make_option(c("-o", "--outdir"), type="character", help="Output directory")
)
opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$counts) || is.null(opt$metadata) || is.null(opt$outdir)) stop("Provide --counts, --metadata, and --outdir")

dir.create(opt$outdir, showWarnings=FALSE, recursive=TRUE)
counts <- read.delim(opt$counts, row.names=1)
meta <- read.delim(opt$metadata)

pca <- prcomp(t(counts), scale.=TRUE)
pca_df <- as.data.frame(pca$x)
pca_df$sample <- rownames(pca_df)
pca_df <- merge(pca_df, meta, by="sample", all.x=TRUE)

p <- ggplot(pca_df, aes(x=PC1, y=PC2, color=condition, label=sample)) +
  geom_point(size=3) + geom_text_repel(size=3) +
  theme_minimal(base_size=12) + labs(title="PCA of Samples", x="PC1", y="PC2")
ggsave(file.path(opt$outdir, "pca_samples.pdf"), p, width=7, height=6)

cat("✅ PCA QC complete. Results saved to", opt$outdir, "\n")
