#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(ggplot2)
  library(ggrepel)
  library(stats)
})

option_list <- list(
  make_option(c("-r","--raw"), type="character", help="Raw counts TSV (genes x samples) with gene_id"),
  make_option(c("-n","--normalized"), type="character", help="Normalized matrix CSV/TSV (genes x samples) to compare"),
  make_option(c("-m","--metadata"), type="character", default=NULL, help="Optional metadata TSV with sample,condition"),
  make_option(c("-o","--outdir"), type="character", help="Output directory for plots")
)
opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$raw) || is.null(opt$normalized) || is.null(opt$outdir)) stop("Provide --raw, --normalized, --outdir")
if (!file.exists(opt$raw) || !file.exists(opt$normalized)) stop("Input files not found")
dir.create(opt$outdir, recursive=TRUE, showWarnings=FALSE)

raw <- read_delim(opt$raw, delim="\t", col_types = cols())
norm <- read_delim(opt$normalized, delim="\t", col_types = cols())

# align columns
if (!all(colnames(raw)[-1] == colnames(norm)[-1])) {
  warning("Sample column names differ between raw and normalized matrices. Trying to match by names.")
}
samples <- intersect(colnames(raw)[-1], colnames(norm)[-1])
if (length(samples) < 2) stop("Need at least 2 matched samples to visualize")

raw_mat <- as.data.frame(raw[, samples])
rownames(raw_mat) <- raw$gene_id
norm_mat <- as.data.frame(norm[, samples])
rownames(norm_mat) <- norm$gene_id

# Boxplots before/after (log2 scale)
pdf(file.path(opt$outdir, "boxplots_pre_post_log2.pdf"), width=10, height=6)
par(mfrow=c(1,2))
boxplot(log2(as.matrix(raw_mat)+1), main="Raw counts (log2+1)", las=2)
boxplot(log2(as.matrix(norm_mat)+1), main="Normalized (log2+1)", las=2)
dev.off()

# PCA before / after
pr_raw <- prcomp(t(log2(as.matrix(raw_mat)+1)), center=TRUE, scale.=TRUE)
pr_norm <- prcomp(t(log2(as.matrix(norm_mat)+1)), center=TRUE, scale.=TRUE)

pdf(file.path(opt$outdir, "pca_pre_post.pdf"), width=12, height=5)
par(mfrow=c(1,2))
plot(pr_raw$x[,1], pr_raw$x[,2], xlab="PC1", ylab="PC2", main="PCA raw")
text(pr_raw$x[,1], pr_raw$x[,2], labels=rownames(pr_raw$x), cex=0.7, pos=3)
plot(pr_norm$x[,1], pr_norm$x[,2], xlab="PC1", ylab="PC2", main="PCA normalized")
text(pr_norm$x[,1], pr_norm$x[,2], labels=rownames(pr_norm$x), cex=0.7, pos=3)
dev.off()

# density plots of sample distributions
pdf(file.path(opt$outdir, "density_pre_post.pdf"), width=10, height=6)
matplot(density(as.numeric(as.matrix(raw_mat)[,1]))$x, density(as.numeric(as.matrix(raw_mat)[,1]))$y, type='n', xlab="log2 expression", ylab="density", main="Density raw vs normalized")
for (i in 1:ncol(raw_mat)) lines(density(log2(as.numeric(as.matrix(raw_mat)[,i])+1))$x, density(log2(as.numeric(as.matrix(raw_mat)[,i])+1))$y, col=i)
for (i in 1:ncol(norm_mat)) lines(density(log2(as.numeric(as.matrix(norm_mat)[,i])+1))$x, density(log2(as.numeric(as.matrix(norm_mat)[,i])+1))$y, col=i, lty=2)
legend("topright", legend=c("raw","normalized"), lty=c(1,2))
dev.off()

cat("✅ Normalization visualizations written to", opt$outdir, "\n")
