#!/usr/bin/env Rscript
# optimize_resolution.R
# Sweep clustering resolutions and produce a summary table + plot of cluster counts per resolution.
#
# Usage:
# Rscript optimize_resolution.R -i input.rds -o results_folder -r "0.2,0.4,0.6,0.8,1.0"
suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(ggplot2)
  library(dplyr)
})

option_list <- list(
  make_option(c("-i","--input"), type="character", help="Input Seurat RDS"),
  make_option(c("-o","--output"), type="character", help="Output folder"),
  make_option(c("-r","--resolutions"), type="character", default="0.2,0.4,0.6,0.8,1.0", help="Comma-separated list")
)
opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$input) || is.null(opt$output)) stop("Provide --input and --output")

if (!file.exists(opt$input)) stop("Input file not found: ", opt$input)
dir.create(opt$output, showWarnings=FALSE, recursive=TRUE)

res <- as.numeric(unlist(strsplit(opt$resolutions,",")))
seu <- readRDS(opt$input)
seu <- FindNeighbors(seu, dims=1:30, verbose=FALSE)
summary_list <- list()
for (r in res) {
  seu_tmp <- FindClusters(seu, resolution=r, verbose=FALSE)
  clusters <- Idents(seu_tmp)
  counts <- table(as.character(clusters))
  summary_list[[as.character(r)]] <- data.frame(resolution=r, n_clusters=length(unique(clusters)), cluster_sizes=I(list(as.numeric(counts))))
}
summary_df <- do.call(rbind, lapply(summary_list, function(x) data.frame(resolution=x$resolution, n_clusters=x$n_clusters)))
write.csv(summary_df, file=file.path(opt$output, "resolution_cluster_counts.csv"), row.names=FALSE)
# Plot
p <- ggplot(summary_df, aes(x=resolution, y=n_clusters)) + geom_line() + geom_point() + theme_minimal() + ggtitle("Clusters per resolution")
ggsave(filename=file.path(opt$output, "resolution_vs_clusters.png"), plot=p, width=6, height=4)
cat("Wrote resolution summary and plot to", opt$output, "\n")
