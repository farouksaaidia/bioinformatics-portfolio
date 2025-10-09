#!/usr/bin/env Rscript
# export_cluster_markers.R
# Identify top marker genes per cluster using Seurat::FindAllMarkers and write results.
#
# Usage:
# Rscript export_cluster_markers.R -i clustered.rds -o markers_folder -c cluster_col -m "wilcox" -p 0.05
suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(dplyr)
})

option_list <- list(
  make_option(c("-i","--input"), type="character", help="Input clustered Seurat RDS"),
  make_option(c("-o","--output"), type="character", help="Output folder for markers"),
  make_option(c("-c","--cluster_col"), type="character", default="seurat_clusters", help="Cluster column in metadata"),
  make_option(c("-m","--method"), type="character", default="wilcox", help="Test method for FindAllMarkers"),
  make_option(c("-p","--pval"), type="numeric", default=0.05, help="Adjusted p-value cutoff")
)
opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$input) || is.null(opt$output)) stop("Provide --input and --output")
if (!file.exists(opt$input)) stop("Input file not found")

dir.create(opt$output, showWarnings = FALSE, recursive = TRUE)
seu <- readRDS(opt$input)

if (!(opt$cluster_col %in% colnames(seu@meta.data))) stop("Cluster column not found in metadata.")

Idents(seu) <- seu@meta.data[[opt$cluster_col]]
markers <- FindAllMarkers(seu, only.pos = TRUE, test.use = opt$method)
markers_filt <- markers %>% filter(p_val_adj <= opt$pval)
write.csv(markers_filt, file=file.path(opt$output, "markers_all_clusters.csv"), row.names = FALSE)

# Also write one CSV per cluster (top N)
topn <- markers_filt %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
cluster_list <- unique(topn$cluster)
for (cl in cluster_list) {
  df <- topn %>% filter(cluster == cl) %>% arrange(desc(avg_log2FC))
  write.csv(df, file = file.path(opt$output, paste0("cluster_", cl, "_top_markers.csv")), row.names = FALSE)
}
cat("Wrote marker tables to", opt$output, "\n")
