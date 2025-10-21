#!/usr/bin/env Rscript
# Plot spatial domains/clusters with optional boundaries and cluster composition summary
suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(glue)
})

option_list <- list(
  make_option(c("-i","--input"), type="character", help="Input Seurat .rds (with cluster/domain metadata)"),
  make_option(c("-c","--cluster_col"), type="character", default="seurat_clusters", help="Metadata column for domain/cluster"),
  make_option(c("-o","--outdir"), type="character", help="Output directory for domain plots"),
  make_option(c("--show_counts"), action="store_true", default=TRUE, help="Save cluster counts summary (default TRUE)")
)
opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$input) || is.null(opt$outdir)) stop("❌ Provide --input and --outdir")
dir.create(opt$outdir, recursive=TRUE, showWarnings=FALSE)

so <- readRDS(opt$input)
if (!(opt$cluster_col %in% colnames(so@meta.data))) stop(glue("❌ Cluster column not found: {opt$cluster_col}"))

# Spatial clustered plot
p_sp <- tryCatch({
  SpatialDimPlot(so, group.by = opt$cluster_col, label = TRUE) + ggtitle("Spatial domains")
}, error = function(e) { ggplot() + ggtitle("Failed to create SpatialDimPlot") })

ggsave(file.path(opt$outdir, "spatial_domains.png"), p_sp, width=8, height=6, dpi=300)

# Cluster composition barplot
if (opt$show_counts) {
  comp <- as.data.frame(table(so@meta.data[[opt$cluster_col]]))
  colnames(comp) <- c("cluster","count")
  p_bar <- ggplot(comp, aes(x = cluster, y = count, fill = cluster)) + geom_bar(stat="identity") +
    ggtitle("Cluster counts") + theme_minimal() + theme(axis.text.x = element_text(angle=45, hjust=1))
  ggsave(file.path(opt$outdir, "cluster_counts.png"), p_bar, width=8, height=4, dpi=300)
}

message("✅ Domain plots saved to ", opt$outdir)
