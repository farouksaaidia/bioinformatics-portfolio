#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(dplyr)
})

option_list <- list(
  make_option(c("-i","--input"), type="character", help="Input clustered Seurat .rds (must contain cluster column)"),
  make_option(c("-c","--cluster_col"), type="character", default="seurat_clusters", help="Cluster/region metadata column (default: seurat_clusters)"),
  make_option(c("-o","--output_dir"), type="character", help="Output directory for marker CSVs"),
  make_option(c("--min_pct"), type="double", default=0.1, help="Min percent expressed (default 0.1)"),
  make_option(c("--logfc"), type="double", default=0.25, help="Log fold-change threshold (default 0.25)")
)
opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$input) || is.null(opt$output_dir)) stop("❌ Provide --input and --output_dir")
dir.create(opt$output_dir, recursive=TRUE, showWarnings=FALSE)

so <- readRDS(opt$input)
if (!(opt$cluster_col %in% colnames(so@meta.data))) {
  stop(paste0("❌ Cluster column not found: ", opt$cluster_col))
}

clusters <- sort(unique(so@meta.data[[opt$cluster_col]]))
message("🧾 Finding markers for ", length(clusters), " clusters using column: ", opt$cluster_col)

for (cl in clusters) {
  message("🔎 Computing markers for cluster: ", cl)
  try({
    markers <- FindMarkers(so, ident.1 = cl, group.by = opt$cluster_col,
                           min.pct = opt$min_pct, logfc.threshold = opt$logfc, only.pos = TRUE)
    if (nrow(markers) > 0) {
      markers$gene <- rownames(markers)
      out <- file.path(opt$output_dir, paste0("markers_cluster_", gsub('[^A-Za-z0-9]', '_', cl), ".csv"))
      write.csv(markers[order(markers$p_val_adj), ], out, row.names=FALSE)
      message("✅ Saved markers to ", out)
    } else {
      message("⚠️ No markers found for cluster: ", cl)
    }
  }, silent = FALSE)
}

message("🎯 Marker detection complete. Outputs in: ", opt$output_dir)
