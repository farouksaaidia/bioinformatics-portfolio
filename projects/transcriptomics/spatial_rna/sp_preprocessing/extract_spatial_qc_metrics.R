#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(igraph)
  library(FNN)
})

option_list <- list(
  make_option(c("-i","--input"), type="character", help="Input Seurat .rds (spatial)"),
  make_option(c("-o","--output_prefix"), type="character", help="Output prefix (CSV + updated .rds will be written)"),
  make_option(c("-k","--k_neighbors"), type="integer", default=6, help="Number of nearest spatial neighbors (default:6)")
)
opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$input) || is.null(opt$output_prefix)) stop("❌ Provide --input and --output_prefix")

so <- readRDS(opt$input)

# Basic per-spot metrics
message("🧾 Computing per-spot QC metrics...")
so@meta.data$spots_total_counts <- Matrix::rowSums(GetAssayData(so, slot="counts"))
so@meta.data$spots_n_features <- Matrix::rowSums(GetAssayData(so, slot="counts") > 0)
if ("percent.mt" %in% colnames(so@meta.data)) {
  so@meta.data$percent_mt <- so@meta.data$percent.mt
} else {
  so@meta.data$percent_mt <- PercentageFeatureSet(so, pattern="^MT-")
}

# Spatial neighbor metrics: require coordinates
coords <- NULL
# try common places for coordinates
if (length(so@images) > 0) {
  img <- so@images[[1]]
  if (!is.null(img@coordinates)) coords <- as.data.frame(img@coordinates)
}
if (is.null(coords) && all(c("imagecol","imagerow") %in% colnames(so@meta.data))) {
  coords <- so@meta.data[, c("imagecol","imagerow")]
}
if (is.null(coords)) {
  message("⚠️ No spatial coordinates found; skipping neighbor metrics")
} else {
  message("🔗 Computing spatial neighbor metrics (k=", opt$k_neighbors, ")")
  coords <- as.matrix(coords)
  knn <- FNN::get.knn(coords, k = opt$k_neighbors)
  # mean neighbor counts and mean neighbor features
  mean_nb_counts <- sapply(seq_len(nrow(coords)), function(i) mean(so@meta.data$spots_total_counts[knn$nn.index[i,]]))
  mean_nb_features <- sapply(seq_len(nrow(coords)), function(i) mean(so@meta.data$spots_n_features[knn$nn.index[i,]]))
  so@meta.data$mean_neighbor_counts <- mean_nb_counts
  so@meta.data$mean_neighbor_features <- mean_nb_features
}

# write outputs
csv_out <- paste0(opt$output_prefix, "_spot_qc_metrics.csv")
rds_out <- paste0(opt$output_prefix, ".rds")
write.csv(so@meta.data, csv_out, row.names = TRUE)
saveRDS(so, rds_out)
message("✅ QC metrics written to ", csv_out, " and updated Seurat saved to ", rds_out)
