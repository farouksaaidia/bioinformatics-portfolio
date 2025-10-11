#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(optparse)
  library(plyr)
})

option_list <- list(
  make_option(c("-i", "--input"), type="character", help="Input clustered Seurat .rds file"),
  make_option(c("-m", "--markers"), type="character", help="CSV of known markers (columns: gene, cell_type)"),
  make_option(c("-o", "--output"), type="character", help="Output annotated Seurat .rds file"),
  make_option(c("-k", "--k_top"), type="integer", default=20, help="Top N markers per cluster to consider (default: 20)")
)
opt <- parse_args(OptionParser(option_list=option_list))

if (is.null(opt$input) || is.null(opt$markers) || is.null(opt$output)) {
  stop("❌ Missing arguments: --input, --markers, and --output are required.")
}

message("📂 Loading clustered Seurat object...")
seurat_obj <- readRDS(opt$input)
markers <- read.csv(opt$markers, stringsAsFactors = FALSE)

if (!all(c("gene", "cell_type") %in% colnames(markers))) {
  stop("❌ Marker CSV must have columns: gene, cell_type")
}

message("🧬 Calculating cluster marker enrichment (FindAllMarkers)...")
cluster_markers <- tryCatch({
  FindAllMarkers(seurat_obj, only.pos = TRUE, logfc.threshold = 0.25)
}, error = function(e) stop("❌ FindAllMarkers failed: ", e$message))

if (nrow(cluster_markers) == 0) stop("❌ No markers found by FindAllMarkers. Check your Seurat object and parameters.")

message("🔎 Matching clusters to reference marker list using top ", opt$k_top, " genes per cluster...")
top_markers <- cluster_markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = opt$k_top, with_ties = FALSE) %>%
  ungroup()

# Join with reference markers and take majority vote per cluster
annotation_table <- top_markers %>%
  inner_join(markers, by = c("gene" = "gene")) %>%
  group_by(cluster, cell_type) %>%
  tally(name = "hits") %>%
  arrange(cluster, desc(hits)) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  select(cluster, cell_type)

if (nrow(annotation_table) == 0) {
  stop("❌ No overlap between detected markers and reference marker list.")
}

message("🧩 Assigning cell-type labels to metadata (key: predicted_cell_type)...")
seurat_obj$predicted_cell_type <- plyr::mapvalues(
  as.character(seurat_obj$seurat_clusters),
  from = as.character(annotation_table$cluster),
  to = as.character(annotation_table$cell_type),
  warn_missing = FALSE
)

# Any clusters not matched will keep their original cluster id as NA -> replace unmapped with "unknown"
seurat_obj$predicted_cell_type[seurat_obj$predicted_cell_type %in% c("") | is.na(seurat_obj$predicted_cell_type)] <- "unknown"

saveRDS(seurat_obj, opt$output)
message("✅ Annotation saved to ", opt$output)
