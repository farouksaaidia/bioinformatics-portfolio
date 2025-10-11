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
  make_option(c("-o", "--output"), type="character", help="Output annotated Seurat .rds file")
)
opt <- parse_args(OptionParser(option_list=option_list))

if (is.null(opt$input) || is.null(opt$markers) || is.null(opt$output)) {
  stop("❌ Missing arguments: --input, --markers, and --output are required.")
}

cat("📂 Loading clustered Seurat object...\n")
seurat_obj <- readRDS(opt$input)
markers <- read.csv(opt$markers, stringsAsFactors = FALSE)

if (!all(c("gene", "cell_type") %in% colnames(markers))) {
  stop("❌ Marker CSV must have columns: gene, cell_type")
}

cat("🧬 Calculating cluster marker enrichment...\n")
cluster_markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, logfc.threshold = 0.25)

cat("🔎 Matching clusters to known marker sets...\n")
annotations <- cluster_markers %>%
  inner_join(markers, by = c("gene" = "gene")) %>%
  group_by(cluster) %>%
  count(cell_type, sort = TRUE) %>%
  slice_max(n, n = 1, with_ties = FALSE) %>%
  select(cluster, cell_type)

if (nrow(annotations) == 0) {
  stop("❌ No overlap between detected markers and reference marker list.")
}

cat("🧩 Assigning cell-type labels to metadata...\n")
seurat_obj$predicted_cell_type <- plyr::mapvalues(
  seurat_obj$seurat_clusters,
  from = annotations$cluster,
  to = annotations$cell_type
)

saveRDS(seurat_obj, opt$output)
cat("✅ Annotation saved to", opt$output, "\n")
