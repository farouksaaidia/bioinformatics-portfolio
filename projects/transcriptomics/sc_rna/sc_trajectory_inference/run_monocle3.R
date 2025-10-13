#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(monocle3)
  library(dplyr)
})

option_list <- list(
  make_option(c("-i","--input"), type="character", help="Input Seurat .rds (normalized, with PCA/UMAP)"),
  make_option(c("-c","--cluster_col"), type="character", default="seurat_clusters", help="Metadata column for clusters (default: seurat_clusters)"),
  make_option(c("-r","--root_cells"), type="character", default=NULL, help="Optional comma-separated cell IDs to use as root for pseudotime"),
  make_option(c("-o","--output"), type="character", help="Output Seurat .rds with pseudotime metadata")
)
opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$input) || is.null(opt$output)) stop("❌ Provide --input and --output")

so <- readRDS(opt$input)

# convert Seurat -> CellDataSet (Monocle3)
cds <- as.cell_data_set(so)

# ensure reducedDims (UMAP) exists; if not, try to compute
if (!("UMAP" %in% reducedDims(cds) %>% names())) {
  if ("umap" %in% names(so@reductions)) {
    reducedDims(cds)$UMAP <- Embeddings(so, "umap")
  } else {
    stop("❌ No UMAP found in Seurat object. Run RunUMAP first.")
  }
}

# preprocess and learn graph (monocle3 expects preprocessed)
cds <- cluster_cells(cds, reduction_method="UMAP")
cds <- learn_graph(cds, use_partition = TRUE)

# order cells by root
if (!is.null(opt$root_cells)) {
  roots <- unlist(strsplit(opt$root_cells, ","))
  cds <- order_cells(cds, root_cells = roots)
} else {
  # try to find root by cluster with minimum pseudotime after picking a root node
  cds <- order_cells(cds)
}

# attach pseudotime to Seurat metadata
ptime <- pseudotime(cds)
so$monocle3_pseudotime <- as.numeric(ptime[Cells(cds)])

# attach graph coordinates for plotting
if ("UMAP" %in% names(so@reductions)) {
  # nothing additional needed; Monocle graph drawn on UMAP
}

saveRDS(so, opt$output)
message("✅ Monocle3 pseudotime computed and saved to ", opt$output)
