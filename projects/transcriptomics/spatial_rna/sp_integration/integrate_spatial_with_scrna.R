#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(dplyr)
})

option_list <- list(
  make_option(c("-r","--scrna"), type="character", help="Single-cell reference Seurat .rds (with cell type labels)"),
  make_option(c("-s","--spatial"), type="character", help="Spatial Seurat .rds to annotate/integrate"),
  make_option(c("-o","--output"), type="character", help="Output annotated spatial Seurat .rds"),
  make_option(c("--transfer_dims"), type="integer", default=30, help="Number of dims for anchor transfer (default 30)")
)
opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$scrna) || is.null(opt$spatial) || is.null(opt$output)) stop("❌ Provide --scrna, --spatial, and --output")

ref <- readRDS(opt$scrna)
qry <- readRDS(opt$spatial)

if (is.null(ref$cell_type)) {
  stop("❌ Reference Seurat must contain cell type labels in metadata column 'cell_type'")
}

# Normalize and find anchors if necessary
if (!("SCT" %in% Assays(ref))) {
  ref <- SCTransform(ref, verbose = FALSE)
}
if (!("SCT" %in% Assays(qry))) {
  qry <- SCTransform(qry, verbose = FALSE)
}

message("🔗 Finding transfer anchors between scRNA and spatial...")
anchors <- FindTransferAnchors(reference = ref, query = qry, normalization.method = "SCT", dims = 1:opt$transfer_dims)

message("📦 Transferring labels and probabilities...")
pred <- TransferData(anchorset = anchors, refdata = ref$cell_type, dims = 1:opt$transfer_dims)
qry <- AddMetaData(qry, metadata = pred)
qry$predicted_cell_type <- pred$predicted.id
qry$prediction_score <- pred$prediction.score.max

message("💾 Saving annotated spatial object ->", opt$output)
saveRDS(qry, opt$output)
message("✅ Integration and label transfer complete.")
