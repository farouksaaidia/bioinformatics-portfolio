#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(Matrix)
})

option_list <- list(
  make_option(c("-i","--input"), type="character", help="Input Seurat .rds (normalized)"),
  make_option(c("-o","--output"), type="character", help="Output Seurat .rds with PCA embeddings"),
  make_option(c("-n","--npcs"), type="integer", default=30, help="Number of PCs to compute (default 30)"),
  make_option(c("--spatial_pca"), action="store_true", default=FALSE, help="Also compute a spatially-informed PCA (append coords)"),
  make_option(c("--coords_cols"), type="character", default="imagecol,imagerow", help="Comma-separated image coord columns to use (default: imagecol,imagerow)")
)
opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$input) || is.null(opt$output)) stop("❌ Provide --input and --output")

so <- readRDS(opt$input)
message("📥 Loaded Seurat object: ", opt$input)

# Ensure we have a default assay with data
if (is.null(DefaultAssay(so))) stop("❌ No default assay set in Seurat object")
message("ℹ️ Default assay: ", DefaultAssay(so))

# Find HVGs and run PCA
message("🧭 Finding variable features and running PCA (npcs=", opt$npcs, ")")
so <- FindVariableFeatures(so, selection.method = "vst", nfeatures = 3000, verbose = FALSE)
so <- ScaleData(so, features = VariableFeatures(so), verbose = FALSE)
so <- RunPCA(so, features = VariableFeatures(so), npcs = opt$npcs, verbose = FALSE)

# Optionally compute a spatially-informed PCA by appending scaled coordinates into the expression matrix
if (opt$spatial_pca) {
  coords_cols <- unlist(strsplit(opt$coords_cols, ","))
  coords_cols <- trimws(coords_cols)
  if (!all(coords_cols %in% colnames(so@meta.data))) {
    warning("⚠️ Some coordinate columns not found in metadata: expected ", paste(coords_cols, collapse=", "))
  } else {
    message("🗺️ Computing spatially-informed PCA by appending scaled coords to HVG matrix")
    # build expression matrix for HVGs (genes x cells)
    expr <- as.matrix(GetAssayData(so, slot = "data")[VariableFeatures(so), ])
    # scale coords and append as synthetic 'features'
    coords <- as.data.frame(so@meta.data[, coords_cols, drop=FALSE])
    coords_scaled <- scale(coords)
    # create a combined matrix: rows = genes + coords_cols, cols = cells
    combined <- rbind(expr, t(coords_scaled))
    # center and scale combined features
    combined_scaled <- t(scale(t(combined)))
    # SVD for PCA
    svd_res <- svd(combined_scaled, nu = opt$npcs, nv = 0)
    pcs <- as.data.frame(svd_res$u %*% diag(svd_res$d[1:opt$npcs]))
    colnames(pcs) <- paste0("spatialPCA_", 1:ncol(pcs))
    rownames(pcs) <- colnames(so)
    # attach as dimensional reduction
    so[["spatialPCA"]] <- CreateDimReducObject(embeddings = as.matrix(pcs), key = "sPCA_", assay = DefaultAssay(so))
    message("✅ Spatially-informed PCA stored as reduction 'spatialPCA'")
  }
}

message("💾 Saving Seurat object with PCA to ", opt$output)
saveRDS(so, file = opt$output)
message("✅ Done.")
