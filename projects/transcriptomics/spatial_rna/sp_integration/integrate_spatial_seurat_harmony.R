#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(harmony)
  library(dplyr)
})

option_list <- list(
  make_option(c("-i","--inputs"), type="character", help="Comma-separated Seurat .rds files (or a merged RDS)"),
  make_option(c("-b","--batch_col"), type="character", default="sample_id", help="Metadata column indicating batch/sample (default: sample_id)"),
  make_option(c("-o","--output"), type="character", help="Output integrated Seurat .rds file"),
  make_option(c("--npcs"), type="integer", default=30, help="Number of PCs for integration (default:30)")
)
opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$inputs) || is.null(opt$output)) stop("❌ Provide --inputs and --output")

files <- unlist(strsplit(opt$inputs, ","))
files <- trimws(files)
if (length(files) == 1) {
  so <- readRDS(files[1])
} else {
  # load and merge with sample_id prefixing
  objs <- lapply(seq_along(files), function(i){
    x <- readRDS(files[i])
    if (is.null(x$sample_id)) x$sample_id <- paste0("sample", i)
    # prefix cell names to avoid duplicates
    cn <- paste0(x$sample_id[1], "_", colnames(x))
    colnames(x) <- cn
    x
  })
  so <- Reduce(function(a,b) merge(a, y=b), objs)
}

if (!(opt$batch_col %in% colnames(so@meta.data))) stop(paste0("❌ batch column '", opt$batch_col, "' not found"))
DefaultAssay(so) <- DefaultAssay(so)

message("🔎 Running SCTransform (if not already present)...")
if (!("SCT" %in% Assays(so))) {
  so <- SCTransform(so, verbose = FALSE, conserve.memory = TRUE)
  DefaultAssay(so) <- "SCT"
}

message("⚙️ Running PCA...")
so <- RunPCA(so, verbose = FALSE, npcs = opt$npcs)

message("🌀 Running Harmony integration on PCs with batch:", opt$batch_col)
so@reductions$pca@cell.embeddings <- HarmonyMatrix(as.matrix(Embeddings(so, "pca")), so@meta.data, vars_use = opt$batch_col)
so <- RunUMAP(so, reduction = "pca", dims = 1:opt$npcs)
so <- FindNeighbors(so, reduction = "pca", dims = 1:opt$npcs)
so <- FindClusters(so, resolution = 0.6)

message("💾 Saving integrated Seurat object to: ", opt$output)
saveRDS(so, file = opt$output)
message("✅ Integration complete.")
