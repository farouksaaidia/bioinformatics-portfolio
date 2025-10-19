#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(scran)
  library(scater)
  library(SingleCellExperiment)
})

option_list <- list(
  make_option(c("-i","--input"), type="character", help="Input Seurat .rds"),
  make_option(c("-o","--output"), type="character", help="Output Seurat .rds (with 'scran' normalized assay)"),
  make_option(c("-p","--pool_size"), type="integer", default=20, help="Pool size for computeSumFactors (default:20)")
)
opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$input) || is.null(opt$output)) stop("❌ Provide --input and --output")

so <- readRDS(opt$input)
message("🔁 Converting Seurat to SingleCellExperiment for scran...")
sce <- as.SingleCellExperiment(so)

# perform quick pre-clustering if needed (scran advice)
if (is.null(names(rowData(sce))) && !("cluster" %in% colnames(colData(sce)))) {
  message("ℹ️ No pre-clusters found. Running quick cluster for pool-based size factors.")
  clusters <- quickCluster(sce)
} else {
  clusters <- colData(sce)$cluster %||% NULL
}

message("🧮 Computing size factors with scran (pooling)...")
sce <- computeSumFactors(sce, clusters = clusters, min.mean = 0.1, sizes = seq(10, opt$pool_size, by=10))
sce <- logNormCounts(sce)
# attach normalized assay back to Seurat
norm_mat <- assay(sce, "logcounts")
so[["scran"]] <- CreateAssayObject(counts = expm1(norm_mat))
DefaultAssay(so) <- "scran"
message("✅ Scran pooling normalization complete. Saving to ", opt$output)
saveRDS(so, opt$output)
