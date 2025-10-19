#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(BayesSpace)
  library(dplyr)
})

option_list <- list(
  make_option(c("-i","--input"), type="character", help="Input Seurat .rds (Visium)"),
  make_option(c("-k","--k"), type="integer", default=7, help="Number of spatial clusters for BayesSpace (default:7)"),
  make_option(c("-r","--nrep"), type="integer", default=1000, help="Number of MCMC iterations (default:1000)"),
  make_option(c("-o","--output"), type="character", help="Output Seurat .rds (BayesSpace results attached)")
)
opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$input) || is.null(opt$output)) stop("❌ Provide --input and --output")

so <- readRDS(opt$input)
message("ℹ️ Converting Seurat -> SingleCellExperiment for BayesSpace")
sce <- as.SingleCellExperiment(so)

# run BayesSpace enhancement (SPOT-level clustering + smoothing)
message("🏁 Running BayesSpace spatial clustering (k=", opt$k, ")")
set.seed(123)
sce <- spatialPreprocess(sce, platform = "Visium", assay.type = "counts", n.PCs = 15, verbose = FALSE)
sce <- spatialCluster(sce, q = opt$k, platform = "Visium", d = 15, init.method = "mclust", save.chain = TRUE, nrep = opt$nrep)

# attach BayesSpace cluster and enhanced expression if available
bayes_clusters <- colData(sce)$spatial.cluster
so$bayespace_cluster <- factor(bayes_clusters)
if ("enhanced" %in% assayNames(sce)) {
  enhanced <- assay(sce, "enhanced")
  so[["bayespace_enhanced"]] <- CreateAssayObject(counts = as.matrix(enhanced))
}

saveRDS(so, opt$output)
message("✅ BayesSpace processing complete. Saved to ", opt$output)
