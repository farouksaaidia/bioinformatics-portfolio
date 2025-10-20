#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(BayesSpace)
  library(Seurat)
})

option_list <- list(
  make_option(c("-i","--input"), type="character", help="Input Seurat .rds (with BayesSpace preprocessing done)"),
  make_option(c("-k","--k"), type="integer", default=7, help="Number of spatial domains"),
  make_option(c("-r","--nrep"), type="integer", default=1000, help="Number of MCMC iterations"),
  make_option(c("-o","--output"), type="character", help="Output Seurat .rds (with BayesSpace domains)")
)
opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$input) || is.null(opt$output)) stop("❌ Provide --input and --output")

so <- readRDS(opt$input)
sce <- as.SingleCellExperiment(so)

message("🧠 Running BayesSpace domain detection (k=", opt$k, ") ...")
sce <- spatialCluster(sce, q=opt$k, platform="Visium", d=15, init.method="mclust", nrep=opt$nrep, save.chain=TRUE)

so$bayes_domains <- as.factor(colData(sce)$spatial.cluster)
saveRDS(so, opt$output)
message("✅ Spatial domains inferred -> ", opt$output)
