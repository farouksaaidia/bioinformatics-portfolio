#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(dplyr)
})

option_list <- list(
  make_option(c("-i","--input"), type="character", help="Input Seurat .rds (with bayespace_enhanced assay or bayespace cluster metadata)"),
  make_option(c("-a","--assay"), type="character", default="bayespace_enhanced", help="Assay to use for PCA (default: bayespace_enhanced)"),
  make_option(c("-n","--npcs"), type="integer", default=30, help="Number of PCs"),
  make_option(c("-o","--output"), type="character", help="Output Seurat .rds with bayespace PCA reduction")
)
opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$input) || is.null(opt$output)) stop("❌ Provide --input and --output")

so <- readRDS(opt$input)
message("📥 Loaded Seurat object: ", opt$input)

if (!(opt$assay %in% Assays(so))) stop(paste0("❌ Assay not found: ", opt$assay))
DefaultAssay(so) <- opt$assay

message("🔎 Scaling and running PCA on assay: ", opt$assay)
so <- ScaleData(so, verbose = FALSE)
so <- RunPCA(so, npcs = opt$npcs, verbose = FALSE)
# rename reduction to bayesPCA to avoid clashes
red_name <- "bayesPCA"
reductions(so)[[red_name]] <- reductions(so)[["pca"]]
so@reductions[["pca"]] <- NULL

message("💾 Saving Seurat object with BayesSpace PCA to ", opt$output)
saveRDS(so, opt$output)
message("✅ Done.")
