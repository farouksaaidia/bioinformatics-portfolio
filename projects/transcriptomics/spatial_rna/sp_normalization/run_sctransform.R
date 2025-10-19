#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(sctransform)
  library(dplyr)
})

option_list <- list(
  make_option(c("-i","--input"), type="character", help="Input Seurat .rds (spatial)"),
  make_option(c("-o","--output"), type="character", help="Output Seurat .rds (normalized)"),
  make_option(c("--vars_to_regress"), type="character", default=NULL, help="Comma-separated metadata columns to regress (e.g., percent_mt,batch)"),
  make_option(c("--assay_name"), type="character", default="SCT", help="Name for SCTransform assay (default: SCT)"),
  make_option(c("--vst_max_features"), type="integer", default=3000, help="n_genes for variable features (default: 3000)")
)
opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$input) || is.null(opt$output)) stop("❌ Provide --input and --output")

so <- readRDS(opt$input)
message("ℹ️ Running SCTransform on ", opt$input)

vars <- if (!is.null(opt$vars_to_regress)) unlist(strsplit(opt$vars_to_regress, ",")) else NULL
so <- tryCatch({
  if (!is.null(vars)) {
    SCTransform(so, assay = DefaultAssay(so), new.assay.name = opt$assay_name, vars.to.regress = vars, verbose = TRUE, conserve.memory = TRUE)
  } else {
    SCTransform(so, assay = DefaultAssay(so), new.assay.name = opt$assay_name, verbose = TRUE, conserve.memory = TRUE)
  }
}, error = function(e) {
  stop("❌ SCTransform failed: ", e$message)
})

# Set default assay to SCT for downstream usage
DefaultAssay(so) <- opt$assay_name
message("✅ SCTransform complete. Saving to ", opt$output)
saveRDS(so, opt$output)
