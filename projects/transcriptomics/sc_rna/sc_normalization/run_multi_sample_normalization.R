#!/usr/bin/env Rscript

# Multi-sample normalization wrapper
suppressPackageStartupMessages({
  library(Seurat)
})

args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 3){
  stop("Usage: Rscript run_multi_sample_normalization.R <method: log/sctransform> <input_folder> <output_folder>")
}

method <- args[1]
input_folder <- args[2]
output_folder <- args[3]

if(!dir.exists(input_folder)) stop("Input folder does not exist")
if(!dir.exists(output_folder)) dir.create(output_folder, recursive=TRUE)

files <- list.files(input_folder, pattern="*.rds$", full.names = TRUE)
if(length(files) == 0) stop("No RDS files found in input folder")

for(f in files){
  out_file <- file.path(output_folder, basename(f))
  seurat_obj <- readRDS(f)
  seurat_obj <- switch(method,
                       log = NormalizeData(seurat_obj, normalization.method="LogNormalize", scale.factor=10000),
                       sctransform = SCTransform(seurat_obj, verbose=FALSE),
                       stop("Invalid method"))
  saveRDS(seurat_obj, out_file)
  cat("Normalized ", f, " -> ", out_file, "\n")
}
