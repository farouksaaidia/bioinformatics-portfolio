#!/usr/bin/env Rscript

# SCTransform normalization for Seurat
suppressPackageStartupMessages({
  library(Seurat)
})

args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 2){
  stop("Usage: Rscript run_sctransform.R <input_rds> <output_rds>")
}

input_file <- args[1]
output_file <- args[2]

if(!file.exists(input_file)) stop("Input file not found: ", input_file)

tryCatch({
  seurat_obj <- readRDS(input_file)
  seurat_obj <- SCTransform(seurat_obj, verbose = FALSE)
  saveRDS(seurat_obj, output_file)
  cat("SCTransform normalization completed: ", output_file, "\n")
}, error=function(e){
  message("Error during SCTransform: ", e$message)
  quit(status=1)
})
