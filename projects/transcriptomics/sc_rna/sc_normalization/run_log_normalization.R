#!/usr/bin/env Rscript

# Log-normalization for scRNA-seq
suppressPackageStartupMessages({
  library(Seurat)
})

args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 2){
  stop("Usage: Rscript run_log_normalization.R <input_rds> <output_rds>")
}

input_file <- args[1]
output_file <- args[2]

if(!file.exists(input_file)){
  stop("Input file does not exist: ", input_file)
}

tryCatch({
  seurat_obj <- readRDS(input_file)
  seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
  saveRDS(seurat_obj, output_file)
  cat("Log-normalization completed: ", output_file, "\n")
}, error=function(e){
  message("Error during log-normalization: ", e$message)
  quit(status=1)
})
