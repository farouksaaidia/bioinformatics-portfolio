#!/usr/bin/env Rscript

# Scaling and regression (Seurat)
suppressPackageStartupMessages({
  library(Seurat)
})

args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 2){
  stop("Usage: Rscript run_scaling.R <input_rds> <output_rds>")
}

input_file <- args[1]
output_file <- args[2]

if(!file.exists(input_file)) stop("Input file not found: ", input_file)

tryCatch({
  seurat_obj <- readRDS(input_file)
  seurat_obj <- ScaleData(seurat_obj, vars.to.regress = c("percent.mt", "nCount_RNA"))
  saveRDS(seurat_obj, output_file)
  cat("Scaling completed: ", output_file, "\n")
}, error=function(e){
  message("Error during scaling: ", e$message)
  quit(status=1)
})
