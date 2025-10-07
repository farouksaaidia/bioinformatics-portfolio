#!/usr/bin/env Rscript

# Scran pooling-based normalization
suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(scran)
})

args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 2){
  stop("Usage: Rscript run_scran_pooling.R <input_sce_rds> <output_sce_rds>")
}

input_file <- args[1]
output_file <- args[2]

if(!file.exists(input_file)) stop("Input file not found: ", input_file)

tryCatch({
  sce <- readRDS(input_file)
  clusters <- quickCluster(sce)
  sce <- computeSumFactors(sce, clusters=clusters)
  sce <- logNormCounts(sce)
  saveRDS(sce, output_file)
  cat("Scran pooling normalization completed: ", output_file, "\n")
}, error=function(e){
  message("Error during Scran normalization: ", e$message)
  quit(status=1)
})
