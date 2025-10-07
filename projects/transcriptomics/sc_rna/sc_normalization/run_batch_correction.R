#!/usr/bin/env Rscript

# Batch correction and integration (Seurat)
suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
})

args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 3){
  stop("Usage: Rscript run_batch_correction.R <input_folder> <output_file> <method: harmony/integration>")
}

input_folder <- args[1]
output_file <- args[2]
method <- args[3]

if(!dir.exists(input_folder)) stop("Input folder does not exist")

files <- list.files(input_folder, pattern="*.rds$", full.names = TRUE)
if(length(files) == 0) stop("No RDS files found in input folder")

seurat_list <- lapply(files, readRDS)
if(method == "harmony"){
  combined <- merge(seurat_list[[1]], y=seurat_list[-1])
  combined <- RunHarmony(combined, group.by.vars="batch")
} else if(method == "integration"){
  anchors <- FindIntegrationAnchors(object.list=seurat_list)
  combined <- IntegrateData(anchors)
} else {
  stop("Invalid method")
}

saveRDS(combined, output_file)
cat("Batch correction completed: ", output_file, "\n")
