#!/usr/bin/env Rscript
# Purpose: Preprocess scRNA-seq with Seurat (normalization, QC, filtering)

suppressMessages(library(Seurat))
suppressMessages(library(optparse))

option_list <- list(
  make_option(c("-i","--input"), type="character", help="Input count matrix or RDS, comma-separated for multiple samples"),
  make_option(c("-o","--output"), type="character", help="Output directory to save processed Seurat objects")
)
opt <- parse_args(OptionParser(option_list=option_list))

inputs <- strsplit(opt$input, ",")[[1]]
dir.create(opt$output, showWarnings = FALSE, recursive = TRUE)

for(sample_file in inputs){
    cat("Processing:", sample_file, "\n")
    seurat_obj <- if(grepl("\\.rds$", sample_file)) readRDS(sample_file) else Read10X(sample_file) %>% CreateSeuratObject()
    seurat_obj <- NormalizeData(seurat_obj)
    seurat_obj <- FindVariableFeatures(seurat_obj)
    seurat_obj <- ScaleData(seurat_obj)
    saveRDS(seurat_obj, file=file.path(opt$output, paste0(basename(sample_file), "_processed.rds")))
}

cat("Preprocessing finished for all samples.\n")
