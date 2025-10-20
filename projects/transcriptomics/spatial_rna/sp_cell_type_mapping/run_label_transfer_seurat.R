#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Seurat)
  library(optparse)
  library(dplyr)
})

option_list <- list(
  make_option(c("-r","--reference"), type="character", help="Reference Seurat object with known cell type labels"),
  make_option(c("-q","--query"), type="character", help="Spatial Seurat object to annotate"),
  make_option(c("-o","--output"), type="character", help="Output annotated Seurat .rds")
)
opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$reference) || is.null(opt$query) || is.null(opt$output)) stop("❌ Provide reference, query, and output files")

cat("📥 Loading reference and query Seurat objects...\n")
ref <- readRDS(opt$reference)
qry <- readRDS(opt$query)

cat("🔗 Finding transfer anchors...\n")
anchors <- FindTransferAnchors(reference=ref, query=qry, dims=1:30)
predictions <- TransferData(anchorset=anchors, refdata=ref$cell_type, dims=1:30)

qry <- AddMetaData(qry, metadata=predictions)
qry$predicted_cell_type <- predictions$predicted.id

saveRDS(qry, file=opt$output)
cat("✅ Label transfer complete ->", opt$output, "\n")
