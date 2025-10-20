#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(optparse)
})

option_list <- list(
  make_option(c("-i","--input"), type="character", help="Annotated Seurat .rds"),
  make_option(c("-o","--output_pdf"), type="character", help="Output PDF report path")
)
opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$input) || is.null(opt$output_pdf)) stop("❌ Provide input and output")

so <- readRDS(opt$input)
pdf(opt$output_pdf, width=12, height=8)
if ("predicted_cell_type" %in% colnames(so@meta.data)) {
  print(SpatialDimPlot(so, group.by="predicted_cell_type", label=TRUE) + ggtitle("Predicted Cell Types (Seurat Transfer)"))
}
if ("tangram_predicted" %in% colnames(so@meta.data)) {
  print(SpatialDimPlot(so, group.by="tangram_predicted") + ggtitle("Tangram Predicted Mapping"))
}
dev.off()
message("✅ Mapping report generated -> ", opt$output_pdf)
