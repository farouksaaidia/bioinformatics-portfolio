#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(optparse)
})

option_list <- list(
  make_option(c("-i", "--input"), type="character", help="Input annotated Seurat .rds file"),
  make_option(c("-o", "--output"), type="character", help="Output directory for plots and tables")
)
opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$input) || is.null(opt$output)) stop("❌ Provide both --input and --output")

if (!file.exists(opt$input)) stop(paste0("❌ Input file not found: ", opt$input))

dir.create(opt$output, showWarnings=FALSE, recursive=TRUE)

message("📂 Loading Seurat object...")
seurat_obj <- tryCatch(readRDS(opt$input),
                       error = function(e) stop(paste0("❌ Failed to read RDS: ", e)))

# Determine which metadata column holds cell type labels
ct_cols <- c("predicted_cell_type", "cell_type", "predicted_labels")
found_col <- intersect(ct_cols, colnames(seurat_obj@meta.data))
if (length(found_col) == 0) {
  stop("❌ No cell type column found in metadata. Expected one of: predicted_cell_type, cell_type, predicted_labels")
}
ct_col <- found_col[1]
message(paste0("ℹ️ Using metadata column for cell types: ", ct_col))

# Basic counts table
counts_tbl <- table(seurat_obj[[ct_col]])
write.csv(as.data.frame(counts_tbl), file = file.path(opt$output, "cell_type_counts.csv"), row.names = TRUE)

# DimPlot (wrap in tryCatch to avoid crashing)
pdf(file.path(opt$output, "annotation_summary.pdf"))
tryCatch({
  p <- DimPlot(seurat_obj, group.by = ct_col, label = TRUE) + ggtitle("Cell Type Annotation")
  print(p)
}, error = function(e) {
  message(paste0("⚠️ Failed to generate DimPlot: ", e))
})
dev.off()

message(paste0("✅ Report generated in ", opt$output))
