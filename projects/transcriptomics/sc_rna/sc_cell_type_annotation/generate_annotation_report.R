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

dir.create(opt$output, showWarnings=FALSE, recursive=TRUE)
seurat_obj <- readRDS(opt$input)

pdf(file.path(opt$output, "annotation_summary.pdf"))
print(DimPlot(seurat_obj, group.by="cell_type", label=TRUE) + ggtitle("Cell Type Annotation"))
dev.off()

write.csv(table(seurat_obj$cell_type), file=file.path(opt$output, "cell_type_counts.csv"))
cat("✅ Report generated in", opt$output, "\n")
