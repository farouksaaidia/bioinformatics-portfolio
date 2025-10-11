#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(optparse)
  library(scales)
})

option_list <- list(
  make_option(c("-i", "--input"), type="character", help="Input annotated Seurat .rds file"),
  make_option(c("-o", "--output"), type="character", help="Output directory for plots and tables"),
  make_option(c("-g", "--group_by"), type="character", default="predicted_cell_type", help="Metadata column to summarize (default: predicted_cell_type)"),
  make_option(c("--umap_key"), type="character", default="umap", help="Reductions key for UMAP (default: umap)")
)
opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$input) || is.null(opt$output)) stop("❌ Provide both --input and --output")

dir.create(opt$output, showWarnings=FALSE, recursive=TRUE)
seurat_obj <- readRDS(opt$input)

if (!(opt$group_by %in% colnames(seurat_obj@meta.data))) {
  stop(paste0("❌ Metadata column '", opt$group_by, "' not found in Seurat object."))
}

# Safe plotting: check if UMAP exists, otherwise attempt to run RunUMAP if possible
if (! ("umap" %in% names(seurat_obj@reductions)) && ! (opt$umap_key %in% names(seurat_obj@reductions))) {
  message("⚠️ UMAP reduction not found. Attempting to compute a quick UMAP (requires PCA).")
  if (! ("pca" %in% names(seurat_obj@reductions))) {
    message("⚠️ PCA not found. Skipping UMAP plot.")
    generate_umap <- FALSE
  } else {
    seurat_obj <- tryCatch({
      RunUMAP(seurat_obj, reduction = "pca", dims = 1:20, verbose = FALSE)
    }, error = function(e) { message("⚠️ UMAP computation failed: ", e$message); seurat_obj })
    generate_umap <- "umap" %in% names(seurat_obj@reductions)
  }
} else {
  generate_umap <- TRUE
}

pdf(file.path(opt$output, "annotation_summary.pdf"), width = 8, height = 6)
if (generate_umap) {
  print(DimPlot(seurat_obj, reduction = opt$umap_key, group.by = opt$group_by, label = TRUE) + ggtitle("Cell Type Annotation"))
} else {
  plot.new()
  title("UMAP not available - see logs")
}
dev.off()

# Barplot of counts
counts_df <- as.data.frame(table(seurat_obj@meta.data[[opt$group_by]]))
colnames(counts_df) <- c("cell_type", "count")
png(filename = file.path(opt$output, "cell_type_counts_barplot.png"), width = 1000, height = 600)
print(ggplot(counts_df, aes(x = reorder(cell_type, -count), y = count)) +
        geom_bar(stat = "identity") +
        coord_flip() +
        xlab("") + ylab("Cell count") +
        ggtitle("Cell counts per annotated cell type") +
        theme_minimal())
dev.off()

write.csv(counts_df, file = file.path(opt$output, "cell_type_counts.csv"), row.names = FALSE)
message("✅ Report generated in ", opt$output)
