#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(ggplot2)
  library(pheatmap)
})

option_list <- list(
  make_option(c("-i","--input"), type="character", help="Input Seurat .rds with pseudotime and DEG outputs in folder"),
  make_option(c("-d","--deg_dir"), type="character", default=NULL, help="Optional directory with pseudotime DE CSVs"),
  make_option(c("-o","--output"), type="character", help="Output PDF report path")
)
opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$input) || is.null(opt$output)) stop("❌ Provide --input and --output")

so <- readRDS(opt$input)
pdf(opt$output, width=10, height=8)
# UMAP colored by available pseudotime columns
ptime_cols <- colnames(so@meta.data)[grepl("pseudotime", colnames(so@meta.data))]
for (col in ptime_cols) {
  if ("umap" %in% names(so@reductions)) {
    p <- FeaturePlot(so, features = col, reduction="umap") + ggtitle(paste("UMAP:", col))
    print(p)
  }
}

# Top pseudotime genes if DEG files present
if (!is.null(opt$deg_dir) && dir.exists(opt$deg_dir)) {
  csvs <- list.files(opt$deg_dir, pattern="\\.csv$", full.names=TRUE)
  for (f in csvs) {
    degs <- read.csv(f)
    top <- head(degs[order(degs$p_val_adj), ], 20)
    genes <- intersect(top$gene, rownames(so))
    if (length(genes) > 0) {
      mat <- as.matrix(GetAssayData(so, slot="data")[genes, ])
      pheatmap(mat, main=paste("Top pseudotime genes:", tools::file_path_sans_ext(basename(f))))
    }
  }
}
dev.off()
message("✅ Trajectory report written to ", opt$output)
