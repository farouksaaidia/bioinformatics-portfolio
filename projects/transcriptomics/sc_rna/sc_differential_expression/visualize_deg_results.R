#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(ggplot2)
  library(Seurat)
  library(pheatmap)
  library(optparse)
})

option_list <- list(
  make_option(c("-i","--input"), type="character", help="Input annotated Seurat .rds file"),
  make_option(c("-d","--deg_dir"), type="character", help="Directory containing DEG CSV files"),
  make_option(c("-o","--output"), type="character", help="Output directory for plots")
)
opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$input) || is.null(opt$deg_dir) || is.null(opt$output)) stop("❌ Missing args")

dir.create(opt$output, showWarnings=FALSE, recursive=TRUE)
so <- readRDS(opt$input)

csvs <- list.files(opt$deg_dir, pattern="\\.csv$", full.names=TRUE)
for (csv in csvs) {
  degs <- read.csv(csv)
  name <- tools::file_path_sans_ext(basename(csv))
  top <- head(degs[order(degs$p_val_adj), ], 20)

  # Volcano plot
  g <- ggplot(degs, aes(x=avg_log2FC, y=-log10(p_val_adj))) +
    geom_point(aes(color=p_val_adj < 0.05), alpha=0.6) +
    ggtitle(paste("Volcano Plot:", name)) +
    theme_minimal()
  ggsave(file.path(opt$output, paste0(name, "_volcano.png")), g)

  # Heatmap of top genes
  if (all(top$gene %in% rownames(so))) {
    m <- GetAssayData(so, slot="data")[top$gene, ]
    pheatmap(m, show_rownames=TRUE, main=paste("Top 20 DEGs:", name),
             filename=file.path(opt$output, paste0(name, "_heatmap.png")))
  }
}
cat("✅ DEG visualization complete. Plots saved to", opt$output, "\n")
