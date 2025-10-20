#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(ggplot2)
  library(pheatmap)
  library(dplyr)
})

option_list <- list(
  make_option(c("-s","--seurat"), type="character", help="Annotated Seurat .rds (with clusters/domains)"),
  make_option(c("-d","--deg_dir"), type="character", help="Directory containing DEG CSV files"),
  make_option(c("-e","--enrich_dir"), type="character", default=NULL, help="Optional enrichment results directory"),
  make_option(c("-o","--output_pdf"), type="character", help="Output PDF report path"),
  make_option(c("--top_n"), type="integer", default=20, help="Top N genes to display per comparison (default 20)")
)
opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$seurat) || is.null(opt$deg_dir) || is.null(opt$output_pdf)) stop("❌ Provide --seurat, --deg_dir, and --output_pdf")

so <- readRDS(opt$seurat)
pdf(opt$output_pdf, width=12, height=10)

csvs <- list.files(opt$deg_dir, pattern="\\.csv$", full.names=TRUE)
for (f in csvs) {
  name <- tools::file_path_sans_ext(basename(f))
  degs <- read.csv(f)
  top <- head(degs[order(degs$p_val_adj), ], opt$top_n)
  genes <- intersect(top$gene, rownames(so))
  if (length(genes) > 0) {
    mat <- as.matrix(GetAssayData(so, slot="data")[genes, ])
    pheatmap(mat, main=paste("Top DEGs:", name))
  } else {
    plot(1,1, type="n", axes=FALSE, xlab="", ylab="")
    text(1,1, paste("No overlapping genes in object for", name))
  }
  # volcano-like quick plot
  if (all(c("avg_log2FC","p_val_adj","gene") %in% colnames(degs))) {
    p <- ggplot(degs, aes(x=avg_log2FC, y=-log10(p_val_adj))) +
      geom_point(aes(size = -log10(p_val_adj)), alpha=0.4) +
      ggtitle(paste("Volcano:", name))
    print(p)
  }
}

dev.off()
message("✅ DEG report generated at ", opt$output_pdf)
