#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(ggplot2)
  library(gridExtra)
})

option_list <- list(
  make_option(c("-d","--deg_dir"), type="character", help="DEG directory"),
  make_option(c("-e","--enrich_dir"), type="character", help="Enrichment results directory"),
  make_option(c("-o","--output"), type="character", help="Output PDF report")
)
opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$deg_dir) || is.null(opt$output)) stop("❌ Provide --deg_dir and --output")

pdf(opt$output, width=10, height=8)
csvs <- list.files(opt$deg_dir, pattern="\\.csv$", full.names=TRUE)
for (f in csvs) {
  name <- tools::file_path_sans_ext(basename(f))
  degs <- read.csv(f)
  p <- ggplot(degs, aes(x=avg_log2FC, y=-log10(p_val_adj))) +
       geom_point(aes(color=p_val_adj<0.05), alpha=0.6) +
       ggtitle(paste("Volcano Plot -", name))
  print(p)
}
dev.off()
cat("✅ DEG summary report generated:", opt$output, "\n")
