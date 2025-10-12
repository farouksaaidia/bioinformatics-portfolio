#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(optparse)
})

option_list <- list(
  make_option(c("-i","--input"), type="character", help="Input annotated Seurat .rds file"),
  make_option(c("-g","--group"), type="character", default="cell_type", help="Grouping variable (default: cell_type)"),
  make_option(c("-o","--output"), type="character", help="Output directory for DEGs")
)
opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$input) || is.null(opt$output)) stop("❌ Provide --input and --output")

dir.create(opt$output, showWarnings=FALSE, recursive=TRUE)
so <- readRDS(opt$input)

cat("🧬 Finding differentially expressed genes by", opt$group, "...\n")
idents <- unique(so@meta.data[[opt$group]])
for (id in idents) {
  degs <- FindMarkers(so, ident.1=id, group.by=opt$group, logfc.threshold=0.25)
  outfile <- file.path(opt$output, paste0("DEGs_", gsub("[^A-Za-z0-9]", "_", id), ".csv"))
  write.csv(degs, outfile)
  cat("✅ Saved DEGs for", id, "→", outfile, "\n")
}
cat("🎯 DEG analysis complete. Results in:", opt$output, "\n")
