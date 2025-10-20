#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(dplyr)
})

option_list <- list(
  make_option(c("-i","--input"), type="character", help="Input Seurat .rds (with domain labels)"),
  make_option(c("-d","--domain_col"), type="character", default="bayes_domains", help="Domain/cluster column to compare (default: bayes_domains)"),
  make_option(c("-o","--output_dir"), type="character", help="Output directory for pairwise DE CSVs"),
  make_option(c("--min_pct"), type="double", default=0.05, help="Min percent expressed"),
  make_option(c("--logfc"), type="double", default=0.25, help="Log fold-change threshold")
)
opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$input) || is.null(opt$output_dir)) stop("❌ Provide --input and --output_dir")
dir.create(opt$output_dir, recursive=TRUE, showWarnings=FALSE)

so <- readRDS(opt$input)
if (!(opt$domain_col %in% colnames(so@meta.data))) stop(paste0("❌ Domain column not found: ", opt$domain_col))

domains <- unique(so@meta.data[[opt$domain_col]])
combos <- combn(domains, 2, simplify = FALSE)
message("🔁 Performing pairwise DE for ", length(combos), " domain pairs")

for (pair in combos) {
  a <- pair[1]; b <- pair[2]
  name <- paste0("DE_", gsub('[^A-Za-z0-9]', '_', a), "_vs_", gsub('[^A-Za-z0-9]', '_', b))
  message("🔎 ", name)
  try({
    degs <- FindMarkers(so, ident.1 = a, ident.2 = b, group.by = opt$domain_col,
                        min.pct = opt$min_pct, logfc.threshold = opt$logfc)
    if (nrow(degs) > 0) {
      degs$gene <- rownames(degs)
      out <- file.path(opt$output_dir, paste0(name, ".csv"))
      write.csv(degs[order(degs$p_val_adj), ], out, row.names=FALSE)
      message("✅ Saved ", out)
    } else {
      message("⚠️ No DEGs for pair: ", name)
    }
  }, silent = FALSE)
}

message("🎯 Domain differential expression complete. Outputs in: ", opt$output_dir)
