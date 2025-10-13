#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(tradeSeq)
  library(splines)
  library(dplyr)
})

option_list <- list(
  make_option(c("-i","--input"), type="character", help="Input Seurat .rds with pseudotime (columns: monocle3_pseudotime, sling_pseudotime, or dpt_pseudotime)"),
  make_option(c("-p","--pseudotime_col"), type="character", default="monocle3_pseudotime", help="Pseudotime column to use"),
  make_option(c("-c","--cluster_col"), type="character", default="seurat_clusters", help="Cluster/lineage column (optional; used for branching tests)"),
  make_option(c("-o","--output"), type="character", help="Output directory for pseudotime DE results")
)
opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$input) || is.null(opt$output)) stop("❌ Provide --input and --output")

dir.create(opt$output, recursive=TRUE, showWarnings=FALSE)
so <- readRDS(opt$input)
if (!(opt$pseudotime_col %in% colnames(so@meta.data))) stop("❌ pseudotime column not found in metadata")

ptime <- so@meta.data[[opt$pseudotime_col]]
# Prepare counts matrix for tradeSeq (requires SingleCellExperiment)
sce <- as.SingleCellExperiment(so)
# remove NA pseudotime cells
valid <- !is.na(ptime)
sce <- sce[, valid]
ptime <- ptime[valid]

# Fit GAMs using tradeSeq; use 6 knots default
k <- 6
counts <- assays(sce)[[1]]
# construct design: pseudotime only (simplest case)
cellWeights <- matrix(1, nrow=1, ncol=ncol(sce))
rownames(cellWeights) <- "W"
rownames(sce) <- rownames(counts)

# Use tradeSeq wrapper: fitGAM
tryCatch({
  w <- fitGAM(counts = counts, pseudotime = matrix(ptime, ncol=1), nknots = k)
  res <- associationTest(w)
  res_df <- as.data.frame(res)
  res_df <- res_df[order(res_df$waldStat, decreasing=TRUE), ]
  write.csv(res_df, file=file.path(opt$output, paste0("pseudotime_deg_", opt$pseudotime_col, ".csv")), row.names = TRUE)
  message("✅ Pseudotime DE results written to ", opt$output)
}, error = function(e) {
  stop("❌ tradeSeq fit failed: ", e$message)
})
