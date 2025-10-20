#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(dplyr)
})

option_list <- list(
  make_option(c("-i","--inputs"), type="character", help="Comma-separated list of input Seurat .rds files to merge"),
  make_option(c("-n","--names"), type="character", default=NULL, help="Optional comma-separated sample names (must match inputs order)"),
  make_option(c("-o","--output"), type="character", help="Output merged Seurat .rds file")
)
opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$inputs) || is.null(opt$output)) stop("❌ Provide --inputs and --output")

inputs <- unlist(strsplit(opt$inputs, ","))
inputs <- trimws(inputs)
if (!all(file.exists(inputs))) stop("❌ One or more input files not found")

sample_names <- NULL
if (!is.null(opt$names)) {
  sample_names <- unlist(strsplit(opt$names, ","))
  sample_names <- trimws(sample_names)
  if (length(sample_names) != length(inputs)) stop("❌ --names length must equal number of inputs")
}

objs <- lapply(seq_along(inputs), function(i) {
  x <- readRDS(inputs[i])
  sn <- if (!is.null(sample_names)) sample_names[i] else paste0("sample", i)
  x$sample_id <- sn
  # ensure unique cell/spot names by prefix
  colnames(x) <- paste0(sn, "_", colnames(x))
  return(x)
})

message("🔗 Merging ", length(objs), " Seurat objects...")
merged <- Reduce(function(a, b) merge(a, y=b), objs)

message("💾 Saving merged object: ", opt$output)
saveRDS(merged, file = opt$output)
message("✅ Done.")
