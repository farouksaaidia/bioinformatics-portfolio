#!/usr/bin/env Rscript
# run_clustering_seurat.R
# Robust graph-based clustering wrapper using Seurat (Louvain or Leiden)
#
# Usage:
#   Rscript run_clustering_seurat.R -i input.rds -o output.rds -m leiden -r 0.5
#   Rscript run_clustering_seurat.R -d input_folder/ -o out_folder/ -m louvain -r 0.4,0.6,0.8
#
# Supports single RDS input or directory of RDS objects (batch).
# Produces: saved Seurat object(s) with clustering metadata 'seurat_clusters' and 'cluster_resolution_<r>' when multiple.
suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(dplyr)
})

option_list <- list(
  make_option(c("-i","--input"), type="character", default=NULL, help="Input Seurat RDS file (single)"),
  make_option(c("-d","--dir"), type="character", default=NULL, help="Input folder containing .rds Seurat objects (batch)"),
  make_option(c("-o","--output"), type="character", help="Output file (single) or output folder (batch)"),
  make_option(c("-m","--method"), type="character", default="leiden", help="Clustering method: 'leiden' or 'louvain' [default: leiden]"),
  make_option(c("-r","--resolutions"), type="character", default="0.5", help="Comma-separated resolution(s), e.g. '0.4,0.6'"),
  make_option(c("-n","--npcs"), type="integer", default=30, help="Number of PCs to use [default:30]"),
  make_option(c("--verbose"), action="store_true", default=FALSE, help="Verbose messages")
)
opt <- parse_args(OptionParser(option_list=option_list))

log <- function(...) {
  msg <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", ...)
  cat(msg, "\n")
}

if (is.null(opt$input) && is.null(opt$dir)) stop("Either --input or --dir must be provided.")
if (is.null(opt$output)) stop("--output must be provided.")

res_list <- as.numeric(unlist(strsplit(opt$resolutions,",")))
if (any(is.na(res_list))) stop("Invalid resolution list provided.")

process_one <- function(infile, outfile_base) {
  log("Reading", infile)
  if (!file.exists(infile)) { stop("Input file not found: ", infile) }
  seurat_obj <- readRDS(infile)
  if (!"pca" %in% names(seurat_obj@reductions)) {
    log("PCA not found. Running RunPCA with variable features.")
    seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(seurat_obj), verbose = FALSE)
  }
  # Build neighbor graph
  seurat_obj <- FindNeighbors(seurat_obj, dims = 1:opt$npcs, verbose = FALSE)
  for (r in res_list) {
    if (tolower(opt$method) == "leiden") {
      seurat_obj <- FindClusters(seurat_obj, resolution = r, algorithm = 4, verbose = FALSE) # 4 = Leiden
    } else {
      seurat_obj <- FindClusters(seurat_obj, resolution = r, algorithm = 1, verbose = FALSE) # 1 = Louvain
    }
    # rename default column to include resolution
    colname <- paste0("cluster_resolution_", r)
    seurat_obj@meta.data[[colname]] <- Idents(seurat_obj)
    if (opt$verbose) log("Assigned clusters for resolution", r)
  }
  # keep last result as seurat_clusters for downstream compatibility
  Idents(seurat_obj) <- seurat_obj@meta.data[[paste0("cluster_resolution_", res_list[length(res_list)])]]
  seurat_obj@meta.data$seurat_clusters <- Idents(seurat_obj)
  saveRDS(seurat_obj, paste0(outfile_base, ".rds"))
  log("Wrote clustered object to", paste0(outfile_base, ".rds"))
}

# Single input
if (!is.null(opt$input)) {
  out_base <- ifelse(grepl("\\.rds$", opt$output, ignore.case=TRUE), sub("\\.rds$","",opt$output), file.path(dirname(opt$output), sub("\\.rds$","", basename(opt$input))))
  process_one(opt$input, out_base)
} else {
  # batch mode
  if (!dir.exists(opt$dir)) stop("Input directory not found: ", opt$dir)
  if (!dir.exists(opt$output)) dir.create(opt$output, recursive = TRUE)
  files <- list.files(opt$dir, pattern="\\.rds$", full.names = TRUE)
  if (length(files) == 0) stop("No .rds files found in input directory.")
  for (f in files) {
    basename_noext <- sub("\\.rds$","", basename(f))
    out_base <- file.path(opt$output, basename_noext)
    process_one(f, out_base)
  }
}

log("All done.")
