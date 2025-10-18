#!/usr/bin/env Rscript
# Load 10x Visium data into a Seurat Spatial object, run basic QC, save canonical .rds
suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(dplyr)
  library(glue)
})

option_list <- list(
  make_option(c("-s","--sample_dir"), type="character", help="10x Space Ranger/Visium sample folder (contains 'filtered_feature_bc_matrix' or 'outs')"),
  make_option(c("-n","--sample_name"), type="character", help="Sample name (used for metadata)"),
  make_option(c("-o","--output"), type="character", help="Output Seurat .rds file (recommended extension: .rds)"),
  make_option(c("--image_assay"), type="character", default="image", help="Name to store histology image assay (default: image)")
)
opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$sample_dir) || is.null(opt$sample_name) || is.null(opt$output)) {
  stop("❌ Required: --sample_dir, --sample_name, and --output")
}

sample_dir <- normalizePath(opt$sample_dir, mustWork = TRUE)
out_file <- opt$output

message(glue("📥 Loading Visium data from: {sample_dir}"))
# Try common 10x layouts: directly pointing to output folder or to raw filtered matrix path
seurat_obj <- tryCatch({
  Load10X_Spatial(data.dir = sample_dir, filename = NULL, assay = "Spatial", image = NULL)
}, error = function(e){
  # Try the 'outs' subfolder common to spaceranger
  alt <- file.path(sample_dir, "outs")
  if (dir.exists(alt)) {
    tryCatch({
      Load10X_Spatial(data.dir = alt, assay = "Spatial")
    }, error = function(e2) stop("❌ Failed to Load10X_Spatial: ", e2$message))
  } else {
    stop("❌ Load10X_Spatial failed and no 'outs' folder found.")
  }
})

# Add sample name
seurat_obj$sample_id <- opt$sample_name

# Basic QC metrics
message("🧪 Computing basic QC metrics...")
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
# rename default spatial counts slots if needed
if (!"Spatial" %in% Assays(seurat_obj)) {
  DefaultAssay(seurat_obj) <- Assays(seurat_obj)[1]
}

# Save canonical Seurat object
message(glue("💾 Saving Seurat object to {out_file}"))
saveRDS(seurat_obj, file = out_file)
message("✅ Done.")
