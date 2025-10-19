#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(SPOTlight)
  library(NMF)
  library(dplyr)
})

option_list <- list(
  make_option(c("-r","--reference"), type="character", help="Reference single-cell Seurat .rds with cell type in metadata (e.g., 'cell_type')"),
  make_option(c("-s","--spatial"), type="character", help="Spatial Seurat .rds to deconvolve"),
  make_option(c("-l","--label_col"), type="character", default="cell_type", help="Reference metadata column with labels"),
  make_option(c("-o","--output"), type="character", help="Output Seurat .rds (with spot-level proportions in metadata or assay)"),
  make_option(c("-n","--n_topics"), type="integer", default=300, help="Number of NMF topics (default:300)")
)
opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$reference) || is.null(opt$spatial) || is.null(opt$output)) stop("❌ Provide --reference, --spatial, and --output")

ref <- readRDS(opt$reference)
sp <- readRDS(opt$spatial)
if (!(opt$label_col %in% colnames(ref@meta.data))) stop(paste0("❌ Label column '", opt$label_col, "' not found in reference"))

message("🔎 Running SPOTlight deconvolution...")
set.seed(123)
spotlight_ls <- spotlight_deconvolution(
  se_sc = ref,
  counts_spatial = GetAssayData(sp, assay = DefaultAssay(sp), slot = "counts"),
  clust_vr = opt$label_col,
  cluster_markers = NULL,
  cl_n = NULL,
  hvg = NULL,
  ntop = opt$n_topics,
  method = "nsNMF",
  verbose = TRUE
)

# spotlight_deconvolution returns a list with 'mat' (cell-type proportions per spot)
proportions <- spotlight_ls$mat
# attach proportions to spatial Seurat as metadata columns (prefix SPOTLIGHT_)
for (ct in colnames(proportions)) {
  colname <- paste0("SPOTLIGHT_", ct)
  sp@meta.data[[colname]] <- proportions[, ct]
}
# also add as assay
sp[["spotlight_props"]] <- CreateAssayObject(data = t(proportions))

saveRDS(sp, opt$output)
message("✅ SPOTlight deconvolution complete. Saved to ", opt$output)
