#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Seurat)
  library(optparse)
  library(dplyr)
})

option_list <- list(
  make_option(c("-r", "--reference"), type="character", help="Reference Seurat .rds with cell type in metadata (e.g. ref$cell_type)"),
  make_option(c("-q", "--query"), type="character", help="Query Seurat .rds to annotate"),
  make_option(c("-o", "--output"), type="character", help="Output annotated Seurat .rds"),
  make_option(c("-l", "--ref_label"), type="character", default="cell_type", help="Metadata column in reference with labels (default: cell_type)"),
  make_option(c("-k", "--k"), type="integer", default=20, help="k.anchor for FindTransferAnchors (default: 20)"),
  make_option(c("-d", "--dims"), type="integer", default=30, help="Number of dims to use in anchors and TransferData (default: 30)")
)
opt <- parse_args(OptionParser(option_list=option_list))

if (is.null(opt$reference) || is.null(opt$query) || is.null(opt$output)) {
  stop("❌ Provide --reference, --query, and --output")
}

message("📥 Loading reference and query Seurat objects...")
ref <- readRDS(opt$reference)
qry <- readRDS(opt$query)

if (!(opt$ref_label %in% colnames(ref@meta.data))) {
  stop(paste0("❌ Reference metadata column '", opt$ref_label, "' not found"))
}

message("🔗 Running FindTransferAnchors (k=", opt$k, ", dims=", opt$dims, ")...")
anchors <- FindTransferAnchors(reference = ref, query = qry, dims = 1:opt$dims, k.anchor = opt$k)

message("🎯 Transferring labels...")
pred <- TransferData(anchorset = anchors, refdata = ref@meta.data[[opt$ref_label]], dims = 1:opt$dims)
qry <- AddMetaData(qry, metadata = pred)

# Standardize metadata key
if ("predicted.id" %in% colnames(pred)) {
  meta_key <- "predicted.id"
} else {
  meta_key <- colnames(pred)[1]
}
names(qry@meta.data)[names(qry@meta.data) == meta_key] <- "predicted_cell_type_transfer"

message("💾 Saving annotated query to ", opt$output)
saveRDS(qry, opt$output)
message("✅ Done.")
