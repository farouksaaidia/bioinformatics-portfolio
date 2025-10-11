#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Seurat)
  library(optparse)
  library(dplyr)
  library(jsonlite)
})

option_list <- list(
  make_option(c("-i","--input"), type="character", help="Input annotated Seurat .rds (expects predicted_cell_type or ensemble_label)"),
  make_option(c("-o","--output"), type="character", help="Output Seurat .rds with ontology mappings"),
  make_option(c("-m","--mapping"), type="character", default=NULL, help="Optional JSON mapping file mapping label->{ontology_id, parent_id, preferred_name}. If NULL, uses a small built-in mapping")
)
opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$input) || is.null(opt$output)) stop("❌ Provide --input and --output")

seurat_obj <- readRDS(opt$input)
meta_key <- ifelse("ensemble_label" %in% colnames(seurat_obj@meta.data), "ensemble_label",
                   ifelse("predicted_cell_type" %in% colnames(seurat_obj@meta.data), "predicted_cell_type", NA))
if (is.na(meta_key)) stop("❌ No known annotation column found (expected ensemble_label or predicted_cell_type)")

# load mapping
if (!is.null(opt$mapping) && file.exists(opt$mapping)) {
  mapping <- fromJSON(opt$mapping, simplifyDataFrame = TRUE)
} else {
  # small example mapping; expand as needed
  mapping <- data.frame(
    label = c("B cell", "T cell", "CD4 T cell", "CD8 T cell", "monocyte"),
    ontology_id = c("CL:0000789", "CL:0000084", "CL:0000624", "CL:0000623", "CL:0000576"),
    parent_id = c("CL:0000236","CL:0000236","CL:0000084","CL:0000084","CL:0000576"),
    preferred_name = c("B cell","T cell","CD4-positive, alpha-beta T cell","CD8-positive, alpha-beta T cell","monocyte"),
    stringsAsFactors = FALSE
  )
}

# join metadata
meta <- seurat_obj@meta.data
meta$label_lookup <- meta[[meta_key]]
map_df <- as.data.frame(mapping)
mapped <- left_join(meta, map_df, by = c("label_lookup" = "label"))

# attach columns
seurat_obj@meta.data$cell_ontology_id <- mapped$ontology_id
seurat_obj@meta.data$cell_ontology_parent <- mapped$parent_id
seurat_obj@meta.data$cell_ontology_preferred_name <- mapped$preferred_name

saveRDS(seurat_obj, opt$output)
message("✅ Ontology mappings added and saved to ", opt$output)
