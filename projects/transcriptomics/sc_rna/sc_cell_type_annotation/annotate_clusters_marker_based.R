#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(optparse)
})

option_list <- list(
  make_option(c("-i", "--input"), type="character", help="Input clustered Seurat .rds file"),
  make_option(c("-m", "--markers"), type="character", help="CSV file of known marker genes (columns: gene, cell_type)"),
  make_option(c("-o", "--output"), type="character", help="Output annotated Seurat .rds file")
)
opt <- parse_args(OptionParser(option_list=option_list))

if (is.null(opt$input) || is.null(opt$markers) || is.null(opt$output)) {
  stop("Missing arguments: provide --input, --markers, and --output")
}

cat("📂 Loading clustered Seurat object...\n")
seurat_obj <- readRDS(opt$input)
markers <- read.csv(opt$markers)

cat("🧬 Annotating clusters based on marker enrichment...\n")
cluster_markers <- FindAllMarkers(seurat_obj, only.pos=TRUE)
annotations <- cluster_markers %>%
  group_by(cluster) %>%
  summarize(cell_type = markers$cell_type[match(gene[1], markers$gene)])

seurat_obj$cell_type <- plyr::mapvalues(seurat_obj$seurat_clusters, 
                                        from=annotations$cluster, to=annotations$cell_type)

saveRDS(seurat_obj, file=opt$output)
cat("✅ Annotation saved to", opt$output, "\n")
