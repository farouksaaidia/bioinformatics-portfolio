#!/usr/bin/env bash
set -euo pipefail
# Export images and embeddings for downstream figure assembly
# Usage: export_visual_assets.sh -i input_seurat.rds -o outdir

usage() {
  echo "Usage: $0 -i <input_seurat.rds> -o <outdir>"
  exit 1
}

INPUT=""
OUTDIR=""

while getopts "i:o:" opt; do
  case $opt in
    i) INPUT=$OPTARG ;;
    o) OUTDIR=$OPTARG ;;
    *) usage ;;
  esac
done

if [[ -z "$INPUT" || -z "$OUTDIR" ]]; then
  usage
fi

mkdir -p "$OUTDIR"

Rscript - "$INPUT" "$OUTDIR" <<'RS'
args <- commandArgs(trailingOnly=TRUE)
input <- args[1]
outdir <- args[2]

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(readr)
  library(patchwork)
})

so <- readRDS(input)

# Export embeddings
if ("umap" %in% names(so@reductions)) {
  write_csv(as.data.frame(Embeddings(so, "umap")), file.path(outdir, "umap_embeddings.csv"))
}
if ("pca" %in% names(so@reductions)) {
  write_csv(as.data.frame(Embeddings(so, "pca")), file.path(outdir, "pca_embeddings.csv"))
}

# Export spatial coordinates if available
if (length(so@images) > 0 && !is.null(so@images[[1]]@coordinates)) {
  coords <- as.data.frame(so@images[[1]]@coordinates)
  write_csv(coords, file.path(outdir, "spatial_coordinates.csv"))
} else if (all(c("imagecol","imagerow") %in% colnames(so@meta.data))) {
  write_csv(so@meta.data[, c("imagecol","imagerow")], file.path(outdir, "spatial_coordinates.csv"))
}

# Export top features as PNGs
features <- head(VariableFeatures(so), 6)
for (f in features) {
  p <- tryCatch(SpatialFeaturePlot(so, features=f) + ggtitle(f),
                error=function(e) ggplot() + ggtitle(paste("Missing", f)))
  ggsave(filename=file.path(outdir, paste0("feature_", f, ".png")),
         plot=p, width=6, height=5, dpi=300)
}

cat("✅ Visual assets exported to", outdir, "\n")
RS

echo "✅ export_visual_assets completed: $OUTDIR"
