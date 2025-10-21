#!/usr/bin/env Rscript
# Plot spatial feature overlays (Seurat Spatial objects). Saves PNG/PDF for each feature.
suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(glue)
})

option_list <- list(
  make_option(c("-i","--input"), type="character", help="Input Seurat .rds (spatial, normalized)"),
  make_option(c("-f","--features"), type="character", default=NULL, help="Comma-separated feature list to plot (default: top variable features)"),
  make_option(c("-o","--outdir"), type="character", help="Output directory for plots"),
  make_option(c("--png"), action="store_true", default=TRUE, help="Save PNG files (default TRUE)"),
  make_option(c("--pdf"), action="store_true", default=FALSE, help="Save a combined PDF")
)
opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$input) || is.null(opt$outdir)) stop("❌ Provide --input and --outdir")
dir.create(opt$outdir, recursive=TRUE, showWarnings=FALSE)

so <- readRDS(opt$input)
DefaultAssay(so) <- DefaultAssay(so)

# select features
features <- NULL
if (!is.null(opt$features)) {
  features <- trimws(unlist(strsplit(opt$features, ",")))
} else {
  v <- VariableFeatures(so)
  if (length(v) >= 10) features <- head(v, 10) else features <- head(rownames(so), 10)
}

plots <- list()
message("🔎 Creating spatial feature plots for: ", paste(features, collapse=", "))
for (feat in features) {
  if (!feat %in% rownames(so)) {
    warning(glue("⚠️ Feature not found in object: {feat} — skipping"))
    next
  }
  p <- tryCatch({
    SpatialFeaturePlot(so, features = feat) + ggtitle(feat)
  }, error = function(e) {
    ggplot() + ggtitle(glue("Failed to plot: {feat}"))
  })
  plots[[feat]] <- p
  if (opt$png) {
    png_file <- file.path(opt$outdir, paste0("spatial_feature_", feat, ".png"))
    ggsave(png_file, plot = p, width = 6, height = 5, dpi = 300)
  }
}

if (opt$pdf && length(plots) > 0) {
  pdf(file.path(opt$outdir, "spatial_features_combined.pdf"), width=8, height=10)
  print(wrap_plots(plots, ncol=2))
  dev.off()
}

message("✅ Plots written to ", opt$outdir)
