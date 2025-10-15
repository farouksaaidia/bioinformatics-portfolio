#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(ggplot2)
  library(patchwork)
})

option_list <- list(
  make_option(c("-i","--input"), type="character", help="Input Seurat .rds with pseudotime columns"),
  make_option(c("-p","--pseudotime"), type="character", default="monocle3_pseudotime", help="Pseudotime metadata column to plot"),
  make_option(c("-g","--genes"), type="character", default=NULL, help="Comma-separated genes to plot along pseudotime"),
  make_option(c("-o","--output"), type="character", help="Output directory for plots")
)
opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$input) || is.null(opt$output)) stop("❌ Provide --input and --output")

dir.create(opt$output, showWarnings=FALSE, recursive=TRUE)
so <- readRDS(opt$input)
if (!(opt$pseudotime %in% colnames(so@meta.data))) stop("❌ pseudotime column not found")

# UMAP colored by pseudotime
if ("umap" %in% names(so@reductions)) {
  p1 <- FeaturePlot(so, features = opt$pseudotime, reduction="umap") + ggtitle(paste("UMAP -", opt$pseudotime))
  ggsave(file.path(opt$output, paste0("umap_pseudotime_", opt$pseudotime, ".png")), p1)
} else {
  message("⚠️ UMAP not found; skipping UMAP pseudotime plot")
}

# gene trends along pseudotime
if (!is.null(opt$genes)) {
  genes <- unlist(strsplit(opt$genes, ","))
  # make a simple scatter + smooth per gene
  df_meta <- so@meta.data
  df_meta$cell <- rownames(df_meta)
  for (g in genes) {
    if (g %in% rownames(so)) {
      expr <- FetchData(so, vars = g)[,1]
      df <- data.frame(pseudotime = df_meta[[opt$pseudotime]], expr = expr)
      p <- ggplot(df, aes(x=pseudotime, y=expr)) + geom_point(alpha=0.3) + geom_smooth() + ggtitle(paste(g, "vs", opt$pseudotime))
      ggsave(file.path(opt$output, paste0("gene_trend_", g, ".png")), p)
    } else {
      message("⚠️ Gene not found in object: ", g)
    }
  }
}
message("✅ Trajectory visualizations saved to ", opt$output)
