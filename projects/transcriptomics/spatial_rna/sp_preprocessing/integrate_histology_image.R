#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(jpeg)
  library(png)
  library(imager)
  library(dplyr)
  library(glue)
})

option_list <- list(
  make_option(c("-i","--input"), type="character", help="Input Seurat .rds with spatial coordinates"),
  make_option(c("-img","--image"), type="character", help="High-resolution histology image (PNG or JPEG)"),
  make_option(c("-o","--output"), type="character", help="Output Seurat .rds (updated)"),
  make_option(c("-r","--radius"), type="integer", default=5, help="Radius in pixels to sample around each spot (default:5)")
)
opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$input) || is.null(opt$image) || is.null(opt$output)) stop("❌ Provide --input, --image, --output")

so <- readRDS(opt$input)
img_path <- normalizePath(opt$image, mustWork = TRUE)
message(glue("📸 Loading image: {img_path}"))
# load image robustly
img <- tryCatch({
  if (grepl("\\.png$", img_path, ignore.case = TRUE)) {
    png::readPNG(img_path)
  } else if (grepl("\\.(jpg|jpeg)$", img_path, ignore.case = TRUE)) {
    jpeg::readJPEG(img_path)
  } else {
    stop("Unsupported image format (use PNG or JPEG)")
  }
}, error = function(e) stop("❌ Failed to read image: ", e$message))

# convert to grayscale intensity if necessary
if (length(dim(img)) == 3) {
  img_gray <- 0.2989 * img[,,1] + 0.5870 * img[,,2] + 0.1140 * img[,,3]
} else {
  img_gray <- img
}
# store image in @images if not already
img_name <- "histology_custom"
if (!img_name %in% names(so@images)) {
  so@images[[img_name]] <- new(Class = "SlideSeq", image = img) # best-effort; user may want to adjust
}

# compute mean intensity within a small radius around each spatial coordinate (if coords exist)
coords <- tryCatch({
  if (!is.null(so@images) && length(so@images) > 0 && !is.null(so@images[[1]]@coordinates)) {
    so@images[[1]]@coordinates
  } else if ("imagecol" %in% colnames(so@meta.data)) {
    so@meta.data[, c("imagecol", "imagerow")]
  } else {
    NULL
  }
}, error = function(e) NULL)

if (!is.null(coords)) {
  message("🔎 Computing image-derived mean intensity per spot...")
  coords <- as.data.frame(coords)
  intensities <- numeric(nrow(coords))
  h <- nrow(img_gray); w <- ncol(img_gray)
  for (i in seq_len(nrow(coords))) {
    x <- round(coords[i, "imagecol"]); y <- round(coords[i, "imagerow"])
    x1 <- max(1, x - opt$radius); x2 <- min(w, x + opt$radius)
    y1 <- max(1, y - opt$radius); y2 <- min(h, y + opt$radius)
    patch <- img_gray[y1:y2, x1:x2]
    intensities[i] <- mean(patch, na.rm = TRUE)
  }
  so@meta.data$histology_mean_intensity <- intensities
} else {
  message("⚠️ Spatial coordinates for image sampling not found; skipping intensity computation.")
}

saveRDS(so, opt$output)
message("✅ Image integrated and object saved to ", opt$output)
