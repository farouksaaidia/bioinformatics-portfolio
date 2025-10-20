#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(spdep)
  library(Matrix)
  library(dplyr)
})

option_list <- list(
  make_option(c("-i","--input"), type="character", help="Input Seurat .rds (spatial)"),
  make_option(c("-g","--genes"), type="character", default=NULL, help="Comma-separated genes to test (default: top variable genes)"),
  make_option(c("-k","--k"), type="integer", default=6, help="k nearest neighbors for spatial weights (default 6)"),
  make_option(c("-o","--output"), type="character", help="Output CSV path for Moran's I results")
)
opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$input) || is.null(opt$output)) stop("❌ Provide --input and --output")
so <- readRDS(opt$input)

# extract coords
coords <- NULL
if (length(so@images) > 0 && !is.null(so@images[[1]]@coordinates)) {
  coords <- as.data.frame(so@images[[1]]@coordinates)
} else if (all(c("imagecol","imagerow") %in% colnames(so@meta.data))) {
  coords <- so@meta.data[, c("imagecol","imagerow")]
} else {
  stop("❌ Spatial coordinates not found in Seurat object")
}

# select genes
genes <- NULL
if (!is.null(opt$genes)) {
  genes <- unlist(strsplit(opt$genes, ","))
  genes <- trimws(genes)
} else {
  # default: top variable features
  if ("RNA" %in% Assays(so)) {
    genes <- head(VariableFeatures(so), 500)
  } else {
    genes <- rownames(so)[1:500]
  }
}
genes <- intersect(genes, rownames(so))
message("ℹ️ Testing Moran's I for ", length(genes), " genes")

# build neighbor list
coords_mat <- as.matrix(coords)
knn <- FNN::get.knn(coords_mat, k = opt$k)
nb_list <- vector("list", nrow(coords_mat))
for (i in seq_len(nrow(coords_mat))) nb_list[[i]] <- knn$nn.index[i,]
lw <- nb2listw(neig2nb(nb_list), style="W", zero.policy=TRUE)

results <- data.frame(gene=character(0), moran_I=double(0), p_value=double(0), stringsAsFactors=FALSE)
for (g in genes) {
  vals <- as.numeric(GetAssayData(so, slot="data")[g, ])
  if (all(is.na(vals)) || sd(vals, na.rm=TRUE) == 0) next
  test <- try(moran.test(vals, lw, zero.policy=TRUE), silent=TRUE)
  if (!inherits(test, "try-error")) {
    results <- rbind(results, data.frame(gene=g, moran_I = test$estimate[["Moran I statistic"]], p_value = test$p.value))
  }
}

if (nrow(results) == 0) stop("❌ No Moran's I results computed")
results <- results %>% arrange(p_value)
write.csv(results, opt$output, row.names=FALSE)
message("✅ Moran's I results written to ", opt$output)
