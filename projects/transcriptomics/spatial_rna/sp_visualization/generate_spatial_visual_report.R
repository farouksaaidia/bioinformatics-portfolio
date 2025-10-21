#!/usr/bin/env Rscript
# Generate a combined visual report (PDF) with key spatial plots and summary statistics.
suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(ggplot2)
  library(rmarkdown)
  library(glue)
  library(patchwork)
})

option_list <- list(
  make_option(c("-i","--input"), type="character", help="Input Seurat .rds (spatial)"),
  make_option(c("-o","--output_pdf"), type="character", help="Output PDF file path"),
  make_option(c("-f","--features"), type="character", default=NULL, help="Comma-separated features for highlight")
)
opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$input) || is.null(opt$output_pdf)) stop("❌ Provide --input and --output_pdf")

so <- readRDS(opt$input)
features <- if (!is.null(opt$features)) trimws(unlist(strsplit(opt$features, ","))) else head(VariableFeatures(so), 6)

# Create a temporary Rmd file dynamically and render it to PDF
rmd <- tempfile(fileext = ".Rmd")

rmd_lines <- c(
  "---",
  "title: 'Spatial Visualization Report'",
  "output: pdf_document",
  "---",
  "",
  "```{r setup, include=FALSE}",
  "library(Seurat)",
  "library(ggplot2)",
  "library(patchwork)",
  paste0("so <- readRDS('", opt$input, "')"),
  "```",
  "",
  "# Summary",
  "",
  "**Number of spots:** `r ncol(so)`  ",
  "**Number of features:** `r nrow(so)`",
  "",
  "# Spatial Domains",
  "```{r domains, echo=FALSE, message=FALSE}",
  "if ('seurat_clusters' %in% colnames(so@meta.data)) {",
  "  print(SpatialDimPlot(so, group.by='seurat_clusters', label=TRUE))",
  "} else {",
  "  print('No seurat_clusters column found')",
  "}",
  "```",
  "",
  "# Highlighted Features",
  "```{r features, echo=FALSE, message=FALSE}",
  paste0("features <- c('", paste(features, collapse=\"','\"), "')"),
  "plots <- lapply(features, function(f) {",
  "  if (f %in% rownames(so)) {",
  "    SpatialFeaturePlot(so, features = f) + ggtitle(f)",
  "  } else {",
  "    ggplot() + ggtitle(paste('Missing', f))",
  "  }",
  "})",
  "print(wrap_plots(plots, ncol=2))",
  "```"
)

writeLines(rmd_lines, rmd)

message('🧾 Rendering report to ', opt$output_pdf)
rmarkdown::render(rmd, output_file = opt$output_pdf, quiet = TRUE)
message('✅ Report written to ', opt$output_pdf)
