#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(rmarkdown)
  library(Seurat)
  library(ggplot2)
})

option_list <- list(
  make_option(c("-i","--input"), type="character", help="Input annotated Seurat .rds (expects predicted_cell_type or ensemble_label)"),
  make_option(c("-o","--output_html"), type="character", help="Output HTML report path"),
  make_option(c("-t","--title"), type="character", default="Annotation Report", help="Report title")
)
opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$input) || is.null(opt$output_html)) stop("❌ Provide --input and --output_html")

seurat_obj_path <- normalizePath(opt$input)
output_html <- normalizePath(opt$output_html, mustWork = FALSE)
rmd_path <- tempfile(fileext = ".Rmd")

# Write RMarkdown document (contains code chunks)
cat(file = rmd_path, "
---
title: '`r params$title`'
output:
  html_document:
    toc: true
    toc_depth: 2
params:
  input_rds: NULL
  title: 'Annotation Report'
---

```{r setup, include=FALSE}
library(Seurat)
library(ggplot2)
library(pheatmap)
library(htmlwidgets)
so <- readRDS(params$input_rds)
ann_col <- ifelse('ensemble_label' %in% colnames(so@meta.data),'ensemble_label','predicted_cell_type')
```

# Interactive UMAP colored by annotation

This section shows an interactive UMAP (if present) colored by the chosen annotation column.

```{r umap_plot, echo=FALSE, message=FALSE, warning=FALSE}
if (!('umap' %in% names(so@reductions))) {
  if ('pca' %in% names(so@reductions)) {
    so <- tryCatch({ RunUMAP(so, reduction='pca', dims=1:20, verbose = FALSE) }, error = function(e) so)
  }
}
p <- DimPlot(so, reduction='umap', group.by=ann_col, label=TRUE, repel=TRUE) + ggtitle('Annotated UMAP')
print(p)
```

# Confusion heatmap (if gold annotation provided)

If the Seurat object metadata contains a `label_gold` column, a confusion heatmap between gold and predicted labels will be produced.

```{r confusion, echo=FALSE, message=FALSE, warning=FALSE}
if ('label_gold' %in% colnames(so@meta.data)) {
  ann <- data.frame(cell_id = colnames(so), label_gold = so@meta.data$label_gold, label_pred = so@meta.data[[ann_col]])
  tbl <- table(ann$label_gold, ann$label_pred)
  if (nrow(tbl) > 0 && ncol(tbl) > 0) {
    pheatmap::pheatmap(as.matrix(tbl), main = 'Confusion: gold vs predicted', fontsize = 10)
  } else {
    cat('Confusion table is empty or malformed.')
  }
} else {
  cat('No gold labels found in metadata; provide gold labels as metadata column \"label_gold\" for confusion heatmap.')
}
```

# Summary tables and downloads

The report includes a simple table of counts per predicted label and a link to download the table.

```{r summary_table, echo=FALSE}
counts_df <- as.data.frame(table(so@meta.data[[ann_col]]))
colnames(counts_df) <- c('label', 'count')
print(counts_df)
```

")

# Render the Rmd to HTML
render(input = rmd_path, output_file = output_html, params = list(input_rds = seurat_obj_path, title = opt$title), quiet = FALSE)
message('✅ Interactive HTML report written to ', output_html)
