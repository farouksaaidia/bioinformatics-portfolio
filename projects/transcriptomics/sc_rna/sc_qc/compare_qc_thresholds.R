#!/usr/bin/env Rscript
# compare_qc_thresholds.R
# Sweeps min_genes and max_pct_mt thresholds and reports retained cell counts.
#
# Usage examples:
# Rscript compare_qc_thresholds.R --input sample.rds --outdir ./qc_thresholds --minGenes "100,200,500" --maxMT "5,10,15"
# Rscript compare_qc_thresholds.R --input counts.csv --outdir ./qc_thresholds --minGenes "100,200" --maxMT "5,10"

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(ggplot2)
  library(reshape2)
})

# Helper logging
log <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), ..., "\n")

# Parse arguments
option_list = list(
  make_option(c("-i", "--input"), type="character", help="Input Seurat RDS file (.rds/.rds) or counts matrix CSV (.csv). REQUIRED"),
  make_option(c("-o", "--outdir"), type="character", default="qc_thresholds", help="Output directory [default %default]"),
  make_option(c("--minGenes"), type="character", default="200,500,1000", help="Comma-separated min genes list [default %default]"),
  make_option(c("--maxMT"), type="character", default="5,10,15", help="Comma-separated max mitochondrial % list [default %default]"),
  make_option(c("-v", "--verbose"), action="store_true", default=FALSE, help="Verbose output")
)
opt = parse_args(OptionParser(option_list=option_list))

if (is.null(opt$input)) {
  stop("ERROR: --input is required. See --help for usage.")
}

# Create output dir
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

# Read input
input_path <- opt$input
log("Input:", input_path)
seu <- NULL

if (grepl("\\.rds$|\\.RDS$|\\.rda$|\\.RDA$", input_path)) {
  log("Reading Seurat object from RDS...")
  tryCatch({
    seu <- readRDS(input_path)
    if (!"Seurat" %in% class(seu)) {
      stop("Input RDS is not a Seurat object.")
    }
  }, error = function(e) { stop("Failed to read Seurat RDS: ", e$message) })
} else if (grepl("\\.csv$|\\.CSV$", input_path)) {
  log("Reading count matrix CSV and building Seurat object...")
  tryCatch({
    counts <- as.matrix(read.csv(input_path, row.names = 1, check.names = FALSE))
    seu <- CreateSeuratObject(counts = counts)
  }, error = function(e) { stop("Failed to read counts CSV: ", e$message) })
} else {
  stop("Unsupported input file type. Provide a Seurat RDS or counts CSV.")
}

# Ensure basic QC columns exist; if not, compute common metrics
if (!"nFeature_RNA" %in% colnames(seu@meta.data)) {
  log("Computing per-cell QC metrics (nFeature_RNA, nCount_RNA, percent.mt) ...")
  # compute percent.mt if MT genes exist
  mt_genes <- grep(pattern = "^MT-|^MT\\.", x = rownames(seu), value = TRUE, ignore.case = TRUE)
  if (length(mt_genes) > 0) {
    seu <- PercentageFeatureSet(seu, pattern = "^MT-|^MT\\.", col.name = "percent.mt")
  } else {
    # create placeholder percent.mt (NA) if no MT genes
    seu@meta.data$percent.mt <- NA
  }
}

# Parse threshold lists
minGenes <- as.integer(unlist(strsplit(opt$minGenes, ",")))
maxMT    <- as.numeric(unlist(strsplit(opt$maxMT, ",")))

log("minGenes:", paste(minGenes, collapse = ", "))
log("maxMT:", paste(maxMT, collapse = ", "))

results <- data.frame()
for (mg in minGenes) {
  for (mm in maxMT) {
    # Build subset expression
    # allow NA percent.mt -> treat NA as pass (i.e., don't filter those cells by mt if NA)
    keep_idx <- rep(TRUE, ncol(seu))
    if ("nFeature_RNA" %in% colnames(seu@meta.data)) {
      keep_idx <- keep_idx & (seu@meta.data$nFeature_RNA > mg)
    }
    if ("percent.mt" %in% colnames(seu@meta.data)) {
      # if percent.mt is NA for some cells, treat NA as FALSE for >mm (i.e. keep them)
      mtvals <- seu@meta.data$percent.mt
      mt_ok <- ifelse(is.na(mtvals), TRUE, mtvals < mm)
      keep_idx <- keep_idx & mt_ok
    }
    retained <- sum(keep_idx, na.rm = TRUE)
    results <- rbind(results, data.frame(minGenes = mg, maxMT = mm, retained = retained))
    if (opt$verbose) log(sprintf("Threshold minGenes=%d maxMT=%g -> retained %d cells", mg, mm, retained))
  }
}

# Save table
out_table <- file.path(opt$outdir, "qc_thresholds_retained_cells.csv")
write.csv(results, out_table, row.names = FALSE)
log("Wrote thresholds table to", out_table)

# Create heatmap
library(viridis)
mat <- reshape2::acast(results, minGenes ~ maxMT, value.var = "retained")
png(file.path(opt$outdir, "qc_thresholds_heatmap.png"), width = 1000, height = 800)
par(mar = c(5,5,4,2))
image(x = as.numeric(colnames(mat)), y = as.numeric(rownames(mat)), z = t(mat[nrow(mat):1, , drop = FALSE]),
      xlab = "max %MT", ylab = "min Genes", main = "Cells retained per threshold", col = viridis::viridis(100))
# add text values
nr <- nrow(mat); nc <- ncol(mat)
xs <- as.numeric(colnames(mat)); ys <- rev(as.numeric(rownames(mat)))
for (i in seq_len(nr)) {
  for (j in seq_len(nc)) {
    text(xs[j], ys[i], labels = mat[nr - i + 1, j], cex = 0.9)
  }
}
dev.off()
log("Wrote heatmap to", file.path(opt$outdir, "qc_thresholds_heatmap.png"))

log("Done.")
