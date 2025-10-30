#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(dplyr)
})

option_list <- list(
  make_option(c("-i","--counts"), type="character", help="Raw integer counts CSV/TSV (genes x samples)"),
  make_option(c("-m","--metadata"), type="character", help="Sample metadata TSV with columns: sample,batch,condition (batch required)"),
  make_option(c("-o","--outdir"), type="character", help="Output directory"),
  make_option(c("--method"), type="character", default="combat_seq", help="Method: combat_seq (default) or removeBatchEffect")
)
opt <- parse_args(OptionParser(option_list=option_list))

if (is.null(opt$counts) || is.null(opt$metadata) || is.null(opt$outdir)) stop("Provide --counts, --metadata, --outdir")
if (!file.exists(opt$counts)) stop("Counts file not found")
if (!file.exists(opt$metadata)) stop("Metadata file not found")
dir.create(opt$outdir, recursive=TRUE, showWarnings=FALSE)

counts_df <- read_delim(opt$counts, delim="\t", col_types = cols())
if (!"gene_id" %in% colnames(counts_df)) stop("Counts must have gene_id column")
mat <- as.matrix(counts_df[ , setdiff(colnames(counts_df), "gene_id")])
rownames(mat) <- counts_df$gene_id
mat <- round(mat)

meta <- read_delim(opt$metadata, delim="\t", col_types = cols())
if (!all(c("sample","batch") %in% colnames(meta))) stop("Metadata must include columns: sample and batch")
batch <- meta$batch[match(colnames(mat), meta$sample)]
if (any(is.na(batch))) stop("Some samples in counts are missing in metadata")

method <- tolower(opt$method)
if (method == "combat_seq") {
  if (!requireNamespace("sva", quietly=TRUE)) {
    stop("Package 'sva' required for ComBat_seq. Install Bioconductor 'sva'.")
  }
  library(sva)
  message("Running ComBat_seq on integer counts...")
  corrected <- ComBat_seq(counts = mat, batch = batch)
  write.csv(as.data.frame(corrected), file.path(opt$outdir, "counts_combatseq_corrected.csv"), row.names=TRUE)
  message("ComBat_seq output saved.")
} else if (method == "removebatcheffect") {
  if (!requireNamespace("limma", quietly=TRUE)) {
    stop("Package 'limma' required for removeBatchEffect fallback.")
  }
  library(limma)
  message("Converting to logCPM and running removeBatchEffect...")
  if (!requireNamespace("edgeR", quietly=TRUE)) stop("edgeR required for logCPM conversion")
  library(edgeR)
  dge <- DGEList(counts=mat)
  logcpm <- cpm(dge, log=TRUE, prior.count=1)
  adj <- removeBatchEffect(logcpm, batch = batch)
  write.csv(as.data.frame(adj), file.path(opt$outdir, "counts_logcpm_removeBatchEffect.csv"), row.names=TRUE)
  message("removeBatchEffect output saved.")
} else {
  stop("Unknown method: choose combat_seq or removeBatchEffect")
}
