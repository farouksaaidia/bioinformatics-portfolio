#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(dplyr)
})

option_list <- list(
  make_option(c("-i","--counts"), type="character", help="Counts matrix TSV/CSV (genes x samples). First column gene_id"),
  make_option(c("-l","--lengths"), type="character", default=NULL, help="Optional gene lengths TSV/CSV with columns: gene_id,length"),
  make_option(c("-m","--methods"), type="character", default="tpm,cpm,tmm", help="Comma-separated methods to compute: tpm,cpm,rpkm,tmm"),
  make_option(c("-o","--outdir"), type="character", help="Output directory")
)
opt <- parse_args(OptionParser(option_list=option_list))

if (is.null(opt$counts) || is.null(opt$outdir)) {
  stop("Provide --counts and --outdir")
}
if (!file.exists(opt$counts)) stop("Counts file not found: ", opt$counts)
dir.create(opt$outdir, recursive=TRUE, showWarnings=FALSE)

methods <- tolower(trimws(unlist(strsplit(opt$methods, ","))))

message("Loading counts: ", opt$counts)
counts_df <- read_delim(opt$counts, delim = "\t", col_types = cols())

if (!"gene_id" %in% colnames(counts_df)) stop("Counts file must have 'gene_id' first column")

gene_ids <- counts_df$gene_id
mat <- as.matrix(counts_df[ , setdiff(colnames(counts_df), "gene_id")])
rownames(mat) <- gene_ids

# Ensure numeric matrix
mat <- apply(mat, 2, as.numeric)
rownames(mat) <- gene_ids

# helper functions
safe_write <- function(x, fname) {
  write.csv(as.data.frame(x), fname, row.names = TRUE)
  message("Wrote: ", fname)
}

if ("cpm" %in% methods) {
  message("Computing CPM...")
  libsize <- colSums(mat, na.rm=TRUE)
  cpm <- sweep(mat, 2, libsize/1e6, "/")
  safe_write(cpm, file.path(opt$outdir, "counts_cpm.csv"))
}

if ("rpkm" %in% methods || "tpm" %in% methods) {
  if (is.null(opt$lengths)) {
    warning("Gene lengths not provided; skipping RPKM/TPM")
  } else {
    if (!file.exists(opt$lengths)) stop("Lengths file not found: ", opt$lengths)
    lengths_df <- read_delim(opt$lengths, delim="\t", col_types = cols())
    if (!all(c("gene_id","length") %in% colnames(lengths_df))) stop("Lengths file must contain columns: gene_id,length")
    lengths <- lengths_df$length[match(gene_ids, lengths_df$gene_id)]
    if (any(is.na(lengths))) stop("Some gene lengths are missing or gene_id mismatch")
    message("Computing RPKM/TPM using provided gene lengths")
    # RPKM
    if ("rpkm" %in% methods) {
      rpkm <- sweep(mat, 1, lengths/1000, "/")
      libsize <- colSums(mat, na.rm=TRUE)
      rpkm <- sweep(rpkm, 2, libsize/1e6, "/")
      safe_write(rpkm, file.path(opt$outdir, "counts_rpkm.csv"))
    }
    # TPM
    if ("tpm" %in% methods) {
      rpk <- sweep(mat, 1, lengths/1000, "/")   # reads per kilobase
      per_million <- sweep(rpk, 2, colSums(rpk, na.rm=TRUE)/1e6, "/")
      safe_write(per_million, file.path(opt$outdir, "counts_tpm.csv"))
    }
  }
}

if ("tmm" %in% methods) {
  message("Computing TMM normalization factors (edgeR)...")
  if (!requireNamespace("edgeR", quietly = TRUE)) {
    stop("Package 'edgeR' required for TMM. Install Bioconductor edgeR.")
  }
  library(edgeR)
  dge <- DGEList(counts = round(mat))
  dge <- calcNormFactors(dge, method="TMM")
  tmm_factors <- dge$samples$norm.factors
  write.csv(data.frame(sample=colnames(mat), tmm_factor=tmm_factors), file.path(opt$outdir, "tmm_factors.csv"), row.names=FALSE)
  # produce TMM-adjusted CPM as well
  tmm_cpm <- cpm(dge, normalized.lib.sizes=TRUE)
  write.csv(as.data.frame(tmm_cpm), file.path(opt$outdir, "counts_cpm_tmm.csv"), row.names=TRUE)
  message("TMM factors and adjusted CPM saved")
}

message("All requested normalization methods complete. Outputs in: ", opt$outdir)
