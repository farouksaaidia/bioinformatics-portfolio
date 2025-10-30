#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
})

option_list <- list(
  make_option(c("-c","--counts"), type="character", help="Counts matrix TSV (rows=transcript, cols=samples). First col must be transcript_id and include length if third column 'length' exists"),
  make_option(c("-l","--lengths"), type="character", default=NULL, help="Optional two-column TSV of transcript_id and length (bp). If absent, script expects 'length' column in counts file."),
  make_option(c("-o","--output"), type="character", help="Output directory for normalized matrices"),
  make_option(c("-m","--method"), type="character", default="TPM", help="Normalization method: TPM, CPM, RPKM (default: TPM)")
)
opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$counts) || is.null(opt$output)) stop("Provide --counts and --output")

if (!file.exists(opt$counts)) stop("Counts file not found: ", opt$counts)
dir.create(opt$output, showWarnings=FALSE, recursive=TRUE)

counts <- read.table(opt$counts, header=TRUE, sep="\t", stringsAsFactors=FALSE, check.names=FALSE)
if (!("transcript_id" %in% colnames(counts))) {
  # try first column as transcript id
  colnames(counts)[1] <- "transcript_id"
}
if (!is.null(opt$lengths)) {
  if (!file.exists(opt$lengths)) stop("Lengths file not found: ", opt$lengths)
  lengths <- read.table(opt$lengths, header=TRUE, sep="\t", stringsAsFactors=FALSE)
  if (!all(c("transcript_id","length") %in% colnames(lengths))) stop("Lengths file must have columns: transcript_id,length")
  counts <- left_join(counts, lengths, by="transcript_id")
} else {
  if (!("length" %in% colnames(counts))) stop("No length column found: provide --lengths or ensure counts file has 'length' column")
}

# extract matrix
mat <- as.matrix(counts[, setdiff(colnames(counts), c("transcript_id","length"))])
rownames(mat) <- counts$transcript_id
len <- counts$length

method <- toupper(opt$method)
if (!method %in% c("TPM","CPM","RPKM")) stop("Method must be one of TPM, CPM, RPKM")

# CPM
cpm <- apply(mat, 2, function(x) (x / sum(x)) * 1e6)

# RPKM: reads per kilobase per million
rpkm <- sweep(mat, 1, len/1000, "/")
rpkm <- apply(rpkm, 2, function(x) (x / sum(x)) * 1e6)

# TPM: normalize by length then scale per million
rpk <- sweep(mat, 1, len/1000, "/")
tpm <- apply(rpk, 2, function(x) (x / sum(x)) * 1e6)

write.table(cpm, file=file.path(opt$output, "counts_CPM.tsv"), sep="\t", quote=FALSE, col.names=NA)
write.table(rpkm, file=file.path(opt$output, "counts_RPKM.tsv"), sep="\t", quote=FALSE, col.names=NA)
write.table(tpm, file=file.path(opt$output, "counts_TPM.tsv"), sep="\t", quote=FALSE, col.names=NA)

cat(sprintf("✅ Normalization done. Outputs: %s\n", opt$output))
