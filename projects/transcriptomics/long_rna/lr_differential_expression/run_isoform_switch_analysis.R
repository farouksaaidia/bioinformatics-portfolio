#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(IsoformSwitchAnalyzeR)
  library(optparse)
})

option_list <- list(
  make_option(c("-a", "--abundance"), type="character", help="Isoform abundance file (TPM or counts)"),
  make_option(c("-m", "--metadata"), type="character", help="Sample metadata (TSV: sample, condition)"),
  make_option(c("-g", "--annotation"), type="character", help="GTF file for annotation"),
  make_option(c("-o", "--outdir"), type="character", help="Output directory")
)
opt <- parse_args(OptionParser(option_list=option_list))

if (is.null(opt$abundance) || is.null(opt$metadata) || is.null(opt$annotation) || is.null(opt$outdir))
  stop("Provide --abundance, --metadata, --annotation, and --outdir")

dir.create(opt$outdir, showWarnings=FALSE, recursive=TRUE)

meta <- read.delim(opt$metadata)
designMatrix <- data.frame(sampleID=meta$sample, condition=meta$condition)

aSwitchList <- importIsoformExpression(opt$abundance)
aSwitchList <- importRdata(aSwitchList, designMatrix=designMatrix, isoformAnnotation=opt$annotation)
aSwitchList <- isoformSwitchTestDEXSeq(aSwitchList)

write.csv(extractSwitchSummary(aSwitchList), file.path(opt$outdir, "isoform_switch_summary.csv"))
cat("✅ IsoformSwitchAnalyzeR complete. Results in", opt$outdir, "\n")
