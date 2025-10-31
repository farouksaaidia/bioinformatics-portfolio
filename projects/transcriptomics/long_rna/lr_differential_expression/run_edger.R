#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(edgeR)
  library(optparse)
})

option_list <- list(
  make_option(c("-c", "--counts"), type="character", help="Raw counts matrix (CSV/TSV)"),
  make_option(c("-m", "--metadata"), type="character", help="Sample metadata file (TSV: sample, condition)"),
  make_option(c("-o", "--outdir"), type="character", help="Output directory")
)
opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$counts) || is.null(opt$metadata) || is.null(opt$outdir))
  stop("Provide --counts, --metadata, and --outdir")

dir.create(opt$outdir, showWarnings=FALSE, recursive=TRUE)

counts <- read.delim(opt$counts, row.names=1)
meta <- read.delim(opt$metadata, header=TRUE)

group <- factor(meta$condition)
y <- DGEList(counts=counts, group=group)
y <- calcNormFactors(y)
y <- estimateDisp(y)
fit <- glmQLFit(y)
qlf <- glmQLFTest(fit)

res <- topTags(qlf, n=Inf)
write.csv(as.data.frame(res), file.path(opt$outdir, "edger_results.csv"))

cat("✅ edgeR analysis complete. Results saved in", opt$outdir, "\n")
