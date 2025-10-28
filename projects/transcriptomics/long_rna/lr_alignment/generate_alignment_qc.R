#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Rsamtools)
  library(ggplot2)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: generate_alignment_qc.R <bam> <outdir>")

bam <- args[1]
out <- args[2]
dir.create(out, showWarnings=FALSE, recursive=TRUE)

idxstats <- idxstatsBam(bam) %>% as.data.frame()
write.csv(idxstats, file.path(out, "alignment_idxstats.csv"), row.names=FALSE)

pdf(file.path(out, "mapping_stats.pdf"))
barplot(idxstats$V2, names.arg=idxstats$V1, las=2, main="Mapped Reads per Chromosome")
dev.off()

cat("✅ Alignment QC complete →", out, "\n")
