#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: generate_isoform_report.R <sqanti_report> <output_dir>")

report <- read.csv(args[1], sep="\t")
out <- args[2]
dir.create(out, showWarnings=FALSE, recursive=TRUE)

write.csv(table(report$structural_category), file.path(out, "category_counts.csv"))

pdf(file.path(out, "isoform_types.pdf"))
barplot(table(report$structural_category), main="Isoform Types (SQANTI3)")
dev.off()

cat("📊 Isoform summary exported →", out, "\n")
