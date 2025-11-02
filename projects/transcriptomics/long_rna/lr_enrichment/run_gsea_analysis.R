#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(optparse)
  library(dplyr)
  library(readr)
})

option_list <- list(
  make_option(c("-i", "--input"), type="character", help="Input DE results file (CSV with gene_id, log2FoldChange, padj)"),
  make_option(c("-o", "--outdir"), type="character", help="Output directory for GSEA results"),
  make_option(c("-m", "--minGSSize"), type="integer", default=10, help="Minimum gene set size"),
  make_option(c("-M", "--maxGSSize"), type="integer", default=500, help="Maximum gene set size")
)
opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$input) || is.null(opt$outdir)) stop("Provide --input and --outdir")

dir.create(opt$outdir, showWarnings=FALSE, recursive=TRUE)
de <- read.csv(opt$input)
de <- de[!is.na(de$log2FoldChange), ]
geneList <- sort(setNames(de$log2FoldChange, de$gene_id), decreasing = TRUE)

gsea_res <- gseGO(geneList = geneList, OrgDb = org.Hs.eg.db, ont = "BP", keyType = "SYMBOL",
                  minGSSize = opt$minGSSize, maxGSSize = opt$maxGSSize, pAdjustMethod = "BH")
write.csv(as.data.frame(gsea_res), file.path(opt$outdir, "gsea_results.csv"))
cat("✅ GSEA complete. Results saved to", opt$outdir, "\n")
