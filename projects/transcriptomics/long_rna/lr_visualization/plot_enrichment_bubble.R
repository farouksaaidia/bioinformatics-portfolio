#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(ggplot2)
  library(optparse)
  library(readr)
  library(dplyr)
})

option_list <- list(
  make_option(c("-i", "--input"), type="character", help="Enrichment results file (GO/KEGG/Reactome)"),
  make_option(c("-o", "--outdir"), type="character", help="Output directory")
)
opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$input) || is.null(opt$outdir)) stop("Provide --input and --outdir")

dir.create(opt$outdir, showWarnings=FALSE, recursive=TRUE)
df <- read.csv(opt$input)
df <- head(df[order(df$p.adjust), ], 30)

p <- ggplot(df, aes(x=GeneRatio, y=reorder(Description, GeneRatio), size=Count, color=-log10(p.adjust))) +
  geom_point(alpha=0.8) +
  scale_color_gradient(low="skyblue", high="red") +
  labs(title="Top Enriched Terms", x="Gene Ratio", y="Term", color="-log10(p.adj)") +
  theme_minimal(base_size=12)
ggsave(file.path(opt$outdir, "enrichment_bubble_plot.pdf"), p, width=8, height=6)
cat("✅ Enrichment bubble plot saved to", opt$outdir, "\n")
