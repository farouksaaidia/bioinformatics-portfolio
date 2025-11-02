#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(ggplot2)
  library(optparse)
  library(dplyr)
})

option_list <- list(
  make_option(c("-i", "--indir"), type="character", help="Directory containing enrichment results (GO, KEGG, Reactome, GSEA)"),
  make_option(c("-o", "--outdir"), type="character", help="Output directory for report")
)
opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$indir) || is.null(opt$outdir)) stop("Provide --indir and --outdir")

dir.create(opt$outdir, showWarnings=FALSE, recursive=TRUE)

files <- list.files(opt$indir, pattern="_results.csv$", full.names=TRUE)
plots <- list()

for (f in files) {
  df <- read.csv(f)
  df <- head(df[order(df$p.adjust), ], 10)
  p <- ggplot(df, aes(x=reorder(Description, -p.adjust), y=-log10(p.adjust))) +
    geom_bar(stat="identity", fill="#3182bd") +
    coord_flip() +
    labs(title=basename(f), x="Pathway / Term", y="-log10(Adjusted p-value)") +
    theme_minimal(base_size=12)
  plots[[f]] <- p
  ggsave(filename=file.path(opt$outdir, paste0(basename(f), "_top10.pdf")), plot=p, width=7, height=5)
}

cat("✅ Enrichment summary report generated at", opt$outdir, "\n")
