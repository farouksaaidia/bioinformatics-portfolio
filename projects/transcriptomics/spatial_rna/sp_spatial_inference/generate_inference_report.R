#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(ggplot2)
  library(igraph)
  library(ggraph)
})

option_list <- list(
  make_option(c("-i","--input"), type="character", help="Merged ligand–receptor CSV"),
  make_option(c("-o","--output_pdf"), type="character", help="Output PDF file")
)
opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$input) || is.null(opt$output_pdf)) stop("❌ Provide input and output")

df <- read.csv(opt$input)
pdf(opt$output_pdf, width=12, height=10)

if (all(c("ligand","receptor") %in% colnames(df))) {
  g <- graph_from_data_frame(df[,c("ligand","receptor")])
  p <- ggraph(g, layout="fr") +
    geom_edge_link(alpha=0.4) +
    geom_node_point(size=3, color="steelblue") +
    geom_node_text(aes(label=name), repel=TRUE, size=3) +
    ggtitle("Ligand–Receptor Network (Consensus)") +
    theme_void()
  print(p)
} else {
  plot(1,1,type="n",axes=FALSE)
  text(1,1,"No valid ligand–receptor pairs found")
}

dev.off()
cat("✅ Inference network report saved ->", opt$output_pdf, "\n")
